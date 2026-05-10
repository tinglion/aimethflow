#!/usr/bin/env python3
"""
temporal_utils.py — Step 17: Temporal Prospective Validation.

Train/test split by sample collection date (default cutoff: 2024-07-01).
Evaluates the same 11-model lineup as Step 15 (5 cfDNA-only + CA-125 alone +
5 cfDNA + CA-125 fusion) under leakage-free conditions: feature selection,
imputation, and CA-125 standardisation are all fitted exclusively on the
training partition.

Outputs:
  temporal_validation_summary.csv      — AUC / Acc / Sens / Spec + 95% bootstrap CI
  temporal_validation_predictions.csv  — per-sample predictions
  ca125_distribution_by_period.{png,csv}  — diagnostic for selection bias
  temporal_delong_pairwise.csv         — formal pairwise DeLong p-values

Extracted from classifier_utils.py to keep that module focused on the main
classify step.
"""

import os
import warnings

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

from utils.logger import logger
from scripts.modules.classifier_utils import (
    _load_feature_matrix,
    _select_features_in_fold,
    _impute_fit,
    _impute_transform,
    _bootstrap_ci,
    _norm_id,
    _build_sid_to_case_map,
    _build_sample_date_map,
    _load_borderline_features,
)


def run_temporal_validation(workdir, result_dir, config, benign_samples, malignant_samples,
                            borderline_samples=None):
    """
    Supplementary temporal split validation (plasma mode only).

    Reads data/plasma_exp.xlsx (columns: case_id, datetime) joined with sample_info_excel
    to map case_id → sample_id.
    Splits at 2024-07-01: datetime < cutoff → train; datetime >= cutoff → validation.
    Trains each classifier once on the train set and evaluates on the validation set
    using the pre-selected panel features (panel_cpg.csv).

    Outputs:
        temporal_validation_predictions.csv — per-sample predictions for validation set
        temporal_validation_summary.csv     — AUC / accuracy / sensitivity / specificity
    """

    tv_cfg = (config or {}).get("temporal_validate", {})
    CUTOFF = pd.Timestamp(tv_cfg.get("cutoff_date", "2024-07-01"))
    TOP_N  = int(tv_cfg.get("top_n_features", 20))

    # 1-2. Build {sample_id: collection_date} via plasma_exp.xlsx + sample_info_excel bridge
    date_map = _build_sample_date_map(config)
    if not date_map:
        return None, None

    # 3. Load full feature matrix (all features retained until feature selection below)
    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    # 4. Apply exclude_samples (consistent with run_classifiers)
    ml_cfg = (config or {}).get("ml", {})
    exclude_set = set(ml_cfg.get("exclude_samples", []))
    if exclude_set:
        keep_mask = [sid not in exclude_set for sid in sample_ids]
        X          = X[[i for i, k in enumerate(keep_mask) if k]]
        y          = y[[i for i, k in enumerate(keep_mask) if k]]
        sample_ids = [sid for sid, k in zip(sample_ids, keep_mask) if k]

    # 5. Partition by collection date — must happen BEFORE feature selection so that
    #    validation samples never influence the feature set (true temporal leakage-free split)
    train_idx, test_idx = [], []
    no_date_count = 0
    for i, sid in enumerate(sample_ids):
        if sid not in date_map:
            no_date_count += 1
            continue
        if date_map[sid] < CUTOFF:
            train_idx.append(i)
        else:
            test_idx.append(i)

    if no_date_count:
        logger.info(f"[temporal] {no_date_count} sample(s) not in plasma_exp.xlsx — excluded")
    logger.info(f"[temporal] Split: {len(train_idx)} train (before {CUTOFF.date()})  "
                f"{len(test_idx)} validation (on/after {CUTOFF.date()})")

    if len(train_idx) < 5:
        logger.warning(f"[temporal] Insufficient train samples ({len(train_idx)}) — skipping")
        return None, None
    if len(test_idx) == 0:
        logger.warning("[temporal] No validation samples — skipping")
        return None, None

    # 6. Feature selection — two modes controlled by temporal_validate.reselect_features:
    #
    #   reselect_features = true  (default, leakage-free):
    #     Run Wilcoxon + LASSO∪RF on training-period samples only.
    #     Validation samples are never seen during feature selection.
    #
    #   reselect_features = false  (legacy, minor leakage):
    #     Use top-N features from feature_selection_freq.csv (computed by classify step on
    #     all samples including the validation-period ones — introduces minor temporal leakage).
    #
    # Blood-contaminated sites are excluded in both modes when panel_cpg_blood_filtered.csv exists.

    # Build blood whitelist (used by both paths below)
    blood_whitelist = None
    blood_filtered_path = os.path.join(result_dir, "panel_cpg_blood_filtered.csv")
    if os.path.exists(blood_filtered_path):
        blood_panel_df = pd.read_csv(blood_filtered_path)
        blood_whitelist = set(blood_panel_df["feature"].tolist())
    else:
        logger.debug(
            "[temporal] panel_cpg_blood_filtered.csv not found — "
            "run --step external_validate to enable blood-site exclusion"
        )

    reselect = tv_cfg.get("reselect_features", True)

    if reselect:
        # Leakage-free: select features from training-period samples only.
        # Disable the internal max_features_clf cap so _select_features_in_fold returns
        # the full ranked candidate set; we apply the cap ourselves after blood filtering,
        # ensuring max_feat clean features rather than fewer due to filter shrinkage.
        train_idx_arr = np.array(train_idx)
        max_feat = ml_cfg.get("max_features_clf", 20)
        _cfg_no_cap = {**config, "ml": {**ml_cfg, "max_features_clf": 0}}
        sel_idx = list(_select_features_in_fold(
            X[train_idx_arr], y[train_idx_arr], region_in_88, _cfg_no_cap
        ))
        if not sel_idx:
            sel_idx = list(range(X.shape[1]))
        # Blood filter: remove contaminated sites, then trim to max_feat from the clean set.
        if blood_whitelist is not None:
            n_before = len(sel_idx)
            sel_idx  = [i for i in sel_idx if feature_names[i] in blood_whitelist]
            n_removed = n_before - len(sel_idx)
            if n_removed:
                logger.info(
                    f"[temporal] Blood-filter: {n_removed} contaminated site(s) removed "
                    f"→ {len(sel_idx)} clean features remaining"
                )
        if max_feat > 0 and len(sel_idx) > max_feat:
            sel_idx = sel_idx[:max_feat]
        if not sel_idx:
            logger.warning("[temporal] No features remain after selection — skipping")
            return None, None
        X = X[:, np.array(sel_idx)]
        feature_names = [feature_names[i] for i in sel_idx]
        logger.info(
            f"[temporal] Leakage-free feature selection on training data: "
            f"{len(sel_idx)} features selected"
        )
    else:
        # Legacy: top-N from classify-step frequency table (minor temporal leakage)
        freq_path = os.path.join(result_dir, "feature_selection_freq.csv")
        if not os.path.exists(freq_path):
            logger.warning(
                f"[temporal] feature_selection_freq.csv not found — run classify first: {freq_path}"
            )
            return None, None
        freq_df = pd.read_csv(freq_path)
        top_features = set(freq_df["feature"].iloc[:TOP_N].tolist())
        if blood_whitelist is not None:
            n_before = len(top_features)
            top_features = top_features & blood_whitelist
            n_removed = n_before - len(top_features)
            if n_removed:
                logger.info(
                    f"[temporal] Blood-filter applied: {n_removed} site(s) removed "
                    f"from top features → {len(top_features)} retained"
                )
        top_idx = [i for i, fn in enumerate(feature_names) if fn in top_features]
        if not top_idx:
            logger.warning("[temporal] None of the top features found in feature matrix — skipping")
            return None, None
        X = X[:, top_idx]
        feature_names = [feature_names[i] for i in top_idx]
        logger.info(
            f"[temporal] Using top {len(top_idx)} features from feature_selection_freq.csv "
            f"(NOTE: computed on all samples — minor temporal leakage)"
        )

    X_train, y_train = X[train_idx], y[train_idx]
    X_test,  y_test  = X[test_idx],  y[test_idx]
    val_sample_ids  = [sample_ids[i] for i in test_idx]
    val_true_labels = ["malignant" if y[i] == 1 else "benign" for i in test_idx]
    test_dates      = [date_map.get(sample_ids[i]) for i in test_idx]

    # 6b. Optionally augment validation set with borderline samples (treated as benign / y=0).
    #     Borderline samples are NOT in the feature matrix; their bedGraph files are read on the fly.
    #     They are never added to the train set — only to the validation set.
    if tv_cfg.get("borderline_as_benign", False):
        X_bl, sids_bl, dates_bl = _load_borderline_features(
            workdir, borderline_samples or [], feature_names, date_map, CUTOFF
        )
        if len(sids_bl) > 0:
            if exclude_set:
                keep_bl  = [i for i, s in enumerate(sids_bl) if s not in exclude_set]
                X_bl     = X_bl[keep_bl]
                sids_bl  = [sids_bl[i]  for i in keep_bl]
                dates_bl = [dates_bl[i] for i in keep_bl]
            if len(sids_bl) > 0:
                X_test          = np.vstack([X_test, X_bl])
                y_test          = np.concatenate([y_test, np.zeros(len(sids_bl), dtype=int)])
                val_sample_ids  = val_sample_ids  + sids_bl
                val_true_labels = val_true_labels + ["borderline"] * len(sids_bl)
                test_dates      = test_dates      + dates_bl
                logger.info(
                    f"[temporal] borderline_as_benign: {len(sids_bl)} borderline sample(s) "
                    f"added to validation as benign (y=0): {sids_bl}"
                )

    # 7. Impute then scale — both fitted on train only to prevent leakage
    X_train_imp, col_medians = _impute_fit(X_train)
    X_test_imp               = _impute_transform(X_test, col_medians)

    scaler     = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train_imp)
    X_test_sc  = scaler.transform(X_test_imp)

    # 7b. Look up CA-125 for train/test samples (for baseline + combined models).
    # Same join pattern as Step 15: sid → case_id → ca125_map.
    include_ca125 = bool(ml_cfg.get("include_ca125_baseline", True))
    ca125_thresh  = float(ml_cfg.get("ca125_threshold", 35.0))
    ca125_train, ca125_test = None, None
    if include_ca125:
        try:
            from methylation.sample_grouping import load_ca125_map
            mode_cfg    = (config or {}).get("plasma", {})
            excel_sheet = mode_cfg.get("excel_sheet", "plasma")
            ca125_map   = load_ca125_map(category=excel_sheet)
            if ca125_map:
                sid_to_case = _build_sid_to_case_map(config)
                def _ca125_for(sid):
                    cid = sid_to_case.get(sid) or sid
                    return ca125_map.get(_norm_id(cid), np.nan)
                train_sids = [sample_ids[i] for i in train_idx]
                ca125_train = np.array([_ca125_for(s) for s in train_sids],     dtype=float)
                ca125_test  = np.array([_ca125_for(s) for s in val_sample_ids], dtype=float)
                n_train_ca = int((~np.isnan(ca125_train)).sum())
                n_test_ca  = int((~np.isnan(ca125_test)).sum())
                logger.info(f"[temporal] CA-125: train {n_train_ca}/{len(ca125_train)}, "
                            f"test {n_test_ca}/{len(ca125_test)} samples have measurements")
                if n_test_ca == 0:
                    logger.warning("[temporal] no test samples have CA-125 — skipping baseline")
                    ca125_train, ca125_test = None, None
            else:
                logger.info("[temporal] CA-125 baseline disabled: no CA-125 column in sample Excel")
        except Exception as exc:
            logger.warning(f"[temporal] CA-125 baseline disabled ({exc})")
            ca125_train, ca125_test = None, None

    # 7. Classifiers (same hyperparameters as run_classifiers)
    n_neg = int(sum(y_train == 0))
    n_pos = int(sum(y_train == 1))
    spw   = n_neg / n_pos if n_pos > 0 else 1.0
    logger.info(f"[temporal] train n_neg={n_neg} n_pos={n_pos}")
    logger.info(f"[temporal] test  n_neg={int(sum(y_test == 0))} n_pos={int(sum(y_test == 1))}")

    classifiers = {
        "LR": LogisticRegression(
            C=ml_cfg.get("lr_C", 1.0),
            class_weight="balanced",
            max_iter=2000,
            random_state=42,
        ),
        "SVM": SVC(
            kernel="rbf",
            C=ml_cfg.get("svm_C", 1.0),
            probability=True,
            class_weight="balanced",
            random_state=42,
        ),
        "RF": RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "XGB": XGBClassifier(
            n_estimators=ml_cfg.get("xgb_n_estimators", 200),
            max_depth=ml_cfg.get("xgb_max_depth", 4),
            learning_rate=ml_cfg.get("xgb_learning_rate", 0.1),
            scale_pos_weight=spw,
            use_label_encoder=False,
            eval_metric="logloss",
            random_state=42,
            n_jobs=-1,
            verbosity=0,
        ),
        "LGBM": LGBMClassifier(
            n_estimators=ml_cfg.get("lgbm_n_estimators", 200),
            max_depth=ml_cfg.get("lgbm_max_depth", 4),
            learning_rate=ml_cfg.get("lgbm_learning_rate", 0.1),
            scale_pos_weight=spw,
            random_state=42,
            n_jobs=-1,
            verbosity=-1,
        ),
    }

    # 8. Train and evaluate
    results = {}

    for clf_name, clf in classifiers.items():
        clf.fit(X_train_sc, y_train)
        y_prob = clf.predict_proba(X_test_sc)[:, 1]
        y_pred = clf.predict(X_test_sc)

        acc  = float(accuracy_score(y_test, y_pred))
        cm   = confusion_matrix(y_test, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
        specificity = tn / (tn + fp) if (tn + fp) > 0 else float("nan")

        n_classes_test = len(set(y_test))
        if n_classes_test < 2:
            logger.warning(f"[temporal] {clf_name}: only one class in validation set — AUC undefined")
            auc = float("nan")
            ci_dict = {"auc": ("", ""), "accuracy": ("", ""),
                       "sensitivity": ("", ""), "specificity": ("", "")}
        else:
            auc = float(roc_auc_score(y_test, y_prob))
            # Bootstrap 1000-resample 95% CI — small n_val (typically ~22) so wide CIs expected
            ci_dict = _bootstrap_ci(y_test, y_prob, y_pred, n_bootstrap=1000)

        ci_lo, ci_hi = ci_dict["auc"]
        ci_str = f" [{ci_lo:.3f}–{ci_hi:.3f}]" if ci_lo != "" else ""
        logger.info(f"[temporal] {clf_name}: AUC={auc:.3f}{ci_str}  Acc={acc:.3f}  "
                    f"Sens={sensitivity:.3f}  Spec={specificity:.3f}")

        results[clf_name] = {
            "y_true": y_test,
            "y_prob": y_prob,
            "y_pred": y_pred,
            "auc": auc,
            "accuracy": acc,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "confusion_matrix": cm,
            "ci": ci_dict,
        }

    # 8b. CA-125 baseline (filter NaN, evaluate on samples with a measurement only)
    if ca125_test is not None:
        valid_ca = ~np.isnan(ca125_test)
        n_valid_ca = int(valid_ca.sum())
        if n_valid_ca >= 5 and len(np.unique(y_test[valid_ca])) >= 2:
            y_true_ca = y_test[valid_ca]
            y_prob_ca = ca125_test[valid_ca]
            y_pred_ca = (y_prob_ca >= ca125_thresh).astype(int)
            cm_ca = confusion_matrix(y_true_ca, y_pred_ca, labels=[0, 1])
            tn_, fp_, fn_, tp_ = cm_ca.ravel()
            sens_ca = tp_ / (tp_ + fn_) if (tp_ + fn_) > 0 else float("nan")
            spec_ca = tn_ / (tn_ + fp_) if (tn_ + fp_) > 0 else float("nan")
            auc_ca  = float(roc_auc_score(y_true_ca, y_prob_ca))
            ci_ca   = _bootstrap_ci(y_true_ca, y_prob_ca, y_pred_ca, n_bootstrap=1000)
            ci_lo, ci_hi = ci_ca["auc"]
            ci_str = f" [{ci_lo:.3f}–{ci_hi:.3f}]" if ci_lo != "" else ""
            logger.info(f"[temporal] CA-125: AUC={auc_ca:.3f}{ci_str}  "
                        f"Acc={accuracy_score(y_true_ca, y_pred_ca):.3f}  "
                        f"Sens={sens_ca:.3f}  Spec={spec_ca:.3f}  (n={n_valid_ca})")
            results["CA-125"] = {
                "y_true":           y_true_ca,
                "y_prob":           y_prob_ca,
                "y_pred":           y_pred_ca,
                "valid_mask":       valid_ca,                  # length = len(val_sample_ids)
                "auc":              auc_ca,
                "accuracy":         float(accuracy_score(y_true_ca, y_pred_ca)),
                "sensitivity":      sens_ca,
                "specificity":      spec_ca,
                "confusion_matrix": cm_ca,
                "ci":               ci_ca,
            }
        else:
            logger.warning(f"[temporal] CA-125 baseline skipped (n_valid={n_valid_ca}, "
                           f"unique_classes={len(np.unique(y_test[valid_ca]))})")

    # 8c. Combined models: cfDNA features + CA-125 column
    if ca125_train is not None and ca125_test is not None:
        train_med = np.nanmedian(ca125_train) if not np.all(np.isnan(ca125_train)) else 0.0
        ca125_train_imp = np.where(np.isnan(ca125_train), train_med, ca125_train)
        ca125_test_imp  = np.where(np.isnan(ca125_test),  train_med, ca125_test)
        ca_scaler = StandardScaler()
        ca125_train_sc = ca_scaler.fit_transform(ca125_train_imp.reshape(-1, 1))
        ca125_test_sc  = ca_scaler.transform(ca125_test_imp.reshape(-1, 1))
        X_train_aug = np.hstack([X_train_sc, ca125_train_sc])
        X_test_aug  = np.hstack([X_test_sc,  ca125_test_sc])

        # Re-instantiate (the classifiers in `classifiers` were already fitted on cfDNA-only).
        combined_classifiers = {
            "LR": LogisticRegression(C=ml_cfg.get("lr_C", 1.0), class_weight="balanced",
                                     max_iter=2000, random_state=42),
            "SVM": SVC(kernel="rbf", C=ml_cfg.get("svm_C", 1.0), probability=True,
                       class_weight="balanced", random_state=42),
            "RF": RandomForestClassifier(n_estimators=ml_cfg.get("rf_n_estimators", 200),
                                         class_weight="balanced", random_state=42, n_jobs=-1),
            "XGB": XGBClassifier(n_estimators=ml_cfg.get("xgb_n_estimators", 200),
                                 max_depth=ml_cfg.get("xgb_max_depth", 4),
                                 learning_rate=ml_cfg.get("xgb_learning_rate", 0.1),
                                 scale_pos_weight=spw, use_label_encoder=False,
                                 eval_metric="logloss", random_state=42, n_jobs=-1, verbosity=0),
            "LGBM": LGBMClassifier(n_estimators=ml_cfg.get("lgbm_n_estimators", 200),
                                   max_depth=ml_cfg.get("lgbm_max_depth", 4),
                                   learning_rate=ml_cfg.get("lgbm_learning_rate", 0.1),
                                   scale_pos_weight=spw, random_state=42, n_jobs=-1, verbosity=-1),
        }

        for base_name, clf in combined_classifiers.items():
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",
                                        message="X does not have valid feature names",
                                        category=UserWarning)
                clf.fit(X_train_aug, y_train)
            y_prob = clf.predict_proba(X_test_aug)[:, 1]
            y_pred = clf.predict(X_test_aug)
            cm_c = confusion_matrix(y_test, y_pred, labels=[0, 1])
            tn_, fp_, fn_, tp_ = cm_c.ravel()
            sens_c = tp_ / (tp_ + fn_) if (tp_ + fn_) > 0 else float("nan")
            spec_c = tn_ / (tn_ + fp_) if (tn_ + fp_) > 0 else float("nan")
            if len(np.unique(y_test)) >= 2:
                auc_c   = float(roc_auc_score(y_test, y_prob))
                ci_c    = _bootstrap_ci(y_test, y_prob, y_pred, n_bootstrap=1000)
            else:
                auc_c   = float("nan")
                ci_c    = {"auc": ("", ""), "accuracy": ("", ""),
                           "sensitivity": ("", ""), "specificity": ("", "")}
            ci_lo, ci_hi = ci_c["auc"]
            ci_str = f" [{ci_lo:.3f}–{ci_hi:.3f}]" if ci_lo != "" else ""
            logger.info(f"[temporal] {base_name}+CA125: AUC={auc_c:.3f}{ci_str}  "
                        f"Acc={accuracy_score(y_test, y_pred):.3f}  "
                        f"Sens={sens_c:.3f}  Spec={spec_c:.3f}")
            results[f"{base_name}+CA125"] = {
                "y_true":           y_test,
                "y_prob":           y_prob,
                "y_pred":           y_pred,
                "auc":              auc_c,
                "accuracy":         float(accuracy_score(y_test, y_pred)),
                "sensitivity":      sens_c,
                "specificity":      spec_c,
                "confusion_matrix": cm_c,
                "ci":               ci_c,
            }

    # 9. Save predictions
    pred_rows = []
    for clf_name, res in results.items():
        if clf_name == "CA-125":
            # CA-125 baseline was evaluated only on samples with a measurement
            valid = res["valid_mask"]
            sids   = [val_sample_ids[i]  for i in range(len(val_sample_ids))  if valid[i]]
            labels = [val_true_labels[i] for i in range(len(val_true_labels)) if valid[i]]
            dates  = [test_dates[i]      for i in range(len(test_dates))      if valid[i]]
        else:
            sids, labels, dates = val_sample_ids, val_true_labels, test_dates
        for i, sid in enumerate(sids):
            pred_rows.append({
                "sample_id":       sid,
                "true_label":      labels[i],
                "classifier":      clf_name,
                "prob_malignant":  round(float(res["y_prob"][i]), 4),
                "predicted":       "malignant" if res["y_pred"][i] == 1 else "benign",
                "correct":         bool(res["y_true"][i] == res["y_pred"][i]),
                "split":           "validation",
                "collection_date": dates[i].strftime("%Y%m%d") if dates[i] is not None else "",
            })
    pred_df = pd.DataFrame(pred_rows)
    pred_path = os.path.join(result_dir, "temporal_validation_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"[temporal] Predictions saved: {pred_path}")

    # 10. Save summary
    def _r(v, n=4):
        try:
            return round(float(v), n)
        except (ValueError, TypeError):
            return ""

    summary_rows = []
    for clf_name, res in results.items():
        # For CA-125 baseline n_validation reflects the matched subset; for cfDNA
        # / combined models it's the full validation cohort (NaN imputed for combined).
        n_val_used = int(len(res["y_true"]))
        ci = res.get("ci", {})
        auc_lo,  auc_hi  = ci.get("auc",         ("", ""))
        acc_lo,  acc_hi  = ci.get("accuracy",    ("", ""))
        sens_lo, sens_hi = ci.get("sensitivity", ("", ""))
        spec_lo, spec_hi = ci.get("specificity", ("", ""))
        summary_rows.append({
            "classifier":   clf_name,
            "n_train":      len(train_idx),
            "n_validation": n_val_used,
            "auc":          _r(res["auc"]),
            "auc_ci95_lo":  _r(auc_lo),
            "auc_ci95_hi":  _r(auc_hi),
            "accuracy":     _r(res["accuracy"]),
            "acc_ci95_lo":  _r(acc_lo),
            "acc_ci95_hi":  _r(acc_hi),
            "sensitivity":  _r(res["sensitivity"]),
            "sens_ci95_lo": _r(sens_lo),
            "sens_ci95_hi": _r(sens_hi),
            "specificity":  _r(res["specificity"]),
            "spec_ci95_lo": _r(spec_lo),
            "spec_ci95_hi": _r(spec_hi),
            "cutoff_date":  str(CUTOFF.date()),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "temporal_validation_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"[temporal] Summary saved: {summary_path}")

    # 11. Console summary table — quick visual to see whether *+CA125 vs base CIs overlap
    logger.info("\n  ┌────────────────────────────────────────────────────────────────────────────┐")
    logger.info("  │  Temporal Validation Summary  (95% bootstrap CI in brackets)              │")
    logger.info("  ├─────────────┬──────────────────────┬───────────┬───────────┬───────────────┤")
    logger.info("  │ Classifier  │        AUC           │  Accuracy │  Sensitiv │  Specificity  │")
    logger.info("  ├─────────────┼──────────────────────┼───────────┼───────────┼───────────────┤")
    for row in summary_rows:
        def _fmt(v, lo, hi):
            if lo != "" and not (isinstance(v, str) and v == ""):
                return f"{v:.3f}[{lo:.3f}-{hi:.3f}]"
            return f"{v}" if isinstance(v, str) else f"{v:.3f}"
        logger.info(
            f"  │ {row['classifier']:<11} │ {_fmt(row['auc'], row['auc_ci95_lo'], row['auc_ci95_hi']):<20} │"
            f" {_fmt(row['accuracy'], row['acc_ci95_lo'], row['acc_ci95_hi']):<9} │"
            f" {_fmt(row['sensitivity'], row['sens_ci95_lo'], row['sens_ci95_hi']):<9} │"
            f"   {_fmt(row['specificity'], row['spec_ci95_lo'], row['spec_ci95_hi']):<13} │"
        )
    logger.info("  └─────────────┴──────────────────────┴───────────┴───────────┴───────────────┘")
    logger.info(f"  Train n={len(train_idx)} (before {CUTOFF.date()}), "
                f"Validation n_total={len(test_idx)} (on/after {CUTOFF.date()})")
    logger.info("  ⚠ With small validation cohorts (n≈22), wide bootstrap CIs are expected.")
    logger.info("    Check whether the cfDNA-only and *+CA125 CIs OVERLAP before claiming the "
                "fusion model improves or degrades vs cfDNA alone — a non-overlapping CI is the "
                "minimum bar for a real difference.")

    # 12. Diagnostic plot — CA-125 distribution by (period × group). Tests the
    # selection-bias hypothesis: if validation-period malignant CA-125 is markedly
    # higher than training-period malignant CA-125, the val cohort is enriched for
    # advanced/clinically-obvious cases.
    if ca125_train is not None and ca125_test is not None:
        try:
            _plot_ca125_distribution_by_period(
                ca125_train, y_train,
                ca125_test,  y_test,
                ca125_thresh,
                CUTOFF.date().isoformat(),
                result_dir,
            )
        except Exception as exc:
            logger.warning(f"[temporal] CA-125 distribution plot failed: {exc}")

    # 13. Pairwise DeLong tests across all 11 models (only on models sharing the
    # same y_test — cfDNA-only and *+CA125 share full val cohort; CA-125 alone may
    # be on a subset, so it's compared separately on its own subset).
    try:
        _delong_pairwise_temporal(results, result_dir)
    except Exception as exc:
        logger.warning(f"[temporal] Pairwise DeLong test failed: {exc}")

    return results, summary_df


def _plot_ca125_distribution_by_period(ca125_train, y_train, ca125_test, y_test,
                                       cutoff_value, cutoff_date_str, result_dir):
    """
    Box + strip plot of CA-125 by (period × group).
    4 boxes: Train-benign, Train-malignant, Val-benign, Val-malignant.
    Log y-axis (CA-125 is right-skewed). Cutoff line at clinical threshold.
    Saved to result_dir/ca125_distribution_by_period.png.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "DejaVu Sans"})

    # Filter NaN; build (period, group, value) tuples
    groups = []
    for arr, y, period in [(ca125_train, y_train, "Train"), (ca125_test, y_test, "Validation")]:
        for v, lab in zip(arr, y):
            if not np.isnan(v):
                groups.append((period, "Malignant" if int(lab) == 1 else "Benign", float(v)))

    if not groups:
        logger.info("[temporal] No CA-125 values for distribution plot — skipping")
        return

    df = pd.DataFrame(groups, columns=["period", "group", "ca125"])
    layout = [("Train", "Benign"), ("Train", "Malignant"),
              ("Validation", "Benign"), ("Validation", "Malignant")]
    data = [df[(df["period"] == p) & (df["group"] == g)]["ca125"].values
            for p, g in layout]
    labels = [f"{p}\n{g}\n(n={len(d)})" for (p, g), d in zip(layout, data)]
    n_total = sum(len(d) for d in data)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    bplot = ax.boxplot(data, labels=labels, widths=0.55, showfliers=False, patch_artist=True)
    palette = ["#4FA3D1", "#D67676", "#4FA3D1", "#D67676"]
    for patch, c in zip(bplot["boxes"], palette):
        patch.set_facecolor(c)
        patch.set_alpha(0.45)

    # Stripplot (jittered points) on top
    rng = np.random.default_rng(42)
    for i, d in enumerate(data, start=1):
        if len(d) == 0:
            continue
        x = rng.normal(loc=i, scale=0.06, size=len(d))
        ax.scatter(x, d, color=palette[i - 1], edgecolor="black",
                   linewidth=0.5, s=22, alpha=0.85, zorder=3)

    ax.axhline(cutoff_value, color="grey", lw=1, ls="--", alpha=0.7,
               label=f"Clinical cutoff = {cutoff_value:.0f} U/mL")
    ax.set_yscale("log")
    ax.set_ylabel("CA-125 (U/mL, log scale)")
    ax.set_title(f"CA-125 distribution by period × group  (cutoff: {cutoff_date_str},  n={n_total})")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.3, axis="y", which="both")

    # Annotate medians
    for i, d in enumerate(data, start=1):
        if len(d) > 0:
            med = float(np.median(d))
            ax.text(i, med, f"  {med:.0f}", va="center", fontsize=8,
                    color="black", fontweight="bold")

    plt.tight_layout()
    out_path = os.path.join(result_dir, "ca125_distribution_by_period.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"[temporal] CA-125 distribution plot saved: {out_path}")

    # Also dump summary statistics
    stats_rows = []
    for (p, g), d in zip(layout, data):
        if len(d) > 0:
            stats_rows.append({
                "period": p, "group": g, "n": len(d),
                "median":   round(float(np.median(d)), 2),
                "mean":     round(float(np.mean(d)),   2),
                "q25":      round(float(np.percentile(d, 25)), 2),
                "q75":      round(float(np.percentile(d, 75)), 2),
                "min":      round(float(np.min(d)),    2),
                "max":      round(float(np.max(d)),    2),
                "frac_above_cutoff": round(float((d >= cutoff_value).mean()), 3),
            })
    stats_path = os.path.join(result_dir, "ca125_distribution_by_period.csv")
    pd.DataFrame(stats_rows).to_csv(stats_path, index=False)
    logger.info(f"[temporal] CA-125 distribution stats: {stats_path}")


def _delong_pairwise_temporal(results, result_dir):
    """
    Pairwise DeLong AUC comparison across all temporal-validation models.

    Models with mismatched sample sets (CA-125 alone, when its valid_mask
    excludes some samples) are compared separately on their own subset.
    Otherwise (typical case: all 22 val samples have CA-125) a single
    pairwise matrix covers all 11 models.
    """
    from scripts.modules.imeta_utils import _fastdelong_cov
    from scipy import stats as _stats

    if not results:
        return

    # Group models by their y_true tuple — models sharing the same evaluation set
    # can be compared in a single fastDeLong matrix.
    groups: dict = {}
    for clf_name, res in results.items():
        key = tuple(np.asarray(res["y_true"], dtype=int).tolist())
        groups.setdefault(key, []).append(clf_name)

    rows = []
    for y_key, clf_list in groups.items():
        if len(clf_list) < 2:
            continue
        y_true = np.array(y_key, dtype=int)
        if len(np.unique(y_true)) < 2:
            continue
        score_matrix = np.vstack([results[c]["y_prob"] for c in clf_list])
        aucs, cov = _fastdelong_cov(y_true, score_matrix)

        for i, name_a in enumerate(clf_list):
            for j, name_b in enumerate(clf_list):
                if j <= i:
                    continue
                delta = aucs[i] - aucs[j]
                var_d = cov[i, i] + cov[j, j] - 2 * cov[i, j]
                if var_d <= 0:
                    z = float("nan")
                    p = float("nan")
                else:
                    se = float(np.sqrt(var_d))
                    z  = float(delta / se)
                    p  = float(2 * _stats.norm.sf(abs(z)))
                rows.append({
                    "classifier_A":  name_a,
                    "classifier_B":  name_b,
                    "n_eval":        int(len(y_true)),
                    "auc_A":         round(float(aucs[i]), 4),
                    "auc_B":         round(float(aucs[j]), 4),
                    "delta_auc":     round(float(delta),   4),
                    "se_delta":      round(float(np.sqrt(max(var_d, 0))), 4) if var_d > 0 else "",
                    "z_stat":        round(z, 3) if not np.isnan(z) else "",
                    "p_delong":      round(p, 4) if not np.isnan(p) else "",
                    "significant":   bool(p < 0.05) if not np.isnan(p) else False,
                })

    if not rows:
        logger.info("[temporal] No pairs for DeLong test — skipping")
        return

    df = pd.DataFrame(rows).sort_values("p_delong")
    out_path = os.path.join(result_dir, "temporal_delong_pairwise.csv")
    df.to_csv(out_path, index=False, encoding="utf-8")
    logger.info(f"[temporal] Pairwise DeLong saved: {out_path}")

    # Highlight significant pairs in the log
    sig = df[df["significant"] == True]
    if len(sig) > 0:
        logger.info(f"  ── Significantly different pairs (p<0.05) — {len(sig)} of {len(df)} ──")
        for _, row in sig.iterrows():
            logger.info(f"    {row['classifier_A']:<12} vs {row['classifier_B']:<12} "
                        f"Δ={row['delta_auc']:+.3f}  p={row['p_delong']:.4f}")
    else:
        logger.info(f"  No significantly different pairs (all p ≥ 0.05, n_pairs={len(df)})")


# ─────────────────────────────────────────────────────────────────────────────
