#!/usr/bin/env python3
"""
fs_benchmark.py — Step 15b: Feature selection ablation study.

Runs CV across four feature-selection modes (wilcoxon / lasso / rf / union) on
the same cohort, classifiers, and CV folds, so the marginal contribution of each
selection strategy can be compared in isolation.

Output: result/feature_selection_benchmark.csv

Extracted from classifier_utils.py to keep that module focused on the main
classify step. All shared helpers (_load_feature_matrix, _build_classifiers,
_select_features_in_fold, _impute_fit/_transform) are imported from
classifier_utils.
"""

import copy
import os
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from utils.logger import logger
from scripts.modules.classifier_utils import (
    _load_feature_matrix,
    _build_classifiers,
    _select_features_in_fold,
    _impute_fit,
    _impute_transform,
)


def run_feature_selection_benchmark(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Ablation study: run CV for each feature_selection_mode and compare AUC.

    Modes tested: wilcoxon | lasso | rf | union
    Reuses the existing feature matrix (feature_matrix.npz) — no re-alignment needed.
    Saves result/feature_selection_benchmark.csv.

    Returns: pd.DataFrame with columns
        [mode, classifier, n_features_mean, auc, accuracy, sensitivity, specificity]
    """
    logger.info("━━ Feature Selection Benchmark (ablation) ━━")

    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )
    if X.shape[1] == 0:
        logger.error("[benchmark] Feature matrix empty — aborting")
        return pd.DataFrame()

    ml_cfg_base = (config or {}).get("ml", {})
    exclude_set = set(ml_cfg_base.get("exclude_samples", []))
    if exclude_set:
        keep = np.array([sid not in exclude_set for sid in sample_ids])
        X, y = X[keep], y[keep]
        sample_ids = [s for s, k in zip(sample_ids, keep) if k]

    n_neg = int((y == 0).sum())
    n_pos = int((y == 1).sum())
    spw   = n_neg / n_pos if n_pos > 0 else 1.0

    cv_folds = ml_cfg_base.get("cv_folds", 5)
    if cv_folds <= 0:
        splitter = LeaveOneOut()
    else:
        splitter = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)

    modes = ["wilcoxon", "lasso", "rf", "union"]
    all_rows = []

    for mode in modes:
        logger.info(f"  [benchmark] mode={mode} ...")
        cfg_mode = copy.deepcopy(config)
        cfg_mode.setdefault("ml", {})["feature_selection_mode"] = mode

        classifiers = _build_classifiers(ml_cfg_base, spw)
        clf_names   = list(classifiers.keys())
        # Note: this benchmark only computes aggregate metrics (AUC/Acc/Sens/Spec)
        # which are permutation-invariant, so appending in fold order is safe.
        # If you ever extend this to write per-sample predictions, switch to the
        # pre-allocated test_idx-indexed pattern used in run_classifiers().
        clf_y_true  = {n: [] for n in clf_names}
        clf_y_prob  = {n: [] for n in clf_names}
        clf_y_pred  = {n: [] for n in clf_names}
        fold_n_feat = []

        for train_idx, test_idx in splitter.split(X, y):
            X_tr_raw, X_te_raw = X[train_idx], X[test_idx]
            y_tr, y_te         = y[train_idx], y[test_idx]

            sel_idx = _select_features_in_fold(X_tr_raw, y_tr, region_in_88, cfg_mode)
            if len(sel_idx) == 0:
                sel_idx = np.arange(X.shape[1])
            fold_n_feat.append(len(sel_idx))

            X_tr = X_tr_raw[:, sel_idx]
            X_te = X_te_raw[:, sel_idx]
            X_tr_imp, meds = _impute_fit(X_tr)
            X_te_imp       = _impute_transform(X_te, meds)
            scaler         = StandardScaler()
            X_tr_sc        = scaler.fit_transform(X_tr_imp)
            X_te_sc        = scaler.transform(X_te_imp)

            for clf_name, clf in classifiers.items():
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        clf.fit(X_tr_sc, y_tr)
                    prob = clf.predict_proba(X_te_sc)[:, 1]
                    pred = clf.predict(X_te_sc)
                    clf_y_true[clf_name].extend(y_te.tolist())
                    clf_y_prob[clf_name].extend(prob.tolist())
                    clf_y_pred[clf_name].extend(pred.tolist())
                except Exception as exc:
                    logger.debug(f"[benchmark] {mode}/{clf_name} fold error: {exc}")

        n_feat_mean = int(round(np.mean(fold_n_feat))) if fold_n_feat else 0

        for clf_name in clf_names:
            y_true = np.array(clf_y_true[clf_name])
            y_prob = np.array(clf_y_prob[clf_name])
            y_pred = np.array(clf_y_pred[clf_name])
            if len(y_true) == 0 or len(np.unique(y_true)) < 2:
                continue
            auc  = roc_auc_score(y_true, y_prob)
            acc  = accuracy_score(y_true, y_pred)
            cm   = confusion_matrix(y_true, y_pred, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel()
            sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            all_rows.append({
                "mode":            mode,
                "classifier":      clf_name,
                "n_features_mean": n_feat_mean,
                "auc":             round(auc,  4),
                "accuracy":        round(acc,  4),
                "sensitivity":     round(sens, 4),
                "specificity":     round(spec, 4),
            })

        mode_aucs = {r["classifier"]: r["auc"] for r in all_rows if r["mode"] == mode}
        logger.info(
            f"    {mode}: AUCs = " +
            "  ".join(f"{clf}={v:.3f}" for clf, v in mode_aucs.items())
        )

    df = pd.DataFrame(all_rows)
    out_path = os.path.join(result_dir, "feature_selection_benchmark.csv")
    df.to_csv(out_path, index=False, encoding="utf-8")
    logger.info(f"  [benchmark] Saved: {out_path}")

    logger.info("\n  Feature Selection Benchmark — Best AUC per mode")
    logger.info("  ┌────────────┬──────────┬─────────────────┬────────┐")
    logger.info("  │ Mode       │ Best clf │ n_features_mean │  AUC   │")
    logger.info("  ├────────────┼──────────┼─────────────────┼────────┤")
    for mode in modes:
        mode_rows = df[df["mode"] == mode]
        if mode_rows.empty:
            continue
        best = mode_rows.loc[mode_rows["auc"].idxmax()]
        logger.info(
            f"  │ {best['mode']:<10} │ {best['classifier']:<8} │"
            f" {int(best['n_features_mean']):>15} │ {best['auc']:.3f}  │"
        )
    logger.info("  └────────────┴──────────┴─────────────────┴────────┘")

    return df
