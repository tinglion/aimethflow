#!/usr/bin/env python3
"""
classifier_utils.py - Diagnostic classifier module (Step 15)

Evaluates five classifiers with configurable cross-validation (ml.cv_folds in config.toml):
  cv_folds = 0  → Leave-One-Out CV (LOO-CV; default for very small cohorts)
  cv_folds = k  → Stratified k-fold CV (default k=5; recommended for n≥50)

  LR    — L2 logistic regression (interpretable, well-suited for small samples)
  SVM   — RBF-kernel SVM (non-linear decision boundary)
  RF    — Random Forest (ensemble method with built-in feature importance)
  XGB   — XGBoost (gradient-boosted trees; scale_pos_weight for class imbalance)
  LGBM  — LightGBM (efficient gradient boosting; scale_pos_weight for class imbalance)

Outputs:
  roc_curve.png              — ROC curve comparison across classifiers
  confusion_matrix.png       — confusion matrix per classifier
  classifier_predictions.csv — per-sample predicted probabilities and labels
  classifier_summary.csv     — AUC / accuracy / sensitivity / specificity summary

Dependencies: scikit-learn, matplotlib, xgboost, lightgbm
"""

import os
import warnings
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

from sklearn.model_selection import LeaveOneOut, StratifiedKFold

from utils.logger import logger


def _load_feature_matrix(workdir, result_dir, benign_samples, malignant_samples, config):
    """
    Load the feature matrix: prefer restoring from feature_matrix.npz (saved by
    feature_analysis step); rebuild from bedGraph if not found (requires workdir).

    Returns:
        X, y, sample_ids, feature_names, region_in_88
    """
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if os.path.exists(fm_path):
        logger.info(f"Loading feature matrix from cache: {fm_path}")
        data = np.load(fm_path, allow_pickle=True)
        X = data["X"]
        y = data["y"]
        sample_ids    = list(data["sample_ids"])
        feature_names = list(data["feature_names"])
        region_in_88  = data["region_in_88"] if "region_in_88" in data else np.zeros(X.shape[1], dtype=bool)
        # Drop borderline rows (y=2) — not used in classifiers or feature selection
        _bm = y != 2
        if not _bm.all():
            X, y = X[_bm], y[_bm]
            sample_ids   = [s for s, k in zip(sample_ids, _bm) if k]
            region_in_88 = region_in_88  # feature-axis, unchanged
    else:
        logger.info("feature_matrix.npz not found — rebuilding from bedGraph...")
        from scripts.modules.feature_utils import build_feature_matrix
        X, y, sample_ids, feature_names, _, region_in_88 = build_feature_matrix(
            workdir, result_dir, benign_samples, malignant_samples, config
        )
    return X, y, sample_ids, feature_names, region_in_88





def _apply_panel_filter(X, feature_names, region_in_88, result_dir):
    """
    Restrict the feature matrix to the deployment panel; prefer the
    blood-filtered variant when available.

    Lookup order:
      1. panel_cpg_blood_filtered.csv  (Step 19 external_validate output)
      2. panel_cpg.csv                  (Step 14 feature_analysis output)
      3. no panel → return X unchanged

    Returns
    -------
    (X_filt, feature_names_filt, region_in_88_filt, mode_label)
    """
    import pandas as pd
    blood_path = os.path.join(result_dir, "panel_cpg_blood_filtered.csv")
    panel_path = os.path.join(result_dir, "panel_cpg.csv")

    if os.path.exists(blood_path):
        chosen_path, mode_label = blood_path, "blood-filtered panel"
    elif os.path.exists(panel_path):
        chosen_path, mode_label = panel_path, "Step 14 panel (no blood filter)"
    else:
        return X, feature_names, region_in_88, "no panel restriction (full matrix)"

    try:
        panel_set = set(pd.read_csv(chosen_path)["feature"].astype(str))
    except Exception as exc:
        logger.warning(f"  Failed to read {os.path.basename(chosen_path)}: {exc} "
                       "— falling back to full matrix")
        return X, feature_names, region_in_88, "no panel restriction (read failed)"

    keep = np.array([fn in panel_set for fn in feature_names])
    if keep.sum() == 0:
        logger.warning(f"  Panel filter dropped all features — falling back to full matrix")
        return X, feature_names, region_in_88, "no panel restriction (empty intersection)"

    return (X[:, keep],
            [fn for fn, k in zip(feature_names, keep) if k],
            region_in_88[keep],
            mode_label)


def _impute(X):
    """Impute NaN values using column medians."""
    X_imp = X.copy()
    col_medians = np.nanmedian(X_imp, axis=0)
    for j in range(X_imp.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j] if not np.isnan(col_medians[j]) else 0.0
    return X_imp


def _impute_fit(X):
    """
    Impute NaN values with training-set column medians; return (imputed matrix, medians).
    Columns that are entirely NaN default to 0.0.
    """
    col_medians = np.nanmedian(X, axis=0)
    col_medians = np.where(np.isnan(col_medians), 0.0, col_medians)
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp, col_medians


def _impute_transform(X, col_medians):
    """Impute test-set NaN values using training-set column medians."""
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp


def _select_features_in_fold(X_train, y_train, region_in_88, config):
    """
    Run feature selection on the training set of a single LOO fold; return selected column indices.

    Pipeline: stratified Wilcoxon pre-filter (relaxed thresholds for 88-region sites,
    standard thresholds elsewhere) → LASSO → RF → union.
    Imputes NaN internally using training-set medians; test samples are never touched.

    Args:
        X_train:      (n_train, n_features)
        y_train:      (n_train,)
        region_in_88: (n_features,) bool — whether each feature is in the 88 candidate regions
        config:       config dict

    Returns:
        selected_indices: np.ndarray of column indices (relative to X_train)
    """
    try:
        from scipy import stats as _stats
        _HAS_SCIPY = True
    except ImportError:
        _HAS_SCIPY = False

    ml_cfg = (config or {}).get("ml", {})
    wilcox_p        = ml_cfg.get("sample_wilcox_p",       0.05)
    wilcox_delta    = ml_cfg.get("sample_delta_beta",     0.10)
    wilcox_p_88     = ml_cfg.get("sample_wilcox_p_88",    0.10)
    wilcox_delta_88 = ml_cfg.get("sample_delta_beta_88",  0.05)
    selection_mode  = ml_cfg.get("feature_selection_mode", "union")

    n_features = X_train.shape[1]

    # 1. impute (training-set medians)
    X_imp, _ = _impute_fit(X_train)

    # 2. stratified Wilcoxon pre-filter
    if _HAS_SCIPY:
        pvals  = np.ones(n_features)
        deltas = np.zeros(n_features)
        for j in range(n_features):
            b_vals = X_imp[y_train == 0, j]
            m_vals = X_imp[y_train == 1, j]
            if len(b_vals) >= 2 and len(m_vals) >= 2:
                _, p = _stats.mannwhitneyu(b_vals, m_vals, alternative="two-sided")
                d    = abs(float(np.mean(m_vals)) - float(np.mean(b_vals)))
            else:
                p, d = 1.0, 0.0
            pvals[j]  = p
            deltas[j] = d

        p_thresh = np.where(region_in_88, wilcox_p_88,    wilcox_p)
        d_thresh = np.where(region_in_88, wilcox_delta_88, wilcox_delta)
        wilcox_mask = (pvals < p_thresh) & (deltas >= d_thresh)

        if wilcox_mask.sum() == 0:   # fallback: relax thresholds
            wilcox_mask = (pvals < 0.20) & (deltas >= 0.05)
        if wilcox_mask.sum() == 0:
            wilcox_mask = np.ones(n_features, dtype=bool)
    else:
        wilcox_mask = np.ones(n_features, dtype=bool)

    wilcox_indices = np.where(wilcox_mask)[0]

    # "wilcoxon" mode: skip LASSO/RF, rank candidates by |Δβ| (if available), cap, return early
    if selection_mode == "wilcoxon":
        if _HAS_SCIPY:
            order = np.argsort(deltas[wilcox_indices])[::-1]
            ranked = wilcox_indices[order]
        else:
            ranked = wilcox_indices   # no delta scores available; keep original order
        max_feat = ml_cfg.get("max_features_clf", 0)
        if max_feat > 0 and len(ranked) > max_feat:
            ranked = ranked[:max_feat]
        return ranked

    X_w = X_imp[:, wilcox_indices]

    # 3. LASSO (in-fold speedup: Cs=10, cv≤3)
    cv_folds = min(ml_cfg.get("lasso_cv_folds", 5), 3)
    scaler = StandardScaler()
    X_ws = scaler.fit_transform(X_w)
    lasso = LogisticRegressionCV(
        Cs=10, cv=cv_folds, penalty="l1", solver="liblinear",
        class_weight="balanced", max_iter=1000, random_state=42,
    )
    lasso.fit(X_ws, y_train)
    lasso_local = np.where(lasso.coef_[0] != 0)[0]

    # 4. RF (100 trees, quick in-fold estimate)
    rf = RandomForestClassifier(
        n_estimators=100, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )
    rf.fit(X_w, y_train)
    rf_local = np.where(rf.feature_importances_ > 0)[0]

    # 5. union, ranked by combined LASSO coef + RF importance, capped at max_features
    if selection_mode == "lasso":
        local_idx = lasso_local
    elif selection_mode == "rf":
        local_idx = rf_local
    else:
        local_idx = np.unique(np.concatenate([lasso_local, rf_local]))

    # combined score (each normalized then summed): higher = higher priority
    lasso_abs  = np.abs(lasso.coef_[0])                  # shape: (n_wilcox,)
    rf_imp     = rf.feature_importances_                  # shape: (n_wilcox,)
    lasso_norm = lasso_abs / (lasso_abs.max() + 1e-10)
    rf_norm    = rf_imp    / (rf_imp.max()    + 1e-10)
    combined   = lasso_norm + rf_norm                     # shape: (n_wilcox,)

    max_feat = ml_cfg.get("max_features_clf", 0)
    if max_feat > 0 and len(local_idx) > max_feat:
        scores    = combined[local_idx]
        top_n     = np.argsort(scores)[-max_feat:]        # top-N by combined score
        local_idx = local_idx[np.sort(top_n)]             # preserve column order

    return wilcox_indices[local_idx]


def _bootstrap_ci(y_true, y_prob, y_pred, n_bootstrap=1000, ci=0.95, seed=42):
    """Bootstrap percentile CI for AUC / accuracy / sensitivity / specificity."""
    from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix as _cm
    rng = np.random.default_rng(seed)
    n = len(y_true)
    alpha = (1 - ci) / 2
    aucs, accs, sens, specs = [], [], [], []
    for _ in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        yt, yp, ypr = y_true[idx], y_pred[idx], y_prob[idx]
        if len(np.unique(yt)) < 2:
            continue
        aucs.append(roc_auc_score(yt, ypr))
        accs.append(accuracy_score(yt, yp))
        tn, fp, fn, tp = _cm(yt, yp, labels=[0, 1]).ravel()
        sens.append(tp / (tp + fn) if (tp + fn) > 0 else 0.0)
        specs.append(tn / (tn + fp) if (tn + fp) > 0 else 0.0)

    def _pct(arr):
        a = np.array(arr)
        return float(np.percentile(a, 100 * alpha)), float(np.percentile(a, 100 * (1 - alpha)))

    return {
        "auc":         _pct(aucs),
        "accuracy":    _pct(accs),
        "sensitivity": _pct(sens),
        "specificity": _pct(specs),
    }


def _save_roc_plot(results, result_dir):
    """Plot and save the classifier ROC comparison figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    matplotlib.rcParams.update({"font.family": "DejaVu Sans"})

    fig, ax = plt.subplots(figsize=(8, 7))
    colors = {
        "LR":     "#2196F3",
        "SVM":    "#4CAF50",
        "RF":     "#FF5722",
        "XGB":    "#9C27B0",
        "LGBM":   "#FF9800",
        "CA-125": "#555555",   # baseline serum biomarker — dashed grey
    }

    def _style(clf_name):
        """Pick (color, linestyle, linewidth) based on classifier kind."""
        if clf_name == "CA-125":
            return colors["CA-125"], "--", 1.5
        if clf_name.endswith("+CA125"):
            base = clf_name.split("+", 1)[0]
            return colors.get(base, "gray"), ":", 2.5     # combined: same colour, dotted, thicker
        return colors.get(clf_name, "gray"), "-", 2

    for clf_name, res in results.items():
        fpr, tpr, _ = roc_curve(res["y_true"], res["y_prob"])
        ci = res.get("ci", {})
        auc_ci = ci.get("auc", (None, None))
        if auc_ci[0] is not None:
            lbl = f"{clf_name} (AUC={res['auc']:.3f} [{auc_ci[0]:.3f}–{auc_ci[1]:.3f}])"
        else:
            lbl = f"{clf_name} (AUC={res['auc']:.3f})"
        c, ls, lw = _style(clf_name)
        ax.plot(fpr, tpr, label=lbl, color=c, lw=lw, linestyle=ls)

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate (1 − Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title("CV ROC Curves (Benign vs Malignant)")
    n_curves = len(results)
    ax.legend(loc="lower right",
              fontsize=8 if n_curves > 6 else 9,
              ncol=2 if n_curves > 7 else 1)
    ax.grid(alpha=0.3)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])

    roc_path = os.path.join(result_dir, "roc_curve.png")
    fig.savefig(roc_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  ROC curve saved: {roc_path}")
    return roc_path


def _save_confusion_matrix_plot(results, result_dir):
    """Plot a confusion matrix for each classifier; wraps to multiple rows beyond 3."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "DejaVu Sans"})

    n_clf = len(results)
    n_cols = min(n_clf, 3)
    n_rows = (n_clf + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 4 * n_rows),
                             squeeze=False)
    axes_flat = [axes[r][c] for r in range(n_rows) for c in range(n_cols)]

    for ax, (clf_name, res) in zip(axes_flat, results.items()):
        cm = res["confusion_matrix"]
        im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
        ax.set_title(f"{clf_name}\nAUC={res['auc']:.3f}")
        tick_labels = ["Benign", "Malignant"]
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(tick_labels)
        ax.set_yticklabels(tick_labels)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]),
                        ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black",
                        fontsize=14)

    # hide unused subplots
    for ax in axes_flat[n_clf:]:
        ax.set_visible(False)

    fig.tight_layout()
    cm_path = os.path.join(result_dir, "confusion_matrix.png")
    fig.savefig(cm_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  Confusion matrix saved: {cm_path}")
    return cm_path


def _build_classifiers(ml_cfg, spw):
    """Instantiate the 5 standard classifiers from config. spw = scale_pos_weight."""
    return {
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


# ──────────────────────────────────────────────────────
# Step 15 entry point
# ──────────────────────────────────────────────────────

def run_classifiers(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Step 15 entry point: Nested CV evaluation of LR / SVM / RF / XGB / LGBM classifiers.

    CV strategy is controlled by ml.cv_folds in config.toml:
      cv_folds = 0  → Leave-One-Out CV (maximally data-efficient, slow)
      cv_folds = k  → Stratified k-fold CV (default 5; faster, recommended for n≥50)

    Feature selection (Wilcoxon + LASSO + RF) is performed independently on each fold's
    training set; test samples are never used during feature selection or scaler fitting —
    no data leakage.

    Args:
        workdir:           working directory (contains bedGraph files; used only when
                           no cached feature matrix exists)
        result_dir:        results directory (reads feature_matrix.npz, writes outputs)
        config:            config dict
        benign_samples:    list of benign sample IDs
        malignant_samples: list of malignant sample IDs

    Returns:
        results:    dict {clf_name: {y_true, y_prob, y_pred, auc, accuracy,
                                    sensitivity, specificity, confusion_matrix}}
        summary_df: pd.DataFrame of per-classifier summary metrics
    """
    ml_cfg_cv = (config or {}).get("ml", {})
    cv_folds = ml_cfg_cv.get("cv_folds", 5)
    cv_name  = "LOO-CV" if cv_folds <= 0 else f"{cv_folds}-fold CV"

    logger.info(f"Starting classifier evaluation ({cv_name})  "
                f"(benign={len(benign_samples)}, malignant={len(malignant_samples)})")

    # 1. load full feature matrix (post-coverage filter, pre-Wilcoxon)
    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    if X.shape[1] == 0:
        logger.error("Feature matrix has no features — skipping classifier evaluation")
        return {}, pd.DataFrame()

    # 1b. outlier detection (skip if Step 13 already ran and wrote outlier_report.csv)
    ml_cfg_pre = (config or {}).get("ml", {})
    outlier_path = os.path.join(result_dir, "outlier_report.csv")
    if os.path.exists(outlier_path):
        logger.info(f"  Outlier report already exists (from Step 13) — skipping re-run: {outlier_path}")
    else:
        from scripts.modules.feature_utils import detect_sample_outliers
        outlier_df, _ = detect_sample_outliers(X, y, sample_ids, config)
        if len(outlier_df) > 0:
            outlier_df.to_csv(outlier_path, index=False)
            logger.info(f"  Outlier report saved: {outlier_path}")

    exclude_set = set(ml_cfg_pre.get("exclude_samples", []))
    if exclude_set:
        keep_mask = np.array([sid not in exclude_set for sid in sample_ids])
        n_excluded = int((~keep_mask).sum())
        logger.info(f"  Excluding {n_excluded} sample(s) per exclude_samples: "
                    f"{[sid for sid in sample_ids if sid in exclude_set]}")
        X          = X[keep_mask]
        y          = y[keep_mask]
        sample_ids = [sid for sid, k in zip(sample_ids, keep_mask) if k]
        region_in_88 = region_in_88   # feature dimension unchanged, no filtering needed

    # ── Panel filter: prefer blood-filtered panel (Step 19) when available ──
    # Restrict the candidate feature pool to the Step 14 panel; if Step 19
    # external_validate has run and written panel_cpg_blood_filtered.csv, that
    # variant is preferred so leukocyte-cfDNA confounded sites are excluded
    # from per-fold selection. Falls back to panel_cpg.csv if the blood-filter
    # is absent, or to the full feature matrix if no panel file exists.
    X, feature_names, region_in_88, panel_mode = _apply_panel_filter(
        X, feature_names, region_in_88, result_dir
    )

    n_in_88 = int(region_in_88.sum())
    logger.info(f"  Feature pool ({panel_mode}): {X.shape[0]} samples × {X.shape[1]} features  "
                f"(in 88-regions: {n_in_88}  outside: {X.shape[1]-n_in_88};  "
                f"feature selection runs independently per fold)")

    # 2. define classifiers
    ml_cfg = (config or {}).get("ml", {})

    n_neg = int(sum(y == 0))
    n_pos = int(sum(y == 1))
    spw   = n_neg / n_pos if n_pos > 0 else 1.0

    classifiers = _build_classifiers(ml_cfg, spw)

    # 2b. CA-125 baseline values aligned with sample_ids (NaN if no measurement).
    # Treated as a pseudo-classifier in the CV loop so AUC is computed on the
    # IDENTICAL test samples as cfDNA models — fair head-to-head per fold.
    include_ca125 = bool(ml_cfg.get("include_ca125_baseline", True))
    ca125_thresh  = float(ml_cfg.get("ca125_threshold", 35.0))   # standard clinical cutoff (U/mL)
    ca125_vals    = None
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
                ca125_vals = np.array([_ca125_for(s) for s in sample_ids], dtype=float)
                n_with_ca = int((~np.isnan(ca125_vals)).sum())
                logger.info(f"  CA-125 baseline: {n_with_ca}/{len(sample_ids)} samples have "
                            f"CA-125 values (clinical cutoff = {ca125_thresh:.1f} U/mL)")
                if n_with_ca == 0:
                    ca125_vals = None
            else:
                logger.info("  CA-125 baseline disabled: no CA-125 column in sample Excel")
        except Exception as exc:
            logger.warning(f"  CA-125 baseline disabled (load failure: {exc})")

    # 3. Nested CV (LOO or stratified k-fold)
    if cv_folds <= 0:
        splitter = LeaveOneOut()
        n_folds  = X.shape[0]
    else:
        splitter = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
        n_folds  = cv_folds
    logger.info(f"  CV strategy: {cv_name}  ({n_folds} folds)")

    clf_names = list(classifiers.keys())
    base_clf_names = list(classifiers.keys())   # for combined models — see below
    if ca125_vals is not None:
        clf_names.append("CA-125")    # treat as pseudo-classifier in storage + aggregation
        # Combined models: each base classifier trained on cfDNA features + CA-125 column.
        # NaN CA-125 imputed with training-fold median (leakage-free).
        clf_names.extend([f"{n}+CA125" for n in base_clf_names])

    # Per-classifier prediction storage, pre-allocated and indexed by ORIGINAL sample
    # position (NOT fold-iteration order). This is critical: StratifiedKFold(shuffle=True)
    # returns test_idx in shuffled order, so appending in fold order would misalign
    # predictions with sample_ids when later writing classifier_predictions.csv.
    n_samples  = X.shape[0]
    clf_y_true = {n: np.full(n_samples, -1, dtype=int)            for n in clf_names}
    clf_y_prob = {n: np.full(n_samples, np.nan, dtype=float)      for n in clf_names}
    clf_y_pred = {n: np.full(n_samples, -1, dtype=int)            for n in clf_names}

    fold_n_features    = []  # features selected per fold
    fold_selected_names = []  # selected feature names per fold (for frequency stats)

    log_interval = max(1, n_folds // 10)

    for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test         = y[train_idx], y[test_idx]

        # 3a. independent feature selection on training set (Wilcoxon + LASSO + RF)
        sel_idx = _select_features_in_fold(X_train_raw, y_train, region_in_88, config)
        if len(sel_idx) == 0:
            sel_idx = np.arange(X.shape[1])   # fallback: all features
        fold_n_features.append(len(sel_idx))
        fold_selected_names.append([feature_names[i] for i in sel_idx])

        X_tr = X_train_raw[:, sel_idx]
        X_te = X_test_raw[:, sel_idx]

        # 3b. impute (training-set medians)
        X_tr_imp, col_meds = _impute_fit(X_tr)
        X_te_imp = _impute_transform(X_te, col_meds)

        # 3c. standardize (fit on training set)
        scaler = StandardScaler()
        X_tr_sc = scaler.fit_transform(X_tr_imp)
        X_te_sc = scaler.transform(X_te_imp)

        # 3d. train and predict with all classifiers on this fold's feature set
        # LightGBM >= 4.x generates internal feature names when fitting numpy arrays;
        # sklearn's check layer emits a UserWarning about missing names — suppress locally.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="X does not have valid feature names",
                category=UserWarning,
            )
            for clf_name, clf in classifiers.items():
                clf.fit(X_tr_sc, y_train)
                prob = clf.predict_proba(X_te_sc)[:, 1]
                pred = clf.predict(X_te_sc)
                # Write back at ORIGINAL sample positions (test_idx) so predictions
                # stay aligned with sample_ids regardless of fold shuffling.
                clf_y_true[clf_name][test_idx] = y_test
                clf_y_prob[clf_name][test_idx] = prob
                clf_y_pred[clf_name][test_idx] = pred

        # CA-125 baseline: no training, just look up the value for each test sample.
        # Same-fold, same-test-set as the classifiers above → fair AUC comparison.
        # NaN samples (no CA-125 measurement) stay NaN and are filtered at aggregation.
        if ca125_vals is not None:
            ca_test = ca125_vals[test_idx]
            clf_y_true["CA-125"][test_idx] = y_test
            clf_y_prob["CA-125"][test_idx] = ca_test
            # Hard prediction for accuracy/sens/spec uses the clinical cutoff.
            ca_pred = np.where(np.isnan(ca_test), -1, (ca_test >= ca125_thresh).astype(int))
            clf_y_pred["CA-125"][test_idx] = ca_pred

            # ── Combined models: cfDNA features + CA-125 as an extra column ──
            # Leakage-free preprocessing: impute training-fold median, then fit a
            # separate StandardScaler on the CA-125 column (its scale is U/mL,
            # very different from beta values which are already standardized).
            ca_train = ca125_vals[train_idx]
            train_med = np.nanmedian(ca_train) if not np.all(np.isnan(ca_train)) else 0.0
            ca_train_imp = np.where(np.isnan(ca_train), train_med, ca_train)
            ca_test_imp  = np.where(np.isnan(ca_test),  train_med, ca_test)

            ca_scaler = StandardScaler()
            ca_train_sc = ca_scaler.fit_transform(ca_train_imp.reshape(-1, 1))
            ca_test_sc  = ca_scaler.transform(ca_test_imp.reshape(-1, 1))

            X_tr_aug = np.hstack([X_tr_sc, ca_train_sc])
            X_te_aug = np.hstack([X_te_sc, ca_test_sc])

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",
                                        message="X does not have valid feature names",
                                        category=UserWarning)
                # Need fresh classifier instances — the ones in `classifiers` were
                # already fitted on the cfDNA-only matrix above.
                combined_clfs = _build_classifiers(ml_cfg, spw)
                for clf_name, clf in combined_clfs.items():
                    clf.fit(X_tr_aug, y_train)
                    prob = clf.predict_proba(X_te_aug)[:, 1]
                    pred = clf.predict(X_te_aug)
                    combined_key = f"{clf_name}+CA125"
                    clf_y_true[combined_key][test_idx] = y_test
                    clf_y_prob[combined_key][test_idx] = prob
                    clf_y_pred[combined_key][test_idx] = pred

        if (fold_idx + 1) % log_interval == 0 or fold_idx == n_folds - 1:
            logger.info(f"    fold {fold_idx+1}/{n_folds}  "
                        f"n_features={len(sel_idx)}")

    # per-fold feature count stats
    fold_arr = np.array(fold_n_features)
    logger.info(f"  Nested CV features per fold: mean={fold_arr.mean():.1f}  "
                f"median={float(np.median(fold_arr)):.1f}  "
                f"min={fold_arr.min()}  max={fold_arr.max()}")

    # feature selection frequency table
    from collections import Counter
    name_counter = Counter()
    for names in fold_selected_names:
        name_counter.update(names)

    feat_to_in88 = {fn: bool(region_in_88[i]) for i, fn in enumerate(feature_names)}
    freq_rows = sorted(
        [{"feature": feat, "in_88_regions": feat_to_in88.get(feat, False),
          "selection_count": cnt, "selection_rate": round(cnt / n_folds, 4)}
         for feat, cnt in name_counter.items()],
        key=lambda r: -r["selection_count"]
    )
    freq_df = pd.DataFrame(freq_rows)
    freq_path = os.path.join(result_dir, "feature_selection_freq.csv")
    freq_df.to_csv(freq_path, index=False)
    logger.info(f"  Feature selection frequency table saved: {freq_path}")

    # stable features (selection rate ≥ 50%)
    stable_thr = 0.50
    stable_rows = [r for r in freq_rows if r["selection_rate"] >= stable_thr]
    logger.info(f"\n  ══ Stable feature panel (rate≥{stable_thr:.0%}, {len(stable_rows)} features) ══")
    for r in stable_rows:
        marker = "★" if r["in_88_regions"] else "○"
        logger.info(f"    {marker} {r['feature']:<30s}  "
                    f"{r['selection_rate']:.0%} ({r['selection_count']}/{n_folds} folds)")
    if not stable_rows:
        logger.info("    (no stable features — consider checking max_features_clf or relaxing Wilcoxon thresholds)")
    logger.info("")

    # 4. aggregate metrics
    # Sanity check: every sample must have been a test sample exactly once.
    # If any slot is still the -1 sentinel, fold coverage is broken — surface loudly.
    # CA-125 is special: NaN-prob slots are samples WITHOUT a CA-125 measurement,
    # not coverage gaps; they're filtered out before metric computation.
    # Both "CA-125" alone and the "*+CA125" combined models may have NaN slots
    # for samples without a CA-125 measurement — those are filtered, not errors.
    ca_dependent = {"CA-125", *(f"{n}+CA125" for n in base_clf_names)}

    results = {}
    for clf_name in clf_names:
        y_true_full = clf_y_true[clf_name]
        y_prob_full = clf_y_prob[clf_name]
        y_pred_full = clf_y_pred[clf_name]

        if clf_name in ca_dependent:
            valid = ~np.isnan(y_prob_full)
            n_valid = int(valid.sum())
            if n_valid < 5 or len(np.unique(y_true_full[valid])) < 2:
                logger.warning(f"  [{clf_name}] only {n_valid} samples with valid CA-125 — skipping")
                continue
            y_true = y_true_full[valid]
            y_prob = y_prob_full[valid]
            y_pred = y_pred_full[valid]
        else:
            unfilled = int((y_true_full == -1).sum())
            if unfilled > 0:
                logger.error(f"  [{clf_name}] {unfilled} sample(s) never appeared as a test sample — "
                             "CV fold coverage is incomplete. Predictions for these samples will be wrong.")
            y_true, y_prob, y_pred = y_true_full, y_prob_full, y_pred_full

        auc = roc_auc_score(y_true, y_prob)
        acc = accuracy_score(y_true, y_pred)
        cm  = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        ci_dict = _bootstrap_ci(y_true, y_prob, y_pred)

        results[clf_name] = {
            "y_true": y_true,                # filtered (used for AUC + ROC plot)
            "y_prob": y_prob,
            "y_pred": y_pred,
            "y_true_full": y_true_full,      # length n_samples (used for per-sample CSV)
            "y_prob_full": y_prob_full,
            "y_pred_full": y_pred_full,
            "n_used": int(len(y_true)),      # how many samples contributed to the metric
            "auc": auc,
            "accuracy": acc,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "confusion_matrix": cm,
            "ci": ci_dict,
        }
        auc_lo, auc_hi = ci_dict["auc"]
        logger.info(f"    {clf_name}: AUC={auc:.3f} [{auc_lo:.3f}–{auc_hi:.3f}]  "
                    f"Acc={acc:.3f}  Sens={sensitivity:.3f}  Spec={specificity:.3f}")

    # 5. save ROC curve
    _save_roc_plot(results, result_dir)

    # 6. save confusion matrix plots
    _save_confusion_matrix_plot(results, result_dir)

    # 7. save per-sample predictions  (uses y_*_full → length matches sample_ids)
    rows = []
    for clf_name, res in results.items():
        full_true = res["y_true_full"]
        full_prob = res["y_prob_full"]
        full_pred = res["y_pred_full"]
        is_ca_dep = clf_name in ca_dependent
        for i, sid in enumerate(sample_ids):
            # Skip CA-125-dependent rows for samples without a CA-125 measurement.
            if is_ca_dep and np.isnan(full_prob[i]):
                continue
            rows.append({
                "sample_id":      sid,
                "true_label":     "malignant" if full_true[i] == 1 else "benign",
                "classifier":     clf_name,
                "prob_malignant": round(float(full_prob[i]), 4),
                "predicted":      "malignant" if full_pred[i] == 1 else "benign",
                "correct":        bool(full_true[i] == full_pred[i]),
            })
    pred_df = pd.DataFrame(rows)
    pred_path = os.path.join(result_dir, "classifier_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"  Per-sample predictions saved: {pred_path}")

    # 8. save summary metrics
    summary_rows = []
    for clf_name, res in results.items():
        cm = res["confusion_matrix"]
        tn, fp, fn, tp = cm.ravel()
        ci = res.get("ci", {})
        auc_lo,  auc_hi  = ci.get("auc",         (None, None))
        acc_lo,  acc_hi  = ci.get("accuracy",     (None, None))
        sens_lo, sens_hi = ci.get("sensitivity",  (None, None))
        spec_lo, spec_hi = ci.get("specificity",  (None, None))
        # n_benign / n_malignant reflect what was actually used. For cfDNA classifiers
        # that's the full cohort; for CA-125 it's only samples with a measurement.
        y_used = res["y_true"]
        summary_rows.append({
            "classifier":        clf_name,
            "auc":               round(res["auc"],         4),
            "auc_ci95_lo":       round(auc_lo,  4) if auc_lo  is not None else "",
            "auc_ci95_hi":       round(auc_hi,  4) if auc_hi  is not None else "",
            "accuracy":          round(res["accuracy"],    4),
            "acc_ci95_lo":       round(acc_lo,  4) if acc_lo  is not None else "",
            "acc_ci95_hi":       round(acc_hi,  4) if acc_hi  is not None else "",
            "sensitivity":       round(res["sensitivity"], 4),
            "sens_ci95_lo":      round(sens_lo, 4) if sens_lo is not None else "",
            "sens_ci95_hi":      round(sens_hi, 4) if sens_hi is not None else "",
            "specificity":       round(res["specificity"], 4),
            "spec_ci95_lo":      round(spec_lo, 4) if spec_lo is not None else "",
            "spec_ci95_hi":      round(spec_hi, 4) if spec_hi is not None else "",
            "TP": int(tp), "FP": int(fp),
            "TN": int(tn), "FN": int(fn),
            "n_used":            int(len(y_used)),
            "n_benign":          int((y_used == 0).sum()),
            "n_malignant":       int((y_used == 1).sum()),
            "n_features_mean":   round(float(fold_arr.mean()), 1) if clf_name != "CA-125" else "",
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "classifier_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  Classifier summary saved: {summary_path}")

    # print summary table
    logger.info("\n  ┌────────────────────────────────────────────────────────────────────────────┐")
    logger.info("  │  Classifier CV Summary  (95% Bootstrap CI in brackets)                    │")
    logger.info("  ├─────────────┬──────────────────────┬───────────┬───────────┬───────────────┤")
    logger.info("  │ Classifier  │        AUC           │  Accuracy │  Sensitiv │  Specificity  │")
    logger.info("  ├─────────────┼──────────────────────┼───────────┼───────────┼───────────────┤")
    for row in summary_rows:
        def _fmt(v, lo, hi):
            if lo != "":
                return f"{v:.3f}[{lo:.3f}-{hi:.3f}]"
            return f"{v:.3f}"
        logger.info(
            f"  │ {row['classifier']:<11} │ {_fmt(row['auc'], row['auc_ci95_lo'], row['auc_ci95_hi']):<20} │"
            f" {_fmt(row['accuracy'], row['acc_ci95_lo'], row['acc_ci95_hi']):<9} │"
            f" {_fmt(row['sensitivity'], row['sens_ci95_lo'], row['sens_ci95_hi']):<9} │"
            f"   {_fmt(row['specificity'], row['spec_ci95_lo'], row['spec_ci95_hi']):<13} │"
        )
    logger.info("  └─────────────┴──────────────────────┴───────────┴───────────┴───────────────┘")

    return results, summary_df


def _norm_id(raw):
    """Normalise Excel case_id: '640849.0' → '640849'; non-numeric strings kept as-is."""
    s = str(raw).strip()
    try:
        return str(int(float(s)))
    except (ValueError, OverflowError):
        return s


def _build_sid_to_case_map(config):
    """
    Build {sample_id: case_id} by reading the sample info Excel.

    sample_id = fq_R1 filename prefix (e.g. 'B10001'), same key used in
    feature_matrix.npz and classifier_predictions.csv.
    case_id   = normalized value of the Excel 'case_id' column
    (Excel often stores integers as floats; _norm_id handles '640849.0' → '640849').

    Returns an empty dict on any failure.
    """
    excel_file = (config or {}).get("data", {}).get("sample_info_excel", "")
    if not excel_file or not os.path.exists(excel_file):
        logger.warning(f"[id-map] sample_info_excel not found: {excel_file!r}")
        return {}
    sid_to_case = {}
    try:
        xl = pd.ExcelFile(excel_file)
        for sheet in xl.sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet)
            if "fq_R1" not in df.columns or "case_id" not in df.columns:
                continue
            for _, row in df.iterrows():
                cid = _norm_id(row.get("case_id", ""))
                r1  = str(row.get("fq_R1", "")).strip()
                sid = r1.split(".")[0]
                if cid and cid != "nan" and sid and sid != "nan":
                    sid_to_case[sid] = cid
    except Exception as exc:
        logger.warning(f"[id-map] Failed to read sample_info_excel: {exc}")
    logger.debug(f"[id-map] {len(sid_to_case)} sample_id → case_id mappings built")
    return sid_to_case


def _build_sample_date_map(config):
    """
    Build {sample_id: collection_date} keyed by fq_R1 prefix.

    Two-step join via _build_sid_to_case_map:
      plasma_exp.xlsx  : case_id → collection_date
      sample_info_excel: sample_id (fq_R1 prefix) → case_id
    """
    # Step 1: plasma_exp.xlsx → {case_id: timestamp}
    from methylation.sample_grouping import _resolve_excel
    exp_filename = (config or {}).get("data", {}).get("plasma_exp_excel", "plasma_exp.xlsx")
    exp_path = _resolve_excel(exp_filename) if exp_filename else None
    if not exp_path or not os.path.exists(exp_path):
        logger.warning(
            f"[temporal] plasma collection-date Excel '{exp_filename}' not found "
            f"(searched data/{exp_filename}, {exp_filename}) — skipping temporal validation"
        )
        return {}

    try:
        df_exp = pd.read_excel(exp_path)
    except Exception as exc:
        logger.warning(f"[temporal] Failed to read {exp_path}: {exc}")
        return {}

    df_exp.columns = df_exp.columns.str.strip().str.lower()
    if "case_id" not in df_exp.columns or "datetime" not in df_exp.columns:
        logger.warning("[temporal] plasma_exp.xlsx missing required columns (case_id, datetime)")
        return {}

    case_date = {}
    for _, row in df_exp.iterrows():
        cid    = _norm_id(row["case_id"])
        dt_raw = str(row["datetime"]).strip()
        ts     = pd.to_datetime(dt_raw, format="%Y%m%d", errors="coerce")
        if pd.isna(ts):
            ts = pd.to_datetime(dt_raw, errors="coerce")
        if not pd.isna(ts):
            case_date[cid] = ts

    if not case_date:
        logger.warning("[temporal] No valid dates parsed from plasma_exp.xlsx")
        return {}

    # Step 2: sample_id → case_id from sample_info_excel
    sid_to_case = _build_sid_to_case_map(config)
    if not sid_to_case:
        return {}

    # Step 3: join
    sample_date_map = {
        sid: case_date[cid]
        for sid, cid in sid_to_case.items()
        if cid in case_date
    }
    logger.info(
        f"[temporal] date map: {len(case_date)} case_ids in plasma_exp.xlsx, "
        f"{len(sid_to_case)} sample_ids in sample_info_excel, "
        f"{len(sample_date_map)} matched → {list(sample_date_map.keys())[:5]}"
    )
    return sample_date_map


def _load_borderline_features(workdir, borderline_samples, feature_names, date_map, cutoff):
    """
    Load per-sample feature vectors for borderline-category plasma samples whose
    collection date falls in the validation period (>= cutoff).

    Returns:
        X_bl    : np.ndarray (n, len(feature_names)), methylation 0-1, NaN = not covered
        sids_bl : list[str]  sample IDs
        dates_bl: list[pd.Timestamp]
    """
    _empty = (np.empty((0, len(feature_names))), [], [])

    if not borderline_samples:
        logger.info("[temporal/borderline] No borderline samples in pipeline state")
        return _empty

    # Filter to validation period (>= cutoff); log those excluded
    test_sids  = [s for s in borderline_samples if s in date_map and date_map[s] >= cutoff]
    n_train_bl = sum(1 for s in borderline_samples if s in date_map and date_map[s] < cutoff)
    n_nodate   = sum(1 for s in borderline_samples if s not in date_map)
    if n_nodate:
        logger.info(f"[temporal/borderline] {n_nodate} borderline sample(s) not in date map — excluded")
    if n_train_bl:
        logger.info(f"[temporal/borderline] {n_train_bl} borderline sample(s) in train period — excluded from validation")
    if not test_sids:
        logger.info("[temporal/borderline] No borderline samples fall in the validation period")
        return _empty

    # Parse feature_names ("chrom:pos") into (chrom, int_pos) lookup keys
    feat_keys = []
    for fn in feature_names:
        parts = fn.rsplit(":", 1)
        if len(parts) == 2:
            try:
                feat_keys.append((parts[0], int(parts[1])))
            except ValueError:
                feat_keys.append(None)
        else:
            feat_keys.append(None)

    # Read each borderline bedGraph and extract values for panel features
    X_rows, sids_out, dates_out = [], [], []
    for sid in test_sids:
        bg_path = os.path.join(workdir, f"{sid}_CpG.bedGraph")
        if not os.path.exists(bg_path):
            logger.warning(f"[temporal/borderline] bedGraph not found: {bg_path} — skipping {sid}")
            continue
        try:
            bg = pd.read_csv(bg_path, sep="\t", header=None,
                             names=["chrom", "start", "end", "meth_pct",
                                    "meth_count", "unmeth_count"])
            bg = bg[pd.to_numeric(bg["start"], errors="coerce").notna()].copy()
            bg["start"] = bg["start"].astype(int)
            site_map = {(row.chrom, row.start): float(row.meth_pct)
                        for row in bg.itertuples(index=False)}
        except Exception as e:
            logger.warning(f"[temporal/borderline] Failed to read bedGraph for {sid}: {e}")
            continue

        row_vals = np.full(len(feature_names), np.nan)
        for j, key in enumerate(feat_keys):
            if key is not None and key in site_map:
                row_vals[j] = site_map[key]
        X_rows.append(row_vals)
        sids_out.append(sid)
        dates_out.append(date_map[sid])

    if not X_rows:
        return _empty

    logger.info(
        f"[temporal/borderline] Added {len(X_rows)} borderline sample(s) to validation benign: {sids_out}"
    )
    return np.array(X_rows), sids_out, dates_out


