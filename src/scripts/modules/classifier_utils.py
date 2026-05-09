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
    feature_select step); rebuild from bedGraph if not found (requires workdir).

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

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {
        "LR":   "#2196F3",
        "SVM":  "#4CAF50",
        "RF":   "#FF5722",
        "XGB":  "#9C27B0",
        "LGBM": "#FF9800",
    }

    for clf_name, res in results.items():
        fpr, tpr, _ = roc_curve(res["y_true"], res["y_prob"])
        ci = res.get("ci", {})
        auc_ci = ci.get("auc", (None, None))
        if auc_ci[0] is not None:
            lbl = f"{clf_name} (AUC={res['auc']:.3f} [{auc_ci[0]:.3f}–{auc_ci[1]:.3f}])"
        else:
            lbl = f"{clf_name} (AUC={res['auc']:.3f})"
        ax.plot(fpr, tpr, label=lbl, color=colors.get(clf_name, "gray"), lw=2)

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate (1 − Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title("LOO-CV ROC Curves (Benign vs Malignant)")
    ax.legend(loc="lower right")
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

    n_in_88 = int(region_in_88.sum())
    logger.info(f"  Full feature matrix: {X.shape[0]} samples × {X.shape[1]} features  "
                f"(in 88-regions: {n_in_88}  outside: {X.shape[1]-n_in_88};  "
                f"feature selection runs independently per fold)")

    # 2. define classifiers
    ml_cfg = (config or {}).get("ml", {})

    n_neg = int(sum(y == 0))
    n_pos = int(sum(y == 1))
    spw   = n_neg / n_pos if n_pos > 0 else 1.0

    classifiers = _build_classifiers(ml_cfg, spw)

    # 3. Nested CV (LOO or stratified k-fold)
    if cv_folds <= 0:
        splitter = LeaveOneOut()
        n_folds  = X.shape[0]
    else:
        splitter = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42)
        n_folds  = cv_folds
    logger.info(f"  CV strategy: {cv_name}  ({n_folds} folds)")

    clf_names = list(classifiers.keys())

    # per-classifier prediction storage
    clf_y_true = {n: [] for n in clf_names}
    clf_y_prob = {n: [] for n in clf_names}
    clf_y_pred = {n: [] for n in clf_names}

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
                clf_y_true[clf_name].extend(y_test.tolist())
                clf_y_prob[clf_name].extend(prob.tolist())
                clf_y_pred[clf_name].extend(pred.tolist())

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
    results = {}
    for clf_name in clf_names:
        y_true = np.array(clf_y_true[clf_name])
        y_prob = np.array(clf_y_prob[clf_name])
        y_pred = np.array(clf_y_pred[clf_name])

        auc = roc_auc_score(y_true, y_prob)
        acc = accuracy_score(y_true, y_pred)
        cm  = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        ci_dict = _bootstrap_ci(y_true, y_prob, y_pred)

        results[clf_name] = {
            "y_true": y_true,
            "y_prob": y_prob,
            "y_pred": y_pred,
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

    # 7. save per-sample predictions
    rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(sample_ids):
            rows.append({
                "sample_id": sid,
                "true_label": "malignant" if res["y_true"][i] == 1 else "benign",
                "classifier": clf_name,
                "prob_malignant": round(float(res["y_prob"][i]), 4),
                "predicted": "malignant" if res["y_pred"][i] == 1 else "benign",
                "correct": bool(res["y_true"][i] == res["y_pred"][i]),
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
            "n_benign":          int(sum(y == 0)),
            "n_malignant":       int(sum(y == 1)),
            "n_features_mean":   round(float(fold_arr.mean()), 1),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "classifier_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  Classifier summary saved: {summary_path}")

    # print summary table
    logger.info("\n  ┌───────────────────────────────────────────────────────────────────────┐")
    logger.info("  │  Classifier Nested LOO-CV Summary  (95% Bootstrap CI in brackets)    │")
    logger.info("  ├────────┬──────────────────────┬───────────┬───────────┬───────────────┤")
    logger.info("  │  Clf   │        AUC           │  Accuracy │  Sensitiv │  Specificity  │")
    logger.info("  ├────────┼──────────────────────┼───────────┼───────────┼───────────────┤")
    for row in summary_rows:
        def _fmt(v, lo, hi):
            if lo != "":
                return f"{v:.3f}[{lo:.3f}-{hi:.3f}]"
            return f"{v:.3f}"
        logger.info(
            f"  │ {row['classifier']:<6} │ {_fmt(row['auc'], row['auc_ci95_lo'], row['auc_ci95_hi']):<20} │"
            f" {_fmt(row['accuracy'], row['acc_ci95_lo'], row['acc_ci95_hi']):<9} │"
            f" {_fmt(row['sensitivity'], row['sens_ci95_lo'], row['sens_ci95_hi']):<9} │"
            f"   {_fmt(row['specificity'], row['spec_ci95_lo'], row['spec_ci95_hi']):<13} │"
        )
    logger.info("  └────────┴──────────────────────┴───────────┴───────────┴───────────────┘")

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


def run_feature_selection_benchmark(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Ablation study: run 5-fold CV for each feature_selection_mode and compare AUC.

    Modes tested: wilcoxon | lasso | rf | union
    Reuses the existing feature matrix (feature_matrix.npz) — no re-alignment needed.
    Saves result/feature_selection_benchmark.csv.

    Returns: pd.DataFrame with columns
        [mode, classifier, n_features_mean, auc, accuracy, sensitivity, specificity]
    """
    import copy

    logger.info("━━ Feature Selection Benchmark (ablation) ━━")

    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )
    if X.shape[1] == 0:
        logger.error("[benchmark] Feature matrix empty — aborting")
        return pd.DataFrame()

    ml_cfg_base  = (config or {}).get("ml", {})
    exclude_set  = set(ml_cfg_base.get("exclude_samples", []))
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

    # Print summary table
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
        else:
            auc = float(roc_auc_score(y_test, y_prob))

        logger.info(f"[temporal] {clf_name}: AUC={auc:.3f}  Acc={acc:.3f}  "
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
        }

    # 9. Save predictions
    pred_rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(val_sample_ids):
            pred_rows.append({
                "sample_id":       sid,
                "true_label":      val_true_labels[i],
                "classifier":      clf_name,
                "prob_malignant":  round(float(res["y_prob"][i]), 4),
                "predicted":       "malignant" if res["y_pred"][i] == 1 else "benign",
                "correct":         bool(res["y_true"][i] == res["y_pred"][i]),
                "split":           "validation",
                "collection_date": test_dates[i].strftime("%Y%m%d") if test_dates[i] is not None else "",
            })
    pred_df = pd.DataFrame(pred_rows)
    pred_path = os.path.join(result_dir, "temporal_validation_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"[temporal] Predictions saved: {pred_path}")

    # 10. Save summary
    summary_rows = []
    for clf_name, res in results.items():
        summary_rows.append({
            "classifier":   clf_name,
            "n_train":      len(train_idx),
            "n_validation": len(test_idx),
            "auc":          round(res["auc"], 4) if not np.isnan(res["auc"]) else "",
            "accuracy":     round(res["accuracy"], 4),
            "sensitivity":  round(res["sensitivity"], 4) if not np.isnan(res["sensitivity"]) else "",
            "specificity":  round(res["specificity"], 4) if not np.isnan(res["specificity"]) else "",
            "cutoff_date":  str(CUTOFF.date()),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "temporal_validation_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"[temporal] Summary saved: {summary_path}")

    return results, summary_df


# ─────────────────────────────────────────────────────────────────────────────
# Step 15c: CA-125 vs classifier AUC comparison
# ─────────────────────────────────────────────────────────────────────────────

def run_ca125_comparison(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Compare CA-125 AUC with classifier AUCs on the same matched sample subset.

    Reads CA-125 values from the sample info Excel (column 'ca125' / 'CA-125').
    If no CA-125 data is available, logs a warning and returns (None, None).

    Outputs
    -------
    result/ca125_comparison.csv  — AUC ± 95% CI for CA-125 and each classifier
                                   (evaluated on the CA-125-available subset)
    result/ca125_values.csv      — per-sample CA-125 value + true label
    """
    from methylation.sample_grouping import load_ca125_map

    mode_cfg    = (config or {}).get("plasma", {})
    excel_sheet = mode_cfg.get("excel_sheet", "plasma")

    # ca125_map is keyed by normalised case_id (may be numeric like "640849")
    ca125_map = load_ca125_map(category=excel_sheet)
    if not ca125_map:
        logger.warning("[ca125] No CA-125 values found in Excel — skipping CA-125 comparison. "
                       "Add a 'ca125' column to the sample info Excel to enable this step.")
        return None, None

    pred_path = os.path.join(result_dir, "classifier_predictions.csv")
    if not os.path.exists(pred_path):
        logger.warning("[ca125] classifier_predictions.csv not found — run --step classify first.")
        return None, None

    pred_df = pd.read_csv(pred_path, encoding="utf-8")
    clfs = [c for c in ["LR", "SVM", "RF", "XGB", "LGBM"] if c in pred_df["classifier"].values]
    if not clfs:
        logger.warning("[ca125] No classifiers found in classifier_predictions.csv")
        return None, None

    # Two-step join: sample_id (fq_R1 prefix) → case_id → ca125 value
    # Same join used by _build_sample_date_map for temporal validation.
    sid_to_case = _build_sid_to_case_map(config)
    logger.info(f"[ca125] sid→case map: {len(sid_to_case)} entries; "
                f"ca125 map: {len(ca125_map)} entries")

    # Reference sample order from the first classifier
    ref = pred_df[pred_df["classifier"] == clfs[0]].sort_values("sample_id").reset_index(drop=True)
    sample_ids = ref["sample_id"].tolist()
    y_true_all = (ref["true_label"] == "malignant").astype(int).values

    def _lookup_ca125(sid):
        case_id = sid_to_case.get(sid) or sid   # fall back to sid itself if no mapping
        return ca125_map.get(_norm_id(case_id), np.nan)

    ca125_vals = np.array([_lookup_ca125(sid) for sid in sample_ids])
    has_ca125  = ~np.isnan(ca125_vals)
    n_total    = len(sample_ids)
    n_ca125    = int(has_ca125.sum())

    if n_ca125 == 0:
        logger.warning(
            "[ca125] No CA-125 values matched any classifier sample ID.\n"
            f"  sample_ids (first 5): {sample_ids[:5]}\n"
            f"  sid_to_case (first 5): {dict(list(sid_to_case.items())[:5])}\n"
            f"  ca125_map keys (first 5): {list(ca125_map.keys())[:5]}\n"
            "  Check that 'case_id' in the Excel links fq_R1 prefixes to the ca125 column."
        )
        return None, None

    logger.info(f"[ca125] {n_ca125}/{n_total} samples have CA-125 data")
    sub_idx   = np.where(has_ca125)[0]
    y_sub     = y_true_all[sub_idx]
    ca125_sub = ca125_vals[sub_idx]

    if len(np.unique(y_sub)) < 2:
        logger.warning("[ca125] CA-125 subset contains only one class — cannot compute AUC")
        return None, None

    rows = []

    # CA-125 AUC with bootstrap CI
    ca125_auc  = roc_auc_score(y_sub, ca125_sub)
    ca125_pred = (ca125_sub >= float(np.median(ca125_sub))).astype(int)
    ci_ca125   = _bootstrap_ci(y_sub, ca125_sub, ca125_pred)
    auc_lo, auc_hi = ci_ca125["auc"]
    rows.append({
        "method":          "CA-125",
        "n_total":         n_total,
        "n_ca125_subset":  n_ca125,
        "auc":             round(ca125_auc, 4),
        "auc_ci95_lo":     round(auc_lo, 4),
        "auc_ci95_hi":     round(auc_hi, 4),
    })

    # Classifier AUCs re-computed on the same subset for fair comparison
    sub_sids = {sample_ids[i] for i in sub_idx}
    for clf in clfs:
        clf_df = (pred_df[pred_df["classifier"] == clf]
                  .sort_values("sample_id")
                  .reset_index(drop=True))
        clf_df = clf_df[clf_df["sample_id"].isin(sub_sids)].set_index("sample_id")
        sid_order = [sample_ids[i] for i in sub_idx]
        clf_df = clf_df.reindex(sid_order)
        probs  = clf_df["prob_malignant"].astype(float).values
        preds  = (probs >= 0.5).astype(int)
        valid  = ~np.isnan(probs)
        if valid.sum() < 5 or len(np.unique(y_sub[valid])) < 2:
            continue
        clf_auc = roc_auc_score(y_sub[valid], probs[valid])
        ci_clf  = _bootstrap_ci(y_sub[valid], probs[valid], preds[valid])
        a_lo, a_hi = ci_clf["auc"]
        rows.append({
            "method":          clf,
            "n_total":         n_total,
            "n_ca125_subset":  int(valid.sum()),
            "auc":             round(clf_auc, 4),
            "auc_ci95_lo":     round(a_lo, 4),
            "auc_ci95_hi":     round(a_hi, 4),
        })

    df_out = pd.DataFrame(rows)
    comp_path = os.path.join(result_dir, "ca125_comparison.csv")
    df_out.to_csv(comp_path, index=False, encoding="utf-8")
    logger.info(f"[ca125] Comparison saved: {comp_path}")

    df_vals = pd.DataFrame({
        "sample_id":  [sample_ids[i] for i in sub_idx],
        "true_label": ["malignant" if y == 1 else "benign" for y in y_sub],
        "ca125":      ca125_sub,
    })
    vals_path = os.path.join(result_dir, "ca125_values.csv")
    df_vals.to_csv(vals_path, index=False, encoding="utf-8")
    logger.info(f"[ca125] Per-sample values saved: {vals_path}")

    logger.info("\n  ┌────────────────────────────────────────────────────────┐")
    logger.info(f"  │  CA-125 vs Classifier AUC  (n_ca125={n_ca125}, 95% CI)         │")
    logger.info("  ├──────────┬─────────────────────────────────────────────┤")
    for _, r in df_out.iterrows():
        logger.info(
            f"  │  {r['method']:<8} │  AUC {r['auc']:.3f}"
            f"  ({r['auc_ci95_lo']:.3f}–{r['auc_ci95_hi']:.3f})                    │"
        )
    logger.info("  └──────────┴─────────────────────────────────────────────┘")

    return df_out, df_vals
