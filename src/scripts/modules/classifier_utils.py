#!/usr/bin/env python3
"""
classifier_utils.py - Diagnostic classifier module (Step 15)

Evaluates five classifiers with leave-one-out cross-validation (LOO-CV):
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
import traceback
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

from sklearn.model_selection import LeaveOneOut

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


# ──────────────────────────────────────────────────────
# Step 15 entry point
# ──────────────────────────────────────────────────────

def run_classifiers(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Step 15 entry point: Nested LOO-CV evaluation of LR / SVM / RF / XGB / LGBM classifiers.

    Feature selection (Wilcoxon + LASSO + RF) is performed independently on each LOO fold
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

    logger.info(f"Starting classifier evaluation (Nested LOO-CV)  "
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

    # 3. Nested LOO-CV
    loo = LeaveOneOut()
    clf_names = list(classifiers.keys())

    # per-classifier prediction storage
    clf_y_true = {n: [] for n in clf_names}
    clf_y_prob = {n: [] for n in clf_names}
    clf_y_pred = {n: [] for n in clf_names}

    fold_n_features    = []  # features selected per fold
    fold_selected_names = []  # selected feature names per fold (for frequency stats)

    n_folds = X.shape[0]
    log_interval = max(1, n_folds // 10)

    for fold_idx, (train_idx, test_idx) in enumerate(loo.split(X)):
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
                prob = clf.predict_proba(X_te_sc)[0, 1]
                pred = clf.predict(X_te_sc)[0]
                clf_y_true[clf_name].append(int(y_test[0]))
                clf_y_prob[clf_name].append(float(prob))
                clf_y_pred[clf_name].append(int(pred))

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


def _build_sample_date_map(config):
    """
    Build {sample_id: collection_date} keyed by fq_R1 prefix (same as feature matrix).

    Two-step join:
      1. plasma_exp.xlsx  : case_id → collection_date
      2. sample_info_excel: fq_R1 prefix (sample_id) → case_id
    Both tables store integer case_ids as Excel floats ("640849.0"); normalise to "640849".

    Returns an empty dict on any failure (caller should check and skip).
    """

    def _norm_id(raw):
        """'640849.0' → '640849'; non-numeric strings kept as-is."""
        s = str(raw).strip()
        try:
            return str(int(float(s)))
        except (ValueError, OverflowError):
            return s

    # Step 1: plasma_exp.xlsx → {case_id: timestamp}
    exp_path = None
    for candidate in ["data/plasma_exp.xlsx", "plasma_exp.xlsx"]:
        if os.path.exists(candidate):
            exp_path = candidate
            break
    if exp_path is None:
        logger.warning("[temporal] plasma_exp.xlsx not found — skipping temporal validation")
        return {}

    try:
        df_exp = pd.read_excel(exp_path)
    except Exception as e:
        logger.warning(f"[temporal] Failed to read {exp_path}: {e}")
        return {}

    df_exp.columns = df_exp.columns.str.strip().str.lower()
    if "case_id" not in df_exp.columns or "datetime" not in df_exp.columns:
        logger.warning("[temporal] plasma_exp.xlsx missing required columns (case_id, datetime)")
        return {}

    case_date = {}
    for _, row in df_exp.iterrows():
        cid = _norm_id(row["case_id"])
        dt_raw = str(row["datetime"]).strip()
        ts = pd.to_datetime(dt_raw, format="%Y%m%d", errors="coerce")
        if pd.isna(ts):
            ts = pd.to_datetime(dt_raw, errors="coerce")
        if not pd.isna(ts):
            case_date[cid] = ts

    if not case_date:
        logger.warning("[temporal] No valid dates parsed from plasma_exp.xlsx")
        return {}

    # Step 2: sample_info_excel → {sample_id (fq_R1 prefix): case_id}
    excel_file = (config or {}).get("data", {}).get("sample_info_excel", "")
    if not excel_file or not os.path.exists(excel_file):
        logger.warning(f"[temporal] sample_info_excel not found: {excel_file!r}")
        return {}

    sid_to_case = {}
    try:
        xl = pd.ExcelFile(excel_file)
        for sheet in xl.sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet)
            # df.columns = df.columns.str.strip()
            if "fq_R1" not in df.columns or "case_id" not in df.columns:
                continue
            for _, row in df.iterrows():
                cid = _norm_id(row.get("case_id", ""))
                r1  = str(row.get("fq_R1", "")).strip()
                sid = r1.split(".")[0]
                if cid and cid != "nan" and sid and sid != "nan":
                    sid_to_case[sid] = cid
    except Exception as e:
        logger.warning(f"[temporal] Failed to read sample_info_excel: {e} {traceback.format_exc()}")
        return {}

    # Step 3: join → {sample_id: date}
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


def run_temporal_validation(workdir, result_dir, config, benign_samples, malignant_samples):
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

    # 3. Load full feature matrix (no panel filter — use freq-selected features below)
    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    # 4. Restrict to top-N features from feature_selection_freq.csv.
    #    If panel_cpg_blood_filtered.csv exists (written by Step 19 external_validate),
    #    additionally exclude blood-contaminated sites (blood_mean_beta >= 0.10).
    #    These were selected per LOO fold on training data in run_classifiers,
    #    so they carry no information from the temporal test set.
    freq_path = os.path.join(result_dir, "feature_selection_freq.csv")
    if not os.path.exists(freq_path):
        logger.warning(f"[temporal] feature_selection_freq.csv not found — run classify first: {freq_path}")
        return None, None
    freq_df = pd.read_csv(freq_path)
    top_features = set(freq_df["feature"].iloc[:TOP_N].tolist())

    # Apply blood-filter whitelist if available (Step 19 must have run first)
    blood_filtered_path = os.path.join(result_dir, "panel_cpg_blood_filtered.csv")
    if os.path.exists(blood_filtered_path):
        blood_panel_df = pd.read_csv(blood_filtered_path)
        blood_whitelist = set(blood_panel_df["feature"].tolist())
        n_before = len(top_features)
        top_features = top_features & blood_whitelist
        n_removed = n_before - len(top_features)
        if n_removed:
            logger.info(
                f"[temporal] Blood-filter applied: {n_removed} blood-contaminated site(s) "
                f"removed from top features → {len(top_features)} features retained"
            )
        else:
            logger.info("[temporal] Blood-filter: no top features were blood-contaminated")
    else:
        logger.debug(
            "[temporal] panel_cpg_blood_filtered.csv not found — "
            "run --step external_validate to enable blood-site exclusion"
        )

    top_idx = [i for i, fn in enumerate(feature_names) if fn in top_features]
    if not top_idx:
        logger.warning("[temporal] None of the top features found in feature matrix — skipping")
        return None, None
    X = X[:, top_idx]
    feature_names = [feature_names[i] for i in top_idx]
    logger.info(f"[temporal] Using top {len(top_idx)} features from feature_selection_freq.csv")

    # 5. Apply exclude_samples (consistent with run_classifiers)
    ml_cfg = (config or {}).get("ml", {})
    exclude_set = set(ml_cfg.get("exclude_samples", []))
    if exclude_set:
        keep_mask = [sid not in exclude_set for sid in sample_ids]
        X          = X[[i for i, k in enumerate(keep_mask) if k]]
        y          = y[[i for i, k in enumerate(keep_mask) if k]]
        sample_ids = [sid for sid, k in zip(sample_ids, keep_mask) if k]

    # 6. Partition by collection date — no test-set info leaks into train
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

    X_train, y_train = X[train_idx], y[train_idx]
    X_test,  y_test  = X[test_idx],  y[test_idx]

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
    test_dates = [date_map.get(sample_ids[i]) for i in test_idx]

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
    val_sample_ids = [sample_ids[i] for i in test_idx]
    pred_rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(val_sample_ids):
            pred_rows.append({
                "sample_id":       sid,
                "true_label":      "malignant" if res["y_true"][i] == 1 else "benign",
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
