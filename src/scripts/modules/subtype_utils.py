#!/usr/bin/env python3
"""
subtype_utils.py - Subtype classifier module (Step 16)

Multi-class LOO-CV on malignant + borderline samples using the "subtype" column
from the sample Excel.  Classifiers: RF and SVM (OvR).

Also runs a dedicated binary classifier: high-grade serous carcinoma (HGSC) vs all other subtypes.

Outputs (saved to result_dir):
    subtype_predictions.csv          -- per-sample predicted subtype + confidence (multi-class)
    subtype_summary.csv              -- accuracy per classifier, data-quality warnings
    subtype_skip.json                -- present only when training is skipped (too few data)
    subtype_confusion_{RF|SVM}.png   -- multi-class confusion matrix
    subtype_roc.png                  -- multi-class OvR ROC curves
    hgsoc_predictions.csv            -- per-sample HGSOC vs Others binary predictions
    hgsoc_summary.csv                -- AUC / accuracy / sensitivity / specificity
    hgsoc_confusion_{RF|SVM}.png     -- binary confusion matrix
    hgsoc_roc.png                    -- binary ROC curve
"""

import os
import json
import numpy as np
import pandas as pd
from collections import Counter

from utils.logger import logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Canonical English subtype labels used by the classifier.
# These five classes correspond to the standard histological subtypes of
# epithelial ovarian carcinoma. Other malignant labels in the input Excel
# (e.g. "double primary cancer", "unknown type") and all benign/borderline
# entries are dropped from subtype classification by _normalize_subtype().
SUBTYPE_CLASSES = [
    "HGSC",
    "Endometrioid",
    "Clear cell",
    "LGSC",
    "Mucinous",
]

# Identity map kept for backward compatibility with existing call sites that
# previously translated Chinese → English. After normalization, labels are
# already in canonical English form, so `SUBTYPE_EN.get(c, c)` is a safe no-op.
SUBTYPE_EN = {k: k for k in SUBTYPE_CLASSES}

# Subtypes that are recognised but explicitly non-classifiable
# (logged once per run; samples are dropped from subtype classification).
_SUBTYPE_SKIP_LITERALS = {
    "unknown", "unknown type", "unknown typ",
    "double primary cancer", "double primary",
    "n/a", "na", "nan", "",
}

# Legacy Chinese exact-match aliases — preserved for back-compatibility with
# older sample-info Excels that haven't been migrated to English yet.
_SUBTYPE_LEGACY_CN = {
    "高级别浆液性癌": "HGSC",
    "子宫内膜样癌":   "Endometrioid",
    "透明细胞癌":     "Clear cell",
    "低级别浆液性癌": "LGSC",
    "黏液性":         "Mucinous",
    "粘液性":         "Mucinous",   # alternate variant of 黏
}


def _normalize_subtype(raw):
    """
    Map a raw subtype string to one of SUBTYPE_CLASSES.

    Recognised inputs (case-insensitive; hyphen/underscore-tolerant):

      ``"high-grade serous carcinoma"`` / ``"high grade serous carcinoma"`` /
      ``"hgsc"``                                                  → ``"HGSC"``
      ``"low-grade serous carcinoma"`` /  ``"low grade serous carcinoma"`` /
      ``"lgsc"``                                                  → ``"LGSC"``
      ``"clear cell carcinoma"`` / ``"clear cell"``               → ``"Clear cell"``
      ``"endometrioid carcinoma"`` / ``"endometrioid"``           → ``"Endometrioid"``
      ``"mucinous carcinoma"`` / ``"mucinous"``                   → ``"Mucinous"``

    Returns ``None`` (sample is skipped) for:
      - any value in ``_SUBTYPE_SKIP_LITERALS`` (unknown / double-primary / NA)
      - plain ``"serous"`` without a high/low-grade qualifier — clinically
        ambiguous between HGSC, LGSC, and serous borderline tumour
      - benign-condition labels that slip through the category filter
        (e.g. ``"corpus luteum cyst"``, ``"fibroma"``, ``"...cystadenoma"``)
      - any other unrecognised string (caller logs an aggregated warning)

    Legacy Chinese labels remain supported.
    """
    if raw is None:
        return None
    s_raw = str(raw).strip()
    if not s_raw:
        return None

    # Legacy Chinese exact match first
    if s_raw in _SUBTYPE_LEGACY_CN:
        return _SUBTYPE_LEGACY_CN[s_raw]

    # Normalise English: case-fold, treat hyphens/underscores as spaces,
    # collapse whitespace.
    s = s_raw.lower().replace("-", " ").replace("_", " ")
    s = " ".join(s.split())

    if s in _SUBTYPE_SKIP_LITERALS:
        return None

    # Carcinoma matches — order matters: more-specific patterns first.
    if "high grade serous" in s:
        return "HGSC"
    if "low grade serous" in s:
        return "LGSC"
    if s == "hgsc":
        return "HGSC"
    if s == "lgsc":
        return "LGSC"
    if "clear cell" in s:
        return "Clear cell"
    if "endometrioid" in s:
        return "Endometrioid"
    if "mucinous carcinoma" in s or s == "mucinous":
        # Plain "mucinous" → Mucinous carcinoma (single canonical bucket).
        # If category-aware borderline handling is needed later, branch here.
        return "Mucinous"

    # Plain "serous" without grade qualifier — clinically ambiguous, drop.
    if s == "serous":
        return None

    # Anything else (cystadenomas, fibromas, simple cysts, hemorrhages, etc.)
    # is silently unrecognised; caller aggregates these into one log entry.
    return None

# Default viability thresholds — overridden by config.toml [ml].
# These mirror the keys documented in config.toml; if the user has not set them,
# we fall back to these conservative defaults rather than failing.
_DEFAULT_MIN_SAMPLES_TOTAL     = 10
_DEFAULT_MIN_SAMPLES_PER_CLASS = 2
_DEFAULT_WARN_SAMPLES_TOTAL    = 20


def _read_subtype_thresholds(config):
    """Read viability thresholds from config.toml [ml] with safe fallbacks."""
    ml = (config or {}).get("ml", {}) if config else {}
    return (
        int(ml.get("subtype_min_samples_total",     _DEFAULT_MIN_SAMPLES_TOTAL)),
        int(ml.get("subtype_min_samples_per_class", _DEFAULT_MIN_SAMPLES_PER_CLASS)),
        int(ml.get("subtype_warn_samples_total",    _DEFAULT_WARN_SAMPLES_TOTAL)),
    )


# ---------------------------------------------------------------------------
# Viability check
# ---------------------------------------------------------------------------

def _check_viability(labels, config=None):
    """
    Decide whether the dataset is large enough for meaningful LOO-CV.

    Thresholds are read from config.toml [ml]:
      - ``subtype_min_samples_total``     (default 10) — hard floor; below this
        Step 16 is **skipped** and a `subtype_skip.json` marker is written.
      - ``subtype_min_samples_per_class`` (default 2)  — classes below this trigger
        a per-class reliability warning.
      - ``subtype_warn_samples_total``    (default 20) — below this (but ≥ min),
        training proceeds with an "indicative-only" warning.

    Returns:
        viable     (bool)       -- False means skip training entirely
        reason     (str)        -- plain-English explanation when not viable
        warnings   (list[str])  -- non-fatal reliability concerns
        thresholds (dict)       -- the actual threshold values used (for logging
                                   and inclusion in `subtype_skip.json`)
    """
    min_total, min_per_class, warn_total = _read_subtype_thresholds(config)
    thresholds = {
        "subtype_min_samples_total":     min_total,
        "subtype_min_samples_per_class": min_per_class,
        "subtype_warn_samples_total":    warn_total,
    }

    counts    = Counter(labels)
    n_total   = len(labels)
    n_classes = len(counts)

    if n_total < min_total:
        return False, (
            f"Only {n_total} labelled sample(s) — minimum {min_total} required for "
            f"reliable multi-class LOO-CV (configurable via "
            f"`ml.subtype_min_samples_total` in config.toml). "
            f"Class distribution: {dict(counts)}. "
            "Collect more malignant samples before running subtype classification, "
            "or lower the threshold if you accept reduced statistical reliability."
        ), [], thresholds

    if n_classes < 2:
        return False, (
            f"Only 1 subtype class present ({list(counts.keys())[0]}, n={n_total}). "
            "At least 2 distinct subtypes are required."
        ), [], thresholds

    warnings = []
    if n_total < warn_total:
        warnings.append(
            f"Small dataset ({n_total} < {warn_total} samples total): LOO-CV "
            "results should be treated as indicative only "
            "(threshold: `ml.subtype_warn_samples_total`)."
        )

    singletons = [c for c, n in counts.items() if n < min_per_class]
    if singletons:
        en = [SUBTYPE_EN.get(c, c) for c in singletons]
        warnings.append(
            f"Class(es) with < {min_per_class} sample(s) — will never appear in "
            f"training for their leave-out fold: {en}. Predictions for these "
            f"classes are unreliable (threshold: `ml.subtype_min_samples_per_class`)."
        )

    rare = [c for c, n in counts.items() if min_per_class <= n < 3]
    if rare:
        en = [SUBTYPE_EN.get(c, c) for c in rare]
        warnings.append(
            f"Class(es) with <3 samples: {en}. "
            "AUC and per-class sensitivity will be near-random for these groups."
        )

    return True, "", warnings, thresholds


def _save_skip_summary(result_dir, reason, warnings=None, thresholds=None,
                       n_total=None, class_distribution=None):
    """Write a skip marker so the report can display the reason and the
    threshold values that triggered the skip."""
    path = os.path.join(result_dir, "subtype_skip.json")
    payload = {
        "skipped": True,
        "reason":   reason,
        "warnings": warnings or [],
    }
    if thresholds is not None:
        payload["thresholds"] = thresholds
    if n_total is not None:
        payload["n_total"] = n_total
    if class_distribution is not None:
        payload["class_distribution"] = class_distribution
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    logger.info(f"  Subtype skip reason saved: {path}")


# ---------------------------------------------------------------------------
# Label reading
# ---------------------------------------------------------------------------

def _build_subtype_map(config):
    """
    Read all sheets from sample_info_excel and return {sample_id: subtype_str}.
    Only retains rows where category is 'malignant' or 'borderline' and subtype is non-empty.
    """
    excel_file = config.get("data", {}).get("sample_info_excel", "")
    if not excel_file or not os.path.exists(excel_file):
        logger.error(f"sample_info_excel not found: {excel_file}")
        return {}

    xl = pd.ExcelFile(excel_file)
    subtype_map = {}
    for sheet in xl.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet)
        logger.info(f"  sheet '{sheet}' columns: {list(df.columns)}")
        if "fq_R1" not in df.columns:
            logger.info("  → no fq_R1 column, skip")
            continue

        if "category" in df.columns:
            df["category"] = df["category"].astype(str).str.strip()
            logger.info(f"  → category distribution: {df['category'].value_counts().to_dict()}")
        if "subtype" not in df.columns:
            logger.warning("  → no 'subtype' column, skip")
            continue
        df["subtype"] = df["subtype"].astype(str).str.strip()
        logger.info(f"  → subtype distribution: {df['subtype'].value_counts().to_dict()}")

        # Aggregate unrecognised subtype labels into a single warning per sheet
        # so the log isn't flooded with per-sample messages.
        unrecognised = Counter()
        skipped_explicit = Counter()

        for _, row in df.iterrows():
            category = str(row.get("category", "")).strip() if "category" in df.columns else "malignant"
            if category not in ("malignant", "borderline"):
                continue
            raw_subtype = str(row.get("subtype", "")).strip()
            if not raw_subtype or raw_subtype.lower() == "nan":
                continue

            canonical = _normalize_subtype(raw_subtype)
            if canonical is None:
                # Distinguish "explicitly skip-listed" (known non-classifiable
                # labels) from "truly unrecognised" (typo / new pathology).
                s_norm = raw_subtype.lower().replace("-", " ").replace("_", " ")
                s_norm = " ".join(s_norm.split())
                if s_norm in _SUBTYPE_SKIP_LITERALS or s_norm == "serous":
                    skipped_explicit[raw_subtype] += 1
                else:
                    unrecognised[raw_subtype] += 1
                continue

            r1  = str(row.get("fq_R1", "")).strip()
            # Canonical sample-id split: matches feature_matrix.npz keys
            # (e.g. 'B10013_R1.fq.gz' → 'B10013'), same scheme used by
            # classifier_utils._build_sid_to_case_map.
            pid = r1.split(".")[0].split("_")[0]
            if pid and pid != "nan":
                subtype_map[pid] = canonical

        if skipped_explicit:
            logger.info(
                f"  → non-classifiable subtypes dropped (expected): "
                f"{dict(skipped_explicit)}"
            )
        if unrecognised:
            logger.warning(
                f"  → unrecognised subtype labels (dropped): "
                f"{dict(unrecognised)} — extend _normalize_subtype() if these "
                f"are legitimate new categories."
            )

    dist = Counter(subtype_map.values())
    logger.info(f"Subtype labels loaded (canonical): {len(subtype_map)} samples")
    for k, v in sorted(dist.items(), key=lambda x: -x[1]):
        logger.info(f"  {k}: {v}")
    return subtype_map


# ---------------------------------------------------------------------------
# Feature matrix helpers (reuse Step 14 cache)
# ---------------------------------------------------------------------------

def _load_feature_matrix(result_dir):
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if not os.path.exists(fm_path):
        raise FileNotFoundError(
            f"feature_matrix.npz not found — run feature_analysis first: {fm_path}"
        )
    data = np.load(fm_path, allow_pickle=True)
    X            = data["X"]
    sample_ids   = list(data["sample_ids"])
    feature_names = list(data["feature_names"])
    logger.info(f"Feature matrix loaded: {X.shape[0]} samples × {X.shape[1]} features")
    return X, sample_ids, feature_names


def _apply_panel_filter(X, feature_names, result_dir):
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if os.path.exists(panel_path):
        panel_set = set(pd.read_csv(panel_path)["feature"])
        indices = [i for i, fn in enumerate(feature_names) if fn in panel_set]
        logger.info(f"Panel filter applied: {len(indices)} features retained")
        return X[:, indices], [feature_names[i] for i in indices]
    return X, feature_names


# Fold-wise imputation helpers (mirrors classifier_utils to avoid data leakage)
def _impute_fit(X):
    col_medians = np.nanmedian(X, axis=0)
    col_medians = np.where(np.isnan(col_medians), 0.0, col_medians)
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp, col_medians


def _impute_transform(X, col_medians):
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def _save_confusion_plot(cm, en_labels, clf_name, result_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"]})

    n = len(en_labels)
    fig, ax = plt.subplots(figsize=(max(6, n), max(5, n)))
    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    ax.set_title(f"{clf_name} — Subtype Confusion Matrix (LOO-CV)")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(en_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(en_labels, fontsize=9)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    for i in range(n):
        for j in range(n):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black",
                    fontsize=11)
    fig.tight_layout()
    path = os.path.join(result_dir, f"subtype_confusion_{clf_name}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  Confusion matrix saved: {path}")


def _save_roc_plot(results, en_labels, result_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc
    from sklearn.preprocessing import label_binarize
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"]})

    n_classes = len(en_labels)
    n_clf = len(results)
    fig, axes = plt.subplots(1, n_clf, figsize=(6 * n_clf, 5))
    if n_clf == 1:
        axes = [axes]
    colors = plt.cm.tab10.colors

    for ax, (clf_name, res) in zip(axes, results.items()):
        y_true_bin = label_binarize(res["y_true"], classes=list(range(n_classes)))
        if y_true_bin.ndim == 1:           # binary edge case
            y_true_bin = np.column_stack([1 - y_true_bin, y_true_bin])
        y_prob = res["y_prob"]
        for i, lbl in enumerate(en_labels):
            if i >= y_true_bin.shape[1] or i >= y_prob.shape[1]:
                continue
            if y_true_bin[:, i].sum() == 0:
                continue
            fpr, tpr, _ = roc_curve(y_true_bin[:, i], y_prob[:, i])
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
                    label=f"{lbl} (AUC={roc_auc:.2f})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_title(f"{clf_name} OvR ROC")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)

    fig.tight_layout()
    path = os.path.join(result_dir, "subtype_roc.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  ROC curves saved: {path}")



# ---------------------------------------------------------------------------
# One-vs-Rest binary classifiers  (each subtype vs all others)
# ---------------------------------------------------------------------------

def _run_ovr_binary(X_raw, labels, sample_ids, observed_classes, result_dir, config):
    """
    For every class in observed_classes run a binary LOO-CV (class vs all others).

    Skips a class if it has fewer than 2 positive or 2 negative samples.

    Outputs saved to result_dir:
        ovr_summary.csv          — AUC / accuracy / sensitivity / specificity per (subtype, clf)
        ovr_predictions.csv      — per-sample binary predictions for every (subtype, clf)
        ovr_roc.png              — combined ROC: one subplot per classifier, one curve per subtype
        ovr_confusion_{clf}_{en}.png — confusion matrix per (clf, subtype)

    Returns list of summary-row dicts (same content as ovr_summary.csv).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import (accuracy_score, confusion_matrix,
                                 roc_curve, auc as sklearn_auc)

    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["DejaVu Sans"]})

    ml_cfg = (config or {}).get("ml", {})
    # OvR per-class threshold — configurable via config.toml [ml].
    # Subtypes failing this check are skipped with a clear log; if NO subtype
    # passes, the whole OvR step writes ovr_skip.json so the report can render
    # a clean "Skipped" panel.
    ovr_min = int(ml_cfg.get("ovr_min_samples_per_class", 2))
    logger.info(
        f"  OvR per-class threshold (config.toml [ml].ovr_min_samples_per_class): "
        f"{ovr_min} (need ≥{ovr_min} positives AND ≥{ovr_min} negatives per subtype)"
    )

    clf_factories = {
        "RF":  lambda: RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced", random_state=42, n_jobs=-1,
        ),
        "SVM": lambda: SVC(
            kernel="rbf", C=ml_cfg.get("svm_C", 1.0),
            probability=True, class_weight="balanced", random_state=42,
        ),
    }
    clf_names   = list(clf_factories.keys())
    loo         = LeaveOneOut()
    summary_all = []   # list of dicts → ovr_summary.csv
    pred_all    = []   # list of dicts → ovr_predictions.csv

    # Track per-class viability outcomes so the report can show "skipped because …"
    trained_classes = []         # list[str]  — subtypes that were trained
    skipped_classes = []         # list[dict] — {subtype, n_positive, n_negative, reason}

    # Per-classifier ROC data collected across all subtypes for the combined plot
    roc_data = {cn: {} for cn in clf_names}   # {clf_name: {subtype_en: (fpr, tpr, auc)}}

    # Exclude the "Others" merge bucket from OvR entirely (no training, no
    # confusion matrix, no ROC, no entry in any OvR output file or report
    # warning). It is a heterogeneous catch-all of below-threshold rare
    # subtypes — a one-vs-rest binary against it has no clean clinical meaning
    # and would clutter the OvR plots. "Others" still appears in the
    # multi-class confusion matrix above where its catch-all nature is
    # appropriate; it is treated as if it were never an OvR candidate.
    ovr_candidates = [c for c in observed_classes if c != "Others"]
    n_others_excluded = len(observed_classes) - len(ovr_candidates)
    if n_others_excluded:
        logger.info(
            f"  OvR: 'Others' merge bucket is excluded by design "
            f"(heterogeneous catch-all of rare subtypes — not a clinically "
            f"meaningful binary target)."
        )

    for cls in ovr_candidates:
        en = SUBTYPE_EN.get(cls, cls)
        y_bin  = np.array([1 if lbl == cls else 0 for lbl in labels])
        n_pos  = int(y_bin.sum())
        n_neg  = len(y_bin) - n_pos

        if n_pos < ovr_min or n_neg < ovr_min:
            reason = (f"need ≥{ovr_min} positive AND ≥{ovr_min} negative "
                      f"(got pos={n_pos}, neg={n_neg}); "
                      f"threshold: ml.ovr_min_samples_per_class")
            logger.warning(f"  OvR [{en}]: skipped — {reason}")
            skipped_classes.append({
                "subtype":    en,
                "n_positive": n_pos,
                "n_negative": n_neg,
                "reason":     reason,
            })
            continue

        trained_classes.append(en)
        logger.info(f"  OvR [{en}]: pos={n_pos}, neg={n_neg}")

        for clf_name, factory in clf_factories.items():
            y_true_list, y_prob_list, y_pred_list = [], [], []

            for train_idx, test_idx in loo.split(X_raw):
                X_tr_imp, col_meds = _impute_fit(X_raw[train_idx])
                X_te_imp           = _impute_transform(X_raw[test_idx], col_meds)
                scaler             = StandardScaler()
                X_tr_sc            = scaler.fit_transform(X_tr_imp)
                X_te_sc            = scaler.transform(X_te_imp)

                clf = factory()
                clf.fit(X_tr_sc, y_bin[train_idx])

                cls_list   = list(clf.classes_)
                prob_raw   = clf.predict_proba(X_te_sc)[0]
                pos_prob   = prob_raw[cls_list.index(1)] if 1 in cls_list else 0.0
                pred       = int(clf.predict(X_te_sc)[0])

                y_true_list.append(int(y_bin[test_idx[0]]))
                y_prob_list.append(pos_prob)
                y_pred_list.append(pred)

            y_true   = np.array(y_true_list)
            y_prob   = np.array(y_prob_list)
            y_pred   = np.array(y_pred_list)

            acc = accuracy_score(y_true, y_pred)
            cm  = confusion_matrix(y_true, y_pred, labels=[0, 1])
            tn, fp, fn, tp = cm.ravel() if cm.shape == (2, 2) else (0, 0, 0, 0)
            sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0
            fpr, tpr, _ = roc_curve(y_true, y_prob)
            auc_val     = sklearn_auc(fpr, tpr)

            logger.info(
                f"    {clf_name}: AUC={auc_val:.3f}  acc={acc:.3f}  "
                f"sens={sensitivity:.3f}  spec={specificity:.3f}"
            )
            roc_data[clf_name][en] = (fpr, tpr, auc_val)

            # Confusion matrix plot
            _save_confusion_plot(cm, ["Other", en], f"OvR_{clf_name}_{en}", result_dir)

            summary_all.append({
                "subtype":      en,
                "classifier":   clf_name,
                "auc":          round(auc_val, 4),
                "accuracy":     round(acc, 4),
                "sensitivity":  round(sensitivity, 4),
                "specificity":  round(specificity, 4),
                "n_positive":   n_pos,
                "n_negative":   n_neg,
            })
            for i, sid in enumerate(sample_ids):
                pred_all.append({
                    "subtype":      en,
                    "classifier":   clf_name,
                    "sample_id":    sid,
                    "true_subtype": SUBTYPE_EN.get(labels[i], labels[i]),
                    "true_label":   "Pos" if y_true[i] == 1 else "Neg",
                    "predicted":    "Pos" if y_pred[i] == 1 else "Neg",
                    "prob_pos":     round(float(y_prob[i]), 4),
                    "correct":      bool(y_true[i] == y_pred[i]),
                })

    # Whole-step skip: NO subtype had enough samples → write ovr_skip.json so
    # the report can render a clean Skipped panel parallel to subtype_skip.json.
    if not summary_all:
        reason = (
            f"No subtype met the OvR per-class minimum "
            f"(ml.ovr_min_samples_per_class = {ovr_min}). "
            f"Per-class details: {skipped_classes}"
        )
        logger.error(f"  OvR binary: SKIP — {reason}")
        skip_path = os.path.join(result_dir, "ovr_skip.json")
        with open(skip_path, "w", encoding="utf-8") as f:
            json.dump({
                "step_skipped":    True,
                "reason":          reason,
                "thresholds":      {"ovr_min_samples_per_class": ovr_min},
                "skipped_classes": skipped_classes,
                "trained_classes": [],
                "n_observed_classes": len(observed_classes),
            }, f, ensure_ascii=False, indent=2)
        logger.info(f"  OvR skip marker saved: {skip_path}")
        return []

    # At least one subtype trained → clear any stale whole-step skip marker
    # left over from a previous all-failed run.
    stale_skip = os.path.join(result_dir, "ovr_skip.json")
    if os.path.exists(stale_skip):
        try:
            os.remove(stale_skip)
            logger.info(f"  Removed stale OvR skip marker: {stale_skip}")
        except OSError as exc:
            logger.warning(f"  Could not remove stale {stale_skip}: {exc}")

    # Partial skip: some subtypes trained, some skipped. Persist the per-class
    # skip list so the report can show which subtypes were excluded and why.
    partial_path = os.path.join(result_dir, "ovr_skipped_classes.json")
    if skipped_classes:
        with open(partial_path, "w", encoding="utf-8") as f:
            json.dump({
                "step_skipped":    False,
                "thresholds":      {"ovr_min_samples_per_class": ovr_min},
                "skipped_classes": skipped_classes,
                "trained_classes": trained_classes,
                "n_observed_classes": len(observed_classes),
            }, f, ensure_ascii=False, indent=2)
        logger.info(
            f"  OvR partial skip: trained={trained_classes}, "
            f"skipped={[c['subtype'] for c in skipped_classes]} "
            f"→ {partial_path}"
        )
    else:
        # All classes passed → clear any stale partial-skip marker too.
        if os.path.exists(partial_path):
            try:
                os.remove(partial_path)
                logger.info(f"  Removed stale OvR partial-skip marker: {partial_path}")
            except OSError as exc:
                logger.warning(f"  Could not remove stale {partial_path}: {exc}")

    # Save CSVs
    pd.DataFrame(summary_all).to_csv(
        os.path.join(result_dir, "ovr_summary.csv"), index=False
    )
    pd.DataFrame(pred_all).to_csv(
        os.path.join(result_dir, "ovr_predictions.csv"), index=False
    )

    # Combined ROC plot: one subplot per classifier, one curve per subtype
    colors = plt.cm.tab10.colors
    n_clf  = len(clf_names)
    fig, axes = plt.subplots(1, n_clf, figsize=(6 * n_clf, 5))
    if n_clf == 1:
        axes = [axes]
    for ax, clf_name in zip(axes, clf_names):
        for i, (en, (fpr, tpr, auc_val)) in enumerate(roc_data[clf_name].items()):
            ax.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
                    label=f"{en} (AUC={auc_val:.2f})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_title(f"{clf_name} — One-vs-Rest ROC")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    roc_path = os.path.join(result_dir, "ovr_roc.png")
    fig.savefig(roc_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  OvR ROC saved: {roc_path}")
    logger.info(f"  OvR binary outputs saved to {result_dir}")
    return summary_all


# ---------------------------------------------------------------------------
# Main entry: Step 16
# ---------------------------------------------------------------------------

def run_subtype_classifier(workdir, result_dir, config):
    """
    Step 16: Multi-class LOO-CV on malignant + borderline samples.
    Also runs HGSOC vs Others binary LOO-CV.

    Improvements over the original:
    - Viability check with clear skip message when data is insufficient
    - LabelEncoder fitted only on observed classes (no empty confusion matrix rows)
    - Fold-wise imputation to prevent data leakage
    - Confidence score (max predicted probability) in predictions CSV
    - Data-quality warnings propagated to summary CSV and report

    Returns:
        results:    dict {clf_name: {y_true, y_prob, y_pred, accuracy, confusion_matrix}}
        summary_df: pd.DataFrame
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

    # 1. Build label map
    subtype_map = _build_subtype_map(config)
    if not subtype_map:
        logger.error("Cannot read subtype labels — skipping subtype classification")
        return {}, pd.DataFrame()

    # 2. Load feature matrix
    try:
        X_all, all_sample_ids, feature_names = _load_feature_matrix(result_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        return {}, pd.DataFrame()

    X_all, feature_names = _apply_panel_filter(X_all, feature_names, result_dir)

    # 3. Filter to samples with known subtype labels
    indices, labels, sample_ids = [], [], []
    for i, sid in enumerate(all_sample_ids):
        subtype = subtype_map.get(sid)
        if subtype is None:
            continue
        if subtype not in SUBTYPE_CLASSES:
            logger.warning(f"  Sample {sid}: unknown subtype '{subtype}' — skipped")
            continue
        indices.append(i)
        labels.append(subtype)
        sample_ids.append(sid)

    n_benign_skipped = len(all_sample_ids) - len(indices) - sum(
        1 for sid in all_sample_ids if subtype_map.get(sid) not in SUBTYPE_CLASSES + [None]
    )
    logger.info(f"  {n_benign_skipped} samples have no subtype label (benign) — skipped")
    logger.debug(f"  feature_matrix sample_ids[:5]: {all_sample_ids[:5]}")
    logger.debug(f"  subtype_map keys[:5]: {list(subtype_map.keys())[:5]}")

    # 4a. Merge rare classes into "Others" (multi-class + OvR consistency)
    #
    # Classes with fewer than `subtype_min_samples_per_class` samples are merged
    # into a single "Others" bucket rather than (a) treated as their own
    # statistically-unreliable class or (b) silently inflating warnings. Applied
    # BEFORE the viability check so every downstream stage — multi-class LOO,
    # confusion matrix, OvR binary, ROC, predictions CSV — operates on the
    # merged labels consistently.
    #
    # If the merged "Others" bucket itself is still below the threshold (i.e.
    # there are very few rare-class samples in total), those samples are dropped
    # entirely rather than creating an unreliable single-sample-collection class.
    OTHERS_LABEL    = "Others"
    _ml_cfg         = (config or {}).get("ml", {})
    min_per_class   = int(_ml_cfg.get("subtype_min_samples_per_class", 2))
    pre_merge_counts = Counter(labels)
    rare_classes    = sorted(c for c, n in pre_merge_counts.items() if n < min_per_class)

    if rare_classes:
        n_rare_total = sum(pre_merge_counts[c] for c in rare_classes)
        rare_summary = ", ".join(
            f"{SUBTYPE_EN.get(c, c)}={pre_merge_counts[c]}" for c in rare_classes
        )
        if n_rare_total < min_per_class:
            # Even after merging, "Others" would be too small → drop entirely.
            n_before    = len(labels)
            keep_mask   = [c not in set(rare_classes) for c in labels]
            indices     = [i for i, k in zip(indices,    keep_mask) if k]
            labels      = [l for l, k in zip(labels,     keep_mask) if k]
            sample_ids  = [s for s, k in zip(sample_ids, keep_mask) if k]
            logger.warning(
                f"  [subtype] Rare classes ({rare_summary}) total only {n_rare_total} "
                f"sample(s) — below `ml.subtype_min_samples_per_class = {min_per_class}` "
                f"even after merging into '{OTHERS_LABEL}'. Dropping "
                f"{n_before - len(labels)} sample(s) entirely to avoid an "
                f"unreliable class."
            )
        else:
            # Relabel rare-class samples to 'Others' (in-place; counts re-derived below)
            labels = [OTHERS_LABEL if c in set(rare_classes) else c for c in labels]
            logger.info(
                f"  [subtype] Merged {len(rare_classes)} rare class(es) into "
                f"'{OTHERS_LABEL}' (per `ml.subtype_min_samples_per_class = {min_per_class}`): "
                f"{rare_summary} → '{OTHERS_LABEL}'={n_rare_total}"
            )

    # 4. Viability check — fail fast with a clear, config-driven message
    viable, reason, warnings, thresholds = _check_viability(labels, config)
    logger.info(
        f"  [subtype] Viability thresholds (config.toml [ml]): "
        f"min_total={thresholds['subtype_min_samples_total']}, "
        f"min_per_class={thresholds['subtype_min_samples_per_class']}, "
        f"warn_total={thresholds['subtype_warn_samples_total']}"
    )
    for w in warnings:
        logger.warning(f"  [subtype] {w}")
    if not viable:
        logger.error(f"  [subtype] SKIP — {reason}")
        _save_skip_summary(
            result_dir, reason, warnings,
            thresholds=thresholds,
            n_total=len(labels),
            class_distribution=dict(Counter(labels)),
        )
        return {}, pd.DataFrame()

    # Training will proceed → clear any stale skip markers from previous runs
    # so the report doesn't show a "Skipped" panel for a cohort that now passes.
    for stale in ("subtype_skip.json", "ovr_skip.json", "ovr_skipped_classes.json"):
        p = os.path.join(result_dir, stale)
        if os.path.exists(p):
            try:
                os.remove(p)
                logger.info(f"  [subtype] Removed stale skip marker: {p}")
            except OSError as exc:
                logger.warning(f"  [subtype] Could not remove stale {p}: {exc}")

    X_raw = X_all[indices]

    # Encode only observed classes (cleaner confusion matrix)
    observed_classes = sorted(set(labels), key=lambda c: SUBTYPE_CLASSES.index(c)
                              if c in SUBTYPE_CLASSES else 99)
    # CRITICAL: do NOT use sklearn's LabelEncoder here — it internally calls
    # np.unique() which re-sorts the classes alphabetically, ignoring the
    # `observed_classes` ordering above. That produced a silent label/integer
    # mismatch (e.g. predictions correctly biased to HGSC were reported as
    # "Clear cell" because of alphabetical re-sort). We do the encoding by hand
    # so the integer ↔ class-name mapping is guaranteed to follow the intended
    # SUBTYPE_CLASSES display order throughout the module.
    class_to_idx  = {c: i for i, c in enumerate(observed_classes)}
    idx_to_class  = {i: c for c, i in class_to_idx.items()}
    y             = np.array([class_to_idx[l] for l in labels], dtype=int)
    n_classes     = len(observed_classes)
    en_labels     = [SUBTYPE_EN.get(c, c) for c in observed_classes]

    counts = Counter(labels)
    logger.info(f"Subtype classification: {len(sample_ids)} samples, {n_classes} class(es)")
    for cls in observed_classes:
        logger.info(f"  {SUBTYPE_EN.get(cls, cls):>14s}: {counts[cls]}")

    # 5. Classifiers
    ml_cfg = (config or {}).get("ml", {})
    classifiers = {
        "RF": RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "SVM": SVC(
            kernel="rbf",
            C=ml_cfg.get("svm_C", 1.0),
            probability=True,
            class_weight="balanced",
            decision_function_shape="ovr",
            random_state=42,
        ),
    }

    # 6. LOO-CV with fold-wise imputation (no leakage)
    loo     = LeaveOneOut()
    results = {}

    for clf_name, clf in classifiers.items():
        logger.info(f"  LOO-CV: {clf_name} ...")
        y_true_list, y_prob_list, y_pred_list = [], [], []

        for train_idx, test_idx in loo.split(X_raw):
            X_tr_raw, X_te_raw = X_raw[train_idx], X_raw[test_idx]
            y_train,  y_test   = y[train_idx],     y[test_idx]

            # Fold-wise imputation — fit on train only
            X_tr_imp, col_meds = _impute_fit(X_tr_raw)
            X_te_imp           = _impute_transform(X_te_raw, col_meds)

            scaler  = StandardScaler()
            X_tr_sc = scaler.fit_transform(X_tr_imp)
            X_te_sc = scaler.transform(X_te_imp)

            clf.fit(X_tr_sc, y_train)

            # Align probabilities to full n_classes vector
            # clf.classes_ may be a subset when some class absent from this fold
            prob_raw  = clf.predict_proba(X_te_sc)[0]
            prob_full = np.zeros(n_classes)
            for j, cls_idx in enumerate(clf.classes_):
                prob_full[int(cls_idx)] = prob_raw[j]

            pred = int(clf.predict(X_te_sc)[0])

            y_true_list.append(int(y_test[0]))
            y_prob_list.append(prob_full)
            y_pred_list.append(pred)

        y_true = np.array(y_true_list)
        y_prob = np.array(y_prob_list)   # (n_samples, n_classes)
        y_pred = np.array(y_pred_list)

        acc = accuracy_score(y_true, y_pred)
        cm  = confusion_matrix(y_true, y_pred, labels=list(range(n_classes)))

        results[clf_name] = {
            "y_true":           y_true,
            "y_prob":           y_prob,
            "y_pred":           y_pred,
            "accuracy":         acc,
            "confusion_matrix": cm,
        }
        logger.info(f"    {clf_name}: accuracy={acc:.3f}")
        logger.info("\n" + classification_report(
            y_true, y_pred,
            labels=list(range(n_classes)),
            target_names=en_labels,
            zero_division=0,
        ))

    # 7. Plots
    for clf_name, res in results.items():
        _save_confusion_plot(res["confusion_matrix"], en_labels, clf_name, result_dir)
    _save_roc_plot(results, en_labels, result_dir)

    # 8. Per-sample predictions (confidence = max predicted probability)
    rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(sample_ids):
            row = {
                "sample_id":         sid,
                "true_subtype":      idx_to_class[int(res["y_true"][i])],
                "classifier":        clf_name,
                "predicted_subtype": idx_to_class[int(res["y_pred"][i])],
                "confidence":        round(float(res["y_prob"][i].max()), 4),
                "correct":           bool(res["y_true"][i] == res["y_pred"][i]),
            }
            for j, cls in enumerate(observed_classes):
                col = f"prob_{SUBTYPE_EN.get(cls, cls[:6])}"
                row[col] = round(float(res["y_prob"][i, j]), 4)
            rows.append(row)
    pred_df   = pd.DataFrame(rows)
    pred_path = os.path.join(result_dir, "subtype_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"  Per-sample predictions saved: {pred_path}")

    # 9. Summary (includes data-quality warnings)
    summary_rows = []
    for clf_name, res in results.items():
        summary_rows.append({
            "classifier":   clf_name,
            "accuracy":     round(res["accuracy"], 4),
            "n_samples":    len(sample_ids),
            "n_subtypes":   n_classes,
            "class_counts": str({SUBTYPE_EN.get(k, k): v
                                 for k, v in sorted(counts.items(), key=lambda x: -x[1])}),
            "warnings":     " | ".join(warnings) if warnings else "",
        })
    summary_df   = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "subtype_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  Subtype summary saved: {summary_path}")

    # 10. One-vs-Rest binary classifiers (each subtype vs others)
    _run_ovr_binary(X_raw, labels, sample_ids, observed_classes, result_dir, config)

    return results, summary_df


def _build_sample_date_map():
    """
    Build {sample_id: collection_date} keyed by fq_R1 prefix (same as feature matrix).

    Two-step join:
      1. plasma_exp.xlsx  : case_id → collection_date
      2. sample_info_excel: fq_R1 prefix (sample_id) → case_id
    Both tables store integer case_ids as Excel floats ("640849.0"); normalise to "640849".

    Returns an empty dict on any failure (caller should check and skip).
    """
    from utils.config import load_config as _load_config
    _cfg = _load_config()

    def _norm_id(raw):
        s = str(raw).strip()
        try:
            return str(int(float(s)))
        except (ValueError, OverflowError):
            return s

    from methylation.sample_grouping import _resolve_excel
    exp_filename = _cfg.get("data", {}).get("plasma_exp_excel", "plasma_exp.xlsx")
    exp_path = _resolve_excel(exp_filename) if exp_filename else None
    if not exp_path or not os.path.exists(exp_path):
        logger.warning(
            f"[subtype-temporal] plasma collection-date Excel '{exp_filename}' not found "
            f"(searched data/{exp_filename}, {exp_filename}) — skipping"
        )
        return {}

    try:
        df_exp = pd.read_excel(exp_path)
    except Exception as e:
        logger.warning(f"[subtype-temporal] Failed to read {exp_path}: {e}")
        return {}

    df_exp.columns = df_exp.columns.str.strip().str.lower()
    if "case_id" not in df_exp.columns or "datetime" not in df_exp.columns:
        logger.warning("[subtype-temporal] plasma_exp.xlsx missing required columns (case_id, datetime)")
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
        logger.warning("[subtype-temporal] No valid dates parsed from plasma_exp.xlsx")
        return {}

    excel_file = _cfg.get("data", {}).get("sample_info_excel", "")
    if not excel_file or not os.path.exists(excel_file):
        logger.warning(f"[subtype-temporal] sample_info_excel not found: {excel_file!r}")
        return {}

    sid_to_case = {}
    try:
        xl = pd.ExcelFile(excel_file)
        for sheet in xl.sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet)
            df.columns = df.columns.str.strip()
            if "fq_R1" not in df.columns or "case_id" not in df.columns:
                continue
            for _, row in df.iterrows():
                cid = _norm_id(row.get("case_id", ""))
                r1  = str(row.get("fq_R1", "")).strip()
                sid = r1.split(".")[0]
                if cid and cid != "nan" and sid and sid != "nan":
                    sid_to_case[sid] = cid
    except Exception as e:
        logger.warning(f"[subtype-temporal] Failed to read sample_info_excel: {e}")
        return {}

    sample_date_map = {
        sid: case_date[cid]
        for sid, cid in sid_to_case.items()
        if cid in case_date
    }
    logger.info(
        f"[subtype-temporal] date map: {len(case_date)} case_ids in plasma_exp.xlsx, "
        f"{len(sid_to_case)} sample_ids in sample_info_excel, "
        f"{len(sample_date_map)} matched → {list(sample_date_map.keys())[:5]}"
    )
    return sample_date_map


def run_subtype_temporal_validation(result_dir, config):
    """
    Supplementary temporal split validation for subtype classification (plasma mode only).

    Reads data/plasma_exp.xlsx (columns: case_id, datetime) joined with sample_info_excel
    to map case_id → sample_id.
    Splits at 2024-07-01: datetime < cutoff → train; datetime >= cutoff → validation.
    Trains RF and SVM once on the train set, evaluates on the validation set.

    Outputs:
        subtype_temporal_predictions.csv — per-sample predictions for validation set
        subtype_temporal_summary.csv     — accuracy per classifier
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

    tv_cfg = (config or {}).get("temporal_validate", {})
    CUTOFF = pd.Timestamp(tv_cfg.get("cutoff_date", "2024-07-01"))

    # 1-2. Build {sample_id: collection_date} via plasma_exp.xlsx + sample_info_excel bridge
    date_map = _build_sample_date_map()
    if not date_map:
        return None, None

    # 3. Load subtype map and feature matrix
    subtype_map = _build_subtype_map(config)
    if not subtype_map:
        logger.warning("[subtype-temporal] No subtype labels found — skipping")
        return None, None

    try:
        X_all, all_ids, feature_names = _load_feature_matrix(result_dir)
    except FileNotFoundError as e:
        logger.warning(f"[subtype-temporal] {e}")
        return None, None

    X_all, feature_names = _apply_panel_filter(X_all, feature_names, result_dir)

    # 4. Filter to samples with both a subtype label and a date
    indices, labels_raw, sample_ids = [], [], []
    for i, sid in enumerate(all_ids):
        sub = subtype_map.get(sid)
        if sub is None or sub not in SUBTYPE_CLASSES:
            continue
        if sid not in date_map:
            continue
        indices.append(i)
        labels_raw.append(sub)
        sample_ids.append(sid)

    if len(indices) == 0:
        logger.warning("[subtype-temporal] No samples with both subtype labels and dates — skipping")
        return None, None

    X_sub = X_all[indices]
    dates = [date_map[sid] for sid in sample_ids]

    # 5. Partition by date
    train_mask = [d < CUTOFF for d in dates]
    test_mask  = [d >= CUTOFF for d in dates]
    train_idx  = [i for i, m in enumerate(train_mask) if m]
    test_idx   = [i for i, m in enumerate(test_mask)  if m]

    logger.info(f"[subtype-temporal] Split: {len(train_idx)} train  {len(test_idx)} validation")

    if len(train_idx) < 5 or len(test_idx) == 0:
        logger.warning("[subtype-temporal] Insufficient samples for temporal split — skipping")
        return None, None

    X_train = X_sub[train_idx]
    X_test  = X_sub[test_idx]
    y_train_raw = [labels_raw[i] for i in train_idx]
    y_test_raw  = [labels_raw[i] for i in test_idx]
    test_ids    = [sample_ids[i]  for i in test_idx]
    test_dates  = [dates[i]       for i in test_idx]

    # Encode labels on train set only; skip test samples with unseen classes
    le = LabelEncoder()
    le.fit(y_train_raw)
    train_classes = set(y_train_raw)
    valid_test_mask = [l in train_classes for l in y_test_raw]
    if not any(valid_test_mask):
        logger.warning("[subtype-temporal] No test samples have classes seen in train — skipping")
        return None, None
    if not all(valid_test_mask):
        n_skip = sum(1 for m in valid_test_mask if not m)
        logger.warning(f"[subtype-temporal] {n_skip} test sample(s) have unseen subtype classes — excluded")

    X_test_filt   = X_test[[i for i, m in enumerate(valid_test_mask) if m]]
    y_test_filt   = [y_test_raw[i]  for i, m in enumerate(valid_test_mask) if m]
    test_ids_filt = [test_ids[i]    for i, m in enumerate(valid_test_mask) if m]
    test_dates_f  = [test_dates[i]  for i, m in enumerate(valid_test_mask) if m]

    y_train_enc = le.transform(y_train_raw)
    y_test_enc  = le.transform(y_test_filt)
    n_classes   = len(le.classes_)
    en_labels   = [SUBTYPE_EN.get(c, c) for c in le.classes_]

    # 6. Impute (train only) then scale
    X_train_imp, col_medians = _impute_fit(X_train)
    X_test_imp               = _impute_transform(X_test_filt, col_medians)
    scaler    = StandardScaler()
    X_train_sc = scaler.fit_transform(X_train_imp)
    X_test_sc  = scaler.transform(X_test_imp)

    ml_cfg = (config or {}).get("ml", {})
    classifiers = {
        "RF":  RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "SVM": SVC(
            kernel="rbf",
            C=ml_cfg.get("svm_C", 1.0),
            probability=True,
            class_weight="balanced",
            random_state=42,
        ),
    }

    # 7. Train and evaluate
    results = {}
    for clf_name, clf in classifiers.items():
        clf.fit(X_train_sc, y_train_enc)
        y_prob = clf.predict_proba(X_test_sc)   # (n_test, n_train_classes)
        # Align probability vector to full n_classes (handles class subset edge case)
        clf_classes = list(clf.classes_)
        y_prob_full = np.zeros((len(y_test_enc), n_classes))
        for j, cls_idx in enumerate(clf_classes):
            y_prob_full[:, cls_idx] = y_prob[:, j]

        y_pred = clf.predict(X_test_sc)
        acc    = float(accuracy_score(y_test_enc, y_pred))
        cm     = confusion_matrix(y_test_enc, y_pred, labels=list(range(n_classes)))
        logger.info(f"[subtype-temporal] {clf_name}: accuracy={acc:.3f}")
        logger.info("\n" + classification_report(
            y_test_enc, y_pred,
            labels=list(range(n_classes)),
            target_names=en_labels,
            zero_division=0,
        ))
        results[clf_name] = {
            "y_true": y_test_enc,
            "y_prob": y_prob_full,
            "y_pred": y_pred,
            "accuracy": acc,
            "confusion_matrix": cm,
        }

    # 8. Save predictions
    pred_rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(test_ids_filt):
            row = {
                "sample_id":         sid,
                "true_subtype":      le.inverse_transform([res["y_true"][i]])[0],
                "classifier":        clf_name,
                "predicted_subtype": le.inverse_transform([res["y_pred"][i]])[0],
                "confidence":        round(float(res["y_prob"][i].max()), 4),
                "correct":           bool(res["y_true"][i] == res["y_pred"][i]),
                "split":             "validation",
                "collection_date":   test_dates_f[i].strftime("%Y%m%d"),
            }
            for j, cls in enumerate(le.classes_):
                col = f"prob_{SUBTYPE_EN.get(cls, cls[:6])}"
                row[col] = round(float(res["y_prob"][i, j]), 4)
            pred_rows.append(row)

    pred_df   = pd.DataFrame(pred_rows)
    pred_path = os.path.join(result_dir, "subtype_temporal_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"[subtype-temporal] Predictions saved: {pred_path}")

    # 9. Save summary
    summary_rows = []
    for clf_name, res in results.items():
        counts = Counter(y_train_raw)
        summary_rows.append({
            "classifier":   clf_name,
            "n_train":      len(train_idx),
            "n_validation": len(test_ids_filt),
            "accuracy":     round(res["accuracy"], 4),
            "n_subtypes":   n_classes,
            "class_counts_train": str({SUBTYPE_EN.get(k, k): v
                                       for k, v in sorted(counts.items(), key=lambda x: -x[1])}),
            "cutoff_date":  str(CUTOFF.date()),
        })
    summary_df   = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "subtype_temporal_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"[subtype-temporal] Summary saved: {summary_path}")

    return results, summary_df
