"""imeta_utils.py — auto-generate manuscript publication figures and statistical tables.

Called during Step 20 (report). All public entry points are optional: missing input
files log a warning and the figure / table is silently skipped.

Module organisation:
    Configuration         module-level constants (model groups, colours, plot params)
    I/O helpers           _read_csv, _read_config, _savefig, _mpl
    Pipeline context      _PipelineContext, _load_context — cohort metadata read once
    DeLong inference      _fastdelong_cov, _delong_formal, _delong_bootstrap
    Operating point       _mcnemar_operating — paired sensitivity / specificity
    Figure builders       _fig3_combined, _fig5_temporal,
                          _fig2b_panel_annotation, _fig6b_blood_dist
    Public API            generate_imeta_figures, build_imeta_section
"""

import os
import sys
import logging
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

# Classifier groups — matches Step 15 eleven-model evaluation
CFDNA_ONLY    = ["LR", "SVM", "RF", "XGB", "LGBM"]
CA125_ALONE   = "CA-125"
FUSION_MODELS = [f"{m}+CA125" for m in CFDNA_ONLY]
ALL_MODELS    = CFDNA_ONLY + [CA125_ALONE] + FUSION_MODELS

# Colourblind-safe palette; fusion models share their cfDNA-only base colour
_CFDNA_COLORS = {
    "LR":   "#4393c3",
    "SVM":  "#d6604d",
    "RF":   "#4dac26",
    "XGB":  "#8073ac",
    "LGBM": "#e08214",
}
CLF_COLORS = {
    **_CFDNA_COLORS,
    CA125_ALONE: "#000000",
    **{f"{m}+CA125": _CFDNA_COLORS[m] for m in CFDNA_ONLY},
}

# Plot parameters
BOOTSTRAP_N     = 1000
DELONG_BOOT_N   = 2000
FIG_DPI         = 300
FIG_FACECOLOR   = "white"
DEFAULT_SEED    = 42

# Thresholds
BLOOD_BETA_THRESHOLD = 0.10


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_csv(path):
    """Return a DataFrame or None if missing/unreadable."""
    try:
        import pandas as pd
        if not os.path.exists(path):
            return None
        return pd.read_csv(path, encoding="utf-8")
    except Exception as exc:
        logger.warning("[imeta] read %s: %s", path, exc)
        return None


def _read_config(config_path: str = "config.toml") -> dict:
    """Read config.toml from CWD; return {} if absent or unparseable."""
    try:
        if sys.version_info >= (3, 11):
            import tomllib
            with open(config_path, "rb") as f:
                return tomllib.load(f)
        else:
            import toml
            return toml.load(config_path)
    except Exception:
        return {}


def _savefig(fig, path: str) -> None:
    import matplotlib.pyplot as plt
    fig.savefig(path, dpi=FIG_DPI, bbox_inches="tight", facecolor=FIG_FACECOLOR)
    plt.close(fig)


def _mpl():
    """Return pyplot with Agg backend; safe to call multiple times."""
    import matplotlib
    try:
        matplotlib.use("Agg")
    except Exception:
        pass
    import matplotlib.pyplot as plt
    return plt


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline context: cohort metadata derived from result CSVs
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class _PipelineContext:
    """Cohort metadata read once and shared across figure builders.

    All numeric defaults are zero so that missing inputs render captions like
    ``n=0`` rather than crashing; callers should still check ``n_total > 0``
    before assuming results are meaningful.
    """
    # Binary classification cohort
    n_total:     int = 0
    n_benign:    int = 0
    n_malignant: int = 0

    # CV mode
    cv_folds:    int = 5         # 0 = LOO-CV; k > 0 = stratified k-fold
    cv_label:    str = "5-fold CV"

    # Panel composition
    panel_size:  int = 0
    in_88_count: int = 0

    # Temporal validation
    tv_train:     int = 0
    tv_total:     int = 0
    tv_benign:    int = 0
    tv_malignant: int = 0
    tv_cutoff:    str = "2024-01-01"

    # External validation (450K matched-probe + blood baseline)
    ext_matched:  int = 0
    ext_pass:     int = 0
    ext_fail:     int = 0

    def cv_caption(self) -> str:
        return (f"{self.cv_label}, n={self.n_total}: "
                f"{self.n_benign} benign, {self.n_malignant} malignant")

    def tv_caption(self) -> str:
        return (f"train n={self.tv_train}, "
                f"validation n={self.tv_total} "
                f"({self.tv_benign} benign / {self.tv_malignant} malignant); "
                f"cutoff {self.tv_cutoff}")


def _load_context(result_dir: str) -> _PipelineContext:
    """Build a context object from result CSVs + config.toml."""
    ctx = _PipelineContext()

    # CV strategy from config (default 5-fold)
    cfg = _read_config()
    cv_folds = int(cfg.get("ml", {}).get("cv_folds", 5))
    ctx.cv_folds = cv_folds
    ctx.cv_label = "LOO-CV" if cv_folds == 0 else f"{cv_folds}-fold CV"

    # Binary classifier cohort
    clf_sum = _read_csv(os.path.join(result_dir, "classifier_summary.csv"))
    if clf_sum is not None and not clf_sum.empty:
        r = clf_sum.iloc[0]
        ctx.n_benign    = int(r.get("n_benign",    0) or 0)
        ctx.n_malignant = int(r.get("n_malignant", 0) or 0)
        ctx.n_total     = int(r.get("n_used", ctx.n_benign + ctx.n_malignant) or
                              ctx.n_benign + ctx.n_malignant)

    # Panel size: panel_cpg.csv lists the ACTUAL selected panel; panel_annotation.csv
    # may include additional candidate sites. Always trust panel_cpg.csv when present.
    panel_cpg = _read_csv(os.path.join(result_dir, "panel_cpg.csv"))
    panel_set = None
    if panel_cpg is not None and not panel_cpg.empty:
        ctx.panel_size = len(panel_cpg)
        if "feature" in panel_cpg.columns:
            panel_set = set(panel_cpg["feature"].astype(str))

    annot = _read_csv(os.path.join(result_dir, "panel_annotation.csv"))
    if annot is not None and not annot.empty:
        # Restrict to the final panel when possible, so the in-88 count reflects
        # the panel actually used downstream rather than the broader annotation table.
        if panel_set is not None and "feature" in annot.columns:
            annot = annot[annot["feature"].astype(str).isin(panel_set)]
        if not annot.empty and "in_88_regions" in annot.columns:
            ctx.in_88_count = int(
                (annot["in_88_regions"].astype(str).str.lower() == "true").sum()
            )
        if ctx.panel_size == 0:
            ctx.panel_size = len(annot)

    # Temporal validation
    tv_sum = _read_csv(os.path.join(result_dir, "temporal_validation_summary.csv"))
    if tv_sum is not None and not tv_sum.empty:
        r = tv_sum.iloc[0]
        ctx.tv_train  = int(r.get("n_train",      0) or 0)
        ctx.tv_total  = int(r.get("n_validation", 0) or 0)
        ctx.tv_cutoff = str(r.get("cutoff_date",  ctx.tv_cutoff))

        # benign / malignant breakdown is not stored in the summary CSV — recover
        # by counting unique sample_ids in the predictions file
        tv_pred = _read_csv(os.path.join(result_dir, "temporal_validation_predictions.csv"))
        if tv_pred is not None and "sample_id" in tv_pred.columns:
            uniq = tv_pred.drop_duplicates(subset="sample_id")
            ctx.tv_benign    = int((uniq["true_label"] == "benign").sum())
            ctx.tv_malignant = int((uniq["true_label"] == "malignant").sum())

    # External validation: 450K matched probes + leukocyte baseline pass.
    # Restrict to final-panel sites — the broader candidate set in
    # external_validation.csv may include probes that did not survive selection.
    ext = _read_csv(os.path.join(result_dir, "external_validation.csv"))
    if ext is not None and "probe_matched" in ext.columns:
        if panel_set is not None and "feature" in ext.columns:
            ext = ext[ext["feature"].astype(str).isin(panel_set)]
        matched = ext[ext["probe_matched"].astype(str).str.lower() == "true"]
        ctx.ext_matched = len(matched)
        if "blood_pass" in matched.columns:
            ctx.ext_pass = int(
                (matched["blood_pass"].astype(str).str.lower() == "true").sum()
            )
        ctx.ext_fail = ctx.ext_matched - ctx.ext_pass

    return ctx


# ─────────────────────────────────────────────────────────────────────────────
# Classifier prediction extraction
# ─────────────────────────────────────────────────────────────────────────────

def _extract_clf_data(df, models: Iterable[str]) -> dict:
    """Group a long-format predictions DataFrame by classifier.

    Returns ``{name: {"y_true": np.ndarray, "y_prob": np.ndarray}}`` for each
    classifier in ``models`` that has rows in the DataFrame.
    """
    data = {}
    for clf in models:
        sub = df[df["classifier"] == clf]
        if sub.empty:
            continue
        data[clf] = {
            "y_true": (sub["true_label"] == "malignant").astype(int).values,
            "y_prob": sub["prob_malignant"].astype(float).values,
        }
    return data


def _build_aligned_score_matrix(df, models: Iterable[str]):
    """Build a (K, n) score matrix with one row per classifier.

    Requires identical sample order across all classifiers; returns
    ``(y_true, score_matrix, model_names)`` or ``(None, None, None)`` on
    alignment failure.
    """
    ref_clf = None
    ref_ids = None
    aligned = {}
    for clf in models:
        sub = df[df["classifier"] == clf].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("sample_id").reset_index(drop=True)
        if ref_ids is None:
            ref_ids = sub["sample_id"].values
            ref_clf = clf
        elif not np.array_equal(sub["sample_id"].values, ref_ids):
            logger.warning(
                "[delong] sample_id mismatch between %s and %s — skipping",
                ref_clf, clf,
            )
            return None, None, None
        aligned[clf] = sub

    names = list(aligned.keys())
    if not names:
        return None, None, None

    y_true = (aligned[names[0]]["true_label"] == "malignant").astype(int).values
    score_matrix = np.vstack(
        [aligned[c]["prob_malignant"].astype(float).values for c in names]
    )
    return y_true, score_matrix, names


# ─────────────────────────────────────────────────────────────────────────────
# DeLong 1988 — formal variance estimation (fastDeLong algorithm)
# ─────────────────────────────────────────────────────────────────────────────

def _fastdelong_cov(y_true, score_matrix):
    """fastDeLong covariance matrix for K classifiers.

    y_true       : 1-D array of 0/1 labels, length n
    score_matrix : (K, n) array — row k is classifier k's probability scores

    Returns
    -------
    aucs : (K,) AUC for each classifier
    cov  : (K, K) covariance matrix of the K AUC estimates
    """
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n_pos, n_neg = len(pos_idx), len(neg_idx)

    K = score_matrix.shape[0]
    V10 = np.empty((K, n_pos))   # placement values for positives
    V01 = np.empty((K, n_neg))   # placement values for negatives

    for k in range(K):
        ps = score_matrix[k, pos_idx]   # (n_pos,)
        ns = score_matrix[k, neg_idx]   # (n_neg,)
        cmp10 = ns[np.newaxis, :] <  ps[:, np.newaxis]
        tie10 = ns[np.newaxis, :] == ps[:, np.newaxis]
        V10[k] = cmp10.mean(axis=1) + 0.5 * tie10.mean(axis=1)
        cmp01 = ps[np.newaxis, :] >  ns[:, np.newaxis]
        tie01 = ps[np.newaxis, :] == ns[:, np.newaxis]
        V01[k] = cmp01.mean(axis=1) + 0.5 * tie01.mean(axis=1)

    aucs = V10.mean(axis=1)

    # Structural-component covariance (DeLong 1988, eq. 7):
    #   Cov(AUC_a, AUC_b) = Cov_sample(V10_a, V10_b)/n_pos
    #                     + Cov_sample(V01_a, V01_b)/n_neg
    S10 = np.cov(V10, ddof=1) if n_pos > 1 else np.zeros((K, K))
    S01 = np.cov(V01, ddof=1) if n_neg > 1 else np.zeros((K, K))
    if K == 1:
        S10 = np.atleast_2d(S10)
        S01 = np.atleast_2d(S01)
    cov = S10 / n_pos + S01 / n_neg

    return aucs, cov


def _delong_formal(clf_pred_path: str, out_csv_path: str,
                   models: Optional[Iterable[str]] = None):
    """Formal pairwise AUC comparison using the fastDeLong variance estimator.

    Writes a normal Z-test (two-sided) per classifier pair.
    """
    import pandas as pd
    from scipy import stats

    df = _read_csv(clf_pred_path)
    if df is None:
        return None

    y_true, score_matrix, clfs = _build_aligned_score_matrix(
        df, models or ALL_MODELS
    )
    if y_true is None or len(clfs) < 2:
        return None

    aucs, cov = _fastdelong_cov(y_true, score_matrix)

    rows = []
    for i, ca in enumerate(clfs):
        for j in range(i + 1, len(clfs)):
            cb = clfs[j]
            d = aucs[i] - aucs[j]
            var = cov[i, i] + cov[j, j] - 2 * cov[i, j]
            se = float(np.sqrt(max(var, 0.0)))
            z = d / se if se > 0 else float("nan")
            p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
            rows.append({
                "classifier_A": ca,
                "classifier_B": cb,
                "auc_A":        round(float(aucs[i]), 4),
                "auc_B":        round(float(aucs[j]), 4),
                "delta_auc":    round(float(d), 4),
                "se_delta":     round(se, 6),
                "z_stat":       round(float(z), 4) if np.isfinite(z) else "NA",
                "p_delong":     round(float(p), 4) if np.isfinite(p) else "NA",
                "significant":  "Yes" if (np.isfinite(p) and p < 0.05) else "No",
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv_path, index=False, encoding="utf-8")
    logger.info(
        "[delong] %d pairs tested across %d models → %s",
        len(rows), len(clfs), os.path.basename(out_csv_path),
    )
    return df_out


def _delong_bootstrap(clf_pred_path: str, out_csv_path: str,
                      n_boot: int = DELONG_BOOT_N,
                      models: Optional[Iterable[str]] = None):
    """Bootstrap pairwise AUC test (paired resampling)."""
    import pandas as pd
    from sklearn.metrics import roc_auc_score

    df = _read_csv(clf_pred_path)
    if df is None:
        return None

    data = _extract_clf_data(df, models or ALL_MODELS)
    clfs = list(data.keys())
    if len(clfs) < 2:
        return None

    n = len(data[clfs[0]]["y_true"])
    rng = np.random.default_rng(DEFAULT_SEED)
    rows = []

    for i, c1 in enumerate(clfs):
        for c2 in clfs[i + 1:]:
            y_true = data[c1]["y_true"]
            diffs = []
            for _ in range(n_boot):
                idx = rng.integers(0, n, n)
                if len(np.unique(y_true[idx])) < 2:
                    continue
                a1 = roc_auc_score(y_true[idx], data[c1]["y_prob"][idx])
                a2 = roc_auc_score(y_true[idx], data[c2]["y_prob"][idx])
                diffs.append(a1 - a2)
            diffs = np.array(diffs)
            auc1 = roc_auc_score(y_true, data[c1]["y_prob"])
            auc2 = roc_auc_score(y_true, data[c2]["y_prob"])
            obs  = auc1 - auc2
            if len(diffs):
                p_two_sided = min(np.mean(diffs <= 0), np.mean(diffs >= 0)) * 2
                ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])
            else:
                p_two_sided, ci_lo, ci_hi = 1.0, float("nan"), float("nan")
            rows.append({
                "classifier_A": c1, "classifier_B": c2,
                "auc_A": round(auc1, 4), "auc_B": round(auc2, 4),
                "delta_auc": round(obs, 4),
                "ci95_lo": round(float(ci_lo), 4),
                "ci95_hi": round(float(ci_hi), 4),
                "p_bootstrap": round(float(p_two_sided), 4),
                "significant": "Yes" if p_two_sided < 0.05 else "No",
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv_path, index=False, encoding="utf-8")
    return df_out


# ─────────────────────────────────────────────────────────────────────────────
# McNemar — paired operating-point comparison (sensitivity / specificity)
# ─────────────────────────────────────────────────────────────────────────────
#
# DeLong tests the ROC curve as a whole and is threshold-free. It therefore
# cannot speak to sensitivity or specificity, which are properties of one
# operating point. Those are paired binary outcomes on the same samples, so the
# correct test is McNemar's on the discordant pairs. The exact binomial form is
# used rather than the chi-square approximation because the discordant counts
# here are small (5-18), where the approximation is unreliable.
#
# Reported for every comparison attempted, significant or not: cherry-picking
# the operating point after seeing the result is precisely the kind of inflation
# this pipeline exists to measure.

def _mcnemar_operating(clf_pred_path: str, out_csv_path: str,
                       models: Optional[Iterable[str]] = None):
    """Exact McNemar tests on the classification operating point.

    Compares CA-125 alone against each fusion model, and each cfDNA-only model
    against its own fusion counterpart, separately within the benign stratum
    (specificity) and the malignant stratum (sensitivity).
    """
    import pandas as pd
    from scipy import stats

    df = _read_csv(clf_pred_path)
    if df is None or "predicted" not in df.columns:
        return None

    wanted = list(models or ALL_MODELS)
    by_clf = {}
    for clf in wanted:
        sub = df[df["classifier"] == clf]
        if sub.empty:
            continue
        by_clf[clf] = {
            r["sample_id"]: (r["true_label"], r["predicted"])
            for _, r in sub.iterrows()
        }
    if CA125_ALONE not in by_clf:
        logger.warning("[mcnemar] %s absent from predictions — skipping",
                       CA125_ALONE)
        return None

    # (reference, comparator) pairs: the incumbent marker against each fusion
    # model, then each cfDNA-only model against the fusion that extends it.
    pairs = [(CA125_ALONE, f) for f in FUSION_MODELS]
    pairs += [(m, f"{m}+CA125") for m in CFDNA_ONLY]

    strata = [("specificity", "benign"), ("sensitivity", "malignant")]

    rows = []
    for ref, cmp_ in pairs:
        if ref not in by_clf or cmp_ not in by_clf:
            continue
        A, B = by_clf[ref], by_clf[cmp_]
        for metric, label in strata:
            ids = [s for s in A if s in B and A[s][0] == label]
            if not ids:
                continue
            # b10: reference correct, comparator wrong; b01: the reverse.
            b10 = sum(1 for s in ids if A[s][1] == label and B[s][1] != label)
            b01 = sum(1 for s in ids if A[s][1] != label and B[s][1] == label)
            n_disc = b10 + b01
            p = (float(stats.binomtest(min(b10, b01), n_disc, 0.5).pvalue)
                 if n_disc else float("nan"))
            rate_ref = sum(1 for s in ids if A[s][1] == label) / len(ids)
            rate_cmp = sum(1 for s in ids if B[s][1] == label) / len(ids)
            rows.append({
                "reference":       ref,
                "comparator":      cmp_,
                "metric":          metric,
                "stratum":         label,
                "n":               len(ids),
                "rate_reference":  round(rate_ref, 4),
                "rate_comparator": round(rate_cmp, 4),
                "delta":           round(rate_cmp - rate_ref, 4),
                "discordant_ref_only": b10,
                "discordant_cmp_only": b01,
                "p_mcnemar_exact": round(p, 4) if np.isfinite(p) else "NA",
                "significant":     "Yes" if (np.isfinite(p) and p < 0.05) else "No",
            })

    if not rows:
        return None

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv_path, index=False, encoding="utf-8")
    logger.info(
        "[mcnemar] %d operating-point tests → %s",
        len(rows), os.path.basename(out_csv_path),
    )
    return df_out


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Binary classification performance (merged: ROC + headline metrics)
# ─────────────────────────────────────────────────────────────────────────────

_METRIC_SPEC = [
    ("auc",         "auc_ci95_lo",  "auc_ci95_hi",  "AUC"),
    ("accuracy",    "acc_ci95_lo",  "acc_ci95_hi",  "Accuracy"),
    ("sensitivity", "sens_ci95_lo", "sens_ci95_hi", "Sensitivity"),
    ("specificity", "spec_ci95_lo", "spec_ci95_hi", "Specificity"),
]

def _fig3_combined(clf_pred_path: str, clf_sum_path: str, out_path: str,
                   ctx: _PipelineContext) -> Optional[bool]:
    """Merged version of the previous Figure 3 (ROC) + Figure 4 (bar chart).

    Left panel — ROC curves for all eleven models on a single axis with line-style
    distinction: cfDNA-only solid, CA-125 alone dashed black, fusion (cfDNA +
    CA-125) dotted in the cfDNA base colour.

    Right panel — grouped bar chart of the four headline metrics (AUC, Accuracy,
    Sensitivity, Specificity) for the same eleven models with 95% bootstrap CIs.
    Fusion bars are drawn with a `//` hatch so they remain distinguishable from
    their cfDNA-only counterpart of the same base colour.
    """
    plt = _mpl()
    pred_df = _read_csv(clf_pred_path)
    sum_df  = _read_csv(clf_sum_path)
    if pred_df is None or sum_df is None:
        return None

    fig, (ax_roc, ax_bar) = plt.subplots(1, 2, figsize=(16, 6.5))

    # ── Left panel: ROC for all 11 models on one axis ──
    from sklearn.metrics import roc_curve, auc as sk_auc
    rng = np.random.default_rng(DEFAULT_SEED)

    def _style(clf):
        if clf == CA125_ALONE:
            return CLF_COLORS[CA125_ALONE], "--", 1.5
        if clf.endswith("+CA125"):
            return CLF_COLORS[clf], ":", 2.2
        return CLF_COLORS[clf], "-", 2.0

    plotted = 0
    for clf in ALL_MODELS:
        sub = pred_df[pred_df["classifier"] == clf]
        if sub.empty:
            continue
        y_true = (sub["true_label"] == "malignant").astype(int).values
        y_prob = sub["prob_malignant"].astype(float).values
        if len(set(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        auc_val = sk_auc(fpr, tpr)
        boot = []
        for _ in range(BOOTSTRAP_N):
            idx = rng.integers(0, len(y_true), len(y_true))
            if len(np.unique(y_true[idx])) < 2:
                continue
            f, t, _ = roc_curve(y_true[idx], y_prob[idx])
            boot.append(sk_auc(f, t))
        ci_lo, ci_hi = (np.percentile(boot, [2.5, 97.5]) if boot
                        else (float("nan"), float("nan")))
        c, ls, lw = _style(clf)
        ax_roc.plot(fpr, tpr, color=c, lw=lw, linestyle=ls,
                    label=f"{clf}  AUC={auc_val:.3f} ({ci_lo:.3f}–{ci_hi:.3f})")
        plotted += 1

    if plotted == 0:
        plt.close(fig)
        return None

    ax_roc.plot([0, 1], [0, 1], "k:", lw=1, alpha=0.4)
    ax_roc.set_xlim([0, 1]); ax_roc.set_ylim([0, 1.02])
    ax_roc.set_xlabel("1 − Specificity (False Positive Rate)", fontsize=11)
    ax_roc.set_ylabel("Sensitivity (True Positive Rate)", fontsize=11)
    ax_roc.set_title("A. ROC curves (all eleven models)",
                     fontsize=11, fontweight="bold")
    ax_roc.legend(loc="lower right", fontsize=7.5, framealpha=0.92, ncol=1)
    ax_roc.grid(alpha=0.25)

    # ── Right panel: grouped bar chart of headline metrics ──
    clfs_present = [c for c in ALL_MODELS if c in sum_df["classifier"].values]
    if not clfs_present:
        plt.close(fig)
        return None

    x = np.arange(len(_METRIC_SPEC))
    n = len(clfs_present)
    width = 0.9 / n
    offsets = np.linspace(-(n - 1) / 2, (n - 1) / 2, n) * width

    for i, clf in enumerate(clfs_present):
        row = sum_df[sum_df["classifier"] == clf].iloc[0]
        vals = [float(row.get(m, 0)) for m, *_ in _METRIC_SPEC]
        lo = [float(row.get(l, v)) for (m, l, h, _), v in zip(_METRIC_SPEC, vals)]
        hi = [float(row.get(h, v)) for (m, l, h, _), v in zip(_METRIC_SPEC, vals)]
        yerr = [[v - lo_i for v, lo_i in zip(vals, lo)],
                [hi_i - v for v, hi_i in zip(vals, hi)]]
        is_fusion = clf.endswith("+CA125")
        ax_bar.bar(x + offsets[i], vals, width, yerr=yerr,
                   color=CLF_COLORS.get(clf, "#777"),
                   alpha=0.85, capsize=1.5,
                   label=clf,
                   hatch="//" if is_fusion else None,
                   edgecolor="#222" if is_fusion else "none",
                   linewidth=0.4,
                   error_kw={"elinewidth": 0.8, "ecolor": "#555"})

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([m[-1] for m in _METRIC_SPEC], fontsize=11)
    ax_bar.set_ylim([0, 1.15])
    ax_bar.axhline(0.5, color="#bbb", lw=1, ls="--", alpha=0.6)
    ax_bar.set_title("B. Headline metrics (95% bootstrap CI)",
                     fontsize=11, fontweight="bold")
    ax_bar.set_ylabel("Score", fontsize=11)
    ax_bar.legend(loc="upper right", fontsize=7.5, ncol=3, framealpha=0.9)
    ax_bar.grid(axis="y", alpha=0.2)

    fig.suptitle(
        f"Figure 3 — Binary Classification Performance "
        f"({ctx.cv_caption()}; 95% bootstrap CI on AUC and metrics)",
        fontsize=12.5, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — Temporal validation: CV-vs-prospective AUCs + best-model scatter
# ─────────────────────────────────────────────────────────────────────────────

def _fig5_temporal(tv_sum_path: str, tv_pred_path: str, clf_sum_path: str,
                   out_path: str, ctx: _PipelineContext) -> Optional[bool]:
    import pandas as pd
    plt = _mpl()

    tv_sum  = _read_csv(tv_sum_path)
    tv_pred = _read_csv(tv_pred_path)
    clf_sum = _read_csv(clf_sum_path)
    if tv_sum is None or tv_pred is None:
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # ── Panel A: side-by-side CV vs prospective AUC bars (all 11 models) ──
    clfs = [c for c in ALL_MODELS if c in tv_sum["classifier"].values]
    x = np.arange(len(clfs))
    w = 0.36

    cv_aucs = []
    if clf_sum is not None:
        for c in clfs:
            row = clf_sum[clf_sum["classifier"] == c]
            cv_aucs.append(float(row["auc"].iloc[0]) if not row.empty else 0.0)
    else:
        cv_aucs = [0.0] * len(clfs)
    tv_aucs = [float(tv_sum[tv_sum["classifier"] == c]["auc"].iloc[0]) for c in clfs]

    b1 = ax1.bar(x - w / 2, cv_aucs, w,
                 label=f"{ctx.cv_label} (n={ctx.n_total})",
                 color="#4393c3", alpha=0.85)
    b2 = ax1.bar(x + w / 2, tv_aucs, w,
                 label=f"Prospective (n={ctx.tv_total})",
                 color="#d6604d", alpha=0.85)
    ax1.axhline(0.5, color="#bbb", lw=1, ls="--", alpha=0.6)
    ax1.set_xticks(x)
    ax1.set_xticklabels(clfs, fontsize=8.5, rotation=30, ha="right")
    ax1.set_ylim([0, 1.08])
    ax1.set_ylabel("AUC", fontsize=11)
    ax1.set_title(f"A. {ctx.cv_label} vs Prospective AUC",
                  fontsize=11, fontweight="bold")
    ax1.legend(fontsize=9, loc="upper left")
    ax1.grid(axis="y", alpha=0.2)
    for bar in list(b1) + list(b2):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012,
                 f"{bar.get_height():.3f}",
                 ha="center", va="bottom", fontsize=6.5, rotation=90)

    # ── Panel B: prediction scatter for the prospective best model ────────
    best_clf = tv_sum.loc[tv_sum["auc"].astype(float).idxmax(), "classifier"]
    sub = tv_pred[tv_pred["classifier"] == best_clf].copy()
    if sub.empty:
        first_clf = tv_pred["classifier"].iloc[0]
        sub = tv_pred[tv_pred["classifier"] == first_clf].copy()
        best_clf = first_clf

    sub["date"] = pd.to_datetime(sub["collection_date"].astype(str),
                                 format="%Y%m%d", errors="coerce")
    sub["prob"] = sub["prob_malignant"].astype(float)
    sub["correct"] = sub["correct"].astype(str).str.lower().isin(["true", "1"])

    for label, color, marker in [("malignant", "#d6604d", "o"),
                                 ("benign",    "#4393c3", "s")]:
        grp = sub[sub["true_label"] == label]
        if not grp.empty:
            ax2.scatter(grp["date"], grp["prob"],
                        c=color, marker=marker, s=70, alpha=0.85, zorder=3,
                        label=f"{label.capitalize()} (n={len(grp)})")

    wrong = sub[~sub["correct"]]
    if not wrong.empty:
        ax2.scatter(wrong["date"], wrong["prob"],
                    marker="x", s=110, c="black", linewidths=2, zorder=5,
                    label=f"Misclassified (n={len(wrong)})")

    ax2.axhline(0.5, color="#888", lw=1.2, ls="--", alpha=0.7,
                label="Decision boundary (0.5)")
    ax2.set_ylim([-0.07, 1.07])
    ax2.set_xlabel("Sample collection date", fontsize=11)
    ax2.set_ylabel(f"P(Malignant) — {best_clf}", fontsize=11)
    ax2.set_title(
        f"B. Prospective Predictions — best model: {best_clf} "
        f"(validation n={ctx.tv_total}: {ctx.tv_benign} benign, "
        f"{ctx.tv_malignant} malignant)",
        fontsize=10.5, fontweight="bold",
    )
    fig.autofmt_xdate(rotation=30)
    ax2.legend(fontsize=8, loc="center left")
    ax2.grid(alpha=0.2)

    fig.suptitle(
        f"Figure 5 — Temporal Prospective Validation ({ctx.tv_caption()})",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4B / 4C / 5B — Multimodal-fusion marginal-value figures
# ─────────────────────────────────────────────────────────────────────────────

# Metric columns shared by classifier_summary.csv and temporal_validation_summary.csv
_FUSION_METRICS = [
    ("auc",         "auc_ci95_lo",  "auc_ci95_hi",  "AUC"),
    ("accuracy",    "acc_ci95_lo",  "acc_ci95_hi",  "Accuracy"),
    ("sensitivity", "sens_ci95_lo", "sens_ci95_hi", "Sensitivity"),
    ("specificity", "spec_ci95_lo", "spec_ci95_hi", "Specificity"),
]


def _row_for(df, clf: str):
    """Return the single row matching ``clf`` or None if absent."""
    sub = df[df["classifier"] == clf]
    return sub.iloc[0] if not sub.empty else None


def _fig4b_fusion_delta(clf_sum_path: str, out_path: str,
                        ctx: _PipelineContext) -> Optional[bool]:
    """Paired fusion-vs-cfDNA-only Δ for each base classifier (4 metrics)."""
    plt = _mpl()
    df = _read_csv(clf_sum_path)
    if df is None:
        return None

    pairs = [(m, f"{m}+CA125") for m in CFDNA_ONLY]
    rows = []
    for base, fused in pairs:
        rb, rf = _row_for(df, base), _row_for(df, fused)
        if rb is None or rf is None:
            continue
        rec = {"base": base}
        for m, _, _, _ in _FUSION_METRICS:
            rec[f"{m}_base"]  = float(rb.get(m, np.nan))
            rec[f"{m}_fused"] = float(rf.get(m, np.nan))
            rec[f"{m}_delta"] = rec[f"{m}_fused"] - rec[f"{m}_base"]
        rows.append(rec)

    if not rows:
        return None

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4), sharey=True)
    bases = [r["base"] for r in rows]
    y = np.arange(len(bases))

    for ax, (m, _, _, label) in zip(axes, _FUSION_METRICS):
        deltas = np.array([r[f"{m}_delta"] for r in rows])
        colors = ["#2a9d8f" if d >= 0 else "#e63946" for d in deltas]
        ax.barh(y, deltas, color=colors, alpha=0.85,
                edgecolor="#333", linewidth=0.4)
        ax.axvline(0, color="#444", lw=0.8)
        for i, (d, b_val, f_val) in enumerate(zip(
                deltas,
                [r[f"{m}_base"]  for r in rows],
                [r[f"{m}_fused"] for r in rows])):
            sign = "+" if d >= 0 else ""
            ax.text(d + (0.003 if d >= 0 else -0.003), i,
                    f"{sign}{d:.3f}\n({b_val:.3f}→{f_val:.3f})",
                    va="center",
                    ha="left" if d >= 0 else "right",
                    fontsize=7.5, color="#222")
        ax.set_yticks(y)
        ax.set_yticklabels(bases, fontsize=9)
        ax.set_title(f"Δ {label}", fontsize=11, fontweight="bold")
        ax.set_xlabel(f"fused − base ({label})", fontsize=9)
        ax.grid(axis="x", alpha=0.25)
        x_pad = max(0.02, np.max(np.abs(deltas)) * 0.45)
        ax.set_xlim(-x_pad, x_pad + max(0, np.max(deltas)) * 0.6)

    axes[0].set_ylabel("Base classifier", fontsize=10)
    fig.suptitle(
        f"Figure 4B — Marginal value of CA-125 fusion ({ctx.cv_caption()}): "
        f"five paired (cfDNA-only, +CA125) comparisons",
        fontsize=12, fontweight="bold", y=1.04,
    )
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


def _fig4c_spec_gain(clf_sum_path: str, out_path: str,
                     ctx: _PipelineContext) -> Optional[bool]:
    """Δ specificity and Δ AUC for each cfDNA-only/fusion model vs CA-125 alone."""
    plt = _mpl()
    df = _read_csv(clf_sum_path)
    if df is None:
        return None

    ca = _row_for(df, CA125_ALONE)
    if ca is None:
        return None
    ca_spec = float(ca.get("specificity", np.nan))
    ca_auc  = float(ca.get("auc",         np.nan))
    if not (np.isfinite(ca_spec) and np.isfinite(ca_auc)):
        return None

    models = CFDNA_ONLY + FUSION_MODELS
    rows = []
    for clf in models:
        r = _row_for(df, clf)
        if r is None:
            continue
        rows.append({
            "clf":      clf,
            "is_fused": clf in FUSION_MODELS,
            "d_spec":   float(r.get("specificity", np.nan)) - ca_spec,
            "d_auc":    float(r.get("auc",         np.nan)) - ca_auc,
            "spec":     float(r.get("specificity", np.nan)),
            "auc":      float(r.get("auc",         np.nan)),
        })
    if not rows:
        return None

    # Order: cfDNA-only first (sorted by Δ spec ascending), then fusion (same order)
    rows.sort(key=lambda r: (r["is_fused"], r["d_spec"]))
    labels = [r["clf"] for r in rows]
    y      = np.arange(len(labels))

    fig, (ax_spec, ax_auc) = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)

    def _draw(ax, key_d: str, key_abs: str, ca_ref: float, x_label: str, title: str):
        deltas = np.array([r[key_d] for r in rows])
        for i, r in enumerate(rows):
            edge   = "#1d3557" if r["is_fused"] else "#bbbbbb"
            color  = CLF_COLORS.get(r["clf"], "#777")
            alpha  = 0.95 if r["is_fused"] else 0.65
            ax.barh(i, r[key_d], color=color, alpha=alpha,
                    edgecolor=edge, linewidth=1.2 if r["is_fused"] else 0.5)
            sign = "+" if r[key_d] >= 0 else ""
            ax.text(r[key_d] + (0.004 if r[key_d] >= 0 else -0.004), i,
                    f"{sign}{r[key_d]:.3f}\n({r[key_abs]:.3f})",
                    va="center",
                    ha="left" if r[key_d] >= 0 else "right",
                    fontsize=7.5, color="#222")
        ax.axvline(0, color="#444", lw=0.8)
        ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel(x_label, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.grid(axis="x", alpha=0.25)
        x_pad = max(0.05, np.max(np.abs(deltas)) * 0.45)
        ax.set_xlim(-x_pad, x_pad + max(0, np.max(deltas)) * 0.5)
        # Reference baseline annotation
        ax.text(0.99, -0.16,
                f"CA-125 baseline = {ca_ref:.3f}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=8.5, color="#555", style="italic")

    _draw(ax_spec, "d_spec", "spec", ca_spec,
          "Δ Specificity vs CA-125 alone",
          "A. Specificity gain over CA-125")
    _draw(ax_auc, "d_auc", "auc", ca_auc,
          "Δ AUC vs CA-125 alone",
          "B. AUC gain over CA-125")

    # Separator between cfDNA-only and fusion groups
    n_cfdna = sum(1 for r in rows if not r["is_fused"])
    if 0 < n_cfdna < len(rows):
        for ax in (ax_spec, ax_auc):
            ax.axhline(n_cfdna - 0.5, color="#1d3557", lw=0.6, ls="--", alpha=0.6)
    ax_spec.set_ylabel("Classifier  (fusion = dark border)", fontsize=10)

    fig.suptitle(
        f"Figure 4C — Per-model gain over CA-125 alone ({ctx.cv_caption()}); "
        f"fusion models bear a dark border",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


def _fig5b_two_regime(clf_sum_path: str, tv_sum_path: str, out_path: str,
                      ctx: _PipelineContext) -> Optional[bool]:
    """Side-by-side AUC and Specificity across CV and Temporal regimes (all 11 models).

    Highlights whether the fusion advantage is preserved when models are evaluated
    on a chronologically held-out cohort — the central claim of the manuscript.
    """
    plt = _mpl()
    cv = _read_csv(clf_sum_path)
    tv = _read_csv(tv_sum_path)
    if cv is None or tv is None:
        return None

    models = [m for m in ALL_MODELS
              if m in cv["classifier"].values and m in tv["classifier"].values]
    if not models:
        return None

    x = np.arange(len(models))
    width = 0.36

    fig, (ax_auc, ax_spec) = plt.subplots(1, 2, figsize=(16, 5.4), sharex=True)

    def _draw(ax, metric: str, lo_col: str, hi_col: str, ylabel: str, title: str):
        cv_v  = np.array([float(_row_for(cv, m).get(metric, np.nan)) for m in models])
        tv_v  = np.array([float(_row_for(tv, m).get(metric, np.nan)) for m in models])
        cv_lo = np.array([float(_row_for(cv, m).get(lo_col, v)) for m, v in zip(models, cv_v)])
        cv_hi = np.array([float(_row_for(cv, m).get(hi_col, v)) for m, v in zip(models, cv_v)])
        tv_lo = np.array([float(_row_for(tv, m).get(lo_col, v)) for m, v in zip(models, tv_v)])
        tv_hi = np.array([float(_row_for(tv, m).get(hi_col, v)) for m, v in zip(models, tv_v)])
        cv_err = [cv_v - cv_lo, cv_hi - cv_v]
        tv_err = [tv_v - tv_lo, tv_hi - tv_v]

        # Distinguish CV from Temporal by fill style; preserve per-model colour
        for i, m in enumerate(models):
            col = CLF_COLORS.get(m, "#777")
            ax.bar(i - width/2, cv_v[i], width,
                   yerr=[[cv_err[0][i]], [cv_err[1][i]]],
                   color=col, alpha=0.85, edgecolor="#222", linewidth=0.4,
                   error_kw={"elinewidth": 0.9, "ecolor": "#555"},
                   label=ctx.cv_label if i == 0 else None)
            ax.bar(i + width/2, tv_v[i], width,
                   yerr=[[tv_err[0][i]], [tv_err[1][i]]],
                   color=col, alpha=0.85, edgecolor="#222", linewidth=0.4,
                   hatch="///",
                   error_kw={"elinewidth": 0.9, "ecolor": "#555"},
                   label="Temporal validation" if i == 0 else None)

        ax.set_xticks(x)
        ax.set_xticklabels(models, rotation=35, ha="right", fontsize=9)
        ax.set_ylim(0, 1.08)
        ax.axhline(0.5, color="#bbb", lw=1, ls="--", alpha=0.5)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="upper left", fontsize=8.5, framealpha=0.9)

    _draw(ax_auc,  "auc",         "auc_ci95_lo",  "auc_ci95_hi",
          "AUC",         "A. AUC: CV vs Temporal")
    _draw(ax_spec, "specificity", "spec_ci95_lo", "spec_ci95_hi",
          "Specificity", "B. Specificity: CV vs Temporal")

    fig.suptitle(
        f"Figure 5B — Cross-regime consistency: {ctx.cv_label} (solid) vs "
        f"Temporal prospective validation (hatched), 95% CI",
        fontsize=12, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4D — DeLong forest plot (15 focused pairwise AUC comparisons)
# ─────────────────────────────────────────────────────────────────────────────

def _lookup_delong_pair(df, ca: str, cb: str):
    """Return (delta_oriented_as_A_minus_B, ci_lo, ci_hi, p_value_or_None, row).

    Handles both orientations: if the table stored (cb, ca) instead of (ca, cb),
    the sign of delta_auc / CI is flipped so the result is always A - B oriented.
    Bootstrap-mode tables include ``ci95_lo`` / ``ci95_hi`` and ``p_bootstrap``.
    Formal-mode tables omit CI columns and expose ``p_delong``.
    """
    if df is None or df.empty:
        return None
    same = df[(df["classifier_A"] == ca) & (df["classifier_B"] == cb)]
    flip = df[(df["classifier_A"] == cb) & (df["classifier_B"] == ca)]
    if not same.empty:
        r = same.iloc[0]
        sign = 1.0
    elif not flip.empty:
        r = flip.iloc[0]
        sign = -1.0
    else:
        return None

    delta = sign * float(r.get("delta_auc", np.nan))
    ci_lo_raw = r.get("ci95_lo", np.nan)
    ci_hi_raw = r.get("ci95_hi", np.nan)
    try:
        ci_lo = float(ci_lo_raw); ci_hi = float(ci_hi_raw)
        if sign < 0:
            ci_lo, ci_hi = -ci_hi, -ci_lo
    except (TypeError, ValueError):
        ci_lo, ci_hi = float("nan"), float("nan")

    p_val = None
    for col in ("p_bootstrap", "p_delong"):
        if col in df.columns:
            try:
                p_val = float(r.get(col))
            except (TypeError, ValueError):
                p_val = None
            break
    return delta, ci_lo, ci_hi, p_val


def _fig_delong_forest(delong_boot_path: str, delong_formal_path: str,
                       out_path: str, ctx: _PipelineContext) -> Optional[bool]:
    """Forest plot of 15 focused pairwise AUC comparisons (3 semantic groups).

    Groups:
      1. Each cfDNA-only vs CA-125 alone (5 rows)
      2. Each fusion vs CA-125 alone (5 rows)
      3. Each fusion vs its cfDNA-only counterpart (5 rows)

    Shows Δ AUC with 95% CI whiskers (from the bootstrap table) and annotates
    each row with bootstrap p and formal DeLong p side by side.
    """
    plt = _mpl()
    df_boot = _read_csv(delong_boot_path)
    df_form = _read_csv(delong_formal_path)
    if df_boot is None and df_form is None:
        return None

    # Build (group_label, label, ca, cb) for each of the 15 rows
    groups = [
        ("cfDNA vs CA-125",
         [(m, m, CA125_ALONE) for m in CFDNA_ONLY]),
        ("Fusion vs CA-125",
         [(f"{m}+CA125", f"{m}+CA125", CA125_ALONE) for m in CFDNA_ONLY]),
        ("Fusion vs base cfDNA",
         [(f"{m}+CA125 vs {m}", f"{m}+CA125", m) for m in CFDNA_ONLY]),
    ]

    rows = []
    for group_idx, (group_label, items) in enumerate(groups):
        for label, ca, cb in items:
            b = _lookup_delong_pair(df_boot, ca, cb) if df_boot is not None else None
            f = _lookup_delong_pair(df_form, ca, cb) if df_form is not None else None
            if b is None and f is None:
                continue
            delta, ci_lo, ci_hi, p_boot = b if b is not None else (
                f[0], float("nan"), float("nan"), None
            )
            p_form = f[3] if f is not None else None
            rows.append({
                "group_idx": group_idx,
                "group_label": group_label,
                "label": label,
                "delta": delta,
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "p_boot": p_boot,
                "p_form": p_form,
            })

    if not rows:
        return None

    # Layout: rows top-to-bottom in declared order; insert blank slots between groups
    y_positions = []
    y = 0.0
    last_group = None
    for r in rows:
        if last_group is not None and r["group_idx"] != last_group:
            y -= 1.0  # blank row between groups
        y -= 1.0
        y_positions.append(y)
        last_group = r["group_idx"]

    fig_height = max(5.0, 0.45 * (abs(y_positions[-1]) + 2))
    fig, ax = plt.subplots(figsize=(11.5, fig_height))

    # Colour by group
    group_colors = {0: "#4393c3", 1: "#8073ac", 2: "#4dac26"}

    # Draw zero reference and CI whiskers + point estimates
    ax.axvline(0.0, color="#444", lw=1.0, ls="--", alpha=0.7, zorder=1)
    for r, y in zip(rows, y_positions):
        col = group_colors[r["group_idx"]]
        d, lo, hi = r["delta"], r["ci_lo"], r["ci_hi"]
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [y, y], color=col, lw=2.0, alpha=0.85, zorder=2)
            ax.plot([lo, lo], [y - 0.18, y + 0.18], color=col, lw=2.0, alpha=0.85)
            ax.plot([hi, hi], [y - 0.18, y + 0.18], color=col, lw=2.0, alpha=0.85)
        ax.plot([d], [y], marker="o", markersize=8, color=col,
                markeredgecolor="black", markeredgewidth=0.6, zorder=3)

    # Compute x-axis range with padding for annotations on the right
    all_vals = []
    for r in rows:
        for v in (r["delta"], r["ci_lo"], r["ci_hi"]):
            if np.isfinite(v):
                all_vals.append(v)
    x_min = min(all_vals + [0.0]) - 0.03
    x_max = max(all_vals + [0.0]) + 0.03
    x_span = x_max - x_min
    x_anno = x_max + 0.04 * x_span
    ax.set_xlim(x_min, x_anno + 0.55 * x_span)

    # Row labels on left axis
    ax.set_yticks(y_positions)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=9.5)

    # Right-side annotations: Δ AUC [CI]   p_boot / p_DeLong   significance star
    def _fmt_p(p):
        if p is None or not np.isfinite(p):
            return "NA"
        return "<.001" if p < 0.001 else f"{p:.3f}"

    for r, y in zip(rows, y_positions):
        d, lo, hi = r["delta"], r["ci_lo"], r["ci_hi"]
        if np.isfinite(lo) and np.isfinite(hi):
            ci_str = f"Δ={d:+.3f}  [{lo:+.3f}, {hi:+.3f}]"
        else:
            ci_str = f"Δ={d:+.3f}"
        p_str = f"p_boot={_fmt_p(r['p_boot'])} | p_DL={_fmt_p(r['p_form'])}"

        sig = ""
        for p in (r["p_boot"], r["p_form"]):
            if p is not None and np.isfinite(p) and p < 0.05:
                sig = "  ✱"
                break

        ax.text(x_anno, y, f"{ci_str}    {p_str}{sig}",
                va="center", ha="left", fontsize=8.5, family="monospace",
                color="#c0392b" if sig else "#333")

    # Group separators + group labels in left margin
    seen = set()
    for r, y in zip(rows, y_positions):
        if r["group_idx"] not in seen:
            ax.text(x_min, y + 0.7,
                    f"Group {r['group_idx'] + 1} — {r['group_label']}",
                    fontsize=10, fontweight="bold",
                    color=group_colors[r["group_idx"]],
                    ha="left", va="bottom")
            seen.add(r["group_idx"])

    ax.set_ylim(min(y_positions) - 1.0, max(y_positions) + 1.5)
    ax.set_xlabel("Δ AUC (A − B)  with 95% bootstrap CI", fontsize=10.5)
    ax.set_title(
        f"Figure 4D — Focused Pairwise AUC Comparisons under {ctx.cv_caption()}\n"
        f"15 contrasts of interest; {DELONG_BOOT_N:,}-resample bootstrap CI + formal DeLong p",
        fontsize=11, fontweight="bold", pad=12,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.18)
    ax.invert_yaxis()  # group 1 on top, group 3 on bottom
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2B — Panel CpG annotation (3-panel pie chart)
# ─────────────────────────────────────────────────────────────────────────────

def _pie(ax, counts: dict, title: str, colors: list) -> None:
    """Helper: draw a single annotated pie chart on the given axis."""
    labels = [f"{k}\n(n={v})" for k, v in counts.items() if v > 0]
    vals   = [v for v in counts.values() if v > 0]
    if not vals:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title(title, fontsize=11, fontweight="bold")
        return
    _, _, autotexts = ax.pie(
        vals, labels=labels, autopct="%1.1f%%", colors=colors,
        startangle=90, pctdistance=0.78,
        textprops={"fontsize": 9},
        wedgeprops={"linewidth": 0.8, "edgecolor": "white"},
    )
    for at in autotexts:
        at.set_fontsize(8)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)


def _fig2b_panel_annotation(annot_path: str, out_path: str,
                            ctx: _PipelineContext) -> Optional[bool]:
    plt = _mpl()
    df = _read_csv(annot_path)
    if df is None:
        return None

    # Restrict to the final selected panel (panel_cpg.csv) when present, so the
    # figure represents the panel actually used downstream rather than every
    # annotated candidate site.
    panel_cpg = _read_csv(os.path.join(os.path.dirname(annot_path), "panel_cpg.csv"))
    if (panel_cpg is not None and not panel_cpg.empty
            and "feature" in panel_cpg.columns and "feature" in df.columns):
        panel_set = set(panel_cpg["feature"].astype(str))
        df = df[df["feature"].astype(str).isin(panel_set)]
        if df.empty:
            return None

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # CpG-island context
    cpg = (df["cpg_context"].value_counts()
           .reindex(["Island", "Shore", "Shelf", "Open sea"], fill_value=0)
           .to_dict())
    _pie(axes[0], cpg, "CpG Island Context",
         ["#2166ac", "#74add1", "#abd9e9", "#ffffbf"])

    # Genomic region type
    gene = df["region_type"].value_counts().reindex(
        ["promoter", "gene_body", "intergenic"], fill_value=0
    )
    _pie(axes[1],
         {"Promoter":   int(gene.get("promoter",   0)),
          "Gene body":  int(gene.get("gene_body",  0)),
          "Intergenic": int(gene.get("intergenic", 0))},
         "Genomic Region Type",
         ["#d73027", "#fee090", "#91cf60"])

    # Overlap with 88 tissue candidate regions
    in88 = ctx.in_88_count or int(
        (df["in_88_regions"].astype(str).str.lower() == "true").sum()
    )
    _pie(axes[2],
         {"In 88 OC regions": in88, "Genome-wide": len(df) - in88},
         "Tissue Candidate Overlap",
         ["#d7191c", "#bdbdbd"])

    fig.suptitle(f"Figure 2B — Panel CpG Annotation (n={len(df)} sites)",
                 fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 6B — Whole-blood methylation baseline for matched 450K probes
# ─────────────────────────────────────────────────────────────────────────────

def _fig6b_blood_dist(ext_val_path: str, out_path: str,
                      ctx: _PipelineContext) -> Optional[bool]:
    plt = _mpl()
    df = _read_csv(ext_val_path)
    if df is None:
        return None

    # Restrict to the final selected panel (consistent with manuscript reporting).
    panel_cpg = _read_csv(os.path.join(os.path.dirname(ext_val_path), "panel_cpg.csv"))
    if (panel_cpg is not None and not panel_cpg.empty
            and "feature" in panel_cpg.columns and "feature" in df.columns):
        panel_set = set(panel_cpg["feature"].astype(str))
        df = df[df["feature"].astype(str).isin(panel_set)]

    matched = df[df["probe_matched"].astype(str).str.lower() == "true"].copy()
    if matched.empty:
        return None

    matched["beta"] = matched["blood_mean_beta"].astype(float)
    matched["pass"] = matched["blood_pass"].astype(str).str.lower() == "true"
    n_pass = int(matched["pass"].sum())
    n_fail = int((~matched["pass"]).sum())

    fig, ax = plt.subplots(figsize=(8, 4.5))

    rng = np.random.default_rng(DEFAULT_SEED)
    jitter = rng.uniform(-0.35, 0.35, len(matched))
    colors = matched["pass"].map({True: "#4dac26", False: "#d7191c"})
    ax.scatter(matched["beta"], jitter, c=colors, s=70, alpha=0.9,
               edgecolors="white", linewidths=0.5, zorder=3)
    ax.axvline(BLOOD_BETA_THRESHOLD, color="black", lw=1.8, ls="--", zorder=4)

    ax.set_ylim([-1, 1])
    ax.set_xlim([-0.02, 1.02])
    ax.set_yticks([])
    # 656 = the sample columns actually present in GSE40279_series_matrix.txt.gz,
    # and the published cohort size in Hannum et al. 2013. Every other reference to
    # this dataset in the codebase, the docs and the manuscript says 656.
    ax.set_xlabel("Whole-blood mean β (GSE40279, n = 656 samples)", fontsize=12)
    panel_n = ctx.panel_size if ctx.panel_size > 0 else "?"
    ax.set_title(
        f"Figure 6B — Blood Methylation Baseline: "
        f"{len(matched)}/{panel_n} panel sites with matched 450K probe\n"
        f"Blood-pass (β < {BLOOD_BETA_THRESHOLD}): {n_pass} sites   |   "
        f"Blood-fail (β ≥ {BLOOD_BETA_THRESHOLD}): {n_fail} sites "
        f"(potential leukocyte cfDNA confounders)",
        fontsize=10, fontweight="bold",
    )

    from matplotlib.patches import Patch
    import matplotlib.lines as mlines
    ax.legend(handles=[
        Patch(color="#4dac26", label=f"Blood-pass (n={n_pass})"),
        Patch(color="#d7191c", label=f"Blood-fail (n={n_fail})"),
        mlines.Line2D([], [], color="black", lw=1.8, ls="--",
                      label=f"Threshold β={BLOOD_BETA_THRESHOLD}"),
    ], fontsize=9, loc="upper right")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Main orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def generate_imeta_figures(result_dir: str, report_dir: str) -> dict:
    """Generate all auto-computable publication figures + DeLong tables.

    Returns ``{key: abs_path}`` for successfully produced outputs.
    """
    imgs_dir = os.path.join(report_dir, "imgs")
    os.makedirs(imgs_dir, exist_ok=True)

    ctx = _load_context(result_dir)

    clf_pred = os.path.join(result_dir, "classifier_predictions.csv")
    clf_sum  = os.path.join(result_dir, "classifier_summary.csv")
    tv_sum   = os.path.join(result_dir, "temporal_validation_summary.csv")
    tv_pred  = os.path.join(result_dir, "temporal_validation_predictions.csv")
    annot    = os.path.join(result_dir, "panel_annotation.csv")
    ext_val  = os.path.join(result_dir, "external_validation.csv")

    figs = [
        ("fig3_combined", _fig3_combined,
         (clf_pred, clf_sum,              os.path.join(imgs_dir, "pub_fig3_combined.png"),  ctx)),
        ("fig4b_delta", _fig4b_fusion_delta,
         (clf_sum,                        os.path.join(imgs_dir, "pub_fig4b_delta.png"),    ctx)),
        ("fig4c_gain",  _fig4c_spec_gain,
         (clf_sum,                        os.path.join(imgs_dir, "pub_fig4c_gain.png"),     ctx)),
        ("fig5_tv",     _fig5_temporal,
         (tv_sum, tv_pred, clf_sum,       os.path.join(imgs_dir, "pub_fig5_temporal.png"),  ctx)),
        ("fig5b_regime", _fig5b_two_regime,
         (clf_sum, tv_sum,                os.path.join(imgs_dir, "pub_fig5b_regime.png"),   ctx)),
        ("fig2b_panel", _fig2b_panel_annotation,
         (annot,                          os.path.join(imgs_dir, "pub_fig2b_panel.png"),    ctx)),
        ("fig6b_blood", _fig6b_blood_dist,
         (ext_val,                        os.path.join(imgs_dir, "pub_fig6b_blood.png"),    ctx)),
    ]

    out = {"_ctx": ctx}
    for key, fn, args in figs:
        out_path = args[-2]
        try:
            if fn(*args):
                out[key] = out_path
                logger.info("[imeta] %s → %s", key, os.path.basename(out_path))
        except Exception as exc:
            logger.warning("[imeta] %s failed: %s", key, exc, exc_info=False)

    # DeLong tables (all 11 models — covers cfDNA vs CA-125 vs fusion comparisons)
    # plus the paired operating-point tests DeLong cannot provide.
    for key, fn, fname in [
        ("delong_csv",        _delong_bootstrap,  "pub_delong_test.csv"),
        ("delong_formal_csv", _delong_formal,     "pub_delong_formal.csv"),
        ("mcnemar_csv",       _mcnemar_operating, "pub_mcnemar_operating.csv"),
    ]:
        path = os.path.join(result_dir, fname)
        try:
            if fn(clf_pred, path) is not None:
                out[key] = path
                logger.info("[imeta] %s → %s", key, fname)
        except Exception as exc:
            logger.warning("[imeta] %s failed: %s", key, exc)

    # Forest plot of 15 focused pairwise AUC comparisons (reads the CSVs just written)
    forest_path = os.path.join(imgs_dir, "pub_fig_delong_forest.png")
    try:
        if _fig_delong_forest(out.get("delong_csv", ""),
                              out.get("delong_formal_csv", ""),
                              forest_path, ctx):
            out["fig_delong_forest"] = forest_path
            logger.info("[imeta] fig_delong_forest → %s",
                        os.path.basename(forest_path))
    except Exception as exc:
        logger.warning("[imeta] fig_delong_forest failed: %s", exc, exc_info=False)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# HTML section builder — called by report.main.generate_report()
# ─────────────────────────────────────────────────────────────────────────────

def _baseline_vs_ca125_table_html(boot_df, form_df) -> str:
    """Render a focused 10-row table: each non-CA-125 model vs the CA-125 baseline.

    Shows AUC_model, AUC_CA-125, Δ AUC (model − CA-125) with 95% CI, paired
    bootstrap p-value, and formal DeLong p-value. Significant p-values
    (α = 0.05) are bolded in red.
    """
    if boot_df is None and form_df is None:
        return ""

    def _fmt_p(p):
        if p is None or (isinstance(p, float) and not np.isfinite(p)):
            return "<span style='color:#aaa'>NA</span>"
        if p < 0.001:
            return "<strong style='color:#c0392b'>&lt;0.001 ✱</strong>"
        if p < 0.05:
            return f"<strong style='color:#c0392b'>{p:.4f} ✱</strong>"
        return f"<span style='color:#555'>{p:.4f}</span>"

    rows_html = ""
    for model in CFDNA_ONLY + FUSION_MODELS:
        b = _lookup_delong_pair(boot_df, model, CA125_ALONE) if boot_df is not None else None
        f = _lookup_delong_pair(form_df, model, CA125_ALONE) if form_df is not None else None
        if b is None and f is None:
            continue

        # Pull AUCs in (model, CA-125) orientation from the same row
        ref_df = boot_df if boot_df is not None else form_df
        same = ref_df[(ref_df["classifier_A"] == model) & (ref_df["classifier_B"] == CA125_ALONE)]
        flip = ref_df[(ref_df["classifier_A"] == CA125_ALONE) & (ref_df["classifier_B"] == model)]
        if not same.empty:
            auc_model, auc_ca125 = float(same.iloc[0]["auc_A"]), float(same.iloc[0]["auc_B"])
        elif not flip.empty:
            auc_model, auc_ca125 = float(flip.iloc[0]["auc_B"]), float(flip.iloc[0]["auc_A"])
        else:
            continue

        delta, ci_lo, ci_hi, p_boot = b if b is not None else (
            f[0], float("nan"), float("nan"), None
        )
        p_form = f[3] if f is not None else None

        if np.isfinite(ci_lo) and np.isfinite(ci_hi):
            ci_txt = f"[{ci_lo:+.4f}, {ci_hi:+.4f}]"
        else:
            ci_txt = "<span style='color:#aaa'>—</span>"

        group_color = "#4dac26" if "+CA125" in model else "#4393c3"
        model_label = (f"<span style='color:{group_color};font-weight:600'>{model}</span>")

        rows_html += (
            f"<tr>"
            f"<td>{model_label}</td>"
            f"<td>{auc_model:.4f}</td>"
            f"<td>{auc_ca125:.4f}</td>"
            f"<td style='font-weight:600'>{delta:+.4f}</td>"
            f"<td style='font-family:monospace;font-size:.9em'>{ci_txt}</td>"
            f"<td>{_fmt_p(p_boot)}</td>"
            f"<td>{_fmt_p(p_form)}</td>"
            f"</tr>\n"
        )

    if not rows_html:
        return ""

    return f"""
<h4 style="margin-top:14px">Each model vs CA-125 baseline (focused contrasts)</h4>
<table style="font-size:.88em">
  <tr>
    <th>Model</th>
    <th>AUC (model)</th>
    <th>AUC (CA-125)</th>
    <th>Δ AUC<br><span style='font-weight:400;color:#888;font-size:.85em'>model − CA-125</span></th>
    <th>95% bootstrap CI of Δ AUC</th>
    <th><em>p</em> (bootstrap, {DELONG_BOOT_N:,})</th>
    <th><em>p</em> (DeLong)</th>
  </tr>
  {rows_html}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Headline contrasts isolated from the full pairwise matrix: every cfDNA-only and fusion
  model compared head-to-head against CA-125 alone (clinical baseline). Positive Δ AUC
  means the model outperforms CA-125. Bolded red <em>p</em>-values are significant at
  α = 0.05. Bootstrap test: paired-resample two-sided; DeLong test: fastDeLong variance
  estimator (DeLong <em>et al.</em>, 1988). Full 55-pair matrix in the collapsible block below.
</p>"""


def _delong_table_html(df_tbl, mode: str) -> str:
    """Render a DeLong (bootstrap or formal) result table as HTML rows."""
    rows = ""
    for _, r in df_tbl.iterrows():
        sig = ("<strong style='color:#c0392b'>Yes ✱</strong>"
               if str(r.get("significant", "No")) == "Yes" else
               "<span style='color:#aaa'>No</span>")
        if mode == "bootstrap":
            ci = f"{float(r['ci95_lo']):+.4f} to {float(r['ci95_hi']):+.4f}"
            rows += (
                f"<tr><td>{r['classifier_A']} vs {r['classifier_B']}</td>"
                f"<td>{float(r['auc_A']):.4f}</td>"
                f"<td>{float(r['auc_B']):.4f}</td>"
                f"<td>{float(r['delta_auc']):+.4f}</td>"
                f"<td>{ci}</td>"
                f"<td>{float(r['p_bootstrap']):.4f}</td>"
                f"<td>{sig}</td></tr>\n"
            )
        else:  # formal
            z, p = r["z_stat"], r["p_delong"]
            z_fmt = f"{float(z):.3f}" if z != "NA" else "NA"
            p_fmt = f"{float(p):.4f}" if p != "NA" else "NA"
            rows += (
                f"<tr><td>{r['classifier_A']} vs {r['classifier_B']}</td>"
                f"<td>{float(r['auc_A']):.4f}</td>"
                f"<td>{float(r['auc_B']):.4f}</td>"
                f"<td>{float(r['delta_auc']):+.4f}</td>"
                f"<td>{float(r['se_delta']):.5f}</td>"
                f"<td>{z_fmt}</td><td>{p_fmt}</td><td>{sig}</td></tr>\n"
            )
    return rows


def build_imeta_section(result_dir: str, imgs_dir: str) -> str:
    """Generate publication figures + DeLong tables and return an HTML fragment.

    ``imgs_dir`` must be the report's ``imgs/`` subdirectory so that the
    rendered ``<img>`` tags use correct relative paths.
    """
    report_dir = os.path.dirname(imgs_dir)
    fig_paths  = generate_imeta_figures(result_dir, report_dir)
    ctx        = fig_paths.get("_ctx") or _load_context(result_dir)

    def _img_block(key: str, alt: str, fig_label: str, caption: str) -> str:
        path = fig_paths.get(key)
        if not path or not os.path.exists(path):
            return (f'<p style="color:#aaa;font-size:.85em">'
                    f'<em>{fig_label}</em>: not generated (input data absent).</p>')
        rel = os.path.relpath(path, report_dir).replace("\\", "/")
        return (
            f'<figure style="margin:0 0 22px 0">'
            f'<img src="{rel}" alt="{alt}" '
            f'style="max-width:100%;border:1px solid #e0e0e0;border-radius:4px">'
            f'<figcaption style="font-size:.83em;color:#555;margin-top:6px">'
            f'<strong>{fig_label}.</strong> {caption}</figcaption></figure>'
        )

    # ── DeLong forest plot — replaces the two oversized 55-row pairwise tables ─
    forest_cap = (
        f"Focused subset of {len(CFDNA_ONLY) * 3} pairwise AUC contrasts grouped "
        f"by comparison semantics. Group 1 (blue): each cfDNA-only model vs CA-125 alone. "
        f"Group 2 (purple): each fusion model vs CA-125 alone. "
        f"Group 3 (green): each fusion model vs its cfDNA-only counterpart. "
        f"Whiskers = 95% bootstrap CI ({DELONG_BOOT_N:,} resamples); each row is annotated "
        f"with both the bootstrap p-value (p_boot) and the formal DeLong Z-test p-value "
        f"(p_DL). ✱ marks rows where either test is significant at α=0.05. "
        f"Full 55-row pairwise tables remain available as "
        f"<code>pub_delong_test.csv</code> and <code>pub_delong_formal.csv</code>."
    )
    # ── Focused table: each model vs the CA-125 baseline ─────────────────
    boot_df_for_focus = _read_csv(fig_paths.get("delong_csv", ""))
    form_df_for_focus = _read_csv(fig_paths.get("delong_formal_csv", ""))
    baseline_tbl_html = _baseline_vs_ca125_table_html(boot_df_for_focus,
                                                      form_df_for_focus)

    delong_html = (
        f'<h3 style="margin-top:18px">Figure 4D — Focused Pairwise AUC Comparisons '
        f'(Forest Plot)</h3>'
        + _img_block("fig_delong_forest",
                     "DeLong forest plot",
                     "Figure 4D",
                     forest_cap)
        + baseline_tbl_html
    )

    # ── Full pairwise DeLong tables (collapsible; complements the forest plot) ─
    # The forest plot above summarises the 15 contrasts of interest; the two
    # tables below give every pairwise comparison across all 11 models. Wrapped
    # in <details> so the long tables stay collapsed by default.
    boot_tbl_html = ""
    dp = fig_paths.get("delong_csv")
    if dp and os.path.exists(dp):
        try:
            rows = _delong_table_html(_read_csv(dp), "bootstrap")
            boot_tbl_html = f"""
<h4 style="margin-top:14px">Bootstrap Pairwise AUC Comparison ({DELONG_BOOT_N:,} resamples)</h4>
<table style="font-size:.86em">
  <tr><th>Comparison</th><th>AUC A</th><th>AUC B</th><th>Δ AUC</th>
      <th>95% CI of Δ AUC</th><th><em>p</em> (bootstrap)</th><th>Significant</th></tr>
  {rows}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Two-sided bootstrap permutation test with paired resampling ({DELONG_BOOT_N:,} resamples).
  All 11 model comparisons reported (cfDNA-only, CA-125, fusion). ✱ <em>p</em> &lt; 0.05.
  Full table: <code>pub_delong_test.csv</code>
</p>"""
        except Exception as exc:
            logger.warning("[imeta] bootstrap HTML: %s", exc)

    formal_tbl_html = ""
    dfp = fig_paths.get("delong_formal_csv")
    if dfp and os.path.exists(dfp):
        try:
            rows = _delong_table_html(_read_csv(dfp), "formal")
            formal_tbl_html = f"""
<h4 style="margin-top:18px">Formal DeLong Test — fastDeLong Variance Estimator</h4>
<table style="font-size:.86em">
  <tr><th>Comparison</th><th>AUC A</th><th>AUC B</th><th>Δ AUC</th>
      <th>SE(Δ AUC)</th><th><em>Z</em></th>
      <th><em>p</em> (DeLong)</th><th>Significant</th></tr>
  {rows}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Two-sided Z-test using the structural-component variance estimator of
  DeLong <em>et al.</em> (1988); placement-value covariance handles correlated
  AUCs jointly. ✱ <em>p</em> &lt; 0.05. Full table: <code>pub_delong_formal.csv</code>
</p>"""
        except Exception as exc:
            logger.warning("[imeta] formal DeLong HTML: %s", exc)

    if boot_tbl_html or formal_tbl_html:
        delong_formal_html = f"""
<details style="margin-top:10px">
  <summary style="cursor:pointer;font-weight:600;color:#34495e">
    Full pairwise DeLong tables (all 55 comparisons) — click to expand
  </summary>
  {boot_tbl_html}
  {formal_tbl_html}
</details>"""
    else:
        delong_formal_html = ""

    # ── Caption text (built once from context) ─────────────────────────────
    cv_n      = ctx.cv_caption()
    tv_n      = ctx.tv_caption()
    panel_n   = ctx.panel_size or "?"
    ext_n     = ctx.ext_matched or "?"

    def _src(csv_names: str, step_label: str) -> str:
        """Render a small evidence-source annotation appended to a figcaption."""
        return (f'<br><span style="color:#888;font-size:.85em">'
                f'<em>Source:</em> <code>{csv_names}</code> · '
                f'<em>Step:</em> {step_label}</span>')

    fig2b_cap = (
        f"Genomic distribution of the {panel_n}-CpG panel: CpG-island context (left), "
        f"gene region type (centre), and overlap with the 88 ovarian-cancer tissue "
        f"candidate regions (right)."
        + _src("panel_annotation.csv", "Step 18 — annotate")
    )
    fig3_cap = (
        f"Combined binary-classification performance for the eleven evaluated models under {cv_n}. "
        f"Panel A (left): ROC curves on a single axis with line-style distinction — "
        f"cfDNA-only classifiers (LR, SVM, RF, XGB, LGBM) drawn solid, CA-125 alone as dashed black baseline, "
        f"and cfDNA + CA-125 fusion models drawn dotted in the cfDNA base colour. "
        f"Panel B (right): grouped bar chart of the four headline metrics (AUC, Accuracy, Sensitivity, Specificity) "
        f"for the same eleven models; fusion bars use a `/` hatch to remain distinguishable from their cfDNA-only "
        f"counterpart of the same base colour. AUC and metric CIs from {BOOTSTRAP_N:,}-resample bootstrap."
        + _src("classifier_predictions.csv + classifier_summary.csv", "Step 15 — classify")
    )
    fig4b_cap = (
        f"Paired marginal value of CA-125 fusion for each base classifier ({cv_n}). "
        f"Each row contrasts a cfDNA-only model with its corresponding CA-125 fusion "
        f"counterpart (LR vs LR+CA125, …); bars show the fused − base delta with the "
        f"underlying absolute scores annotated. Teal = positive lift, red = regression."
        + _src("classifier_summary.csv", "Step 15 — classify")
    )
    fig4c_cap = (
        f"Per-model gain relative to the CA-125 alone baseline ({cv_n}). "
        f"Panel A: Δ Specificity (the headline benign-rejection metric); "
        f"Panel B: Δ AUC. Fusion models are drawn with a dark border to separate them "
        f"from cfDNA-only models. Bars are sorted by Δ Specificity within each group."
        + _src("classifier_summary.csv", "Step 15 — classify")
    )
    # Append source annotation to the existing DeLong forest-plot caption.
    forest_cap = forest_cap + _src(
        "pub_delong_test.csv, pub_delong_formal.csv",
        "Step 20 — derived from classifier_predictions.csv (Step 15)",
    )
    fig5_cap = (
        f"Left: side-by-side AUC for {ctx.cv_label} vs prospective validation "
        f"({tv_n}). Right: per-sample predicted malignancy probability for the "
        f"prospective best model, plotted against sample collection date. "
        f"✕ = misclassified sample."
        + _src("temporal_validation_summary.csv, temporal_validation_predictions.csv",
               "Step 17 — temporal_validate")
    )
    fig5b_cap = (
        f"Cross-regime consistency of AUC (left) and Specificity (right) across the "
        f"eleven evaluated models. Solid bars = {ctx.cv_label} ({cv_n}); hatched bars = "
        f"temporal prospective validation ({tv_n}). 95% bootstrap CI. Bars sharing the "
        f"same colour represent the same classifier in the two regimes."
        + _src("classifier_summary.csv + temporal_validation_summary.csv",
               "Steps 15 — classify  +  17 — temporal_validate")
    )
    fig6b_cap = (
        f"Whole-blood mean β-value (GSE40279) for the {ext_n} panel CpG sites with "
        f"matched Illumina 450K probes. Sites with blood mean β ≥ {BLOOD_BETA_THRESHOLD} "
        f"({ctx.ext_fail} sites, red) are potential leukocyte-DNA confounders in plasma "
        f"cfDNA liquid biopsy and are candidates for panel refinement."
        + _src("external_validation.csv", "Step 19 — external_validate")
    )

    # Rebuild the DeLong forest-plot block with the updated, source-annotated caption.
    delong_html = (
        f'<h3 style="margin-top:18px">Figure 4D — Focused Pairwise AUC Comparisons '
        f'(Forest Plot)</h3>'
        + _img_block("fig_delong_forest",
                     "DeLong forest plot",
                     "Figure 4D",
                     forest_cap)
        + baseline_tbl_html
    )

    # ── Section index: figure → source CSV → producing pipeline step ──────
    section_index_html = """
<div style="background:#f8f9fa;border-left:3px solid #6c757d;padding:10px 14px;margin:10px 0 18px 0;font-size:.86em">
  <strong>Section index</strong> — figures appear below in pipeline-step order:
  <table style="width:100%;margin-top:6px;font-size:.95em">
    <tr><th style="text-align:left">Figure</th><th style="text-align:left">Source CSV</th><th style="text-align:left">Pipeline step</th></tr>
    <tr><td>Figure 2B</td><td><code>panel_annotation.csv</code></td><td>Step 18 — annotate</td></tr>
    <tr><td>Figure 3</td><td><code>classifier_predictions.csv</code> + <code>classifier_summary.csv</code></td><td>Step 15 — classify</td></tr>
    <tr><td>Figure 4B</td><td><code>classifier_summary.csv</code></td><td>Step 15 — classify</td></tr>
    <tr><td>Figure 4C</td><td><code>classifier_summary.csv</code></td><td>Step 15 — classify</td></tr>
    <tr><td>Figure 4D + full DeLong tables</td><td><code>pub_delong_test.csv</code>, <code>pub_delong_formal.csv</code></td><td>Step 20 — derived from classifier_predictions.csv (Step 15)</td></tr>
    <tr><td>Figure 5</td><td><code>temporal_validation_summary.csv</code>, <code>temporal_validation_predictions.csv</code></td><td>Step 17 — temporal_validate</td></tr>
    <tr><td>Figure 5B</td><td><code>classifier_summary.csv</code> + <code>temporal_validation_summary.csv</code></td><td>Steps 15 + 17</td></tr>
    <tr><td>Figure 6B</td><td><code>external_validation.csv</code></td><td>Step 19 — external_validate</td></tr>
  </table>
</div>"""

    return f"""
<!-- Publication figures — auto-generated by Step 20 -->
<div class="section">
<h2>Statistical &amp; Publication Figures (manuscript prep)</h2>
<p style="color:#666;font-size:.85em">
  Auto-generated from result CSV files during Step 20 (report). {FIG_DPI} DPI PNG files:
  <code>report/imgs/pub_*.png</code>.
  See <code>docs/GPB/01_transfer_checklist.md</code> for the full submission checklist.
</p>
{section_index_html}

<h3>Figure 2B — Panel CpG Site Annotation</h3>
{_img_block("fig2b_panel", "Panel Annotation", "Figure 2B", fig2b_cap)}

<h3 style="margin-top:22px">Figure 3 — Binary Classification Performance (ROC + headline metrics, eleven-model comparison)</h3>
{_img_block("fig3_combined", "Binary Classification Performance", "Figure 3", fig3_cap)}

<h3 style="margin-top:22px">Figure 4B — Marginal Value of CA-125 Fusion (paired Δ per base classifier)</h3>
{_img_block("fig4b_delta", "Fusion vs cfDNA-only paired delta", "Figure 4B", fig4b_cap)}

<h3 style="margin-top:22px">Figure 4C — Per-model Gain over CA-125 Alone</h3>
{_img_block("fig4c_gain", "Gain over CA-125 alone", "Figure 4C", fig4c_cap)}

{delong_html}
{delong_formal_html}

<h3 style="margin-top:22px">Figure 5 — Temporal Prospective Validation</h3>
{_img_block("fig5_tv", "Temporal Validation", "Figure 5", fig5_cap)}

<h3 style="margin-top:22px">Figure 5B — Cross-regime Consistency (CV vs Temporal)</h3>
{_img_block("fig5b_regime", "CV vs Temporal cross-regime", "Figure 5B", fig5b_cap)}

<h3 style="margin-top:22px">Figure 6B — Whole-blood Methylation Baseline</h3>
{_img_block("fig6b_blood", "Blood Methylation Baseline", "Figure 6B", fig6b_cap)}
</div>"""
