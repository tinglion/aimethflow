"""
imeta_utils.py — auto-generate iMeta publication figures and statistical tables.

Called during Step 20 (report). All functions are optional: missing input files
cause a log warning and the figure is silently skipped.
"""

import os
import logging
import numpy as np

logger = logging.getLogger(__name__)

# Colorblind-safe palette per classifier
_CLF_COLORS = {"LR": "#4393c3", "SVM": "#d6604d", "RF": "#4dac26", "XGB": "#8073ac", "LGBM": "#e08214"}
_CLF_ORDER  = ["LR", "SVM", "RF", "XGB", "LGBM"]

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _read_csv(path):
    try:
        import pandas as pd
        if not os.path.exists(path):
            return None
        return pd.read_csv(path, encoding="utf-8")
    except Exception as exc:
        logger.warning("[imeta] read %s: %s", path, exc)
        return None


def _savefig(fig, path):
    import matplotlib.pyplot as plt
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
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
# Figure 3: Multi-classifier ROC curves
# ─────────────────────────────────────────────────────────────────────────────

def _fig3_roc(clf_pred_path, out_path, n_boot=1000):
    from sklearn.metrics import roc_curve, auc as sk_auc
    plt = _mpl()

    df = _read_csv(clf_pred_path)
    if df is None:
        return None

    fig, ax = plt.subplots(figsize=(7, 6))
    rng = np.random.default_rng(42)
    any_plotted = False

    for clf in _CLF_ORDER:
        sub = df[df["classifier"] == clf]
        if sub.empty:
            continue
        y_true = (sub["true_label"] == "malignant").astype(int).values
        y_prob = sub["prob_malignant"].astype(float).values

        fpr, tpr, _ = roc_curve(y_true, y_prob)
        auc_val = sk_auc(fpr, tpr)

        boot_aucs = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(y_true), len(y_true))
            if len(np.unique(y_true[idx])) < 2:
                continue
            f, t, _ = roc_curve(y_true[idx], y_prob[idx])
            boot_aucs.append(sk_auc(f, t))
        ci_lo, ci_hi = (np.percentile(boot_aucs, [2.5, 97.5]) if boot_aucs
                        else (float("nan"), float("nan")))

        ax.plot(fpr, tpr, color=_CLF_COLORS[clf], lw=2,
                label=f"{clf}  AUC={auc_val:.3f} (95%CI {ci_lo:.3f}–{ci_hi:.3f})")
        any_plotted = True

    if not any_plotted:
        plt.close(fig)
        return None

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.4, label="Random (AUC=0.500)")
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1.02])
    ax.set_xlabel("1 − Specificity (False Positive Rate)", fontsize=12)
    ax.set_ylabel("Sensitivity (True Positive Rate)", fontsize=12)
    ax.set_title("Binary Classification ROC — LOO-CV", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# DeLong-equivalent: bootstrap pairwise AUC comparison table
# ─────────────────────────────────────────────────────────────────────────────

def _delong_bootstrap(clf_pred_path, out_csv_path, n_boot=2000):
    import pandas as pd
    from sklearn.metrics import roc_auc_score

    df = _read_csv(clf_pred_path)
    if df is None:
        return None

    data = {}
    for clf in _CLF_ORDER:
        sub = df[df["classifier"] == clf]
        if sub.empty:
            continue
        data[clf] = {
            "y_true": (sub["true_label"] == "malignant").astype(int).values,
            "y_prob": sub["prob_malignant"].astype(float).values,
        }

    clfs = list(data.keys())
    if len(clfs) < 2:
        return None

    n = len(data[clfs[0]]["y_true"])
    rng = np.random.default_rng(42)
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
            p    = min(np.mean(diffs <= 0), np.mean(diffs >= 0)) * 2 if len(diffs) else 1.0
            ci_lo, ci_hi = (np.percentile(diffs, [2.5, 97.5]) if len(diffs) else (float("nan"), float("nan")))
            rows.append({
                "classifier_A": c1, "classifier_B": c2,
                "auc_A": round(auc1, 4), "auc_B": round(auc2, 4),
                "delta_auc": round(obs, 4),
                "ci95_lo": round(ci_lo, 4), "ci95_hi": round(ci_hi, 4),
                "p_bootstrap": round(p, 4),
                "significant": "Yes" if p < 0.05 else "No",
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv_path, index=False, encoding="utf-8")
    return df_out


# ─────────────────────────────────────────────────────────────────────────────
# DeLong 1988 — formal variance estimation (fastDeLong algorithm)
# ─────────────────────────────────────────────────────────────────────────────

def _fastdelong_cov(y_true, score_matrix):
    """
    fastDeLong covariance matrix for K classifiers.

    y_true      : 1-D array of 0/1 labels, length n
    score_matrix: (K, n) array — row k is classifier k's probability scores

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
        # V10[k, i]: fraction of negatives strictly below ps[i] + 0.5 × ties
        cmp10 = ns[np.newaxis, :] < ps[:, np.newaxis]    # (n_pos, n_neg)
        tie10 = ns[np.newaxis, :] == ps[:, np.newaxis]
        V10[k] = cmp10.mean(axis=1) + 0.5 * tie10.mean(axis=1)
        # V01[k, j]: fraction of positives strictly above ns[j] + 0.5 × ties
        cmp01 = ps[np.newaxis, :] > ns[:, np.newaxis]    # (n_neg, n_pos)
        tie01 = ps[np.newaxis, :] == ns[:, np.newaxis]
        V01[k] = cmp01.mean(axis=1) + 0.5 * tie01.mean(axis=1)

    aucs = V10.mean(axis=1)   # AUC_k = mean placement value

    # Structural-component covariance (DeLong 1988, eq. 7)
    # Cov(AUC_a, AUC_b) = Cov_sample(V10_a, V10_b)/n_pos + Cov_sample(V01_a, V01_b)/n_neg
    S10 = np.cov(V10, ddof=1) if n_pos > 1 else np.zeros((K, K))  # (K, K)
    S01 = np.cov(V01, ddof=1) if n_neg > 1 else np.zeros((K, K))  # (K, K)
    # np.cov returns scalar when K==1; force (1,1) shape
    if K == 1:
        S10 = np.atleast_2d(S10)
        S01 = np.atleast_2d(S01)
    cov = S10 / n_pos + S01 / n_neg

    return aucs, cov


def _delong_formal(clf_pred_path, out_csv_path):
    """
    Formal pairwise AUC comparison using the fastDeLong variance estimator.
    Produces a normal Z-test (two-sided) for each classifier pair.
    Saves result/imeta_delong_formal.csv.
    """
    import pandas as pd
    from scipy import stats

    df = _read_csv(clf_pred_path)
    if df is None:
        return None

    # Build (K, n) score matrix — requires identical sample order per classifier
    ref_clf = None
    ref_ids = None
    data = {}
    for clf in _CLF_ORDER:
        sub = df[df["classifier"] == clf].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("sample_id").reset_index(drop=True)
        if ref_ids is None:
            ref_ids = sub["sample_id"].values
            ref_clf = clf
        elif not np.array_equal(sub["sample_id"].values, ref_ids):
            logger.warning(
                "[delong] sample_id mismatch between %s and %s — skipping DeLong", ref_clf, clf
            )
            return None
        data[clf] = {
            "y_true": (sub["true_label"] == "malignant").astype(int).values,
            "y_prob": sub["prob_malignant"].astype(float).values,
        }

    clfs = list(data.keys())
    if len(clfs) < 2:
        return None

    y_true = data[clfs[0]]["y_true"]
    score_matrix = np.vstack([data[c]["y_prob"] for c in clfs])   # (K, n)

    aucs, cov = _fastdelong_cov(y_true, score_matrix)

    rows = []
    for i, ca in enumerate(clfs):
        for j in range(i + 1, len(clfs)):
            cb  = clfs[j]
            d   = aucs[i] - aucs[j]
            var = cov[i, i] + cov[j, j] - 2 * cov[i, j]
            se  = np.sqrt(max(var, 0.0))
            z   = d / se if se > 0 else float("nan")
            p   = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else float("nan")
            rows.append({
                "classifier_A": ca,
                "classifier_B": cb,
                "auc_A":        round(float(aucs[i]), 4),
                "auc_B":        round(float(aucs[j]), 4),
                "delta_auc":    round(float(d),  4),
                "se_delta":     round(float(se), 6),
                "z_stat":       round(float(z),  4) if np.isfinite(z) else "NA",
                "p_delong":     round(float(p),  4) if np.isfinite(p) else "NA",
                "significant":  "Yes" if (np.isfinite(p) and p < 0.05) else "No",
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv_path, index=False, encoding="utf-8")
    logger.info(
        "[delong] %d pairs tested; formal DeLong saved → %s",
        len(rows), os.path.basename(out_csv_path),
    )
    return df_out


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4B: Classifier metric comparison bar chart
# ─────────────────────────────────────────────────────────────────────────────

def _fig4_clf_comparison(clf_sum_path, out_path):
    plt = _mpl()

    df = _read_csv(clf_sum_path)
    if df is None:
        return None

    metrics = [
        ("auc",         "auc_ci95_lo",   "auc_ci95_hi",   "AUC"),
        ("accuracy",    "acc_ci95_lo",    "acc_ci95_hi",   "Accuracy"),
        ("sensitivity", "sens_ci95_lo",   "sens_ci95_hi",  "Sensitivity"),
        ("specificity", "spec_ci95_lo",   "spec_ci95_hi",  "Specificity"),
    ]
    clfs = [c for c in _CLF_ORDER if c in df["classifier"].values]
    if not clfs:
        return None

    x      = np.arange(len(metrics))
    width  = 0.14
    n      = len(clfs)
    offsets = np.linspace(-(n - 1) / 2, (n - 1) / 2, n) * width

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, clf in enumerate(clfs):
        row = df[df["classifier"] == clf].iloc[0]
        vals    = [float(row.get(m, 0))  for m, *_ in metrics]
        cis_lo  = [float(row.get(lo, v)) for (m, lo, hi, _), v in zip(metrics, vals)]
        cis_hi  = [float(row.get(hi, v)) for (m, lo, hi, _), v in zip(metrics, vals)]
        yerr = [[v - lo for v, lo in zip(vals, cis_lo)],
                [hi - v for v, hi in zip(vals, cis_hi)]]
        bars = ax.bar(x + offsets[i], vals, width, yerr=yerr,
                      color=_CLF_COLORS[clf], alpha=0.85, capsize=3.5, label=clf,
                      error_kw={"elinewidth": 1.2, "ecolor": "#555"})
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(yerr[1]) + 0.01,
                    f"{bar.get_height():.2f}",
                    ha="center", va="bottom", fontsize=6.5, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels([m[-1] for m in metrics], fontsize=12)
    ax.set_ylim([0, 1.15])
    ax.axhline(0.5, color="#bbb", lw=1, ls="--", alpha=0.7)
    ax.set_ylabel("Score", fontsize=12)
    ax.set_title("Classifier Performance — LOO-CV (n=116, 95% CI bootstrap)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper right", fontsize=9, ncol=n)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5: Temporal validation — AUC comparison + prediction scatter
# ─────────────────────────────────────────────────────────────────────────────

def _fig5_temporal(tv_sum_path, tv_pred_path, clf_sum_path, out_path):
    import pandas as pd
    plt = _mpl()

    tv_sum  = _read_csv(tv_sum_path)
    tv_pred = _read_csv(tv_pred_path)
    clf_sum = _read_csv(clf_sum_path)
    if tv_sum is None or tv_pred is None:
        return None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # ── Panel A: LOO-CV vs prospective AUC bar chart ──────────────────────
    clfs = [c for c in _CLF_ORDER if c in tv_sum["classifier"].values]
    x = np.arange(len(clfs))
    w = 0.32

    loo_aucs = []
    if clf_sum is not None:
        for c in clfs:
            row = clf_sum[clf_sum["classifier"] == c]
            loo_aucs.append(float(row["auc"].iloc[0]) if not row.empty else 0.0)
    tv_aucs = [float(tv_sum[tv_sum["classifier"] == c]["auc"].iloc[0]) for c in clfs]

    b1 = ax1.bar(x - w / 2, loo_aucs, w, label="LOO-CV (n=116)", color="#4393c3", alpha=0.85)
    b2 = ax1.bar(x + w / 2, tv_aucs,  w, label="Prospective (n=20)", color="#d6604d", alpha=0.85)
    ax1.axhline(0.5, color="#bbb", lw=1, ls="--", alpha=0.7)
    ax1.set_xticks(x); ax1.set_xticklabels(clfs, fontsize=10)
    ax1.set_ylim([0, 1.08]); ax1.set_ylabel("AUC", fontsize=11)
    ax1.set_title("LOO-CV vs Prospective Validation AUC", fontsize=11, fontweight="bold")
    ax1.legend(fontsize=9); ax1.grid(axis="y", alpha=0.2)
    for bar in list(b1) + list(b2):
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.012,
                 f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=7.5)

    # ── Panel B: prediction scatter (LGBM, best prospective performer) ────
    best_clf = tv_sum.loc[tv_sum["auc"].astype(float).idxmax(), "classifier"]
    sub = tv_pred[tv_pred["classifier"] == best_clf].copy()
    if sub.empty:
        sub = tv_pred[tv_pred["classifier"] == tv_pred["classifier"].iloc[0]].copy()

    sub["date"] = pd.to_datetime(sub["collection_date"].astype(str), format="%Y%m%d", errors="coerce")
    sub["prob"] = sub["prob_malignant"].astype(float)
    sub["correct"] = sub["correct"].astype(str).str.lower().isin(["true", "1"])

    for label, color, marker in [("malignant", "#d6604d", "o"), ("benign", "#4393c3", "s")]:
        grp = sub[sub["true_label"] == label]
        ax2.scatter(grp["date"], grp["prob"],
                    c=color, marker=marker, s=80, alpha=0.85, zorder=3,
                    label=f"{label.capitalize()} (n={len(grp)})")

    wrong = sub[~sub["correct"]]
    if not wrong.empty:
        ax2.scatter(wrong["date"], wrong["prob"],
                    marker="x", s=130, c="black", linewidths=2, zorder=5,
                    label="Misclassified")

    ax2.axhline(0.5, color="#888", lw=1.2, ls="--", alpha=0.7, label="Decision boundary (0.5)")
    ax2.set_ylim([-0.07, 1.07])
    ax2.set_xlabel("Sample collection date", fontsize=11)
    ax2.set_ylabel(f"P(Malignant) — {best_clf}", fontsize=11)
    ax2.set_title(f"Prospective Validation Predictions ({best_clf}, n=20)\n"
                  f"Note: 4 benign / 16 malignant — specificity estimate provisional",
                  fontsize=10, fontweight="bold")
    fig.autofmt_xdate(rotation=30)
    ax2.legend(fontsize=8.5, loc="center left"); ax2.grid(alpha=0.2)

    fig.suptitle("Figure 5 — Temporal Prospective Validation", fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2B: Panel CpG annotation pie charts
# ─────────────────────────────────────────────────────────────────────────────

def _fig2b_panel_annotation(annot_path, out_path):
    plt = _mpl()

    df = _read_csv(annot_path)
    if df is None:
        return None

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    def _pie(ax, counts, title, colors):
        labels = [f"{k}\n(n={v})" for k, v in counts.items() if v > 0]
        vals   = [v for v in counts.values() if v > 0]
        if not vals:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title, fontsize=11, fontweight="bold")
            return
        wedges, texts, autotexts = ax.pie(
            vals, labels=labels, autopct="%1.1f%%", colors=colors,
            startangle=90, pctdistance=0.78,
            textprops={"fontsize": 9}, wedgeprops={"linewidth": 0.8, "edgecolor": "white"},
        )
        for at in autotexts:
            at.set_fontsize(8)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=10)

    cpg = df["cpg_context"].value_counts().reindex(["Island", "Shore", "Shelf", "Open sea"], fill_value=0).to_dict()
    _pie(axes[0], cpg, "CpG Island Context",
         ["#2166ac", "#74add1", "#abd9e9", "#ffffbf"])

    gene = df["region_type"].value_counts().reindex(["promoter", "gene_body", "intergenic"], fill_value=0)
    _pie(axes[1],
         {"Promoter": int(gene.get("promoter", 0)),
          "Gene body": int(gene.get("gene_body", 0)),
          "Intergenic": int(gene.get("intergenic", 0))},
         "Genomic Region Type",
         ["#d73027", "#fee090", "#91cf60"])

    in88 = int((df["in_88_regions"].astype(str).str.lower() == "true").sum())
    _pie(axes[2],
         {"In 88 OC regions": in88, "Genome-wide": len(df) - in88},
         "Tissue Candidate Overlap",
         ["#d7191c", "#bdbdbd"])

    fig.suptitle(f"Figure 2B — Panel CpG Annotation (n={len(df)} sites)",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Figure 6B: Whole-blood methylation baseline (cfDNA specificity)
# ─────────────────────────────────────────────────────────────────────────────

def _fig6b_blood_dist(ext_val_path, out_path):
    plt = _mpl()

    df = _read_csv(ext_val_path)
    if df is None:
        return None

    matched = df[df["probe_matched"].astype(str).str.lower() == "true"].copy()
    if matched.empty:
        return None

    matched["beta"] = matched["blood_mean_beta"].astype(float)
    matched["pass"] = matched["blood_pass"].astype(str).str.lower() == "true"
    threshold = 0.10
    n_pass = int(matched["pass"].sum())
    n_fail = int((~matched["pass"]).sum())

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Dot strip chart with jitter
    rng = np.random.default_rng(42)
    jitter = rng.uniform(-0.35, 0.35, len(matched))
    colors = matched["pass"].map({True: "#4dac26", False: "#d7191c"})
    ax.scatter(matched["beta"], jitter, c=colors, s=70, alpha=0.9,
               edgecolors="white", linewidths=0.5, zorder=3)

    ax.axvline(threshold, color="black", lw=1.8, ls="--", zorder=4,
               label=f"Threshold β = {threshold}")
    ax.set_ylim([-1, 1])
    ax.set_xlim([-0.02, 1.02])
    ax.set_yticks([])
    ax.set_xlabel("Whole-blood mean β (GSE40279, n=729 samples)", fontsize=12)
    ax.set_title(
        f"Figure 6B — Blood Methylation Baseline: {len(matched)}/328 panel sites with matched 450K probe\n"
        f"Blood-pass (β < {threshold}): {n_pass} sites   |   Blood-fail (β ≥ {threshold}): {n_fail} sites "
        f"(potential leukocyte cfDNA confounders)",
        fontsize=10, fontweight="bold",
    )

    from matplotlib.patches import Patch
    import matplotlib.lines as mlines
    ax.legend(handles=[
        Patch(color="#4dac26", label=f"Blood-pass (n={n_pass})"),
        Patch(color="#d7191c", label=f"Blood-fail (n={n_fail})"),
        mlines.Line2D([], [], color="black", lw=1.8, ls="--", label=f"Threshold β={threshold}"),
    ], fontsize=9, loc="upper right")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    _savefig(fig, out_path)
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Main orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def generate_imeta_figures(result_dir, report_dir):
    """
    Generate all auto-computable iMeta figures + DeLong table.
    Returns dict {key: abs_path} for successfully produced outputs.
    """
    imgs_dir = os.path.join(report_dir, "imgs")
    os.makedirs(imgs_dir, exist_ok=True)

    clf_pred = os.path.join(result_dir, "classifier_predictions.csv")
    clf_sum  = os.path.join(result_dir, "classifier_summary.csv")
    tv_sum   = os.path.join(result_dir, "temporal_validation_summary.csv")
    tv_pred  = os.path.join(result_dir, "temporal_validation_predictions.csv")
    annot    = os.path.join(result_dir, "panel_annotation.csv")
    ext_val  = os.path.join(result_dir, "external_validation.csv")

    tasks = [
        ("fig3_roc",    _fig3_roc,              clf_pred,                                     os.path.join(imgs_dir, "imeta_fig3_roc.png")),
        ("fig4_clf",    _fig4_clf_comparison,   clf_sum,                                      os.path.join(imgs_dir, "imeta_fig4_clf.png")),
        ("fig5_tv",     _fig5_temporal,         tv_sum, tv_pred, clf_sum,                     os.path.join(imgs_dir, "imeta_fig5_temporal.png")),
        ("fig2b_panel", _fig2b_panel_annotation,annot,                                        os.path.join(imgs_dir, "imeta_fig2b_panel.png")),
        ("fig6b_blood", _fig6b_blood_dist,       ext_val,                                     os.path.join(imgs_dir, "imeta_fig6b_blood.png")),
    ]

    out = {}
    for key, fn, *args in tasks:
        out_path = args[-1]
        in_args  = args[:-1]
        try:
            result = fn(*in_args, out_path)
            if result is not None:
                out[key] = out_path
                logger.info("[imeta] %s → %s", key, os.path.basename(out_path))
        except Exception as exc:
            logger.warning("[imeta] %s failed: %s", key, exc)

    # Bootstrap pairwise AUC table
    delong_path = os.path.join(result_dir, "imeta_delong_test.csv")
    try:
        df_dl = _delong_bootstrap(clf_pred, delong_path)
        if df_dl is not None:
            out["delong_csv"] = delong_path
            logger.info("[imeta] bootstrap AUC table → %s", os.path.basename(delong_path))
    except Exception as exc:
        logger.warning("[imeta] bootstrap AUC table failed: %s", exc)

    # Formal DeLong test (fastDeLong variance estimator)
    delong_formal_path = os.path.join(result_dir, "imeta_delong_formal.csv")
    try:
        df_dl_formal = _delong_formal(clf_pred, delong_formal_path)
        if df_dl_formal is not None:
            out["delong_formal_csv"] = delong_formal_path
    except Exception as exc:
        logger.warning("[imeta] formal DeLong test failed: %s", exc)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# HTML section builder — called by report_utils.generate_report()
# ─────────────────────────────────────────────────────────────────────────────

def build_imeta_section(result_dir, imgs_dir):
    """
    Generate iMeta figures and return an HTML fragment for the Step 20 report.
    imgs_dir must be the report's imgs/ subdirectory (for correct relative paths).
    """
    report_dir = os.path.dirname(imgs_dir)
    fig_paths  = generate_imeta_figures(result_dir, report_dir)

    def _img_block(key, alt, fig_label, caption):
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
            f'<strong>{fig_label}.</strong> {caption}'
            f'</figcaption></figure>'
        )

    import pandas as pd

    def _auc_table_html(df_tbl, mode):
        """Build HTML for either bootstrap or formal DeLong table."""
        rows = ""
        for _, r in df_tbl.iterrows():
            sig = ("<strong style='color:#c0392b'>Yes ✱</strong>"
                   if str(r.get("significant", "No")) == "Yes" else
                   "<span style='color:#aaa'>No</span>")
            if mode == "bootstrap":
                p_val = f"{float(r['p_bootstrap']):.4f}"
                ci    = f"{float(r['ci95_lo']):+.4f} to {float(r['ci95_hi']):+.4f}"
                rows += (
                    f"<tr><td>{r['classifier_A']} vs {r['classifier_B']}</td>"
                    f"<td>{float(r['auc_A']):.4f}</td><td>{float(r['auc_B']):.4f}</td>"
                    f"<td>{float(r['delta_auc']):+.4f}</td>"
                    f"<td>{ci}</td>"
                    f"<td>{p_val}</td><td>{sig}</td></tr>\n"
                )
            else:
                z_s  = r['z_stat']
                p_s  = r['p_delong']
                z_fmt = f"{float(z_s):.3f}" if z_s != "NA" else "NA"
                p_fmt = f"{float(p_s):.4f}" if p_s != "NA" else "NA"
                rows += (
                    f"<tr><td>{r['classifier_A']} vs {r['classifier_B']}</td>"
                    f"<td>{float(r['auc_A']):.4f}</td><td>{float(r['auc_B']):.4f}</td>"
                    f"<td>{float(r['delta_auc']):+.4f}</td>"
                    f"<td>{float(r['se_delta']):.5f}</td>"
                    f"<td>{z_fmt}</td><td>{p_fmt}</td><td>{sig}</td></tr>\n"
                )
        return rows

    # Bootstrap table
    delong_html = ""
    dp = fig_paths.get("delong_csv")
    if dp and os.path.exists(dp):
        try:
            dl = pd.read_csv(dp, encoding="utf-8")
            rows = _auc_table_html(dl, "bootstrap")
            delong_html = f"""
<h3 style="margin-top:18px">Bootstrap Pairwise AUC Comparison (2,000 resamples)</h3>
<table style="font-size:.88em">
  <tr><th>Comparison</th><th>AUC A</th><th>AUC B</th><th>Δ AUC</th>
      <th>95% CI of Δ AUC</th><th><em>p</em> (bootstrap)</th><th>Significant</th></tr>
  {rows}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Two-sided bootstrap permutation test; 2,000 resamples with replacement.
  ✱ <em>p</em> &lt; 0.05. Full table: <code>imeta_delong_test.csv</code>
</p>"""
        except Exception as exc:
            logger.warning("[imeta] bootstrap HTML: %s", exc)

    # Formal DeLong table
    delong_formal_html = ""
    dfp = fig_paths.get("delong_formal_csv")
    if dfp and os.path.exists(dfp):
        try:
            dl_f = pd.read_csv(dfp, encoding="utf-8")
            rows_f = _auc_table_html(dl_f, "formal")
            delong_formal_html = f"""
<h3 style="margin-top:18px">Formal DeLong Test — fastDeLong Variance Estimator</h3>
<table style="font-size:.88em">
  <tr><th>Comparison</th><th>AUC A</th><th>AUC B</th><th>Δ AUC</th>
      <th>SE(Δ AUC)</th><th><em>Z</em></th>
      <th><em>p</em> (DeLong)</th><th>Significant</th></tr>
  {rows_f}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Two-sided Z-test using the structural-component variance estimator of DeLong <em>et al.</em> (1988).
  SE estimated from placement-value covariance; correlated AUCs handled jointly.
  ✱ <em>p</em> &lt; 0.05. Full table: <code>imeta_delong_formal.csv</code>
</p>"""
        except Exception as exc:
            logger.warning("[imeta] formal DeLong HTML: %s", exc)

    return f"""
<!-- iMeta Publication Figures — auto-generated by Step 20 -->
<div class="section">
<h2>Statistical &amp; Publication Figures (iMeta prep)</h2>
<p style="color:#666;font-size:.85em">
  Auto-generated from result CSV files during Step 20 (report). 300 DPI PNG files:
  <code>report/imgs/imeta_*.png</code>.
  See <code>docs/iMeta_todo.md</code> for the full submission checklist.
</p>

<h3>Figure 3 — Binary Classification ROC Curves</h3>
{_img_block("fig3_roc", "ROC Curves", "Figure 3",
    "ROC curves for LR, SVM, RF, XGB, and LGBM evaluated via leave-one-out cross-validation "
    "(n=116: 56 benign, 60 malignant). AUC and 95% CIs from 1,000-resample bootstrap.")}
{delong_html}
{delong_formal_html}

<h3 style="margin-top:22px">Figure 4B — Classifier Performance Comparison</h3>
{_img_block("fig4_clf", "Classifier Comparison", "Figure 4B",
    "AUC, accuracy, sensitivity, and specificity for each classifier (LOO-CV, n=116). "
    "Error bars = 95% bootstrap CI (1,000 resamples).")}

<h3 style="margin-top:22px">Figure 5 — Temporal Prospective Validation</h3>
{_img_block("fig5_tv", "Temporal Validation", "Figure 5",
    "Left: side-by-side AUC for LOO-CV (n=116) vs prospective validation (n=20; cutoff 2024-07-01). "
    "Right: per-sample predicted malignancy probability for the best prospective classifier, "
    "plotted against sample collection date. ✕ = misclassified sample. "
    "<strong>Note:</strong> validation set contains only 4 benign samples — specificity is provisional.")}

<h3 style="margin-top:22px">Figure 2B — Panel CpG Site Annotation</h3>
{_img_block("fig2b_panel", "Panel Annotation", "Figure 2B",
    "Genomic distribution of the 328-CpG panel: "
    "CpG island context (left), gene region type (centre), "
    "and overlap with the 88 ovarian cancer tissue candidate regions (right).")}

<h3 style="margin-top:22px">Figure 6B — Whole-blood Methylation Baseline</h3>
{_img_block("fig6b_blood", "Blood Methylation Baseline", "Figure 6B",
    "Whole-blood mean β-value (GSE40279) for the 25 panel CpG sites with matched Illumina 450K probes. "
    "Sites with blood mean β ≥ 0.10 (n=17, red) are potential leukocyte-DNA confounders "
    "in plasma cfDNA liquid biopsy and are candidates for panel refinement.")}
</div>"""
