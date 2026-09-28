"""build_figure4_combined.py — combine Figure 4A (5 cfDNA-only classifiers)
and Figure 4B (5 CA-125 fusion models) into one two-panel Figure 4, matching
the manuscript's own framing ("(A) ... (B) the same four arms for the five
cfDNA + CA-125 fusion models"). Reuses the exact arm order, colors and
annotation logic as leakage_demo.py::_save_plot so the combined figure is
pixel-identical in style to the two standalone images it replaces.

Run once from the project root:
    python src/scripts/build_figure4_combined.py
"""

import os

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SOURCE_CSV = os.path.join(PROJECT_ROOT, "result", "20260406", "leakage_demo.csv")
OUT_PATH = os.path.join(PROJECT_ROOT, "docs", "GPB", "resources", "Figure4_combined.png")

ARMS = ["in_fold", "naive_selection", "naive_full", "panel_then_in_fold"]
REFERENCE_ARM = "in_fold"
CA125_SUFFIX = "+CA125"
ARM_COLORS = {"in_fold": "#4dac26", "naive_selection": "#e08214",
              "naive_full": "#d6604d", "panel_then_in_fold": "#5e8fc7"}


def _panel(ax, df, clfs, title):
    arms = [a for a in ARMS if a in set(df["arm"])]
    x = np.arange(len(clfs))
    width = 0.8 / max(len(arms), 1)

    for i, arm in enumerate(arms):
        sub = df[df["arm"] == arm].set_index("classifier")
        vals = np.array([sub.loc[c, "auc"] if c in sub.index else np.nan for c in clfs])
        lo = np.array([sub.loc[c, "auc_ci_low"] if c in sub.index else np.nan for c in clfs])
        hi = np.array([sub.loc[c, "auc_ci_high"] if c in sub.index else np.nan for c in clfs])
        err = np.vstack([np.clip(vals - lo, 0, None), np.clip(hi - vals, 0, None)])
        ax.bar(x + i * width - 0.4 + width / 2, vals, width,
               yerr=err, capsize=2.5,
               error_kw={"lw": 0.9, "ecolor": "#333333"},
               label=arm, color=ARM_COLORS.get(arm, "#888888"), edgecolor="white")

    ref = df[df["arm"] == REFERENCE_ARM].set_index("classifier")
    for j, c in enumerate(clfs):
        if c not in ref.index:
            continue
        base = ref.loc[c, "auc"]
        ax.text(x[j] - 0.4 + width / 2, base + 0.015, f"{base:.3f}",
                ha="center", fontsize=8)
        others = df[(df["arm"] != REFERENCE_ARM) & (df["classifier"] == c)]
        if not others.empty:
            top = others.loc[others["auc"].idxmax()]
            delta = top["auc"] - base
            ceiling = max(others["auc_ci_high"].max(), top["auc"])
            ax.text(x[j] + width, ceiling + 0.025, f"{delta:+.3f}",
                    ha="center", fontsize=8,
                    color="#b2182b" if delta > 0.02 else "#666666")

    ax.set_xticks(x)
    ax.set_xticklabels(clfs)
    ax.set_ylabel("AUC (out-of-fold, 95% bootstrap CI)")
    ax.set_ylim(0.4, 1.09)
    ax.axhline(0.5, color="#999999", lw=0.8, ls="--")
    ax.set_title(title, fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.read_csv(SOURCE_CSV)
    is_fusion = df["classifier"].str.endswith(CA125_SUFFIX)

    cfdna_df = df[~is_fusion & (df["arm"] != "ca125_only")]
    fusion_df = df[is_fusion]

    cfdna_clfs = list(dict.fromkeys(cfdna_df["classifier"]))
    fusion_clfs = list(dict.fromkeys(fusion_df["classifier"]))

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2))

    _panel(axes[0], cfdna_df, cfdna_clfs, "A — Standalone cfDNA classifiers")
    _panel(axes[1], fusion_df, fusion_clfs, "B — cfDNA + CA-125 fusion models")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=len(ARMS),
               loc="lower center", bbox_to_anchor=(0.5, -0.04))

    fig.suptitle("Quantification of feature-selection leakage under four protocols",
                 fontsize=13, y=1.02)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
