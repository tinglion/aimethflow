"""build_figure_s5.py — combine the two external-cohort leakage demonstrations
(OvaPrint, GSE19711) into one two-panel Figure S5, sharing the same axes,
colors and layout as Figure 4A (result/20260406/leakage_demo.png).

Run once from the project root:
    python src/scripts/build_figure_s5.py
"""

import os

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OVAPRINT_CSV = os.path.join(PROJECT_ROOT, "result", "external", "ovaprint", "leakage_demo.csv")
GSE19711_CSV = os.path.join(PROJECT_ROOT, "result", "external", "gse19711", "leakage_demo.csv")
OUT_PATH = os.path.join(PROJECT_ROOT, "docs", "GPB", "resources", "figure_s5.png")

# Same palette and y-limits as Figure 4A (leakage_demo.py::_save_plot).
ARM_COLORS = {"in_fold": "#4dac26", "naive_selection": "#e08214", "naive_full": "#d6604d"}
ARMS = ["in_fold", "naive_selection", "naive_full"]
Y_LIM = (0.4, 1.09)
REFERENCE_ARM = "in_fold"


def _panel(ax, df, title, n_label):
    clfs = list(dict.fromkeys(df["classifier"]))
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
    ax.set_ylim(*Y_LIM)
    ax.axhline(0.5, color="#999999", lw=0.8, ls="--")
    ax.set_title(f"{title} ({n_label})", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ova = pd.read_csv(OVAPRINT_CSV)
    gse = pd.read_csv(GSE19711_CSV)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    _panel(axes[0], ova, "OvaPrint", "n = 265")
    _panel(axes[1], gse, "GSE19711", "n = 405")

    axes[0].set_ylabel("AUC (out-of-fold, 95% bootstrap CI)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=len(ARMS),
               loc="lower center", bbox_to_anchor=(0.5, -0.02))

    fig.suptitle("Feature-selection leakage demonstration on external cohorts",
                 fontsize=13, y=1.02)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {OUT_PATH}")


if __name__ == "__main__":
    main()
