#!/usr/bin/env python3
"""
Figure 2A — Feature selection funnel diagram.

Reads counts directly from result files so numbers stay in sync with each run.

Usage:
    python src/scripts/plot_feature_funnel.py [result_dir]
    result_dir defaults to the most recent result/ subdirectory.
"""

import argparse
import os
import sys
from typing import Optional

import matplotlib
if not matplotlib.is_interactive():
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import pandas as pd


# ── colour palette ───────────────────────────────────────────────────────────
BOX_COLORS = [
    "#2B6CB0",  # 0 — total
    "#2F855A",  # 1 — wilcoxon
    "#6B46C1",  # 2 — lasso
    "#6B46C1",  # 3 — rf
    "#B7791F",  # 4 — union
    "#C53030",  # 5 — panel (final)
]
TEXT_COLOR  = "white"
ARROW_COLOR = "#555555"
FONT_MAIN   = 10
FONT_SUB    = 8.5


# ── helper ───────────────────────────────────────────────────────────────────
def _latest_result_dir(base="result"):
    dirs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]
    dirs = [d for d in dirs if os.path.isfile(os.path.join(base, d, "feature_importance.csv"))]
    if not dirs:
        raise FileNotFoundError(f"No result dir with feature_importance.csv under {base}/")
    return os.path.join(base, sorted(dirs)[-1])


def _load_counts(result_dir: str) -> dict:
    fi_path   = os.path.join(result_dir, "feature_importance.csv")
    npz_path  = os.path.join(result_dir, "feature_matrix.npz")
    panel_path = os.path.join(result_dir, "panel_cpg.csv")

    if not os.path.exists(fi_path):
        sys.exit(f"[error] feature_importance.csv not found: {fi_path}")

    fi = pd.read_csv(fi_path)

    # total CpG sites from feature_matrix if available, else from wilcoxon file
    total = None
    if os.path.exists(npz_path):
        d = np.load(npz_path, allow_pickle=True)
        total = int(d["X"].shape[1])

    wilcox_n = len(fi)
    lasso_n  = int(fi["lasso_selected"].sum())
    rf_n     = int(fi["rf_selected"].sum())
    union_n  = int((fi["lasso_selected"] | fi["rf_selected"]).sum())
    panel_n  = (int(fi["panel_selected"].sum())
                if "panel_selected" in fi.columns
                else len(pd.read_csv(panel_path)))

    if total is None:
        # fallback: total assumed proportional (wilcox_n / pass_rate typical ~0.45)
        total = wilcox_n
        print("[warn] feature_matrix.npz not found; showing Wilcoxon input as top node")

    return dict(total=total, wilcox=wilcox_n,
                lasso=lasso_n, rf=rf_n, union=union_n, panel=panel_n)


# ── drawing ──────────────────────────────────────────────────────────────────
def _fmt(n: int) -> str:
    return f"{n:,}"


def _draw_box(ax, cx, cy, w, h, color, label, sublabel, fontsize_main=FONT_MAIN):
    box = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.04", linewidth=1.2,
        edgecolor="white", facecolor=color, zorder=3,
    )
    ax.add_patch(box)
    ax.text(cx, cy + 0.06, label,
            ha="center", va="center", fontsize=fontsize_main,
            fontweight="bold", color=TEXT_COLOR, zorder=4)
    if sublabel:
        ax.text(cx, cy - 0.14, sublabel,
                ha="center", va="center", fontsize=FONT_SUB,
                color=TEXT_COLOR, alpha=0.88, zorder=4)


def _draw_arrow(ax, x, y_top, y_bot, label=""):
    ax.annotate(
        "", xy=(x, y_bot + 0.01), xytext=(x, y_top - 0.01),
        arrowprops=dict(arrowstyle="-|>", color=ARROW_COLOR,
                        lw=1.4, mutation_scale=14),
        zorder=2,
    )
    if label:
        ax.text(x + 0.26, (y_top + y_bot) / 2, label,
                ha="left", va="center", fontsize=7.5,
                color=ARROW_COLOR, style="italic")


# ── main ─────────────────────────────────────────────────────────────────────
def plot_funnel(result_dir: str, out_path: Optional[str] = None):
    c = _load_counts(result_dir)

    # layout constants
    W, H   = 2.6, 0.42    # box width, height  (in data units)
    CX     = 0.0           # horizontal centre
    GAP    = 0.22          # arrow gap between boxes
    STEP   = H + GAP + 0.18

    # Y positions (top to bottom)
    rows = [
        # (y_centre,  colour_idx, label,            sublabel)
        (0.0,          0, f"{_fmt(c['total'])} CpG sites",      "all candidate sites"),
        (-STEP,        1, f"{_fmt(c['wilcox'])} sites",          "Wilcoxon pre-filter (FDR < 0.05, |Δβ| > 0.1)"),
        (-STEP*2,      2, f"LASSO: {_fmt(c['lasso'])}",          "L1 logistic regression"),
        (-STEP*2 - (H + GAP*0.5), 3, f"RF: {_fmt(c['rf'])}",    "Random Forest importance"),
        (-STEP*3,      4, f"Union: {_fmt(c['union'])} sites",    "LASSO ∪ RF"),
        (-STEP*4,      5, f"{_fmt(c['panel'])} panel CpGs",      "≥10× depth in ≥80% samples  ★ final panel"),
    ]

    fig_h = STEP * 4 + H * 2 + 1.2
    fig, ax = plt.subplots(figsize=(4.2, fig_h))
    ax.set_xlim(-1.8, 1.8)
    ax.set_ylim(-STEP * 4 - H - 0.6, H + 0.4)
    ax.axis("off")

    y_prev = None
    for idx, (yc, cidx, lab, sub) in enumerate(rows):
        # split LASSO/RF row: draw side-by-side, slightly narrower
        is_split = idx in (2, 3)
        w = W * 0.72 if is_split else W
        offset = -W * 0.38 if idx == 2 else (W * 0.38 if idx == 3 else 0.0)
        _draw_box(ax, CX + offset, yc, w, H,
                  BOX_COLORS[cidx], lab, sub,
                  fontsize_main=FONT_MAIN - (0.5 if is_split else 0))

        # arrows
        if y_prev is not None and idx not in (3,):
            # from previous box bottom to this box top
            src_y = y_prev - H / 2
            dst_y = yc + H / 2
            if idx == 4:
                # two arrows: from lasso box and rf box converging
                lasso_y = rows[2][0] - H / 2
                rf_y    = rows[3][0] - H / 2
                lasso_cx = CX - W * 0.38
                rf_cx    = CX + W * 0.38
                for sx, sy in [(lasso_cx, lasso_y), (rf_cx, rf_y)]:
                    ax.annotate(
                        "", xy=(CX, dst_y + 0.01), xytext=(sx, sy - 0.01),
                        arrowprops=dict(arrowstyle="-|>", color=ARROW_COLOR,
                                        lw=1.2, mutation_scale=12,
                                        connectionstyle="arc3,rad=0.0"),
                        zorder=2,
                    )
            else:
                src_cx = CX + (CX - W * 0.38 if idx == 3 else 0.0)
                _draw_arrow(ax, CX + offset, src_y, dst_y)

        if idx == 1:
            y_prev = yc   # wilcoxon → feed both lasso and rf
        elif idx not in (2, 3):
            y_prev = yc

    # % retained annotation on right side
    pct_wilcox = c["wilcox"] / c["total"] * 100
    pct_union  = c["union"]  / c["total"] * 100
    pct_panel  = c["panel"]  / c["total"] * 100
    for yc, pct, color in [
        (rows[1][0], pct_wilcox, BOX_COLORS[1]),
        (rows[4][0], pct_union,  BOX_COLORS[4]),
        (rows[5][0], pct_panel,  BOX_COLORS[5]),
    ]:
        ax.text(W / 2 + 0.22, yc, f"{pct:.1f}%",
                ha="left", va="center", fontsize=7.5,
                color=color, fontweight="bold")

    ax.set_title("Feature Selection Funnel", fontsize=12,
                 fontweight="bold", pad=8)

    plt.tight_layout()
    if out_path is None:
        out_path = os.path.join(result_dir, "feature_selection_funnel.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] saved → {out_path}")
    print(f"       {c['total']:,} → {c['wilcox']:,} → LASSO {c['lasso']} + RF {c['rf']} "
          f"→ union {c['union']} → panel {c['panel']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", nargs="?", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    rdir = args.result_dir or _latest_result_dir()
    plot_funnel(rdir, args.out)
