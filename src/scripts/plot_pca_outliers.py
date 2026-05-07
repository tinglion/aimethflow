#!/usr/bin/env python3
"""
Figure S2 — PCA scatter plot with outlier annotation.

Usage:
    python src/scripts/plot_pca_outliers.py [result_dir]
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
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler


# ── colour / style ──────────────────────────────────────────────────────────
COLOR = {"benign": "#4878CF", "malignant": "#D65F5F"}
MARKER = {"normal": "o", "suspect": "X"}
SIZE   = {"normal": 50,  "suspect": 120}
ALPHA  = {"normal": 0.75, "suspect": 1.0}

LABEL_OFFSET = (0.18, 0.18)       # data-unit nudge for outlier text


# ── helpers ─────────────────────────────────────────────────────────────────
def _latest_result_dir(base="result"):
    dirs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]
    dirs = [d for d in dirs if os.path.isfile(os.path.join(base, d, "feature_matrix.npz"))]
    if not dirs:
        raise FileNotFoundError(f"No result subdirectory with feature_matrix.npz under {base}/")
    return os.path.join(base, sorted(dirs)[-1])


# ── main ─────────────────────────────────────────────────────────────────────
def plot_pca_outliers(result_dir: str, out_path: Optional[str] = None):
    npz_path     = os.path.join(result_dir, "feature_matrix.npz")
    outlier_path = os.path.join(result_dir, "outlier_report.csv")

    if not os.path.exists(npz_path):
        sys.exit(f"[error] feature_matrix.npz not found: {npz_path}")

    # ── load feature matrix ──────────────────────────────────────────────────
    data       = np.load(npz_path, allow_pickle=True)
    X          = data["X"].astype(float)
    y          = data["y"]
    sample_ids = list(data["sample_ids"])

    # impute + scale before PCA
    X_imp  = SimpleImputer(strategy="median").fit_transform(X)
    X_scal = StandardScaler().fit_transform(X_imp)

    pca   = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_scal)
    var   = pca.explained_variance_ratio_ * 100

    # ── load outlier flags ───────────────────────────────────────────────────
    suspect_ids: set[str] = set()
    if os.path.exists(outlier_path):
        out_df = pd.read_csv(outlier_path)
        suspect_ids = set(out_df.loc[out_df["suspect"] == True, "sample_id"])

    # group lookup
    group_of = {sid: ("malignant" if int(label) == 1 else "benign")
                for sid, label in zip(sample_ids, y)}

    # ── plot ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(6, 5))

    for i, sid in enumerate(sample_ids):
        grp     = group_of[sid]
        is_out  = sid in suspect_ids
        key     = "suspect" if is_out else "normal"
        ax.scatter(
            coords[i, 0], coords[i, 1],
            c=COLOR[grp], marker=MARKER[key],
            s=SIZE[key], alpha=ALPHA[key],
            edgecolors="black" if is_out else "none", linewidths=0.8,
            zorder=3 if is_out else 2,
        )
        if is_out:
            ax.annotate(
                sid,
                xy=(coords[i, 0], coords[i, 1]),
                xytext=(coords[i, 0] + LABEL_OFFSET[0], coords[i, 1] + LABEL_OFFSET[1]),
                fontsize=7, color="black",
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.6),
            )

    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)", fontsize=11)
    ax.set_title("PCA of CpG Methylation Features", fontsize=12, fontweight="bold")
    ax.axhline(0, color="grey", lw=0.5, ls="--", zorder=0)
    ax.axvline(0, color="grey", lw=0.5, ls="--", zorder=0)

    legend_handles = [
        mpatches.Patch(color=COLOR["benign"],    label=f"Benign (n={sum(y==0)})"),
        mpatches.Patch(color=COLOR["malignant"], label=f"Malignant (n={sum(y==1)})"),
        plt.scatter([], [], marker="X", c="grey", s=80, edgecolors="black",
                    linewidths=0.8, label=f"Outlier excluded (n={len(suspect_ids)})"),
    ]
    ax.legend(handles=legend_handles, fontsize=9, framealpha=0.9)

    plt.tight_layout()

    if out_path is None:
        out_path = os.path.join(result_dir, "pca_outliers.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[done] saved → {out_path}")
    print(f"       PC1={var[0]:.1f}%  PC2={var[1]:.1f}%  "
          f"benign={sum(y==0)}  malignant={sum(y==1)}  outliers={len(suspect_ids)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", nargs="?", default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    rdir = args.result_dir or _latest_result_dir()
    plot_pca_outliers(rdir, args.out)
