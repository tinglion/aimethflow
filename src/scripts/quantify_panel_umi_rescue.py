#!/usr/bin/env python3
"""
Quantify what fraction of the final diagnostic-panel CpG sites would have
been lost under coordinate-only deduplication.

Approach
--------
Instead of re-running the entire pipeline with ``umi_tools dedup --ignore-umi``
(which would take ~70 hours on the 118-sample cohort), this script applies
the per-sample BAM-level rescue rates from ``quantify_umi_rescue.py``
(``docs/figures/umi_rescue_summary.csv``) directly to the per-CpG depth
matrix and tests, for each of the final panel CpG sites, whether the
≥10× coverage ≥80% sample threshold (Step 14 final hard filter) would still
have been met:

    D_coord_only[s, i] = D_umi_aware[s, i] × (1 − rescue_pct[s] / 100)

A panel CpG is "UMI-rescue-dependent" iff, in the simulated coord-only depth,
fewer than ``min_10x_coverage_rate`` (default 0.80) of samples retain ≥10×
coverage at that site. Total runtime: < 1 minute.

Why this is a defensible lower bound
------------------------------------
- Rescue is NOT uniformly distributed across genomic positions; it
  concentrates at high-density positions (which is exactly where our panel
  sites tend to be — they passed the ≥10× QC). The uniform scaling therefore
  *under-estimates* panel-site impact.
- The check only models the Step 14 Stage 5 hard filter; under reduced
  depth, β estimates would also become noisier and could destabilize the
  Wilcoxon / LASSO / RF selection upstream. That additional source of
  rescue-dependence is not modelled here.
- Therefore the reported count is a strict lower bound on the panel sites
  whose discovery genuinely required UMI-aware deduplication.

Outputs (written to ``<out_dir>/``)
-----------------------------------
- ``umi_rescue_panel_attribution.csv``    — Per-site original-vs-counterfactual
- ``SUMMARY.md``                          — Human-readable summary
"""

import argparse
import logging
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _load_config(config_path: str) -> dict:
    if sys.version_info >= (3, 11):
        import tomllib
        with open(config_path, "rb") as f:
            return tomllib.load(f)
    import toml
    return toml.load(config_path)


def _load_rescue_rates(rescue_csv: str) -> dict:
    df = pd.read_csv(rescue_csv)
    return dict(zip(df["sample"].astype(str), df["rescue_pct"].astype(float)))


def analyse_panel_depth(
    feature_matrix_path: str,
    panel_csv: str,
    rescue_rates: dict,
    min_10x_rate: float = 0.80,
    annotation_csv: Optional[str] = None,
) -> tuple:
    """
    For each panel CpG, simulate coord-only depth and check the ≥10× rate.

    Returns (summary_dict, attribution_DataFrame).
    """
    data = np.load(feature_matrix_path, allow_pickle=True)
    sample_ids   = [str(s) for s in data["sample_ids"]]
    feature_names = [str(f) for f in data["feature_names"]]
    D            = data["D"]  # (n_samples, n_features)

    # Map panel CpG → column index in the feature matrix
    panel_df       = pd.read_csv(panel_csv)
    panel_features = [str(f) for f in panel_df["feature"].tolist()]
    feat_to_col    = {fn: i for i, fn in enumerate(feature_names)}
    panel_cols     = [feat_to_col[f] for f in panel_features if f in feat_to_col]
    panel_in_matrix = [feature_names[c] for c in panel_cols]
    missing         = [f for f in panel_features if f not in feat_to_col]
    if missing:
        logger.warning(f"  {len(missing)} panel CpGs not found in feature_matrix "
                       f"(likely the matrix is from a different run): {missing[:5]} ...")

    # Per-sample scaling factor (samples without measured rescue → 1.0)
    scaling = np.array([
        1.0 - rescue_rates.get(sid, 0.0) / 100.0 for sid in sample_ids
    ], dtype=float)
    scaling = np.maximum(scaling, 0.0)

    # Annotation lookup (optional)
    ann_lookup = {}
    if annotation_csv and os.path.exists(annotation_csv):
        ann_df = pd.read_csv(annotation_csv)
        ann_lookup = {
            str(r["feature"]): r.to_dict() for _, r in ann_df.iterrows()
        }

    # Per-CpG attribution
    rows = []
    for feat, col in zip(panel_in_matrix, panel_cols):
        depths        = D[:, col]                                  # original
        depths_coord  = depths * scaling                           # simulated
        n_total       = (depths > 0).sum()                         # samples with any read
        n_10x_umi     = (depths >= 10).sum()
        n_10x_coord   = (depths_coord >= 10).sum()
        rate_umi      = n_10x_umi / len(depths)
        rate_coord    = n_10x_coord / len(depths)
        passes_umi    = rate_umi  >= min_10x_rate
        passes_coord  = rate_coord >= min_10x_rate
        rescue_dep    = passes_umi and not passes_coord
        ann = ann_lookup.get(feat, {})
        rows.append({
            "feature":               feat,
            "n_samples":             len(depths),
            "median_depth_umi":      float(np.median(depths)),
            "median_depth_coord":    float(np.median(depths_coord)),
            "n_10x_umi":             int(n_10x_umi),
            "n_10x_coord":           int(n_10x_coord),
            "rate_10x_umi":          round(float(rate_umi),   4),
            "rate_10x_coord":        round(float(rate_coord), 4),
            "passes_umi":            bool(passes_umi),
            "passes_coord_only":     bool(passes_coord),
            "umi_rescue_dependent":  bool(rescue_dep),
            "nearest_gene":          ann.get("nearest_gene", ""),
            "in_88_regions":         ann.get("in_88_regions", ""),
            "selection_rate":        ann.get("selection_rate", ""),
            "ewas_n_hits":           ann.get("ewas_n_hits", ""),
        })

    df = pd.DataFrame(rows)

    summary = {
        "n_panel_total":            len(panel_features),
        "n_panel_in_matrix":        len(panel_in_matrix),
        "n_panel_missing":          len(missing),
        "n_passes_umi":             int(df["passes_umi"].sum()),
        "n_passes_coord_only":      int(df["passes_coord_only"].sum()),
        "n_umi_rescue_dependent":   int(df["umi_rescue_dependent"].sum()),
        "fraction_rescue_dependent": (
            df["umi_rescue_dependent"].sum() / max(1, len(df))
        ),
        "min_10x_rate_threshold":   min_10x_rate,
    }
    return summary, df


def write_summary_md(summary: dict, attribution: pd.DataFrame,
                     rescue_csv_path: str, out_path: str) -> None:
    lines = [
        "# UMI-aware deduplication: counterfactual panel attribution",
        "",
        "## Approach",
        "",
        f"For each of the **{summary['n_panel_total']}** CpG sites in the final "
        f"UMI-aware diagnostic panel, the per-sample per-CpG depth was scaled "
        f"down by `(1 − rescue_pct / 100)` using the BAM-level per-sample "
        f"rescue rates from `{os.path.basename(rescue_csv_path)}`. A panel CpG "
        f"is flagged **UMI-rescue-dependent** iff it passes the "
        f"≥10× coverage ≥{summary['min_10x_rate_threshold']:.0%} sample threshold "
        f"under the real UMI-aware depth but FAILS under the simulated "
        f"coord-only depth.",
        "",
        "## Results",
        "",
        "| Quantity | Value |",
        "|---|---|",
        f"| Final UMI-aware panel size | {summary['n_panel_total']} |",
        f"| Sites present in feature_matrix (matched) | {summary['n_panel_in_matrix']} |",
        f"| Sites passing ≥10× ≥80% under UMI-aware depth | {summary['n_passes_umi']} |",
        f"| Sites passing ≥10× ≥80% under simulated coord-only depth | {summary['n_passes_coord_only']} |",
        f"| **UMI-rescue-dependent (lost in coord-only)** | **{summary['n_umi_rescue_dependent']}** |",
        f"| Fraction of panel that is rescue-dependent | "
        f"**{summary['fraction_rescue_dependent']:.1%}** |",
        "",
    ]

    rescue_dep = attribution[attribution["umi_rescue_dependent"]]
    if not rescue_dep.empty:
        lines += [
            "## UMI-rescue-dependent panel CpG sites",
            "",
            "These panel CpGs pass the ≥10× ≥80%-sample QC under the "
            "UMI-aware depth but fail under simulated coord-only depth — "
            "i.e., they would not have entered the final diagnostic panel "
            "without UMI rescue.",
            "",
            "| Feature | Nearest gene | In 88 OC regions | Selection rate | "
            "≥10× rate UMI | ≥10× rate coord | EWAS hits |",
            "|---|---|---|---|---|---|---|",
        ]
        for _, row in rescue_dep.iterrows():
            lines.append(
                f"| `{row['feature']}` | "
                f"{row.get('nearest_gene', '?')} | "
                f"{row.get('in_88_regions', '?')} | "
                f"{row.get('selection_rate', '?')} | "
                f"{row['rate_10x_umi']:.3f} | "
                f"{row['rate_10x_coord']:.3f} | "
                f"{row.get('ewas_n_hits', '?')} |"
            )
        lines.append("")
    else:
        lines += [
            "## UMI-rescue-dependent panel CpG sites",
            "",
            "**None of the panel CpGs fail the ≥10× ≥80% threshold under "
            "the simulated coord-only depth.** This indicates the panel "
            "sites have sufficient depth headroom that the BAM-level UMI "
            "rescue does not gate panel inclusion via this hard filter; "
            "the value of UMI-aware dedup lies elsewhere (β-estimate "
            "precision, Wilcoxon power, LASSO/RF stability) — none of "
            "which are modelled by this strict depth-only check.",
            "",
        ]

    lines += [
        "## Caveats — why this is a strict LOWER BOUND",
        "",
        "1. **Non-uniform rescue.** BAM-level rescue concentrates at "
        "high-density genomic positions, which is exactly where panel "
        "sites tend to be (they passed the ≥10× QC). Applying a uniform "
        "per-sample scaling factor therefore *under-estimates* the depth "
        "loss at panel sites; the true rescue-dependence is likely higher.",
        "2. **Depth-QC only.** This check models only the final Stage 5 "
        "hard filter. Under reduced depth, β estimates become noisier, "
        "which would also destabilize Wilcoxon / LASSO / RF rankings "
        "upstream. That additional rescue dependence is not captured here.",
        "3. **Held-constant β.** β values are not perturbed in the simulation.",
        "",
        "## Reproducibility",
        "",
        "```bash",
        "python src/scripts/quantify_panel_umi_rescue.py \\",
        "    --result-dir result/<workdir> \\",
        "    --rescue-csv docs/figures/umi_rescue_summary.csv \\",
        "    --config    config.toml",
        "```",
        "",
    ]

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--result-dir", required=True,
                        help="Result dir with feature_matrix.npz + panel_cpg_blood_filtered.csv")
    parser.add_argument("--rescue-csv", default="docs/figures/umi_rescue_summary.csv")
    parser.add_argument("--config",     default="config.toml")
    parser.add_argument("--panel-csv",  default=None,
                        help="Default: panel_cpg_blood_filtered.csv (or panel_cpg.csv as fallback)")
    parser.add_argument("--out-dir",    default=None,
                        help="Default: <result-dir>/umi_rescue_attribution")
    args = parser.parse_args()

    if args.out_dir is None:
        args.out_dir = os.path.join(args.result_dir, "umi_rescue_attribution")
    os.makedirs(args.out_dir, exist_ok=True)

    if args.panel_csv is None:
        candidate = os.path.join(args.result_dir, "panel_cpg_blood_filtered.csv")
        if os.path.exists(candidate):
            args.panel_csv = candidate
        else:
            args.panel_csv = os.path.join(args.result_dir, "panel_cpg.csv")

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    config       = _load_config(args.config)
    min_10x_rate = float(config.get("ml", {}).get("min_10x_coverage_rate", 0.80))
    logger.info(f"Using ≥10× sample-rate threshold = {min_10x_rate:.2f} (from config)")

    rescue_rates = _load_rescue_rates(args.rescue_csv)
    rate_values  = list(rescue_rates.values())
    logger.info(
        f"Loaded rescue rates for {len(rescue_rates)} samples: "
        f"median {np.median(rate_values):.2f}%, range "
        f"{min(rate_values):.2f}–{max(rate_values):.2f}%"
    )

    feature_matrix = os.path.join(args.result_dir, "feature_matrix.npz")
    annotation     = os.path.join(args.result_dir, "panel_annotation.csv")

    summary, attribution = analyse_panel_depth(
        feature_matrix, args.panel_csv, rescue_rates,
        min_10x_rate=min_10x_rate,
        annotation_csv=annotation,
    )

    attribution_csv = os.path.join(args.out_dir, "umi_rescue_panel_attribution.csv")
    attribution.to_csv(attribution_csv, index=False)

    summary_md = os.path.join(args.out_dir, "SUMMARY.md")
    write_summary_md(summary, attribution, args.rescue_csv, summary_md)

    print()
    print("============================================================")
    print(" UMI-rescue panel attribution — summary")
    print("============================================================")
    print(f"  Panel (input):                       {args.panel_csv}")
    print(f"  Final panel size:                  {summary['n_panel_total']:>4} CpGs")
    print(f"  Matched in feature_matrix:         {summary['n_panel_in_matrix']:>4} CpGs")
    print(f"  Passes ≥10× ≥{int(min_10x_rate*100)}% under UMI-aware:    "
          f"{summary['n_passes_umi']:>4} CpGs")
    print(f"  Passes ≥10× ≥{int(min_10x_rate*100)}% under coord-only:   "
          f"{summary['n_passes_coord_only']:>4} CpGs")
    print(f"  UMI-rescue-dependent (lower bound): {summary['n_umi_rescue_dependent']:>4} CpGs"
          f"  ({summary['fraction_rescue_dependent']:.1%})")
    print()
    print(f"  Summary:       {summary_md}")
    print(f"  Per-site CSV:  {attribution_csv}")
    print("============================================================")


if __name__ == "__main__":
    main()
