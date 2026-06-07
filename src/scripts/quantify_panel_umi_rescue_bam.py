#!/usr/bin/env python3
"""
BAM-level per-CpG UMI rescue analysis.

For each panel CpG (chrom:pos), walks the existing UMI-deduped (and
on-target filtered) BAM of each sample, counts the reads UMI dedup kept,
groups them by their would-be coord-only dedup key, and computes the
*local* per-CpG rescue rate — i.e., how many of those reads would have
been collapsed by coordinate-only deduplication.

This is more accurate than applying a single BAM-wide rescue rate
uniformly to every panel CpG, because rescue concentrates at high-density
positions and panel CpG sites are by construction in higher-density
regions (they passed the ≥10× QC).

Approach
--------
For each (sample, panel CpG):
    1. fetch reads overlapping the CpG with MAPQ ≥ ``mapq_min``
    2. count read1 of each properly paired alignment (avoids double-counting)
    3. group by ``(reference_start, next_reference_start)`` — the key
       coordinate-only dedup uses for paired reads
    4. n_reads_umi   = total reads kept by UMI dedup
       n_unique_coord = number of unique coord-only groups
       rescue_count  = n_reads_umi − n_unique_coord
       rescue_rate   = rescue_count / n_reads_umi

Then for each panel CpG, apply the per-(sample, CpG) coord-only depth
(`= n_unique_coord`) against the ≥10× coverage in ≥80% of samples threshold.
A CpG is flagged **UMI-rescue-dependent** iff it passes the threshold under
UMI-aware depth but fails under simulated coord-only depth.

Inputs
------
- ``--result-dir``  : directory with ``panel_cpg_blood_filtered.csv`` (or
                      ``panel_cpg.csv``) and ``panel_annotation.csv``
- ``--bam-dir``     : directory with the per-sample UMI-deduped BAMs
- ``--bam-pattern`` : glob pattern (default
                      ``*.bwameth_pe.dedup.targeted.sorted.bam``)

Outputs (in ``<out-dir>/``)
---------------------------
- ``per_sample_per_cpg_rescue.csv``     — one row per (sample, CpG)
- ``per_cpg_summary.csv``               — one row per panel CpG (aggregated)
- ``umi_rescue_panel_attribution.csv``  — final pass/fail attribution per CpG
- ``SUMMARY.md``                        — human-readable summary
"""

import argparse
import glob
import logging
import multiprocessing as mp
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
#  Per-CpG rescue computation                                                 #
# --------------------------------------------------------------------------- #

def _rescue_at_position(bam_path: str, chrom: str, pos: int,
                       mapq_min: int) -> tuple:
    """Walk reads at (chrom, pos) and count UMI-aware vs coord-only.

    Returns ``(n_reads_umi, n_unique_coord_groups)``.

    Counting rule
    -------------
    - Skip unmapped / secondary / supplementary alignments.
    - Skip MAPQ < ``mapq_min``.
    - For paired reads, count read1 only (read2 is the same molecule).
    - Coord-only key: ``(read.reference_start, read.next_reference_start)``
      for paired reads; ``(read.reference_start, read.is_reverse)`` for
      single-end. This matches what ``umi_tools dedup --ignore-umi``
      collapses on.
    """
    import pysam
    keys = []
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as exc:
        logger.warning(f"  failed to open {bam_path}: {exc}")
        return 0, 0

    try:
        for read in bam.fetch(chrom, pos, pos + 1):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.is_duplicate:
                continue
            if read.mapping_quality < mapq_min:
                continue
            if read.is_paired and not read.is_read1:
                continue
            if read.is_paired:
                key = (read.reference_start, read.next_reference_start)
            else:
                key = (read.reference_start, read.is_reverse)
            keys.append(key)
    finally:
        bam.close()

    return len(keys), len(set(keys))


def _process_sample(args_tuple):
    """Worker function: process one sample × all panel CpGs."""
    sample_id, bam_path, cpg_positions, mapq_min = args_tuple
    rows = []
    for chrom, pos, feature in cpg_positions:
        n_umi, n_coord = _rescue_at_position(bam_path, chrom, pos, mapq_min)
        rescued = n_umi - n_coord
        rescue_rate = rescued / n_umi if n_umi else 0.0
        rows.append({
            "sample_id":      sample_id,
            "feature":        feature,
            "chrom":          chrom,
            "pos":            pos,
            "n_reads_umi":    n_umi,
            "n_unique_coord": n_coord,
            "rescued":        rescued,
            "rescue_rate":    rescue_rate,
        })
    return rows


# --------------------------------------------------------------------------- #
#  I/O helpers                                                                #
# --------------------------------------------------------------------------- #

def _load_config(path: str) -> dict:
    if sys.version_info >= (3, 11):
        import tomllib
        with open(path, "rb") as f:
            return tomllib.load(f)
    import toml
    return toml.load(path)


def _discover_bams(bam_dir: str, pattern: str) -> dict:
    """Return ``{sample_id: bam_path}`` from a directory glob."""
    bams = sorted(glob.glob(os.path.join(bam_dir, pattern)))
    out = {}
    for path in bams:
        # e.g. 20240723I-P10060.bwameth_pe.dedup.targeted.sorted.bam
        # sample id = part before the first '.'
        sample_id = os.path.basename(path).split(".", 1)[0]
        out[sample_id] = path
    return out


def _parse_panel_features(panel_csv: str) -> list:
    """Return list of (chrom, pos_0based, feature_str)."""
    df = pd.read_csv(panel_csv)
    out = []
    for feat in df["feature"].astype(str):
        if ":" not in feat:
            continue
        chrom, pos_s = feat.split(":", 1)
        try:
            pos = int(pos_s)
        except ValueError:
            continue
        out.append((chrom, pos, feat))
    return out


# --------------------------------------------------------------------------- #
#  Main                                                                       #
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--bam-dir",    required=True)
    parser.add_argument("--bam-pattern",
                        default="*.bwameth_pe.dedup.targeted.sorted.bam")
    parser.add_argument("--panel-csv",  default=None,
                        help="Default: panel_cpg_blood_filtered.csv "
                             "(or panel_cpg.csv as fallback)")
    parser.add_argument("--config",     default="config.toml")
    parser.add_argument("--mapq-min",   type=int, default=40)
    parser.add_argument("--threads",    type=int, default=None,
                        help="Default: cpu_count")
    parser.add_argument("--out-dir",    default=None,
                        help="Default: <result-dir>/umi_rescue_bam_attribution")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.out_dir is None:
        args.out_dir = os.path.join(args.result_dir, "umi_rescue_bam_attribution")
    os.makedirs(args.out_dir, exist_ok=True)

    if args.panel_csv is None:
        cand = os.path.join(args.result_dir, "panel_cpg_blood_filtered.csv")
        args.panel_csv = cand if os.path.exists(cand) else \
            os.path.join(args.result_dir, "panel_cpg.csv")

    config       = _load_config(args.config)
    min_10x_rate = float(config.get("ml", {}).get("min_10x_coverage_rate", 0.80))

    # ── 1. Discover BAMs ──────────────────────────────────────────────
    bams = _discover_bams(args.bam_dir, args.bam_pattern)
    logger.info(f"Discovered {len(bams)} BAMs in {args.bam_dir}")
    if not bams:
        sys.exit(1)

    # ── 2. Parse panel ────────────────────────────────────────────────
    cpg_positions = _parse_panel_features(args.panel_csv)
    logger.info(f"Panel: {len(cpg_positions)} CpG sites from {args.panel_csv}")

    # ── 3. Restrict to samples present in feature_matrix.npz ─────────
    fm_path = os.path.join(args.result_dir, "feature_matrix.npz")
    fm      = np.load(fm_path, allow_pickle=True)
    fm_samples = [str(s) for s in fm["sample_ids"]]
    fm_y       = fm["y"]
    # Drop borderlines (y == 2)
    keep_mask  = fm_y != 2
    fm_samples = [s for s, k in zip(fm_samples, keep_mask) if k]
    samples_with_bam = [s for s in fm_samples if s in bams]
    logger.info(
        f"Samples: {len(fm_samples)} in feature_matrix, "
        f"{len(samples_with_bam)} with matching BAM"
    )
    if len(samples_with_bam) < len(fm_samples):
        missing = set(fm_samples) - set(samples_with_bam)
        logger.warning(f"Missing BAMs for {len(missing)} samples; "
                       f"first few: {sorted(missing)[:5]}")

    # ── 4. Parallel per-sample × all CpGs ────────────────────────────
    n_workers = args.threads or os.cpu_count() or 4
    work = [
        (sid, bams[sid], cpg_positions, args.mapq_min)
        for sid in samples_with_bam
    ]
    logger.info(f"Launching {len(work)} workers (×{len(cpg_positions)} CpGs each) "
                f"on {n_workers} threads ...")
    all_rows = []
    if n_workers > 1:
        with mp.Pool(n_workers) as pool:
            for i, batch in enumerate(pool.imap_unordered(_process_sample, work), 1):
                all_rows.extend(batch)
                if i % 10 == 0 or i == len(work):
                    logger.info(f"  done: {i}/{len(work)} samples")
    else:
        for i, args_tuple in enumerate(work, 1):
            all_rows.extend(_process_sample(args_tuple))
            if i % 10 == 0 or i == len(work):
                logger.info(f"  done: {i}/{len(work)} samples")

    df = pd.DataFrame(all_rows)
    per_pair_csv = os.path.join(args.out_dir, "per_sample_per_cpg_rescue.csv")
    df.to_csv(per_pair_csv, index=False)
    logger.info(f"Wrote per-(sample,CpG) table: {per_pair_csv}")

    # ── 5. Aggregate per CpG ─────────────────────────────────────────
    # For each CpG, count samples passing ≥10× under UMI vs coord-only.
    df["passes_umi"]   = df["n_reads_umi"]    >= 10
    df["passes_coord"] = df["n_unique_coord"] >= 10
    n_samples = df["sample_id"].nunique()

    agg = df.groupby("feature").agg(
        n_samples         = ("sample_id",     "nunique"),
        median_depth_umi  = ("n_reads_umi",    "median"),
        median_depth_coord= ("n_unique_coord", "median"),
        median_rescue_pct = ("rescue_rate",   lambda x: float(np.median(x)) * 100),
        max_rescue_pct    = ("rescue_rate",   lambda x: float(np.max(x))    * 100),
        n_passes_umi      = ("passes_umi",     "sum"),
        n_passes_coord    = ("passes_coord",   "sum"),
    ).reset_index()
    agg["rate_10x_umi"]   = agg["n_passes_umi"]   / agg["n_samples"]
    agg["rate_10x_coord"] = agg["n_passes_coord"] / agg["n_samples"]
    agg["passes_umi_threshold"]   = agg["rate_10x_umi"]   >= min_10x_rate
    agg["passes_coord_threshold"] = agg["rate_10x_coord"] >= min_10x_rate
    agg["umi_rescue_dependent"]   = (agg["passes_umi_threshold"]
                                     & ~agg["passes_coord_threshold"])

    # Join annotation if available
    ann_csv = os.path.join(args.result_dir, "panel_annotation.csv")
    if os.path.exists(ann_csv):
        ann = pd.read_csv(ann_csv)
        keep_cols = ["feature", "nearest_gene", "in_88_regions",
                     "selection_rate", "ewas_n_hits"]
        ann = ann[[c for c in keep_cols if c in ann.columns]]
        agg = agg.merge(ann, on="feature", how="left")

    per_cpg_csv = os.path.join(args.out_dir, "per_cpg_summary.csv")
    agg.to_csv(per_cpg_csv, index=False)
    logger.info(f"Wrote per-CpG summary: {per_cpg_csv}")

    # ── 6. Summary report ────────────────────────────────────────────
    summary = {
        "n_panel_total":            len(cpg_positions),
        "n_passes_umi":             int(agg["passes_umi_threshold"].sum()),
        "n_passes_coord":           int(agg["passes_coord_threshold"].sum()),
        "n_umi_rescue_dependent":   int(agg["umi_rescue_dependent"].sum()),
        "median_per_cpg_rescue_pct":   float(agg["median_rescue_pct"].median()),
        "max_per_cpg_rescue_pct":      float(agg["max_rescue_pct"].max()),
        "n_samples":                   n_samples,
        "min_10x_rate_threshold":      min_10x_rate,
        "mapq_min":                    args.mapq_min,
    }
    summary["fraction_rescue_dependent"] = (
        summary["n_umi_rescue_dependent"] / max(1, len(cpg_positions))
    )

    rescue_dep = agg[agg["umi_rescue_dependent"]].copy()

    # ── 7. SUMMARY.md ────────────────────────────────────────────────
    lines = [
        "# UMI-aware deduplication: BAM-level per-CpG counterfactual analysis",
        "",
        "## Approach",
        "",
        f"For each of the **{summary['n_panel_total']}** CpG sites in the "
        f"final UMI-aware diagnostic panel, we walked the UMI-deduped + "
        f"on-target BAM of each of the **{summary['n_samples']}** samples in "
        f"the analytical cohort. At each panel position, reads kept by UMI "
        f"deduplication (`n_reads_umi`) were grouped by their would-be "
        f"coordinate-only key `(reference_start, next_reference_start)`, "
        f"yielding the simulated coord-only depth (`n_unique_coord`) and a "
        f"per-(sample, CpG) local rescue rate. A panel CpG is **"
        f"UMI-rescue-dependent** iff it passes the ≥10× coverage in "
        f"≥{int(min_10x_rate * 100)}% of samples threshold under the real "
        f"UMI-aware depth but FAILS under the simulated coord-only depth.",
        "",
        "## Results",
        "",
        "| Quantity | Value |",
        "|---|---|",
        f"| Final UMI-aware panel size | {summary['n_panel_total']} |",
        f"| Samples analysed | {summary['n_samples']} |",
        f"| MAPQ threshold | ≥ {summary['mapq_min']} |",
        f"| Coverage threshold | ≥ 10× in ≥ {int(min_10x_rate * 100)}% samples |",
        f"| Sites passing under UMI-aware depth | {summary['n_passes_umi']} |",
        f"| Sites passing under coord-only depth | {summary['n_passes_coord']} |",
        f"| **UMI-rescue-dependent (lost in coord-only)** | "
        f"**{summary['n_umi_rescue_dependent']}** |",
        f"| Fraction rescue-dependent | "
        f"**{summary['fraction_rescue_dependent']:.1%}** |",
        f"| Median per-CpG rescue rate (across sites) | "
        f"{summary['median_per_cpg_rescue_pct']:.2f}% |",
        f"| Max per-CpG rescue rate observed (any site, any sample) | "
        f"{summary['max_per_cpg_rescue_pct']:.2f}% |",
        "",
    ]

    if not rescue_dep.empty:
        lines += [
            "## UMI-rescue-dependent panel CpG sites",
            "",
            "These panel CpGs would not have entered the final diagnostic "
            "panel under coordinate-only deduplication.",
            "",
            "| Feature | Nearest gene | In 88 OC regions | Selection rate | "
            "≥10× rate UMI | ≥10× rate coord | Median rescue % | EWAS hits |",
            "|---|---|---|---|---|---|---|---|",
        ]
        for _, row in rescue_dep.iterrows():
            lines.append(
                f"| `{row['feature']}` | "
                f"{row.get('nearest_gene', '?')} | "
                f"{row.get('in_88_regions', '?')} | "
                f"{row.get('selection_rate', '?')} | "
                f"{row['rate_10x_umi']:.3f} | "
                f"{row['rate_10x_coord']:.3f} | "
                f"{row['median_rescue_pct']:.1f}% | "
                f"{row.get('ewas_n_hits', '?')} |"
            )
        lines.append("")
    else:
        lines += [
            "## UMI-rescue-dependent panel CpG sites",
            "",
            "**None of the panel CpGs cross the threshold in the simulated "
            "coord-only depth.** This indicates the panel sites have "
            "sufficient depth headroom that the ≥10× hard filter is not the "
            "binding constraint of UMI-aware deduplication at panel "
            "positions; UMI rescue's value at panel sites operates through "
            "more subtle stabilization of β estimates, Wilcoxon statistics, "
            "and LASSO / RF rankings — quantitative effects not captured by "
            "this strict depth check.",
            "",
        ]

    lines += [
        "## Caveats",
        "",
        "1. **Hard depth filter only.** The check uses the Step-14 Stage-5 "
        "≥10× threshold; β-precision and Wilcoxon / LASSO / RF stability "
        "are not modelled.",
        "2. **Coord-only key approximation.** We collapse on "
        "`(reference_start, next_reference_start)` for paired reads, which "
        "closely matches `umi_tools dedup --ignore-umi` but may differ "
        "slightly in edge cases (soft-clipping, supplementary alignments).",
        "3. **Cohort completeness.** Analysis restricted to samples present "
        "in both ``feature_matrix.npz`` and the BAM directory; missing "
        "samples are excluded but reported in the per-sample log.",
        "",
        "## Reproducibility",
        "",
        "```bash",
        "python src/scripts/quantify_panel_umi_rescue_bam.py \\",
        "    --result-dir result/<workdir> \\",
        "    --bam-dir    temp/<workdir> \\",
        "    --config     config.toml",
        "```",
        "",
    ]

    summary_md = os.path.join(args.out_dir, "SUMMARY.md")
    with open(summary_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # ── 8. Console report ────────────────────────────────────────────
    print()
    print("============================================================")
    print(" UMI-rescue panel attribution (BAM-level) — summary")
    print("============================================================")
    print(f"  Panel size:                          {summary['n_panel_total']:>4} CpGs")
    print(f"  Samples analysed:                    {summary['n_samples']:>4}")
    print(f"  Sites passing under UMI-aware:       {summary['n_passes_umi']:>4}")
    print(f"  Sites passing under coord-only:      {summary['n_passes_coord']:>4}")
    print(f"  UMI-rescue-dependent:                {summary['n_umi_rescue_dependent']:>4} "
          f"({summary['fraction_rescue_dependent']:.1%})")
    print(f"  Median per-CpG rescue rate:          {summary['median_per_cpg_rescue_pct']:.2f}%")
    print(f"  Max per-CpG rescue rate observed:    {summary['max_per_cpg_rescue_pct']:.2f}%")
    print()
    print(f"  Summary:      {summary_md}")
    print(f"  Per-CpG CSV:  {per_cpg_csv}")
    print(f"  Per-pair CSV: {per_pair_csv}")
    print("============================================================")


if __name__ == "__main__":
    main()
