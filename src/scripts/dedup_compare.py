#!/usr/bin/env python3
"""
UMI-aware vs coordinate-only deduplication: comparison figures for the manuscript.

Generates (in --output-dir, default = docs/figures):
  - dedup_compare_fragment_size.png : per-sample fragment-size overlay (pre / coord / umi)
  - dedup_compare_summary.csv       : retention statistics per (sample, strategy)
  - dedup_compare_beta_scatter.png  : panel-CpG beta correlation (optional)
  - dedup_compare_beta_stats.csv    : per-sample Pearson r / RMSE on the panel (optional)

The beta-correlation section runs only if [data].probe_panel_bed and
[extraction].reference_genome resolve to existing files (typically needs WSL/Docker).

Usage:
  python src/scripts/dedup_compare.py \
      temp/20260327_134610/B10002.bwameth_pe.sorted.bam \
      temp/20260327_134610/B10197.bwameth_pe.sorted.bam \
      --output-dir docs/figures
"""

import argparse
import csv
import logging
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import toml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logger = logging.getLogger("dedup_compare")

SIZE_BIN_MAX = 600
PEAK_BROAD   = (100, 220)
PEAK_NARROW  = (147, 167)


def _resolve_tool(name: str, config: dict) -> str:
    bin_dir = (config.get("tools") or {}).get("conda_env_bin")
    if bin_dir and os.path.exists(os.path.join(bin_dir, name)):
        return os.path.join(bin_dir, name)
    return name


def _umi_tools(config: dict) -> str:
    bin_dir = (config.get("tools") or {}).get("conda_env_bin", "")
    primary = os.path.join(bin_dir, "umi_tools") if bin_dir else "umi_tools"
    if os.path.exists(primary):
        return primary
    fallback = (config.get("tools") or {}).get("umi_tools_fallback")
    if fallback and os.path.exists(fallback):
        return fallback
    return "umi_tools"


def run_dedup(sorted_bam: str, out_bam: str, ignore_umi: bool, config: dict) -> str:
    if os.path.exists(out_bam) and os.path.exists(out_bam + ".bai"):
        logger.info(f"    reuse: {out_bam}")
        return out_bam
    log_path = out_bam.replace(".bam", ".log")
    cmd = [_umi_tools(config), "dedup", "--paired",
           "--stdin", sorted_bam, "--stdout", out_bam, "--log", log_path]
    if ignore_umi:
        cmd.append("--ignore-umi")
    logger.info(f"    run:   {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    subprocess.run([_resolve_tool("samtools", config), "index", out_bam], check=True)
    return out_bam


def extract_tlen(bam: str, samtools: str, max_reads: int = 400_000) -> np.ndarray:
    """Read TLEN of properly-paired primary alignments, positive end only."""
    cmd = [samtools, "view", "-f", "0x2", "-F", "0x904", bam]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    sizes: list = []
    assert proc.stdout is not None
    for line in proc.stdout:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            tlen = int(parts[8])
        except ValueError:
            continue
        if tlen > 0:
            sizes.append(tlen)
            if len(sizes) >= max_reads:
                proc.stdout.close()
                proc.wait()
                break
    else:
        proc.wait()
    return np.asarray(sizes, dtype=np.int32)


def summarise(arr: np.ndarray, sample: str, strategy: str) -> dict:
    arr_in = arr[(arr >= 0) & (arr <= SIZE_BIN_MAX)]
    n = len(arr_in)
    if n == 0:
        return {"sample": sample, "strategy": strategy, "n_reads": 0,
                "median_bp": 0, "pct_100_220": 0.0, "pct_147_167": 0.0}
    pct_broad  = float(((arr_in >= PEAK_BROAD[0])  & (arr_in <= PEAK_BROAD[1])).sum())  / n * 100
    pct_narrow = float(((arr_in >= PEAK_NARROW[0]) & (arr_in <= PEAK_NARROW[1])).sum()) / n * 100
    return {"sample": sample, "strategy": strategy, "n_reads": n,
            "median_bp": int(np.median(arr_in)),
            "pct_100_220": round(pct_broad, 2),
            "pct_147_167": round(pct_narrow, 2)}


def plot_fragment_overlay(per_sample: dict, out_png: str) -> None:
    n = len(per_sample)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 4.5), sharey=True)
    if n == 1:
        axes = [axes]
    bins = np.arange(0, SIZE_BIN_MAX + 2, 2)
    centres = (bins[:-1] + bins[1:]) / 2
    colours = {"pre-dedup": "#7f8c8d", "coord-only": "#e63946", "umi-aware": "#2a9d8f"}
    legend_done = False
    for ax, (sample, sizes) in zip(axes, per_sample.items()):
        for label in ("pre-dedup", "coord-only", "umi-aware"):
            arr = sizes[label]
            arr = arr[(arr >= 0) & (arr <= SIZE_BIN_MAX)]
            if len(arr) == 0:
                continue
            hist, _ = np.histogram(arr, bins=bins, density=True)
            ax.plot(centres, hist, color=colours[label], linewidth=1.8,
                    label=label if not legend_done else None)
        ax.axvspan(PEAK_BROAD[0], PEAK_BROAD[1], color="#f4a261", alpha=0.10,
                   label="100-220 bp cfDNA zone" if not legend_done else None)
        ax.axvline(167, color="#f4a261", linestyle="--", linewidth=0.8,
                   label="167 bp mono-nucleosome" if not legend_done else None)
        ax.set_xlim(0, 400)
        ax.set_xlabel("Fragment size (bp)")
        ax.set_title(sample)
        if ax is axes[0]:
            ax.set_ylabel("Density")
        legend_done = True
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=9, loc="upper right", bbox_to_anchor=(0.99, 0.97))
    fig.suptitle("Effect of deduplication strategy on cfDNA fragment-size distribution",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def extract_panel_beta(bam: str, ref_fa: str, panel_bed: str, out_dir: str,
                       config: dict, tag: str) -> dict:
    """Run MethylDackel extract on a BAM restricted to panel regions, return {(chrom, pos): beta}."""
    methyldackel = _resolve_tool("MethylDackel", config)
    samtools     = _resolve_tool("samtools", config)
    bedtools     = _resolve_tool("bedtools", config)

    sample = Path(bam).name.split(".")[0]
    target_bam = os.path.join(out_dir, f"{sample}.{tag}.targeted.bam")
    sorted_bam = os.path.join(out_dir, f"{sample}.{tag}.targeted.sorted.bam")
    if not (os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai")):
        with open(target_bam, "wb") as fout:
            subprocess.run([bedtools, "intersect", "-a", bam, "-b", panel_bed, "-wa", "-u"],
                           stdout=fout, check=True)
        subprocess.run([samtools, "sort", "-o", sorted_bam, target_bam], check=True)
        subprocess.run([samtools, "index", sorted_bam], check=True)
        os.remove(target_bam)

    prefix = os.path.join(out_dir, f"{sample}.{tag}")
    bg_path = f"{prefix}_CpG.bedGraph"
    if not os.path.exists(bg_path):
        subprocess.run([methyldackel, "extract", "-q", "20", "-p", "5",
                        "-o", prefix, ref_fa, sorted_bam], check=True)

    beta = {}
    with open(bg_path) as f:
        for line in f:
            if line.startswith("track"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, start, _, pct, meth, unmeth = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
            depth = int(meth) + int(unmeth)
            if depth < 10:
                continue
            beta[(chrom, int(start))] = float(pct) / 100.0
    return beta


def plot_beta_correlation(samples: list, workdir: str, out_dir: str,
                          panel_bed: str, ref_fa: str, config: dict) -> None:
    """Compute panel-CpG β for both dedup BAMs, plot scatter + Pearson r."""
    n = len(samples)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), sharex=True, sharey=True)
    if n == 1:
        axes = [axes]
    stats_rows: list = []
    for ax, sample in zip(axes, samples):
        coord_bam = os.path.join(workdir, f"{sample}.coord_dedup.bam")
        umi_bam   = os.path.join(workdir, f"{sample}.umi_dedup.bam")
        logger.info(f"  • β extract: {sample} (coord-only)")
        beta_c = extract_panel_beta(coord_bam, ref_fa, panel_bed, workdir, config, "coord")
        logger.info(f"  • β extract: {sample} (umi-aware)")
        beta_u = extract_panel_beta(umi_bam,   ref_fa, panel_bed, workdir, config, "umi")
        common = sorted(set(beta_c) & set(beta_u))
        if not common:
            ax.set_title(f"{sample} (no shared CpG)")
            continue
        x = np.array([beta_c[k] for k in common])
        y = np.array([beta_u[k] for k in common])
        r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else float("nan")
        rmse = float(np.sqrt(np.mean((x - y) ** 2)))
        ax.scatter(x, y, s=8, alpha=0.4, color="#2a9d8f")
        ax.plot([0, 1], [0, 1], "--", color="#7f8c8d", linewidth=0.8)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("β (coordinate-only dedup)")
        if ax is axes[0]:
            ax.set_ylabel("β (UMI-aware dedup)")
        ax.set_title(f"{sample}\nn={len(common)} CpG  r={r:.3f}  RMSE={rmse:.3f}")
        stats_rows.append({"sample": sample, "n_cpg": len(common),
                           "pearson_r": round(r, 4), "rmse_beta": round(rmse, 4)})
    out_png = os.path.join(out_dir, "dedup_compare_beta_scatter.png")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    logger.info(f"β scatter saved → {out_png}")

    out_csv = os.path.join(out_dir, "dedup_compare_beta_stats.csv")
    if stats_rows:
        with open(out_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(stats_rows[0].keys()))
            w.writeheader()
            w.writerows(stats_rows)
        logger.info(f"β stats saved  → {out_csv}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bams", nargs="+", help="One or more pre-dedup *.sorted.bam files")
    ap.add_argument("--config",     default="config.toml")
    ap.add_argument("--output-dir", default="docs/figures")
    ap.add_argument("--workdir",    default="temp/dedup_compare",
                    help="Scratch dir for the two dedup BAMs (reused across runs)")
    ap.add_argument("--no-beta", action="store_true",
                    help="Skip the β-value correlation section even if reference + panel are accessible")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if not os.path.exists(args.config):
        logger.error(f"config not found: {args.config}")
        return 2
    config = toml.load(args.config)
    samtools = _resolve_tool("samtools", config)
    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    rows: list = []
    per_sample: dict = {}

    for sorted_bam in args.bams:
        if not os.path.exists(sorted_bam):
            logger.warning(f"skip missing BAM: {sorted_bam}")
            continue
        sample = Path(sorted_bam).name.split(".")[0]
        logger.info(f"━━ {sample} ━━")
        coord_bam = os.path.join(args.workdir, f"{sample}.coord_dedup.bam")
        umi_bam   = os.path.join(args.workdir, f"{sample}.umi_dedup.bam")
        logger.info("  • coord-only dedup (--ignore-umi)")
        run_dedup(sorted_bam, coord_bam, ignore_umi=True,  config=config)
        logger.info("  • UMI-aware dedup")
        run_dedup(sorted_bam, umi_bam,   ignore_umi=False, config=config)
        logger.info("  • TLEN extraction")
        pre   = extract_tlen(sorted_bam, samtools)
        coord = extract_tlen(coord_bam,  samtools)
        umi   = extract_tlen(umi_bam,    samtools)
        per_sample[sample] = {"pre-dedup": pre, "coord-only": coord, "umi-aware": umi}
        for arr, label in ((pre, "pre-dedup"), (coord, "coord-only"), (umi, "umi-aware")):
            row = summarise(arr, sample, label)
            rows.append(row)
            logger.info(f"    {label:11s} n={row['n_reads']:7d}  "
                        f"median={row['median_bp']:3d}bp  "
                        f"100-220bp={row['pct_100_220']:5.1f}%  "
                        f"147-167bp={row['pct_147_167']:5.1f}%")

    if not per_sample:
        logger.error("no usable BAMs — nothing produced")
        return 1

    out_png = os.path.join(args.output_dir, "dedup_compare_fragment_size.png")
    plot_fragment_overlay(per_sample, out_png)
    logger.info(f"fragment overlay → {out_png}")

    out_csv = os.path.join(args.output_dir, "dedup_compare_summary.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample", "strategy", "n_reads",
                                          "median_bp", "pct_100_220", "pct_147_167"])
        w.writeheader()
        w.writerows(rows)
    logger.info(f"summary CSV      → {out_csv}")

    # β-value correlation (optional)
    if args.no_beta:
        logger.info("β-value correlation skipped (--no-beta)")
        return 0
    panel_bed = (config.get("data") or {}).get("probe_panel_bed")
    ref_fa    = (config.get("extraction") or {}).get("reference_genome")
    if not (panel_bed and os.path.exists(panel_bed)):
        logger.info(f"β-value correlation skipped (panel BED not accessible: {panel_bed!r})")
        return 0
    if not (ref_fa and os.path.exists(ref_fa)):
        logger.info(f"β-value correlation skipped (reference genome not accessible: {ref_fa!r})")
        return 0
    logger.info("━━ β-value correlation (panel CpGs) ━━")
    try:
        plot_beta_correlation(list(per_sample.keys()), args.workdir, args.output_dir,
                              panel_bed, ref_fa, config)
    except subprocess.CalledProcessError as exc:
        logger.error(f"β-value correlation failed: {exc}")
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
