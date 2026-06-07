#!/usr/bin/env python3
"""
Quantify and visualize UMI rescue from ``umi_tools dedup`` logs.

Workflow
--------
1. Scan one or more directories for ``*.dedup.stats`` log files.
2. Parse each log for input pairs, UMI-output reads,
   coord-only-equivalent positions, and per-position UMI cardinality.
3. Compute per-sample rescue metrics:
       rescue_reads = umi_out_reads − n_positions
       rescue_pct   = rescue_reads / umi_out_reads × 100
4. Write a CSV summary and produce a two-panel diagnostic figure.

The figure has two panels:

  (a) Scatter of rescue % vs input pairs (log x-axis), coloured by
      mean unique UMIs per position. The median rescue rate is shown
      as a horizontal reference line and the per-cohort statistics
      are rendered in a corner text box (no floating annotations that
      overlap the data).

  (b) Scatter of UMI-aware dup rate vs coord-only dup rate — one point
      per sample. Points sit *below* the y = x identity line by an
      amount equal to the rescue magnitude (pp). Scales gracefully to
      any cohort size (handles 19 *or* 119+ samples without per-sample
      bars).

Outputs
-------
- ``--out-csv`` (default: ``docs/figures/umi_rescue_summary.csv``)
- ``--out-png`` (default: ``docs/figures/umi_rescue_vs_depth.png``)

Cohort filtering
----------------
By default the script counts every library that has a ``*.dedup.stats``
log, including borderline samples and any QC failures. To restrict the
analysis to the benign + malignant cohort that drives the downstream
classification experiment, pass ``--sample-info <excel>`` (the project
sample info Excel that lists ``fq_R1`` and ``category`` columns); the
script will keep only samples whose ``category`` matches ``--categories``
(default ``benign,malignant``). If ``--sample-info`` is omitted, the
script falls back to ``[data].sample_info_excel`` from ``config.toml``
(via ``--config``); if neither is available no filtering is applied.

Usage
-----

  # Extract + plot, all 146 libraries (no filter)
  python src/scripts/quantify_umi_rescue.py temp/<workdir>

  # Filter to 120 benign + malignant samples via the project sample Excel
  python src/scripts/quantify_umi_rescue.py temp/<workdir> \\
      --sample-info data/ovarian_fq_series_samples_v4.xlsx

  # Plot only (re-generate PNG from existing CSV; filter still respected)
  python src/scripts/quantify_umi_rescue.py --plot-only \\
      --sample-info data/ovarian_fq_series_samples_v4.xlsx

  # CSV only (skip the plot)
  python src/scripts/quantify_umi_rescue.py temp/<workdir> --no-plot
"""

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────
#  Part 1 — Parse umi_tools dedup stats logs
# ─────────────────────────────────────────────────────────────────────────

PAT_INPUT = re.compile(r"Input Reads:\s*(\d+)")
PAT_OUT   = re.compile(r"Number of reads out:\s*(\d+)")
PAT_POS   = re.compile(r"Total number of positions deduplicated:\s*(\d+)")
PAT_MEAN  = re.compile(r"Mean number of unique UMIs per position:\s*([\d.]+)")
PAT_MAX   = re.compile(r"Max\. number of unique UMIs per position:\s*(\d+)")


def parse_log(path):
    """Parse one ``*.dedup.stats`` file. Returns a per-sample dict or None."""
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    m_in   = PAT_INPUT.search(text)
    m_out  = PAT_OUT.search(text)
    m_pos  = PAT_POS.search(text)
    m_mean = PAT_MEAN.search(text)
    m_max  = PAT_MAX.search(text)
    if not (m_in and m_out and m_pos):
        return None
    sample = Path(path).name.split(".")[0]
    n_in   = int(m_in.group(1))
    n_out  = int(m_out.group(1))
    n_pos  = int(m_pos.group(1))
    return {
        "sample":             sample,
        "input_pairs":        n_in,
        "umi_out_reads":      n_out,
        "n_positions":        n_pos,
        "mean_umi_per_pos":   float(m_mean.group(1)) if m_mean else None,
        "max_umi_per_pos":    int(m_max.group(1))    if m_max  else None,
        "dup_rate_umi_pct":   round((1 - n_out / n_in) * 100, 2),
        "dup_rate_coord_pct": round((1 - n_pos / n_in) * 100, 2),
        "rescue_reads":       n_out - n_pos,
        "rescue_pct":         round((n_out - n_pos) / n_out * 100, 3),
    }


def extract(dirs):
    """Scan ``dirs`` for ``*.dedup.stats`` and parse each file once."""
    rows, seen = [], set()
    for d in dirs:
        for p in sorted(glob.glob(os.path.join(d, "*.dedup.stats"))):
            row = parse_log(p)
            if not row:
                continue
            if row["sample"] in seen:
                continue
            seen.add(row["sample"])
            rows.append(row)
    return rows


# ─────────────────────────────────────────────────────────────────────────
#  Part 1.5 — Sample-cohort filtering (benign + malignant only)
# ─────────────────────────────────────────────────────────────────────────

# Stripped suffixes when converting an Excel ``fq_R1`` value to a sample ID.
_FQ_SUFFIX_RE = re.compile(r"\.R[12]\.(fq|fastq)(\.gz)?$", re.IGNORECASE)


def _sample_id_from_fq(fq_filename):
    """Convert ``B10001.R1.fq.gz`` → ``B10001`` (matches ``parse_log``'s logic)."""
    base = os.path.basename(str(fq_filename))
    stripped = _FQ_SUFFIX_RE.sub("", base)
    # Final guard: take part before the first dot (matches parse_log)
    return stripped.split(".")[0]


def load_included_samples(excel_path, categories):
    """Return the set of sample IDs from the Excel that match ``categories``.

    Reads every sheet of the project sample-info Excel and keeps rows whose
    ``category`` column (case-insensitive, whitespace-stripped) matches one
    of ``categories``. Sample IDs are derived from the ``fq_R1`` column by
    stripping the ``.R[12].fq[.gz]`` suffix, mirroring the prefix that
    ``parse_log()`` extracts from ``*.dedup.stats`` filenames.

    Returns an empty set if the Excel cannot be parsed or has no matching
    rows.
    """
    try:
        import pandas as pd  # type: ignore
    except ImportError:
        print("ERROR: pandas required for --sample-info filtering; "
              "install with `pip install pandas openpyxl`", file=sys.stderr)
        return set()

    if not os.path.exists(excel_path):
        print(f"WARN: --sample-info file not found: {excel_path}",
              file=sys.stderr)
        return set()

    cats = {c.strip().lower() for c in categories if c and c.strip()}
    samples = set()
    try:
        xl = pd.ExcelFile(excel_path)
    except Exception as exc:
        print(f"WARN: failed to open {excel_path}: {exc}", file=sys.stderr)
        return set()

    for sheet in xl.sheet_names:
        try:
            df = xl.parse(sheet)
        except Exception:
            continue
        if "category" not in df.columns or "fq_R1" not in df.columns:
            continue
        cat_series = df["category"].astype(str).str.strip().str.lower()
        mask = cat_series.isin(cats)
        for fq in df.loc[mask, "fq_R1"].dropna():
            sid = _sample_id_from_fq(fq)
            if sid:
                samples.add(sid)
    return samples


def _read_sample_excel_from_config(config_path):
    """Look up ``[data].sample_info_excel`` from ``config_path``; None if absent.

    Resolves the path the same way the pipeline does: bare basename in the
    config gets joined with ``[data].data_dir`` (default ``data``). Returns the
    first candidate that actually exists, else the joined path (so the caller's
    "not found" warning surfaces a useful path).
    """
    if not config_path or not os.path.exists(config_path):
        return None
    try:
        if sys.version_info >= (3, 11):
            import tomllib  # type: ignore
            with open(config_path, "rb") as f:
                cfg = tomllib.load(f)
        else:
            import toml  # type: ignore
            cfg = toml.load(config_path)
    except Exception:
        return None

    data = cfg.get("data", {}) or {}
    excel = data.get("sample_info_excel")
    if not excel:
        return None

    candidates = [excel]
    data_dir = data.get("data_dir") or "data"
    if not os.path.isabs(excel):
        candidates.append(os.path.join(data_dir, excel))

    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[-1]


def apply_sample_filter(rows, excel_path, categories):
    """Filter ``rows`` to samples matching ``categories``. Returns the filtered
    list (or the original list if no Excel was provided or matched)."""
    if not excel_path:
        return rows
    included = load_included_samples(excel_path, categories)
    if not included:
        print(f"WARN: no included samples loaded from {excel_path} "
              f"(category set: {categories}); skipping filter",
              file=sys.stderr)
        return rows

    before = len(rows)
    kept = [r for r in rows if r["sample"] in included]
    matched = {r["sample"] for r in kept}
    missing = included - matched
    print(f"Cohort filter: {len(kept)}/{before} samples retained "
          f"(categories: {','.join(categories)}; "
          f"{len(missing)} Excel-listed samples had no dedup.stats and were skipped)",
          file=sys.stderr)
    return kept


# ─────────────────────────────────────────────────────────────────────────
#  Part 2 — CSV I/O
# ─────────────────────────────────────────────────────────────────────────

CSV_FIELDS = [
    "sample", "input_pairs", "umi_out_reads", "n_positions",
    "mean_umi_per_pos", "max_umi_per_pos",
    "dup_rate_umi_pct", "dup_rate_coord_pct",
    "rescue_reads", "rescue_pct",
]


def write_csv(rows, out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)


def load_csv(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                rows.append({
                    "sample":             r["sample"],
                    "input_pairs":        int(r["input_pairs"]),
                    "umi_out_reads":      int(r["umi_out_reads"]),
                    "n_positions":        int(r["n_positions"]),
                    "mean_umi_per_pos":   float(r["mean_umi_per_pos"]) if r.get("mean_umi_per_pos") else 0.0,
                    "max_umi_per_pos":    int(r["max_umi_per_pos"])    if r.get("max_umi_per_pos")  else 0,
                    "dup_rate_umi_pct":   float(r["dup_rate_umi_pct"]),
                    "dup_rate_coord_pct": float(r["dup_rate_coord_pct"]),
                    "rescue_reads":       int(r["rescue_reads"]),
                    "rescue_pct":         float(r["rescue_pct"]),
                })
            except (KeyError, ValueError):
                continue
    rows.sort(key=lambda x: x["input_pairs"])
    return rows


def print_aggregate(rows):
    tot_in   = sum(r["input_pairs"]   for r in rows)
    tot_out  = sum(r["umi_out_reads"] for r in rows)
    tot_pos  = sum(r["n_positions"]   for r in rows)
    tot_resc = sum(r["rescue_reads"]  for r in rows)
    mean_umi = sum((r["mean_umi_per_pos"] or 0) for r in rows) / len(rows)
    max_max  = max((r["max_umi_per_pos"] or 0) for r in rows)
    per_sample = sorted(r["rescue_pct"] for r in rows)
    median = per_sample[len(per_sample) // 2]
    lo, hi = per_sample[0], per_sample[-1]

    print(f"samples parsed:              {len(rows)}")
    print(f"aggregate input pairs:       {tot_in:,}")
    print(f"aggregate UMI-output:        {tot_out:,}  ({(1-tot_out/tot_in)*100:.2f}% dup rate)")
    print(f"aggregate coord-output:      {tot_pos:,}  ({(1-tot_pos/tot_in)*100:.2f}% dup rate)")
    print(f"aggregate rescued reads:     {tot_resc:,}  ({tot_resc/tot_out*100:.3f}% of UMI output)")
    print(f"per-sample rescue:           median {median:.2f}%   range {lo:.2f}%–{hi:.2f}%")
    print(f"mean of per-sample 'mean UMIs per position': {mean_umi:.3f}")
    print(f"max  of per-sample 'max  UMIs per position': {max_max}")


# ─────────────────────────────────────────────────────────────────────────
#  Part 3 — Plot
# ─────────────────────────────────────────────────────────────────────────

def make_plot(rows, out_path):
    """Two-panel diagnostic figure of UMI rescue (see module docstring)."""
    n = len(rows)
    x_depth = np.array([r["input_pairs"] / 1e6      for r in rows])
    y_resc  = np.array([r["rescue_pct"]              for r in rows])
    c_umi   = np.array([r["mean_umi_per_pos"] or 0.0 for r in rows])
    dup_umi = np.array([r["dup_rate_umi_pct"]        for r in rows])
    dup_crd = np.array([r["dup_rate_coord_pct"]      for r in rows])

    median_rescue = float(np.median(y_resc))
    diff          = dup_crd - dup_umi   # pp below the diagonal (per sample)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    # ─── Panel A: rescue % vs sequencing depth ───
    sc1 = ax1.scatter(x_depth, y_resc,
                      c=c_umi, s=70, cmap="viridis",
                      edgecolor="black", linewidth=0.5, alpha=0.85, zorder=3)
    cb1 = plt.colorbar(sc1, ax=ax1)
    cb1.set_label("Mean unique UMIs per position", fontsize=10)

    ax1.set_xscale("log")
    ax1.set_xlabel("Input read pairs (millions, log scale)", fontsize=11)
    ax1.set_ylabel("Reads rescued by UMI awareness (%)", fontsize=11)
    ax1.set_title("(a) Rescue rate vs sequencing depth", fontsize=12)
    ax1.grid(alpha=0.3, zorder=0)

    # Median reference line — clean, no overlapping floating annotations
    ax1.axhline(median_rescue, color="#c1121f", linestyle="--",
                linewidth=1.4, alpha=0.7, zorder=2,
                label=f"Median = {median_rescue:.2f}%")
    ax1.legend(loc="upper right", fontsize=9, framealpha=0.9)

    # Corner stats box (always out of the data region)
    stats_a = (
        f"n = {n} libraries\n"
        f"Median rescue: {median_rescue:.2f}%\n"
        f"Range: {y_resc.min():.2f}% – {y_resc.max():.2f}%"
    )
    ax1.text(0.02, 0.98, stats_a,
             transform=ax1.transAxes, fontsize=9,
             verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.4",
                       facecolor="white", edgecolor="#666", alpha=0.92))

    # ─── Panel B: dup_umi vs dup_coord scatter ───
    sc2 = ax2.scatter(dup_crd, dup_umi,
                      c=c_umi, s=70, cmap="viridis",
                      edgecolor="black", linewidth=0.5, alpha=0.85, zorder=3)
    cb2 = plt.colorbar(sc2, ax=ax2)
    cb2.set_label("Mean unique UMIs per position", fontsize=10)

    # y = x identity line — points below it = UMI rescued reads coord-only would discard
    lo = min(dup_umi.min(), dup_crd.min()) - 1
    hi = max(dup_umi.max(), dup_crd.max()) + 1
    ax2.plot([lo, hi], [lo, hi], color="#444",
             linestyle=":", linewidth=1.2, alpha=0.6, zorder=2,
             label="y = x  (no rescue)")
    ax2.set_xlim(lo, hi)
    ax2.set_ylim(lo, hi)

    ax2.set_xlabel("Coord-only duplicate rate (%)", fontsize=11)
    ax2.set_ylabel("UMI-aware duplicate rate (%)",  fontsize=11)
    ax2.set_title("(b) Coord-only would over-collapse — distance below diagonal = rescue",
                  fontsize=12)
    ax2.grid(alpha=0.3, zorder=0)
    ax2.legend(loc="lower right", fontsize=9, framealpha=0.9)

    stats_b = (
        f"n = {n} libraries\n"
        f"Median Δ (coord − UMI): {float(np.median(diff)):.2f} pp\n"
        f"Range Δ: {diff.min():.2f} – {diff.max():.2f} pp"
    )
    ax2.text(0.02, 0.98, stats_b,
             transform=ax2.transAxes, fontsize=9,
             verticalalignment="top",
             bbox=dict(boxstyle="round,pad=0.4",
                       facecolor="white", edgecolor="#666", alpha=0.92))

    fig.suptitle(
        f"UMI-aware deduplication rescues a median of {median_rescue:.2f}% "
        f"of cfDNA reads that coord-only deduplication would discard "
        f"(n = {n} libraries)",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────
#  Part 4 — CLI dispatch
# ─────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Quantify and visualize UMI rescue from umi_tools dedup logs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("dirs", nargs="*",
                    help="Directories to scan for *.dedup.stats "
                         "(omit when using --plot-only)")
    ap.add_argument("--out-csv", default="docs/figures/umi_rescue_summary.csv",
                    help="Per-sample CSV output")
    ap.add_argument("--out-png", default="docs/figures/umi_rescue_vs_depth.png",
                    help="Two-panel diagnostic PNG output")
    ap.add_argument("--csv", default=None,
                    help="Existing CSV to read for --plot-only "
                         "(default: same as --out-csv)")
    ap.add_argument("--no-plot",   action="store_true",
                    help="Skip the PNG plot (CSV only)")
    ap.add_argument("--plot-only", action="store_true",
                    help="Skip extraction; only re-plot from an existing CSV")
    ap.add_argument("--sample-info", default=None,
                    help="Path to the project sample info Excel "
                         "(columns 'fq_R1' + 'category'); when given, restrict "
                         "the analysis to samples whose category matches "
                         "--categories")
    ap.add_argument("--categories", default="benign,malignant",
                    help="Comma-separated category names to keep when "
                         "--sample-info is supplied (default: benign,malignant)")
    ap.add_argument("--config", default="config.toml",
                    help="Fallback config.toml to read [data].sample_info_excel "
                         "if --sample-info is not given (default: config.toml)")
    args = ap.parse_args()

    # Resolve cohort filter: explicit flag wins, otherwise fall back to config.toml
    sample_excel = args.sample_info or _read_sample_excel_from_config(args.config)
    categories   = [c.strip() for c in args.categories.split(",") if c.strip()]

    if args.plot_only:
        csv_path = args.csv or args.out_csv
        if not os.path.exists(csv_path):
            print(f"ERROR: {csv_path} not found (need CSV for --plot-only)",
                  file=sys.stderr)
            return 1
        rows = load_csv(csv_path)
        if not rows:
            print(f"ERROR: no rows parsed from {csv_path}", file=sys.stderr)
            return 1
        rows = apply_sample_filter(rows, sample_excel, categories)
        if not rows:
            print("ERROR: cohort filter removed every row", file=sys.stderr)
            return 1
        make_plot(rows, args.out_png)
        print(f"Plot saved → {args.out_png}")
        return 0

    if not args.dirs:
        print("ERROR: need at least one directory to scan, or --plot-only",
              file=sys.stderr)
        return 1

    rows = extract(args.dirs)
    if not rows:
        print("ERROR: no parseable *.dedup.stats files found", file=sys.stderr)
        return 1

    rows = apply_sample_filter(rows, sample_excel, categories)
    if not rows:
        print("ERROR: cohort filter removed every row", file=sys.stderr)
        return 1

    write_csv(rows, args.out_csv)
    print(f"CSV saved  → {args.out_csv}")
    print_aggregate(rows)

    if not args.no_plot:
        make_plot(rows, args.out_png)
        print(f"Plot saved → {args.out_png}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
