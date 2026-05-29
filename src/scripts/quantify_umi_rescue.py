#!/usr/bin/env python3
"""
Quantify the UMI 'rescue' that coordinate-only dedup would have lost.

Reads every umi_tools dedup log (*.dedup.stats) under a directory and extracts:
  - Input read pairs
  - Output reads (after UMI dedup)
  - Total positions deduplicated  (== output of hypothetical coord-only dedup)
  - Mean unique UMIs per position
  - Max  unique UMIs per position

From these, two derived metrics that quantify UMI value WITHOUT re-running anything:
  - rescue_reads      = output_reads - n_positions               (reads UMI kept that coord-only would discard)
  - rescue_pct        = rescue_reads / output_reads * 100        (% of UMI output owed to UMI-awareness)
  - dup_rate_umi      = 1 - output_reads / input_reads
  - dup_rate_coord    = 1 - n_positions  / input_reads           (what coord-only would have called duplicates)
  - false_kill_pct    = (n_positions - output_reads) / output_reads * 100   (== -rescue_pct, kept for clarity)

Writes a CSV summary and prints a one-line aggregate.

Usage:
  python src/scripts/quantify_umi_rescue.py temp/20260328_075240 temp/20260403 \
      --out docs/figures/umi_rescue_summary.csv
"""

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path


PAT_INPUT     = re.compile(r"Input Reads:\s*(\d+)")
PAT_OUT       = re.compile(r"Number of reads out:\s*(\d+)")
PAT_POS       = re.compile(r"Total number of positions deduplicated:\s*(\d+)")
PAT_MEAN      = re.compile(r"Mean number of unique UMIs per position:\s*([\d.]+)")
PAT_MAX       = re.compile(r"Max\. number of unique UMIs per position:\s*(\d+)")


def parse_log(path: str) -> dict | None:
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
        "sample":        sample,
        "input_pairs":   n_in,
        "umi_out_reads": n_out,
        "n_positions":   n_pos,
        "mean_umi_per_pos": float(m_mean.group(1)) if m_mean else None,
        "max_umi_per_pos":  int(m_max.group(1))    if m_max  else None,
        # derived
        "dup_rate_umi_pct":    round((1 - n_out / n_in) * 100, 2),
        "dup_rate_coord_pct":  round((1 - n_pos / n_in) * 100, 2),
        "rescue_reads":        n_out - n_pos,
        "rescue_pct":          round((n_out - n_pos) / n_out * 100, 3),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+", help="Directories to scan for *.dedup.stats")
    ap.add_argument("--out", default="docs/figures/umi_rescue_summary.csv")
    args = ap.parse_args()

    rows: list = []
    seen: set = set()
    for d in args.dirs:
        for p in sorted(glob.glob(os.path.join(d, "*.dedup.stats"))):
            row = parse_log(p)
            if not row:
                continue
            key = row["sample"]
            if key in seen:
                continue   # keep first occurrence
            seen.add(key)
            rows.append(row)

    if not rows:
        print("no parseable .dedup.stats found", file=sys.stderr)
        return 1

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fields = ["sample", "input_pairs", "umi_out_reads", "n_positions",
              "mean_umi_per_pos", "max_umi_per_pos",
              "dup_rate_umi_pct", "dup_rate_coord_pct",
              "rescue_reads", "rescue_pct"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # Aggregate
    tot_in   = sum(r["input_pairs"]   for r in rows)
    tot_out  = sum(r["umi_out_reads"] for r in rows)
    tot_pos  = sum(r["n_positions"]   for r in rows)
    tot_resc = sum(r["rescue_reads"]  for r in rows)
    mean_umi = sum((r["mean_umi_per_pos"] or 0) for r in rows) / len(rows)
    max_max  = max((r["max_umi_per_pos"] or 0) for r in rows)

    print(f"samples parsed:           {len(rows)}")
    print(f"aggregate input pairs:    {tot_in:,}")
    print(f"aggregate UMI-output:     {tot_out:,}  ({(1-tot_out/tot_in)*100:.2f}% dup rate)")
    print(f"aggregate coord-output:   {tot_pos:,}  ({(1-tot_pos/tot_in)*100:.2f}% dup rate)")
    print(f"aggregate rescued reads:  {tot_resc:,}  ({tot_resc/tot_out*100:.3f}% of UMI output)")
    print(f"mean of per-sample 'mean UMIs per position': {mean_umi:.3f}")
    print(f"max  of per-sample 'max  UMIs per position': {max_max}")
    print(f"\nCSV written to: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
