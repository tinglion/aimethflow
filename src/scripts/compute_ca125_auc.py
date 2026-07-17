#!/usr/bin/env python3
"""
compute_ca125_auc.py — standalone CA-125 AUC for benign vs malignant.

Reads (category, ca125) pairs directly from the sample info Excel — no pipeline
workdir, no classifier_predictions.csv, no model training required. CA-125 is a
serum biomarker, so AUC is computed directly from the raw values vs the
binary clinical label.

Usage
-----
    # default: read sample_info_excel from config.toml
    python src/scripts/compute_ca125_auc.py

    # explicit Excel + sheet
    python src/scripts/compute_ca125_auc.py --excel data/ovarian_fq_series_samples_v4.xlsx --sheet plasma

    # save result to a different CSV
    python src/scripts/compute_ca125_auc.py --out result/standalone_ca125_auc.csv

Outputs
-------
    Console: AUC, 95% bootstrap CI, n_benign, n_malignant
    CSV    : single-row table with the same numbers (for record-keeping)
"""

import argparse
import csv
import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


# ── Column-name aliases (case-insensitive, leading/trailing whitespace stripped) ─
_CA125_ALIASES    = {"ca125", "ca-125", "ca_125", "ca 125", "ca125(u/ml)", "ca-125(u/ml)"}
_CATEGORY_ALIASES = {"category", "group", "label"}
_CASE_ID_ALIASES  = {"case_id", "caseid", "patient_id"}
_FQ_R1_ALIASES    = {"fq_r1", "fastq_r1", "r1", "fq_r1_path"}


def _sid_from_fq_r1(raw) -> str:
    """Derive sample_id from a fq_R1 filename: 'B10013_R1.fq.gz' → 'B10013'."""
    s = str(raw).strip() if raw is not None else ""
    if not s or s.lower() == "nan":
        return ""
    return s.split(".")[0].split("_")[0]


def _resolve_col(df: pd.DataFrame, aliases: set) -> Optional[str]:
    """Return the original column name in df whose lowercase form matches any alias."""
    col_map = {c.strip().lower(): c for c in df.columns}
    return next((col_map[a] for a in aliases if a in col_map), None)


def _load_excel_path(arg_excel: Optional[str]) -> str:
    """Resolve the sample-info Excel path: --excel arg → config.toml → conventional defaults."""
    if arg_excel:
        return arg_excel
    try:
        # Project's load_config; falls back to config.toml in repo root
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
        from utils.config import load_config  # type: ignore
        cfg = load_config()
        from methylation.sample_grouping import _resolve_excel  # type: ignore
        return _resolve_excel(cfg.get("data", {}).get("sample_info_excel", "ovarian_fq_series_samples_v4.xlsx"))
    except Exception:
        # Last-resort literal lookup
        for p in ["data/ovarian_fq_series_samples_v4.xlsx",
                 "ovarian_fq_series_samples_v4.xlsx"]:
            if os.path.exists(p):
                return p
        return "ovarian_fq_series_samples_v4.xlsx"


def _bootstrap_auc_ci(y_true, y_score, n_boot: int = 1000, ci: float = 0.95, seed: int = 42):
    """Threshold-free percentile bootstrap CI for AUC."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    alpha = (1 - ci) / 2
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        aucs.append(roc_auc_score(yt, ys))
    if not aucs:
        return float("nan"), float("nan")
    a = np.array(aucs)
    return (float(np.percentile(a, 100 * alpha)),
            float(np.percentile(a, 100 * (1 - alpha))))


def collect_ca125_pairs(excel_path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    """
    Walk the Excel and return rows (case_id, category, ca125) for all sheets
    matching `sheet` (default: every sheet that has ca125 + category columns).
    """
    if not os.path.exists(excel_path):
        sys.exit(f"[error] Excel not found: {excel_path}")

    xl = pd.ExcelFile(excel_path)
    out = []
    for sheet_name in xl.sheet_names:
        if sheet and sheet_name != sheet:
            continue
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
        ca125_col   = _resolve_col(df, _CA125_ALIASES)
        cat_col     = _resolve_col(df, _CATEGORY_ALIASES)
        case_id_col = _resolve_col(df, _CASE_ID_ALIASES)
        fq_r1_col   = _resolve_col(df, _FQ_R1_ALIASES)
        if ca125_col is None:
            print(f"[skip] sheet '{sheet_name}': no CA-125 column", file=sys.stderr)
            continue
        if cat_col is None:
            print(f"[skip] sheet '{sheet_name}': no category column", file=sys.stderr)
            continue
        for _, row in df.iterrows():
            cat = str(row.get(cat_col, "")).strip().lower()
            val = row.get(ca125_col)
            cid = row.get(case_id_col) if case_id_col else ""
            sid = _sid_from_fq_r1(row.get(fq_r1_col)) if fq_r1_col else ""
            if pd.isna(val) or val == "" or cat in ("", "nan"):
                continue
            try:
                v = float(val)
            except (TypeError, ValueError):
                continue
            out.append({"sheet": sheet_name, "case_id": str(cid).strip(),
                        "sample_id": sid, "category": cat, "ca125": v})
    return pd.DataFrame(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("--excel", default=None,
                        help="Path to sample info Excel. Default: from config.toml.")
    parser.add_argument("--sheet", default=None,
                        help="Restrict to a specific sheet (e.g. 'plasma'). Default: all sheets with ca125+category.")
    parser.add_argument("--out", default=None,
                        help="CSV path to save the AUC summary. Default: print only.")
    parser.add_argument("--values-out", default=None,
                        help="Optional CSV path for per-sample (case_id, category, ca125) records.")
    parser.add_argument("--include-borderline", action="store_true",
                        help="Treat borderline-category samples as malignant (default: drop them).")
    parser.add_argument("--n-bootstrap", type=int, default=1000,
                        help="Bootstrap iterations for the AUC 95%% CI (default: 1000).")
    parser.add_argument("--exclude", nargs="*", default=None, metavar="SAMPLE_ID",
                        help="Sample IDs (fq_R1 prefix, e.g. B10013) to exclude. "
                             "Default: ml.exclude_samples from config.toml. "
                             "Pass `--exclude` with no values to disable exclusion.")
    args = parser.parse_args()

    excel_path = _load_excel_path(args.excel)
    print(f"[info] Reading {excel_path}" + (f" (sheet '{args.sheet}')" if args.sheet else ""))

    df = collect_ca125_pairs(excel_path, args.sheet)
    if df.empty:
        sys.exit("[error] No (category, ca125) pairs found. "
                 "Verify the Excel has both a 'ca125' column and a 'category' column.")

    print(f"[info] Loaded {len(df)} rows with CA-125 values across "
          f"{df['sheet'].nunique()} sheet(s): {sorted(df['sheet'].unique())}")

    # Resolve outlier exclusion list
    if args.exclude is None:
        # No --exclude flag → use ml.exclude_samples from config
        try:
            sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
            from utils.config import load_config  # type: ignore
            cfg = load_config()
            exclude_set = set(cfg.get("ml", {}).get("exclude_samples", []) or [])
            if exclude_set:
                print(f"[info] Outlier exclusion (from config.toml ml.exclude_samples): "
                      f"{sorted(exclude_set)}")
        except Exception as exc:
            print(f"[warn] Could not load ml.exclude_samples from config ({exc}) — no exclusion applied.")
            exclude_set = set()
    else:
        # --exclude given (possibly empty)
        exclude_set = set(args.exclude)
        if exclude_set:
            print(f"[info] Outlier exclusion (from --exclude): {sorted(exclude_set)}")
        else:
            print("[info] Outlier exclusion disabled (--exclude with no values)")

    if exclude_set:
        before = len(df)
        present_before = set(df["sample_id"]) & exclude_set
        df = df[~df["sample_id"].isin(exclude_set)].copy()
        n_dropped_outliers = before - len(df)
        print(f"[info] Excluded {n_dropped_outliers} outlier row(s) by sample_id "
              f"(matched: {sorted(present_before)})")
        not_found = exclude_set - present_before
        if not_found:
            print(f"[warn] Exclude IDs not present in this Excel/sheet: {sorted(not_found)}")

    # Map category → binary label
    cat_to_label = {"benign": 0, "control": 0, "healthy": 0,
                    "malignant": 1, "cancer": 1, "tumour": 1, "tumor": 1}
    if args.include_borderline:
        cat_to_label["borderline"] = 1

    df["label"] = df["category"].map(cat_to_label)
    df_drop = df[df["label"].isna()].copy()
    df = df.dropna(subset=["label"]).copy()
    df["label"] = df["label"].astype(int)

    n_benign    = int((df["label"] == 0).sum())
    n_malignant = int((df["label"] == 1).sum())
    n_dropped   = len(df_drop)

    if n_dropped > 0:
        dropped_cats = df_drop["category"].value_counts().to_dict()
        print(f"[info] Dropped {n_dropped} row(s) with non-binary categories: {dropped_cats} "
              "(use --include-borderline to fold borderline → malignant)")

    if n_benign < 2 or n_malignant < 2:
        sys.exit(f"[error] Insufficient data: {n_benign} benign, {n_malignant} malignant — need ≥2 of each.")

    y_true   = df["label"].to_numpy()
    y_score  = df["ca125"].to_numpy()
    df.to_csv("temp/pure_ca125_values.csv", index=False)
    auc      = roc_auc_score(y_true, y_score)
    auc_lo, auc_hi = _bootstrap_auc_ci(y_true, y_score, n_boot=args.n_bootstrap)

    print()
    print(f"  CA-125 AUC : {auc:.4f}")
    print(f"  95% CI     : [{auc_lo:.4f} , {auc_hi:.4f}]   (bootstrap, n={args.n_bootstrap})")
    print(f"  N benign   : {n_benign}")
    print(f"  N malignant: {n_malignant}")
    print(f"  Median     : benign {np.median(y_score[y_true==0]):.1f} U/mL  |  "
          f"malignant {np.median(y_score[y_true==1]):.1f} U/mL")
    print()

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["method", "n_benign", "n_malignant",
                        "auc", "auc_ci95_lo", "auc_ci95_hi",
                        "n_bootstrap", "excel", "sheet"])
            w.writerow(["CA-125", n_benign, n_malignant,
                        round(auc, 4), round(auc_lo, 4), round(auc_hi, 4),
                        args.n_bootstrap, excel_path, args.sheet or "ALL"])
        print(f"[done] AUC summary saved → {args.out}")

    if args.values_out:
        os.makedirs(os.path.dirname(args.values_out) or ".", exist_ok=True)
        df[["sheet", "case_id", "category", "ca125", "label"]].to_csv(
            args.values_out, index=False, encoding="utf-8"
        )
        print(f"[done] Per-sample values saved → {args.values_out}")


if __name__ == "__main__":
    main()
