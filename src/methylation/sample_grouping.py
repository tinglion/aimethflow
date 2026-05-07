#!/usr/bin/env python3
"""
sample_grouping.py - Sample grouping module

Handles sample grouping logic, including index-based grouping and
Excel-file-based grouping using ovarian_fq_series_samples_v3.xlsx.

Excel schema (per sheet):
  fq_R1, fq_R2, age, case_id, pathology, category, subtype[, stage]
  category values: benign | malignant | borderline
"""

import pandas as pd
import os
import re
from utils.logger import logger
from utils.config import load_config

config = load_config()

_DEFAULT_EXCEL = "ovarian_fq_series_samples_v3.xlsx"

# Candidate search paths for the Excel file (relative to project root or CWD)
_EXCEL_SEARCH_DIRS = ["", "data", "data/"]


def _resolve_excel(filename):
    """Return the first existing path for filename, trying data/ prefix if needed."""
    for prefix in _EXCEL_SEARCH_DIRS:
        candidate = os.path.join(prefix, filename) if prefix else filename
        if os.path.exists(candidate):
            return candidate
    # last resort: absolute path as-is
    return filename


def _to_posix_path(path):
    """Convert Windows path (e.g. F:\\meth\\...) to WSL POSIX path (/mnt/f/meth/...)."""
    if not path:
        return path
    m = re.match(r"^([A-Za-z]):[/\\](.*)", path)
    if m:
        drive = m.group(1).lower()
        rest  = m.group(2).replace("\\", "/")
        return f"/mnt/{drive}/{rest}"
    return path.replace("\\", "/")


def group_methylation_files_v1(methylation_files):
    """Split methylation files into two groups by mod-2 index (v1 fallback)."""
    sample1 = [f for i, f in enumerate(methylation_files) if i % 2 == 0]
    sample2 = [f for i, f in enumerate(methylation_files) if i % 2 == 1]
    return sample1, sample2


def group_methylation_files_v2(methylation_files, category="plasma"):
    """
    Split methylation files into (benign, malignant, borderline) using the sample info Excel.
    Falls back to group_methylation_files_v1 on any error (borderline = [] in that case).
    """
    excel_file = _resolve_excel(config.get("data", {}).get("sample_info_excel", _DEFAULT_EXCEL))
    if not os.path.exists(excel_file):
        logger.error(f"Sample info Excel not found: {excel_file}")
        benign_files, malignant_files = group_methylation_files_v1(methylation_files)
        return benign_files, malignant_files, []

    try:
        sample_info = parse_sample_info(excel_file, category)
        if not sample_info:
            logger.error("Failed to parse Excel, falling back to v1 grouping")
            benign_files, malignant_files = group_methylation_files_v1(methylation_files)
            return benign_files, malignant_files, []

        benign_fq     = set(sample_info["benign"]["R1"])
        malignant_fq  = set(sample_info["malignant"]["R1"])
        borderline_fq = set(sample_info.get("borderline", {}).get("R1", []))
        logger.info(
            f"Excel: {len(benign_fq)} benign, {len(malignant_fq)} malignant, "
            f"{len(borderline_fq)} borderline samples"
        )

        benign_files     = []
        malignant_files  = []
        borderline_files = []

        logger.info(f"Grouping {len(methylation_files)} methylation files by Excel")
        for meth_file in methylation_files:
            sample_id = os.path.basename(meth_file)
            for ext in [".txt", ".bedGraph", ".cov", ".gz"]:
                sample_id = sample_id.replace(ext, "")
            sample_id = sample_id.split(".")[0].split("_")[0]

            def _in_set(sid, fq_set):
                for r1 in fq_set:
                    r1_base   = os.path.splitext(os.path.splitext(r1)[0])[0]
                    r1_sample = r1_base.split(".")[0].split("_")[0]
                    if r1_sample == sid:
                        return True
                return False

            if _in_set(sample_id, benign_fq):
                benign_files.append(meth_file)
            elif _in_set(sample_id, malignant_fq):
                malignant_files.append(meth_file)
            elif _in_set(sample_id, borderline_fq):
                borderline_files.append(meth_file)

        logger.info(
            f"Grouping result from {excel_file}:\n"
            f"  benign={len(benign_files)}: {benign_files}\n"
            f"  malignant={len(malignant_files)}: {malignant_files}\n"
            f"  borderline={len(borderline_files)}: {borderline_files}"
        )
        return benign_files, malignant_files, borderline_files

    except Exception as e:
        logger.error(f"Error grouping by {excel_file}: {e}")
        benign_files, malignant_files = group_methylation_files_v1(methylation_files)
        return benign_files, malignant_files, []


def parse_sample_info(excel_file=None, category="plasma"):
    """
    Parse sample info Excel and return {benign: {R1, R2}, malignant: {R1, R2}}.

    Expects columns: fq_R1, fq_R2, category (benign/malignant/borderline), subtype.
    Borderline samples are excluded from both groups.
    """
    try:
        if excel_file is None:
            excel_file = _resolve_excel(config.get("data", {}).get("sample_info_excel", _DEFAULT_EXCEL))

        xl = pd.ExcelFile(excel_file)
        sheet_names = xl.sheet_names
        logger.info(f"Excel sheets: {sheet_names}")

        all_benign_R1     = []
        all_benign_R2     = []
        all_malignant_R1  = []
        all_malignant_R2  = []
        all_borderline_R1 = []
        all_borderline_R2 = []

        for sheet_name in sheet_names:
            if category and sheet_name != category:
                logger.info(f"Skipping sheet '{sheet_name}' (want '{category}')")
                continue

            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            logger.info(f"Sheet '{sheet_name}': {len(df)} rows")

            if "category" not in df.columns:
                logger.warning(f"Sheet '{sheet_name}' has no 'category' column, skipping")
                continue

            df["category"] = df["category"].astype(str).str.strip()
            if "subtype" in df.columns:
                df["subtype"] = df["subtype"].astype(str).str.strip()

            borderline_df = df[df["category"] == "borderline"]
            if not borderline_df.empty:
                logger.info(f"  borderline: {len(borderline_df)} samples")

            benign_df    = df[df["category"] == "benign"]
            malignant_df = df[df["category"] == "malignant"]

            logger.info(f"  benign={len(benign_df)}  malignant={len(malignant_df)}  "
                        f"borderline={len(borderline_df)}")

            if "fq_R1" in df.columns:
                all_benign_R1.extend(benign_df["fq_R1"].dropna().tolist())
                all_malignant_R1.extend(malignant_df["fq_R1"].dropna().tolist())
                all_borderline_R1.extend(borderline_df["fq_R1"].dropna().tolist())
            if "fq_R2" in df.columns:
                all_benign_R2.extend(benign_df["fq_R2"].dropna().tolist())
                all_malignant_R2.extend(malignant_df["fq_R2"].dropna().tolist())
                all_borderline_R2.extend(borderline_df["fq_R2"].dropna().tolist())

        all_benign_R1     = list(set(all_benign_R1))
        all_benign_R2     = list(set(all_benign_R2))
        all_malignant_R1  = list(set(all_malignant_R1))
        all_malignant_R2  = list(set(all_malignant_R2))
        all_borderline_R1 = list(set(all_borderline_R1))
        all_borderline_R2 = list(set(all_borderline_R2))

        logger.info(f"Parsed {excel_file}: {len(all_benign_R1)} benign, "
                    f"{len(all_malignant_R1)} malignant, {len(all_borderline_R1)} borderline")
        return {
            "benign":     {"R1": all_benign_R1,     "R2": all_benign_R2},
            "malignant":  {"R1": all_malignant_R1,  "R2": all_malignant_R2},
            "borderline": {"R1": all_borderline_R1, "R2": all_borderline_R2},
        }

    except Exception as e:
        logger.error(f"Error parsing Excel: {e}")
        return None


def resolve_patient_files(csv_path, excel_file=None, data_dir=None):
    """
    Resolve R1/R2 FastQ paths for a list of patient IDs from a CSV file.

    Args:
        csv_path:   One patient ID per line (no header).
        excel_file: Path to sample info Excel. Defaults to config value.
        data_dir:   Directory containing FastQ files. Defaults to config plasma_dir.

    Returns:
        list: FastQ paths (R1 and R2 interleaved).
    """
    with open(csv_path, encoding="utf-8") as f:
        patient_ids = [line.strip() for line in f if line.strip()]
    logger.info(f"Loaded {len(patient_ids)} patient IDs from {csv_path}")

    if excel_file is None:
        excel_file = config.get("data", {}).get("sample_info_excel", _DEFAULT_EXCEL)

    if data_dir is None:
        data_dir = config.get("data", {}).get("plasma_dir", "")
    data_dir = _to_posix_path(data_dir)

    xl = pd.ExcelFile(excel_file)
    pid_to_files = {}
    for sheet_name in xl.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        if "fq_R1" not in df.columns:
            continue
        if "category" in df.columns:
            df["category"] = df["category"].astype(str).str.strip()
        if "subtype" in df.columns:
            df["subtype"] = df["subtype"].astype(str).str.strip()
        for _, row in df.iterrows():
            r1 = str(row["fq_R1"]).strip() if pd.notna(row.get("fq_R1")) else ""
            r2 = str(row.get("fq_R2", "")).strip() if pd.notna(row.get("fq_R2")) else ""
            if not r1 or r1 == "nan":
                continue
            pid      = r1.split(".")[0]
            category = str(row.get("category", "")).strip() if "category" in df.columns else ""
            subtype  = str(row.get("subtype",  "")).strip() if "subtype"  in df.columns else ""
            if pid not in pid_to_files:
                pid_to_files[pid] = {"R1": r1, "R2": r2,
                                     "category": category, "subtype": subtype}

    logger.info(f"Excel: {len(pid_to_files)} sample mappings")

    fastq_files       = []
    missing           = []
    skipped_borderline = []
    for pid in patient_ids:
        if pid not in pid_to_files:
            logger.warning(f"  Patient ID {pid} not found in Excel")
            missing.append(pid)
            continue
        info = pid_to_files[pid]
        if info.get("category") == "borderline":
            skipped_borderline.append(pid)
            continue
        r1_path = os.path.join(data_dir, info["R1"])
        r2_name = info["R2"]
        r2_path = os.path.join(data_dir, r2_name) if r2_name else ""
        if not os.path.exists(r1_path):
            logger.warning(f"  R1 not found: {r1_path}")
            missing.append(pid)
            continue
        if r2_name and not os.path.exists(r2_path):
            logger.warning(f"  R2 not found: {r2_path}")
            missing.append(pid)
            continue
        fastq_files.append(r1_path)
        if r2_name:
            fastq_files.append(r2_path)

    found = len(patient_ids) - len(missing) - len(skipped_borderline)
    logger.info(f"Resolved {found}/{len(patient_ids)} patient IDs to FastQ files")
    if skipped_borderline:
        logger.info(f"  Skipped borderline ({len(skipped_borderline)}): {skipped_borderline}")
    if missing:
        logger.warning(f"  Missing/not-found IDs: {missing}")

    return fastq_files


def load_ca125_map(excel_file=None, category="plasma"):
    """
    Return {normalized_case_id: float} for samples with a valid CA-125 value.

    Keys are normalised the same way as classifier_utils._norm_id:
    '640849.0' → '640849', plain strings kept as-is.

    Accepts column names: ca125, CA125, CA-125, ca_125, ca 125 (case-insensitive).
    Skips rows where the value is NaN, empty, or not numeric.
    Returns an empty dict if the column is absent — callers handle this gracefully.
    """
    _ALIASES = {"ca125", "ca-125", "ca_125", "ca 125", "ca125(u/ml)", "ca-125(u/ml)"}

    def _norm(raw):
        s = str(raw).strip()
        try:
            return str(int(float(s)))
        except (ValueError, OverflowError):
            return s

    if excel_file is None:
        excel_file = _resolve_excel(config.get("data", {}).get("sample_info_excel", _DEFAULT_EXCEL))
    if not os.path.exists(excel_file):
        logger.warning(f"[ca125] Excel not found: {excel_file}")
        return {}

    result = {}
    try:
        xl = pd.ExcelFile(excel_file)
        for sheet_name in xl.sheet_names:
            if category and sheet_name != category:
                continue
            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            col_map = {c.strip().lower(): c for c in df.columns}
            ca125_col = next((col_map[a] for a in _ALIASES if a in col_map), None)
            if ca125_col is None:
                logger.debug(f"[ca125] Sheet '{sheet_name}': no CA-125 column — skipping")
                continue
            if "case_id" not in df.columns:
                logger.warning(f"[ca125] Sheet '{sheet_name}': no 'case_id' column — cannot match")
                continue
            for _, row in df.iterrows():
                case_id = _norm(row["case_id"]) if pd.notna(row.get("case_id")) else ""
                val     = row[ca125_col]
                if case_id and case_id != "nan" and pd.notna(val):
                    try:
                        result[case_id] = float(val)
                    except (ValueError, TypeError):
                        pass
        logger.info(f"[ca125] Loaded {len(result)} CA-125 values from {excel_file}")
    except Exception as exc:
        logger.warning(f"[ca125] Failed to load CA-125 values: {exc}")
    return result


if __name__ == "__main__":
    excel_file = _resolve_excel(config.get("data", {}).get("sample_info_excel", _DEFAULT_EXCEL))
    print(parse_sample_info(excel_file))
