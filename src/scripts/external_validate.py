#!/usr/bin/env python3
"""
external_validate.py  —  P3 Public-dataset cross-validation (Step 19, optional)

Validates panel CpG sites against two public 450K reference datasets:

  1. GSE40279  (whole blood, 656 samples)
     → blood_mean_beta, blood_pct_high: confirm panel sites are unmethylated in blood
       (cfDNA markers should have blood_mean_beta < 0.1 to avoid blood-cell confounding)

  2. TCGA-OV HumanMethylation450 (ovarian tumour tissue, ~590 samples)
     → tcga_mean_beta, presence of signal in OC tissue

Coordinate mapping pipeline
----------------------------
  panel sites  (hg38, chr:pos)
        │
        ▼  liftOver hg38→hg19 (data/external/hg38ToHg19.over.chain.gz)
  panel sites  (hg19)
        │
        ▼  match within ±2 bp
  Illumina 450K cg probes  (GPL13534, from data/external/; see README prerequisites)
        │
        ▼  extract rows from GSE40279 / TCGA-OV
  β values  (0–1)

Usage
-----
  python src/scripts/external_validate.py \\
      --result-dir result/20260406 \\
      --data-dir   data

Outputs
-------
  result/<ts>/external_validation.csv
"""

import os
import sys
import gzip
import bisect
import logging
import argparse
import urllib.request
import urllib.error
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────
#  Constants
# ─────────────────────────────────────────────────────────────

_GPL13534_URL   = ("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/"
                   "GPL13534/soft/GPL13534_family.soft.gz")
# Probe manifest candidates (checked in order; first found wins)
_PROBE_MANIFEST_CANDIDATES = [
    "GPL13534_platform.txt",          # GEO → Data table → Download full table (~100 MB)
    "GPL13534_platform.txt.gz",
    "HM450.hg19.manifest.tsv.gz",    # Zhou lab (if accessible)
    "HM450.hg19.manifest.tsv",
    "GPL13534_family.soft.gz",        # NCBI family SOFT — do NOT auto-download (214 GB)
]
_ZHOU_HM450_URL  = ("https://zhouserver.research.chop.edu/InfiniumAnnotation/"
                    "20180909/HM450/HM450.hg19.manifest.tsv.gz")
_ZHOU_HM450_FNAME = "HM450.hg19.manifest.tsv.gz"

# Blood baseline threshold: sites with blood mean β ≥ this are flagged
BLOOD_THRESHOLD = 0.10

# Coordinate match window (bp): how close hg19 liftover position must be to a 450K probe
MATCH_WINDOW_BP = 2


# ─────────────────────────────────────────────────────────────
#  GPL13534 download & parsing
# ─────────────────────────────────────────────────────────────

def _download_gpl13534(data_dir: str) -> str:
    """
    Locate or download the Illumina 450K probe coordinate file.

    Checks _PROBE_MANIFEST_CANDIDATES in order; for the Zhou lab file only,
    attempts auto-download if not present.  Returns local path, or "" if none found.
    """
    os.makedirs(data_dir, exist_ok=True)

    for fname in _PROBE_MANIFEST_CANDIDATES:
        path = os.path.join(data_dir, fname)
        if os.path.exists(path):
            logger.info(f"  Probe manifest found: {path}")
            return path

        # Auto-download only the small Zhou lab file
        if fname == _ZHOU_HM450_FNAME:
            logger.info(f"  Trying to download Zhou lab HM450 hg19 manifest (~10 MB) ...")
            logger.info(f"  URL: {_ZHOU_HM450_URL}")
            try:
                req = urllib.request.Request(_ZHOU_HM450_URL,
                                             headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(req, timeout=120) as resp:
                    raw = resp.read()
                with open(path, "wb") as f:
                    f.write(raw)
                logger.info(f"  Downloaded: {len(raw) // 1024 // 1024} MB → {path}")
                return path
            except Exception as exc:
                logger.warning(f"  Auto-download failed: {exc}")

    logger.error("  No probe manifest found. Probe-to-cg mapping will be skipped.")
    logger.error("  To enable it, place one of the following in  %s :", data_dir)
    logger.error("    GPL13534_platform.txt  — NCBI GEO platform table (recommended)")
    logger.error("      1. Open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534")
    logger.error("         Click 'Data table' → 'Download full table...' → save as GPL13534_platform.txt")
    logger.error("    HM450.hg19.manifest.tsv.gz  — Zhou lab hg19 manifest (~10 MB)")
    return ""


def _build_probe_map(gpl_path: str) -> dict:
    """
    Parse a probe coordinate file → {(chrom_hg19, pos_hg19_0based): cg_id}.

    Auto-detects format by peeking at the first non-empty line:
    • Zhou lab TSV  (HM450.hg19.manifest.tsv.gz)
        Header contains 'CpG_chrm' / 'CpG_beg'; coordinates are 0-based.
    • NCBI GEO SOFT (GPL13534_family.soft.gz)
        Data table between '!platform_table_begin' and '!platform_table_end';
        MAPINFO column is 1-based (converted to 0-based here).
    """
    if not gpl_path or not os.path.exists(gpl_path):
        return {}

    logger.info(f"  Parsing probe coordinates: {os.path.basename(gpl_path)}")
    probe_map: dict = {}
    opener = gzip.open if gpl_path.endswith(".gz") else open

    # ── peek to detect format ─────────────────────────────────────────────
    fmt = "unknown"
    try:
        with opener(gpl_path, "rt", encoding="latin-1", errors="replace") as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped:
                    continue
                if "CpG_chrm" in stripped or "CpG_beg" in stripped:
                    fmt = "zhou"
                elif stripped.startswith("!platform_table_begin") or stripped.startswith("^PLATFORM"):
                    fmt = "soft"
                elif stripped.startswith("#") or stripped.startswith("!") or stripped.startswith("^"):
                    continue   # metadata / comment line — keep peeking
                elif "\t" in stripped and stripped.upper().startswith("ID"):
                    fmt = "geo_tsv"
                break
    except Exception as exc:
        logger.error(f"  Cannot read {gpl_path}: {exc}")
        return {}

    if fmt == "unknown":
        fmt = "soft"   # best-effort fallback

    logger.info(f"  Detected format: {fmt}")

    # ── Zhou lab TSV ──────────────────────────────────────────────────────
    if fmt == "zhou":
        try:
            with opener(gpl_path, "rt", encoding="utf-8", errors="replace") as fh:
                header = None
                id_idx = chr_idx = beg_idx = None
                for line in fh:
                    line = line.rstrip("\n\r")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if header is None:
                        header = parts
                        hup = [h.upper() for h in header]
                        id_idx  = next((i for i, h in enumerate(hup)
                                        if h in ("PROBE_ID", "ILMNID", "NAME")), 0)
                        chr_idx = next((i for i, h in enumerate(hup)
                                        if "CPG_CHRM" in h or h == "CHR"), None)
                        beg_idx = next((i for i, h in enumerate(hup)
                                        if "CPG_BEG" in h or "MAPINFO" in h), None)
                        if chr_idx is None or beg_idx is None:
                            logger.error(f"  Cannot find chr/pos columns. Header: {header[:10]}")
                            return {}
                        is_1based = "MAPINFO" in hup[beg_idx]
                        continue
                    if len(parts) <= max(id_idx, chr_idx, beg_idx):
                        continue
                    cg_id = parts[id_idx].strip()
                    if not cg_id.startswith("cg"):
                        continue
                    chrom_raw = parts[chr_idx].strip()
                    beg_raw   = parts[beg_idx].strip()
                    if not chrom_raw or not beg_raw or beg_raw in ("NA", ""):
                        continue
                    try:
                        pos = int(float(beg_raw))
                        if is_1based:
                            pos -= 1          # 1-based → 0-based
                        chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
                        probe_map[(chrom, pos)] = cg_id
                    except (ValueError, TypeError):
                        continue
        except Exception as exc:
            logger.error(f"  Zhou TSV parse error: {exc}")
            return {}

    # ── NCBI GEO SOFT ─────────────────────────────────────────────────────
    elif fmt == "soft":
        try:
            with opener(gpl_path, "rt", encoding="latin-1", errors="replace") as fh:
                in_table = False
                header   = None
                id_idx = chr_idx = mapinfo_idx = None
                for line in fh:
                    line = line.rstrip("\n\r")
                    if line == "!platform_table_begin":
                        in_table = True
                        continue
                    if line == "!platform_table_end":
                        break
                    # Some GEO downloads omit !platform_table_begin; detect header row directly
                    if not in_table and line.startswith("ID\t") and "CHR" in line and "MAPINFO" in line:
                        in_table = True
                        # fall through — parse this line as the header
                    if not in_table or not line:
                        continue
                    parts = line.split("\t")
                    if header is None:
                        header = [h.strip().upper() for h in parts]
                        id_idx      = 0
                        chr_idx     = next((i for i, h in enumerate(header) if h == "CHR"), None)
                        mapinfo_idx = next((i for i, h in enumerate(header) if h == "MAPINFO"), None)
                        if chr_idx is None or mapinfo_idx is None:
                            logger.error(f"  CHR/MAPINFO not found. Header: {header[:15]}")
                            return {}
                        continue
                    if len(parts) <= max(id_idx, chr_idx, mapinfo_idx):
                        continue
                    cg_id     = parts[id_idx].strip()
                    chrom_raw = parts[chr_idx].strip()
                    minfo_raw = parts[mapinfo_idx].strip()
                    if not cg_id.startswith("cg") or not chrom_raw or not minfo_raw:
                        continue
                    try:
                        pos   = int(float(minfo_raw)) - 1   # 1-based → 0-based
                        chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
                        probe_map[(chrom, pos)] = cg_id
                    except (ValueError, TypeError):
                        continue
        except Exception as exc:
            logger.error(f"  SOFT parse error: {exc}")
            return {}

    # ── GEO plain-text Data-table (GPL13534_platform.txt) ─────────────────
    elif fmt == "geo_tsv":
        try:
            with opener(gpl_path, "rt", encoding="latin-1", errors="replace") as fh:
                header = None
                id_idx = chr_idx = mapinfo_idx = None
                for line in fh:
                    line = line.rstrip("\n\r")
                    if not line or line.startswith("#") or line.startswith("!") or line.startswith("^"):
                        continue
                    parts = line.split("\t")
                    if header is None:
                        header = [h.strip().upper() for h in parts]
                        id_idx      = 0
                        chr_idx     = next((i for i, h in enumerate(header) if h == "CHR"), None)
                        mapinfo_idx = next((i for i, h in enumerate(header) if h == "MAPINFO"), None)
                        if chr_idx is None or mapinfo_idx is None:
                            logger.error(f"  CHR/MAPINFO not found. Header: {header[:15]}")
                            return {}
                        continue
                    if len(parts) <= max(id_idx, chr_idx, mapinfo_idx):
                        continue
                    cg_id     = parts[id_idx].strip()
                    chrom_raw = parts[chr_idx].strip()
                    minfo_raw = parts[mapinfo_idx].strip()
                    if not cg_id.startswith("cg") or not chrom_raw or not minfo_raw or minfo_raw == "NA":
                        continue
                    try:
                        pos   = int(float(minfo_raw)) - 1   # MAPINFO is 1-based
                        chrom = chrom_raw if chrom_raw.startswith("chr") else f"chr{chrom_raw}"
                        probe_map[(chrom, pos)] = cg_id
                    except (ValueError, TypeError):
                        continue
        except Exception as exc:
            logger.error(f"  GEO TSV parse error: {exc}")
            return {}

    logger.info(f"  Loaded {len(probe_map):,} probe positions")
    return probe_map


# ─────────────────────────────────────────────────────────────
#  Chain file liftover (hg38 → hg19)
# ─────────────────────────────────────────────────────────────

def _parse_chain(chain_gz: str) -> dict:
    """
    Parse UCSC liftOver chain file (hg38ToHg19.over.chain.gz).
    Returns {src_chrom: sorted list of (src_start, src_end, tgt_chrom, offset)}
    where  hg19_pos = hg38_pos + offset  (valid only within [src_start, src_end)).
    Only '+' → '+' strand blocks are included (covers >99% of autosomal CpGs).
    """
    logger.info(f"  Parsing chain file: {chain_gz}")
    blocks = defaultdict(list)

    try:
        with gzip.open(chain_gz, "rt") as f:
            src_chrom = tgt_chrom = None
            cur_src = cur_tgt = 0
            tgt_plus = True

            for line in f:
                line = line.rstrip("\n")
                if not line:
                    src_chrom = None
                    continue

                if line.startswith("chain"):
                    p = line.split()
                    src_chrom  = p[2]
                    tgt_chrom  = p[7]
                    tgt_plus   = (p[9] == "+")
                    cur_src    = int(p[5])
                    tgt_size   = int(p[8])
                    tgt_start  = int(p[10])
                    tgt_end    = int(p[11])
                    if tgt_plus:
                        cur_tgt = tgt_start
                    else:
                        # minus strand: positions decrease
                        cur_tgt = tgt_size - tgt_end
                    continue

                if src_chrom is None:
                    continue

                p = line.split()
                if not p:
                    continue
                size = int(p[0])

                if tgt_plus:
                    offset = cur_tgt - cur_src
                    blocks[src_chrom].append((cur_src, cur_src + size, tgt_chrom, offset))

                if len(p) == 3:
                    cur_src += size + int(p[1])
                    if tgt_plus:
                        cur_tgt += size + int(p[2])
                    else:
                        cur_tgt -= size + int(p[2])

    except Exception as exc:
        logger.error(f"  Chain parse error: {exc}")
        return {}

    # Sort each chromosome's blocks by src_start for binary search
    for chrom in blocks:
        blocks[chrom].sort(key=lambda x: x[0])

    total = sum(len(v) for v in blocks.values())
    logger.info(f"  Loaded {total:,} chain blocks")
    return dict(blocks)


def _liftover(chrom: str, pos_hg38: int, blocks: dict):
    """Convert a single hg38 0-based position to hg19. Returns (chrom, pos) or None."""
    chrom_blocks = blocks.get(chrom)
    if not chrom_blocks:
        return None

    starts = [b[0] for b in chrom_blocks]
    idx = bisect.bisect_right(starts, pos_hg38) - 1
    if idx < 0 or idx >= len(chrom_blocks):
        return None

    src_start, src_end, tgt_chrom, offset = chrom_blocks[idx]
    if src_start <= pos_hg38 < src_end:
        return (tgt_chrom, pos_hg38 + offset)
    return None


# ─────────────────────────────────────────────────────────────
#  Panel → 450K probe matching
# ─────────────────────────────────────────────────────────────

def _match_probes(panel_features: list, chain_blocks: dict, probe_map: dict,
                  window: int = MATCH_WINDOW_BP) -> dict:
    """
    For each panel feature (hg38 chr:pos), liftover to hg19 and find
    the nearest 450K probe within ±window bp.

    Returns {feature: cg_id or None}
    """
    logger.info(f"  Matching {len(panel_features)} panel sites to 450K probes (window ±{window} bp) ...")

    # Build sorted probe list per chrom for range lookup
    probes_by_chrom = defaultdict(list)
    for (chrom, pos), cg_id in probe_map.items():
        probes_by_chrom[chrom].append((pos, cg_id))
    for chrom in probes_by_chrom:
        probes_by_chrom[chrom].sort()

    matches = {}
    n_lifted = n_matched = 0

    for feat in panel_features:
        chrom38, pos_str = feat.split(":", 1)
        pos38 = int(pos_str)

        hg19 = _liftover(chrom38, pos38, chain_blocks)
        if hg19 is None:
            matches[feat] = None
            continue
        n_lifted += 1

        chrom19, pos19 = hg19
        probe_list = probes_by_chrom.get(chrom19, [])
        if not probe_list:
            matches[feat] = None
            continue

        probe_pos_arr = [p[0] for p in probe_list]
        idx = bisect.bisect_left(probe_pos_arr, pos19 - window)

        best_cg, best_dist = None, window + 1
        while idx < len(probe_list) and probe_list[idx][0] <= pos19 + window:
            dist = abs(probe_list[idx][0] - pos19)
            if dist < best_dist:
                best_dist = dist
                best_cg   = probe_list[idx][1]
            idx += 1

        matches[feat] = best_cg
        if best_cg:
            n_matched += 1

    logger.info(f"  Liftover success: {n_lifted}/{len(panel_features)}  "
                f"Probe match (≤{window} bp): {n_matched}/{len(panel_features)}")
    return matches


# ─────────────────────────────────────────────────────────────
#  GEO series matrix parser
# ─────────────────────────────────────────────────────────────

def _parse_series_matrix(matrix_gz: str, probes_needed: set) -> pd.DataFrame:
    """
    Stream-parse a GEO series_matrix.txt.gz and return a DataFrame
    (probe × sample) of β values for only the probes in `probes_needed`.
    Memory-efficient: skips rows not in probes_needed.
    """
    logger.info(f"  Parsing {os.path.basename(matrix_gz)} for {len(probes_needed)} probes ...")

    rows = {}
    sample_ids = None
    in_table = False

    opener = gzip.open if matrix_gz.endswith(".gz") else open
    with opener(matrix_gz, "rt", encoding="latin-1", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin") or \
               line.startswith("!Series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end") or \
               line.startswith("!Series_matrix_table_end"):
                break
            if not in_table:
                continue

            parts = line.split("\t")
            if sample_ids is None:
                # Header row: "ID_REF" then sample IDs
                sample_ids = [p.strip('"') for p in parts[1:]]
                continue

            probe = parts[0].strip('"')
            if probe not in probes_needed:
                continue

            try:
                betas = [float(p) if p.strip() not in ("", "NA", "null", "NULL", "NaN")
                         else np.nan
                         for p in parts[1:]]
                rows[probe] = betas
            except ValueError:
                continue

    if not rows or sample_ids is None:
        return pd.DataFrame()

    df = pd.DataFrame(rows, index=sample_ids).T   # index = probe, columns = sample
    logger.info(f"  Loaded {len(df)} probes × {len(df.columns)} samples")
    return df


# ─────────────────────────────────────────────────────────────
#  TCGA Xena file parser
# ─────────────────────────────────────────────────────────────

def _parse_xena_matrix(xena_gz: str, probes_needed: set) -> pd.DataFrame:
    """
    Parse UCSC Xena methylation matrix (probe × sample, tab-delimited, gzipped).
    First column is probe ID (cg*), first row is sample IDs.
    """
    logger.info(f"  Parsing {os.path.basename(xena_gz)} for {len(probes_needed)} probes ...")

    rows = {}
    sample_ids = None

    opener = gzip.open if xena_gz.endswith(".gz") else open
    with opener(xena_gz, "rt", encoding="latin-1", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            parts = line.split("\t")
            if sample_ids is None:
                # Header: 'sample' then TCGA-XX-... IDs
                sample_ids = parts[1:]
                continue
            probe = parts[0]
            if probe not in probes_needed:
                continue
            try:
                betas = [float(p) if p.strip() not in ("", "NA", "NaN", "null")
                         else np.nan
                         for p in parts[1:]]
                rows[probe] = betas
            except ValueError:
                continue

    if not rows or sample_ids is None:
        return pd.DataFrame()

    df = pd.DataFrame(rows, index=sample_ids).T
    logger.info(f"  Loaded {len(df)} probes × {len(df.columns)} samples")
    return df


# ─────────────────────────────────────────────────────────────
#  Main validation function
# ─────────────────────────────────────────────────────────────

def compute_blood_filtered_panel(result_dir: str, external_dir: str,
                                 threshold: float = BLOOD_THRESHOLD) -> dict:
    """Step 14a — standalone blood-baseline filter.

    Reads ``panel_cpg.csv``, lifts each panel site hg38 → hg19, matches to a
    GPL13534 cg probe, looks up GSE40279 (n=656 whole-blood samples) mean β,
    and writes ``panel_cpg_blood_filtered.csv`` with sites whose
    ``blood_mean_beta ≥ threshold`` removed. Sites without a matched 450K
    probe are retained as-is (no blood-baseline evidence available).

    Uses the same parse helpers as Step 19 ``run_external_validation`` so the
    output is bit-identical regardless of which step writes it.

    Returns a summary dict (empty if a critical input is missing).
    """
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if not os.path.exists(panel_path):
        logger.error(f"[blood_filter] panel_cpg.csv not found: {panel_path} "
                     "— run --step feature_analysis first")
        return {}

    panel_df       = pd.read_csv(panel_path)
    panel_features = list(panel_df["feature"].astype(str))
    n_total        = len(panel_features)
    logger.info(f"[blood_filter] panel sites: {n_total}")

    # ── 1. GPL13534 probe map (hg19) ────────────────────────────────────
    gpl_cache = _download_gpl13534(external_dir)
    probe_map = _build_probe_map(gpl_cache) if gpl_cache else {}
    if not probe_map:
        logger.warning("[blood_filter] No GPL13534 probe map — cannot evaluate "
                       "blood baseline; panel_cpg_blood_filtered.csv NOT written")
        return {"n_total": n_total, "n_matched": 0, "n_removed": 0,
                "n_retained": n_total}

    # ── 2. Liftover hg38 → hg19, match to cg probes ─────────────────────
    chain_path = os.path.join(external_dir, "hg38ToHg19.over.chain.gz")
    if not os.path.exists(chain_path):
        logger.warning(f"[blood_filter] chain file not found: {chain_path} "
                       "— skipping (no filtered panel written)")
        return {"n_total": n_total, "n_matched": 0, "n_removed": 0,
                "n_retained": n_total}

    chain_blocks    = _parse_chain(chain_path)
    feature_to_probe = _match_probes(panel_features, chain_blocks, probe_map)
    matched_probes  = {cg for cg in feature_to_probe.values() if cg}
    n_matched       = len(matched_probes)
    logger.info(f"[blood_filter] panel sites matched to a 450K probe: {n_matched}/{n_total}")

    if not matched_probes:
        logger.warning("[blood_filter] no panel sites matched to a 450K probe "
                       "— filtered panel not written")
        return {"n_total": n_total, "n_matched": 0, "n_removed": 0,
                "n_retained": n_total}

    # ── 3. Blood baseline lookup (GSE40279) ─────────────────────────────
    gse_path = os.path.join(external_dir, "GSE40279_series_matrix.txt.gz")
    if not os.path.exists(gse_path):
        logger.warning(f"[blood_filter] GSE40279 not found: {gse_path} "
                       "— filtered panel not written")
        return {"n_total": n_total, "n_matched": n_matched,
                "n_removed": 0, "n_retained": n_total}

    blood_df = _parse_series_matrix(gse_path, matched_probes)
    if blood_df.empty:
        logger.warning("[blood_filter] GSE40279 returned no matching probes "
                       "— filtered panel not written")
        return {"n_total": n_total, "n_matched": n_matched,
                "n_removed": 0, "n_retained": n_total}

    # ── 4. Compute blood_mean_beta per site, identify fails ─────────────
    blood_fail_sites = set()
    blood_rows       = []
    for feat in panel_features:
        cg = feature_to_probe.get(feat)
        if cg and cg in blood_df.index:
            betas = blood_df.loc[cg].dropna().values
            if len(betas) > 0:
                mean_beta = round(float(np.mean(betas)), 4)
                passed    = mean_beta < threshold
                if not passed:
                    blood_fail_sites.add(feat)
                blood_rows.append({"feature": feat, "cg_probe": cg,
                                   "blood_mean_beta": mean_beta,
                                   "blood_pass": passed})
                continue
        blood_rows.append({"feature": feat, "cg_probe": cg or "",
                           "blood_mean_beta": None, "blood_pass": None})

    # ── 5. Write filtered panel + summary ───────────────────────────────
    filtered_df    = panel_df[~panel_df["feature"].astype(str).isin(blood_fail_sites)].copy()
    filtered_path  = os.path.join(result_dir, "panel_cpg_blood_filtered.csv")
    filtered_df.to_csv(filtered_path, index=False)

    summary_df   = pd.DataFrame(blood_rows)
    summary_path = os.path.join(result_dir, "panel_cpg_blood_summary.csv")
    summary_df.to_csv(summary_path, index=False)

    n_removed  = len(blood_fail_sites)
    n_retained = len(filtered_df)
    logger.info(
        f"[blood_filter] threshold={threshold}: {n_removed} site(s) removed, "
        f"{n_retained}/{n_total} retained → {os.path.basename(filtered_path)}"
    )
    if blood_fail_sites:
        logger.info(f"[blood_filter] removed sites: {sorted(blood_fail_sites)}")

    return {
        "n_total":    n_total,
        "n_matched":  n_matched,
        "n_removed":  n_removed,
        "n_retained": n_retained,
        "threshold":  threshold,
    }


def run_external_validation(result_dir: str, data_dir: str, external_dir: str):
    """
    Full P3 cross-validation pipeline.
    result_dir:   pipeline result directory (contains panel_cpg.csv)
    data_dir:     project data/ dir (for GPL13534 cache)
    external_dir: directory with downloaded public datasets
    """
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if not os.path.exists(panel_path):
        logger.error(f"panel_cpg.csv not found: {panel_path}  (run feature_analysis first)")
        return

    panel_df       = pd.read_csv(panel_path)
    panel_features = list(panel_df["feature"])
    logger.info(f"Panel sites: {len(panel_features)}")

    # ── 1. GPL13534 probe map (hg19) ────────────────────────────────────────
    gpl_cache = _download_gpl13534(external_dir)
    probe_map = _build_probe_map(gpl_cache)
    if not probe_map:
        logger.warning("No probe map built — probe_matched column will be all False")

    # ── 2. Liftover panel hg38 → hg19, match to cg probes ──────────────────
    chain_path = os.path.join(external_dir, "hg38ToHg19.over.chain.gz")
    if not os.path.exists(chain_path):
        logger.error(f"Chain file not found: {chain_path}")
        return

    chain_blocks = _parse_chain(chain_path)
    feature_to_probe = _match_probes(panel_features, chain_blocks, probe_map)

    matched_probes = {cg for cg in feature_to_probe.values() if cg}
    n_matched = len(matched_probes)
    logger.info(f"Unique 450K probes matched: {n_matched}")
    if not matched_probes:
        logger.warning("No panel sites matched to 450K probes — blood/TCGA lookup will be skipped")

    # ── 3. Blood baseline (GSE40279) ────────────────────────────────────────
    gse_path = os.path.join(external_dir, "GSE40279_series_matrix.txt.gz")
    blood_df = pd.DataFrame()
    if os.path.exists(gse_path):
        blood_df = _parse_series_matrix(gse_path, matched_probes)
    else:
        logger.warning(f"GSE40279 not found: {gse_path} — blood baseline skipped")

    # ── 4. TCGA-OV tissue ───────────────────────────────────────────────────
    tcga_path = os.path.join(external_dir, "HumanMethylation450.gz")
    tcga_df = pd.DataFrame()
    if os.path.exists(tcga_path):
        tcga_df = _parse_xena_matrix(tcga_path, matched_probes)
        n_tcga = len(tcga_df.columns) if not tcga_df.empty else 0
        if n_tcga < 50:
            logger.warning(
                f"TCGA-OV file has only {n_tcga} samples — likely a preview. "
                "Re-download the complete file from Xena browser."
            )
            tcga_df = pd.DataFrame()
    else:
        logger.warning(f"TCGA-OV file not found: {tcga_path} — tissue cross-check skipped")

    # ── 5. Build per-site result table ──────────────────────────────────────
    rows = []
    for feat in panel_features:
        cg = feature_to_probe.get(feat)
        row = {
            "feature":         feat,
            "cg_probe":        cg or "",
            "probe_matched":   cg is not None,
        }

        # Blood baseline
        if cg and not blood_df.empty and cg in blood_df.index:
            betas_blood = blood_df.loc[cg].dropna().values
            if len(betas_blood) > 0:
                row["blood_mean_beta"]   = round(float(np.mean(betas_blood)), 4)
                row["blood_median_beta"] = round(float(np.median(betas_blood)), 4)
                row["blood_pct_high"]    = round(float((betas_blood >= BLOOD_THRESHOLD).mean()), 4)
                row["blood_pass"]        = bool(row["blood_mean_beta"] < BLOOD_THRESHOLD)
            else:
                row.update({"blood_mean_beta": None, "blood_median_beta": None,
                            "blood_pct_high": None, "blood_pass": None})
        else:
            row.update({"blood_mean_beta": None, "blood_median_beta": None,
                        "blood_pct_high": None, "blood_pass": None})

        # TCGA tissue
        if cg and not tcga_df.empty and cg in tcga_df.index:
            betas_tcga = tcga_df.loc[cg].dropna().values
            if len(betas_tcga) > 0:
                row["tcga_mean_beta"]   = round(float(np.mean(betas_tcga)), 4)
                row["tcga_median_beta"] = round(float(np.median(betas_tcga)), 4)
                row["tcga_n_samples"]   = int(len(betas_tcga))
            else:
                row.update({"tcga_mean_beta": None, "tcga_median_beta": None, "tcga_n_samples": 0})
        else:
            row.update({"tcga_mean_beta": None, "tcga_median_beta": None, "tcga_n_samples": 0})

        rows.append(row)

    result = pd.DataFrame(rows)

    # ── 6. Summary ──────────────────────────────────────────────────────────
    n_total      = len(result)
    n_probe_hit  = int(result["probe_matched"].sum())
    n_blood_pass = int(result["blood_pass"].sum()) if "blood_pass" in result else 0
    n_blood_na   = int(result["blood_pass"].isna().sum()) if "blood_pass" in result else n_total

    logger.info("=" * 60)
    logger.info(f"External Validation Summary")
    logger.info(f"  Panel sites total:          {n_total}")
    logger.info(f"  Matched to 450K probe:      {n_probe_hit} ({n_probe_hit/n_total:.1%})")
    if n_blood_na < n_total:
        logger.info(f"  Blood baseline pass (<0.10): {n_blood_pass} "
                    f"({n_blood_pass/(n_total - n_blood_na):.1%} of matched)")
        n_blood_fail = n_probe_hit - n_blood_pass - n_blood_na
        if n_blood_fail > 0:
            fail_sites = result[result["blood_pass"] == False]["feature"].tolist()
            logger.warning(f"  Blood baseline FAIL sites ({n_blood_fail}):")
            for s in fail_sites[:10]:
                bmean = result.loc[result["feature"] == s, "blood_mean_beta"].values[0]
                logger.warning(f"    {s}  blood_mean={bmean:.3f}")
            if len(fail_sites) > 10:
                logger.warning(f"    ... and {len(fail_sites)-10} more (see CSV)")
    logger.info("=" * 60)

    # ── 7. Save external_validation.csv ─────────────────────────────────────
    out_path = os.path.join(result_dir, "external_validation.csv")
    result.to_csv(out_path, index=False)
    logger.info(f"Saved: {out_path}")

    # ── 8. Write blood-filtered panel_cpg_blood_filtered.csv ─────────────────
    #  Probe-matched sites with blood_mean_beta >= BLOOD_THRESHOLD are removed.
    #  Un-matched sites (no 450K probe) are retained as-is (no blood evidence).
    #  Classifiers (Steps 15 / 17) will auto-detect and prefer the filtered file.
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if os.path.exists(panel_path):
        blood_fail_sites = set(
            result.loc[
                (result["probe_matched"] == True) & (result["blood_pass"] == False),
                "feature",
            ]
        )
        panel_df = pd.read_csv(panel_path)
        filtered_df = panel_df[~panel_df["feature"].isin(blood_fail_sites)].copy()
        filtered_path = os.path.join(result_dir, "panel_cpg_blood_filtered.csv")
        filtered_df.to_csv(filtered_path, index=False)

        n_removed  = len(blood_fail_sites)
        n_retained = len(filtered_df)
        n_total_panel = len(panel_df)
        logger.info(
            f"Blood-filter: {n_removed} site(s) removed (blood_mean_beta >= {BLOOD_THRESHOLD}) "
            f"→ {n_retained}/{n_total_panel} sites retained in panel_cpg_blood_filtered.csv"
        )
        if n_removed:
            logger.warning(
                "Re-run --step classify (and --step temporal_validate) to apply the "
                "blood-filtered panel to classifier training."
            )

        # Summary CSV for reporting
        summary_path = os.path.join(result_dir, "blood_filter_summary.csv")
        pd.DataFrame([{
            "panel_total":         n_total_panel,
            "probe_matched":       n_probe_hit,
            "blood_fail_removed":  n_removed,
            "panel_retained":      n_retained,
            "blood_threshold":     BLOOD_THRESHOLD,
        }]).to_csv(summary_path, index=False)
        logger.info(f"Saved: {summary_path}")
    else:
        logger.warning(
            f"panel_cpg.csv not found at {panel_path} — blood-filtered panel not written. "
            "Run --step feature_analysis first."
        )

    return result


# ─────────────────────────────────────────────────────────────
#  CLI entry point
# ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="P3 external validation (blood + TCGA-OV)")
    ap.add_argument("--result-dir", required=True,
                    help="Pipeline result dir containing panel_cpg.csv")
    ap.add_argument("--data-dir",   default="data",
                    help="Project data/ directory (GPL13534 cached here)")
    ap.add_argument("--external-dir", default="data/external",
                    help="Directory with downloaded public datasets")
    args = ap.parse_args()

    run_external_validation(
        result_dir   = args.result_dir,
        data_dir     = args.data_dir,
        external_dir = args.external_dir,
    )


if __name__ == "__main__":
    main()
