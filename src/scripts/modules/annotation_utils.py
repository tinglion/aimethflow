#!/usr/bin/env python3
"""
annotation_utils.py - CpG site annotation module (Step 18)

Annotates each site in panel_cpg.csv with:
  1. Nearest gene + region type    (Ensembl REST API — no local files required)
  2. CpG island context            (UCSC BED auto-downloaded to data/ on first run)
  3. EWAS Catalog associations     (REST API, configurable window around each site)
  4. 88-candidate-region membership (config candidate_regions file)
  5. LOO-CV selection frequency    (from feature_selection_freq.csv if available)

Coordinate convention
---------------------
Features are stored as chr:start where start is the 0-based BED coordinate
produced by MethylDackel.  All external APIs use 1-based coordinates, so
every API call uses  pos_1b = pos_0b + 1.

Outputs (to result_dir/)
------------------------
  panel_annotation.csv     — full per-site annotation table
  annotation_cache.json    — cached API responses (re-runs skip network calls)
"""

import os
import json
import time
import gzip
import urllib.request
import urllib.error
import logging
import numpy as np
import pandas as pd
from collections import defaultdict

from utils.logger import logger

# ---------------------------------------------------------------------------
# External endpoints
# ---------------------------------------------------------------------------
_ENSEMBL_BASE  = "https://rest.ensembl.org"
_EWAS_API      = "https://www.ewascatalog.org/api/"
_EWAS_ATLAS    = "https://ngdc.cncb.ac.cn/ewas/api/study"
_PUBMED_SEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_CPGI_URL      = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/"
                  "cpgIslandExt.txt.gz")
_CPGI_FNAME    = "cpg_islands_hg38.bed"

# Keywords indicating plasma / liquid biopsy / early-detection literature
_PLASMA_KEYWORDS = frozenset([
    "plasma", "serum", "cfDNA", "cell-free", "liquid biopsy",
    "circulating", "ctDNA", "early detection", "early-stage",
    "screening", "noninvasive", "non-invasive",
])

# ---------------------------------------------------------------------------
# HTTP helper with retry + rate limiting
# ---------------------------------------------------------------------------

def _get_json(url, params=None, timeout=15, retries=3, delay=0.25):
    """GET JSON with retry and back-off.  Returns parsed object or None."""
    import urllib.parse
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            time.sleep(delay)
            req = urllib.request.Request(
                url,
                headers={"Accept": "application/json",
                         "Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 429:              # rate-limited
                time.sleep(2 ** (attempt + 1))
            elif e.code >= 500:
                time.sleep(1)
            else:
                return None               # 4xx (not rate limit) → give up
        except Exception:
            if attempt < retries - 1:
                time.sleep(1)
    return None


# ---------------------------------------------------------------------------
# CpG island BED
# ---------------------------------------------------------------------------

def _ensure_cpgi_bed(data_dir):
    """
    Return path to hg38 CpG island BED, auto-downloading from UCSC if absent.
    Returns None when download fails.
    """
    bed_path = os.path.join(data_dir, _CPGI_FNAME)
    if os.path.exists(bed_path):
        return bed_path

    logger.info(f"  Downloading CpG island BED from UCSC → {bed_path}")
    try:
        req = urllib.request.Request(
            _CPGI_URL,
            headers={"User-Agent": "Mozilla/5.0"},
        )
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.read()

        lines = []
        with gzip.open(__import__("io").BytesIO(raw), "rt") as gz:
            for line in gz:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    # UCSC cpgIslandExt: bin chrom chromStart chromEnd name ...
                    lines.append(f"{parts[1]}\t{parts[2]}\t{parts[3]}\t"
                                 f"{parts[4] if len(parts) > 4 else '.'}\n")

        os.makedirs(data_dir, exist_ok=True)
        with open(bed_path, "w") as fout:
            fout.writelines(lines)
        logger.info(f"  Downloaded {len(lines)} CpG islands → {bed_path}")
        return bed_path

    except Exception as exc:
        logger.warning(f"  CpG island download failed: {exc} — context will show 'N/A'")
        return None


def _load_cpgi_intervals(bed_path):
    """Load CpG island intervals into {chrom: sorted[(start, end)]}."""
    ivs = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return ivs
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.strip().split("\t")
            if len(p) >= 3:
                try:
                    ivs[p[0]].append((int(p[1]), int(p[2])))
                except ValueError:
                    pass
    for ch in ivs:
        ivs[ch].sort()
    return dict(ivs)


def _cpg_context(chrom, pos_0b, intervals):
    """
    Island / Shore (≤2 kb) / Shelf (2-4 kb) / Open sea.
    pos_0b is 0-based BED coordinate.
    """
    import bisect
    islands = intervals.get(chrom, [])
    if not islands:
        return "Open sea"

    starts = [s for s, _ in islands]
    idx = bisect.bisect_right(starts, pos_0b) - 1
    candidates = []
    if 0 <= idx < len(islands):
        candidates.append(islands[idx])
    if idx + 1 < len(islands):
        candidates.append(islands[idx + 1])

    min_dist = float("inf")
    for s, e in candidates:
        if s <= pos_0b <= e:
            return "Island"
        min_dist = min(min_dist, min(abs(pos_0b - s), abs(pos_0b - e)))

    if min_dist <= 2000:
        return "Shore"
    if min_dist <= 4000:
        return "Shelf"
    return "Open sea"


# ---------------------------------------------------------------------------
# Ensembl gene annotation
# ---------------------------------------------------------------------------

def _ensembl_nearest_gene(chrom, pos_1b, window=2000, delay=0.25):
    """
    Nearest gene via Ensembl REST overlap endpoint.
    pos_1b is 1-based (= bedGraph start + 1).
    Returns dict with gene_name, gene_id, biotype, distance_bp, region_type.
    """
    chrom_e = chrom.lstrip("chr")   # Ensembl uses "10" not "chr10"
    start   = max(1, pos_1b - window)
    end     = pos_1b + window
    url     = f"{_ENSEMBL_BASE}/overlap/region/human/{chrom_e}:{start}-{end}"

    data = _get_json(url, params={"feature": "gene"}, delay=delay)

    _empty = dict(gene_name="", gene_id="", biotype="",
                  distance_bp=None, region_type="unknown")
    if not data:
        return _empty

    best, best_dist = None, float("inf")
    for g in data:
        gs, ge = g.get("start", 0), g.get("end", 0)
        dist   = 0 if gs <= pos_1b <= ge else min(abs(pos_1b - gs),
                                                   abs(pos_1b - ge))
        if dist < best_dist:
            best_dist, best = dist, g

    if not best:
        return {**_empty, "region_type": "intergenic"}

    gs, ge, strand = (best.get("start", 0), best.get("end", 0),
                      best.get("strand", 1))
    prom_s = (max(1, gs - 2000) if strand == 1 else ge - 500)
    prom_e = (gs + 500           if strand == 1 else ge + 2000)

    if prom_s <= pos_1b <= prom_e:
        rtype = "promoter"
    elif gs <= pos_1b <= ge:
        rtype = "gene_body"
    else:
        rtype = "intergenic"

    return dict(
        gene_name   = best.get("external_name", best.get("id", "")),
        gene_id     = best.get("id", ""),
        biotype     = best.get("biotype", ""),
        distance_bp = 0 if gs <= pos_1b <= ge else best_dist,
        region_type = rtype,
    )


# ---------------------------------------------------------------------------
# EWAS Catalog
# ---------------------------------------------------------------------------

def _ewas_query(chrom, pos_1b, window=500, delay=0.25):
    """
    Query EWAS Catalog for CpG associations within ±window bp.
    Returns list of hit dicts.
    """
    chrom_e = chrom.lstrip("chr")
    data    = _get_json(_EWAS_API, params=dict(
        chromosome     = chrom_e,
        position_lower = pos_1b - window,
        position_upper = pos_1b + window,
        format         = "json",
    ), delay=delay)

    if data is None:
        return []
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get("results", data.get("data", []))
    return []


# ---------------------------------------------------------------------------
# EWAS Atlas (CNCB NGDC)
# ---------------------------------------------------------------------------

def _ewas_atlas_query(chrom, pos_1b, window=500, delay=0.25):
    """
    Query EWAS Atlas (CNCB NGDC) for associations within ±window bp.
    Returns list of hit dicts; empty list on failure.

    Endpoint: https://ngdc.cncb.ac.cn/ewas/api/study
    Positional params: chrom (no "chr"), bp_start, bp_end
    """
    chrom_e = chrom.lstrip("chr")
    data = _get_json(
        _EWAS_ATLAS,
        params=dict(
            chrom     = chrom_e,
            bp_start  = pos_1b - window,
            bp_end    = pos_1b + window,
            page_size = 100,
        ),
        delay=delay,
    )
    if data is None:
        return []
    # Response may be a list, or {"data": [...], "total": N}
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return data.get("data", data.get("results", []))
    return []


# ---------------------------------------------------------------------------
# PubMed E-utilities
# ---------------------------------------------------------------------------

def _pubmed_search(gene_name, ncbi_api_key="", max_hits=20, delay=0.25):
    """
    Search PubMed for gene_name + 'plasma methylation cancer'.
    Returns (count: int, pmids: list[str]).
    Requires gene_name to be non-empty; returns (0, []) otherwise.
    """
    if not gene_name or gene_name.strip() == "":
        return 0, []

    query = f'"{gene_name}" plasma methylation cancer'
    params = {
        "db":      "pubmed",
        "term":    query,
        "retmax":  max_hits,
        "retmode": "json",
    }
    if ncbi_api_key:
        params["api_key"] = ncbi_api_key

    data = _get_json(_PUBMED_SEARCH, params=params, delay=delay)
    if not data:
        return 0, []

    res   = data.get("esearchresult", {})
    count = int(res.get("count", 0))
    pmids = res.get("idlist", [])
    return count, pmids


# ---------------------------------------------------------------------------
# Plasma-evidence flag
# ---------------------------------------------------------------------------

def _has_plasma_evidence(traits_combined: str) -> bool:
    """
    Return True if any plasma/liquid-biopsy keyword appears in traits text
    (case-insensitive).  traits_combined is the concatenated trait strings
    from EWAS Catalog + EWAS Atlas.
    """
    text = traits_combined.lower()
    return any(kw.lower() in text for kw in _PLASMA_KEYWORDS)


# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------

def _load_cache(path):
    if os.path.exists(path):
        try:
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass
    return {}


def _save_cache(cache, path):
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(cache, f, ensure_ascii=False, indent=2)
    except Exception as exc:
        logger.warning(f"  Failed to write annotation cache: {exc}")


# ---------------------------------------------------------------------------
# Site loader
# ---------------------------------------------------------------------------

def _load_sites(result_dir):
    """
    Load panel_cpg.csv; optionally merge LOO-CV selection frequency.
    Returns DataFrame with at least a 'feature' column.
    """
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if not os.path.exists(panel_path):
        raise FileNotFoundError(
            f"panel_cpg.csv not found — run feature_analysis first: {panel_path}"
        )
    df = pd.read_csv(panel_path)
    if "feature" not in df.columns:
        raise ValueError("panel_cpg.csv must have a 'feature' column (chr:pos format)")

    freq_path = os.path.join(result_dir, "feature_selection_freq.csv")
    if os.path.exists(freq_path):
        freq = pd.read_csv(freq_path)
        keep = [c for c in ("feature", "selection_rate", "selection_count")
                if c in freq.columns]
        if "feature" in keep:
            df = df.merge(freq[keep], on="feature", how="left")

    return df


# ---------------------------------------------------------------------------
# Supplementary Table S2 export — mirror the annotation table into the
# manuscript docs/ tree as a formatted xlsx
# ---------------------------------------------------------------------------

def _resolve_project_path(rel_path: str) -> str:
    """Resolve ``rel_path`` against the project root.

    annotation_utils.py lives at src/scripts/modules/annotation_utils.py, so
    the project root is four ``dirname`` calls up. Returns the full path even
    when the parent directory does not yet exist; the caller is responsible
    for ``os.makedirs`` if appropriate.
    """
    here = os.path.abspath(__file__)
    project_root = here
    for _ in range(4):
        project_root = os.path.dirname(project_root)
    return os.path.join(project_root, rel_path)


def _save_annotation_as_xlsx(df, xlsx_path: str) -> bool:
    """Write the annotation DataFrame as a formatted xlsx.

    Rows are sorted by ``selection_rate`` descending so the most
    diagnostically informative CpGs appear first. The header row is bold +
    dark fill; numeric columns are right-aligned, boolean columns
    center-aligned, text columns left-aligned with wrap; ``selection_rate``
    is rendered as a percentage and ``ewas_min_p`` in scientific notation;
    alternating rows receive a subtle zebra fill; the UCSC URL becomes a
    clickable hyperlink; the header row is frozen and an auto-filter is
    enabled. Returns True on success, False if openpyxl is unavailable.
    """
    try:
        import openpyxl
        from openpyxl.styles import Font, Alignment, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError:
        logger.warning("  openpyxl not installed — Supplementary Table S2 "
                       "xlsx not written")
        return False

    # ── Sort: most diagnostically informative CpGs first ────────────────
    df = df.copy()
    if "selection_rate" in df.columns:
        df = df.sort_values("selection_rate", ascending=False,
                            na_position="last").reset_index(drop=True)

    # ── Column-type classification (drives alignment + number formats) ──
    numeric_cols = {"pos_0based", "selection_rate", "selection_count",
                    "distance_bp", "ewas_n_hits", "ewas_n_studies",
                    "ewas_min_p", "atlas_n_hits", "atlas_n_studies",
                    "pubmed_n_hits"}
    bool_cols    = {"in_88_regions", "has_plasma_evidence"}
    url_cols     = {"ucsc_url"}

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Panel CpG Annotation"

    headers = list(df.columns)
    ws.append(headers)

    # ── Header row styling ──────────────────────────────────────────────
    header_font  = Font(bold=True, color="FFFFFF", size=11)
    header_fill  = PatternFill(start_color="34495E", end_color="34495E",
                               fill_type="solid")
    header_align = Alignment(horizontal="center", vertical="center",
                             wrap_text=True)
    for c in range(1, len(headers) + 1):
        cell = ws.cell(row=1, column=c)
        cell.font      = header_font
        cell.fill      = header_fill
        cell.alignment = header_align

    # ── Data row styling ────────────────────────────────────────────────
    left_align   = Alignment(horizontal="left",   vertical="top", wrap_text=True)
    right_align  = Alignment(horizontal="right",  vertical="top")
    center_align = Alignment(horizontal="center", vertical="top")
    zebra_fill   = PatternFill(start_color="F4F6F8", end_color="F4F6F8",
                               fill_type="solid")
    link_font    = Font(color="0563C1", underline="single")

    for r_idx, row in enumerate(df.itertuples(index=False), start=2):
        for c_idx, (header, value) in enumerate(zip(headers, row), start=1):
            cell = ws.cell(row=r_idx, column=c_idx, value=value)

            if header in numeric_cols:
                cell.alignment = right_align
            elif header in bool_cols:
                cell.alignment = center_align
            else:
                cell.alignment = left_align

            if header == "selection_rate":
                cell.number_format = "0.00%"
            elif header == "ewas_min_p" and isinstance(value, (int, float)):
                cell.number_format = "0.00E+00"

            if header in url_cols and isinstance(value, str) and value.startswith("http"):
                cell.hyperlink = value
                cell.font      = link_font

            if r_idx % 2 == 0:
                cell.fill = zebra_fill

    # ── Column widths — measured from actual content (sampled rows) ─────
    sample_upper = min(ws.max_row + 1, 60)
    for c_idx, header in enumerate(headers, start=1):
        max_len = len(str(header))
        for r in range(2, sample_upper):
            v = ws.cell(r, c_idx).value
            if v is None:
                continue
            text = str(v)
            if header in url_cols:
                text = "ucsc.edu link"   # shorten URL footprint
            max_len = max(max_len, len(text))
        ws.column_dimensions[get_column_letter(c_idx)].width = min(45,
                                                                   max(10,
                                                                       max_len + 2))

    ws.auto_filter.ref       = ws.dimensions
    ws.freeze_panes          = "A2"
    ws.row_dimensions[1].height = 32

    os.makedirs(os.path.dirname(xlsx_path), exist_ok=True)
    wb.save(xlsx_path)
    return True


# ---------------------------------------------------------------------------
# Main entry: Step 18
# ---------------------------------------------------------------------------

def run_annotation(result_dir, config):
    """
    Step 17: Annotate panel CpG sites with gene context, CpG island status,
    EWAS Catalog + EWAS Atlas literature evidence, and PubMed hit count.

    Responses are cached in annotation_cache.json so re-runs skip network calls.

    Returns:
        ann_df:  pd.DataFrame  (saved as panel_annotation.csv)
        summary: dict          (aggregate counts for report)
    """
    ann_cfg       = (config or {}).get("annotation", {})
    ewas_win      = ann_cfg.get("ewas_window",   500)
    gene_win      = ann_cfg.get("gene_window",   2000)
    api_delay     = ann_cfg.get("api_delay",     0.25)
    do_atlas      = ann_cfg.get("ewas_atlas",    True)
    do_pubmed     = ann_cfg.get("pubmed_search", True)
    ncbi_key      = ann_cfg.get("ncbi_api_key",  "")
    pubmed_max    = ann_cfg.get("pubmed_max_hits", 20)

    # ── Load sites ──────────────────────────────────────────────────────────
    try:
        df = _load_sites(result_dir)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(str(exc))
        return pd.DataFrame(), {}

    sites   = list(df["feature"])
    n_total = len(sites)
    logger.info(f"Annotating {n_total} panel sites ...")

    # ── 88-region membership ─────────────────────────────────────────────────
    candidate_regions = (config or {}).get("data", {}).get("candidate_regions", "")
    try:
        from scripts.modules.feature_utils import _check_region_membership
        in_88 = _check_region_membership(sites, candidate_regions)
    except Exception:
        in_88 = np.zeros(n_total, dtype=bool)

    # ── CpG island BED ───────────────────────────────────────────────────────
    excel_path = (config or {}).get("data", {}).get("sample_info_excel", "data/x.xlsx")
    data_dir   = os.path.dirname(excel_path) if os.path.dirname(excel_path) else "data"
    cpgi_path  = ann_cfg.get("cpgi_bed", "") or _ensure_cpgi_bed(data_dir)
    cpgi_ivs   = _load_cpgi_intervals(cpgi_path)
    logger.info(f"  CpG island BED: {cpgi_path or 'unavailable (context = N/A)'}")

    # ── Cache ────────────────────────────────────────────────────────────────
    cache_path = os.path.join(result_dir, "annotation_cache.json")
    cache      = _load_cache(cache_path)
    n_cached   = sum(1 for s in sites if s in cache)
    logger.info(f"  API cache: {n_cached}/{n_total} sites already cached")

    # ── Annotate ─────────────────────────────────────────────────────────────
    rows = []
    for i, site in enumerate(sites):
        chrom, pos_str = site.split(":", 1)
        pos_0b = int(pos_str)
        pos_1b = pos_0b + 1          # convert 0-based BED → 1-based for APIs

        if (i + 1) % 5 == 0 or i == n_total - 1:
            logger.info(f"  [{i+1}/{n_total}] {site}")

        entry = cache.get(site, {})

        # Gene annotation
        if "gene" not in entry:
            entry["gene"] = _ensembl_nearest_gene(
                chrom, pos_1b, window=gene_win, delay=api_delay)
            cache[site] = entry

        # EWAS Catalog hits
        if "ewas" not in entry:
            entry["ewas"] = _ewas_query(
                chrom, pos_1b, window=ewas_win, delay=api_delay)
            cache[site] = entry

        # EWAS Atlas hits (CNCB)
        if do_atlas and "ewas_atlas" not in entry:
            entry["ewas_atlas"] = _ewas_atlas_query(
                chrom, pos_1b, window=ewas_win, delay=api_delay)
            cache[site] = entry

        gene           = entry["gene"]
        ewas_hits      = entry["ewas"]
        atlas_hits     = entry.get("ewas_atlas", []) if do_atlas else []

        # PubMed search (keyed on gene name to avoid per-site redundancy)
        gene_name = gene.get("gene_name", "")
        pubmed_key = f"pubmed:{gene_name}"
        if do_pubmed and gene_name and pubmed_key not in cache:
            pm_count, pm_ids = _pubmed_search(
                gene_name, ncbi_api_key=ncbi_key,
                max_hits=pubmed_max, delay=api_delay)
            cache[pubmed_key] = {"count": pm_count, "pmids": pm_ids}
        pm_entry  = cache.get(pubmed_key, {}) if do_pubmed and gene_name else {}
        pm_count  = pm_entry.get("count",  0)
        pm_pmids  = pm_entry.get("pmids",  [])

        cpg_ctx = (_cpg_context(chrom, pos_0b, cpgi_ivs)
                   if cpgi_ivs else "N/A")

        # EWAS Catalog fields
        traits  = list({h.get("trait", "")  for h in ewas_hits if h.get("trait")})
        pmids   = list({str(h.get("pmid","") or h.get("study_id",""))
                        for h in ewas_hits if h.get("pmid") or h.get("study_id")})
        p_vals  = [float(h["p_value"]) for h in ewas_hits
                   if h.get("p_value") not in (None, "")]
        min_p   = min(p_vals) if p_vals else None

        # EWAS Atlas fields
        atlas_traits = list({h.get("trait", "") or h.get("phenotype", "")
                             for h in atlas_hits
                             if h.get("trait") or h.get("phenotype")})
        atlas_pmids  = list({str(h.get("pmid", "") or h.get("pubmed_id", ""))
                             for h in atlas_hits
                             if h.get("pmid") or h.get("pubmed_id")})

        # Plasma evidence: keyword scan across both EWAS sources
        all_traits_text = " ".join(traits + atlas_traits)
        plasma_flag = _has_plasma_evidence(all_traits_text)

        ucsc_url = (f"https://genome.ucsc.edu/cgi-bin/hgTracks?"
                    f"db=hg38&position={chrom}%3A{pos_1b-500}-{pos_1b+500}")

        row = {
            "feature":        site,
            "chrom":          chrom,
            "pos_0based":     pos_0b,
            "in_88_regions":  bool(in_88[i]),
        }
        if "selection_rate" in df.columns:
            mask = df["feature"] == site
            row["selection_rate"]  = (df.loc[mask, "selection_rate"].iloc[0]
                                      if mask.any() else None)
            row["selection_count"] = (df.loc[mask, "selection_count"].iloc[0]
                                      if "selection_count" in df.columns and mask.any()
                                      else None)
        row.update({
            "nearest_gene":      gene.get("gene_name", ""),
            "gene_id":           gene.get("gene_id",   ""),
            "biotype":           gene.get("biotype",   ""),
            "distance_bp":       gene.get("distance_bp"),
            "region_type":       gene.get("region_type", ""),
            "cpg_context":       cpg_ctx,
            # EWAS Catalog
            "ewas_n_hits":       len(ewas_hits),
            "ewas_n_studies":    len(pmids),
            "ewas_traits":       "; ".join(traits[:8]),
            "ewas_min_p":        round(min_p, 8) if min_p is not None else None,
            # EWAS Atlas
            "atlas_n_hits":      len(atlas_hits),
            "atlas_n_studies":   len(atlas_pmids),
            "atlas_traits":      "; ".join(atlas_traits[:8]),
            # Plasma evidence flag
            "has_plasma_evidence": plasma_flag,
            # PubMed
            "pubmed_n_hits":     pm_count,
            "pubmed_pmids":      "; ".join(pm_pmids[:5]),
            "pubmed_query":      f'"{gene_name}" plasma methylation cancer' if gene_name else "",
            # UCSC browser link
            "ucsc_url":          ucsc_url,
        })
        rows.append(row)

    # Persist cache after all annotations
    _save_cache(cache, cache_path)

    ann_df   = pd.DataFrame(rows)
    out_path = os.path.join(result_dir, "panel_annotation.csv")
    ann_df.to_csv(out_path, index=False)
    logger.info(f"  Annotation saved: {out_path}")

    # Mirror the table into the manuscript docs tree as Supplementary Table S2
    # (formatted xlsx). Failures are non-fatal — the pipeline still produces
    # the canonical CSV above.
    try:
        s2_xlsx = _resolve_project_path("docs/iMeta_0.2/Table_S2_annotation.xlsx")
        if _save_annotation_as_xlsx(ann_df, s2_xlsx):
            logger.info(f"  Supplementary Table S2 saved: {s2_xlsx}")
    except Exception as exc:
        logger.warning(f"  Failed to write Supplementary Table S2 xlsx: {exc}")

    # ── Summary ──────────────────────────────────────────────────────────────
    def _cnt(col, val):
        return int((ann_df[col] == val).sum()) if col in ann_df.columns else 0

    summary = {
        "n_sites":              n_total,
        "n_in_88":              int(ann_df["in_88_regions"].sum()),
        "n_with_gene":          int((ann_df["nearest_gene"] != "").sum()),
        "n_promoter":           _cnt("region_type", "promoter"),
        "n_gene_body":          _cnt("region_type", "gene_body"),
        "n_intergenic":         _cnt("region_type", "intergenic"),
        "n_island":             _cnt("cpg_context", "Island"),
        "n_shore":              _cnt("cpg_context", "Shore"),
        "n_shelf":              _cnt("cpg_context", "Shelf"),
        "n_open_sea":           _cnt("cpg_context", "Open sea"),
        "n_with_ewas":          int((ann_df["ewas_n_hits"] > 0).sum()),
        "total_ewas_hits":      int(ann_df["ewas_n_hits"].sum()),
        "ewas_window_bp":       ewas_win,
        "n_with_atlas":         int((ann_df["atlas_n_hits"] > 0).sum()) if "atlas_n_hits" in ann_df else 0,
        "n_with_plasma_evidence": int(ann_df["has_plasma_evidence"].sum()) if "has_plasma_evidence" in ann_df else 0,
        "n_with_pubmed":        int((ann_df["pubmed_n_hits"] > 0).sum()) if "pubmed_n_hits" in ann_df else 0,
    }

    logger.info(
        f"  Summary: in_88={summary['n_in_88']}  "
        f"promoter={summary['n_promoter']}  gene_body={summary['n_gene_body']}  "
        f"intergenic={summary['n_intergenic']}  "
        f"island={summary['n_island']}  shore={summary['n_shore']}  "
        f"EWAS_hits={summary['n_with_ewas']}/{n_total}  "
        f"Atlas_hits={summary['n_with_atlas']}/{n_total}  "
        f"plasma_evidence={summary['n_with_plasma_evidence']}/{n_total}  "
        f"PubMed_hits={summary['n_with_pubmed']}/{n_total}"
    )

    return ann_df, summary
