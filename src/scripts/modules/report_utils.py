#!/usr/bin/env python3
"""
report_utils.py - Report generation module

Report structure:
  1.  Sample Overview (all samples — fastp QC)
  2.  Alignment & Deduplication
  3.  M-bias Check (inline SVG)
  4.  Sample Grouping (Benign vs Malignant)
  5.  Methylation Extraction & Candidate Regions
  6.  Methylation Distribution Visualization
  7.  Differential Methylation Analysis (DMC)
  8.  Sample Quality & Outlier Detection [Step 13, optional]
  9.  Feature Selection (LASSO + RF) [Step 14, optional]
  10. Diagnostic Classifiers (LR / SVM / RF / XGB / LGBM, LOO-CV) [Step 15, optional]
  11. Subtype Classifier (multi-class LOO-CV) [Step 16, optional]
  12. Temporal Split Validation (plasma only) [Step 17, optional]
  13. Feature Site Annotation (Ensembl + EWAS + CpG island) [Step 18, optional]
  14. Public Dataset Cross-validation (blood baseline + TCGA-OV) [Step 19, optional]
  15. Results Interpretation & Quality Assessment
  16. Analysis Pipeline
"""

import os
import re
import json
import glob
import shutil
import logging
from datetime import datetime

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────
#  Data-reading helpers
# ─────────────────────────────────────────────

def _parse_fastp_json(path):
    try:
        with open(path, encoding="utf-8") as f:
            d = json.load(f)
        sm = d.get("summary", {})
        bf = sm.get("before_filtering", {})
        af = sm.get("after_filtering", {})
        dup_rate = d.get("duplication", {}).get("rate", None)
        return {
            "total_reads_before": bf.get("total_reads", 0),
            "total_bases_before": bf.get("total_bases", 0),
            "q30_before":         bf.get("q30_rate", 0),
            "total_reads_after":  af.get("total_reads", 0),
            "total_bases_after":  af.get("total_bases", 0),
            "q30_after":          af.get("q30_rate", 0),
            "passed_filter_reads": d.get("filtering_result", {}).get("passed_filter_reads", 0),
            "low_quality_reads":   d.get("filtering_result", {}).get("low_quality_reads", 0),
            "too_short_reads":     d.get("filtering_result", {}).get("too_short_reads", 0),
            "dup_rate": dup_rate,
        }
    except Exception as e:
        logger.warning(f"Cannot parse fastp JSON {path}: {e}")
        return {}


def _parse_dedup_stats(path):
    result = {}
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
        m = re.search(r"Input Reads:\s*([\d,]+)", content)
        if m:
            result["input_reads"] = int(m.group(1).replace(",", ""))
        m = re.search(r"Read pairs:\s*([\d,]+)", content)
        if m:
            result["input_pairs"] = int(m.group(1).replace(",", ""))
        m = re.search(r"Number of reads out:\s*([\d,]+)", content)
        if m:
            result["output_reads"] = int(m.group(1).replace(",", ""))
        if "input_pairs" in result and "output_reads" in result and result["input_pairs"] > 0:
            result["dup_rate"] = 1.0 - result["output_reads"] / result["input_pairs"]
    except Exception as e:
        logger.warning(f"Cannot parse dedup stats {path}: {e}")
    return result


def _count_bedgraph_sites(path):
    try:
        count = 0
        meth_sum = 0.0
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    count += 1
                    try:
                        meth_sum += float(parts[3])
                    except ValueError:
                        pass
        avg_meth = meth_sum / count if count > 0 else 0
        return count, avg_meth
    except Exception:
        return 0, 0.0


def _parse_diff_results(path):
    rows = []
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 5:
                    rows.append({
                        "chrom": parts[0],
                        "pos":   int(parts[1]),
                        "meth1": float(parts[2]),
                        "meth2": float(parts[3]),
                        "diff":  float(parts[4]),
                    })
    except Exception as e:
        logger.warning(f"Cannot parse diff results {path}: {e}")
    return rows


def _fmt_num(n):
    if isinstance(n, int):
        return f"{n:,}"
    return str(n)

def _fmt_pct(r, digits=1):
    if r is None:
        return "N/A"
    return f"{r * 100:.{digits}f}%"

def _fmt_mb(b):
    return f"{b / 1024**2:.0f} MB" if b else "—"

def _fmt_gb(b):
    return f"{b / 1024**3:.2f} GB" if b else "—"

def _safe_pct(num, denom, digits=1):
    """Return percentage string; returns '—' when denominator is 0."""
    if not denom:
        return "—"
    return f"{num / denom * 100:.{digits}f}%"


def _copy_img(src_path, imgs_dir):
    """Copy *src_path* to *imgs_dir* and return a relative URL ``imgs/<fname>``.
    Returns None if the source file does not exist."""
    if not src_path or not os.path.exists(src_path):
        return None
    os.makedirs(imgs_dir, exist_ok=True)
    fname = os.path.basename(src_path)
    shutil.copy2(src_path, os.path.join(imgs_dir, fname))
    return f"imgs/{fname}"


def _img_tag(src_rel, alt):
    if not src_rel:
        return f"<p style='color:#aaa;font-size:.85em'>{alt} not yet generated — run --step visualize first</p>"
    return f"<img src='{src_rel}' alt='{alt}' style='max-width:100%;border-radius:6px;margin-top:6px'/>"


def _read_csv_safe(path):
    """Read CSV and return list[dict]; returns None if file missing or pandas unavailable."""
    if not _HAS_PANDAS or not path or not os.path.exists(path):
        return None
    try:
        import pandas as pd
        df = pd.read_csv(path)
        return df.to_dict(orient="records")
    except Exception:
        return None


def _pagination_widget(tbl_id: str, ctrl_id: str, n_total: int, page_size: int = 20):
    """Return (controls_html, js_html) for a paginated <table id=tbl_id>.

    Each <tr> must carry a data-row='index' attribute (0-based).
    The table must have a <tbody>.
    """
    n_pages = max(1, (n_total + page_size - 1) // page_size)
    controls = f"""
<div id="{ctrl_id}_bar" style="display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin:8px 0;font-size:.83em">
  <button id="{ctrl_id}_prev"
          style="padding:2px 10px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px"
          onclick="window.{ctrl_id}_show(window.{ctrl_id}_cur - 1)">&#8592; Prev</button>
  <span id="{ctrl_id}_pages" style="display:inline-flex;flex-wrap:wrap;gap:2px"></span>
  <button id="{ctrl_id}_next"
          style="padding:2px 10px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px"
          onclick="window.{ctrl_id}_show(window.{ctrl_id}_cur + 1)">Next &#8594;</button>
  <span id="{ctrl_id}_info" style="color:#666;margin-left:6px"></span>
</div>"""
    js = f"""
<script>
(function() {{
  var PAGE  = {page_size};
  var TOTAL = {n_total};
  var PAGES = {n_pages};
  var cur   = 0;
  function show(p) {{
    p = Math.max(0, Math.min(p, PAGES - 1));
    cur = p;
    window.{ctrl_id}_cur = cur;
    var rows = document.querySelectorAll('#{tbl_id} tbody tr');
    var start = cur * PAGE;
    rows.forEach(function(tr) {{
      var i = parseInt(tr.getAttribute('data-row'));
      tr.style.display = (i >= start && i < start + PAGE) ? '' : 'none';
    }});
    var info = document.getElementById('{ctrl_id}_info');
    if (info) info.textContent = 'Showing ' + (start+1) + '\u2013' + Math.min(start+PAGE, TOTAL) + ' of ' + TOTAL;
    var prev = document.getElementById('{ctrl_id}_prev');
    var next = document.getElementById('{ctrl_id}_next');
    if (prev) prev.disabled = (cur === 0);
    if (next) next.disabled = (cur === PAGES - 1);
    document.querySelectorAll('.{ctrl_id}_pg').forEach(function(btn) {{
      var active = parseInt(btn.getAttribute('data-p')) === cur;
      btn.style.fontWeight  = active ? 'bold'    : 'normal';
      btn.style.background  = active ? '#e3f0ff' : '#fff';
      btn.style.borderColor = active ? '#1a6fbb' : '#ccc';
    }});
  }}
  window.{ctrl_id}_show = show;
  window.{ctrl_id}_cur  = 0;
  document.addEventListener('DOMContentLoaded', function() {{
    var pgArea = document.getElementById('{ctrl_id}_pages');
    if (pgArea) {{
      for (var p = 0; p < PAGES; p++) {{
        (function(pp) {{
          var b = document.createElement('button');
          b.textContent = pp + 1;
          b.setAttribute('data-p', pp);
          b.className = '{ctrl_id}_pg';
          b.style.cssText = 'margin:0 1px;padding:2px 7px;cursor:pointer;border:1px solid #ccc;background:#fff;border-radius:3px';
          b.onclick = function() {{ show(pp); }};
          pgArea.appendChild(b);
        }})(p);
      }}
    }}
    show(0);
  }});
}})();
</script>"""
    return controls, js


def _build_outlier_section(result_dir, imgs_dir=None):
    """
    Read outlier_report.csv and build the HTML section + JSON data.
    Returns (html_str, data_dict). Both empty/None when file is absent.
    The sample table uses JS pagination (PAGE_SIZE rows per page).
    """
    path = os.path.join(result_dir, "outlier_report.csv")
    records = _read_csv_safe(path)
    if records is None:
        return (
            '<div class="section"><h2>Step 13. Sample Quality &amp; Outlier Detection</h2>'
            '<p style="color:#aaa;font-size:.9em">Step 13 not yet run — '
            'execute <code>--step outlier_detect</code> to populate this section.</p></div>',
            None,
        )

    PAGE_SIZE = 20
    n_total   = len(records)
    n_suspect = sum(1 for r in records if str(r.get("suspect", "")).lower() == "true")

    # Build one <tr> per record; embed row index as data-row attribute for JS paging
    rows_html = ""
    for idx, r in enumerate(records):
        is_suspect = str(r.get("suspect", "")).lower() == "true"
        bg = "background:#fff3cd" if is_suspect else ""
        badge = ("<span style='color:#c0392b;font-weight:bold'>&#9888; suspect</span>"
                 if is_suspect else "<span style='color:#aaa'>ok</span>")
        flag_cov = "&#10003;" if str(r.get("flag_coverage",  "")).lower() == "true" else "&mdash;"
        flag_pca = "&#10003;" if str(r.get("flag_pca",       "")).lower() == "true" else "&mdash;"
        flag_if  = "&#10003;" if str(r.get("flag_isolation", "")).lower() == "true" else "&mdash;"
        try:
            pca_z = f"{float(r.get('pca_z_score', '')):.2f}"
        except (TypeError, ValueError):
            pca_z = "&mdash;"
        try:
            nan_frac = f"{float(r.get('nan_frac', ''))*100:.1f}%"
        except (TypeError, ValueError):
            nan_frac = str(r.get("nan_frac", ""))
        rows_html += (
            f"<tr data-row='{idx}' style='{bg}'>"
            f"<td>{r.get('sample_id','')}</td>"
            f"<td>{r.get('group','')}</td>"
            f"<td style='text-align:right'>{nan_frac}</td>"
            f"<td style='text-align:center'>{flag_cov}</td>"
            f"<td style='text-align:right'>{pca_z}</td>"
            f"<td style='text-align:center'>{flag_pca}</td>"
            f"<td style='text-align:center'>{flag_if}</td>"
            f"<td style='text-align:center'>{r.get('n_flags',0)}</td>"
            f"<td style='text-align:center'>{badge}</td>"
            f"</tr>\n"
        )

    if n_suspect > 0:
        suspect_note = (
            f"<p style='color:#c0392b;font-size:.85em;margin-top:8px'>"
            f"&#9888; <strong>{n_suspect}</strong> suspect sample(s) detected. "
            f"After reviewing <code>outlier_report.csv</code>, add confirmed IDs to "
            f"<code>config.toml</code> <code>[ml] exclude_samples = [...]</code> "
            f"and re-run the <code>classify</code> step.</p>"
        )
    else:
        suspect_note = ("<p style='color:#22863a;font-size:.85em;margin-top:8px'>"
                        "&#10003; No suspect samples detected.</p>")

    # PCA scatter plot (Figure S2)
    pca_img_src = None
    if imgs_dir:
        pca_img_src = _copy_img(os.path.join(result_dir, "pca_outliers.png"), imgs_dir)
    pca_img_html = (
        f"<div style='margin:12px 0'><h3 style='font-size:.95em;margin-bottom:4px'>"
        f"Figure S2 — PCA Scatter (PC1 vs PC2)</h3>"
        f"{_img_tag(pca_img_src, 'PCA outlier scatter')}</div>"
        if pca_img_src else
        "<p style='color:#aaa;font-size:.85em;margin:8px 0'>"
        "PCA scatter plot not yet generated — run <code>--step outlier_detect</code>.</p>"
    )

    # Unique JS identifier so multiple reports on the same page don't conflict
    tbl_id  = "outlier_tbl"
    ctrl_id = "outlier_ctrl"
    controls, js = _pagination_widget(tbl_id, ctrl_id, n_total, PAGE_SIZE)

    html = f"""
<div class="section">
<h2>Step 13. Sample Quality &amp; Outlier Detection</h2>
<p style="color:#666;font-size:.85em">
  Three complementary methods flag unusual samples (flagged by &ge;2 methods = suspect):
  (1) <strong>Coverage</strong> &mdash; fraction of CpG sites with no read depth;
  (2) <strong>PCA within-group distance z-score</strong>;
  (3) <strong>Isolation Forest</strong> per-group anomaly score.
  Suspect samples are not automatically removed; use <code>exclude_samples</code> in
  <code>config.toml [ml]</code> to exclude confirmed outliers before re-running
  <code>classify</code>.
</p>
<table>
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Total samples analysed</td><td>{n_total}</td></tr>
  <tr><td>Suspect (&ge;2 flags)</td>
      <td><strong style="color:{'#c0392b' if n_suspect else '#22863a'}">{n_suspect}</strong></td></tr>
</table>
{suspect_note}
{pca_img_html}
<details style="margin-top:12px" {'open' if n_suspect > 0 else ''}>
  <summary style="cursor:pointer;color:#1a6fbb;font-size:.88em">
    Show all {n_total} samples &mdash; paginated, {PAGE_SIZE} per page (yellow = suspect)
  </summary>
  <div style="overflow-x:auto;margin-top:8px">
  {controls}
  <table id="{tbl_id}" style="font-size:.82em">
    <thead>
      <tr>
        <th>Sample</th><th>Group</th>
        <th>NaN frac</th><th>Flag (cov)</th>
        <th>PCA z</th><th>Flag (PCA)</th>
        <th>Flag (IF)</th><th># Flags</th><th>Status</th>
      </tr>
    </thead>
    <tbody>
    {rows_html}
    </tbody>
  </table>
  </div>
</details>
<p style="color:#666;font-size:.83em;margin-top:6px">
  Full report: <code>outlier_report.csv</code> &nbsp;|&nbsp;
  Thresholds: <code>outlier_coverage_min</code>, <code>outlier_pca_z</code> in
  <code>config.toml [ml]</code>
</p>
{js}
</div>"""

    data = {
        "n_samples": n_total,
        "n_suspect": n_suspect,
        "suspect_ids": [r.get("sample_id") for r in records
                        if str(r.get("suspect","")).lower() == "true"],
    }
    return html, data


def _build_annotation_section(result_dir):
    """
    Read panel_annotation.csv and build the feature annotation HTML section.
    Returns (html_str, data_dict).  Both empty/None when file is absent.
    Paginated table: PAGE_SIZE rows per page (same JS pattern as outlier section).
    """
    path    = os.path.join(result_dir, "panel_annotation.csv")
    records = _read_csv_safe(path)
    if records is None:
        return (
            '<div class="section"><h2>Step 18. Feature Site Annotation</h2>'
            '<p style="color:#aaa;font-size:.9em">Step 18 not yet run — '
            'execute <code>--step annotate</code> to populate this section.</p></div>',
            None,
        )

    PAGE_SIZE = 20
    n_total   = len(records)

    # ── Summary counts ────────────────────────────────────────────────────
    def _cnt(col, val):
        return sum(1 for r in records if str(r.get(col, "")).strip() == val)

    n_promoter  = _cnt("region_type", "promoter")
    n_gene_body = _cnt("region_type", "gene_body")
    n_intergenic= _cnt("region_type", "intergenic")
    n_island    = _cnt("cpg_context", "Island")
    n_shore     = _cnt("cpg_context", "Shore")
    n_shelf     = _cnt("cpg_context", "Shelf")
    n_open_sea  = _cnt("cpg_context", "Open sea")
    n_in_88     = sum(1 for r in records
                      if str(r.get("in_88_regions","")).lower() == "true")
    n_ewas      = sum(1 for r in records if int(r.get("ewas_n_hits", 0) or 0) > 0)
    total_hits  = sum(int(r.get("ewas_n_hits", 0) or 0) for r in records)

    # Sort by LOO-CV selection rate descending; NaN → treated as 0 (not in freq log)
    def _sort_key(r):
        try:
            import math
            v = float(r.get("selection_rate", 0) or 0)
            return 0.0 if math.isnan(v) else -v
        except (TypeError, ValueError):
            return 0.0
    records = sorted(records, key=_sort_key)

    # ── Table rows ────────────────────────────────────────────────────────
    region_colors = {
        "promoter":   "#e8f5e9",
        "gene_body":  "#e3f2fd",
        "intergenic": "",
    }
    cpg_badges = {
        "Island":   ("CpG island",  "#1b5e20", "#c8e6c9"),
        "Shore":    ("Shore",        "#1a237e", "#c5cae9"),
        "Shelf":    ("Shelf",        "#4a148c", "#e1bee7"),
        "Open sea": ("Open sea",     "#555",    "#eeeeee"),
        "N/A":      ("N/A",          "#555",    "#eeeeee"),
    }
    rows_html = ""
    def _s(val, default=""):
        """Convert a potentially NaN/float cell value to a clean string."""
        v = str(val).strip() if val is not None else ""
        return default if v.lower() in ("nan", "none", "") else v

    for idx, r in enumerate(records):
        in88   = _s(r.get("in_88_regions")).lower() == "true"
        site   = _s(r.get("feature"))
        ucsc   = _s(r.get("ucsc_url"))
        site_link = (f"<a href='{ucsc}' target='_blank' "
                     f"style='font-family:monospace;color:#1a6fbb'>{site}</a>"
                     if ucsc else f"<code>{site}</code>")
        marker = "&#9733;" if in88 else "&#9675;"   # ★ / ○
        rate   = r.get("selection_rate", "")
        try:
            import math
            v = float(r.get("selection_rate", ""))
            rate = "&mdash;" if math.isnan(v) else f"{v*100:.1f}%"
        except (TypeError, ValueError):
            rate = "&mdash;"
        gene     = _s(r.get("nearest_gene")) or "&mdash;"
        dist     = r.get("distance_bp", "")
        try:
            dist = f"{int(float(dist))} bp" if float(dist) > 0 else "overlap"
        except (TypeError, ValueError):
            dist = "&mdash;"
        rtype  = _s(r.get("region_type"))
        bg     = region_colors.get(rtype, "")
        cpg    = _s(r.get("cpg_context"), "N/A")
        lbl, fc, bc = cpg_badges.get(cpg, cpg_badges["N/A"])
        cpg_badge = (f"<span style='background:{bc};color:{fc};"
                     f"padding:1px 5px;border-radius:3px;font-size:.8em'>{lbl}</span>")
        n_hits = int(r.get("ewas_n_hits", 0) or 0)
        traits = _s(r.get("ewas_traits"))
        if n_hits > 0 and traits:
            safe_traits = traits.replace("'", "&apos;")
            ewas_cell = (f"<td style='text-align:center'>"
                         f"<span title='{safe_traits}' "
                         f"style='cursor:help;color:#c0392b;font-weight:bold'>"
                         f"{n_hits} &#9432;</span></td>")
        elif n_hits > 0:
            ewas_cell = f"<td style='text-align:center'>{n_hits}</td>"
        else:
            ewas_cell = "<td style='text-align:center;color:#aaa'>&mdash;</td>"
        biotype = _s(r.get("biotype"))
        rows_html += (
            f"<tr data-row='{idx}' style='background:{bg}'>"
            f"<td style='font-size:.8em'>{site_link}</td>"
            f"<td style='text-align:center'>{marker}</td>"
            f"<td style='text-align:center'>{rate}</td>"
            f"<td><strong>{gene}</strong>"
            f"{'<br><span style=color:#888;font-size:.78em>' + biotype + '</span>' if biotype else ''}"
            f"</td>"
            f"<td style='text-align:right;font-size:.85em'>{dist}</td>"
            f"<td style='font-size:.85em'>{rtype or '&mdash;'}</td>"
            f"<td>{cpg_badge}</td>"
            f"{ewas_cell}"
            f"</tr>\n"
        )

    # ── Pagination JS (same pattern as outlier section) ───────────────────
    TBL  = "annot_tbl"
    CTRL = "annot_ctrl"
    controls, js = _pagination_widget(TBL, CTRL, n_total, PAGE_SIZE)

    html = f"""
<div class="section">
<h2>Step 18. Feature Site Annotation</h2>
<p style="color:#666;font-size:.85em">
  Each panel CpG site is annotated with the nearest gene (Ensembl REST API),
  CpG island context (UCSC hg38), and associations from the EWAS Catalog
  (&plusmn;500 bp window).
  &#9733; = within 88 candidate tissue regions; &#9675; = outside.
  Hover over EWAS hit counts to see associated traits.
  Coordinates are 0-based BED (MethylDackel); UCSC links use 1-based display.
</p>
<table style="margin-bottom:10px">
  <tr>
    <th colspan="2">Region type</th>
    <th colspan="4">CpG island context</th>
    <th colspan="2">EWAS evidence</th>
  </tr>
  <tr>
    <td>Promoter</td><td><strong>{n_promoter}</strong></td>
    <td>Island</td><td><strong>{n_island}</strong></td>
    <td>Shore</td><td><strong>{n_shore}</strong></td>
    <td>Sites with &ge;1 hit</td><td><strong style="color:#c0392b">{n_ewas}</strong></td>
  </tr>
  <tr>
    <td>Gene body</td><td><strong>{n_gene_body}</strong></td>
    <td>Shelf</td><td><strong>{n_shelf}</strong></td>
    <td>Open sea</td><td><strong>{n_open_sea}</strong></td>
    <td>Total EWAS hits</td><td><strong>{total_hits}</strong></td>
  </tr>
  <tr>
    <td>Intergenic</td><td><strong>{n_intergenic}</strong></td>
    <td colspan="4"></td>
    <td>In 88 candidate regions</td><td><strong style="color:#1b5e20">{n_in_88}</strong></td>
  </tr>
</table>
{controls}
<div style="overflow-x:auto">
<table id="{TBL}" style="font-size:.82em">
  <thead>
    <tr>
      <th>Site (UCSC &rarr;)</th>
      <th>88?</th>
      <th>Sel. rate (%)</th>
      <th>Nearest gene</th>
      <th>Distance</th>
      <th>Region</th>
      <th>CpG context</th>
      <th>EWAS hits</th>
    </tr>
  </thead>
  <tbody>
  {rows_html}
  </tbody>
</table>
</div>
<p style="color:#666;font-size:.83em;margin-top:6px">
  Full table: <code>panel_annotation.csv</code> &nbsp;|&nbsp;
  API responses cached in <code>annotation_cache.json</code> (re-runs are instant) &nbsp;|&nbsp;
  Configurable in <code>config.toml [annotation]</code>
</p>
{js}
</div>"""

    data = {
        "n_sites":        n_total,
        "n_in_88":        n_in_88,
        "n_promoter":     n_promoter,
        "n_gene_body":    n_gene_body,
        "n_intergenic":   n_intergenic,
        "n_island":       n_island,
        "n_with_ewas":    n_ewas,
        "total_ewas_hits": total_hits,
    }
    return html, data


def _build_external_validation_section(result_dir):
    """
    Read external_validation.csv (produced by external_validate.py) and build
    the public-dataset cross-validation HTML section.
    Returns (html_str, data_dict). Both empty/None when file is absent.
    Joins selection_rate from panel_annotation.csv (if present).
    Table is sortable by clicking column headers; default sort: selection_rate desc.
    """
    path    = os.path.join(result_dir, "external_validation.csv")
    records = _read_csv_safe(path)
    if records is None:
        return (
            '<div class="section"><h2>Step 19. Public Dataset Cross-validation (optional)</h2>'
            '<p style="color:#aaa;font-size:.9em">Step 19 not yet run — '
            'execute <code>--step external_validate</code> to populate this section.</p></div>',
            None,
        )

    # ── Join selection_rate from panel_annotation.csv ────────────────────────
    annot_path = os.path.join(result_dir, "panel_annotation.csv")
    annot_records = _read_csv_safe(annot_path)
    sel_rate_map = {}
    if annot_records:
        import math
        for ar in annot_records:
            feat = ar.get("feature", "")
            try:
                v = float(ar.get("selection_rate", "") or "")
                sel_rate_map[feat] = 0.0 if math.isnan(v) else v
            except (TypeError, ValueError):
                sel_rate_map[feat] = 0.0
    for r in records:
        r["selection_rate"] = sel_rate_map.get(r.get("feature", ""), None)

    # ── Default sort: selection_rate descending (NaN → bottom) ───────────────
    def _sort_key(r):
        v = r.get("selection_rate")
        if v is None:
            return -1.0
        try:
            return -float(v)
        except (TypeError, ValueError):
            return -1.0
    records = sorted(records, key=_sort_key)

    n_total      = len(records)
    n_matched    = sum(1 for r in records if str(r.get("probe_matched", "")).lower() == "true")
    n_blood_pass = sum(1 for r in records
                       if str(r.get("blood_pass",     "")).lower() == "true"
                       and str(r.get("probe_matched", "")).lower() == "true")
    n_blood_fail = sum(1 for r in records
                       if str(r.get("blood_pass",     "")).lower() == "false"
                       and str(r.get("probe_matched", "")).lower() == "true")
    n_tcga = sum(1 for r in records
                 if str(r.get("tcga_n_samples", "0")).strip() not in ("0", "nan", "", "None"))

    def _fv(val, fmt=".3f"):
        try:
            return f"{float(val):{fmt}}"
        except (TypeError, ValueError):
            return "&mdash;"

    PAGE_SIZE = 20
    TBL  = "extval_tbl"
    CTRL = "extval_ctrl"

    # ── Table rows ─────────────────────────────────────────────────────────────
    # Embed numeric sort keys as data-* for the JS sorter
    rows_html = ""
    for idx, r in enumerate(records):
        matched = str(r.get("probe_matched", "")).lower() == "true"
        bpass   = str(r.get("blood_pass",   "")).lower() == "true"
        bg = "background:#fff3cd" if (matched and not bpass) else ""
        if not matched:
            blood_badge = "<span style='color:#aaa'>no probe</span>"
        elif bpass:
            blood_badge = "<span style='color:#22863a;font-weight:bold'>&#10003; pass</span>"
        else:
            blood_badge = "<span style='color:#c0392b;font-weight:bold'>&#9888; fail</span>"

        # numeric values for sort; use -1 as sentinel for missing
        sr_raw = r.get("selection_rate")
        try:
            sr_num = float(sr_raw) if sr_raw is not None else -1.0
            sr_disp = f"{sr_num*100:.1f}%" if sr_num >= 0 else "&mdash;"
        except (TypeError, ValueError):
            sr_num, sr_disp = -1.0, "&mdash;"

        def _num(v, default=-1.0):
            try:
                return float(v)
            except (TypeError, ValueError):
                return default

        rows_html += (
            f"<tr data-row='{idx}' style='{bg}'"
            f" data-sel='{sr_num}'"
            f" data-site='{str(r.get('feature','') or '').lower()}'"
            f" data-cg='{str(r.get('cg_probe','') or '').lower()}'"
            f" data-bm='{_num(r.get('blood_mean_beta'))}'"
            f" data-bmed='{_num(r.get('blood_median_beta'))}'"
            f" data-bp='{'1' if bpass else ('0' if matched else '-1')}'"
            f" data-tm='{_num(r.get('tcga_mean_beta'))}'"
            f" data-tn='{_num(r.get('tcga_n_samples',0), 0)}'>"
            f"<td style='font-family:monospace;font-size:.82em'>{r.get('feature', '')}</td>"
            f"<td style='font-family:monospace;font-size:.82em'>{r.get('cg_probe', '') or '&mdash;'}</td>"
            f"<td style='text-align:right'><strong>{sr_disp}</strong></td>"
            f"<td style='text-align:right'>{_fv(r.get('blood_mean_beta'))}</td>"
            f"<td style='text-align:right'>{_fv(r.get('blood_median_beta'))}</td>"
            f"<td style='text-align:center'>{blood_badge}</td>"
            f"<td style='text-align:right'>{_fv(r.get('tcga_mean_beta'))}</td>"
            f"<td style='text-align:right'>{r.get('tcga_n_samples', '') or '&mdash;'}</td>"
            f"</tr>\n"
        )

    controls_extval, js_pagination = _pagination_widget(TBL, CTRL, n_total, PAGE_SIZE)

    # ── Blood-filter panel summary (written by external_validate after panel_cpg.csv exists) ─
    bf_path = os.path.join(result_dir, "blood_filter_summary.csv")
    bf_records = _read_csv_safe(bf_path)
    if bf_records:
        bfr = bf_records[0]
        try:
            bf_removed  = int(float(bfr.get("blood_fail_removed", 0)))
            bf_retained = int(float(bfr.get("panel_retained", 0)))
            bf_total    = int(float(bfr.get("panel_total", 0)))
            bf_thresh   = float(bfr.get("blood_threshold", 0.10))
        except (TypeError, ValueError):
            bf_removed = bf_retained = bf_total = 0
            bf_thresh = 0.10
        if bf_removed > 0:
            blood_filter_note = (
                f"<div style='background:#fff3cd;border:1px solid #f0ad4e;border-radius:4px;"
                f"padding:8px 12px;margin-top:10px;font-size:.87em'>"
                f"<strong>&#9888; Blood-filter applied (Step 19 &rarr; re-run Step 15):</strong> "
                f"<strong>{bf_removed}</strong> site(s) with blood mean &beta; &ge; {bf_thresh} "
                f"removed from <code>panel_cpg.csv</code> &rarr; "
                f"<code>panel_cpg_blood_filtered.csv</code> "
                f"(<strong>{bf_retained}/{bf_total}</strong> sites retained). "
                f"The blood-filtered panel is automatically used by <strong>Step 17</strong> "
                f"(temporal validation). To apply it to <strong>Step 15</strong> LOO-CV results, "
                f"re-run: <code>--step classify</code>."
                f"</div>"
            )
        else:
            blood_filter_note = (
                "<p style='color:#22863a;font-size:.85em;margin-top:8px'>"
                "&#10003; Blood-filter: no sites removed — all probe-matched sites pass "
                f"blood baseline threshold (&beta; &lt; {bf_thresh}).</p>"
            )
    else:
        blood_filter_note = ""

    if n_blood_fail > 0:
        blood_note = (
            f"<p style='color:#c0392b;font-size:.85em;margin-top:8px'>"
            f"&#9888; <strong>{n_blood_fail}</strong> site(s) fail blood baseline check "
            f"(blood_mean_beta &ge; 0.10, highlighted yellow). "
            f"These sites are potentially confounded by leukocyte DNA in plasma cfDNA — "
            f"they have been removed from the classifier panel "
            f"(see <code>panel_cpg_blood_filtered.csv</code>).</p>"
        )
    else:
        blood_note = ("<p style='color:#22863a;font-size:.85em;margin-top:8px'>"
                      "&#10003; All matched sites pass blood baseline check (mean &beta; &lt; 0.10).</p>")

    # ── Sort JS ─────────────────────────────────────────────────────────────────
    # Columns: 0=site(str), 1=cg(str), 2=sel_rate(num desc), 3=blood_mean(num),
    #          4=blood_med(num), 5=blood_pass(num desc), 6=tcga_mean(num), 7=tcga_n(num)
    js_sort = f"""
<script>
(function() {{
  var TBL   = document.getElementById ? null : null;   // resolved at runtime
  var sortState = {{ col: 2, asc: false }};   // default: sel_rate descending
  var DATA_ATTRS = ['data-site','data-cg','data-sel','data-bm','data-bmed','data-bp','data-tm','data-tn'];
  var STR_COLS   = new Set([0, 1]);

  function sortTable(col) {{
    TBL = document.getElementById('{TBL}');
    if (!TBL) return;
    var asc = (sortState.col === col) ? !sortState.asc : (col === 2 ? false : true);
    sortState = {{ col: col, asc: asc }};

    var tbody = TBL.querySelector('tbody');
    var rows  = Array.from(tbody.querySelectorAll('tr'));
    rows.sort(function(a, b) {{
      var av = a.getAttribute(DATA_ATTRS[col]);
      var bv = b.getAttribute(DATA_ATTRS[col]);
      var cmp;
      if (STR_COLS.has(col)) {{
        cmp = (av || '').localeCompare(bv || '');
      }} else {{
        var an = parseFloat(av), bn = parseFloat(bv);
        var aM = isNaN(an) || an < -0.5, bM = isNaN(bn) || bn < -0.5;
        if (aM && bM) cmp = 0;
        else if (aM)  cmp = 1;   // missing always last
        else if (bM)  cmp = -1;
        else          cmp = an - bn;
      }}
      return asc ? cmp : -cmp;
    }});
    // Re-assign data-row after sort so pagination still works
    rows.forEach(function(tr, i) {{
      tr.setAttribute('data-row', i);
      tbody.appendChild(tr);
    }});
    // Update header arrows
    document.querySelectorAll('#{TBL} thead th').forEach(function(th, i) {{
      th.style.cursor = 'pointer';
      var arrow = th.querySelector('.sort-arrow');
      if (arrow) {{
        arrow.textContent = (i === col) ? (asc ? ' ▲' : ' ▼') : ' ⇅';
      }}
    }});
    // Reset pagination to first page
    if (window.{CTRL}_show) window.{CTRL}_show(0);
  }}

  document.addEventListener('DOMContentLoaded', function() {{
    var ths = document.querySelectorAll('#{TBL} thead th');
    ths.forEach(function(th, i) {{
      var span = document.createElement('span');
      span.className = 'sort-arrow';
      span.style.cssText = 'font-size:.75em;margin-left:3px;color:#888';
      span.textContent = (i === 2) ? ' ▼' : ' ⇅';
      th.appendChild(span);
      th.style.cursor = 'pointer';
      th.onclick = function() {{ sortTable(i); }};
    }});
  }});
}})();
</script>"""

    html = f"""
<div class="section">
<h2>Step 19. Public Dataset Cross-validation (optional)</h2>
<p style="color:#666;font-size:.85em">
  Panel sites mapped to Illumina 450K probes via hg38&rarr;hg19 liftOver (&plusmn;2 bp) and
  cross-referenced against two public reference datasets:
  <strong>GSE40279</strong> (whole blood, n=656) for blood baseline, and
  <strong>TCGA-OV HumanMethylation450</strong> (&sim;590 tumour samples) for tumour signal confirmation.
  Blood baseline threshold: mean &beta; &lt; 0.10.
  Table sorted by LOO-CV feature selection rate by default; click any header to re-sort.
</p>
<table style="margin-bottom:10px">
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Panel sites</td><td>{n_total}</td></tr>
  <tr><td>Sites matched to 450K probe</td><td>{n_matched} ({_safe_pct(n_matched, n_total)})</td></tr>
  <tr><td>Blood baseline pass (mean &beta; &lt; 0.10)</td>
      <td><strong style="color:#22863a">{n_blood_pass}</strong></td></tr>
  <tr><td>Blood baseline fail (mean &beta; &ge; 0.10)</td>
      <td><strong style="color:{'#c0392b' if n_blood_fail else '#22863a'}">{n_blood_fail}</strong></td></tr>
  <tr><td>Sites with TCGA-OV coverage</td><td>{n_tcga}</td></tr>
</table>
{blood_note}
{blood_filter_note}
{controls_extval}
<div style="overflow-x:auto;margin-top:4px">
<table id="{TBL}" style="font-size:.82em">
  <thead>
    <tr>
      <th>Panel site</th>
      <th>cg probe</th>
      <th>Sel. rate</th>
      <th>Blood mean &beta;</th>
      <th>Blood median &beta;</th>
      <th>Blood baseline</th>
      <th>TCGA-OV mean &beta;</th>
      <th>TCGA-OV n</th>
    </tr>
  </thead>
  <tbody>
  {rows_html}
  </tbody>
</table>
</div>
<p style="color:#666;font-size:.83em;margin-top:6px">
  Yellow = blood baseline fail. Sel. rate = LOO-CV feature selection frequency (from panel_annotation.csv).
  Full data: <code>external_validation.csv</code> &nbsp;|&nbsp;
  Coordinate mapping: panel (hg38) &rarr; liftOver &rarr; hg19 &rarr; GPL13534 (&plusmn;2 bp).
</p>
{js_pagination}
{js_sort}
</div>"""

    data = {
        "n_total":      n_total,
        "n_matched":    n_matched,
        "n_blood_pass": n_blood_pass,
        "n_blood_fail": n_blood_fail,
        "n_tcga":       n_tcga,
    }
    return html, data


def _build_ovr_subsection(ovr_records, roc_src):
    """Return HTML fragment for the One-vs-Rest binary classifier subsection."""
    if not ovr_records:
        return ""

    # Group rows by subtype
    from collections import OrderedDict
    groups = OrderedDict()
    for r in ovr_records:
        st = r.get("subtype", "")
        groups.setdefault(st, []).append(r)

    rows = ""
    for subtype, recs in groups.items():
        first = True
        for r in recs:
            try:
                auc  = f"{float(r.get('auc',  0)):.3f}"
                acc  = f"{float(r.get('accuracy',   0)):.3f}"
                sens = f"{float(r.get('sensitivity', 0)):.3f}"
                spec = f"{float(r.get('specificity', 0)):.3f}"
            except (TypeError, ValueError):
                auc = acc = sens = spec = "&mdash;"
            st_cell = f"<td rowspan='{len(recs)}'>{subtype}</td>" if first else ""
            rows += (
                f"<tr>{st_cell}<td>{r.get('classifier','')}</td>"
                f"<td>{auc}</td><td>{acc}</td><td>{sens}</td><td>{spec}</td>"
                f"<td>{r.get('n_positive','')}</td><td>{r.get('n_negative','')}</td></tr>\n"
            )
            first = False

    return f"""
<hr style="margin:20px 0;border:none;border-top:1px solid #e0e0e0">
<h3 style="margin-top:0">11b. One-vs-Rest Binary Classifiers (LOO-CV)</h3>
<p style="color:#666;font-size:.85em">
  For each subtype, binary re-labelling: that subtype → positive, all others → negative.
  Same panel features and LOO-CV protocol as the multi-class task above.
  Sensitivity = subtype recall; Specificity = other-class recall.
</p>
<table>
  <tr><th>Subtype</th><th>Classifier</th><th>AUC</th><th>Accuracy</th>
      <th>Sensitivity</th><th>Specificity</th><th>n Positive</th><th>n Negative</th></tr>
  {rows}
</table>
<div style="margin-top:14px"><h3>ROC Curves (all subtypes)</h3>{_img_tag(roc_src, 'OvR ROC')}</div>
<p style="color:#666;font-size:.83em">
  Outputs: <code>ovr_summary.csv</code>, <code>ovr_predictions.csv</code>,
  per-subtype confusion matrices: <code>subtype_confusion_OvR_&lt;clf&gt;_&lt;subtype&gt;.png</code>
</p>"""


def _build_classifier_section(result_dir, imgs_dir):
    """
    Read Step 13-16 output files and build HTML fragments + JSON data.
    All sections are optional: returns empty string / None when files are absent.
    Returns (sec8_html, sec_outlier_html, sec10_html, sec11_html, sec12_html, clf_data_dict)
    """
    # ── Section 8: Feature Selection (Step 14) ──
    panel_path    = os.path.join(result_dir, "panel_cpg.csv")
    feat_imp_path = os.path.join(result_dir, "feature_importance.csv")
    feat_records  = _read_csv_safe(feat_imp_path)
    panel_records = _read_csv_safe(panel_path)

    if panel_records is not None and feat_records is not None:
        n_panel   = len(panel_records)
        n_total   = len(feat_records)
        lasso_sel = sum(1 for r in feat_records if r.get("lasso_selected"))
        rf_sel    = sum(1 for r in feat_records if r.get("rf_selected"))

        # ── Wilcoxon sample-level filter stats ──
        wilcox_records = _read_csv_safe(os.path.join(result_dir, "wilcoxon_stats.csv"))
        if wilcox_records:
            n_candidates = len(wilcox_records)
            n_wpass   = sum(1 for r in wilcox_records if r.get("pass_both"))
            n_fail_p  = sum(1 for r in wilcox_records if not r.get("pass_p") and r.get("pass_delta"))
            n_fail_d  = sum(1 for r in wilcox_records if r.get("pass_p") and not r.get("pass_delta"))
            n_fail_b  = sum(1 for r in wilcox_records if not r.get("pass_p") and not r.get("pass_delta"))
            all_deltas = [float(r.get("sample_delta", 0)) for r in wilcox_records]
            all_pvals  = [float(r.get("wilcox_pval",  1)) for r in wilcox_records]
            delta_mean   = sum(all_deltas) / len(all_deltas) if all_deltas else 0
            delta_max    = max(all_deltas) if all_deltas else 0
            delta_median = sorted(all_deltas)[len(all_deltas) // 2] if all_deltas else 0
            n_p05  = sum(1 for p in all_pvals if p < 0.05)
            n_p10  = sum(1 for p in all_pvals if p < 0.10)
            n_p20  = sum(1 for p in all_pvals if p < 0.20)

            # All-sites detail table (collapsed, sorted by delta desc)
            WTBL  = "wilcox_tbl"
            WCTRL = "wilcox_ctrl"
            all_rows = "".join(
                f"<tr data-row='{idx}' style='{'background:#e8f5e9' if r.get('pass_both') else ''}'>"
                f"<td>{r['feature']}</td>"
                f"<td style='text-align:right'>{float(r.get('benign_mean',0)):.3f}</td>"
                f"<td style='text-align:right'>{float(r.get('malignant_mean',0)):.3f}</td>"
                f"<td style='text-align:right'><strong>{float(r.get('sample_delta',0)):.3f}</strong></td>"
                f"<td style='text-align:right'>{float(r.get('wilcox_pval',1)):.4f}</td>"
                f"<td style='text-align:center'>{'✓' if r.get('pass_p') else '—'}</td>"
                f"<td style='text-align:center'>{'✓' if r.get('pass_delta') else '—'}</td>"
                f"<td style='text-align:center'>{'<strong style=color:#22863a>Pass</strong>' if r.get('pass_both') else '<span style=color:#aaa>Fail</span>'}</td>"
                f"</tr>"
                for idx, r in enumerate(wilcox_records)   # already sorted by delta desc
            )
            wcontrols, wjs = _pagination_widget(WTBL, WCTRL, n_candidates)

            wilcox_section = f"""
<h3 style="margin-top:20px">Sample-level Wilcoxon Pre-filter</h3>
<p style="color:#666;font-size:.83em">
  The diff_methylation.txt is derived from <em>merged group bedGraphs</em>, which average out
  within-group noise and can inflate apparent |Δβ|. This step runs a Mann-Whitney U test on
  individual sample methylation values to retain only sites with genuine per-sample separation.
</p>
<table>
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Sites entering filter (after coverage filter)</td><td>{n_candidates}</td></tr>
  <tr><td>Pass both (p &lt; threshold <em>and</em> |Δβ| ≥ threshold)</td>
      <td><strong style="color:#22863a">{n_wpass}</strong></td></tr>
  <tr><td>Fail — p not significant, |Δβ| OK</td><td>{n_fail_p}</td></tr>
  <tr><td>Fail — p significant, |Δβ| too small</td><td>{n_fail_d}</td></tr>
  <tr><td>Fail — both criteria missed</td><td>{n_fail_b}</td></tr>
</table>
<table style="margin-top:10px">
  <tr><th>Sample-level |Δβ| stats</th><th>Value</th></tr>
  <tr><td>Mean</td><td>{delta_mean:.3f}</td></tr>
  <tr><td>Median</td><td>{delta_median:.3f}</td></tr>
  <tr><td>Max</td><td>{delta_max:.3f}</td></tr>
</table>
<table style="margin-top:10px">
  <tr><th>Wilcoxon p-value distribution (of {n_candidates} sites)</th><th>Count</th></tr>
  <tr><td>p &lt; 0.05</td><td>{n_p05}</td></tr>
  <tr><td>p &lt; 0.10</td><td>{n_p10}</td></tr>
  <tr><td>p &lt; 0.20</td><td>{n_p20}</td></tr>
</table>
<details style="margin-top:10px">
  <summary style="cursor:pointer;color:#1a6fbb;font-size:.88em">
    Show all {n_candidates} candidate sites (sorted by sample |Δβ|, green = passed)
  </summary>
  <div style="overflow-x:auto;margin-top:8px">
  {wcontrols}
  <table id="{WTBL}" style="font-size:.82em">
    <thead><tr><th>Site</th><th>Benign mean</th><th>Malignant mean</th>
        <th>|Δβ| (sample)</th><th>p-value</th>
        <th>p OK</th><th>|Δβ| OK</th><th>Result</th></tr></thead>
    <tbody>{all_rows}</tbody>
  </table>
  </div>
</details>
<p style="color:#666;font-size:.83em;margin-top:6px">
  Full stats: <code>wilcoxon_stats.csv</code> &nbsp;|&nbsp;
  Thresholds configurable in <code>config.toml</code> [ml] sample_wilcox_p / sample_delta_beta
</p>
{wjs}"""
            wilcox_data = {
                "n_candidates": n_candidates, "n_pass": n_wpass,
                "n_fail_p": n_fail_p, "n_fail_delta": n_fail_d, "n_fail_both": n_fail_b,
                "delta_mean": round(delta_mean, 3), "delta_max": round(delta_max, 3),
            }
        else:
            wilcox_section = ""
            wilcox_data    = None

        # ── depth stats for panel sites ──
        depth_records = _read_csv_safe(os.path.join(result_dir, "panel_cpg_depth.csv"))
        if depth_records:
            depths       = [r.get("mean_depth", 0) for r in depth_records]
            medians      = [r.get("median_depth", 0) for r in depth_records]
            cov_rates    = [r.get("coverage_rate", 0) for r in depth_records]
            overall_mean = sum(depths) / len(depths) if depths else 0
            overall_med  = sorted(medians)[len(medians) // 2] if medians else 0
            mean_cov     = sum(cov_rates) / len(cov_rates) if cov_rates else 0
            n_lt5   = sum(1 for d in depths if d < 5)
            n_lt10  = sum(1 for d in depths if d < 10)
            n_ge20  = sum(1 for d in depths if d >= 20)

            depth_warn = ""
            if n_lt10 > 0:
                pct = n_lt10 / n_panel * 100
                depth_warn = (f"<p style='color:#c0392b;font-size:.83em'>"
                              f"&#9888; {n_lt10} sites ({pct:.0f}%) have mean depth &lt; 10x "
                              f"— low-depth sites may reduce classifier reliability.</p>")

            # top-10 lowest depth sites
            top_low = sorted(depth_records, key=lambda r: r.get("mean_depth", 0))[:10]
            low_rows = "".join(
                f"<tr><td>{r['feature']}</td>"
                f"<td>{r.get('mean_depth',0):.1f}</td>"
                f"<td>{r.get('median_depth',0):.1f}</td>"
                f"<td>{r.get('coverage_rate',0)*100:.0f}%</td>"
                f"<td>{r.get('min_depth',0)}</td>"
                f"<td>{r.get('max_depth',0)}</td></tr>"
                for r in top_low
            )
            depth_section = f"""
<h3 style="margin-top:18px">Panel CpG Site Depth (classifier features)</h3>
<table>
  <tr><th>Metric</th><th>Value</th></tr>
  <tr><td>Panel sites</td><td>{n_panel}</td></tr>
  <tr><td>Mean depth across all sites &amp; samples</td><td>{overall_mean:.1f}×</td></tr>
  <tr><td>Median depth across all sites &amp; samples</td><td>{overall_med:.1f}×</td></tr>
  <tr><td>Mean sample coverage rate per site</td><td>{mean_cov*100:.1f}%</td></tr>
  <tr><td>Sites with mean depth &lt; 5×</td><td>{n_lt5} ({n_lt5/n_panel*100:.0f}%)</td></tr>
  <tr><td>Sites with mean depth &lt; 10×</td><td>{n_lt10} ({n_lt10/n_panel*100:.0f}%)</td></tr>
  <tr><td>Sites with mean depth ≥ 20×</td><td>{n_ge20} ({n_ge20/n_panel*100:.0f}%)</td></tr>
</table>
{depth_warn}
<details style="margin-top:10px">
  <summary style="cursor:pointer;color:#1a6fbb;font-size:.88em">Show 10 lowest-depth panel sites</summary>
  <table style="margin-top:8px">
    <tr><th>Site</th><th>Mean depth</th><th>Median depth</th><th>Coverage rate</th><th>Min</th><th>Max</th></tr>
    {low_rows}
  </table>
</details>
<p style="color:#666;font-size:.83em;margin-top:6px">Depth = meth_reads + unmeth_reads per site per sample (MethylDackel extract);
full per-site stats: <code>panel_cpg_depth.csv</code></p>"""
            depth_data = {"mean_depth": round(overall_mean, 1),
                          "median_depth": round(overall_med, 1),
                          "mean_coverage_rate": round(mean_cov, 4),
                          "n_lt5x": n_lt5, "n_lt10x": n_lt10, "n_ge20x": n_ge20}
        else:
            depth_section = ""
            depth_data = None

        funnel_src = _copy_img(os.path.join(result_dir, "feature_selection_funnel.png"), imgs_dir)
        funnel_html = (
            f"<div style='margin:12px 0'><h3 style='font-size:.95em;margin-bottom:4px'>"
            f"Figure 2A — Feature Selection Funnel</h3>"
            f"{_img_tag(funnel_src, 'Feature selection funnel')}</div>"
            if funnel_src else ""
        )
        sec9_html = f"""
<div class="section">
<h2>Step 14. Feature Selection — LASSO + RF</h2>
<p style="color:#666;font-size:.85em">
  Pipeline: group-level DMC (merged bedGraph) →
  coverage filter (≥50% samples) →
  <strong>sample-level Wilcoxon filter</strong> →
  LASSO ∪ RF → final panel
</p>
<table>
  <tr><th>Stage</th><th>Sites</th></tr>
  <tr><td>Candidate sites entering Wilcoxon filter (after ≥50% coverage filter)</td>
      <td>{n_candidates if wilcox_records else n_total}</td></tr>
  <tr><td>Passed sample-level Wilcoxon (p &lt; threshold <em>and</em> |Δβ| ≥ threshold)</td>
      <td>{n_wpass if wilcox_records else "—"}</td></tr>
  <tr><td>LASSO selected (from Wilcoxon-passed)</td><td>{lasso_sel}</td></tr>
  <tr><td>RF selected (from Wilcoxon-passed)</td><td>{rf_sel}</td></tr>
  <tr><td>Final panel (union, for classifier)</td><td><strong>{n_panel}</strong></td></tr>
</table>
<p style="color:#666;font-size:.83em">Method: LASSO (LogisticRegressionCV, L1) ∪ Random Forest (Gini importance &gt; 0);
outputs: <code>panel_cpg.csv</code>, <code>feature_importance.csv</code></p>
{funnel_html}
{wilcox_section}
{depth_section}
</div>"""
        sec9_data = {"n_total": n_total, "lasso_selected": lasso_sel,
                     "rf_selected": rf_sel, "panel_size": n_panel,
                     "wilcoxon": wilcox_data, "depth": depth_data}
    else:
        sec9_html = ('<div class="section"><h2>Step 14. Feature Selection — LASSO + RF</h2>'
                     '<p style="color:#aaa;font-size:.9em">Step 14 not yet run — '
                     'execute <code>--step feature_select</code> to populate this section.</p></div>')
        sec9_data = None

    # ── Section 10: Diagnostic Classifiers (Step 15) ──
    clf_sum_path  = os.path.join(result_dir, "classifier_summary.csv")
    roc_src       = _copy_img(os.path.join(result_dir, "roc_curve.png"), imgs_dir)
    cm_src        = _copy_img(os.path.join(result_dir, "confusion_matrix.png"), imgs_dir)
    clf_records   = _read_csv_safe(clf_sum_path)

    if clf_records is not None:
        clf_rows = ""
        for r in clf_records:
            clf_rows += (
                f"<tr><td>{r.get('classifier','')}</td>"
                f"<td>{r.get('auc', 0):.3f}</td>"
                f"<td>{r.get('accuracy', 0):.3f}</td>"
                f"<td>{r.get('sensitivity', 0):.3f}</td>"
                f"<td>{r.get('specificity', 0):.3f}</td>"
                f"<td>{r.get('n_features', '')}</td></tr>\n"
            )
        roc_note   = _img_tag(roc_src,  "ROC Curve")
        cm_note    = _img_tag(cm_src,   "Confusion Matrix")

        # Feature selection ablation table (Step 15b — optional)
        bench_records = _read_csv_safe(os.path.join(result_dir, "feature_selection_benchmark.csv"))
        if bench_records:
            _MODES  = ["wilcoxon", "lasso", "rf", "union"]
            _CLFS   = ["LR", "SVM", "RF", "XGB", "LGBM"]
            _MODE_LABEL = {
                "wilcoxon": "Wilcoxon-only",
                "lasso":    "LASSO-only",
                "rf":       "RF-only",
                "union":    "LASSO ∪ RF (AIMethFlow)",
            }
            # Build lookup: (mode, clf) → row
            bench_map = {(r.get("mode",""), r.get("classifier","")): r for r in bench_records}
            # Header: Classifier | Wilcoxon | LASSO | RF | LASSO∪RF
            bench_header = "<tr><th>Classifier</th>" + "".join(
                f"<th>{_MODE_LABEL.get(m, m)}<br/><small>n_feat</small></th>"
                for m in _MODES
            ) + "</tr>\n"
            bench_body = ""
            for clf in _CLFS:
                bench_body += f"<tr><td><strong>{clf}</strong></td>"
                best_auc = max(
                    float(bench_map.get((m, clf), {}).get("auc", 0)) for m in _MODES
                )
                for m in _MODES:
                    r = bench_map.get((m, clf), {})
                    auc_v  = float(r.get("auc", 0))
                    nf     = r.get("n_features_mean", "—")
                    bold   = " font-weight:700;color:#1a6fbb" if abs(auc_v - best_auc) < 1e-6 else ""
                    bench_body += (
                        f"<td style='text-align:center{bold}'>"
                        f"{auc_v:.3f}<br/><small style='color:#888'>{nf}</small></td>"
                    )
                bench_body += "</tr>\n"
            bench_html = f"""
<h3 style="margin-top:20px">Step 15b. Feature Selection Ablation (wilcoxon / lasso / rf / union)</h3>
<p style="color:#666;font-size:.83em">
  Same CV strategy, same classifiers. AUC shown; bold = best per row.
  <em>n_feat</em> = mean features selected per fold.
  Full table: <code>feature_selection_benchmark.csv</code>
</p>
<table style="font-size:.88em">
  {bench_header}
  {bench_body}
</table>"""
        else:
            bench_html = (
                "<p style='color:#aaa;font-size:.85em;margin-top:10px'>"
                "Feature selection ablation not yet run — "
                "execute <code>--step fs_benchmark</code> to populate.</p>"
            )

        # ── CA-125 comparison block (Step 15c — optional) ──────────────────────
        ca125_records = _read_csv_safe(os.path.join(result_dir, "ca125_comparison.csv"))
        if ca125_records:
            ca125_rows = ""
            for r in ca125_records:
                is_ca125 = r.get("method") == "CA-125"
                style    = " style='background:#eaf4fb'" if is_ca125 else ""
                ca125_rows += (
                    f"<tr{style}>"
                    f"<td><strong>{r.get('method','')}</strong></td>"
                    f"<td style='text-align:right'>{r.get('n_ca125_subset','')}</td>"
                    f"<td style='text-align:right'><strong>{float(r.get('auc',0)):.3f}</strong></td>"
                    f"<td style='text-align:right'>{float(r.get('auc_ci95_lo',0)):.3f}"
                    f" – {float(r.get('auc_ci95_hi',0)):.3f}</td>"
                    f"</tr>\n"
                )
            n_sub = ca125_records[0].get("n_ca125_subset", "?")
            n_tot = ca125_records[0].get("n_total", "?")
            ca125_html = f"""
<h3 style="margin-top:20px">Step 15c. CA-125 vs Classifier AUC
  <span style="font-size:.82em;font-weight:normal;color:#666">
    (n = {n_sub} with CA-125 data / {n_tot} total; 95% bootstrap CI)
  </span>
</h3>
<table style="font-size:.88em">
  <tr><th>Method</th><th>n subset</th><th>AUC</th><th>95% CI</th></tr>
  {ca125_rows}
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Classifiers re-evaluated on the CA-125-available subset for direct comparison.
  CA-125 binary threshold = median value in subset.
  Full data: <code>ca125_comparison.csv</code>, <code>ca125_values.csv</code>
</p>"""
        else:
            ca125_html = (
                "<p style='color:#aaa;font-size:.85em;margin-top:10px'>"
                "CA-125 comparison not available — add a <code>ca125</code> column to the "
                "sample info Excel, then run <code>--step ca125_compare</code>.</p>"
            )

        sec10_html = f"""
<div class="section">
<h2>Step 15. Diagnostic Classifiers — LOO-CV, Benign vs Malignant</h2>
<table>
  <tr><th>Classifier</th><th>AUC</th><th>Accuracy</th><th>Sensitivity</th><th>Specificity</th><th>Features</th></tr>
  {clf_rows}
</table>
<div class="two-col" style="margin-top:14px">
  <div><h3>ROC Curves</h3>{roc_note}</div>
  <div><h3>Confusion Matrices</h3>{cm_note}</div>
</div>
<p style="color:#666;font-size:.83em">Classifiers: LR / SVM / RF / XGB / LGBM evaluated via Leave-One-Out CV (LOO-CV);
StandardScaler fitted within each training fold to prevent data leakage;
outputs: <code>classifier_summary.csv</code>, <code>classifier_predictions.csv</code></p>
{bench_html}
{ca125_html}
</div>"""
        sec10_data = clf_records
    else:
        sec10_html = ('<div class="section"><h2>Step 15. Diagnostic Classifiers — LOO-CV, Benign vs Malignant</h2>'
                      '<p style="color:#aaa;font-size:.9em">Step 15 not yet run — '
                      'execute <code>--step classify</code> to populate this section.</p></div>')
        sec10_data = None

    # ── Section 11: Subtype Classifier (Step 16) ──
    sub_sum_path   = os.path.join(result_dir, "subtype_summary.csv")
    sub_roc_src    = _copy_img(os.path.join(result_dir, "subtype_roc.png"), imgs_dir)
    rf_cm_src      = _copy_img(os.path.join(result_dir, "subtype_confusion_RF.png"), imgs_dir)
    svm_cm_src     = _copy_img(os.path.join(result_dir, "subtype_confusion_SVM.png"), imgs_dir)
    hgsoc_sum_path = os.path.join(result_dir, "ovr_summary.csv")
    hgsoc_roc_src  = _copy_img(os.path.join(result_dir, "ovr_roc.png"), imgs_dir)
    sub_records    = _read_csv_safe(sub_sum_path)
    hgsoc_records  = _read_csv_safe(hgsoc_sum_path)

    # Check for a skip marker written when data was insufficient
    skip_path = os.path.join(result_dir, "subtype_skip.json")
    skip_info = None
    if os.path.exists(skip_path):
        try:
            import json as _json
            with open(skip_path, encoding="utf-8") as _f:
                skip_info = _json.load(_f)
        except Exception:
            pass

    if skip_info and skip_info.get("skipped"):
        skip_reason = skip_info.get("reason", "")
        skip_warns  = skip_info.get("warnings", [])
        warn_html   = "".join(f"<li>{w}</li>" for w in skip_warns) if skip_warns else ""
        sec11_html  = f"""
<div class="section">
<h2>Step 16. Subtype Classifier — Skipped</h2>
<p style="color:#c0392b;font-size:.9em">
  &#9888; Subtype classification was skipped because the dataset does not meet
  the minimum requirements for reliable LOO-CV training.
</p>
<table>
  <tr><th>Reason</th><td>{skip_reason}</td></tr>
</table>
{"<ul style='color:#856404;font-size:.85em;margin-top:8px'>" + warn_html + "</ul>" if warn_html else ""}
<p style="color:#666;font-size:.83em;margin-top:8px">
  Collect additional malignant/borderline samples with subtype labels and re-run
  <code>--step subtype_classify</code> to enable this section.
</p>
</div>"""
        sec11_data = None

    elif sub_records is not None:
        sub_rows = ""
        for r in sub_records:
            # Per-row data-quality warning indicator
            row_warn = r.get("warnings", "")
            warn_cell = (f"<td style='color:#856404;font-size:.82em'>{row_warn}</td>"
                         if row_warn else "<td></td>")
            sub_rows += (
                f"<tr><td>{r.get('classifier','')}</td>"
                f"<td>{float(r.get('accuracy', 0)):.3f}</td>"
                f"<td>{r.get('n_samples', '')}</td>"
                f"<td>{r.get('n_subtypes', '')}</td>"
                f"<td style='font-size:.82em;color:#555'>{r.get('class_counts','')}</td>"
                f"{warn_cell}</tr>\n"
            )
        # Collect any warnings to show above the table
        all_warns = list({r.get("warnings", "") for r in sub_records if r.get("warnings")})
        warn_banner = ""
        if all_warns:
            items = "".join(f"<li>{w}</li>" for w in all_warns)
            warn_banner = (
                f"<div style='background:#fff3cd;border:1px solid #ffc107;"
                f"border-radius:4px;padding:8px 12px;margin-bottom:10px;"
                f"font-size:.85em;color:#856404'>"
                f"&#9888; Data quality notes:<ul style='margin:4px 0 0 16px'>{items}</ul></div>"
            )
        sec11_html = f"""
<div class="section">
<h2>Step 16. Subtype Classifier — LOO-CV, Malignant + Borderline, Multi-class</h2>
<p style="color:#666;font-size:.85em">
  Subtypes: HGSC / Endometrioid / Clear cell / LGSC / Mucinous.
  Strategy: One-vs-Rest (OvR); evaluation: Leave-One-Out CV.
  Fold-wise imputation prevents data leakage.
</p>
{warn_banner}
<table>
  <tr><th>Classifier</th><th>Accuracy</th><th>Samples</th><th>Subtypes</th>
      <th>Class counts</th><th>Data warnings</th></tr>
  {sub_rows}
</table>
<div class="two-col" style="margin-top:14px">
  <div><h3>RF Confusion Matrix</h3>{_img_tag(rf_cm_src, 'RF Confusion Matrix')}</div>
  <div><h3>SVM Confusion Matrix</h3>{_img_tag(svm_cm_src, 'SVM Confusion Matrix')}</div>
</div>
<div style="margin-top:14px"><h3>OvR ROC Curves</h3>{_img_tag(sub_roc_src, 'Subtype ROC')}</div>
<p style="color:#666;font-size:.83em">
  Outputs: <code>subtype_summary.csv</code>, <code>subtype_predictions.csv</code>
  (includes per-sample confidence = max predicted probability)
</p>
{_build_ovr_subsection(hgsoc_records, hgsoc_roc_src)}
</div>"""
        sec11_data = sub_records

    else:
        sec11_html = ('<div class="section"><h2>Step 16. Subtype Classifier — LOO-CV, Malignant + Borderline</h2>'
                      '<p style="color:#aaa;font-size:.9em">Step 16 not yet run — '
                      'execute <code>--step subtype_classify</code> to populate this section.</p></div>')
        sec11_data = None

    clf_data = {
        "feature_selection":  sec9_data,
        "outlier_detection":  None,          # filled below
        "classifier":         sec10_data,
        "subtype_classifier": sec11_data,
    }

    # ── Outlier detection section (Step 13) ──
    sec_outlier_html, outlier_data = _build_outlier_section(result_dir, imgs_dir)
    clf_data["outlier_detection"] = outlier_data

    # ── Section 12: Temporal Split Validation (Step 17, plasma only) ──
    tv_sum_path  = os.path.join(result_dir, "temporal_validation_summary.csv")
    tv_pred_path = os.path.join(result_dir, "temporal_validation_predictions.csv")
    tv_records   = _read_csv_safe(tv_sum_path)

    if tv_records is not None:
        tv_rows = ""
        for r in tv_records:
            tv_rows += (
                f"<tr><td>{r.get('classifier','')}</td>"
                f"<td>{r.get('n_train','')}</td>"
                f"<td>{r.get('n_validation','')}</td>"
                f"<td>{r.get('auc','')}</td>"
                f"<td>{r.get('accuracy','')}</td>"
                f"<td>{r.get('sensitivity','')}</td>"
                f"<td>{r.get('specificity','')}</td></tr>\n"
            )
        cutoff = tv_records[0].get("cutoff_date", "2024-07-01") if tv_records else "2024-07-01"
        sec12_html = f"""
<div class="section">
<h2>Step 17. Temporal Split Validation (plasma mode)</h2>
<p style="color:#666;font-size:.85em">
  Prospective simulation: samples collected before <strong>{cutoff}</strong> form the train set;
  samples collected on/after that date form the held-out validation set.
  Features: top-N sites by LOO-CV selection frequency (<code>feature_selection_freq.csv</code>),
  fitted on train only — no information from the validation set leaks in.
</p>
<table>
  <tr><th>Classifier</th><th>Train (n)</th><th>Validation (n)</th>
      <th>AUC</th><th>Accuracy</th><th>Sensitivity</th><th>Specificity</th></tr>
  {tv_rows}
</table>
<p style="color:#666;font-size:.83em">
  Outputs: <code>temporal_validation_summary.csv</code>,
  <code>temporal_validation_predictions.csv</code>
</p>
</div>"""
        clf_data["temporal_validation"] = tv_records
    else:
        sec12_html = ('<div class="section"><h2>Step 17. Temporal Split Validation (plasma mode)</h2>'
                      '<p style="color:#aaa;font-size:.9em">Step 17 not yet run — '
                      'execute <code>--step temporal_validate</code> (plasma mode only) to populate this section.</p></div>')
        clf_data["temporal_validation"] = None

    return sec9_html, sec_outlier_html, sec10_html, sec11_html, sec12_html, clf_data


# ─────────────────────────────────────────────
#  Main report generation function
# ─────────────────────────────────────────────

def generate_report(
    temp_dir,
    result_dir,
    sample1_methylation,
    sample2_methylation,
    candidate_sample1,
    candidate_sample2,
    benign_samples=None,
    malignant_samples=None,
):
    """Generate full analysis report (HTML + JSON) inside result_dir/report/."""
    logger.info("Generating analysis report")

    benign_samples    = benign_samples or []
    malignant_samples = malignant_samples or []

    # ── Report output directory (self-contained, easy to zip) ────────
    timestamp = os.path.basename(os.path.normpath(result_dir))
    report_dir = os.path.join(result_dir, f"{timestamp}_report")
    os.makedirs(report_dir, exist_ok=True)

    # ── Image output directory (inside report/) ──────────────────────
    imgs_dir = os.path.join(report_dir, "imgs")
    os.makedirs(imgs_dir, exist_ok=True)

    # ── Discover all samples via fastp JSON ──────────────────────
    all_fastp_files = sorted(glob.glob(os.path.join(temp_dir, "*.fastp.json")))
    all_samples     = [os.path.basename(f).replace(".fastp.json", "") for f in all_fastp_files]

    # ── Per-sample QC stats ──────────────────────────────────────
    fastp  = {s: _parse_fastp_json(os.path.join(temp_dir, f"{s}.fastp.json")) for s in all_samples}
    dedup  = {s: _parse_dedup_stats(os.path.join(temp_dir, f"{s}.bwameth_pe.dedup.stats")) for s in all_samples}

    def _bam(s, suffix):
        p = os.path.join(temp_dir, f"{s}{suffix}")
        return os.path.getsize(p) if os.path.exists(p) else 0

    bam_sorted  = {s: _bam(s, ".bwameth_pe.sorted.bam")                for s in all_samples}
    bam_dedup   = {s: _bam(s, ".bwameth_pe.dedup.bam")                 for s in all_samples}
    bam_target  = {s: _bam(s, ".bwameth_pe.dedup.targeted.sorted.bam") for s in all_samples}

    # ── Per-sample CpG site counts ───────────────────────────────
    cpg_per_sample = {}
    for s in all_samples:
        bg = os.path.join(temp_dir, f"{s}_CpG.bedGraph")
        cnt, avg = _count_bedgraph_sites(bg) if os.path.exists(bg) else (0, 0.0)
        cpg_per_sample[s] = {"count": cnt, "avg": avg}

    # ── Group-level CpG & candidate region stats ─────────────────
    n_cpg1,  avg_meth1  = _count_bedgraph_sites(sample1_methylation)  if sample1_methylation  else (0, 0.0)
    n_cpg2,  avg_meth2  = _count_bedgraph_sites(sample2_methylation)  if sample2_methylation  else (0, 0.0)
    n_cand1, avg_cand1  = _count_bedgraph_sites(candidate_sample1)    if candidate_sample1    else (0, 0.0)
    n_cand2, avg_cand2  = _count_bedgraph_sites(candidate_sample2)    if candidate_sample2    else (0, 0.0)

    def _full_meth_cnt(path):
        try:
            c = 0
            with open(path, encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("track"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        try:
                            if float(parts[3]) >= 100.0:
                                c += 1
                        except ValueError:
                            pass
            return c
        except Exception:
            return 0

    full1 = _full_meth_cnt(candidate_sample1) if candidate_sample1 else 0
    full2 = _full_meth_cnt(candidate_sample2) if candidate_sample2 else 0

    # ── Differential methylation results ─────────────────────────
    diff_file = os.path.join(result_dir, "diff_methylation.txt")
    diff_rows = _parse_diff_results(diff_file)
    n_dmc = len(diff_rows)
    chrom_dist = {}
    s1_higher = s2_higher = 0
    diffs = []
    for r in diff_rows:
        chrom_dist[r["chrom"]] = chrom_dist.get(r["chrom"], 0) + 1
        diffs.append(r["diff"])
        if r["meth1"] > r["meth2"]:
            s1_higher += 1
        else:
            s2_higher += 1
    mean_diff = sum(diffs) / len(diffs) if diffs else 0
    max_diff  = max(diffs) if diffs else 0

    # ── Visualization images ──────────────────────────────────────
    density_src = _copy_img(os.path.join(result_dir, "density_plot.png"), imgs_dir)
    scatter_src = _copy_img(os.path.join(result_dir, "scatter_plot.png"), imgs_dir)
    heatmap_src = _copy_img(os.path.join(result_dir, "heatmap_plot.png"), imgs_dir)

    # ── M-bias: copy SVGs to imgs/, build JS-paginated viewer ────
    mbias_samples = []
    for svg_path in sorted(glob.glob(os.path.join(temp_dir, "*.mbias_OT.svg"))):
        sname  = os.path.basename(svg_path).replace(".mbias_OT.svg", "")
        ot_src = _copy_img(svg_path, imgs_dir)
        ob_src = _copy_img(svg_path.replace("_OT.svg", "_OB.svg"), imgs_dir)
        if ot_src:
            mbias_samples.append({"name": sname, "ot": ot_src, "ob": ob_src})

    if mbias_samples:
        n_mbias = len(mbias_samples)
        slides_html = ""
        for i, s in enumerate(mbias_samples):
            flex_style = "display:flex;gap:10px;flex-wrap:wrap"
            hide_style = "display:none"
            slide_style = flex_style if i == 0 else hide_style
            ot_img = (f"<div style='flex:1;min-width:260px;text-align:center'>"
                      f"<p style='font-size:.8em;color:#555;margin:0 0 4px'>OT strand (Watson)</p>"
                      f"<img src='{s['ot']}' style='max-width:100%'/></div>")
            ob_img = ""
            if s["ob"]:
                ob_img = (f"<div style='flex:1;min-width:260px;text-align:center'>"
                          f"<p style='font-size:.8em;color:#555;margin:0 0 4px'>OB strand (Crick)</p>"
                          f"<img src='{s['ob']}' style='max-width:100%'/></div>")
            slides_html += (
                f"<div class='mbias-slide' style='{slide_style}'>"
                f"{ot_img}{ob_img}</div>"
            )
        names_json = json.dumps([s["name"] for s in mbias_samples])
        btn_style = ("cursor:pointer;padding:5px 14px;border:1px solid #1a6fbb;"
                     "border-radius:4px;background:#fff;color:#1a6fbb;font-size:.88em")
        mbias_body = f"""
<p style='color:#666;font-size:.83em;margin:0 0 10px'>
  Total: <strong>{n_mbias}</strong> sample(s) — OT (Watson) + OB (Crick) per sample.
</p>
<div style='display:flex;align-items:center;gap:14px;margin-bottom:12px;flex-wrap:wrap'>
  <button onclick='mbiasPrev()' style='{btn_style}'>&#9664; Previous</button>
  <span id='mbias-counter' style='font-weight:600;color:#1a6fbb;min-width:60px;text-align:center'>
    1 / {n_mbias}
  </span>
  <button onclick='mbiasNext()' style='{btn_style}'>Next &#9654;</button>
  <span id='mbias-name' style='color:#444;font-size:.9em;margin-left:6px'>
    {mbias_samples[0]["name"]}
  </span>
</div>
{slides_html}
<script>
(function(){{
  var idx=0, total={n_mbias};
  var names={names_json};
  function show(i){{
    var slides=document.querySelectorAll('.mbias-slide');
    slides.forEach(function(el,j){{el.style.display=j===i?'flex':'none';}});
    document.getElementById('mbias-counter').textContent=(i+1)+' / '+total;
    document.getElementById('mbias-name').textContent=names[i];
    idx=i;
  }}
  window.mbiasPrev=function(){{show((idx-1+total)%total);}};
  window.mbiasNext=function(){{show((idx+1)%total);}};
  document.addEventListener('keydown',function(e){{
    if(e.key==='ArrowLeft')window.mbiasPrev();
    if(e.key==='ArrowRight')window.mbiasNext();
  }});
}})();
</script>"""
    else:
        mbias_body = "<p style='color:#aaa'>M-bias SVG not yet generated — run --step mbias first</p>"

    # ── Fragment size distribution ────────────────────────────────
    frag_img_src = _copy_img(os.path.join(result_dir, "fragment_size_distribution.png"), imgs_dir)
    frag_stats_records = _read_csv_safe(os.path.join(result_dir, "fragment_size_stats.csv"))

    if frag_img_src:
        frag_img_html = (
            f"<img src='{frag_img_src}' style='max-width:100%;margin:10px 0'/>"
        )
    else:
        frag_img_html = ("<p style='color:#aaa;font-size:.9em'>"
                         "fragment_size_distribution.png not found — "
                         "run --step mbias to generate.</p>")

    if frag_stats_records:
        n_frag      = len(frag_stats_records)
        page_size   = 10
        n_pages     = max(1, (n_frag + page_size - 1) // page_size)
        first_end   = min(page_size, n_frag)
        rows_js     = json.dumps([
            {"sample": r.get("sample", ""),
             "n":      r.get("n_reads", ""),
             "med":    r.get("median", ""),
             "mean":   r.get("mean", ""),
             "pct":    r.get("pct_100_220", "")}
            for r in frag_stats_records
        ])
        btn_style = ("cursor:pointer;padding:5px 14px;border:1px solid #1a6fbb;"
                     "border-radius:4px;background:#fff;color:#1a6fbb;font-size:.88em")
        frag_stats_html = f"""
<div style='margin-top:10px'>
<div style='display:flex;align-items:center;gap:14px;margin-bottom:8px;flex-wrap:wrap'>
  <button onclick='fragPrev()' style='{btn_style}'>&#9664; Previous</button>
  <span id='frag-counter' style='font-weight:600;color:#1a6fbb;min-width:100px;text-align:center'>
    1–{first_end} / {n_frag}
  </span>
  <button onclick='fragNext()' style='{btn_style}'>Next &#9654;</button>
  <span style='color:#888;font-size:.83em;margin-left:4px'>({n_pages} page(s), {page_size} rows/page)</span>
</div>
<table id='frag-table' style='font-size:.83em'>
  <tr>
    <th>Sample</th><th>Reads sampled</th>
    <th>Median (bp)</th><th>Mean (bp)</th><th>% in 100–220 bp</th>
  </tr>
  <tbody id='frag-tbody'></tbody>
</table>
</div>
<script>
(function(){{
  var data={rows_js};
  var PAGE={page_size}, page=0, total=data.length;
  var pages=Math.ceil(total/PAGE);
  function render(){{
    var tbody=document.getElementById('frag-tbody');
    var start=page*PAGE, end=Math.min(start+PAGE,total);
    var html='';
    for(var i=start;i<end;i++){{
      var r=data[i];
      html+='<tr>'
        +'<td>'+r.sample+'</td>'
        +'<td style="text-align:right">'+r.n+'</td>'
        +'<td style="text-align:right">'+r.med+'</td>'
        +'<td style="text-align:right">'+r.mean+'</td>'
        +'<td style="text-align:right">'+r.pct+'%</td>'
        +'</tr>';
    }}
    tbody.innerHTML=html;
    document.getElementById('frag-counter').textContent=
      (start+1)+'–'+end+' / '+total;
  }}
  window.fragPrev=function(){{if(page>0){{page--;render();}}}};
  window.fragNext=function(){{if(page<pages-1){{page++;render();}}}};
  render();
}})();
</script>"""
    else:
        frag_stats_html = ""

    fragment_body = frag_img_html + frag_stats_html

    # ── Group labels ──────────────────────────────────────────────
    g1_label = "Benign"
    g2_label = "Malignant"
    g1_samples_str = ", ".join(benign_samples)    if benign_samples    else "—"
    g2_samples_str = ", ".join(malignant_samples) if malignant_samples else "—"

    # ── Results interpretation ────────────────────────────────────
    interp_items = []

    n_benign    = len(benign_samples)
    n_malignant = len(malignant_samples)
    interp_items.append(
        f"<li>Analysis includes <strong>{len(all_samples)}</strong> samples: "
        f"<strong>{n_benign}</strong> benign and <strong>{n_malignant}</strong> malignant.</li>"
    )

    q30_vals = [fastp[s].get("q30_after", 0) for s in all_samples if fastp[s]]
    avg_q30  = sum(q30_vals) / len(q30_vals) if q30_vals else 0
    q30_eval = ("excellent (≥85%)" if avg_q30 >= 0.85
                else "acceptable (70–85%)" if avg_q30 >= 0.70
                else "low (<70%) — consider re-checking sequencing quality")
    interp_items.append(
        f"<li>Mean Q30 rate: <strong>{avg_q30*100:.1f}%</strong> — data quality is {q30_eval}.</li>"
    )

    dup_vals = [dedup[s].get("dup_rate", None) for s in all_samples if dedup[s].get("dup_rate") is not None]
    if dup_vals:
        avg_dup = sum(dup_vals) / len(dup_vals)
        dup_eval = ("high (>60%) — check PCR cycle count during library preparation" if avg_dup > 0.60
                    else "normal (30–60%)" if avg_dup >= 0.30
                    else "low (<30%) — high data utilization rate")
        interp_items.append(
            f"<li>Mean UMI deduplication rate: <strong>{avg_dup*100:.1f}%</strong> — {dup_eval}.</li>"
        )

    ont_vals = []
    for s in all_samples:
        if bam_dedup[s] and bam_target[s]:
            ont_vals.append(bam_target[s] / bam_dedup[s])
    if ont_vals:
        avg_ont = sum(ont_vals) / len(ont_vals)
        ont_eval = ("good (>30%)" if avg_ont > 0.30
                    else "low (<30%) — possible off-target reads or poor capture efficiency")
        interp_items.append(
            f"<li>Mean on-target rate (file-size estimate): <strong>{avg_ont*100:.1f}%</strong> — {ont_eval}.</li>"
        )

    if n_dmc == 0:
        interp_items.append(
            "<li>No differentially methylated CpGs found at |Δβ| ≥ 0.2. "
            "The two groups may have similar methylation profiles, or sample size may be insufficient.</li>"
        )
    elif n_dmc < 10:
        interp_items.append(
            f"<li>Found <strong>{n_dmc}</strong> DMC sites — relatively few; interpret in clinical context.</li>"
        )
    elif n_dmc < 50:
        interp_items.append(
            f"<li>Found <strong>{n_dmc}</strong> DMC sites — moderate methylation differences between groups.</li>"
        )
    else:
        interp_items.append(
            f"<li>Found <strong>{n_dmc}</strong> DMC sites — strong methylation differences, "
            f"good discriminatory potential.</li>"
        )

    if n_dmc > 0:
        direction = g2_label if s2_higher > s1_higher else g1_label
        interp_items.append(
            f"<li><strong>{direction}</strong> shows higher overall methylation "
            f"({max(s1_higher, s2_higher)}/{n_dmc} sites); "
            f"mean |Δβ| = {mean_diff:.3f}, max |Δβ| = {max_diff:.3f}.</li>"
        )
        if max_diff >= 0.5:
            interp_items.append(
                "<li>Sites with |Δβ| ≥ 0.5 are present — prioritise these as candidate biomarkers for validation.</li>"
            )

    if n_cand1 > 0 or n_cand2 > 0:
        cov_note = ("sufficient (recommend proceeding to model training)"
                    if min(n_cand1, n_cand2) > 100
                    else "limited (consider lowering depth threshold)")
        interp_items.append(
            f"<li>Within the 88 candidate regions: Benign group has <strong>{_fmt_num(n_cand1)}</strong> CpG sites, "
            f"Malignant group has <strong>{_fmt_num(n_cand2)}</strong>; coverage is {cov_note}.</li>"
        )

    interp_html = "<ul style='line-height:1.9'>" + "".join(interp_items) + "</ul>"

    # ── Build HTML ────────────────────────────────────────────────

    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Optional ML sections (Steps 13–17; empty string if files absent)
    sec8_html, sec9_html, sec10_html, sec11_html, sec12_html, clf_data = _build_classifier_section(result_dir, imgs_dir)
    # Section 13: annotation (Step 18; silent if panel_annotation.csv absent)
    sec13_html, annot_data = _build_annotation_section(result_dir)
    clf_data["annotation"] = annot_data
    # Section 14: public-dataset cross-validation (Step 19; silent if external_validation.csv absent)
    sec14_html, ext_val_data = _build_external_validation_section(result_dir)
    clf_data["external_validation"] = ext_val_data

    # iMeta publication figures (auto-generated; silently skipped if input CSVs absent)
    try:
        from scripts.modules.imeta_utils import build_imeta_section as _build_imeta
        sec_imeta_html = _build_imeta(result_dir, imgs_dir)
    except Exception as _imeta_exc:
        logger.warning("iMeta figure generation skipped: %s", _imeta_exc)
        sec_imeta_html = ""

    # QC table rows
    qc_rows = ""
    for s in all_samples:
        fp = fastp.get(s, {})
        pairs_in   = fp.get("total_reads_before", 0) // 2
        pairs_out  = fp.get("total_reads_after",  0) // 2
        q30        = fp.get("q30_after", 0)
        pass_rate  = pairs_out / pairs_in if pairs_in else 0
        grp_badge  = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                      f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                      "<span class='badge badge-grey'>Ungrouped</span>")
        qc_rows += (
            f"<tr><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(pairs_in)}</td><td>{_fmt_num(pairs_out)}</td>"
            f"<td>{_fmt_pct(pass_rate)}</td><td>{_fmt_pct(q30)}</td></tr>\n"
        )

    # Alignment / dedup / on-target table rows
    aln_rows = ""
    for s in all_samples:
        dp = dedup.get(s, {})
        inp  = dp.get("input_pairs", 0)
        out  = dp.get("output_reads", 0) // 2
        dr   = dp.get("dup_rate", None)
        ont_r = (bam_target[s] / bam_dedup[s]) if bam_dedup[s] else None
        aln_rows += (
            f"<tr><td>{s}</td>"
            f"<td>{_fmt_gb(bam_sorted[s])}</td>"
            f"<td>{_fmt_num(inp)}</td>"
            f"<td>{_fmt_num(out)}</td>"
            f"<td>{_fmt_pct(dr)}</td>"
            f"<td>{_fmt_mb(bam_target[s])}</td>"
            f"<td>{_fmt_pct(ont_r)}</td></tr>\n"
        )

    # Per-sample CpG table rows
    cpg_rows = ""
    for s in all_samples:
        c = cpg_per_sample.get(s, {})
        grp_badge = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                     f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                     "<span class='badge badge-grey'>Ungrouped</span>")
        cpg_rows += (
            f"<tr><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(c.get('count', 0))}</td>"
            f"<td>{c.get('avg', 0):.1f}%</td></tr>\n"
        )

    # Chromosome distribution table
    chrom_rows = "".join(
        f"<tr><td>{c}</td><td>{n}</td><td>{n/n_dmc*100:.1f}%</td></tr>"
        for c, n in sorted(chrom_dist.items())
    ) if n_dmc > 0 else "<tr><td colspan='3'>No DMC sites</td></tr>"

    # Top-20 DMC sites
    top_diff = sorted(diff_rows, key=lambda x: -x["diff"])[:20]
    top_diff_rows = "".join(
        f"<tr><td>{r['chrom']}</td><td>{r['pos']:,}</td>"
        f"<td>{r['meth1']:.1f}%</td><td>{r['meth2']:.1f}%</td>"
        f"<td class='{'up' if r['meth1']>r['meth2'] else 'down'}'>{r['diff']:.3f}</td></tr>"
        for r in top_diff
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8"/>
<title>Plasma Methylation Analysis Report</title>
<style>
  * {{ box-sizing: border-box; }}
  body {{ font-family: Arial, sans-serif; margin: 0; background: #f4f6f9; color: #333; }}
  .wrap {{ max-width: 1200px; margin: 0 auto; padding: 20px; }}
  h1 {{ background: linear-gradient(135deg,#1a6fbb,#0d4f9c); color: #fff; padding: 22px 30px; margin: 0 0 24px; border-radius: 8px; font-size: 1.5em; }}
  h1 span {{ font-size: .75em; font-weight: normal; opacity: .85; }}
  h2 {{ color: #1a6fbb; border-left: 4px solid #1a6fbb; padding-left: 10px; margin-top: 30px; font-size: 1.1em; }}
  h3 {{ font-size: .95em; color: #444; margin-top: 16px; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fill,minmax(180px,1fr)); gap: 14px; margin: 14px 0; }}
  .card {{ background: #fff; border-radius: 8px; padding: 16px 18px; box-shadow: 0 1px 4px rgba(0,0,0,.08); }}
  .card .num {{ font-size: 1.7em; font-weight: bold; color: #1a6fbb; }}
  .card .lbl {{ font-size: .82em; color: #666; margin-top: 4px; }}
  table {{ width: 100%; border-collapse: collapse; background: #fff; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 4px rgba(0,0,0,.08); margin: 10px 0; font-size: .88em; }}
  th {{ background: #1a6fbb; color: #fff; padding: 9px 12px; text-align: left; white-space: nowrap; }}
  td {{ padding: 8px 12px; border-bottom: 1px solid #eee; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f0f6ff; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: .78em; font-weight: bold; }}
  .badge-blue  {{ background: #dbeafe; color: #1d4ed8; }}
  .badge-red   {{ background: #fee2e2; color: #b91c1c; }}
  .badge-grey  {{ background: #f3f4f6; color: #6b7280; }}
  .up   {{ color: #c0392b; font-weight: bold; }}
  .down {{ color: #2980b9; font-weight: bold; }}
  .section {{ background: #fff; border-radius: 8px; padding: 20px 24px; box-shadow: 0 1px 4px rgba(0,0,0,.08); margin-bottom: 18px; }}
  .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }}
  .group-box {{ border-radius: 8px; padding: 14px 18px; }}
  .group-benign    {{ background: #eff6ff; border: 1px solid #bfdbfe; }}
  .group-malignant {{ background: #fef2f2; border: 1px solid #fecaca; }}
  .funnel {{ display: flex; align-items: center; gap: 0; margin: 14px 0; flex-wrap: wrap; }}
  .funnel-step {{ background: #1a6fbb; color: #fff; padding: 9px 16px; font-size: .82em; text-align: center; min-width: 100px; }}
  .funnel-step:nth-child(odd) {{ background: #2980b9; }}
  .funnel-arrow {{ color: #1a6fbb; font-size: 1.3em; padding: 0 3px; }}
  .interp {{ background: #f0fdf4; border-left: 4px solid #22c55e; padding: 14px 18px; border-radius: 0 8px 8px 0; margin: 10px 0; }}
  .meta {{ color: #888; font-size: .82em; text-align: right; margin-top: 8px; }}
  @media(max-width:700px){{ .two-col {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<div class="wrap">

<h1>Plasma cfDNA Targeted Methylation Analysis Report
  <br/><span>Analysis date: {analysis_date} &nbsp;|&nbsp;
  {len(all_samples)} samples ({g1_label}: {len(benign_samples)} / {g2_label}: {len(malignant_samples)})</span>
</h1>

<!-- Key metric cards -->
<div class="cards">
  <div class="card"><div class="num">{len(all_samples)}</div><div class="lbl">Total Samples</div></div>
  <div class="card"><div class="num">{len(benign_samples)}</div><div class="lbl">Benign</div></div>
  <div class="card"><div class="num">{len(malignant_samples)}</div><div class="lbl">Malignant</div></div>
  <div class="card"><div class="num">{_fmt_num(n_cand1 + n_cand2)}</div><div class="lbl">Candidate Region CpGs</div></div>
  <div class="card"><div class="num">{n_dmc}</div><div class="lbl">DMC Sites</div></div>
</div>

<!-- 1. Data Quality Overview -->
<div class="section">
<h2>Step 4. Data Quality Overview — fastp QC</h2>
<table>
  <tr><th>Sample</th><th>Group</th><th>Raw Read Pairs</th><th>Filtered Pairs</th><th>Pass Rate</th><th>Q30 (post-QC)</th></tr>
  {qc_rows}
</table>
<p style="color:#666;font-size:.83em">Parameters: UMI extraction (R2 5' 6 N) + trim_front1=12 + trim_front2=12 + trim_tail2=5 + min_length=35 + Q20</p>
</div>

<!-- 2. Alignment & Deduplication -->
<div class="section">
<h2>Steps 5–7. Alignment, Deduplication &amp; On-target Filter</h2>
<div style="overflow-x:auto">
<table>
  <tr><th>Sample</th><th>Sorted BAM</th><th>Input Pairs (pre-dedup)</th><th>Dedup Pairs</th><th>Dup Rate</th><th>On-target BAM</th><th>On-target Rate</th></tr>
  {aln_rows}
</table>
</div>
<p style="color:#666;font-size:.83em">Dedup: umi_tools dedup --paired; On-target: bedtools intersect panel_B10308.bed; BAM size of 0 means the file was auto-cleaned (expected).</p>
</div>

<!-- 3. M-bias Check -->
<div class="section">
<h2>Step 8. M-bias Check — MethylDackel mbias</h2>
<p style="color:#666;font-size:.85em">OT = Watson-strand methylation (~99% of reads). If the R2 3' curve drops, increase trim_tail2 and re-run from --step extract.</p>
{mbias_body}
<h3 style="margin-top:20px">cfDNA Fragment Size Distribution</h3>
<p style="color:#666;font-size:.85em">
  Insert sizes from paired-end BAM (TLEN, samtools view -f 0x2 -F 0x904).
  Expected mono-nucleosomal cfDNA peak: 147–167 bp.
  Shaded zone 100–220 bp; dashed line at 167 bp.
</p>
{fragment_body}
</div>

<!-- 4. Sample Grouping -->
<div class="section">
<h2>Sample Grouping</h2>
<div class="two-col">
  <div class="group-box group-benign">
    <h3 style="color:#1d4ed8;margin-top:0">&#x1F535; {g1_label} ({len(benign_samples)} samples)</h3>
    <p style="font-size:.9em">{g1_samples_str}</p>
  </div>
  <div class="group-box group-malignant">
    <h3 style="color:#b91c1c;margin-top:0">&#x1F534; {g2_label} ({len(malignant_samples)} samples)</h3>
    <p style="font-size:.9em">{g2_samples_str}</p>
  </div>
</div>
<p style="color:#666;font-size:.83em;margin-top:10px">Source: data/ovarian_fq_series_samples_v3.xlsx</p>
</div>

<!-- 5. Methylation Extraction -->
<div class="section">
<h2>Step 9. Methylation Extraction — MethylDackel extract</h2>
<table>
  <tr><th>Sample</th><th>Group</th><th>CpG Sites Detected</th><th>Mean Methylation</th></tr>
  {cpg_rows}
</table>
<table style="margin-top:16px">
  <tr><th>Metric</th><th>{g1_label} (merged)</th><th>{g2_label} (merged)</th></tr>
  <tr><td>Total CpG sites detected</td><td>{_fmt_num(n_cpg1)}</td><td>{_fmt_num(n_cpg2)}</td></tr>
  <tr><td>Mean methylation level</td><td>{avg_meth1:.1f}%</td><td>{avg_meth2:.1f}%</td></tr>
  <tr><td>Candidate region CpG sites</td><td>{_fmt_num(n_cand1)}</td><td>{_fmt_num(n_cand2)}</td></tr>
  <tr><td>Candidate region mean methylation</td><td>{avg_cand1:.1f}%</td><td>{avg_cand2:.1f}%</td></tr>
  <tr><td>Fully methylated sites (100%) in candidate regions</td>
      <td>{full1} ({_safe_pct(full1, n_cand1)})</td>
      <td>{full2} ({_safe_pct(full2, n_cand2)})</td></tr>
</table>
<p style="color:#666;font-size:.83em">Candidate regions: OC_tissue_88_regions.hg38 (88 ovarian cancer differential regions); extraction: --minDepth 1 --mergeContext</p>
</div>

<!-- 6. Methylation Visualization -->
<div class="section">
<h2>Step 12. Methylation Distribution Visualization</h2>
<div class="two-col">
  <div><h3>Density Plot</h3>{_img_tag(density_src, 'Density Plot')}</div>
  <div><h3>Position Scatter Plot</h3>{_img_tag(scatter_src, 'Position Scatter Plot')}</div>
</div>
<div style="margin-top:14px"><h3>Methylation Heatmap</h3>{_img_tag(heatmap_src, 'Methylation Heatmap')}</div>
</div>

<!-- 7. Differential Methylation Analysis -->
<div class="section">
<h2>Step 11. Differential Methylation Analysis — DMC, |Δβ| ≥ 0.2</h2>
<div class="two-col">
  <div>
    <table>
      <tr><th>Statistic</th><th>Value</th></tr>
      <tr><td>Total DMC sites</td><td><strong>{n_dmc}</strong></td></tr>
      <tr><td>{g1_label} methylation higher</td><td>{s1_higher} sites</td></tr>
      <tr><td>{g2_label} methylation higher</td><td>{s2_higher} sites</td></tr>
      <tr><td>Mean |Δβ|</td><td>{mean_diff:.3f}</td></tr>
      <tr><td>Max |Δβ|</td><td>{max_diff:.3f}</td></tr>
    </table>
  </div>
  <div>
    <table>
      <tr><th>Chromosome</th><th>DMC Count</th><th>Fraction</th></tr>
      {chrom_rows}
    </table>
  </div>
</div>
<h3 style="margin-top:20px">Top 20 DMC sites by |Δβ|</h3>
<table>
  <tr><th>Chromosome</th><th>Position</th><th>{g1_label} β</th><th>{g2_label} β</th><th>|Δβ|</th></tr>
  {top_diff_rows}
</table>
<p style="color:#666;font-size:.83em">Red = {g1_label} higher; Blue = {g2_label} higher</p>
</div>

{sec8_html}

{sec9_html}

{sec10_html}

{sec11_html}

{sec12_html}

{sec13_html}

{sec14_html}

{sec_imeta_html}

<!-- Step 20. Results Interpretation -->
<div class="section">
<h2>Step 20. Results Interpretation &amp; Quality Assessment</h2>
<div class="interp">
{interp_html}
</div>
</div>

<!-- 16. Analysis Pipeline -->
<div class="section">
<h2>Pipeline Overview</h2>
<div class="funnel">
  <div class="funnel-step">fastp QC<br/>UMI extract</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">bwameth<br/>alignment</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">UMI<br/>dedup</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">On-target<br/>filter</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">M-bias<br/>check</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">MethylDackel<br/>extract</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Sample<br/>grouping</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Candidate<br/>regions</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Diff.<br/>analysis</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Feature<br/>selection</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">ML<br/>classifiers</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Temporal<br/>validation</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">Report</div>
</div>
<table style="margin-top:14px">
  <tr><th>Step</th><th>Tool</th><th>Key Parameters</th></tr>
  <tr><td>Preprocessing (Step 4)</td><td>fastp</td><td>trim_front=12/12, trim_tail2=5, min_length=35, Q20, UMI R2 5' 6 N</td></tr>
  <tr><td>Alignment (Step 5)</td><td>bwameth + bwa-mem</td><td>-t 8, two-step (view + sort -m 8 G), reference hg38</td></tr>
  <tr><td>Deduplication (Step 6)</td><td>umi_tools dedup</td><td>--paired, UMI at R2 5' 6 N</td></tr>
  <tr><td>On-target filtering (Step 7)</td><td>bedtools intersect</td><td>panel_B10308.bed</td></tr>
  <tr><td>M-bias check (Step 8)</td><td>MethylDackel mbias</td><td>confirms trim_tail2 parameter</td></tr>
  <tr><td>Methylation extraction (Step 9)</td><td>MethylDackel extract</td><td>--minDepth 1 --mergeContext → CpG bedGraph; merge per group (Benign/Malignant)</td></tr>
  <tr><td>Filtering (Step 10)</td><td>custom script + bedtools</td><td>repeat masking, depth filter, focus on 88 candidate regions</td></tr>
  <tr><td>Differential analysis (Step 11)</td><td>custom script</td><td>threshold |Δβ| ≥ 0.2, coverage statistics (Wilcoxon)</td></tr>
  <tr><td>Visualization (Step 12)</td><td>matplotlib / seaborn</td><td>density, scatter, heatmap; candidate regions merged bedGraph</td></tr>
  <tr><td>Sample QC &amp; outlier detection (Step 13)</td><td>scikit-learn / scipy</td><td>coverage fraction, PCA z-score, Isolation Forest; ≥2 flags = suspect</td></tr>
  <tr><td>Feature selection (Step 14)</td><td>scikit-learn</td><td>LASSO (LogisticRegressionCV, L1) ∪ RF (Gini importance > 0); outlier samples excluded from selection</td></tr>
  <tr><td>Diagnostic classifiers (Step 15)</td><td>scikit-learn / XGBoost / LightGBM</td><td>LR / SVM / RF / XGB / LGBM, LOO-CV, StandardScaler per fold</td></tr>
  <tr><td>Subtype classifier (Step 16)</td><td>scikit-learn</td><td>RF + SVM OvR, LOO-CV, malignant + borderline samples</td></tr>
  <tr><td>Temporal validation (Step 17)</td><td>scikit-learn</td><td>Train ← before cutoff_date; Validation ← on/after; top-N freq features; no test-set leakage</td></tr>
  <tr><td>Site annotation (Step 18)</td><td>Ensembl REST + EWAS Catalog + UCSC CpGi</td><td>gene overlap, CpG island context, EWAS associations</td></tr>
  <tr><td>External validation (Step 19)</td><td>custom script</td><td>GSE40279 blood baseline + TCGA-OV tumour cross-reference</td></tr>
</table>
</div>

<p class="meta">Working dir: {temp_dir} &nbsp;|&nbsp; Result dir: {result_dir}<br/>Report generated: {analysis_date}</p>
</div>
</body>
</html>
"""

    # ── JSON ────────────────────────────────────────────────────
    report_data = {
        "analysis_date":   analysis_date,
        "samples":         all_samples,
        "benign_samples":  benign_samples,
        "malignant_samples": malignant_samples,
        "fastp":           fastp,
        "dedup":           dedup,
        "cpg_per_sample":  cpg_per_sample,
        "group_merged": {
            g1_label: {"cpg_total": n_cpg1, "avg_meth": avg_meth1, "candidate": n_cand1, "full_meth": full1},
            g2_label: {"cpg_total": n_cpg2, "avg_meth": avg_meth2, "candidate": n_cand2, "full_meth": full2},
        },
        "diff": {
            "n_dmc": n_dmc, "g1_higher": s1_higher, "g2_higher": s2_higher,
            "mean_diff": mean_diff, "max_diff": max_diff, "chrom_dist": chrom_dist,
        },
        "ml": clf_data,
    }
    json_report = os.path.join(report_dir, "analysis_report.json")
    with open(json_report, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2, ensure_ascii=False)

    html_report = os.path.join(report_dir, "analysis_report.html")
    with open(html_report, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info(f"Analysis report saved: {html_report}")
    logger.info(f"  Report directory (ready to zip): {report_dir}")
    return html_report


def _sample_name(methylation_path):
    """Extract sample name from methylation file path."""
    if not methylation_path:
        return "unknown"
    base = os.path.basename(methylation_path)
    base = re.sub(r"(candidate_filtered_|filtered_)", "", base)
    base = re.sub(r"_CpG\.bedGraph$", "", base)
    return base
