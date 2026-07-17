#!/usr/bin/env python3
"""report.annotation — Step 18 panel-CpG site annotation section.

Reads `panel_annotation.csv` (Ensembl + EWAS Catalog + EWAS Atlas + PubMed)
and renders a summary table + paginated full table.
"""

import os
import logging

from scripts.modules.report.helpers import _read_csv_safe, _pagination_widget

logger = logging.getLogger(__name__)


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

