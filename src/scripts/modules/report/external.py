#!/usr/bin/env python3
"""report.external — Step 19 public-dataset cross-validation section
(blood baseline GSE40279 + TCGA-OV concordance).
"""

import os
import logging

from scripts.modules.report.helpers import _read_csv_safe, _pagination_widget, _safe_pct

logger = logging.getLogger(__name__)


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

