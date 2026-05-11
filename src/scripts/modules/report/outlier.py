#!/usr/bin/env python3
"""report.outlier — Step 13 sample-quality / outlier-detection section.

Reads `outlier_report.csv` + (optionally) embeds the PCA scatter plot,
and renders a paginated per-sample table.
"""

import os
import logging

from scripts.modules.report.helpers import (
    _read_csv_safe,
    _pagination_widget,
    _copy_img,
    _img_tag,
)

logger = logging.getLogger(__name__)


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

