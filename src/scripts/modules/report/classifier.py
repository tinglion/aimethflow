#!/usr/bin/env python3
"""report.classifier — Steps 14 (feature selection) through 17 (temporal validation).

Builds the largest collection of report fragments:
  - Section 9:  Feature Selection (Step 14) — `sec9_html`
  - Section 10: Diagnostic Classifiers (Step 15) — `sec10_html` + fs_benchmark + CA-125
  - Section 11: Subtype Classifier (Step 16) — `sec11_html` (incl. OvR subsection)
  - Section 12: Temporal Validation (Step 17) — `sec12_html`
The outlier section (Step 13) is built separately via `report.outlier`.
"""

import os
import json
import logging

from scripts.modules.report.helpers import (
    _read_csv_safe,
    _pagination_widget,
    _copy_img,
    _img_tag,
)
from scripts.modules.report.outlier import _build_outlier_section

logger = logging.getLogger(__name__)


def _build_ovr_subsection(ovr_records, roc_src, result_dir=None):
    """
    Return HTML fragment for the One-vs-Rest binary classifier subsection.

    When `result_dir` is given, also looks for:
      - `ovr_skip.json`            → renders a full Skipped panel (no records present)
      - `ovr_skipped_classes.json` → appends a footer listing per-class skips
        when the OvR step ran but excluded one or more subtypes
    """
    # Whole-step skip path — no records, but a skip marker explains why.
    if not ovr_records and result_dir:
        skip_path = os.path.join(result_dir, "ovr_skip.json")
        if os.path.exists(skip_path):
            try:
                import json as _json
                with open(skip_path, encoding="utf-8") as _f:
                    skip = _json.load(_f)
            except Exception:
                skip = {}
            thresholds = skip.get("thresholds", {})
            ovr_min    = thresholds.get("ovr_min_samples_per_class", "?")
            skipped    = skip.get("skipped_classes", []) or []
            n_obs      = skip.get("n_observed_classes", "?")
            sk_rows = "".join(
                f"<tr><td>{s.get('subtype','?')}</td>"
                f"<td style='text-align:right'>{s.get('n_positive','?')}</td>"
                f"<td style='text-align:right'>{s.get('n_negative','?')}</td></tr>"
                for s in skipped
            )
            return f"""
<hr style="margin:20px 0;border:none;border-top:1px solid #e0e0e0">
<h3 style="margin-top:0">11b. One-vs-Rest Binary Classifiers — Skipped</h3>
<p style="color:#c0392b;font-size:.95em;font-weight:bold">
  &#9888; OvR binary classification was skipped — no subtype met the configured
  per-class minimum of <code>{ovr_min}</code> positive AND <code>{ovr_min}</code> negative samples.
</p>
<p style="color:#666;font-size:.85em">
  Of {n_obs} observed subtype(s), {len(skipped)} failed the threshold:
</p>
<table style="font-size:.85em">
  <thead><tr><th>Subtype</th><th style='text-align:right'>n positive</th><th style='text-align:right'>n negative</th></tr></thead>
  <tbody>{sk_rows}</tbody>
</table>
<p style="color:#666;font-size:.83em;margin-top:8px">
  To enable this section: collect more samples for the under-represented subtypes,
  or lower the threshold in <code>config.toml [ml]</code>
  (e.g. <code>ovr_min_samples_per_class = 2</code>).
</p>"""
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

    # Partial-skip footer: when some subtypes were excluded from OvR, surface
    # the list so the table-row absence isn't silent.
    partial_html = ""
    if result_dir:
        partial_path = os.path.join(result_dir, "ovr_skipped_classes.json")
        if os.path.exists(partial_path):
            try:
                import json as _json
                with open(partial_path, encoding="utf-8") as _f:
                    sk = _json.load(_f)
            except Exception:
                sk = {}
            skipped = sk.get("skipped_classes", []) or []
            ovr_min = sk.get("thresholds", {}).get("ovr_min_samples_per_class", "?")
            if skipped:
                rows_sk = "".join(
                    f"<tr><td>{s.get('subtype','?')}</td>"
                    f"<td style='text-align:right'>{s.get('n_positive','?')}</td>"
                    f"<td style='text-align:right'>{s.get('n_negative','?')}</td></tr>"
                    for s in skipped
                )
                partial_html = f"""
<details style="margin-top:10px;font-size:.85em" open>
  <summary style="cursor:pointer;color:#856404;font-weight:bold">
    &#9888; {len(skipped)} subtype(s) excluded from OvR (per-class minimum
    <code>ml.ovr_min_samples_per_class = {ovr_min}</code> not met)
  </summary>
  <table style="font-size:.85em;margin-top:6px">
    <thead><tr><th>Subtype</th><th style='text-align:right'>n positive</th><th style='text-align:right'>n negative</th></tr></thead>
    <tbody>{rows_sk}</tbody>
  </table>
</details>"""

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
{partial_html}
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
        thresholds  = skip_info.get("thresholds", {}) or {}
        n_total     = skip_info.get("n_total")
        cls_dist    = skip_info.get("class_distribution", {}) or {}
        warn_html   = "".join(f"<li>{w}</li>" for w in skip_warns) if skip_warns else ""

        # Optional threshold table — only rendered when the new skip-marker
        # schema populated it (older skip.json files still display gracefully).
        threshold_rows = ""
        if thresholds:
            for key, val in thresholds.items():
                threshold_rows += (
                    f"<tr><td><code>{key}</code></td>"
                    f"<td style='text-align:right'>{val}</td></tr>"
                )
        threshold_html = (
            f"<h3 style='margin-top:14px;font-size:.95em'>Configured thresholds "
            f"(config.toml <code>[ml]</code>)</h3>"
            f"<table style='font-size:.85em'>"
            f"<thead><tr><th>Key</th><th style='text-align:right'>Current value</th></tr></thead>"
            f"<tbody>{threshold_rows}</tbody></table>"
            if threshold_rows else ""
        )

        # Class distribution table — shows exactly which subtypes were seen and
        # how many samples each had, so the user can decide whether to lower
        # the threshold or collect more samples.
        cls_rows = "".join(
            f"<tr><td>{c}</td><td style='text-align:right'>{n}</td></tr>"
            for c, n in sorted(cls_dist.items(), key=lambda kv: -kv[1])
        )
        cls_html = (
            f"<h3 style='margin-top:14px;font-size:.95em'>Cohort summary "
            f"(n = {n_total if n_total is not None else 'unknown'})</h3>"
            f"<table style='font-size:.85em'>"
            f"<thead><tr><th>Subtype (canonical)</th>"
            f"<th style='text-align:right'>Count</th></tr></thead>"
            f"<tbody>{cls_rows}</tbody></table>"
            if cls_rows else ""
        )

        sec11_html  = f"""
<div class="section">
<h2>Step 16. Subtype Classifier — Skipped</h2>
<p style="color:#c0392b;font-size:.95em;font-weight:bold">
  &#9888; Subtype classification was skipped because the cohort does not meet
  the configured minimum requirements for reliable multi-class LOO-CV.
</p>
<table>
  <tr><th>Reason</th><td>{skip_reason}</td></tr>
</table>
{"<ul style='color:#856404;font-size:.85em;margin-top:8px'>" + warn_html + "</ul>" if warn_html else ""}
{cls_html}
{threshold_html}
<p style="color:#666;font-size:.83em;margin-top:10px">
  <strong>Two ways to enable this section:</strong>
  (1) collect additional malignant samples with subtype labels and re-run
  <code>--step subtype_classify</code>, or
  (2) lower the threshold in <code>config.toml [ml]</code>
  (e.g. <code>subtype_min_samples_total = 5</code>) if you accept reduced
  statistical reliability for exploratory analysis.
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
{_build_ovr_subsection(hgsoc_records, hgsoc_roc_src, result_dir=result_dir)}
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
        # Build a quick lookup of full-cohort AUC by classifier for the Δ column in the subset table.
        full_auc_by_clf = {}
        for r in tv_records:
            try:
                full_auc_by_clf[r.get("classifier", "")] = float(r.get("auc", ""))
            except (TypeError, ValueError):
                pass
            tv_rows += (
                f"<tr><td>{r.get('classifier','')}</td>"
                f"<td>{r.get('n_train','')}</td>"
                f"<td>{r.get('n_validation','')}</td>"
                f"<td>{r.get('auc','')}</td>"
                f"<td>{r.get('accuracy','')}</td>"
                f"<td>{r.get('sensitivity','')}</td>"
                f"<td>{r.get('specificity','')}</td></tr>\n"
            )
        cutoff = tv_records[0].get("cutoff_date", "2024-01-01") if tv_records else "2024-01-01"

        # ── Low-CA-125 subset block (optional; only present when Step 17 was
        # run with low_ca125_threshold > 0 and the subset was viable).
        low_ca_path    = os.path.join(result_dir, "temporal_validation_summary_lowca125.csv")
        low_ca_records = _read_csv_safe(low_ca_path)
        if low_ca_records:
            first = low_ca_records[0]
            low_thresh   = first.get("low_ca125_threshold", "100")
            n_used       = first.get("n_validation_used", "")
            n_benign_s   = first.get("n_benign", "")
            n_malign_s   = first.get("n_malignant", "")
            low_rows = ""
            for r in low_ca_records:
                clf  = r.get("classifier", "")
                try:
                    auc_sub = float(r.get("auc", ""))
                except (TypeError, ValueError):
                    auc_sub = None
                full_auc = full_auc_by_clf.get(clf)
                if auc_sub is not None and full_auc is not None:
                    delta = auc_sub - full_auc
                    delta_color = "#c0392b" if delta <= -0.10 else ("#e67e22" if delta <= -0.04 else "#22863a")
                    delta_str = f"<span style='color:{delta_color};font-weight:bold'>{delta:+.3f}</span>"
                else:
                    delta_str = "—"
                full_auc_cell = f"{full_auc:.4f}" if full_auc is not None else "—"
                low_rows += (
                    f"<tr><td>{clf}</td>"
                    f"<td>{r.get('auc','')}</td>"
                    f"<td>{r.get('auc_ci95_lo','')}–{r.get('auc_ci95_hi','')}</td>"
                    f"<td>{full_auc_cell}</td>"
                    f"<td>{delta_str}</td>"
                    f"<td>{r.get('accuracy','')}</td>"
                    f"<td>{r.get('sensitivity','')}</td>"
                    f"<td>{r.get('specificity','')}</td></tr>\n"
                )
            low_ca_html = f"""
<h3 style="margin-top:18px">Low-CA-125 Subgroup (CA-125 &lt; {low_thresh} U/mL)</h3>
<p style="color:#666;font-size:.85em">
  Subset where the standard serum marker is least informative — the clinically most relevant zone for cfDNA.
  <strong>n = {n_used}</strong> ({n_benign_s} benign, {n_malign_s} malignant). Δ column compares each model's
  subset AUC against its full-validation AUC (above). <span style="color:#22863a">Green</span> = robust
  (Δ &gt; −0.04), <span style="color:#e67e22">orange</span> = moderate drop, <span style="color:#c0392b">red</span> = ≥ 0.10 AUC loss.
</p>
<table>
  <tr><th>Classifier</th><th>AUC (subset)</th><th>95% CI</th><th>AUC (full)</th>
      <th>Δ vs full</th><th>Accuracy</th><th>Sensitivity</th><th>Specificity</th></tr>
  {low_rows}
</table>
<p style="color:#666;font-size:.83em">
  Clinical interpretation: CA-125 alone typically drops sharply in this zone (it's defined by low CA-125),
  while cfDNA-only and linear-fusion (LR+CA-125) models tend to remain stable — the cfDNA panel
  effectively rescues diagnostic performance where CA-125 is uninformative.
  Output: <code>temporal_validation_summary_lowca125.csv</code>; threshold:
  <code>temporal_validate.low_ca125_threshold</code> in config.toml.
</p>"""
            clf_data["temporal_validation_lowca125"] = low_ca_records
        else:
            low_ca_html = ""
            clf_data["temporal_validation_lowca125"] = None

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
{low_ca_html}
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
