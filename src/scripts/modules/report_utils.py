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
  8.  Feature Selection (LASSO + RF) [optional]
  9.  Diagnostic Classifiers (LR / SVM / RF / XGB / LGBM, LOO-CV) [optional]
  10. Subtype Classifier (multi-class LOO-CV) [optional]
  11. Results Interpretation & Quality Assessment
  12. Analysis Pipeline
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


def _build_classifier_section(result_dir, imgs_dir):
    """
    Read Step 14-16 output files and build HTML fragments + JSON data.
    All three sections are optional: returns empty string / None when files are absent.
    Returns (sec9_html, sec10_html, sec11_html, clf_data_dict)
    """
    # ── Section 9: Feature Selection ──
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
            all_rows = "".join(
                f"<tr style='{'background:#e8f5e9' if r.get('pass_both') else ''}'>"
                f"<td>{r['feature']}</td>"
                f"<td style='text-align:right'>{float(r.get('benign_mean',0)):.3f}</td>"
                f"<td style='text-align:right'>{float(r.get('malignant_mean',0)):.3f}</td>"
                f"<td style='text-align:right'><strong>{float(r.get('sample_delta',0)):.3f}</strong></td>"
                f"<td style='text-align:right'>{float(r.get('wilcox_pval',1)):.4f}</td>"
                f"<td style='text-align:center'>{'✓' if r.get('pass_p') else '—'}</td>"
                f"<td style='text-align:center'>{'✓' if r.get('pass_delta') else '—'}</td>"
                f"<td style='text-align:center'>{'<strong style=color:#22863a>Pass</strong>' if r.get('pass_both') else '<span style=color:#aaa>Fail</span>'}</td>"
                f"</tr>"
                for r in wilcox_records   # already sorted by delta desc
            )

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
  <table style="font-size:.82em">
    <tr><th>Site</th><th>Benign mean</th><th>Malignant mean</th>
        <th>|Δβ| (sample)</th><th>p-value</th>
        <th>p OK</th><th>|Δβ| OK</th><th>Result</th></tr>
    {all_rows}
  </table>
  </div>
</details>
<p style="color:#666;font-size:.83em;margin-top:6px">
  Full stats: <code>wilcoxon_stats.csv</code> &nbsp;|&nbsp;
  Thresholds configurable in <code>config.toml</code> [ml] sample_wilcox_p / sample_delta_beta
</p>"""
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

        sec9_html = f"""
<div class="section">
<h2>9. Feature Selection (LASSO + RF)</h2>
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
{wilcox_section}
{depth_section}
</div>"""
        sec9_data = {"n_total": n_total, "lasso_selected": lasso_sel,
                     "rf_selected": rf_sel, "panel_size": n_panel,
                     "wilcoxon": wilcox_data, "depth": depth_data}
    else:
        sec9_html = ""
        sec9_data = None

    # ── Section 10: Diagnostic Classifiers ──
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
        sec10_html = f"""
<div class="section">
<h2>10. Diagnostic Classifiers (LOO-CV — Benign vs Malignant)</h2>
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
</div>"""
        sec10_data = clf_records
    else:
        sec10_html = ""
        sec10_data = None

    # ── Section 11: Subtype Classifier ──
    sub_sum_path = os.path.join(result_dir, "subtype_summary.csv")
    sub_roc_src  = _copy_img(os.path.join(result_dir, "subtype_roc.png"), imgs_dir)
    rf_cm_src    = _copy_img(os.path.join(result_dir, "subtype_confusion_RF.png"), imgs_dir)
    svm_cm_src   = _copy_img(os.path.join(result_dir, "subtype_confusion_SVM.png"), imgs_dir)
    sub_records  = _read_csv_safe(sub_sum_path)

    if sub_records is not None:
        sub_rows = ""
        for r in sub_records:
            sub_rows += (
                f"<tr><td>{r.get('classifier','')}</td>"
                f"<td>{r.get('accuracy', 0):.3f}</td>"
                f"<td>{r.get('n_samples', '')}</td>"
                f"<td>{r.get('n_subtypes', '')}</td></tr>\n"
            )
        sec11_html = f"""
<div class="section">
<h2>11. Subtype Classifier (LOO-CV — Malignant + Borderline, Multi-class)</h2>
<p style="color:#666;font-size:.85em">Subtypes: High-grade serous / Endometrioid / Clear cell / Low-grade serous / Mucinous</p>
<table>
  <tr><th>Classifier</th><th>Accuracy</th><th>Samples</th><th>Subtypes</th></tr>
  {sub_rows}
</table>
<div class="two-col" style="margin-top:14px">
  <div><h3>RF Confusion Matrix</h3>{_img_tag(rf_cm_src, 'RF Confusion Matrix')}</div>
  <div><h3>SVM Confusion Matrix</h3>{_img_tag(svm_cm_src, 'SVM Confusion Matrix')}</div>
</div>
<div style="margin-top:14px"><h3>OvR ROC Curves</h3>{_img_tag(sub_roc_src, 'Subtype ROC')}</div>
<p style="color:#666;font-size:.83em">Strategy: One-vs-Rest (OvR); evaluation: LOO-CV;
outputs: <code>subtype_summary.csv</code>, <code>subtype_predictions.csv</code></p>
</div>"""
        sec11_data = sub_records
    else:
        sec11_html = ""
        sec11_data = None

    clf_data = {
        "feature_selection":  sec9_data,
        "classifier":         sec10_data,
        "subtype_classifier": sec11_data,
    }
    return sec9_html, sec10_html, sec11_html, clf_data


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
    """Generate full analysis report (HTML + JSON)."""
    logger.info("开始生成分析报告")

    benign_samples    = benign_samples or []
    malignant_samples = malignant_samples or []

    # ── Image output directory ────────────────────────────────────
    imgs_dir = os.path.join(result_dir, "imgs")
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

    # Optional classifier sections (Steps 14-16; empty string if files absent)
    sec9_html, sec10_html, sec11_html, clf_data = _build_classifier_section(result_dir, imgs_dir)

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
  <div class="card"><div class="num">{mean_diff:.3f}</div><div class="lbl">Mean |Δβ|</div></div>
</div>

<!-- 1. Data Quality Overview -->
<div class="section">
<h2>1. Data Quality Overview (fastp QC)</h2>
<table>
  <tr><th>Sample</th><th>Group</th><th>Raw Read Pairs</th><th>Filtered Pairs</th><th>Pass Rate</th><th>Q30 (post-QC)</th></tr>
  {qc_rows}
</table>
<p style="color:#666;font-size:.83em">Parameters: UMI extraction (R2 5' 6 N) + trim_front1=12 + trim_front2=12 + trim_tail2=5 + min_length=35 + Q20</p>
</div>

<!-- 2. Alignment & Deduplication -->
<div class="section">
<h2>2. Alignment (bwameth) &amp; UMI Deduplication</h2>
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
<h2>3. M-bias Check (MethylDackel mbias)</h2>
<p style="color:#666;font-size:.85em">OT = Watson-strand methylation (~99% of reads). If the R2 3' curve drops, increase trim_tail2 and re-run from --step extract.</p>
{mbias_body}
</div>

<!-- 4. Sample Grouping -->
<div class="section">
<h2>4. Sample Grouping</h2>
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
<p style="color:#666;font-size:.83em;margin-top:10px">Source: data/卵巢癌fq编号与样本信息对应表.xlsx</p>
</div>

<!-- 5. Methylation Extraction -->
<div class="section">
<h2>5. Methylation Extraction (MethylDackel extract)</h2>
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
<h2>6. Methylation Distribution Visualization (candidate regions, merged)</h2>
<div class="two-col">
  <div><h3>Density Plot</h3>{_img_tag(density_src, 'Density Plot')}</div>
  <div><h3>Position Scatter Plot</h3>{_img_tag(scatter_src, 'Position Scatter Plot')}</div>
</div>
<div style="margin-top:14px"><h3>Methylation Heatmap</h3>{_img_tag(heatmap_src, 'Methylation Heatmap')}</div>
</div>

<!-- 7. Differential Methylation Analysis -->
<div class="section">
<h2>7. Differential Methylation Analysis (DMC, |Δβ| ≥ 0.2)</h2>
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

{sec9_html}

{sec10_html}

{sec11_html}

<!-- 12. Results Interpretation -->
<div class="section">
<h2>12. Results Interpretation &amp; Quality Assessment</h2>
<div class="interp">
{interp_html}
</div>
</div>

<!-- 13. Analysis Pipeline -->
<div class="section">
<h2>13. Analysis Pipeline</h2>
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
  <div class="funnel-step">Report</div>
</div>
<table style="margin-top:14px">
  <tr><th>Step</th><th>Tool</th><th>Key Parameters</th></tr>
  <tr><td>Preprocessing</td><td>fastp</td><td>trim_front=12/12, trim_tail2=5, min_length=35, Q20, UMI R2 5' 6 N</td></tr>
  <tr><td>Alignment</td><td>bwameth + bwa-mem</td><td>-t 8, two-step (view + sort -m 8 G), reference hg38</td></tr>
  <tr><td>Deduplication</td><td>umi_tools dedup</td><td>--paired, UMI at R2 5' 6 N</td></tr>
  <tr><td>On-target filtering</td><td>bedtools intersect</td><td>panel_B10308.bed</td></tr>
  <tr><td>M-bias check</td><td>MethylDackel mbias</td><td>confirms trim_tail2 parameter</td></tr>
  <tr><td>Methylation extraction</td><td>MethylDackel extract</td><td>--minDepth 1 --mergeContext → CpG bedGraph</td></tr>
  <tr><td>Group merging</td><td>custom script</td><td>split by Excel into Benign/Malignant; merge bedGraphs per group</td></tr>
  <tr><td>Filtering</td><td>custom script + bedtools</td><td>repeat masking, depth filter, focus on 88 candidate regions</td></tr>
  <tr><td>Differential analysis</td><td>custom script</td><td>threshold |Δβ| ≥ 0.2, coverage statistics (Wilcoxon)</td></tr>
  <tr><td>Feature selection</td><td>scikit-learn</td><td>LASSO (LogisticRegressionCV, L1) ∪ RF (Gini importance > 0)</td></tr>
  <tr><td>Diagnostic classifiers</td><td>scikit-learn / XGBoost / LightGBM</td><td>LR / SVM / RF / XGB / LGBM, LOO-CV, StandardScaler per fold</td></tr>
  <tr><td>Subtype classifier</td><td>scikit-learn</td><td>RF + SVM OvR, LOO-CV, malignant + borderline samples</td></tr>
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
    json_report = os.path.join(result_dir, "analysis_report.json")
    with open(json_report, "w", encoding="utf-8") as f:
        json.dump(report_data, f, indent=2, ensure_ascii=False)

    html_report = os.path.join(result_dir, "analysis_report.html")
    with open(html_report, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info(f"分析报告生成完成: {html_report}")
    return html_report


def _sample_name(methylation_path):
    """Extract sample name from methylation file path."""
    if not methylation_path:
        return "unknown"
    base = os.path.basename(methylation_path)
    base = re.sub(r"(candidate_filtered_|filtered_)", "", base)
    base = re.sub(r"_CpG\.bedGraph$", "", base)
    return base
