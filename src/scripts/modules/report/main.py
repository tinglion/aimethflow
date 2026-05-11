#!/usr/bin/env python3
"""report.main — top-level `generate_report` entry point + the giant HTML f-string.

Aggregates every section builder in the report package (helpers, progress,
outlier, classifier, annotation, external) and writes the final HTML + JSON
output to `result/<wd>/<wd>_report/`.
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

# Section builders + shared helpers
from scripts.modules.report.helpers import (
    _parse_fastp_json, _parse_dedup_stats, _count_bedgraph_sites,
    _parse_diff_results,
    _fmt_num, _fmt_pct, _fmt_mb, _fmt_gb, _safe_pct,
    _copy_img, _img_tag, _read_csv_safe, _pagination_widget,
)
from scripts.modules.report.progress    import _build_pipeline_progress_section
from scripts.modules.report.outlier     import _build_outlier_section
from scripts.modules.report.classifier  import (
    _build_classifier_section, _build_ovr_subsection,
)
from scripts.modules.report.annotation  import _build_annotation_section
from scripts.modules.report.external    import _build_external_validation_section

logger = logging.getLogger(__name__)


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

    # Fallback: if no fastp.json files (workdir cleaned, regenerated standalone, or
    # different temp_dir), use grouped sample IDs from pipeline state. Without this
    # the "Total Samples" card shows 0 even when classifiers ran on 116 samples.
    grouped_ids = sorted(set((benign_samples or []) + (malignant_samples or [])))
    if not all_samples and grouped_ids:
        logger.info(f"  No fastp.json found in {temp_dir} — falling back to "
                    f"{len(grouped_ids)} grouped sample IDs from pipeline state")
        all_samples = grouped_ids

    # Total-samples count for the header card: union of fastp-discovered and
    # grouped IDs. Robust whether or not the workdir has been cleaned.
    n_total_samples = len(set(all_samples) | set(grouped_ids))

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
        f"<li>Analysis includes <strong>{n_total_samples}</strong> samples: "
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

    # QC table rows (paginated — see _pagination_widget below)
    qc_rows = ""
    for idx, s in enumerate(all_samples):
        fp = fastp.get(s, {})
        pairs_in   = fp.get("total_reads_before", 0) // 2
        pairs_out  = fp.get("total_reads_after",  0) // 2
        q30        = fp.get("q30_after", 0)
        pass_rate  = pairs_out / pairs_in if pairs_in else 0
        grp_badge  = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                      f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                      "<span class='badge badge-grey'>Ungrouped</span>")
        qc_rows += (
            f"<tr data-row='{idx}'><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(pairs_in)}</td><td>{_fmt_num(pairs_out)}</td>"
            f"<td>{_fmt_pct(pass_rate)}</td><td>{_fmt_pct(q30)}</td></tr>\n"
        )
    QC_PAGE_SIZE = 10
    qc_controls, qc_js = _pagination_widget("qc_tbl", "qc_ctrl",
                                            len(all_samples), QC_PAGE_SIZE)

    # Alignment / dedup / on-target table rows (paginated)
    aln_rows = ""
    for idx, s in enumerate(all_samples):
        dp = dedup.get(s, {})
        inp  = dp.get("input_pairs", 0)
        out  = dp.get("output_reads", 0) // 2
        dr   = dp.get("dup_rate", None)
        ont_r = (bam_target[s] / bam_dedup[s]) if bam_dedup[s] else None
        aln_rows += (
            f"<tr data-row='{idx}'><td>{s}</td>"
            f"<td>{_fmt_gb(bam_sorted[s])}</td>"
            f"<td>{_fmt_num(inp)}</td>"
            f"<td>{_fmt_num(out)}</td>"
            f"<td>{_fmt_pct(dr)}</td>"
            f"<td>{_fmt_mb(bam_target[s])}</td>"
            f"<td>{_fmt_pct(ont_r)}</td></tr>\n"
        )
    aln_controls, aln_js = _pagination_widget("aln_tbl", "aln_ctrl",
                                              len(all_samples), QC_PAGE_SIZE)

    # Per-sample CpG table rows (paginated)
    cpg_rows = ""
    for idx, s in enumerate(all_samples):
        c = cpg_per_sample.get(s, {})
        grp_badge = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                     f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                     "<span class='badge badge-grey'>Ungrouped</span>")
        cpg_rows += (
            f"<tr data-row='{idx}'><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(c.get('count', 0))}</td>"
            f"<td>{c.get('avg', 0):.1f}%</td></tr>\n"
        )
    cpg_controls, cpg_js = _pagination_widget("cpg_tbl", "cpg_ctrl",
                                              len(all_samples), QC_PAGE_SIZE)

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

    # Pipeline progress (Steps 1–20 status from pipeline_state.json)
    progress_html = _build_pipeline_progress_section(temp_dir)

    # ── Trim parameters actually applied (for Step 12 visualization context) ──
    #
    # mbias_auto_trim.json is written by --step mbias and consumed by --step extract
    # when nOT/nOB are set to "auto" in config.toml. Surfacing the resolved skip
    # values in the report makes it explicit which positional filtering produced
    # the methylation distribution shown below.
    trim_params_html = ""
    try:
        _trim_path = os.path.join(temp_dir, "mbias_auto_trim.json")
        if os.path.exists(_trim_path):
            with open(_trim_path, encoding="utf-8") as _f:
                _trim_info = json.load(_f)
            _n_ot     = _trim_info.get("nOT", "—")
            _n_ob     = _trim_info.get("nOB", "—")
            _rlen     = _trim_info.get("read_len", "—")
            _n_smp    = _trim_info.get("n_samples", "—")
            _pct      = _trim_info.get("percentile_used", "—")
            trim_params_html = f"""
<details open style="margin:6px 0 14px;border:1px solid #ddd;border-radius:6px;padding:10px 14px;background:#fafbfc">
  <summary style="cursor:pointer;font-weight:bold;color:#1a6fbb;font-size:.92em">
    Trim parameters actually applied (data-derived; auto-resolved from M-bias)
  </summary>
  <table style="font-size:.85em;margin-top:8px">
    <thead><tr><th>Stage / parameter</th><th>Value applied</th><th>Notes</th></tr></thead>
    <tbody>
      <tr><td>fastp <code>trim_front1</code>, <code>trim_front2</code>, <code>trim_tail2</code></td>
          <td><code>12, 12, 5</code></td>
          <td>R1/R2 5′ adapter + R2 3′ poly-C (Step 4; values from <code>config.toml [plasma]</code>)</td></tr>
      <tr><td>MethylDackel <code>--nOT</code> skip counts</td>
          <td><code>{_n_ot}</code></td>
          <td>OT strand 5′,3′,3′,5′ — auto-resolved from M-bias</td></tr>
      <tr><td>MethylDackel <code>--nOB</code> skip counts</td>
          <td><code>{_n_ob}</code></td>
          <td>OB strand 5′,3′,3′,5′ — auto-resolved from M-bias</td></tr>
      <tr><td>Detected read length</td>
          <td><code>{_rlen}</code> bp</td>
          <td>used to convert MethylDackel inclusion positions → skip counts</td></tr>
      <tr><td>Aggregation across samples</td>
          <td>{_pct}th percentile of {_n_smp} sample(s)</td>
          <td><code>mbias_trim_percentile</code> in <code>config.toml [plasma]</code></td></tr>
    </tbody>
  </table>
  <p style="color:#666;font-size:.8em;margin:6px 0 0">
    Full per-sample distribution: <code>{temp_dir}/mbias_auto_trim.json</code>.
    The methylation distribution shown below is the outcome of applying these positional skips.
  </p>
</details>"""
        else:
            trim_params_html = (
                "<p style='color:#856404;background:#fff3cd;padding:8px 12px;"
                "border-radius:6px;font-size:.85em;margin:6px 0 14px'>"
                "&#9888; <code>mbias_auto_trim.json</code> not found — "
                "<code>--step mbias</code> has not been run, or "
                "<code>nOT</code>/<code>nOB</code> were set to explicit numeric values "
                "in <code>config.toml [plasma]</code> (bypassing the auto-trim path).</p>"
            )
    except Exception as _exc:
        logger.warning(f"  [report] Could not read mbias_auto_trim.json: {_exc}")
        trim_params_html = ""

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
  {n_total_samples} samples ({g1_label}: {len(benign_samples)} / {g2_label}: {len(malignant_samples)})</span>
</h1>

<!-- Key metric cards -->
<div class="cards">
  <div class="card"><div class="num">{n_total_samples}</div><div class="lbl">Total Samples</div></div>
  <div class="card"><div class="num">{len(benign_samples)}</div><div class="lbl">Benign</div></div>
  <div class="card"><div class="num">{len(malignant_samples)}</div><div class="lbl">Malignant</div></div>
  <div class="card"><div class="num">{_fmt_num(n_cand1 + n_cand2)}</div><div class="lbl">Candidate Region CpGs</div></div>
  <div class="card"><div class="num">{n_dmc}</div><div class="lbl">DMC Sites</div></div>
</div>

{progress_html}

<!-- 1. Data Quality Overview -->
<div class="section">
<h2>Step 4. Data Quality Overview — fastp QC</h2>
<div style="overflow-x:auto">
{qc_controls}
<table id="qc_tbl">
  <thead><tr><th>Sample</th><th>Group</th><th>Raw Read Pairs</th><th>Filtered Pairs</th><th>Pass Rate</th><th>Q30 (post-QC)</th></tr></thead>
  <tbody>
  {qc_rows}
  </tbody>
</table>
</div>
<p style="color:#666;font-size:.83em">Parameters: UMI extraction (R2 5' 6 N) + trim_front1=12 + trim_front2=12 + trim_tail2=5 + min_length=35 + Q20</p>
{qc_js}
</div>

<!-- 2. Alignment & Deduplication -->
<div class="section">
<h2>Steps 5–7. Alignment, Deduplication &amp; On-target Filter</h2>
<div style="overflow-x:auto">
{aln_controls}
<table id="aln_tbl">
  <thead><tr><th>Sample</th><th>Sorted BAM</th><th>Input Pairs (pre-dedup)</th><th>Dedup Pairs</th><th>Dup Rate</th><th>On-target BAM</th><th>On-target Rate</th></tr></thead>
  <tbody>
  {aln_rows}
  </tbody>
</table>
</div>
<p style="color:#666;font-size:.83em">Dedup: umi_tools dedup --paired; On-target: bedtools intersect panel_B10308.bed; BAM size of 0 means the file was auto-cleaned (expected).</p>
{aln_js}
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
<div style="overflow-x:auto">
{cpg_controls}
<table id="cpg_tbl">
  <thead><tr><th>Sample</th><th>Group</th><th>CpG Sites Detected</th><th>Mean Methylation</th></tr></thead>
  <tbody>
  {cpg_rows}
  </tbody>
</table>
{cpg_js}
</div>
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

{trim_params_html}

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
