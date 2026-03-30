#!/usr/bin/env python3
"""
report_utils.py - 报告生成模块

报告结构：
  1. 样本概览（所有样本 fastp QC）
  2. 比对与去重（所有样本）
  3. M-bias 检查（SVG 内嵌）
  4. On-target 过滤（所有样本）
  5. 样本分组（良性 vs 恶性）
  6. 甲基化提取与候选区域
  7. 甲基化分布可视化
  8. 差异甲基化分析（DMC）
  9. 结果解读与质量评价
 10. 分析流程
"""

import os
import re
import json
import glob
import base64
import logging
from datetime import datetime

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────
#  数据读取辅助函数
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
        logger.warning(f"无法解析 fastp JSON {path}: {e}")
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
        # umi_tools --paired 模式下：Input Reads = Read pairs（单位都是 pair）
        # Number of reads out 也是 pairs 数，因此去重率 = 1 - out_pairs / in_pairs
        if "input_pairs" in result and "output_reads" in result and result["input_pairs"] > 0:
            result["dup_rate"] = 1.0 - result["output_reads"] / result["input_pairs"]
    except Exception as e:
        logger.warning(f"无法解析 dedup stats {path}: {e}")
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
        logger.warning(f"无法解析 diff 结果 {path}: {e}")
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
    """安全百分比：分母为 0 时返回 '—'。"""
    if not denom:
        return "—"
    return f"{num / denom * 100:.{digits}f}%"


def _embed_image_b64(path):
    if not path or not os.path.exists(path):
        return None
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    mime = {"png": "image/png", "jpg": "image/jpeg", "jpeg": "image/jpeg"}.get(ext, "image/png")
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:{mime};base64,{data}"


def _embed_svg(path):
    if not path or not os.path.exists(path):
        return ""
    try:
        with open(path, encoding="utf-8") as f:
            content = f.read()
        content = re.sub(r"<\?xml[^?]*\?>", "", content).strip()
        return content
    except Exception:
        return ""


def _img_tag(b64, alt):
    if not b64:
        return f"<p style='color:#aaa;font-size:.85em'>{alt} 尚未生成，请先运行 --step visualize</p>"
    return f"<img src='{b64}' alt='{alt}' style='max-width:100%;border-radius:6px;margin-top:6px'/>"


# ─────────────────────────────────────────────
#  主报告生成函数
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
    """生成完整分析报告（HTML + JSON）。"""
    logger.info("开始生成分析报告")

    benign_samples    = benign_samples or []
    malignant_samples = malignant_samples or []

    # ── 自动发现所有样本（按 fastp JSON 枚举） ──────────────
    all_fastp_files = sorted(glob.glob(os.path.join(temp_dir, "*.fastp.json")))
    all_samples     = [os.path.basename(f).replace(".fastp.json", "") for f in all_fastp_files]

    # ── 每样本 QC 统计 ──────────────────────────────────────
    fastp  = {s: _parse_fastp_json(os.path.join(temp_dir, f"{s}.fastp.json")) for s in all_samples}
    dedup  = {s: _parse_dedup_stats(os.path.join(temp_dir, f"{s}.bwameth_pe.dedup.stats")) for s in all_samples}

    def _bam(s, suffix):
        p = os.path.join(temp_dir, f"{s}{suffix}")
        return os.path.getsize(p) if os.path.exists(p) else 0

    bam_sorted  = {s: _bam(s, ".bwameth_pe.sorted.bam")                for s in all_samples}
    bam_dedup   = {s: _bam(s, ".bwameth_pe.dedup.bam")                 for s in all_samples}
    bam_target  = {s: _bam(s, ".bwameth_pe.dedup.targeted.sorted.bam") for s in all_samples}

    # ── 每样本 CpG 位点统计（个人 bedGraph） ────────────────
    cpg_per_sample = {}
    for s in all_samples:
        bg = os.path.join(temp_dir, f"{s}_CpG.bedGraph")
        cnt, avg = _count_bedgraph_sites(bg) if os.path.exists(bg) else (0, 0.0)
        cpg_per_sample[s] = {"count": cnt, "avg": avg}

    # ── 分组级 CpG 与候选区域统计 ───────────────────────────
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

    # ── 差异甲基化结果 ───────────────────────────────────────
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

    # ── 可视化图表 ────────────────────────────────────────────
    density_b64 = _embed_image_b64(os.path.join(result_dir, "density_plot.png"))
    scatter_b64 = _embed_image_b64(os.path.join(result_dir, "scatter_plot.png"))
    heatmap_b64 = _embed_image_b64(os.path.join(result_dir, "heatmap_plot.png"))

    # ── M-bias SVG ────────────────────────────────────────────
    mbias_blocks = []
    for svg_path in sorted(glob.glob(os.path.join(temp_dir, "*.mbias_OT.svg"))):
        sname = os.path.basename(svg_path).replace(".mbias_OT.svg", "")
        ot = _embed_svg(svg_path)
        ob = _embed_svg(svg_path.replace("_OT.svg", "_OB.svg"))
        if ot:
            ot_b = f"<div style='flex:1;min-width:0'><p style='font-size:.8em;color:#555;text-align:center'>OT链（Watson）</p>{ot}</div>"
            ob_b = (f"<div style='flex:1;min-width:0'><p style='font-size:.8em;color:#555;text-align:center'>OB链（Crick）</p>{ob}</div>" if ob else "")
            mbias_blocks.append(
                f"<h3 style='font-size:.95em;color:#555;margin:16px 0 6px'>{sname}</h3>"
                f"<div style='display:flex;gap:10px;flex-wrap:wrap'>{ot_b}{ob_b}</div>"
            )
    mbias_body = "\n".join(mbias_blocks) if mbias_blocks else "<p style='color:#aaa'>M-bias SVG 尚未生成，请先运行 --step mbias</p>"

    # ── 分组标签 ────────────────────────────────────────────
    g1_label = "良性组（Benign）"
    g2_label = "恶性组（Malignant）"
    g1_samples_str = "、".join(benign_samples)    if benign_samples    else "—"
    g2_samples_str = "、".join(malignant_samples) if malignant_samples else "—"

    # ── 结果解读 ────────────────────────────────────────────
    interp_items = []

    # 样本规模
    n_benign    = len(benign_samples)
    n_malignant = len(malignant_samples)
    interp_items.append(
        f"<li>本次分析共 <strong>{len(all_samples)}</strong> 个样本，"
        f"其中良性 <strong>{n_benign}</strong> 例、恶性 <strong>{n_malignant}</strong> 例。</li>"
    )

    # 数据质量
    q30_vals = [fastp[s].get("q30_after", 0) for s in all_samples if fastp[s]]
    avg_q30  = sum(q30_vals) / len(q30_vals) if q30_vals else 0
    q30_eval = "优良（≥85%）" if avg_q30 >= 0.85 else "尚可（70–85%）" if avg_q30 >= 0.70 else "偏低（<70%），建议复检测序质量"
    interp_items.append(f"<li>平均 Q30 率 <strong>{avg_q30*100:.1f}%</strong>，数据质量{q30_eval}。</li>")

    # 去重率
    dup_vals = [dedup[s].get("dup_rate", None) for s in all_samples if dedup[s].get("dup_rate") is not None]
    if dup_vals:
        avg_dup = sum(dup_vals) / len(dup_vals)
        dup_eval = ("偏高（>60%），建议检查建库过程 PCR 循环数" if avg_dup > 0.60
                    else "正常（30–60%）" if avg_dup >= 0.30
                    else "偏低（<30%），数据利用率高")
        interp_items.append(f"<li>平均 UMI 去重率 <strong>{avg_dup*100:.1f}%</strong>，{dup_eval}。</li>")

    # on-target
    ont_vals = []
    for s in all_samples:
        if bam_dedup[s] and bam_target[s]:
            ont_vals.append(bam_target[s] / bam_dedup[s])
    if ont_vals:
        avg_ont = sum(ont_vals) / len(ont_vals)
        ont_eval = "良好（>30%）" if avg_ont > 0.30 else "偏低（<30%），可能存在脱靶或捕获效率问题"
        interp_items.append(f"<li>平均 On-target 率（文件大小估算）<strong>{avg_ont*100:.1f}%</strong>，{ont_eval}。</li>")

    # DMC 结果
    if n_dmc == 0:
        interp_items.append("<li>未发现达到阈值（|Δβ|≥0.2）的差异甲基化位点，两组甲基化特征相似，或样本量不足导致统计效力不足。</li>")
    elif n_dmc < 10:
        interp_items.append(f"<li>发现 <strong>{n_dmc}</strong> 个差异甲基化位点，数量较少，建议结合临床意义判断。</li>")
    elif n_dmc < 50:
        interp_items.append(f"<li>发现 <strong>{n_dmc}</strong> 个差异甲基化位点，两组间存在一定甲基化差异。</li>")
    else:
        interp_items.append(f"<li>发现 <strong>{n_dmc}</strong> 个差异甲基化位点，两组甲基化特征差异明显，具有较好的区分潜力。</li>")

    if n_dmc > 0:
        direction = g2_label if s2_higher > s1_higher else g1_label
        interp_items.append(
            f"<li><strong>{direction}</strong> 整体甲基化水平偏高"
            f"（{max(s1_higher, s2_higher)}/{n_dmc} 个位点）；"
            f"平均 |Δβ| = {mean_diff:.3f}，最大 |Δβ| = {max_diff:.3f}。</li>"
        )
        if max_diff >= 0.5:
            interp_items.append("<li>存在 |Δβ|≥0.5 的强差异位点，可作为候选生物标志物优先验证。</li>")

    # 候选区域覆盖
    if n_cand1 > 0 or n_cand2 > 0:
        interp_items.append(
            f"<li>88 个候选区域内，良性组检出 <strong>{_fmt_num(n_cand1)}</strong> 个 CpG 位点，"
            f"恶性组 <strong>{_fmt_num(n_cand2)}</strong> 个；覆盖率"
            f"{'充分（建议继续后续模型训练）' if min(n_cand1, n_cand2) > 100 else '有限（可考虑降低深度阈值）'}。</li>"
        )

    interp_html = "<ul style='line-height:1.9'>" + "".join(interp_items) + "</ul>"

    # ── 构建 HTML ── ──────────────────────────────────────────

    analysis_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # 所有样本 QC 表格行
    qc_rows = ""
    for s in all_samples:
        fp = fastp.get(s, {})
        pairs_in   = fp.get("total_reads_before", 0) // 2
        pairs_out  = fp.get("total_reads_after",  0) // 2
        q30        = fp.get("q30_after", 0)
        pass_rate  = pairs_out / pairs_in if pairs_in else 0
        grp_badge  = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                      f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                      "<span class='badge badge-grey'>未分组</span>")
        qc_rows += (
            f"<tr><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(pairs_in)}</td><td>{_fmt_num(pairs_out)}</td>"
            f"<td>{_fmt_pct(pass_rate)}</td><td>{_fmt_pct(q30)}</td></tr>\n"
        )

    # 所有样本比对/去重/on-target 表格行
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

    # 每样本 CpG 表格行
    cpg_rows = ""
    for s in all_samples:
        c = cpg_per_sample.get(s, {})
        grp_badge = (f"<span class='badge badge-blue'>{g1_label}</span>"  if s in benign_samples else
                     f"<span class='badge badge-red'>{g2_label}</span>"   if s in malignant_samples else
                     "<span class='badge badge-grey'>未分组</span>")
        cpg_rows += (
            f"<tr><td>{s}</td><td>{grp_badge}</td>"
            f"<td>{_fmt_num(c.get('count', 0))}</td>"
            f"<td>{c.get('avg', 0):.1f}%</td></tr>\n"
        )

    # 染色体分布表
    chrom_rows = "".join(
        f"<tr><td>{c}</td><td>{n}</td><td>{n/n_dmc*100:.1f}%</td></tr>"
        for c, n in sorted(chrom_dist.items())
    ) if n_dmc > 0 else "<tr><td colspan='3'>无差异位点</td></tr>"

    # 前 20 差异位点
    top_diff = sorted(diff_rows, key=lambda x: -x["diff"])[:20]
    top_diff_rows = "".join(
        f"<tr><td>{r['chrom']}</td><td>{r['pos']:,}</td>"
        f"<td>{r['meth1']:.1f}%</td><td>{r['meth2']:.1f}%</td>"
        f"<td class='{'up' if r['meth1']>r['meth2'] else 'down'}'>{r['diff']:.3f}</td></tr>"
        for r in top_diff
    )

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8"/>
<title>血浆甲基化分析报告</title>
<style>
  * {{ box-sizing: border-box; }}
  body {{ font-family: "Microsoft YaHei", Arial, sans-serif; margin: 0; background: #f4f6f9; color: #333; }}
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

<h1>血浆 cfDNA 靶向甲基化分析报告
  <br/><span>分析日期：{analysis_date} &nbsp;|&nbsp;
  {len(all_samples)} 个样本（良性 {len(benign_samples)} 例 / 恶性 {len(malignant_samples)} 例）</span>
</h1>

<!-- 关键指标卡片 -->
<div class="cards">
  <div class="card"><div class="num">{len(all_samples)}</div><div class="lbl">样本总数</div></div>
  <div class="card"><div class="num">{len(benign_samples)}</div><div class="lbl">良性样本</div></div>
  <div class="card"><div class="num">{len(malignant_samples)}</div><div class="lbl">恶性样本</div></div>
  <div class="card"><div class="num">{_fmt_num(n_cand1 + n_cand2)}</div><div class="lbl">候选区域总 CpG</div></div>
  <div class="card"><div class="num">{n_dmc}</div><div class="lbl">差异甲基化位点 (DMC)</div></div>
  <div class="card"><div class="num">{mean_diff:.3f}</div><div class="lbl">平均甲基化差异 |Δβ|</div></div>
</div>

<!-- 1. 数据质量概览 -->
<div class="section">
<h2>1. 数据质量概览（fastp QC）</h2>
<table>
  <tr><th>样本</th><th>分组</th><th>原始 reads 对数</th><th>质控后对数</th><th>通过率</th><th>Q30（质控后）</th></tr>
  {qc_rows}
</table>
<p style="color:#666;font-size:.83em">参数：UMI提取（R2 5' 6N）+ trim_front1=12 + trim_front2=12 + trim_tail2=5 + min_length=35 + Q20</p>
</div>

<!-- 2. 比对与去重 -->
<div class="section">
<h2>2. 比对（bwameth）与 UMI 去重</h2>
<table>
  <tr><th>样本</th><th>比对 BAM 大小</th><th>输入 reads 对数（去重前）</th><th>去重后 pairs</th><th>去重率</th><th>On-target BAM</th><th>On-target 率</th></tr>
  {aln_rows}
</table>
<p style="color:#666;font-size:.83em">去重：umi_tools dedup --paired；On-target：bedtools intersect panel_B10308.bed；BAM 大小为 0 表示已被自动清理（正常）。</p>
</div>

<!-- 3. M-bias 检查 -->
<div class="section">
<h2>3. M-bias 检查（MethylDackel mbias）</h2>
<p style="color:#666;font-size:.85em">OT = 正链甲基化（Watson，占 ~99%）；若 R2 3' 端曲线下滑，需调大 trim_tail2 重跑。</p>
{mbias_body}
</div>

<!-- 4. 样本分组 -->
<div class="section">
<h2>4. 样本分组</h2>
<div class="two-col">
  <div class="group-box group-benign">
    <h3 style="color:#1d4ed8;margin-top:0">🔵 {g1_label}（{len(benign_samples)} 例）</h3>
    <p style="font-size:.9em">{g1_samples_str}</p>
  </div>
  <div class="group-box group-malignant">
    <h3 style="color:#b91c1c;margin-top:0">🔴 {g2_label}（{len(malignant_samples)} 例）</h3>
    <p style="font-size:.9em">{g2_samples_str}</p>
  </div>
</div>
<p style="color:#666;font-size:.83em;margin-top:10px">分组依据：data/卵巢癌fq编号与样本信息对应表.xlsx</p>
</div>

<!-- 5. On-target 与甲基化提取 -->
<div class="section">
<h2>5. 甲基化提取（MethylDackel extract）</h2>
<table>
  <tr><th>样本</th><th>分组</th><th>检出 CpG 位点数</th><th>平均甲基化水平</th></tr>
  {cpg_rows}
</table>
<table style="margin-top:16px">
  <tr><th>指标</th><th>{g1_label}（合并）</th><th>{g2_label}（合并）</th></tr>
  <tr><td>检出 CpG 位点数（全部）</td><td>{_fmt_num(n_cpg1)}</td><td>{_fmt_num(n_cpg2)}</td></tr>
  <tr><td>平均甲基化水平</td><td>{avg_meth1:.1f}%</td><td>{avg_meth2:.1f}%</td></tr>
  <tr><td>候选区域 CpG 位点数</td><td>{_fmt_num(n_cand1)}</td><td>{_fmt_num(n_cand2)}</td></tr>
  <tr><td>候选区域平均甲基化</td><td>{avg_cand1:.1f}%</td><td>{avg_cand2:.1f}%</td></tr>
  <tr><td>候选区域 100% 甲基化位点</td>
      <td>{full1} ({_safe_pct(full1, n_cand1)})</td>
      <td>{full2} ({_safe_pct(full2, n_cand2)})</td></tr>
</table>
<p style="color:#666;font-size:.83em">候选区域：OC_tissue_88_regions.hg38（88 个卵巢癌差异区域）；提取参数：--minDepth 1 --mergeContext</p>
</div>

<!-- 6. 可视化分布图 -->
<div class="section">
<h2>6. 甲基化分布可视化（候选区域合并）</h2>
<div class="two-col">
  <div><h3>密度分布图</h3>{_img_tag(density_b64, '密度分布图')}</div>
  <div><h3>位置散点图</h3>{_img_tag(scatter_b64, '位置散点图')}</div>
</div>
<div style="margin-top:14px"><h3>甲基化热图</h3>{_img_tag(heatmap_b64, '甲基化热图')}</div>
</div>

<!-- 7. 差异甲基化分析 -->
<div class="section">
<h2>7. 差异甲基化分析（DMC，|Δβ|≥0.2）</h2>
<div class="two-col">
  <div>
    <table>
      <tr><th>统计量</th><th>数值</th></tr>
      <tr><td>差异甲基化位点总数</td><td><strong>{n_dmc}</strong></td></tr>
      <tr><td>{g1_label} 甲基化更高</td><td>{s1_higher} 个位点</td></tr>
      <tr><td>{g2_label} 甲基化更高</td><td>{s2_higher} 个位点</td></tr>
      <tr><td>平均 |Δβ|</td><td>{mean_diff:.3f}</td></tr>
      <tr><td>最大 |Δβ|</td><td>{max_diff:.3f}</td></tr>
    </table>
  </div>
  <div>
    <table>
      <tr><th>染色体</th><th>DMC 数</th><th>占比</th></tr>
      {chrom_rows}
    </table>
  </div>
</div>
<h3 style="margin-top:20px">|Δβ| 最大前 20 个差异位点</h3>
<table>
  <tr><th>染色体</th><th>位置</th><th>{g1_label} β</th><th>{g2_label} β</th><th>|Δβ|</th></tr>
  {top_diff_rows}
</table>
<p style="color:#666;font-size:.83em">红色 = {g1_label} 更高；蓝色 = {g2_label} 更高</p>
</div>

<!-- 8. 结果解读 -->
<div class="section">
<h2>8. 结果解读与质量评价</h2>
<div class="interp">
{interp_html}
</div>
</div>

<!-- 9. 分析流程 -->
<div class="section">
<h2>9. 分析流程</h2>
<div class="funnel">
  <div class="funnel-step">fastp QC<br/>UMI 提取</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">bwameth<br/>比对</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">UMI<br/>去重</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">On-target<br/>过滤</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">M-bias<br/>检查</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">MethylDackel<br/>提取</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">样本<br/>分组</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">候选区域<br/>过滤</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">差异<br/>分析</div><div class="funnel-arrow">→</div>
  <div class="funnel-step">可视化<br/>报告</div>
</div>
<table style="margin-top:14px">
  <tr><th>步骤</th><th>工具</th><th>关键参数</th></tr>
  <tr><td>预处理</td><td>fastp</td><td>trim_front=12/12, trim_tail2=5, min_length=35, Q20, UMI R2 5' 6N</td></tr>
  <tr><td>比对</td><td>bwameth + bwa mem</td><td>-t 8，两步法（view + sort -m 8G），参考 hg38</td></tr>
  <tr><td>去重</td><td>umi_tools dedup</td><td>--paired，UMI 位于 R2 5' 6N</td></tr>
  <tr><td>靶向过滤</td><td>bedtools intersect</td><td>panel_B10308.bed</td></tr>
  <tr><td>M-bias 检查</td><td>MethylDackel mbias</td><td>确认 trim_tail2 参数</td></tr>
  <tr><td>甲基化提取</td><td>MethylDackel extract</td><td>--minDepth 1 --mergeContext → CpG bedGraph</td></tr>
  <tr><td>分组合并</td><td>自研脚本</td><td>按 Excel 表分为良性/恶性，各组 bedGraph 合并</td></tr>
  <tr><td>过滤</td><td>自研脚本 + bedtools</td><td>重复序列过滤、深度过滤、88 候选区域聚焦</td></tr>
  <tr><td>差异分析</td><td>自研脚本</td><td>阈值 |Δβ|≥0.2，覆盖深度统计</td></tr>
</table>
</div>

<p class="meta">工作目录: {temp_dir} &nbsp;|&nbsp; 结果目录: {result_dir}<br/>报告生成时间: {analysis_date}</p>
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
    """从甲基化文件路径提取样本名。"""
    if not methylation_path:
        return "unknown"
    base = os.path.basename(methylation_path)
    base = re.sub(r"(candidate_filtered_|filtered_)", "", base)
    base = re.sub(r"_CpG\.bedGraph$", "", base)
    return base
