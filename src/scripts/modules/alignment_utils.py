#!/usr/bin/env python3
"""
alignment_utils.py - 序列比对模块

支持两种比对工具（通过 config.toml [plasma] aligner 切换）：
  bismark  — 精准，慢（3-4h/30M pairs）。首次验证、需要OT/OB方向性报告时用。
  bwameth  — 快3-5倍（40-60min/30M pairs）。100+样本批量分析推荐。

甲基化提取：
  bismark  → bismark_methylation_extractor（.bismark.cov）
  bwameth  → MethylDackel extract（_CpG.bedGraph）
"""

import os
import shutil
import subprocess
import traceback

from utils.logger import logger


def _rm(path):
    """安全删除单个文件，失败时只记录 warning 不抛出。"""
    try:
        if os.path.exists(path):
            os.remove(path)
            logger.debug(f"  [清理] 删除中间文件: {path}")
    except OSError as e:
        logger.warning(f"  [清理] 删除失败: {path} — {e}")


def _resolve_tool(name, config):
    """
    解析工具可执行文件的绝对路径。
    优先顺序：
      1. config [tools] conda_env_bin/{name}（直接给出绝对路径，最可靠）
      2. shutil.which(name)（已在 PATH 中）
      3. 原样返回 name（让 subprocess 自己报错）
    """
    conda_bin = (config or {}).get("tools", {}).get("conda_env_bin", "")
    if conda_bin:
        candidate = os.path.join(conda_bin, name)
        if os.path.exists(candidate):
            return candidate
    found = shutil.which(name)
    if found:
        return found
    return name


# ──────────────────────────────────────────────
# 内部工具：文件名解析
# ──────────────────────────────────────────────

def _find_paired_fastq(fastq_files):
    """
    从文件列表中找出R1/R2配对，兼容多种命名约定。
    优先匹配 .clean.R1/.clean.R2（fastp输出），也兼容 _val_1/_val_2（TrimGalore输出）。

    Returns:
        list of (r1, r2_or_None)
    """
    pairs = []
    r2_seen = set()
    for f in fastq_files:
        if f in r2_seen:
            continue
        r2 = None
        for r1_pat, r2_pat in [
            (".clean.R1.", ".clean.R2."),
            (".R1.", ".R2."),
            ("_R1.", "_R2."),
            ("_R1_", "_R2_"),
            ("_val_1", "_val_2"),
        ]:
            if r1_pat in f:
                candidate = f.replace(r1_pat, r2_pat)
                if candidate in fastq_files:
                    r2 = candidate
                    r2_seen.add(candidate)
                break
        pairs.append((f, r2))
    return pairs


def _sample_base_from_r1(r1_file):
    """从R1文件名提取样本基础名（去除所有已知后缀）。"""
    base = os.path.basename(r1_file)
    for pat in [".clean.R1.fq.gz", ".R1.fq.gz", "_R1.fq.gz",
                ".clean.R1.fastq.gz", ".R1.fastq.gz", "_val_1.fq.gz", "_trimmed.fq.gz"]:
        if base.endswith(pat):
            return base[:-len(pat)]
    return base.split(".")[0]


# ──────────────────────────────────────────────
# bismark 比对（精准模式）
# ──────────────────────────────────────────────

def _align_bismark_pair(r1_file, r2_file, sample_base, temp_dir, config, non_directional=False):
    """
    bismark --bowtie2 双端比对。
    速度：~3-4h / 30M pairs（8核）。
    优势：输出 OT/OB/CTOT/CTOB 链比例报告，用于验证文库方向性。
    """
    bismark_index = os.path.dirname(config.get("extraction", {}).get("reference_genome", ""))
    cmd = [
        _resolve_tool("bismark", config), "--bowtie2", "--bam", "--multicore", "4",
        bismark_index, "-1", r1_file, "-2", r2_file,
        "--output_dir", temp_dir,
    ]
    if non_directional:
        cmd.append("--non_directional")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        # 在 temp_dir 找 bismark 输出的 BAM
        for fname in sorted(os.listdir(temp_dir)):
            if sample_base in fname and fname.endswith(".bam") and "bismark_bt2" in fname and "sorted" not in fname and "temp" not in fname:
                bam_path = os.path.join(temp_dir, fname)
                sorted_bam = bam_path.replace(".bam", ".sorted.bam")
                subprocess.run([_resolve_tool("samtools", config), "sort", "-o", sorted_bam, bam_path],
                                check=True, capture_output=True)
                subprocess.run([_resolve_tool("samtools", config), "index", sorted_bam],
                                check=True, capture_output=True)
                logger.info(f"bismark比对+排序完成: {sorted_bam}")
                _rm(bam_path)  # 删除未排序中间 BAM
                # 打印方向性报告路径（供用户确认 OT/OB 比例）
                for fname2 in os.listdir(temp_dir):
                    if sample_base in fname2 and "PE_report" in fname2:
                        logger.info(f"  ↳ 方向性报告（请检查OT/OB比例）: {os.path.join(temp_dir, fname2)}")
                return sorted_bam
        logger.error(f"bismark未生成预期BAM文件: {sample_base}")
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"bismark比对失败 ({r1_file}): {e.stderr[-500:]}")
        return None


# ──────────────────────────────────────────────
# bwameth 比对（快速模式，推荐批量）
# ──────────────────────────────────────────────

def _align_bwameth_pair(r1_file, r2_file, sample_base, temp_dir, config, non_directional=False):  # noqa: ARG001
    """
    bwameth 双端比对（BWA-MEM + C→T转换）。
    速度：~40-60min / 30M pairs（8核），比 bismark 快 3-5 倍。
    注意：不输出 OT/OB 链比例报告，建议先用 bismark 验证方向性后再切换。

    索引要求：运行前须先执行：
        bwameth.py index /data/meth/bwa/hg38.fa
    """
    reference_genome = config.get("extraction", {}).get("reference_genome", "")
    sorted_bam = os.path.join(temp_dir, f"{sample_base}.bwameth_pe.sorted.bam")

    # 线程数从 config [threads] 读取
    threads_cfg = config.get("threads", {})
    bwameth_t   = str(threads_cfg.get("bwameth", 8))
    sort_t      = str(threads_cfg.get("samtools_sort", 4))
    sort_mem    = str(threads_cfg.get("samtools_sort_mem", "8G"))

    # 幂等检查：已存在则跳过
    if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
        logger.info(f"bwameth BAM已存在，跳过: {sorted_bam}")
        return sorted_bam

    cmd_align   = ["bwameth.py", "--reference", reference_genome, "-t", bwameth_t,
                    r1_file, r2_file]
    unsorted_bam = sorted_bam.replace(".sorted.bam", ".unsorted.bam")

    align_log = sorted_bam.replace(".sorted.bam", ".align.log")
    try:
        # 步骤 1：bwameth → samtools view（流式 SAM→BAM 转换）
        # stderr 写到文件，避免 PIPE 缓冲满后 bwameth 阻塞 stdout，造成死锁
        cmd_view   = ["samtools", "view", "-bS", "-o", unsorted_bam, "-"]
        logger.info(f"  [1/2] bwameth 比对中（输出至未排序 BAM）")
        with open(align_log, "w") as f_log:
            logger.info(f"    bwameth 命令: {' '.join(cmd_align)}")
            proc_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=f_log)
            logger.info(f"    samtools view 命令: {' '.join(cmd_view)}")
            proc_view  = subprocess.Popen(cmd_view, stdin=proc_align.stdout,
                                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            proc_align.stdout.close()
            proc_view.wait()
            proc_align.wait()

        if proc_view.returncode != 0:
            raise RuntimeError(f"samtools view 失败，退出码: {proc_view.returncode}")
        if proc_align.returncode != 0:
            raise RuntimeError(f"bwameth 失败，退出码: {proc_align.returncode}，详见: {align_log}")

        # 步骤 2：从文件排序（无管道，bwameth 已退出，彻底消除死锁）
        logger.info(f"  [2/2] samtools sort 排序中...")
        subprocess.run(
            ["samtools", "sort", "-@", sort_t, "-m", sort_mem, "-o", sorted_bam, unsorted_bam],
            check=True, capture_output=True
        )
        _rm(unsorted_bam)

        subprocess.run(["samtools", "index", sorted_bam], check=True, capture_output=True)
        _rm(align_log)
        logger.info(f"bwameth比对+排序完成: {sorted_bam}")
        return sorted_bam
    except Exception as e:
        logger.error(f"bwameth比对失败 ({r1_file}): {e}")
        traceback.print_exc()
        _rm(unsorted_bam)
        return None


# ──────────────────────────────────────────────
# 公共接口：align_methylation
# ──────────────────────────────────────────────

def align_methylation(fastq_files, temp_dir, config, plasma=False):
    """
    甲基化序列比对。

    比对工具由 config.toml [plasma] aligner 决定：
      aligner = "bismark"  — 慢但给出方向性报告
      aligner = "bwameth"  — 快3-5倍，批量推荐（需先建 c2t 索引）

    Args:
        fastq_files: 预处理后的 FastQ 文件列表
        temp_dir:    临时目录
        config:      配置对象
        plasma:      是否为血浆模式

    Returns:
        list: 坐标排序 + 建索引的 BAM 文件路径列表
    """
    aligner = config.get("plasma", {}).get("aligner", "bismark") if plasma else "bismark"
    non_directional = config.get("plasma", {}).get("non_directional", False)
    logger.info(f"开始甲基化序列比对 [工具: {aligner}, "
                f"{'血浆' if plasma else '标准'}模式"
                f"{', non-directional' if plasma and non_directional else ''}]")

    bam_files = []
    paired_files = _find_paired_fastq(fastq_files)

    for r1_file, r2_file in paired_files:
        logger.info(f"处理文件: {r1_file}")
        if r2_file is None:
            logger.warning(f"未找到R2配对文件，跳过单端文件: {r1_file}")
            continue
        sample_base = _sample_base_from_r1(r1_file)
        logger.info(f"  比对样本: {sample_base}  ({r1_file})")

        if aligner == "bwameth":
            bam = _align_bwameth_pair(r1_file, r2_file, sample_base,
                                      temp_dir, config, non_directional)
        else:
            bam = _align_bismark_pair(r1_file, r2_file, sample_base,
                                      temp_dir, config, non_directional)
        if bam:
            bam_files.append(bam)
        logger.info(f"finished {r1_file}")

    logger.info(f"比对完成，生成 {len(bam_files)} 个BAM文件")
    return bam_files


# ──────────────────────────────────────────────
# 去重
# ──────────────────────────────────────────────

def deduplicate_bam(bam_files, temp_dir, plasma=False, config=None):
    """
    PCR 重复去除。
    血浆模式：umi_tools dedup（UMI感知，避免误删cfDNA真实片段）。
    标准模式：deduplicate_bismark（坐标去重）。
    """
    logger.info(f"开始PCR重复去除 ({'umi_tools血浆模式' if plasma else 'bismark标准模式'})")
    delete_prev = (config or {}).get("cleanup", {}).get("delete_after_dedup", False)
    deduplicated_files = []

    for bam_file in bam_files:
        logger.info(f"处理文件: {bam_file}")
        if plasma:
            dedup_bam = bam_file.replace(".sorted.bam", ".dedup.bam")
            if dedup_bam == bam_file:
                dedup_bam = bam_file.replace(".bam", ".dedup.bam")
            stats_file = dedup_bam.replace(".bam", ".stats")
            cmd = ["umi_tools", "dedup", "--paired",
                   "--stdin", bam_file, "--stdout", dedup_bam, "--log", stats_file]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                subprocess.run(["samtools", "index", dedup_bam],
                                check=True, capture_output=True)
                deduplicated_files.append(dedup_bam)
                logger.info(f"umi_tools去重完成: {dedup_bam}")
                if delete_prev:
                    _rm(bam_file)
                    _rm(bam_file + ".bai")
            except subprocess.CalledProcessError as e:
                logger.error(f"umi_tools dedup失败: {e.stderr[-300:]}")
                logger.warning("回退到 deduplicate_bismark（不推荐用于血浆数据）")
                deduplicated_files.extend(_deduplicate_bismark([bam_file], temp_dir))
        else:
            deduplicated_files.extend(_deduplicate_bismark([bam_file], temp_dir))

    logger.info(f"去重完成，生成 {len(deduplicated_files)} 个去重BAM文件")
    return deduplicated_files


def _deduplicate_bismark(bam_files, temp_dir):
    """内部函数：deduplicate_bismark 坐标去重（组织模式）。"""
    result = []
    for bam_file in bam_files:
        cmd = ["deduplicate_bismark", "--bam", bam_file, "--output_dir", temp_dir]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            dedup_file = bam_file.replace(".bam", ".deduplicated.bam")
            if os.path.exists(dedup_file):
                result.append(dedup_file)
                logger.info(f"deduplicate_bismark去重完成: {dedup_file}")
            else:
                logger.error(f"deduplicate_bismark未生成预期文件: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"deduplicate_bismark失败: {e.stderr[-300:]}")
    return result


# ──────────────────────────────────────────────
# On-target 过滤
# ──────────────────────────────────────────────

def extract_target_regions(bam_files, temp_dir, config):
    """bedtools intersect：只保留 panel.bed 覆盖区域的 reads。"""
    logger.info("开始靶向区域提取")
    target_files = []
    delete_prev = (config or {}).get("cleanup", {}).get("delete_after_target", False)

    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if not probe_panel or not os.path.exists(probe_panel):
        logger.warning("探针区域文件不存在，跳过靶向区域提取")
        return bam_files

    for bam_file in bam_files:
        logger.info(f"处理文件: {bam_file}")
        output_file = os.path.join(
            temp_dir, os.path.basename(bam_file).replace(".bam", ".targeted.bam")
        )
        cmd = ["bedtools", "intersect", "-abam", bam_file,
               "-b", probe_panel, "-wa", "-header"]
        try:
            with open(output_file, "wb") as fout:
                subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.PIPE)
            sorted_file = output_file.replace(".bam", ".sorted.bam")
            subprocess.run(["samtools", "sort", output_file, "-o", sorted_file],
                            check=True, capture_output=True)
            subprocess.run(["samtools", "index", sorted_file],
                            check=True, capture_output=True)
            if os.path.exists(sorted_file + ".bai"):
                target_files.append(sorted_file)
                logger.info(f"靶向区域提取完成: {sorted_file}")
                _rm(output_file)       # 删除未排序中间 targeted.bam
                if delete_prev:
                    _rm(bam_file)
                    _rm(bam_file + ".bai")
            else:
                logger.error(f"靶向区域提取失败: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"bedtools intersect失败: {e.stderr}")

    logger.info(f"靶向区域提取完成，生成 {len(target_files)} 个靶向BAM文件")
    return target_files


# ──────────────────────────────────────────────
# M-bias 检查
# ──────────────────────────────────────────────

def check_mbias(bam_files, temp_dir, reference_genome, config=None):
    """
    MethylDackel mbias：检查 reads 各位置甲基化率是否存在系统偏差。
    生成 SVG 图像，供分析人员判断是否需要调整 config.toml trim_tail2 参数。
    """
    logger.info("开始M-bias检查")
    n_threads = str((config or {}).get("threads", {}).get("methyldackel", 4))
    report_files = []
    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split(".")[0]
        output_prefix = os.path.join(temp_dir, f"{sample_name}.mbias")
        cmd = ["MethylDackel", "mbias", "--nOT", "0,0,0,0", "-@", n_threads,
               reference_genome, bam_file, output_prefix]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            svg_file = f"{output_prefix}_OT.svg"
            if os.path.exists(svg_file):
                report_files.append(svg_file)
            logger.info(
                f"M-bias检查完成: {sample_name}\n"
                f"  → SVG 文件: {output_prefix}_OT.svg / _OB.svg\n"
                f"  → 若R2 3'端甲基化率断崖下降，调整 config.toml [plasma] trim_tail2 后重跑"
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"MethylDackel mbias失败 ({sample_name}): {e.stderr[-300:]}")
    return report_files


# ──────────────────────────────────────────────
# 甲基化提取：MethylDackel（配合 bwameth 使用）
# ──────────────────────────────────────────────

def extract_methylation_methyldackel(bam_files, temp_dir, config):
    """
    使用 MethylDackel extract 提取 CpG 甲基化信息。
    配合 bwameth 比对结果使用；输出为 bedGraph 格式。

    输出格式（每行一个CpG位点）：
        chrom  start  end  methylation%  meth_reads  total_reads

    Args:
        bam_files: on-target BAM 文件列表（已排序建索引）
        temp_dir:  临时目录
        config:    配置对象

    Returns:
        list: bedGraph 文件路径列表（{sample}_CpG.bedGraph）
    """
    logger.info("开始甲基化信息提取 [MethylDackel extract]")
    reference_genome = config.get("extraction", {}).get("reference_genome", "")
    if not reference_genome or not os.path.exists(reference_genome):
        logger.error(f"参考基因组不存在: {reference_genome}")
        return []

    result_files = []
    n_threads = str(config.get("threads", {}).get("methyldackel", 4))
    plasma_cfg = config.get("plasma", {})
    nOT = plasma_cfg.get("nOT", "0,0,0,0")
    nOB = plasma_cfg.get("nOB", "0,0,0,0")
    logger.info(f"  M-bias校正参数: --nOT {nOT}  --nOB {nOB}")
    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split(".")[0]
        output_prefix = os.path.join(temp_dir, sample_name)
        output_file = f"{output_prefix}_CpG.bedGraph"

        # 幂等检查：已提取则跳过
        if os.path.exists(output_file):
            logger.info(f"MethylDackel bedGraph已存在，跳过: {output_file}")
            result_files.append(output_file)
            continue

        cmd = [
            _resolve_tool("MethylDackel", config), "extract",
            "--minDepth", "1",
            "--mergeContext",
            "--nOT", nOT,
            "--nOB", nOB,
            "-@", n_threads,
            reference_genome, bam_file,
            "-o", output_prefix,
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if os.path.exists(output_file):
                result_files.append(output_file)
                logger.info(f"MethylDackel提取完成: {output_file}")
            else:
                logger.error(f"MethylDackel未生成预期输出: {output_prefix}_CpG.bedGraph")
        except subprocess.CalledProcessError as e:
            logger.error(f"MethylDackel extract失败 ({sample_name}): {e.stderr[-300:]}")

    logger.info(f"甲基化提取完成，生成 {len(result_files)} 个bedGraph文件")
    return result_files
