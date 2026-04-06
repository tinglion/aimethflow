#!/usr/bin/env python3
"""
fastq_utils.py - FastQ文件处理模块

该模块包含FastQ文件质量分析和预处理相关的函数
"""

import os
import subprocess
import random
import gzip
from Bio import SeqIO
from utils.logger import logger


def analyze_quality(fastq_files):
    """
    分析FastQ文件质量

    Args:
        fastq_files: FastQ文件列表
    """
    logger.info("开始分析FastQ文件质量")
    from fq_inspect import inspect_fastq
    
    for fastq_file in fastq_files:
        logger.info(f"分析文件: {fastq_file}")
        try:
            inspect_fastq(fastq_file)
            logger.info(f"成功分析文件: {fastq_file}")
        except Exception as e:
            error_msg = f"分析文件 {fastq_file} 时出错: {str(e)}"
            logger.error(error_msg)


def get_fastq_stats(fastq_file, sample_size=None, random_seed=None, config=None):
    """
    获取FastQ文件的基本统计信息（使用随机抽样）

    Args:
        fastq_file: FastQ文件路径
        sample_size: 抽样大小，默认为None，使用配置文件中的值
        random_seed: 随机种子，默认为None，使用配置文件中的值
        config: 配置对象

    Returns:
        dict: 包含统计信息的字典
    """

    # 使用配置文件中的值或默认值
    if sample_size is None and config:
        sample_size = config.get("fastq", {}).get("sample_size", 1000)
    else:
        sample_size = 1000
    
    if random_seed is None and config:
        random_seed = config.get("fastq", {}).get("random_seed", 42)
    else:
        random_seed = 42

    stats = {"total_records": 0, "total_bases": 0, "avg_length": 0, "avg_quality": 0}

    try:
        # 设置随机种子
        random.seed(random_seed)

        # 快速估计文件中的总序列数
        # 对于gzip文件，我们读取前100000字节来估计序列密度
        total_sequences = 0
        try:
            # 获取文件大小
            file_size = os.path.getsize(fastq_file)

            # 读取文件的前一部分来估计序列密度
            with gzip.open(fastq_file, "rt") as handle:
                # 读取前100000字节
                sample_data = handle.read(100000)
                # 计算样本中的序列数（以@开头的行）
                sample_sequences = sample_data.count("\n@")

                # 估计总序列数
                if sample_sequences > 0:
                    # 计算样本数据的大小（字节）
                    # 估计总序列数
                    total_sequences = int(
                        (file_size / len(sample_data.encode("utf-8"))) * sample_sequences
                    )
                else:
                    # 如果样本中没有序列，使用默认值
                    total_sequences = 1000000
        except Exception as e:
            logger.warning(f"无法估计序列数，使用默认值: {str(e)}")
            total_sequences = 1000000

        # 计算抽样概率
        if total_sequences <= sample_size:
            # 如果文件序列数小于等于抽样大小，全部读取
            sample_prob = 1.0
        else:
            sample_prob = sample_size / total_sequences

        logger.info(f"sample {fastq_file} p={sample_prob:.8f}={sample_size}/{total_sequences}")
        with gzip.open(fastq_file, "rt") as handle:
            lengths = []
            quality_scores = []

            for record in SeqIO.parse(handle, "fastq"):
                # 随机抽样
                if random.random() <= sample_prob:
                    stats["total_records"] += 1
                    length = len(record)
                    stats["total_bases"] += length
                    lengths.append(length)

                    if record.letter_annotations.get("phred_quality"):
                        avg_quality = sum(record.letter_annotations["phred_quality"]) / length
                        quality_scores.append(avg_quality)

                # 如果已经达到抽样大小，停止
                if stats["total_records"] >= sample_size:
                    break

            if lengths:
                stats["avg_length"] = stats["total_bases"] / stats["total_records"]

            if quality_scores:
                stats["avg_quality"] = sum(quality_scores) / len(quality_scores)

        # 添加总序列数到统计信息中
        stats["original_total_records"] = total_sequences

    except Exception as e:
        logger.error(f"获取FastQ统计信息时出错: {str(e)}")

    return stats


def find_r2_for_r1(r1_file):
    """
    根据R1文件路径找到对应的R2文件，兼容 .R1. 和 _R1 两种命名方式。

    Returns:
        str or None: R2文件路径，如果未找到返回None
    """
    for r1_pat, r2_pat in [(".R1.", ".R2."), ("_R1.", "_R2."), ("_R1_", "_R2_"),
                            (".r1.", ".r2."), ("_r1.", "_r2.")]:
        if r1_pat in r1_file:
            r2_file = r1_file.replace(r1_pat, r2_pat)
            if os.path.exists(r2_file):
                return r2_file
    return None


def is_r1_file(fastq_file):
    """判断该文件是否是R1文件（避免重复处理R2）。"""
    for pat in [".R1.", "_R1.", "_R1_", ".r1.", "_r1."]:
        if pat in fastq_file:
            return True
    return False


def preprocess_fastq(fastq_files, temp_dir, plasma=False, config=None):
    """
    数据预处理（质量控制和接头去除）。
    使用 fastp 替代 TrimGalore，支持血浆模式的 UMI 提取。
    幂等：若输出 .clean.R1.fq.gz 已存在则跳过该样本。

    Args:
        fastq_files: FastQ文件列表
        temp_dir: 临时目录
        plasma: 是否启用血浆模式（UMI提取 + 5'端低复杂度切除 + 3'端poly-C切除）
        config: 配置对象

    Returns:
        list: 预处理后的FastQ文件列表（R1在前，R2在后，配对顺序）
    """
    logger.info(f"开始数据预处理 ({'血浆模式' if plasma else '标准模式'})")
    processed_files = []

    # 读取血浆配置
    plasma_cfg = (config or {}).get("plasma", {})
    trim_front1 = plasma_cfg.get("trim_front1", 12)
    trim_front2 = plasma_cfg.get("trim_front2", 12)
    trim_tail2  = plasma_cfg.get("trim_tail2", 5)
    min_length  = plasma_cfg.get("min_length", 35)
    min_quality = plasma_cfg.get("min_quality", 20)
    umi_length  = plasma_cfg.get("umi_length", 6)

    # 线程数从 config [threads] 读取
    n_threads = str((config or {}).get("threads", {}).get("fastp", 8))

    # 找出所有R1文件（避免重复处理R2）
    r1_files = [f for f in fastq_files if is_r1_file(f)]
    # 补充没有R1/R2标记的文件（单端）
    unpaired_files = [f for f in fastq_files if not is_r1_file(f)
                      and not any(p in f for p in [".R2.", "_R2.", "_R2_", ".r2.", "_r2."])]

    for r1_file in r1_files:
        r2_file = find_r2_for_r1(r1_file)
        sample_name = os.path.basename(r1_file)
        # 移除所有常见后缀得到样本名
        for pat in [".R1.fq.gz", ".R1.fastq.gz", "_R1.fq.gz", "_R1.fastq.gz",
                    ".r1.fq.gz", ".r1.fastq.gz"]:
            if sample_name.endswith(pat):
                sample_name = sample_name[:-len(pat)]
                break
        else:
            sample_name = sample_name.split(".")[0]

        out_r1 = os.path.join(temp_dir, f"{sample_name}.clean.R1.fq.gz")
        out_r2 = os.path.join(temp_dir, f"{sample_name}.clean.R2.fq.gz")
        json_report = os.path.join(temp_dir, f"{sample_name}.fastp.json")
        html_report = os.path.join(temp_dir, f"{sample_name}.fastp.html")

        # ── 幂等检查：输出文件已存在则跳过 ──
        if r2_file and os.path.exists(out_r1) and os.path.exists(out_r2):
            logger.info(f"fastp输出已存在，跳过: {sample_name}")
            processed_files.extend([out_r1, out_r2])
            continue
        if not r2_file and os.path.exists(out_r1):
            logger.info(f"fastp输出已存在，跳过: {sample_name}")
            processed_files.append(out_r1)
            continue

        if r2_file:
            cmd = [
                "fastp",
                "-i", r1_file, "-I", r2_file,
                "-o", out_r1, "-O", out_r2,
                "-q", str(min_quality),
                "--length_required", str(min_length),
                "--json", json_report,
                "--html", html_report,
                "--thread", n_threads,
            ]
            logger.info(f"正在处理双端样本({n_threads}): {sample_name} (R1: {r1_file}, R2: {r2_file})")
            if plasma:
                cmd += [
                    "--umi", "--umi_loc", "read2",
                    "--umi_len", str(umi_length),
                    "--umi_prefix", "UMI",
                    "--trim_front1", str(trim_front1),
                    "--trim_front2", str(trim_front2),
                    "--trim_tail2", str(trim_tail2),
                ]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"fastp双端处理完成: {sample_name}")
                if os.path.exists(out_r1) and os.path.exists(out_r2):
                    processed_files.extend([out_r1, out_r2])
                else:
                    logger.error(f"fastp未生成预期输出文件: {sample_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"fastp处理失败 ({sample_name}): {e.stderr}")
        else:
            # 单端
            cmd = [
                "fastp",
                "-i", r1_file,
                "-o", out_r1,
                "-q", str(min_quality),
                "--length_required", str(min_length),
                "--json", json_report,
                "--html", html_report,
                "--thread", n_threads,
            ]
            if plasma:
                cmd += ["--trim_front1", str(trim_front1)]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"fastp单端处理完成: {sample_name}")
                if os.path.exists(out_r1):
                    processed_files.append(out_r1)
                else:
                    logger.error(f"fastp未生成预期输出文件: {sample_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"fastp处理失败 ({sample_name}): {e.stderr}")

    for fastq_file in unpaired_files:
        sample_name = os.path.basename(fastq_file).split(".")[0]
        out_r1 = os.path.join(temp_dir, f"{sample_name}.clean.fq.gz")
        json_report = os.path.join(temp_dir, f"{sample_name}.fastp.json")
        html_report = os.path.join(temp_dir, f"{sample_name}.fastp.html")
        if os.path.exists(out_r1):
            logger.info(f"fastp输出已存在，跳过: {sample_name}")
            processed_files.append(out_r1)
            continue
        cmd = [
            "fastp", "-i", fastq_file, "-o", out_r1,
            "-q", str(min_quality), "--length_required", str(min_length),
            "--json", json_report, "--html", html_report, "--thread", n_threads,
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if os.path.exists(out_r1):
                processed_files.append(out_r1)
        except subprocess.CalledProcessError as e:
            logger.error(f"fastp处理失败 ({fastq_file}): {e.stderr}")

    logger.info(f"数据预处理完成，处理了 {len(processed_files)} 个文件")
    return processed_files
