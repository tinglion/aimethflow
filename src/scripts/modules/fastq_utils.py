#!/usr/bin/env python3
"""
fastq_utils.py - FastQ文件处理模块

该模块包含FastQ文件质量分析和预处理相关的函数
"""

import os
import logging
import subprocess
import random
import gzip
from Bio import SeqIO

logger = logging.getLogger(__name__)


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


def preprocess_fastq(fastq_files, temp_dir):
    """
    数据预处理（质量控制和接头去除）

    Args:
        fastq_files: FastQ文件列表
        temp_dir: 临时目录

    Returns:
        list: 预处理后的FastQ文件列表
    """
    logger.info("开始数据预处理")
    processed_files = []

    for fastq_file in fastq_files:
        logger.info(f"处理文件: {fastq_file}")

        # 检查是否是双端测序文件
        is_paired = False
        if "_R1" in fastq_file and os.path.exists(fastq_file.replace("_R1", "_R2")):
            is_paired = True
            r1_file = fastq_file
            r2_file = fastq_file.replace("_R1", "_R2")
            output_prefix = os.path.join(temp_dir, os.path.basename(fastq_file).replace("_R1", ""))

            # 运行TrimGalore进行双端处理
            cmd = [
                "trim_galore",
                "--paired",
                "--quality",
                "20",
                "--phred33",
                "--illumina",
                "--length",
                "35",
                "--clip_R1",
                "5",
                "--clip_R2",
                "5",
                r1_file,
                r2_file,
                "--output_dir",
                temp_dir,
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                processed_r1 = os.path.join(
                    temp_dir,
                    os.path.basename(r1_file)
                    .replace(".fq.gz", "_val_1.fq.gz")
                    .replace(".fastq.gz", "_val_1.fq.gz"),
                )
                processed_r2 = os.path.join(
                    temp_dir,
                    os.path.basename(r2_file)
                    .replace(".fq.gz", "_val_2.fq.gz")
                    .replace(".fastq.gz", "_val_2.fq.gz"),
                )
                if os.path.exists(processed_r1) and os.path.exists(processed_r2):
                    processed_files.extend([processed_r1, processed_r2])
                    logger.info(f"双端文件处理完成: {processed_r1}, {processed_r2}")
                else:
                    logger.error(f"双端文件处理失败: {fastq_file}")
            except subprocess.CalledProcessError as e:
                logger.error(f"运行TrimGalore失败: {e.stderr}")
        else:
            # 单端处理
            output_prefix = os.path.join(temp_dir, os.path.basename(fastq_file))
            cmd = [
                "trim_galore",
                "--quality",
                "20",
                "--phred33",
                "--illumina",
                "--length",
                "35",
                "--clip_R1",
                "5",
                fastq_file,
                "--output_dir",
                temp_dir,
            ]

            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                processed_file = os.path.join(
                    temp_dir,
                    os.path.basename(fastq_file)
                    .replace(".fq.gz", "_trimmed.fq.gz")
                    .replace(".fastq.gz", "_trimmed.fq.gz"),
                )
                if os.path.exists(processed_file):
                    processed_files.append(processed_file)
                    logger.info(f"单端文件处理完成: {processed_file}")
                else:
                    logger.error(f"单端文件处理失败: {fastq_file}")
            except subprocess.CalledProcessError as e:
                logger.error(f"运行TrimGalore失败: {e.stderr}")

    logger.info(f"数据预处理完成，处理了 {len(processed_files)} 个文件")
    return processed_files
