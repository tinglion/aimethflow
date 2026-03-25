#!/usr/bin/env python3
"""
alignment_utils.py - 序列比对模块

该模块包含序列比对和PCR重复去除相关的函数
"""

import os
import logging
import subprocess

logger = logging.getLogger(__name__)


def align_methylation(fastq_files, temp_dir, config):
    """
    甲基化序列比对

    Args:
        fastq_files: 预处理后的FastQ文件列表
        temp_dir: 临时目录
        config: 配置对象

    Returns:
        list: 比对后的BAM文件列表
    """
    logger.info("开始甲基化序列比对")
    bam_files = []

    # 从配置文件中读取参考基因组路径
    reference_genome = config.get("extraction", {}).get("reference_genome")
    bismark_index = os.path.dirname(reference_genome)

    # 检查是否是双端测序文件
    paired_files = []
    i = 0
    while i < len(fastq_files):
        if (
            i + 1 < len(fastq_files)
            and "_val_1" in fastq_files[i]
            and "_val_2" in fastq_files[i + 1]
        ):
            paired_files.append((fastq_files[i], fastq_files[i + 1]))
            i += 2
        else:
            paired_files.append((fastq_files[i], None))
            i += 1

    for r1_file, r2_file in paired_files:
        logger.info(f"比对文件: {r1_file}")
        output_prefix = os.path.join(
            temp_dir,
            os.path.basename(r1_file).replace("_val_1.fq.gz", "").replace("_trimmed.fq.gz", ""),
        )

        if r2_file:
            # 双端比对
            cmd = [
                "bismark",
                "--bowtie2",
                "--bam",
                "--multicore",
                "4",
                bismark_index,
                "-1",
                r1_file,
                "-2",
                r2_file,
                "--output_dir",
                temp_dir,
            ]
        else:
            # 单端比对
            cmd = [
                "bismark",
                "--bowtie2",
                "--bam",
                "--multicore",
                "4",
                bismark_index,
                r1_file,
                "--output_dir",
                temp_dir,
            ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            # 查找生成的BAM文件
            for file in os.listdir(temp_dir):
                if os.path.basename(r1_file).split(".")[0] in file and file.endswith(".bam"):
                    bam_file = os.path.join(temp_dir, file)
                    bam_files.append(bam_file)
                    logger.info(f"比对完成: {bam_file}")
                    break
            else:
                logger.error(f"比对失败，未找到生成的BAM文件: {r1_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"运行Bismark失败: {e.stderr}")

    logger.info(f"甲基化序列比对完成，生成了 {len(bam_files)} 个BAM文件")
    return bam_files


def deduplicate_bam(bam_files, temp_dir):
    """
    PCR重复去除

    Args:
        bam_files: BAM文件列表
        temp_dir: 临时目录

    Returns:
        list: 去重后的BAM文件列表
    """
    logger.info("开始PCR重复去除")
    deduplicated_files = []

    for bam_file in bam_files:
        logger.info(f"处理文件: {bam_file}")

        # 运行deduplicate_bismark
        cmd = ["deduplicate_bismark", "--bam", bam_file, "--output_dir", temp_dir]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            # 查找生成的去重BAM文件
            deduplicated_file = bam_file.replace(".bam", ".deduplicated.bam")
            if os.path.exists(deduplicated_file):
                deduplicated_files.append(deduplicated_file)
                logger.info(f"去重完成: {deduplicated_file}")
            else:
                logger.error(f"去重失败，未找到生成的文件: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"运行deduplicate_bismark失败: {e.stderr}")

    logger.info(f"PCR重复去除完成，生成了 {len(deduplicated_files)} 个去重BAM文件")
    return deduplicated_files


def extract_target_regions(bam_files, temp_dir, config):
    """
    靶向区域提取

    Args:
        bam_files: BAM文件列表
        temp_dir: 临时目录
        config: 配置对象

    Returns:
        list: 提取靶向区域后的BAM文件列表
    """
    logger.info("开始靶向区域提取")
    target_files = []

    # 从配置文件中读取探针区域文件
    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if not probe_panel or not os.path.exists(probe_panel):
        logger.warning("探针区域文件不存在，跳过靶向区域提取")
        return bam_files

    for bam_file in bam_files:
        logger.info(f"处理文件: {bam_file}")
        output_file = os.path.join(
            temp_dir, os.path.basename(bam_file).replace(".bam", ".targeted.bam")
        )

        # 运行bedtools intersect
        cmd = [
            "bedtools",
            "intersect",
            "-a",
            bam_file,
            "-b",
            probe_panel,
            "-wa",
            "-header",
            "-o",
            output_file,
        ]

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            # 排序和建索引
            sorted_file = output_file.replace(".bam", ".sorted.bam")
            index_file = sorted_file + ".bai"

            subprocess.run(["samtools", "sort", output_file, "-o", sorted_file], check=True)
            subprocess.run(["samtools", "index", sorted_file], check=True)

            if os.path.exists(sorted_file) and os.path.exists(index_file):
                target_files.append(sorted_file)
                logger.info(f"靶向区域提取完成: {sorted_file}")
            else:
                logger.error(f"靶向区域提取失败: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"运行bedtools失败: {e.stderr}")

    logger.info(f"靶向区域提取完成，生成了 {len(target_files)} 个靶向BAM文件")
    return target_files
