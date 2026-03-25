#!/usr/bin/env python3
"""
run_analysis.py - 执行完整的甲基化分析流程

该脚本用于执行完整的甲基化分析流程，包括：
1. 创建带时间戳的子目录
2. 分析FastQ文件质量（可选）
3. 提取甲基化信息（可选）
4. 进行差异甲基化分析（可选）
5. 可视化结果（可选）

使用方法:
    python src/scripts/run_analysis.py [选项] <文件/文件夹...>

参数:
    文件/文件夹: 可以是一个或多个FastQ文件或包含FastQ文件的文件夹

选项:
    --skip-quality: 跳过FastQ文件质量分析
    --skip-extraction: 跳过甲基化信息提取
    --skip-diff: 跳过差异甲基化分析
    --skip-visualization: 跳过结果可视化
"""

import os
import sys
import time
import argparse
import glob
import logging
import random

import gzip
from Bio import SeqIO

# 导入各个脚本的函数
from utils.logger import setup_logger
from utils.config import load_config, get_output_dirs, create_output_dirs
from fq_inspect import inspect_fastq
from methylation.methylation_extractor import extract_methylation
from methylation.bam_generator import generate_bam
from methylation.diff_methylation import run_diff_analysis
from methylation.methylation_visualizer import visualize_methylation
from methylation.filter_methylation import filter_methylation
from methylation.candidate_regions import filter_to_candidate_regions
from methylation.sample_grouping import group_methylation_files_v1, group_methylation_files_v2

# 加载配置文件
config = load_config()

# 生成日志文件路径
log_dir = config["directories"]["log_dir"]
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/analysis_{time.strftime('%Y%m%d_%H%M%S')}.log"
logger = setup_logger(log_file=log_file, console_level=logging.INFO, file_level=logging.DEBUG)


def find_fastq_files(paths):
    """
    查找所有FastQ文件

    Args:
        paths: 文件或文件夹路径列表
        logger: 日志记录器

    Returns:
        list: FastQ文件路径列表
    """
    fastq_files = []
    for path in paths:
        if os.path.isfile(path):
            # 如果是文件，检查是否是FastQ文件
            if (
                path.endswith(".fq.gz")
                or path.endswith(".fastq.gz")
                or path.endswith(".fq")
                or path.endswith(".fastq")
            ):
                fastq_files.append(os.path.abspath(path))
        elif os.path.isdir(path):
            # 如果是文件夹，递归查找FastQ文件
            for ext in ["*.fq.gz", "*.fastq.gz", "*.fq", "*.fastq"]:
                fastq_files.extend(glob.glob(os.path.join(path, "**", ext), recursive=True))
    # 去重并排序
    fastq_files = sorted(list(set(fastq_files)))
    return fastq_files


def analyze_quality(fastq_files):
    """
    分析FastQ文件质量

    Args:
        fastq_files: FastQ文件列表
        logger: 日志记录器
    """
    logger.info("开始分析FastQ文件质量")
    for fastq_file in fastq_files:
        logger.info(f"分析文件: {fastq_file}")
        try:
            inspect_fastq(fastq_file)
            logger.info(f"成功分析文件: {fastq_file}")
        except Exception as e:
            error_msg = f"分析文件 {fastq_file} 时出错: {str(e)}"
            logger.error(error_msg)


def get_fastq_stats(fastq_file, sample_size=None, random_seed=None):
    """
    获取FastQ文件的基本统计信息（使用随机抽样）

    Args:
        fastq_file: FastQ文件路径
        sample_size: 抽样大小，默认为None，使用配置文件中的值
        random_seed: 随机种子，默认为None，使用配置文件中的值

    Returns:
        dict: 包含统计信息的字典
    """

    # 使用配置文件中的值或默认值
    if sample_size is None:
        sample_size = config.get("fastq", {}).get("sample_size", 1000)
    if random_seed is None:
        random_seed = config.get("fastq", {}).get("random_seed", 42)

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


def extract_methylation_info(temp_dir, fastq_files, group_method="v2"):
    """
    提取甲基化信息（真实提取）

    Args:
        temp_dir: 临时目录
        fastq_files: FastQ文件列表
        group_method: 分组方法，可选值："v1"（基于索引）或 "v2"（基于Excel文件）

    Returns:
        tuple: (sample1_methylation, sample2_methylation) 文件路径
    """
    logger.info("开始提取甲基化信息")

    # 从配置文件中读取参考基因组路径
    reference_genome = config.get("extraction", {}).get("reference_genome")

    # 为每个文件生成BAM文件（源文件名后缀改为.bam）
    logger.info("开始为每个FastQ文件生成BAM文件")
    all_bam_files = []

    # 从配置文件中读取线程数，默认为4
    threads = config.get("extraction", {}).get("threads", 4)

    # 为所有文件生成BAM文件
    if fastq_files:
        logger.info(f"为 {len(fastq_files)} 个文件生成BAM文件")
        all_bam_files = generate_bam(fastq_files, reference_genome, threads=threads)

    # 提取所有BAM文件的甲基化信息
    all_methylation_files = []
    if all_bam_files:
        logger.info(f"为 {len(all_bam_files)} 个BAM文件提取甲基化信息")
        for i, bam_file in enumerate(all_bam_files):
            # 为每个BAM文件生成单独的甲基化输出文件，使用同文件名
            individual_output = bam_file.replace(".bam", ".txt")
            logger.info(f"为BAM文件 {bam_file} 提取甲基化信息")
            extract_methylation(bam_file, individual_output, reference_genome=reference_genome)
            # 验证提取结果
            validate_methylation_output(individual_output)
            all_methylation_files.append(individual_output)

    # 生成甲基化数据文件
    sample1_methylation = os.path.join(temp_dir, "sample1_methylation.txt")
    sample2_methylation = os.path.join(temp_dir, "sample2_methylation.txt")

    # 根据分组方法选择使用哪个分组函数
    if group_method == "v2":
        logger.info("使用基于Excel文件的分组方法（v2版本）")
        sample1_methylation_files, sample2_methylation_files = group_methylation_files_v2(
            all_methylation_files, fastq_files
        )
    else:
        logger.info("使用基于索引的分组方法（v1版本）")
        sample1_methylation_files, sample2_methylation_files = group_methylation_files_v1(
            all_methylation_files
        )
    logger.info(f"分组完成\n样本1: {sample1_methylation_files}\n样本2: {sample2_methylation_files}")

    # 合并样本1的所有甲基化信息
    if sample1_methylation_files:
        logger.info(f"合并样本1的 {len(sample1_methylation_files)} 个甲基化文件")
        merge_methylation_files(sample1_methylation_files, sample1_methylation)
    logger.info(f"成功生成样本1甲基化信息: {sample1_methylation}")

    # 合并样本2的所有甲基化信息
    if sample2_methylation_files:
        logger.info(f"合并样本2的 {len(sample2_methylation_files)} 个甲基化文件")
        merge_methylation_files(sample2_methylation_files, sample2_methylation)
    logger.info(f"成功生成样本2甲基化信息: {sample2_methylation}")

    return sample1_methylation, sample2_methylation


def validate_methylation_output(output_file):
    """
    验证甲基化提取结果

    Args:
        output_file: 甲基化提取输出文件路径
    """
    logger.info(f"验证甲基化提取结果: {output_file}")

    try:
        # 检查文件是否存在
        if not os.path.exists(output_file):
            logger.error(f"错误: 输出文件 {output_file} 不存在")
            return False

        # 检查文件大小
        file_size = os.path.getsize(output_file)
        if file_size == 0:
            logger.error(f"错误: 输出文件 {output_file} 为空")
            return False

        # 检查文件格式
        with open(output_file, "r", encoding="utf-8") as f:
            # 检查表头
            header = f.readline()
            if not header.startswith("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数"):
                logger.warning(f"警告: 输出文件 {output_file} 表头格式可能不正确")

            # 检查数据行
            line_count = 0
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) != 6:
                        logger.warning(
                            f"警告: 输出文件 {output_file} 第 {line_count+2} 行格式不正确"
                        )
                    line_count += 1

        logger.info(f"验证完成: {output_file} 包含 {line_count} 个甲基化位点")
        return True
    except Exception as e:
        logger.error(f"验证过程中出错: {str(e)}")
        return False


def generate_mock_methylation_data(output_file, min_level, max_level):
    """
    生成模拟的甲基化数据

    Args:
        output_file: 输出文件路径
        min_level: 最小甲基化水平
        max_level: 最大甲基化水平
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")
        num_sites = 1000
        range_level = max_level - min_level
        for i in range(1, num_sites + 1):
            chrom = "chr1"
            pos = 100000 + i * 100
            strand = "+" if i % 2 == 0 else "-"
            # 生成指定范围内的甲基化水平
            meth_level = round(min_level + (i % 50) / 50 * range_level, 2)
            meth_count = int(meth_level * 20)
            unmeth_count = 20 - meth_count
            f.write(f"{chrom}\t{pos}\t{strand}\t{meth_level}\t{meth_count}\t{unmeth_count}\n")


def merge_methylation_files(input_files, output_file):
    """
    合并多个甲基化文件的信息

    Args:
        input_files: 输入甲基化文件列表
        output_file: 输出文件路径
    """
    logger.info(f"合并 {len(input_files)} 个甲基化文件到 {output_file}")

    # 打开输出文件
    with open(output_file, "w", encoding="utf-8") as outfile:
        # 写入表头
        outfile.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")

        # 遍历所有输入文件
        for input_file in input_files:
            logger.info(f"处理文件: {input_file}")
            with open(input_file, "r", encoding="utf-8") as infile:
                # 跳过表头和 track 行
                header_skipped = False
                for line in infile:
                    if line.startswith("#") or line.startswith("track"):
                        if not header_skipped and line.startswith("#"):
                            header_skipped = True
                        continue
                    # 写入数据行
                    outfile.write(line)

    logger.info(f"合并完成: {output_file}")


def filter_methylation_data(methylation_file, temp_dir):
    """
    过滤甲基化数据

    Args:
        methylation_file: 甲基化文件路径
        temp_dir: 临时目录

    Returns:
        str: 过滤后的甲基化文件路径
    """
    logger.info(f"开始过滤甲基化数据: {methylation_file}")

    # 从配置文件中读取过滤相关配置
    repeat_bed = config.get("filtering", {}).get("repeat_bed")
    promoter_bed = config.get("filtering", {}).get("promoter_bed")
    depth_min = config.get("filtering", {}).get("min_depth", 5)
    depth_max = config.get("filtering", {}).get("max_depth", 50)

    # 调用过滤函数
    try:
        # 构建过滤后的输出文件路径
        base_name = os.path.basename(methylation_file)
        filtered_file = os.path.join(temp_dir, f"filtered_{base_name}")
        if os.path.exists(methylation_file):
            # 进行常规过滤
            filter_methylation(
                input_file=methylation_file,
                output_file=filtered_file,
                output_dir=temp_dir,
                repeat_bed=repeat_bed,
                promoter_bed=promoter_bed,
                depth_min=depth_min,
                depth_max=depth_max
            )
            logger.info(f"成功过滤甲基化数据: {filtered_file}")
            return filtered_file
        else:
            logger.warning(f"未找到对应的 counts 文件，跳过过滤步骤: {methylation_file}")
            return methylation_file
    except Exception as e:
        error_msg = f"过滤甲基化数据时出错: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def filter_to_candidate_regions_wrapper(methylation_file, temp_dir):
    """
    过滤甲基化数据到候选区域

    Args:
        methylation_file: 甲基化文件路径
        temp_dir: 临时目录

    Returns:
        str: 过滤到候选区域后的甲基化文件路径
    """
    logger.info(f"开始过滤甲基化数据到候选区域: {methylation_file}")
    
    # 候选区域文件路径
    candidate_regions = "data/OC_tissue_88_regions.hg38.txt"

    try:
        # 构建过滤后的输出文件路径
        base_name = os.path.basename(methylation_file)
        candidate_output = os.path.join(temp_dir, f"candidate_{base_name}")
        
        if os.path.exists(methylation_file) and os.path.exists(candidate_regions):
            if filter_to_candidate_regions(
                input_file=methylation_file,
                output_file=candidate_output,
                candidate_regions=candidate_regions,
                output_dir=temp_dir
            ):
                logger.info(f"成功过滤到候选区域: {candidate_output}")
                return candidate_output
            else:
                logger.warning("过滤到候选区域失败，使用原始文件")
                return methylation_file
        else:
            if not os.path.exists(candidate_regions):
                logger.warning(f"候选区域文件 {candidate_regions} 不存在，跳过过滤到候选区域步骤")
            return methylation_file
    except Exception as e:
        error_msg = f"过滤到候选区域时出错: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def perform_diff_analysis(sample1_methylation, sample2_methylation, result_dir):
    """
    进行差异甲基化分析

    Args:
        sample1_methylation: 样本1甲基化文件路径
        sample2_methylation: 样本2甲基化文件路径
        result_dir: 结果目录
        logger: 日志记录器
    """
    logger.info("开始进行差异甲基化分析")
    diff_result = os.path.join(result_dir, "diff_methylation.txt")

    try:
        run_diff_analysis(sample1_methylation, sample2_methylation, diff_result, threshold=0.2)
        logger.info("成功进行差异甲基化分析")
    except Exception as e:
        error_msg = f"进行差异甲基化分析时出错: {str(e)}"
        logger.error(error_msg)


def visualize_results(sample1_methylation, result_dir):
    """
    可视化结果

    Args:
        sample1_methylation: 样本1甲基化文件路径
        result_dir: 结果目录
        logger: 日志记录器
    """
    logger.info("开始可视化结果")
    density_plot = os.path.join(result_dir, "density_plot.png")
    scatter_plot = os.path.join(result_dir, "scatter_plot.png")
    heatmap_plot = os.path.join(result_dir, "heatmap_plot.png")

    try:
        visualize_methylation(sample1_methylation, density_plot, plot_type="density")
        logger.info("成功生成密度图")
    except Exception as e:
        error_msg = f"生成密度图时出错: {str(e)}"
        logger.error(error_msg)

    try:
        visualize_methylation(sample1_methylation, scatter_plot, plot_type="scatter")
        logger.info("成功生成散点图")
    except Exception as e:
        error_msg = f"生成散点图时出错: {str(e)}"
        logger.error(error_msg)

    try:
        visualize_methylation(sample1_methylation, heatmap_plot, plot_type="heatmap")
        logger.info("成功生成热图")
    except Exception as e:
        error_msg = f"生成热图时出错: {str(e)}"
        logger.error(error_msg)


def main():
    """
    主函数
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="执行完整的甲基化分析流程")
    parser.add_argument("files", nargs="+", help="FastQ文件或文件夹路径列表")

    # 添加模块开关参数
    parser.add_argument(
        "--skip-quality", action="store_true", default=True, help="跳过FastQ文件质量分析"
    )
    parser.add_argument("--skip-extraction", action="store_true", help="跳过甲基化信息提取")
    parser.add_argument("--skip-diff", action="store_true", help="跳过差异甲基化分析")
    parser.add_argument("--skip-visualization", action="store_true", help="跳过结果可视化")
    parser.add_argument(
        "--group-method",
        choices=["v1", "v2"],
        default="v2",
        help="分组方法: v1（基于索引）或 v2（基于Excel文件）",
    )

    args = parser.parse_args()

    # 查找所有FastQ文件
    logger.info("开始查找FastQ文件")
    fastq_files = find_fastq_files(args.files)
    if not fastq_files:
        logger.error("未找到任何FastQ文件")
        return

    logger.info(f"找到 {len(fastq_files)} 个FastQ文件:")
    for f in fastq_files:
        logger.info(f"  - {f}")

    # 创建带时间戳的子目录
    output_dirs = get_output_dirs(config)
    create_output_dirs(output_dirs)
    temp_dir = output_dirs["temp"]
    result_dir = output_dirs["result"]

    logger.info(f"创建目录:\n- 临时目录: {temp_dir}\n" f"- 结果目录: {result_dir}\n")

    # 1. 分析FastQ文件质量
    if not args.skip_quality:
        print("\n1. 分析FastQ文件质量...")
        analyze_quality(fastq_files)
    else:
        logger.info("1. 跳过FastQ文件质量分析")

    # 2. 提取甲基化信息
    sample1_methylation = None
    sample2_methylation = None
    if not args.skip_extraction:
        print("\n2. 提取甲基化信息...")
        sample1_methylation, sample2_methylation = extract_methylation_info(
            temp_dir, fastq_files, group_method=args.group_method
        )
    else:
        logger.info("2. 跳过甲基化信息提取")

    # 2.5 过滤甲基化数据
    filtered_sample1 = None
    filtered_sample2 = None
    if sample1_methylation or sample2_methylation:
        print("\n2.5 过滤甲基化数据...")
        if sample1_methylation:
            filtered_sample1 = filter_methylation_data(sample1_methylation, temp_dir)
        if sample2_methylation:
            filtered_sample2 = filter_methylation_data(sample2_methylation, temp_dir)
    else:
        logger.info("2.5 跳过甲基化数据过滤")

    # 2.6 过滤到候选区域
    candidate_sample1 = None
    candidate_sample2 = None
    if filtered_sample1 or filtered_sample2:
        print("\n2.6 过滤到候选区域...")
        if filtered_sample1:
            candidate_sample1 = filter_to_candidate_regions_wrapper(filtered_sample1, temp_dir)
        if filtered_sample2:
            candidate_sample2 = filter_to_candidate_regions_wrapper(filtered_sample2, temp_dir)
    else:
        logger.info("2.6 跳过过滤到候选区域步骤")

    # 3. 进行差异甲基化分析
    if not args.skip_diff and candidate_sample1 and candidate_sample2:
        print("\n3. 进行差异甲基化分析...")
        perform_diff_analysis(candidate_sample1, candidate_sample2, result_dir)
    else:
        logger.info("3. 跳过差异甲基化分析")

    # 4. 可视化结果
    if not args.skip_visualization and (candidate_sample1 or filtered_sample1 or sample1_methylation):
        print("\n4. 可视化结果...")
        # 使用候选区域的样本数据进行可视化，如果没有则使用过滤数据，如果没有过滤数据则使用原始数据
        if candidate_sample1:
            visualization_sample = candidate_sample1
        elif filtered_sample1:
            visualization_sample = filtered_sample1
        else:
            visualization_sample = sample1_methylation
        visualize_results(visualization_sample, result_dir)
    else:
        logger.info("4. 跳过结果可视化")

    print("\n分析完成！")
    print(f"结果保存在: {result_dir}")
    logger.info("分析完成")


if __name__ == "__main__":
    main()
