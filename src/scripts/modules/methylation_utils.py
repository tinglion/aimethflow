#!/usr/bin/env python3
"""
methylation_utils.py - 甲基化处理模块

该模块包含甲基化信息提取和处理相关的函数
"""

import os
import logging

logger = logging.getLogger(__name__)


def extract_methylation_info(temp_dir, fastq_files, group_method="v2", config=None, bam_files=None):
    """
    提取甲基化信息。

    Args:
        temp_dir: 临时目录
        fastq_files: FastQ文件列表（当 bam_files 为 None 时，用于生成BAM）
        group_method: 分组方法，"v1"（索引）或 "v2"（Excel）
        config: 配置对象
        bam_files: 预计算的BAM文件列表（血浆模式下由比对→去重→on-target步骤生成）。
                   提供此参数时跳过内部BAM生成（避免使用普通bwa mem进行bs-seq比对）。

    Returns:
        tuple: (sample1_methylation, sample2_methylation) 文件路径
    """
    logger.info("开始提取甲基化信息")

    reference_genome = config.get("extraction", {}).get("reference_genome")
    threads = config.get("extraction", {}).get("threads", 4)

    all_bam_files = []
    if bam_files:
        # 使用外部提供的BAM（血浆模式：已经过bismark比对+umi去重+on-target过滤）
        logger.info(f"使用预计算的BAM文件（共{len(bam_files)}个），跳过内部比对步骤")
        all_bam_files = bam_files
    elif fastq_files:
        # 回退：从FastQ重新生成BAM（仅用于非甲基化测序数据或测试）
        logger.warning("未提供预计算BAM，将使用 bwa mem 生成BAM（注意：此方式不适用于BS-seq数据！）")
        from methylation.bam_generator import generate_bam
        all_bam_files = generate_bam(fastq_files, reference_genome, threads=threads)

    # 提取所有BAM文件的甲基化信息
    all_methylation_files = []
    if all_bam_files:
        logger.info(f"为 {len(all_bam_files)} 个BAM文件提取甲基化信息")
        from methylation.methylation_extractor import extract_methylation
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
        from methylation.sample_grouping import group_methylation_files_v2
        sample1_methylation_files, sample2_methylation_files = group_methylation_files_v2(
            all_methylation_files
        )
    else:
        logger.info("使用基于索引的分组方法（v1版本）")
        from methylation.sample_grouping import group_methylation_files_v1
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


def filter_methylation_data(methylation_file, temp_dir, config):
    """
    过滤甲基化数据

    Args:
        methylation_file: 甲基化文件路径
        temp_dir: 临时目录
        config: 配置对象

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
            from methylation.filter_methylation import filter_methylation
            filter_methylation(
                input_file=methylation_file,
                output_file=filtered_file,
                output_dir=temp_dir,
                repeat_bed=repeat_bed,
                promoter_bed=promoter_bed,
                depth_min=depth_min,
                depth_max=depth_max,
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


def filter_to_candidate_regions_wrapper(methylation_file, temp_dir, config):
    """
    过滤甲基化数据到候选区域

    Args:
        methylation_file: 甲基化文件路径
        temp_dir: 临时目录
        config: 配置对象

    Returns:
        str: 过滤到候选区域后的甲基化文件路径
    """
    logger.info(f"开始过滤甲基化数据到候选区域: {methylation_file}")

    # 候选区域文件路径
    candidate_regions = config.get("data", {}).get("candidate_regions")

    try:
        # 构建过滤后的输出文件路径
        base_name = os.path.basename(methylation_file)
        candidate_output = os.path.join(temp_dir, f"candidate_{base_name}")

        if os.path.exists(methylation_file) and os.path.exists(candidate_regions):
            from methylation.candidate_regions import filter_to_candidate_regions
            if filter_to_candidate_regions(
                input_file=methylation_file,
                output_file=candidate_output,
                candidate_regions=candidate_regions,
                output_dir=temp_dir,
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
    """
    logger.info("开始进行差异甲基化分析")
    diff_result = os.path.join(result_dir, "diff_methylation.txt")

    try:
        from methylation.diff_methylation import run_diff_analysis
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
    """
    logger.info("开始可视化结果")
    density_plot = os.path.join(result_dir, "density_plot.png")
    scatter_plot = os.path.join(result_dir, "scatter_plot.png")
    heatmap_plot = os.path.join(result_dir, "heatmap_plot.png")

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, density_plot, plot_type="density")
        logger.info("成功生成密度图")
    except Exception as e:
        error_msg = f"生成密度图时出错: {str(e)}"
        logger.error(error_msg)

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, scatter_plot, plot_type="scatter")
        logger.info("成功生成散点图")
    except Exception as e:
        error_msg = f"生成散点图时出错: {str(e)}"
        logger.error(error_msg)

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, heatmap_plot, plot_type="heatmap")
        logger.info("成功生成热图")
    except Exception as e:
        error_msg = f"生成热图时出错: {str(e)}"
        logger.error(error_msg)
