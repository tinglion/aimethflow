#!/usr/bin/env python3
"""
methylation_utils.py - Methylation processing module

Functions for methylation extraction and processing.
"""

import os
import logging

logger = logging.getLogger(__name__)


def extract_methylation_info(temp_dir, fastq_files, group_method="v2", config=None, bam_files=None, category="plasma"):
    """
    提取甲基化信息并按样本分组
    
    Args:
        temp_dir (str): 工作目录，用于存储临时文件和输出结果
        fastq_files (list): FastQ 文件列表（当未提供 bam_files 时用于生成 BAM 文件）
        group_method (str, optional): 分组方法，可选 "v1"（基于索引）或 "v2"（基于 Excel）。默认为 "v2"
        config (dict, optional): 配置字典，包含参考基因组和线程数等参数。默认为 None
        bam_files (list, optional): 预计算的 BAM 文件列表（血浆模式：由比对→去重→靶向捕获步骤生成）。
                                     当提供时，跳过内部 BAM 生成（避免对 BS-seq 数据运行普通 bwa mem）。默认为 None
        category (str, optional): 分组类别，例如 "plasma" 或 "tissue"。默认为 "plasma"
    
    Returns:
        tuple: 包含三个元素的元组 (sample1_methylation, sample2_methylation, borderline_methylation)，
               分别为样本组1、样本组2和边界样本的甲基化文件路径
    """
    """
    Extract methylation information.

    Args:
        temp_dir:     working directory
        fastq_files:  FastQ file list (used to generate BAM when bam_files is None)
        group_method: grouping method — "v1" (index-based) or "v2" (Excel-based)
        config:       config dict
        bam_files:    pre-computed BAM list (plasma mode: produced by align → dedup →
                      on-target steps). When provided, internal BAM generation is skipped
                      (avoids running plain bwa mem on BS-seq data).
        category:     category for grouping (e.g., "plasma" or "tissue")

    Returns:
        tuple: (sample1_methylation, sample2_methylation) file paths
    """
    logger.info("Starting methylation extraction")

    reference_genome = config.get("extraction", {}).get("reference_genome")
    threads = config.get("extraction", {}).get("threads", 4)

    all_bam_files = []
    if bam_files:
        # use externally supplied BAMs (plasma mode: already bismark-aligned, UMI-deduped, on-target)
        logger.info(f"Using {len(bam_files)} pre-computed BAM file(s) — skipping internal alignment")
        all_bam_files = bam_files
    elif fastq_files:
        # fallback: generate BAM from FastQ (only for non-BS-seq data or testing)
        logger.warning("No pre-computed BAMs provided — generating BAM with bwa mem (NOTE: not suitable for BS-seq data!)")
        from methylation.bam_generator import generate_bam
        all_bam_files = generate_bam(fastq_files, reference_genome, threads=threads)

    all_methylation_files = []
    if all_bam_files:
        logger.info(f"Extracting methylation from {len(all_bam_files)} BAM file(s)")
        from methylation.methylation_extractor import extract_methylation
        for i, bam_file in enumerate(all_bam_files):
            individual_output = bam_file.replace(".bam", ".txt")
            logger.info(f"extracting methylation: {bam_file}")
            extract_methylation(bam_file, individual_output, reference_genome=reference_genome)
            validate_methylation_output(individual_output)
            all_methylation_files.append(individual_output)

    sample1_methylation = os.path.join(temp_dir, "sample1_methylation.txt")
    sample2_methylation = os.path.join(temp_dir, "sample2_methylation.txt")
    borderline_methylation = os.path.join(temp_dir, "borderline_methylation.txt")

    if group_method == "v2":
        logger.info("using Excel-based grouping (v2)")
        from methylation.sample_grouping import group_methylation_files_v2
        sample1_methylation_files, sample2_methylation_files, borderline_methylation_files = group_methylation_files_v2(
            all_methylation_files, category=category
        )
    else:
        logger.info("using index-based grouping (v1)")
        from methylation.sample_grouping import group_methylation_files_v1
        sample1_methylation_files, sample2_methylation_files = group_methylation_files_v1(
            all_methylation_files
        )
        borderline_methylation_files = []
    logger.info(f"Grouping done\ngroup1: {sample1_methylation_files}\ngroup2: {sample2_methylation_files}")

    if sample1_methylation_files:
        logger.info(f"merging {len(sample1_methylation_files)} methylation file(s) for group 1")
        merge_methylation_files(sample1_methylation_files, sample1_methylation)
    logger.info(f"group 1 methylation ready: {sample1_methylation}")

    if sample2_methylation_files:
        logger.info(f"merging {len(sample2_methylation_files)} methylation file(s) for group 2")
        merge_methylation_files(sample2_methylation_files, sample2_methylation)
    logger.info(f"group 2 methylation ready: {sample2_methylation}")

    if borderline_methylation_files:
        logger.info(f"merging {len(borderline_methylation_files)} methylation file(s) for borderline samples")
        merge_methylation_files(borderline_methylation_files, borderline_methylation)
        logger.info(f"borderline samples methylation ready: {borderline_methylation}")

    return sample1_methylation, sample2_methylation, borderline_methylation


def validate_methylation_output(output_file):
    """
    Validate methylation extraction output.

    Args:
        output_file: path to the methylation extraction output file
    """
    logger.info(f"Validating methylation output: {output_file}")

    try:
        if not os.path.exists(output_file):
            logger.error(f"output file not found: {output_file}")
            return False

        file_size = os.path.getsize(output_file)
        if file_size == 0:
            logger.error(f"output file is empty: {output_file}")
            return False

        with open(output_file, "r", encoding="utf-8") as f:
            header = f.readline()
            if not header.startswith("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数"):
                logger.warning(f"unexpected header in {output_file}")

            line_count = 0
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) != 6:
                        logger.warning(f"malformed line {line_count+2} in {output_file}")
                    line_count += 1

        logger.info(f"Validation done: {output_file} has {line_count} methylation sites")
        return True
    except Exception as e:
        logger.error(f"validation error: {str(e)}")
        return False


def generate_mock_methylation_data(output_file, min_level, max_level):
    """
    Generate mock methylation data for testing.

    Args:
        output_file: output file path
        min_level:   minimum methylation level
        max_level:   maximum methylation level
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")
        num_sites = 1000
        range_level = max_level - min_level
        for i in range(1, num_sites + 1):
            chrom = "chr1"
            pos = 100000 + i * 100
            strand = "+" if i % 2 == 0 else "-"
            meth_level = round(min_level + (i % 50) / 50 * range_level, 2)
            meth_count = int(meth_level * 20)
            unmeth_count = 20 - meth_count
            f.write(f"{chrom}\t{pos}\t{strand}\t{meth_level}\t{meth_count}\t{unmeth_count}\n")


def merge_methylation_files(input_files, output_file):
    """
    Merge multiple methylation files.

    Args:
        input_files: list of input methylation file paths
        output_file: output file path
    """
    logger.info(f"Merging {len(input_files)} methylation file(s) → {output_file}")

    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")

        for input_file in input_files:
            logger.info(f"  merging: {input_file}")
            with open(input_file, "r", encoding="utf-8") as infile:
                # skip header and track lines
                header_skipped = False
                for line in infile:
                    if line.startswith("#") or line.startswith("track"):
                        if not header_skipped and line.startswith("#"):
                            header_skipped = True
                        continue
                    outfile.write(line)

    logger.info(f"Merge done: {output_file}")


def filter_methylation_data(methylation_file, temp_dir, config):
    """
    Filter methylation data.

    Args:
        methylation_file: methylation file path
        temp_dir:         working directory
        config:           config dict

    Returns:
        str: path to the filtered methylation file
    """
    logger.info(f"Filtering methylation data: {methylation_file}")

    # load filtering config
    repeat_bed = config.get("filtering", {}).get("repeat_bed")
    promoter_bed = config.get("filtering", {}).get("promoter_bed")
    depth_min = config.get("filtering", {}).get("min_depth", 5)
    depth_max = config.get("filtering", {}).get("max_depth", 50)

    try:
        base_name = os.path.basename(methylation_file)
        filtered_file = os.path.join(temp_dir, f"filtered_{base_name}")
        if os.path.exists(methylation_file):
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
            logger.info(f"Filtering done: {filtered_file}")
            return filtered_file
        else:
            logger.warning(f"input file not found, skipping filter: {methylation_file}")
            return methylation_file
    except Exception as e:
        error_msg = f"error filtering methylation data: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def filter_to_candidate_regions_wrapper(methylation_file, temp_dir, config):
    """
    Filter methylation data to candidate regions.

    Args:
        methylation_file: methylation file path
        temp_dir:         working directory
        config:           config dict

    Returns:
        str: path to the candidate-region-filtered methylation file
    """
    logger.info(f"Filtering to candidate regions: {methylation_file}")

    candidate_regions = config.get("data", {}).get("candidate_regions")

    try:
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
                logger.info(f"Candidate-region filter done: {candidate_output}")
                return candidate_output
            else:
                logger.warning("Candidate-region filter failed, using original file")
                return methylation_file
        else:
            if not os.path.exists(candidate_regions):
                logger.warning(f"Candidate regions file not found, skipping: {candidate_regions}")
            return methylation_file
    except Exception as e:
        error_msg = f"error filtering to candidate regions: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def perform_diff_analysis(sample1_methylation, sample2_methylation, result_dir):
    """
    Run differential methylation analysis.

    Args:
        sample1_methylation: group 1 methylation file path
        sample2_methylation: group 2 methylation file path
        result_dir:          results directory
    """
    logger.info("Starting differential methylation analysis")
    diff_result = os.path.join(result_dir, "diff_methylation.txt")

    try:
        from methylation.diff_methylation import run_diff_analysis
        run_diff_analysis(sample1_methylation, sample2_methylation, diff_result, threshold=0.2)
        logger.info("Differential methylation analysis done")
    except Exception as e:
        error_msg = f"differential methylation analysis error: {str(e)}"
        logger.error(error_msg)


def visualize_results(sample1_methylation, result_dir):
    """
    Visualize analysis results.

    Args:
        sample1_methylation: group 1 methylation file path
        result_dir:          results directory
    """
    logger.info("Starting result visualization")
    density_plot = os.path.join(result_dir, "density_plot.png")
    scatter_plot = os.path.join(result_dir, "scatter_plot.png")
    heatmap_plot = os.path.join(result_dir, "heatmap_plot.png")

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, density_plot, plot_type="density")
        logger.info("Density plot saved")
    except Exception as e:
        error_msg = f"error generating density plot: {str(e)}"
        logger.error(error_msg)

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, scatter_plot, plot_type="scatter")
        logger.info("Scatter plot saved")
    except Exception as e:
        error_msg = f"error generating scatter plot: {str(e)}"
        logger.error(error_msg)

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, heatmap_plot, plot_type="heatmap")
        logger.info("Heatmap saved")
    except Exception as e:
        error_msg = f"error generating heatmap: {str(e)}"
        logger.error(error_msg)
