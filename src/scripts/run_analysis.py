#!/usr/bin/env python3
"""
run_analysis.py - 执行完整的甲基化分析流程

该脚本用于执行完整的甲基化分析流程，包括：
1. 初始化与配置
2. 文件检测与准备
3. 分析FastQ文件质量（可选）
4. 数据预处理（质量控制和接头去除）
5. 甲基化序列比对
6. PCR重复去除
7. 靶向区域提取
8. 提取甲基化信息
9. 过滤甲基化数据（包括候选区域过滤）
10. 进行差异甲基化分析（可选）
11. 可视化结果（可选）
12. 生成分析报告
13. 清理

使用方法:
    python src/scripts/run_analysis.py [选项] <文件/文件夹...>

参数:
    文件/文件夹: 可以是一个或多个FastQ文件或包含FastQ文件的文件夹

选项:
    --skip-quality: 跳过FastQ文件质量分析
    --skip-preprocessing: 跳过数据预处理
    --skip-alignment: 跳过序列比对
    --skip-deduplication: 跳过PCR重复去除
    --skip-target-extraction: 跳过靶向区域提取
    --skip-extraction: 跳过甲基化信息提取
    --skip-filtering: 跳过甲基化数据过滤
    --skip-candidate-regions: 跳过候选区域过滤
    --skip-diff: 跳过差异甲基化分析
    --skip-visualization: 跳过结果可视化
    --skip-report: 跳过报告生成
    --clean-temp: 分析完成后清理临时文件
    --config FILE: 使用指定的配置文件
    --verbose: 启用详细输出
"""

import os
import sys
import time
import argparse
import logging
import shutil

# 导入各个模块的函数
from utils.logger import setup_logger
from utils.config import load_config, get_output_dirs, create_output_dirs
from scripts.modules import (
    find_files,
    analyze_quality,
    preprocess_fastq,
    align_methylation,
    deduplicate_bam,
    extract_target_regions,
    extract_methylation_info,
    merge_methylation_files,
    filter_methylation_data,
    filter_to_candidate_regions_wrapper,
    perform_diff_analysis,
    visualize_results,
    generate_report,
)

# 加载配置文件
config = load_config()

# 生成日志文件路径
log_dir = config["directories"]["log_dir"]
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/analysis_{time.strftime('%Y%m%d_%H%M%S')}.log"
logger = setup_logger(log_file=log_file, console_level=logging.INFO, file_level=logging.DEBUG)


def check_presets():
    """
    检查分析前的预设信息

    Returns:
        bool: 预设信息检查是否通过
    """
    logger.info("开始检查预设信息")

    # 1. 检查探针区域文件
    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if probe_panel:
        if not os.path.exists(probe_panel):
            logger.error(f"探针区域文件不存在: {probe_panel}")
            logger.warning("将跳过靶向区域提取步骤")
        else:
            logger.info(f"探针区域文件: {probe_panel}")
    else:
        logger.warning("未配置探针区域文件，跳过靶向区域提取")

    # 2. 检查参考基因组
    reference_genome = config.get("extraction", {}).get("reference_genome")
    if not reference_genome or not os.path.exists(reference_genome):
        logger.error(f"参考基因组不存在: {reference_genome}")
        return False
    logger.info(f"参考基因组: {reference_genome}")

    # 3. 检查Bismark索引
    bismark_index = os.path.dirname(reference_genome)
    if not os.path.exists(os.path.join(bismark_index, "CT_conversion", "Genome", "genome.2.bt2")):
        logger.warning("Bismark索引未找到，可能需要运行bismark_genome_preparation")
        logger.warning("将跳过比对步骤")

    # 4. 检查必要的工具
    required_tools = ["trim_galore", "bismark", "deduplicate_bismark", "bedtools"]
    missing_tools = []
    for tool in required_tools:
        try:
            import subprocess

            subprocess.run([tool, "--version"], capture_output=True, check=True)
            logger.info(f"工具 {tool} 已安装")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning(f"工具 {tool} 未安装")
            missing_tools.append(tool)

    if missing_tools:
        logger.warning(f"以下工具未安装: {', '.join(missing_tools)}")
        logger.warning("将跳过需要这些工具的步骤")

    logger.info("预设信息检查完成")
    return True


def main():
    """
    主函数
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="执行完整的甲基化分析流程")
    parser.add_argument("files", nargs="+", help="FastQ文件或文件夹路径列表")

    # 添加模块开关参数
    parser.add_argument(
        "--skip-quality",
        action="store_true",
        default=True,
        help="跳过FastQ文件质量分析",
    )
    parser.add_argument(
        "--skip-preprocessing",
        action="store_true",
        default=True,
        help="跳过数据预处理",
    )
    parser.add_argument("--skip-alignment", action="store_true", default=True, help="跳过序列比对")
    parser.add_argument(
        "--skip-deduplication",
        action="store_true",
        default=True,
        help="跳过PCR重复去除",
    )
    parser.add_argument(
        "--skip-target-extraction",
        action="store_true",
        default=True,
        help="跳过靶向区域提取",
    )
    parser.add_argument("--skip-extraction", action="store_true", help="跳过甲基化信息提取")
    parser.add_argument("--skip-filtering", action="store_true", help="跳过甲基化数据过滤")
    parser.add_argument("--skip-candidate-regions", action="store_true", help="跳过候选区域过滤")
    parser.add_argument("--skip-diff", action="store_true", help="跳过差异甲基化分析")
    parser.add_argument("--skip-visualization", action="store_true", help="跳过结果可视化")
    parser.add_argument("--skip-report", action="store_true", help="跳过报告生成")
    parser.add_argument(
        "--clean-temp",
        action="store_true",
        default=False,
        help="分析完成后清理临时文件",
    )
    parser.add_argument("--config", type=str, help="使用指定的配置文件")
    parser.add_argument("--verbose", action="store_true", default=False, help="启用详细输出")
    parser.add_argument(
        "--group-method",
        choices=["v1", "v2"],
        default="v2",
        help="分组方法: v1（基于索引）或 v2（基于Excel文件）",
    )

    args = parser.parse_args()

    # 处理自定义配置文件
    global config
    if args.config and os.path.exists(args.config):
        logger.info(f"使用自定义配置文件: {args.config}")
        from utils.config import load_config as load_custom_config

        config = load_custom_config(args.config)

    # 处理详细输出
    if args.verbose:
        # 重新设置日志级别为DEBUG
        logger.setLevel(logging.DEBUG)
        console_handler = next((h for h in logger.handlers if h.name == "console"), None)
        if console_handler:
            console_handler.setLevel(logging.DEBUG)

    # 检查预设信息
    if not check_presets():
        logger.error("预设信息检查失败，分析终止")
        return

    # 查找所有支持的文件类型
    logger.info("开始查找支持的文件类型")
    found_files = find_files(args.files)

    # 统计文件数量
    total_files = sum(len(files) for files in found_files.values())
    if total_files == 0:
        logger.error("未找到任何支持的文件")
        return

    # 打印找到的文件
    logger.info(f"找到 {total_files} 个文件:")
    for file_type, files in found_files.items():
        if files:
            logger.info(f"  - {file_type}: {len(files)} 个文件")
            for f in files:
                logger.info(f"    - {f}")

    # 确定输入文件类型
    fastq_files = found_files["fastq"]
    bam_files = found_files["bam"]
    methylation_files = found_files["methylation"]

    # 创建带时间戳的子目录
    output_dirs = get_output_dirs(config)
    create_output_dirs(output_dirs)
    temp_dir = output_dirs["temp"]
    result_dir = output_dirs["result"]

    logger.info(f"创建目录:\n- 临时目录: {temp_dir}\n" f"- 结果目录: {result_dir}\n")

    # 1. 初始化与配置
    logger.info("1. 初始化与配置完成")

    # 2. 文件检测与准备
    logger.info("2. 文件检测与准备完成")

    # 3. 分析FastQ文件质量
    if not args.skip_quality:
        print("\n3. 分析FastQ文件质量...")
        analyze_quality(fastq_files)
    else:
        logger.info("3. 跳过FastQ文件质量分析")

    # 4. 数据预处理
    processed_files = fastq_files
    if not args.skip_preprocessing and fastq_files:
        print("\n4. 数据预处理...")
        processed_files = preprocess_fastq(fastq_files, temp_dir)
        if not processed_files:
            logger.error("数据预处理失败，使用原始文件继续")
            processed_files = fastq_files
    else:
        logger.info("4. 跳过数据预处理")

    # 5. 甲基化序列比对
    if not bam_files:
        if not args.skip_alignment and processed_files:
            print("\n5. 甲基化序列比对...")
            bam_files = align_methylation(processed_files, temp_dir, config)
        else:
            logger.info("5. 跳过甲基化序列比对")
    else:
        logger.info("5. 使用提供的BAM文件")

    # 6. PCR重复去除
    deduplicated_files = bam_files
    if not args.skip_deduplication and bam_files:
        print("\n6. PCR重复去除...")
        deduplicated_files = deduplicate_bam(bam_files, temp_dir)
    else:
        logger.info("6. 跳过PCR重复去除")

    # 7. 靶向区域提取
    target_files = deduplicated_files
    if not args.skip_target_extraction and deduplicated_files:
        print("\n7. 靶向区域提取...")
        target_files = extract_target_regions(deduplicated_files, temp_dir, config)
    else:
        logger.info("7. 跳过靶向区域提取")

    # 8. 提取甲基化信息
    sample1_methylation = None
    sample2_methylation = None
    if not methylation_files:
        if not args.skip_extraction:
            print("\n8. 提取甲基化信息...")
            sample1_methylation, sample2_methylation = extract_methylation_info(
                temp_dir, fastq_files, group_method=args.group_method, config=config
            )
        else:
            logger.info("8. 跳过甲基化信息提取")
    else:
        logger.info("8. 使用提供的甲基化数据文件")
        # 根据文件数量和分组方法处理甲基化文件
        if len(methylation_files) >= 2:
            # 如果提供了多个甲基化文件，根据分组方法进行分组
            from methylation.sample_grouping import (
                group_methylation_files_v1,
                group_methylation_files_v2,
            )

            if args.group_method == "v2":
                sample1_methylation_files, sample2_methylation_files = group_methylation_files_v2(
                    methylation_files, fastq_files
                )
            else:
                sample1_methylation_files, sample2_methylation_files = group_methylation_files_v1(
                    methylation_files
                )

            # 合并样本1的所有甲基化信息
            sample1_methylation = os.path.join(temp_dir, "sample1_methylation.txt")
            if sample1_methylation_files:
                merge_methylation_files(sample1_methylation_files, sample1_methylation)

            # 合并样本2的所有甲基化信息
            sample2_methylation = os.path.join(temp_dir, "sample2_methylation.txt")
            if sample2_methylation_files:
                merge_methylation_files(sample2_methylation_files, sample2_methylation)
        elif len(methylation_files) == 1:
            # 如果只提供了一个甲基化文件，同时用作样本1和样本2
            sample1_methylation = methylation_files[0]
            sample2_methylation = methylation_files[0]
        logger.info(f"使用提供的甲基化文件: {methylation_files}")

    # 9. 过滤甲基化数据
    filtered_sample1 = None
    filtered_sample2 = None
    candidate_sample1 = None
    candidate_sample2 = None
    if not args.skip_filtering and (sample1_methylation or sample2_methylation):
        print("\n9. 过滤甲基化数据...")
        if sample1_methylation:
            filtered_sample1 = filter_methylation_data(sample1_methylation, temp_dir, config)
        if sample2_methylation:
            filtered_sample2 = filter_methylation_data(sample2_methylation, temp_dir, config)
        
        # 过滤到候选区域
        if not args.skip_candidate_regions:
            print("\n9. 过滤到候选区域...")
            if filtered_sample1:
                candidate_sample1 = filter_to_candidate_regions_wrapper(filtered_sample1, temp_dir, config)
            if filtered_sample2:
                candidate_sample2 = filter_to_candidate_regions_wrapper(filtered_sample2, temp_dir, config)
    else:
        logger.info("9. 跳过甲基化数据过滤")

    # 10. 进行差异甲基化分析
    if not args.skip_diff and candidate_sample1 and candidate_sample2:
        print("\n10. 进行差异甲基化分析...")
        perform_diff_analysis(candidate_sample1, candidate_sample2, result_dir)
    else:
        logger.info("10. 跳过差异甲基化分析")

    # 11. 可视化结果
    if not args.skip_visualization and (
        candidate_sample1 or filtered_sample1 or sample1_methylation
    ):
        print("\n11. 可视化结果...")
        # 使用候选区域的样本数据进行可视化，如果没有则使用过滤数据，如果没有过滤数据则使用原始数据
        if candidate_sample1:
            visualization_sample = candidate_sample1
        elif filtered_sample1:
            visualization_sample = filtered_sample1
        else:
            visualization_sample = sample1_methylation
        visualize_results(visualization_sample, result_dir)
    else:
        logger.info("11. 跳过结果可视化")

    # 12. 生成分析报告
    if not args.skip_report:
        print("\n12. 生成分析报告...")
        generate_report(
            temp_dir,
            result_dir,
            sample1_methylation,
            sample2_methylation,
            candidate_sample1,
            candidate_sample2,
        )
    else:
        logger.info("12. 跳过报告生成")

    # 13. 清理
    if args.clean_temp:
        print("\n13. 清理临时文件...")
        try:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
                logger.info(f"13. 临时目录已清理: {temp_dir}")
                print("临时文件清理完成")
        except Exception as e:
            logger.error(f"清理临时文件时出错: {str(e)}")
            print("临时文件清理失败")
    else:
        logger.info("13. 跳过清理步骤")

    print("\n分析完成！")
    print(f"结果保存在: {result_dir}")
    logger.info("分析完成")


if __name__ == "__main__":
    main()
