#!/usr/bin/env python3
"""
candidate_regions.py - 候选区域处理模块

该模块用于处理候选区域的过滤和分析，将甲基化数据过滤到指定的候选区域中。
"""

import os
import subprocess
from utils.logger import logger


def intersect_with_regions(cpg_bed, region_bed, output_bed):
    """
    提取特定区域的CpG位点

    Args:
        cpg_bed: CpG位点BED文件路径
        region_bed: 区域BED文件路径
        output_bed: 输出BED文件路径
    """
    try:
        logger.info(f"提取特定区域的CpG位点")
        subprocess.run(
            ['bedtools', 'intersect', '-a', cpg_bed, '-b', region_bed],
            stdout=open(output_bed, 'w'),
            check=True
        )
        return True
    except Exception as e:
        logger.error(f"提取特定区域失败: {str(e)}")
        return False


def extract_methylation_ratio(input_file, filtered_bed, output_ratio):
    """
    提取过滤后位点的甲基化比例，并保存为标准BED六列格式

    Args:
        input_file: bedGraph文件路径
        filtered_bed: 过滤后的BED文件路径
        output_ratio: 输出文件路径
    """
    try:
        logger.info(f"提取甲基化比例并保存为标准BED六列格式")
        # 读取过滤后的位点
        filtered_sites = set()
        with open(filtered_bed, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 3:
                    key = f"{parts[0]}:{parts[1]}:{parts[2]}"
                    filtered_sites.add(key)
        
        # 提取甲基化比例并保存为标准BED六列格式
        with open(input_file, 'r', encoding='utf-8') as infile, open(output_ratio, 'w', encoding='utf-8') as outfile:
            # 写入表头
            outfile.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")
            
            for line in infile:
                if not line.startswith('track') and not line.startswith('#'):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        chrom = parts[0]
                        start = parts[1]
                        end = parts[2]
                        key = f"{chrom}:{start}:{end}"
                        if key in filtered_sites:
                            # 检查文件格式
                            if len(parts) == 4:
                                # 4列格式：chrom, start, end, methylation_count
                                # 假设非甲基化计数为0，甲基化比例为100%
                                try:
                                    meth_count = int(parts[3])
                                    unmeth_count = 0
                                    total_count = meth_count + unmeth_count
                                    if total_count > 0:
                                        meth_level = round(meth_count / total_count, 2)
                                    else:
                                        meth_level = 0.0
                                except ValueError:
                                    continue
                            elif len(parts) >= 6:
                                # 6列格式：chrom, start, end, methylation_percentage, meth_count, unmeth_count
                                try:
                                    meth_count = int(parts[4])
                                    unmeth_count = int(parts[5])
                                    total_count = meth_count + unmeth_count
                                    if total_count > 0:
                                        meth_level = round(meth_count / total_count, 2)
                                    else:
                                        meth_level = 0.0
                                except ValueError:
                                    continue
                            else:
                                # 其他格式，跳过
                                continue
                            # 输出标准BED六列格式
                            outfile.write(f"{chrom}\t{start}\t+\t{meth_level}\t{meth_count}\t{unmeth_count}\n")
        return True
    except Exception as e:
        logger.error(f"提取甲基化比例失败: {str(e)}")
        return False


def calculate_percentage(input_file):
    """
    计算100%位点的占比

    Args:
        input_file: 甲基化数据文件路径（标准BED六列格式）

    Returns:
        tuple: (total_sites, hundred_percent_sites, percentage)
    """
    try:
        total_sites = 0
        hundred_percent_sites = 0
        with open(input_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) >= 4:
                    total_sites += 1
                    try:
                        meth_level = float(parts[3])
                        if meth_level == 1.0:
                            hundred_percent_sites += 1
                    except ValueError:
                        pass
        if total_sites > 0:
            percentage = (hundred_percent_sites / total_sites) * 100
        else:
            percentage = 0.0
        return total_sites, hundred_percent_sites, percentage
    except Exception as e:
        logger.error(f"计算百分比失败: {str(e)}")
        return 0, 0, 0.0


def filter_to_candidate_regions(input_file, output_file, candidate_regions, output_dir="temp"):
    """
    将甲基化数据过滤到候选区域

    Args:
        input_file: MethylDackel生成的bedGraph文件路径
        output_file: 输出过滤后的文件路径
        candidate_regions: 候选区域文件路径
        output_dir: 输出目录，默认为temp
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 获取输入文件的基本名称（不含扩展名），用于生成唯一的临时文件名
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    # 步骤1：将bedGraph文件转换为BED格式
    cpg_sites_bed = os.path.join(output_dir, f"{base_name}_cpg_sites.bed")
    if not convert_to_bed(input_file, cpg_sites_bed):
        logger.error("转换bedGraph文件为BED格式失败")
        return False

    # 步骤2：提取候选区域的CpG位点
    cpg_candidate_bed = os.path.join(output_dir, f"{base_name}_cpg_candidate.bed")
    if not intersect_with_regions(cpg_sites_bed, candidate_regions, cpg_candidate_bed):
        logger.error("提取候选区域的CpG位点失败")
        return False

    # 步骤3：提取甲基化比例并保存为标准BED六列格式
    if not extract_methylation_ratio(input_file, cpg_candidate_bed, output_file):
        logger.error("提取甲基化比例失败")
        return False

    # 步骤4：计算统计信息
    total, hundred, percent = calculate_percentage(output_file)
    logger.info(f"100%位点数量：{hundred}")
    logger.info(f"总位点数：{total}")
    logger.info(f"100%位点占比：{percent:.2f}%")

    logger.info(f"候选区域过滤结果已保存到: {output_file}")
    return True


def convert_to_bed(input_file, output_bed):
    """
    将bedGraph文件转换为BED格式

    Args:
        input_file: bedGraph文件路径
        output_bed: 输出BED文件路径
    """
    try:
        logger.info(f"转换为BED格式: {input_file} -> {output_bed}")
        with open(input_file, 'r', encoding='utf-8') as infile, open(output_bed, 'w', encoding='utf-8') as outfile:
            for line in infile:
                if not line.startswith('track') and not line.startswith('#'):
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        outfile.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")
        return True
    except Exception as e:
        logger.error(f"转换为BED格式失败: {str(e)}")
        return False
