#!/usr/bin/env python3
"""
methylation_extractor.py - 从比对结果中提取甲基化信息的脚本

该脚本用于从Bismark等工具的比对结果中提取甲基化信息。

使用方法:
    python methylation_extractor.py <input_file> <output_file>

参数:
    input_file: 比对结果文件路径（.bam格式）
    output_file: 输出文件路径（.bedmethyl或.txt格式）
"""

import sys
import argparse
import os
import time
import subprocess
import shutil

# 导入配置模块
from utils.config import load_config
from utils.logger import logger


def extract_methylation(input_file, output_file, reference_genome=None, threads=8):
    """
    从比对结果中提取甲基化信息

    Args:
        input_file: 比对结果文件路径
        output_file: 输出文件路径
        reference_genome: 参考基因组文件路径，默认为None，使用配置文件中的值
        threads: 线程数，默认为8
    """
    # 加载配置文件
    config = load_config()
    # 如果没有提供参考基因组路径，使用配置文件中的值
    if reference_genome is None:
        reference_genome = config.get("extraction", {}).get("reference_genome")
    try:
        # 检查输入文件是否存在
        if not os.path.exists(input_file):
            logger.error(f"错误: 输入文件 '{input_file}' 不存在")
            sys.exit(1)

        # 检查文件格式
        if not input_file.endswith(".bam"):
            logger.error("错误: 输入文件必须是.bam格式")
            sys.exit(1)

        # 检查参考基因组文件是否存在
        if not os.path.exists(reference_genome):
            logger.error(f"错误: 参考基因组文件 '{reference_genome}' 不存在")
            sys.exit(1)

        logger.info(f"从 {input_file} 提取甲基化信息...")
        logger.info(f"输出到 {output_file}")
        logger.info(f"使用参考基因组: {reference_genome}")
        start_time = time.time()

        # 使用samtools和methyldackel提取甲基化信息
        # 首先使用samtools对BAM文件进行索引
        logger.info("正在为BAM文件创建索引...")
        index_cmd = f"samtools index {input_file}"
        logger.info(f"执行命令: {index_cmd}")
        subprocess.run(index_cmd, shell=True, check=True)

        # 使用methyldackel提取甲基化信息
        logger.info("正在使用MethylDackel提取甲基化信息...")  

        # 0.6.1
        methyldackel_cmd = f"MethylDackel extract {reference_genome} {input_file} -c CHG -c CHH -o {output_file} -@ {threads}"
        # 1.0.0
        # methyldackel_cmd = f"MethylDackel extract --CHG --CHH --genome {reference_genome} -o {output_file} {input_file}"
        
        logger.info(f"执行命令: {methyldackel_cmd}")
        subprocess.run(methyldackel_cmd, shell=True, check=True)

        # 查找MethylDackel生成的实际文件
        # MethylDackel会生成多个文件，我们需要找到包含CpG信息的文件
        cpg_file = f"{output_file}_CpG.counts.bedGraph"
        if os.path.exists(cpg_file):
            logger.info(f"找到MethylDackel输出文件: {cpg_file}")
            # 转换MethylDackel输出格式为我们需要的格式
            logger.info("正在转换输出格式...")
            convert_output(cpg_file, backup=True)
            # 将转换后的文件重命名为期望的输出文件名
            os.rename(cpg_file, output_file)
            logger.info(f"已将文件重命名为: {output_file}")

        logger.info(f"甲基化信息提取完成 time={time.time() - start_time:.2f}")

    except Exception as e:
        logger.error(f"错误: {str(e)}")
        sys.exit(1)


def convert_output(input_file, backup=False):
    """
    转换MethylDackel输出格式为标准格式

    Args:
        input_file: MethylDackel输出文件路径
        backup: 是否备份原文件，默认为False
    """
    temp_file = input_file + ".tmp"
    with open(input_file, "r", encoding="utf-8") as infile, open(temp_file, "w", encoding="utf-8") as outfile:
        # 写入表头
        outfile.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")
        
        for line in infile:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split()
            if len(parts) >= 4:
                chrom = parts[0]
                # 对于bedGraph格式，使用起始位置作为位置
                pos = parts[1]
                # 假设链为'+'，因为bedGraph格式通常不包含链信息
                strand = "+"
                # 尝试解析甲基化计数
                try:
                    # 检查是否为bedGraph格式（4列）
                    if len(parts) == 4:
                        # bedGraph格式: chrom, start, end, value
                        # 这里假设value是甲基化计数
                        meth_count = int(parts[3])
                        # 非甲基化计数设为0（根据实际情况调整）
                        unmeth_count = 0
                    else:
                        # 原始格式: chrom, pos, strand, meth_count, unmeth_count
                        meth_count = int(parts[3])
                        unmeth_count = int(parts[4]) if len(parts) >= 5 else 0
                    
                    # 计算甲基化水平
                    total_count = meth_count + unmeth_count
                    if total_count > 0:
                        meth_level = round(meth_count / total_count, 2)
                    else:
                        meth_level = 0.0
                except ValueError:
                    # 如果解析失败，使用默认值
                    meth_count = 0
                    unmeth_count = 0
                    meth_level = 0.0
                
                # 写入转换后的格式
                outfile.write(f"{chrom}\t{pos}\t{strand}\t{meth_level}\t{meth_count}\t{unmeth_count}\n")
    
    # 备份原文件（如果需要）
    if backup and os.path.exists(input_file):
        backup_file = input_file + ".bak"
        os.rename(input_file, backup_file)
        logger.info(f"已备份原文件到: {backup_file}")
    
    # 替换原文件
    os.replace(temp_file, input_file)


def align_fastq(fastq_files, output_bam, reference_genome):
    """
    将FastQ文件比对到参考基因组，生成BAM文件

    Args:
        fastq_files: FastQ文件列表
        output_bam: 输出BAM文件路径
        reference_genome: 参考基因组文件路径
    """
    try:
        logger.info(f"开始比对FastQ文件到参考基因组...")
        logger.info(f"参考基因组: {reference_genome}")
        logger.info(f"输出BAM文件: {output_bam}")
        start_time = time.time()

        # 检查参考基因组文件是否存在
        if not os.path.exists(reference_genome):
            logger.error(f"错误: 参考基因组文件 '{reference_genome}' 不存在")
            sys.exit(1)

        # 检查参考基因组索引是否存在
        index_files = [
            f"{reference_genome}.amb",
            f"{reference_genome}.ann",
            f"{reference_genome}.bwt",
            f"{reference_genome}.pac",
            f"{reference_genome}.sa"
        ]
        index_exists = all(os.path.exists(idx_file) for idx_file in index_files)

        if not index_exists:
            logger.info("参考基因组索引不存在，正在构建索引...")
            # 使用bwa构建索引
            index_cmd = f"bwa index {reference_genome}"
            logger.info(f"执行命令: {index_cmd}")
            subprocess.run(index_cmd, shell=True, check=True)

        # 准备FastQ文件参数
        fastq_args = " ".join(fastq_files)

        # 使用bwa mem进行比对
        logger.info("正在使用bwa mem进行比对...")
        sam_file = output_bam.replace(".bam", ".sam")
        align_cmd = f"bwa mem -M {reference_genome} {fastq_args} > {sam_file}"
        logger.info(f"执行命令: {align_cmd}")
        subprocess.run(align_cmd, shell=True, check=True)

        # 使用samtools将SAM转换为BAM并排序
        logger.info("正在将SAM转换为BAM并排序...")
        sort_cmd = f"samtools sort -o {output_bam} {sam_file}"
        logger.info(f"执行命令: {sort_cmd}")
        subprocess.run(sort_cmd, shell=True, check=True)

        # 清理临时SAM文件
        if os.path.exists(sam_file):
            os.remove(sam_file)

        logger.info(f"比对完成，生成BAM文件: {output_bam}")
        logger.info(f"比对耗时: {time.time() - start_time:.2f} 秒")

        return True

    except Exception as e:
        logger.error(f"比对过程中出错: {str(e)}")
        return False


if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="从比对结果中提取甲基化信息")
    parser.add_argument("input_file", help="比对结果文件路径（.bam格式）")
    parser.add_argument("output_file", help="输出文件路径（.bedmethyl或.txt格式）")

    args = parser.parse_args()

    # 执行提取
    extract_methylation(args.input_file, args.output_file)
