#!/usr/bin/env python3
"""
bam_generator.py - 从FastQ文件生成BAM文件的脚本

该脚本用于将FastQ文件比对到参考基因组，生成BAM文件。
"""

import os
import subprocess
import time
import argparse
from concurrent.futures import ThreadPoolExecutor


def check_and_build_index(reference_genome):
    """
    检查参考基因组索引是否存在，不存在则构建

    Args:
        reference_genome: 参考基因组文件路径

    Returns:
        bool: 索引是否成功构建
    """
    # 检查参考基因组文件是否存在
    if not os.path.exists(reference_genome):
        print(f"错误: 参考基因组文件 '{reference_genome}' 不存在")
        return False

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
        print("参考基因组索引不存在，正在构建索引...")
        # 使用bwa构建索引
        index_cmd = f"bwa index {reference_genome}"
        print(f"执行命令: {index_cmd}")
        try:
            subprocess.run(index_cmd, shell=True, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"构建索引失败: {str(e)}")
            return False
    return True


def process_single_file(fastq_file, reference_genome, output_dir=None, threads=1):
    """
    处理单个FastQ文件，生成对应的BAM文件

    Args:
        fastq_file: FastQ文件路径
        reference_genome: 参考基因组文件路径
        output_dir: 输出目录，如果为None则在源文件所在目录生成

    Returns:
        str: 生成的BAM文件路径，如果失败则返回None
    """
    # 生成BAM文件名
    if output_dir:
        base_name = os.path.basename(os.path.splitext(os.path.splitext(fastq_file)[0])[0])
        bam_file = os.path.join(output_dir, f"{base_name}.bam")
    else:
        base_name = os.path.splitext(os.path.splitext(fastq_file)[0])[0]
        bam_file = f"{base_name}.bam"
    
    # 检查BAM文件是否已经存在
    if os.path.exists(bam_file):
        print(f"BAM文件已存在: {bam_file}")
        return bam_file
    
    print(f"为 {fastq_file} 生成BAM文件...")
    print(f"输出BAM文件: {bam_file}")
    print(f"使用参考基因组: {reference_genome}")
    start_time = time.time()
    
    # 检查参考基因组索引
    if not check_and_build_index(reference_genome):
        return None
    
    # 使用bwa mem进行比对
    print("正在使用bwa mem进行比对...")
    sam_file = bam_file.replace(".bam", ".sam")
    align_cmd = f"bwa mem -M -t {threads} {reference_genome} {fastq_file} > {sam_file}"
    print(f"执行命令: {align_cmd}")
    try:
        subprocess.run(align_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"比对失败: {str(e)}")
        if os.path.exists(sam_file):
            os.remove(sam_file)
        return None
    
    # 使用samtools将SAM转换为BAM并排序
    print("正在将SAM转换为BAM并排序...")
    sort_cmd = f"samtools sort  -@ {threads} -m 4G -o {bam_file} {sam_file}"
    print(f"执行命令: {sort_cmd}")
    try:
        subprocess.run(sort_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"转换失败: {str(e)}")
        if os.path.exists(sam_file):
            os.remove(sam_file)
        return None
    
    # 清理临时SAM文件
    if os.path.exists(sam_file):
        os.remove(sam_file)
    
    print(f"BAM文件生成完成: {bam_file}")
    print(f"耗时: {time.time() - start_time:.2f} 秒")
    return bam_file


def generate_bam(fastq_files, reference_genome,parallel=1, threads=1, output_dir=None):
    """
    为每个FastQ文件生成对应的BAM文件

    Args:
        fastq_files: FastQ文件列表
        reference_genome: 参考基因组文件路径
        parallel: 并发处理数，默认为1
        threads: 并发线程数，默认为1
        output_dir: 输出目录，如果为None则在源文件所在目录生成

    Returns:
        list: 生成的BAM文件路径列表
    """
    bam_files = []
    
    # 首先检查并构建参考基因组索引
    if not check_and_build_index(reference_genome):
        return bam_files
    
    # 使用线程池进行并发处理
    if parallel > 1:
        print(f"使用 {parallel} 个线程进行并发处理")
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            # 提交所有任务
            futures = [executor.submit(process_single_file, fastq_file, reference_genome, output_dir, threads) for fastq_file in fastq_files]
            # 收集结果
            for future in futures:
                result = future.result()
                if result:
                    bam_files.append(result)
    else:
        # 单线程处理
        for fastq_file in fastq_files:
            result = process_single_file(fastq_file, reference_genome, output_dir, threads)
            if result:
                bam_files.append(result)
    
    return bam_files

# python src/methylation/bam_generator.py /data/meth/bwa/hg38.fa /mnt/f/meth/
if __name__ == "__main__":
    # 对命令行参数输入的 文件/多层级目录下的文件 进行提取
    parser = argparse.ArgumentParser(description="批量生成FastQ文件对应的BAM文件")
    parser.add_argument("reference_genome", help="参考基因组文件路径")
    parser.add_argument("fastq_files", nargs="+", help="FastQ文件路径或目录路径")
    parser.add_argument("--parallel", type=int, default=1, help="并发处理数，默认为1")
    parser.add_argument("--threads", type=int, default=1, help="并发线程数，默认为1")
    args = parser.parse_args()

    file_list = []
    for fastq_file in args.fastq_files:
        if os.path.isdir(fastq_file):
            # 递归查找目录下的FastQ文件
            for root, dirs, files in os.walk(fastq_file):
                for file in files:
                    if file.endswith((".fq.gz", ".fastq.gz", ".fq", ".fastq")):
                        file_list.append(os.path.join(root, file))
        else:
            file_list.append(fastq_file)
    
    # 过滤出存在的文件
    file_list = [f for f in file_list if os.path.exists(f)]
    
    if not file_list:
        print("未找到有效的FastQ文件")
        exit(1)
    
    print(f"找到 {len(file_list)} 个FastQ文件")
    generate_bam(file_list, args.reference_genome, threads=args.threads)
    print("所有BAM文件生成完成")

