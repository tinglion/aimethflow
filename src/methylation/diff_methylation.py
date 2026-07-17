#!/usr/bin/env python3
"""
diff_methylation.py - 差异甲基化分析脚本

该脚本用于比较两个样本的甲基化数据，找出差异甲基化位点。

使用方法:
    python diff_methylation.py <sample1_file> <sample2_file> <output_file> [--threshold <value>]

参数:
    sample1_file: 第一个样本的甲基化数据文件
    sample2_file: 第二个样本的甲基化数据文件
    output_file: 差异甲基化分析结果输出文件
    --threshold: 甲基化差异阈值，默认为0.2
"""

import sys
import argparse
import os
import time

from utils.logger import logger


def load_methylation_data(file_path):
    """
    加载甲基化数据
    
    Args:
        file_path: 甲基化数据文件路径
    
    Returns:
        dict: 甲基化数据字典，键为(染色体, 位置)，值为甲基化水平
    """
    data = {}
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or line.startswith('track'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                pos = int(parts[1])
                # 处理不同格式的文件
                if len(parts) == 4:
                    # 4列格式：chrom, start, end, methylation_percentage
                    meth_level = float(parts[3]) / 100  # 转换为0-1范围
                elif len(parts) == 6:
                    # 6列格式：chrom, start, end, methylation_percentage, meth_count, unmeth_count
                    try:
                        meth_count = int(parts[4])
                        unmeth_count = int(parts[5])
                        total_count = meth_count + unmeth_count
                        if total_count > 0:
                            meth_level = meth_count / total_count
                        else:
                            meth_level = 0.0
                    except ValueError:
                        # 如果解析失败，尝试使用第4列作为甲基化水平
                        meth_level = float(parts[3]) / 100
                else:
                    # 其他格式，使用第4列作为甲基化水平
                    meth_level = float(parts[3])
                data[(chrom, pos)] = meth_level
    
    return data


def analyze_diff_methylation(sample1_data, sample2_data, threshold=0.2):
    """
    分析差异甲基化位点
    
    Args:
        sample1_data: 第一个样本的甲基化数据
        sample2_data: 第二个样本的甲基化数据
        threshold: 甲基化差异阈值
    
    Returns:
        list: 差异甲基化位点列表
    """
    diff_sites = []
    
    # 获取所有共同的位点
    common_sites = set(sample1_data.keys()) & set(sample2_data.keys())
    logger.info(f"共同的位点数量: {len(common_sites)}")
    
    for site in common_sites:
        meth1 = sample1_data[site]
        meth2 = sample2_data[site]
        diff = abs(meth1 - meth2)
        
        if diff >= threshold:
            diff_sites.append((site[0], site[1], meth1, meth2, diff))
    
    return diff_sites


def run_diff_analysis(sample1_file, sample2_file, output_file, threshold=0.2):
    """
    运行差异甲基化分析
    
    Args:
        sample1_file: 第一个样本的甲基化数据文件
        sample2_file: 第二个样本的甲基化数据文件
        output_file: 输出文件路径
        threshold: 甲基化差异阈值
    """
    start_time = time.time()
    try:
        # 检查输入文件是否存在
        if not os.path.exists(sample1_file):
            print(f"错误: 样本1文件 '{sample1_file}' 不存在")
            sys.exit(1)
        
        if not os.path.exists(sample2_file):
            print(f"错误: 样本2文件 '{sample2_file}' 不存在")
            sys.exit(1)
        
        print(f"加载样本1数据: {sample1_file}")
        sample1_data = load_methylation_data(sample1_file)
        
        print(f"加载样本2数据: {sample2_file}")
        sample2_data = load_methylation_data(sample2_file)
        
        print("分析差异甲基化位点...")
        diff_sites = analyze_diff_methylation(sample1_data, sample2_data, threshold)
        
        print(f"找到 {len(diff_sites)} 个差异甲基化位点")
        
        # 输出结果
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("# 染色体\t位置\t样本1甲基化水平\t样本2甲基化水平\t差异\n")
            for site in diff_sites:
                chrom, pos, meth1, meth2, diff = site
                f.write(f"{chrom}\t{pos}\t{meth1}\t{meth2}\t{diff}\n")
        
        print(f"结果已输出到: {output_file} time={time.time() - start_time:.2f}")
        
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='差异甲基化分析')
    parser.add_argument('sample1_file', help='第一个样本的甲基化数据文件')
    parser.add_argument('sample2_file', help='第二个样本的甲基化数据文件')
    parser.add_argument('output_file', help='差异甲基化分析结果输出文件')
    parser.add_argument('--threshold', type=float, default=0.2, 
                        help='甲基化差异阈值，默认为0.2')
    
    args = parser.parse_args()
    
    # 执行分析
    run_diff_analysis(args.sample1_file, args.sample2_file, args.output_file, args.threshold)
