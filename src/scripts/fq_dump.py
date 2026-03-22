#!/usr/bin/env python3
"""
fq_inspect.py - 用于观察.fq.gz文件数据的脚本

该脚本可以读取并分析FastQ格式的压缩文件，显示文件的基本信息和序列数据。
支持分块加载，避免处理大文件时的内存问题。

使用方法:
    python fq_inspect.py <file_path> [--start <position>] [--count <number>]

参数:
    file_path: .fq.gz文件的路径
    --start: 可选，读取开始位置，默认为0
    --count: 可选，读取数量，默认为10
"""

import sys
import gzip
from Bio import SeqIO
import argparse


def inspect_fastq(file_path, start=0, count=10):
    """
    分析FastQ文件并显示相关信息
    
    Args:
        file_path: FastQ文件路径
        start: 读取开始位置
        count: 读取数量
    """
    try:
        # 打开压缩文件
        with gzip.open(file_path, 'rt') as handle:
            # 分块读取序列
            records = []
            total_records = 0
            
            for record in SeqIO.parse(handle, 'fastq'):
                total_records += 1
                # 跳过直到开始位置
                if total_records <= start:
                    continue
                # 收集指定数量的序列
                if len(records) < count:
                    records.append(record)
                else:
                    break
            
        # 显示基本信息
        print(f"文件: {file_path}")
        print(f"总序列数量: {total_records}")
        print(f"读取位置: 从 {start} 开始，读取 {len(records)} 条序列")
        
        if records:
            # 计算序列长度
            lengths = [len(record) for record in records]
            min_len = min(lengths)
            max_len = max(lengths)
            avg_len = sum(lengths) / len(lengths)
            
            print(f"序列长度范围: {min_len} - {max_len}")
            print(f"平均序列长度: {avg_len:.2f}")
            
            # 显示读取的序列信息
            print("\n读取的序列信息:")
            for i, record in enumerate(records):
                print(f"\n序列 {start + i + 1}:")
                print(f"ID: {record.id}")
                print(f"描述: {record.description}")
                print(f"序列: {record.seq}")
                print(f"质量分数: {record.letter_annotations['phred_quality']}")
                print(f"序列长度: {len(record)}")
        else:
            print("没有读取到序列")
            
    except FileNotFoundError:
        print(f"错误: 文件 '{file_path}' 不存在")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='分析FastQ文件')
    parser.add_argument('file_path', help='.fq.gz文件的路径')
    parser.add_argument('--start', type=int, default=0, 
                        help='读取开始位置，默认为0')
    parser.add_argument('--count', type=int, default=10, 
                        help='读取数量，默认为10')
    
    args = parser.parse_args()
    
    # 执行分析
    inspect_fastq(args.file_path, args.start, args.count)
