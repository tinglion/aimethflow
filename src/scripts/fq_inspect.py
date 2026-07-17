#!/usr/bin/env python3
"""
fq_inspect.py - 快速查看.fq.gz文件统计信息的脚本

该脚本可以快速分析FastQ格式的压缩文件，显示文件的统计信息。

使用方法:
    python fq_inspect.py <file_path>

参数:
    file_path: .fq.gz文件的路径
"""

import sys
import gzip
from Bio import SeqIO
import argparse
import time


def inspect_fastq(file_path):
    """
    快速分析FastQ文件并显示统计信息

    Args:
        file_path: FastQ文件路径
    """
    try:
        start_time = time.time()

        # 打开压缩文件
        with gzip.open(file_path, "rt") as handle:
            # 初始化统计变量
            total_records = 0
            total_bases = 0
            lengths = []
            quality_scores = []

            # 分块读取序列
            for record in SeqIO.parse(handle, "fastq"):
                total_records += 1
                length = len(record)
                total_bases += length
                lengths.append(length)
                if total_records % 1000000 ==0:
                    print(
                        f"read record={total_records} time={time.time() - start_time:.2f}"
                    )

                # 计算平均质量分数
                if record.letter_annotations.get("phred_quality"):
                    avg_quality = (
                        sum(record.letter_annotations["phred_quality"]) / length
                    )
                    quality_scores.append(avg_quality)

        # 显示统计信息
        print(f"文件: {file_path} time={time.time() - start_time:.2f}")
        print(f"总序列数量: {total_records}")
        print(f"总碱基数: {total_bases}")

        if lengths:
            # 计算长度统计
            min_len = min(lengths)
            max_len = max(lengths)
            avg_len = total_bases / total_records

            print(f"序列长度范围: {min_len} - {max_len}")
            print(f"平均序列长度: {avg_len:.2f}")

            # 计算质量分数统计
            if quality_scores:
                avg_quality = sum(quality_scores) / len(quality_scores)
                print(f"平均质量分数: {avg_quality:.2f}")
        else:
            print("文件中没有序列")

    except FileNotFoundError:
        print(f"错误: 文件 '{file_path}' 不存在")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="快速查看FastQ文件统计信息")
    parser.add_argument("file_path", help=".fq.gz文件的路径")

    args = parser.parse_args()

    # 执行分析
    inspect_fastq(args.file_path)
