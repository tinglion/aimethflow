#!/usr/bin/env python3
"""
file_utils.py - 文件工具模块

该模块包含文件检测和查找相关的函数
"""

import os
import glob
import logging

logger = logging.getLogger(__name__)

def detect_file_type(file_path):
    """
    检测文件类型

    Args:
        file_path: 文件路径

    Returns:
        str: 文件类型，可能的值：fastq, bam, methylation, unknown
    """
    if (
        file_path.endswith(".fq.gz")
        or file_path.endswith(".fastq.gz")
        or file_path.endswith(".fq")
        or file_path.endswith(".fastq")
    ):
        return "fastq"
    elif file_path.endswith(".bam"):
        return "bam"
    elif (
        file_path.endswith(".txt")
        or file_path.endswith(".bedGraph")
        or file_path.endswith(".cov")
        or "methylation" in file_path
    ):
        return "methylation"
    else:
        return "unknown"


def find_files(paths):
    """
    查找所有支持的文件类型

    Args:
        paths: 文件或文件夹路径列表

    Returns:
        dict: 包含不同类型文件的字典
    """
    files = {"fastq": [], "bam": [], "methylation": [], "unknown": []}

    for path in paths:
        if os.path.isfile(path):
            # 如果是文件，检测文件类型
            file_type = detect_file_type(path)
            files[file_type].append(os.path.abspath(path))
        elif os.path.isdir(path):
            # 如果是文件夹，递归查找所有支持的文件
            for ext in [
                "*.fq.gz",
                "*.fastq.gz",
                "*.fq",
                "*.fastq",
                "*.bam",
                "*.txt",
                "*.bedGraph",
                "*.cov",
            ]:
                found_files = glob.glob(os.path.join(path, "**", ext), recursive=True)
                for file_path in found_files:
                    file_type = detect_file_type(file_path)
                    files[file_type].append(os.path.abspath(file_path))

    # 去重并排序
    for file_type in files:
        files[file_type] = sorted(list(set(files[file_type])))

    return files
