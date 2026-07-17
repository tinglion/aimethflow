#!/usr/bin/env python3
"""
file_utils.py - File utility module

Functions for file type detection and discovery.
"""

import os
import glob
import logging

logger = logging.getLogger(__name__)

def detect_file_type(file_path):
    """
    Detect the type of a file.

    Args:
        file_path: path to the file

    Returns:
        str: one of 'fastq', 'bam', 'methylation', 'unknown'
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
    Find all supported file types under the given paths.

    Args:
        paths: list of file or directory paths

    Returns:
        dict mapping type name to list of absolute file paths
    """
    files = {"fastq": [], "bam": [], "methylation": [], "unknown": []}

    for path in paths:
        if os.path.isfile(path):
            file_type = detect_file_type(path)
            files[file_type].append(os.path.abspath(path))
        elif os.path.isdir(path):
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

    # deduplicate and sort
    for file_type in files:
        files[file_type] = sorted(list(set(files[file_type])))

    return files
