#!/usr/bin/env python3
"""
modules package init
"""

from .file_utils import detect_file_type, find_files
from .fastq_utils import analyze_quality, get_fastq_stats, preprocess_fastq
from .alignment_utils import (align_methylation, deduplicate_bam, extract_target_regions,
                              check_mbias, extract_methylation_methyldackel)
from .methylation_utils import (
    extract_methylation_info,
    validate_methylation_output,
    generate_mock_methylation_data,
    merge_methylation_files,
    filter_methylation_data,
    filter_to_candidate_regions_wrapper,
    perform_diff_analysis,
    visualize_results
)
from .report_utils import generate_report
from .feature_utils import build_feature_matrix, select_features
from .classifier_utils import run_classifiers

__all__ = [
    "detect_file_type",
    "find_files",
    "analyze_quality",
    "get_fastq_stats",
    "preprocess_fastq",
    "align_methylation",
    "deduplicate_bam",
    "extract_target_regions",
    "check_mbias",
    "extract_methylation_methyldackel",
    "extract_methylation_info",
    "validate_methylation_output",
    "generate_mock_methylation_data",
    "merge_methylation_files",
    "filter_methylation_data",
    "filter_to_candidate_regions_wrapper",
    "perform_diff_analysis",
    "visualize_results",
    "generate_report",
    "build_feature_matrix",
    "select_features",
    "run_classifiers",
]
