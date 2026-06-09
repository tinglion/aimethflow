#!/usr/bin/env python3
"""
modules package init
"""

from .file_utils import detect_file_type, find_files
from .fastq_utils import analyze_quality, get_fastq_stats, preprocess_fastq
from .alignment_utils import (align_methylation, deduplicate_bam, extract_target_regions,
                              check_mbias, extract_methylation_methyldackel,
                              plot_fragment_size)
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

# NOTE: `feature_utils`, `classifier_utils`, and (sibling) `external_validate`
# are intentionally NOT re-exported here. They live in the private companion
# package `aimethflow-core` during the pre-publication phase. Callers that
# need them must import the full path (e.g. `from scripts.modules.feature_utils
# import build_feature_matrix`), which triggers the stub shim's transparent
# substitution. Top-level re-exports here would force every user of the
# `scripts.modules` package to depend on `aimethflow-core`, breaking the
# public-only access path for Steps 1-13, 16-18, 20.

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
    "plot_fragment_size",
    "extract_methylation_info",
    "validate_methylation_output",
    "generate_mock_methylation_data",
    "merge_methylation_files",
    "filter_methylation_data",
    "filter_to_candidate_regions_wrapper",
    "perform_diff_analysis",
    "visualize_results",
    "generate_report",
]
