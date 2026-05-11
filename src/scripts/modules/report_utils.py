"""report_utils.py — backward-compatibility shim.

The implementation moved to the `scripts.modules.report` package (May 2026)
to split the former 2,455-line monolith into focused sub-modules. Existing
callers that do ``from scripts.modules.report_utils import generate_report``
continue to work through this re-export.

New code should import directly from `scripts.modules.report`.
"""

from scripts.modules.report import generate_report  # noqa: F401

# Re-export private helpers used by external callers (e.g. one-off scripts
# that imported underscore-prefixed names directly).
from scripts.modules.report import (  # noqa: F401
    _parse_fastp_json, _parse_dedup_stats, _count_bedgraph_sites,
    _parse_diff_results,
    _fmt_num, _fmt_pct, _fmt_mb, _fmt_gb, _safe_pct,
    _copy_img, _img_tag, _read_csv_safe, _pagination_widget,
    _build_pipeline_progress_section,
    _build_outlier_section,
    _build_classifier_section,
    _build_ovr_subsection,
    _build_annotation_section,
    _build_external_validation_section,
)

__all__ = ["generate_report"]
