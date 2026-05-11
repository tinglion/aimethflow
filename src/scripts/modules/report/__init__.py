"""scripts.modules.report — HTML/JSON report generation package.

This package replaces the former 2,455-line `report_utils.py` module with
focused sub-modules:

    helpers     — parsers, formatters, image/CSV helpers, pagination widget
    progress    — Step 1–20 pipeline progress chip grid
    outlier     — Step 13 sample-quality / outlier-detection section
    classifier  — Steps 14–17 (feature selection, classifiers, subtype, temporal)
    annotation  — Step 18 panel-CpG site annotation section
    external    — Step 19 public-dataset cross-validation section
    main        — top-level `generate_report` entry point + final HTML assembly

The public entry point is `generate_report(...)`. Existing call sites that
import from `scripts.modules.report_utils` continue to work via the shim
left in place at that module path.
"""

from scripts.modules.report.main import generate_report

# Re-export selected private helpers so the shim's `from ... import *` keeps
# working for any code that imported underscore-prefixed names directly.
from scripts.modules.report.helpers import (  # noqa: F401
    _parse_fastp_json, _parse_dedup_stats, _count_bedgraph_sites,
    _parse_diff_results,
    _fmt_num, _fmt_pct, _fmt_mb, _fmt_gb, _safe_pct,
    _copy_img, _img_tag, _read_csv_safe, _pagination_widget,
)
from scripts.modules.report.progress    import _build_pipeline_progress_section  # noqa: F401
from scripts.modules.report.outlier     import _build_outlier_section            # noqa: F401
from scripts.modules.report.classifier  import (                                  # noqa: F401
    _build_classifier_section, _build_ovr_subsection,
)
from scripts.modules.report.annotation  import _build_annotation_section         # noqa: F401
from scripts.modules.report.external    import _build_external_validation_section  # noqa: F401

__all__ = ["generate_report"]
