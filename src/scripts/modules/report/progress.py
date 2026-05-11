#!/usr/bin/env python3
"""report.progress — Pipeline Progress chip grid (Step 1–20 status panel)."""

import os
import json
import logging

logger = logging.getLogger(__name__)


def _build_pipeline_progress_section(temp_dir):
    """
    Render a compact 22-step progress grid based on pipeline_state.json/completed_steps.
    Returns an HTML <div class="section">...</div> block.
    """
    # Match PIPELINE_STEPS in run_analysis.py — keep numbers/labels in sync.
    STEPS = [
        ("data_check",        "1",   "Validate config"),
        ("ref_index",         "2",   "bwameth hg38 index"),
        ("panel_derive",      "3",   "Derive panel BED"),
        ("preprocess",        "4",   "fastp QC + UMI"),
        ("align",             "5",   "bwameth align"),
        ("dedup",             "6",   "UMI dedup"),
        ("target",            "7",   "On-target filter"),
        ("mbias",             "8",   "M-bias check"),
        ("extract",           "9",   "Methylation extract"),
        ("filter",            "10",  "Depth + repeat filter"),
        ("diff",              "11",  "Differential methylation"),
        ("visualize",         "12",  "Plots"),
        ("outlier_detect",    "13",  "Outlier detection"),
        ("feature_select",    "14",  "Feature selection"),
        ("classify",          "15",  "Classifiers (LOO/CV)"),
        ("fs_benchmark",      "15b", "FS ablation"),
        ("ca125_compare",     "15c", "CA-125 vs cfDNA"),
        ("subtype_classify",  "16",  "Subtype classifier"),
        ("temporal_validate", "17",  "Temporal split val."),
        ("annotate",          "18",  "Site annotation"),
        ("external_validate", "19",  "External validation"),
        ("report",            "20",  "Report"),
    ]

    state_path = os.path.join(temp_dir, "pipeline_state.json")
    completed = set()
    if os.path.exists(state_path):
        try:
            with open(state_path, encoding="utf-8") as f:
                completed = set(json.load(f).get("completed_steps", []))
        except Exception:
            pass

    n_done = sum(1 for k, *_ in STEPS if k in completed)
    chips = []
    for key, num, label in STEPS:
        is_done = key in completed
        bg   = "#dcfce7" if is_done else "#f3f4f6"
        col  = "#166534" if is_done else "#6b7280"
        bord = "#86efac" if is_done else "#e5e7eb"
        mark = "&#10003;" if is_done else "&middot;"
        chips.append(
            f"<div title='{key}' style='display:inline-block;background:{bg};color:{col};"
            f"border:1px solid {bord};border-radius:6px;padding:4px 8px;margin:2px;"
            f"font-size:.8em;min-width:96px'>"
            f"<span style='font-weight:bold;margin-right:4px'>{mark}</span>"
            f"<span style='font-family:monospace'>{num}</span> "
            f"<span style='color:{col}'>{label}</span></div>"
        )

    html = (
        '<div class="section">'
        '<h2 style="margin-top:0">Pipeline Progress</h2>'
        f'<p style="color:#666;font-size:.85em;margin:0 0 8px">'
        f'<strong>{n_done} / {len(STEPS)}</strong> steps complete '
        f'(state file: <code>{state_path}</code>). '
        f'Green = done; grey = not yet run. '
        f'Hover for the step key.</p>'
        '<div style="line-height:1.8">'
        + "".join(chips) +
        '</div></div>'
    )
    return html

