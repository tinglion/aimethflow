#!/usr/bin/env python3
"""Region-level Δβ direction concordance for the committed panel (E4).

Tests the *panel* rather than the framework: the committed panel's hg38 sites are
lifted to hg19, matched against the Illumina EPIC manifest, and their plasma
malignant-vs-benign direction is compared with the tumour-vs-normal direction in an
independent EPIC tissue cohort.

This is direction concordance, **not** panel transfer — see
``scripts.modules.panel_concordance`` for the full statement of what the result
does and does not support.

Usage:

    python src/scripts/run_panel_concordance.py \\
        --result-dir result/20260406 \\
        --series-matrix data/external/GSE133556/GSE133556_series_matrix.txt.gz

Defaults resolve every input from the repository layout, so the bare form works
after ``download_external_data.sh`` has run.
"""

import os
import sys
import argparse
import logging

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from scripts.modules.panel_concordance import (   # noqa: E402
    run_panel_concordance,
    _find_epic_manifest,
    MATCH_WINDOWS,
)


def main():
    ap = argparse.ArgumentParser(
        description="Δβ direction concordance between the committed panel's plasma "
                    "effects and an independent EPIC tissue cohort.")
    ap.add_argument("--result-dir", default="result/20260406",
                    help="Directory holding panel_cpg.csv and wilcoxon_stats.csv")
    ap.add_argument("--external-dir", default="data/external",
                    help="Directory holding the liftOver chain and EPIC manifest")
    ap.add_argument("--series-matrix",
                    default="data/external/GSE133556/GSE133556_series_matrix.txt.gz",
                    help="GEO series matrix for the EPIC tissue cohort")
    ap.add_argument("--case-label", default="High grade serous ovarian cancer",
                    help="Value of the 'tissue:' characteristic marking cases")
    ap.add_argument("--control-label", default="Normal Fallopian tube control",
                    help="Value of the 'tissue:' characteristic marking controls")
    ap.add_argument("--out", default=None,
                    help="Output directory (default: <result-dir>/panel_concordance)")
    ap.add_argument("--windows", default=",".join(str(w) for w in MATCH_WINDOWS),
                    help="Comma-separated matching windows in bp")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    log = logging.getLogger(__name__)

    panel_csv = os.path.join(args.result_dir, "panel_cpg.csv")
    wilcox_csv = os.path.join(args.result_dir, "wilcoxon_stats.csv")
    chain_gz = os.path.join(args.external_dir, "hg38ToHg19.over.chain.gz")
    out_dir = args.out or os.path.join(args.result_dir, "panel_concordance")

    missing = [p for p in (panel_csv, wilcox_csv, chain_gz, args.series_matrix)
               if not os.path.exists(p)]
    if missing:
        for p in missing:
            log.error(f"Missing input: {p}")
        log.error("See download_external_data.sh for the external files.")
        return 1

    epic = _find_epic_manifest(args.external_dir)
    if not epic:
        log.error(f"No EPIC manifest found in {args.external_dir}/.")
        log.error("Download the GPL21145 platform table:")
        log.error("  curl -L -o data/external/GPL21145_platform.txt \\")
        log.error("    'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
                  "?targ=self&view=data&form=text&acc=GPL21145'")
        return 1

    windows = [int(w) for w in args.windows.split(",") if w.strip()]

    run_panel_concordance(
        panel_csv=panel_csv,
        wilcoxon_csv=wilcox_csv,
        series_matrix=args.series_matrix,
        chain_gz=chain_gz,
        epic_manifest=epic,
        out_dir=out_dir,
        case_label=args.case_label,
        control_label=args.control_label,
        windows=windows,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
