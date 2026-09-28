#!/usr/bin/env python3
"""
Command-line runner for the feature-selection leakage demonstration
(``scripts.modules.leakage_demo``).

Operates on any result directory that already contains a ``feature_matrix.npz``
— the internal cohort produced by Steps 1-13, or an external cohort produced by
``scripts/prepare_external_matrix.py``. Because the matrix is the entry point,
no FASTQ, alignment or annotation step is required.

Examples
--------
::

    # internal demonstration cohort
    python src/scripts/run_leakage_demo.py --result-dir result/20260406

    # external cohort (OvaPrint plasma cfDNA)
    python src/scripts/run_leakage_demo.py \\
        --result-dir result/external/ovaprint \\
        --cohort-label ovaprint --no-exclude-samples
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from utils.config import load_config          # noqa: E402
from utils.logger import logger               # noqa: E402
from scripts.modules.leakage_demo import run_leakage_demo   # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Quantify AUC inflation from leaky feature selection.")
    ap.add_argument("--result-dir", required=True,
                    help="directory containing feature_matrix.npz; outputs are "
                         "written alongside it")
    ap.add_argument("--workdir", default=None,
                    help="working directory used only if the matrix has to be "
                         "rebuilt from bedGraph (default: --result-dir)")
    ap.add_argument("--cohort-label", default=None,
                    help="label used in output tables and the figure title "
                         "(default: basename of --result-dir)")
    ap.add_argument("--cv-folds", type=int, default=None,
                    help="override [ml] cv_folds; 0 or less selects LOO")
    ap.add_argument("--n-boot", type=int, default=1000,
                    help="bootstrap resamples for the AUC CI (default: 1000)")
    ap.add_argument("--no-exclude-samples", action="store_true",
                    help="ignore [ml] exclude_samples — those IDs belong to the "
                         "internal cohort and must not be applied to external data")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.result_dir):
        logger.error(f"result dir not found: {args.result_dir}")
        return 1

    config = load_config()
    config.setdefault("ml", {})
    if args.cv_folds is not None:
        config["ml"]["cv_folds"] = args.cv_folds
    if args.no_exclude_samples:
        config["ml"]["exclude_samples"] = []

    df, ddf = run_leakage_demo(
        workdir=args.workdir or args.result_dir,
        result_dir=args.result_dir,
        config=config,
        cohort_label=args.cohort_label,
        n_boot=args.n_boot,
    )
    return 0 if not df.empty else 1


if __name__ == "__main__":
    sys.exit(main())
