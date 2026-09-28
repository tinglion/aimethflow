#!/usr/bin/env python3
"""
Command-line runner for T2 — re-run the AIMethFlow ML core on an external
cohort (``docs/GPB/02_external_validation_plan.md``).

This is a thin wrapper over ``scripts.modules.external_cohort``, which is also
what pipeline Step 19b (``--step external_cohort``) calls. Use this script to
run a single cohort ad hoc; use the pipeline step to run every cohort listed in
``config.toml [external_cohort] dirs`` and fold the results into the HTML
report. Both paths execute identical code.

See the module docstring for what is overridden and why.

Examples
--------
::

    # primary: framework re-run, full feature pool, in-fold selection
    python src/scripts/run_external_cohort.py --result-dir result/external/ovaprint

    # plus a de-novo panel derived on this cohort
    python src/scripts/run_external_cohort.py \\
        --result-dir result/external/ovaprint --discovery

    # the same work as a pipeline step, for every configured cohort
    run-analysis --mode plasma --step external_cohort --workdir <workdir>
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from utils.logger import logger                                  # noqa: E402
from scripts.modules.external_cohort import run_one_cohort       # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Re-run the AIMethFlow ML core on an external cohort (T2).")
    ap.add_argument("--result-dir", required=True,
                    help="directory containing feature_matrix.npz, e.g. "
                         "result/external/ovaprint; outputs are written alongside it")
    ap.add_argument("--workdir", default=None,
                    help="working directory, used only if the matrix has to be "
                         "rebuilt from bedGraph (default: --result-dir)")
    ap.add_argument("--cohort-label", default=None,
                    help="label used in the log and summary rows "
                         "(default: basename of --result-dir)")
    ap.add_argument("--cv-folds", type=int, default=None,
                    help="override [ml] cv_folds; 0 or less selects LOO")
    ap.add_argument("--discovery", action="store_true",
                    help="also derive a panel de novo on this cohort (Step 14) in "
                         "a panel_discovery/ subdirectory and re-run Step 15 in it")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.result_dir):
        logger.error(f"result dir not found: {args.result_dir}")
        return 1

    rows = run_one_cohort(
        args.result_dir,
        workdir=args.workdir,
        cv_folds=args.cv_folds,
        discovery=args.discovery,
        label=args.cohort_label,
    )
    return 0 if rows else 1


if __name__ == "__main__":
    sys.exit(main())
