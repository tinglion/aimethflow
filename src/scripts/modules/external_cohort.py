#!/usr/bin/env python3
"""
external_cohort.py — Step 19b: re-run the AIMethFlow ML core on one or more
independent external cohorts.

Step 19 (``external_validate``) cross-references the panel against public
*reference* datasets (GSE40279 whole blood, TCGA-OV). It is a panel-level
provenance check, not an independent cohort. This module is the independent
cohort: it invokes the **production** Step 15 routine
(``classifier_utils.run_classifiers``) unmodified on an external feature matrix,
so the numbers it produces are directly comparable to the internal cohort's
``classifier_summary.csv``.

What it validates is the *framework*, not the 60-CpG panel. No public cohort
shares our probe set, so the panel itself cannot be transferred; see
``docs/GPB/02_external_validation_plan.md`` §0.

Three cohort-specific settings have to be overridden, and only these three:

``[ml] include_ca125_baseline``
    Forced off. CA-125 is read from the internal sample Excel; external cohorts
    have no CA-125 column, so the fusion arm is simply not run and the benchmark
    reduces to the five cfDNA-only classifiers.

``[ml] exclude_samples``
    Cleared. Those sample IDs were resolved by Step 13 against the internal
    cohort and would silently match nothing — or, worse, coincidentally match —
    in an external one.

Feature space
    Step 15 restricts the matrix to ``panel_cpg*.csv`` when one is present in
    the result directory. An external result directory produced by
    ``prepare_external_matrix.py`` has none, so the run uses the full matrix
    with selection confined to each training fold — the leakage-free protocol,
    no panel prior.

Extra arms
    ``[[external_cohort.arms]]`` attaches further result directories to a
    cohort. Each one is a second build of the *same* deposit made by
    ``prepare_external_matrix.py`` with a different sample selection; Step 15 is
    re-run in it and the rows are reported beside the primary arm under their
    own ``arm`` name. Two are defined for OvaPrint:

    ``hgsoc_matched`` (``--subset hgsoc_matched``)
        The original publication's own model set — 82 HGSOC versus 100 benign
        masses. The primary arm keeps every malignant histotype, which is the
        harder all-comers question; this arm exists so the reproduction can be
        checked against the published figures on the definition those figures
        were computed under.

    ``with_healthy_donors`` (``--include-normal``)
        The deposit's 26 healthy donors folded into the benign class. Excluded
        from the primary arm because every AIMethFlow participant presented with
        an adnexal mass, so easy negatives would raise AUC for a reason
        unrelated to the method.

    Both are sensitivity arms and neither is ever the headline number.

The manifest's ``batch`` column is carried through from the deposit's
provenance field but is **deliberately not analysed**. The OvaPrint deposit
ships that field already anonymised (``SOURCE_1`` … ``SOURCE_8``), so nothing
can be established about what distinguishes the groups, and any within-group
result would invite an interpretation the data cannot support.

``discovery`` additionally exercises the *discovery* half of the pipeline:
Step 14 (``select_features``) is run on a copy of the matrix in a
``panel_discovery/`` subdirectory to derive a panel de novo from that cohort,
and Step 15 is then re-run inside it. Comparing the two runs reproduces, on
independent data, the panel-definition effect quantified in
``docs/GPB/03_leakage_demonstration.md`` §4. The copy exists so the primary
result directory keeps its unrestricted feature space and its original
``feature_matrix.npz`` (Step 14 rewrites that file, dropping borderline rows).

Outputs
-------
Per external cohort, in that cohort's own directory:

- ``classifier_summary.csv``, ``classifier_predictions.csv``, ``roc_curve.png``, …
  — the ordinary Step 15 artefacts.
- ``panel_discovery/`` — the de-novo panel run, when ``discovery`` is enabled.

In the **main** result directory:

- ``external_cohort_summary.csv`` — every cohort × arm × classifier in one
  table, with the internal cohort's own rows included so the HTML report can
  render them side by side.
"""

import os
import shutil

import pandas as pd

from utils.logger import logger

# The five cfDNA-only classifiers. The internal cohort also carries CA-125 and
# the *+CA125 fusion models; no external cohort can, so the side-by-side table
# is restricted to the models both cohorts actually have.
CORE_CLASSIFIERS = ["LR", "SVM", "RF", "XGB", "LGBM"]

SUMMARY_FILE = "external_cohort_summary.csv"

# Column order of external_cohort_summary.csv.
SUMMARY_COLUMNS = [
    "cohort", "arm", "classifier", "n_samples", "n_features_mean",
    "auc", "auc_ci_low", "auc_ci_high", "accuracy", "sensitivity",
    "specificity", "delta_auc_vs_full", "delong_p",
]


def external_config(config=None, cv_folds=None):
    """Return a copy of *config* with every internal-cohort-specific setting neutralised.

    Pass the pipeline's own config when calling from Step 19b; omit it to load
    config.toml, which is what the standalone CLI does.
    """
    if config is None:
        from utils.config import load_config
        config = load_config()
    # Shallow copy plus a fresh [ml] dict: the caller's config must not be
    # mutated, because the pipeline keeps using it after this step.
    config = dict(config)
    ml = dict(config.get("ml", {}))
    # CA-125 is read from the internal sample Excel and keyed by internal sample
    # IDs, so no external sample can resolve to a value. GSE19711 does deposit
    # CA-125, but there is no path from it into this code, and its screening-era
    # values are not the preoperative triage measurement the internal cohort has.
    ml["include_ca125_baseline"] = False
    # These IDs belong to the internal cohort only.
    ml["exclude_samples"] = []
    if cv_folds is not None:
        ml["cv_folds"] = cv_folds
    config["ml"] = ml
    return config


def _run_step15(result_dir, workdir, config, banner):
    """Run the production Step 15 and return its summary DataFrame."""
    from scripts.modules.classifier_utils import run_classifiers

    logger.info(f"━━ {banner} ━━")
    _, summary_df = run_classifiers(workdir, result_dir, config, [], [])
    if summary_df is None or summary_df.empty:
        logger.error(f"  no classifier results produced in {result_dir}")
        return None
    logger.info(f"  → {os.path.join(result_dir, 'classifier_summary.csv')}")
    return summary_df


def _run_discovery(result_dir, workdir, config):
    """Derive a panel de novo on this cohort, then re-run Step 15 inside it.

    Runs in a subdirectory so the primary result directory keeps its
    unrestricted feature space: Step 14 writes panel_cpg.csv (which Step 15
    would then auto-apply) and rewrites feature_matrix.npz.
    """
    from scripts.modules.feature_utils import select_features

    disc_dir = os.path.join(result_dir, "panel_discovery")
    os.makedirs(disc_dir, exist_ok=True)
    src_npz = os.path.join(result_dir, "feature_matrix.npz")
    dst_npz = os.path.join(disc_dir, "feature_matrix.npz")
    if not os.path.exists(src_npz):
        logger.error(f"  feature_matrix.npz not found in {result_dir}")
        return None
    shutil.copy2(src_npz, dst_npz)

    logger.info("━━ Step 14 on the external cohort: de-novo panel definition ━━")
    panel_features, _ = select_features(workdir, disc_dir, config, [], [])
    if not panel_features:
        logger.error("  Step 14 selected no features — discovery run aborted")
        return None
    logger.info(f"  → {len(panel_features)} panel features → "
                f"{os.path.join(disc_dir, 'panel_cpg.csv')}")

    return _run_step15(disc_dir, workdir, config,
                       "Step 15 restricted to the de-novo panel")


def _paired_probs(result_dir):
    """{classifier: {sample_id: (label, prob)}} from classifier_predictions.csv."""
    path = os.path.join(result_dir, "classifier_predictions.csv")
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path)
    out = {}
    for clf, sub in df.groupby("classifier"):
        out[clf] = {
            str(r.sample_id): (1 if r.true_label == "malignant" else 0,
                               float(r.prob_malignant))
            for r in sub.itertuples()
        }
    return out


def compare_panel_effect(primary_dir, panel_dir):
    """Unrestricted vs panel-restricted AUC, paired sample-by-sample (§4 of the T1 doc).

    Both runs score the same samples over the same fold partition
    (StratifiedKFold depends only on y and n), so the AUCs are paired and the
    difference is testable with DeLong rather than only described.

    Logs the comparison table and returns
    ``{classifier: {"auc_full", "auc_panel", "delta_auc", "delong_p"}}``.
    """
    import numpy as np
    from scipy import stats as _stats
    from scripts.modules.pub_figures import _fastdelong_cov

    a_all, b_all = _paired_probs(primary_dir), _paired_probs(panel_dir)
    # groupby returns classifiers alphabetically; restore the benchmark's order.
    shared_clf = sorted(
        (c for c in a_all if c in b_all),
        key=lambda c: (CORE_CLASSIFIERS.index(c) if c in CORE_CLASSIFIERS
                       else len(CORE_CLASSIFIERS)))
    if not shared_clf:
        logger.warning("  comparison skipped — no per-sample predictions to pair")
        return {}

    out = {}
    logger.info("\n  Effect of a cohort-defined panel on this external cohort")
    logger.info("  ┌──────────┬───────────────┬───────────────┬───────────┬──────────┐")
    logger.info("  │ Clf      │ full matrix   │ de-novo panel │ ΔAUC      │ DeLong p │")
    logger.info("  ├──────────┼───────────────┼───────────────┼───────────┼──────────┤")
    for clf in shared_clf:
        a, b = a_all[clf], b_all[clf]
        sids = [s for s in a if s in b]
        if len(sids) < 5:
            continue
        y = np.array([a[s][0] for s in sids])
        pa = np.array([a[s][1] for s in sids])
        pb = np.array([b[s][1] for s in sids])
        if len(np.unique(y)) < 2:
            continue
        aucs, cov = _fastdelong_cov(y, np.vstack([pb, pa]))   # panel first
        var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
        if var > 0:
            z = (aucs[0] - aucs[1]) / np.sqrt(var)
            p_val = 2.0 * (1.0 - _stats.norm.cdf(abs(z)))
            p_txt = f"{p_val:.2e}"
        else:
            p_val, p_txt = float("nan"), "n/a"
        out[clf] = {
            "auc_full": float(aucs[1]),
            "auc_panel": float(aucs[0]),
            "delta_auc": float(aucs[0] - aucs[1]),
            "delong_p": p_val,
        }
        logger.info(f"  │ {clf:<8} │ {aucs[1]:>13.4f} │ {aucs[0]:>13.4f} │ "
                    f"{aucs[0] - aucs[1]:>+9.4f} │ {p_txt:>8} │")
    logger.info("  └──────────┴───────────────┴───────────────┴───────────┴──────────┘")
    return out


def _summary_rows(result_dir, cohort, arm, core_only=True):
    """Translate a classifier_summary.csv into external_cohort_summary.csv rows."""
    path = os.path.join(result_dir, "classifier_summary.csv")
    if not os.path.exists(path):
        return []
    df = pd.read_csv(path)
    rows = []
    for r in df.itertuples():
        clf = str(r.classifier)
        if core_only and clf not in CORE_CLASSIFIERS:
            continue
        rows.append({
            "cohort": cohort,
            "arm": arm,
            "classifier": clf,
            "n_samples": getattr(r, "n_used", None),
            "n_features_mean": getattr(r, "n_features_mean", None),
            "auc": getattr(r, "auc", None),
            "auc_ci_low": getattr(r, "auc_ci95_lo", None),
            "auc_ci_high": getattr(r, "auc_ci95_hi", None),
            "accuracy": getattr(r, "accuracy", None),
            "sensitivity": getattr(r, "sensitivity", None),
            "specificity": getattr(r, "specificity", None),
            "delta_auc_vs_full": None,
            "delong_p": None,
        })
    # Preserve the benchmark's classifier order rather than the file's.
    rows.sort(key=lambda d: (CORE_CLASSIFIERS.index(d["classifier"])
                             if d["classifier"] in CORE_CLASSIFIERS
                             else len(CORE_CLASSIFIERS)))
    return rows


def run_one_cohort(cohort_dir, workdir=None, config=None, cv_folds=None,
                   discovery=False, label=None):
    """Re-run Step 15 (and optionally a de-novo Step 14) on one external cohort.

    Returns the list of summary rows for this cohort, or an empty list if the
    run could not proceed.
    """
    if not os.path.isdir(cohort_dir):
        logger.error(f"  external cohort dir not found: {cohort_dir}")
        return []
    if not os.path.exists(os.path.join(cohort_dir, "feature_matrix.npz")):
        logger.error(f"  feature_matrix.npz not found in {cohort_dir} — run "
                     "src/scripts/prepare_external_matrix.py first")
        return []

    label = label or os.path.basename(os.path.normpath(cohort_dir))
    workdir = workdir or cohort_dir
    config = external_config(config, cv_folds)

    primary_df = _run_step15(
        cohort_dir, workdir, config,
        f"Step 15 on {label}: full feature pool, in-fold selection")
    if primary_df is None:
        return []

    rows = _summary_rows(cohort_dir, label, "full_matrix")

    if discovery:
        panel_df = _run_discovery(cohort_dir, workdir, config)
        if panel_df is not None:
            disc_dir = os.path.join(cohort_dir, "panel_discovery")
            panel_rows = _summary_rows(disc_dir, label, "denovo_panel")
            stats = compare_panel_effect(cohort_dir, disc_dir)
            for row in panel_rows:
                st = stats.get(row["classifier"])
                if st:
                    row["delta_auc_vs_full"] = round(st["delta_auc"], 4)
                    row["delong_p"] = st["delong_p"]
            rows.extend(panel_rows)

    return rows


def run_extra_arm(arm_dir, cohort, arm, config=None, cv_folds=None):
    """Re-run Step 15 on an alternative build of an already-configured cohort.

    An extra arm is the same deposit prepared by ``prepare_external_matrix.py``
    under a different sample selection — a published sub-cohort definition, or
    healthy donors folded in. It answers a different question from the primary
    arm rather than a better-posed version of the same one, so it is reported
    under its own ``arm`` name and is never the headline number.
    """
    if not os.path.isdir(arm_dir):
        logger.error(f"  arm dir not found: {arm_dir}")
        return []
    if not os.path.exists(os.path.join(arm_dir, "feature_matrix.npz")):
        logger.error(f"  feature_matrix.npz not found in {arm_dir} — build it with "
                     "prepare_external_matrix.py (--subset / --include-normal)")
        return []

    config = external_config(config, cv_folds)
    summary_df = _run_step15(
        arm_dir, arm_dir, config,
        f"Step 15 on {cohort}: {arm} (sensitivity arm)")
    if summary_df is None:
        return []
    return _summary_rows(arm_dir, cohort, arm)


def run_external_cohorts(workdir, result_dir, config):
    """Step 19b entry point — run every cohort in ``[external_cohort] dirs``.

    Writes ``external_cohort_summary.csv`` into *result_dir*, including the
    internal cohort's own rows so the report can render a side-by-side table.
    Returns the summary DataFrame (empty when nothing ran).
    """
    ext_cfg = (config or {}).get("external_cohort", {}) or {}
    dirs = list(ext_cfg.get("dirs", []) or [])
    labels = list(ext_cfg.get("labels", []) or [])
    extra_arms = list(ext_cfg.get("arms", []) or [])
    discovery = bool(ext_cfg.get("discovery", False))

    if not dirs:
        logger.warning(
            "  [external_cohort] dirs is empty — nothing to run. Prepare a cohort "
            "with src/scripts/prepare_external_matrix.py, then list its result "
            "directory in config.toml [external_cohort] dirs.")
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    # Resolve display labels up front: extra arms are attached by label, and an
    # arm naming a cohort that is not configured would otherwise vanish silently.
    resolved = [(labels[i] if i < len(labels) and labels[i] else
                 os.path.basename(os.path.normpath(d)))
                for i, d in enumerate(dirs)]
    orphans = sorted({(a.get("cohort") or "") for a in extra_arms} - set(resolved))
    if orphans:
        logger.warning(f"  [external_cohort.arms] entries name cohorts that are not "
                       f"configured and will not run: {orphans}. Expected one of "
                       f"{resolved}.")

    logger.info(f"  {len(dirs)} external cohort(s) configured; "
                f"{len(extra_arms)} extra arm(s); "
                f"de-novo panel discovery {'ON' if discovery else 'OFF'}")

    rows = []
    for cohort_dir, label in zip(dirs, resolved):
        cohort_rows = run_one_cohort(
            cohort_dir, config=config, discovery=discovery, label=label)
        if not cohort_rows:
            logger.error(f"  cohort skipped: {cohort_dir}")
            continue
        rows.extend(cohort_rows)

        # Sensitivity arms attached to this cohort, in configured order.
        for spec in extra_arms:
            if (spec.get("cohort") or "") != label:
                continue
            arm = spec.get("arm") or os.path.basename(
                os.path.normpath(spec.get("dir", "")))
            rows.extend(run_extra_arm(spec.get("dir", ""), label, arm, config))

    if not rows:
        logger.error("  no external cohort produced results — summary not written")
        return pd.DataFrame(columns=SUMMARY_COLUMNS)

    # The internal cohort's own numbers, for the side-by-side comparison. Its
    # arm depends on whether Step 14 has committed a panel that Step 15 applied.
    has_panel = any(os.path.exists(os.path.join(result_dir, f))
                    for f in ("panel_cpg_blood_filtered.csv", "panel_cpg.csv"))
    internal_label = os.path.basename(os.path.normpath(result_dir))
    internal_rows = _summary_rows(
        result_dir, f"{internal_label} (internal)",
        "panel_restricted" if has_panel else "full_matrix")
    if internal_rows:
        rows = internal_rows + rows
    else:
        logger.warning("  no internal classifier_summary.csv — the summary will "
                       "carry external cohorts only. Run --step classify first.")

    df = pd.DataFrame(rows, columns=SUMMARY_COLUMNS)
    out_csv = os.path.join(result_dir, SUMMARY_FILE)
    df.to_csv(out_csv, index=False, encoding="utf-8")
    logger.info(f"  [external_cohort] Saved: {out_csv} ({len(df)} rows)")
    return df
