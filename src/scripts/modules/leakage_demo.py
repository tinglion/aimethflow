#!/usr/bin/env python3
"""
leakage_demo.py — quantify the AUC inflation caused by information leakage
during feature selection.

AIMethFlow's central methodological claim is that selecting features on the
full cohort before cross-validation inflates apparent performance, and that
confining selection to each training fold prevents it. This module turns that
claim from a citation into a measurement by running the *same* classifiers over
the *same* cross-validation partitions under three protocols that differ only
in what the selection / imputation / scaling steps are allowed to see.

Arms
----
``in_fold``
    The AIMethFlow protocol. Mann-Whitney pre-filter, LASSO, RF and their union
    are fitted inside each training fold; imputation medians and the standard
    scaler are fitted on the training fold and applied to the held-out fold.

``naive_selection``
    Feature selection is fitted once on the complete matrix (labels included),
    then frozen; imputation and scaling remain in-fold. This isolates the
    selection component, which is the specific error described in refs 15-18.

``naive_full``
    Selection, imputation and scaling are all fitted on the complete matrix —
    the textbook worst case, and the protocol an unwary reimplementation would
    produce.

``panel_then_in_fold``
    What Step 15 actually does: restrict the matrix to the deployment panel
    (``panel_cpg_blood_filtered.csv``, else ``panel_cpg.csv``) and *then* run
    the in-fold protocol inside it. The panel is defined by Step 14 from the
    labels of the whole cohort, so this arm measures how much of the reported
    performance depends on that upstream, cohort-level selection step. Run only
    when a panel file exists — external cohorts have none.

Because every arm reuses one set of folds, each sample receives exactly one
out-of-fold prediction per arm, so the arms are paired sample-by-sample and the
AUC differences can be tested with DeLong rather than only described.

CA-125
------
When the cohort has CA-125 measurements (``[ml] include_ca125_baseline``), each
arm additionally runs the ``<classifier>+CA125`` fusion models that Step 15
reports, built exactly as ``classifier_utils.run_classifiers`` builds them:
the CA-125 column is median-imputed and standardised on the training fold only
and appended to the already-scaled cfDNA block. Only the cfDNA half of a fusion
model passes through feature selection, so the fusion models inherit whichever
protocol their arm defines and are directly comparable across arms.

CA-125 *alone* is a directly measured clinical assay — no training, no feature
selection — so it is identical under every protocol. It is therefore emitted
once under the pseudo-arm ``ca125_only`` rather than repeated per arm, and it is
excluded from the DeLong table, where it could only ever compare against itself.

Outputs
-------
- ``<result_dir>/leakage_demo.csv``         — per arm x classifier metrics + bootstrap CI
- ``<result_dir>/leakage_demo_delong.csv``  — paired DeLong tests, each naive arm vs in_fold,
  with a paired bootstrap CI on the AUC difference itself (``delta_ci_low/high``): a
  non-significant p says the data cannot resolve the difference, only the interval says how
  large a difference is still compatible with it, and both are needed to call an arm robust
- ``<result_dir>/leakage_demo_ca125_gain.csv`` — within each arm, what each half of a fusion
  model contributes: fusion vs its own cfDNA classifier, and fusion vs CA-125 alone
- ``<result_dir>/leakage_demo.png``         — grouped AUC comparison, cfDNA-only models
- ``<result_dir>/leakage_demo_ca125.png``   — same for the CA-125 fusion models (when present)

Runs on any cohort that has a ``feature_matrix.npz``, including external ones
prepared by ``scripts/prepare_external_matrix.py``.
"""

import os
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, confusion_matrix, accuracy_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from utils.logger import logger
from scripts.modules.classifier_utils import (
    _load_feature_matrix,
    _apply_panel_filter,
    _build_classifiers,
    _select_features_in_fold,
    _impute_fit,
    _impute_transform,
)
from scripts.modules.imeta_utils import _fastdelong_cov

# arm -> what each preprocessing stage is allowed to see
#   panel          : restrict to the cohort-defined deployment panel first
#   select_full    : fit feature selection on the complete matrix
#   transform_full : fit imputation medians and the scaler on the complete matrix
ARMS = {
    "in_fold":            {"panel": False, "select_full": False, "transform_full": False},
    "naive_selection":    {"panel": False, "select_full": True,  "transform_full": False},
    "naive_full":         {"panel": False, "select_full": True,  "transform_full": True},
    "panel_then_in_fold": {"panel": True,  "select_full": False, "transform_full": False},
}
REFERENCE_ARM = "in_fold"

# CA-125 alone belongs to no protocol — see the module docstring.
CA125_ARM = "ca125_only"
CA125_SUFFIX = "+CA125"


def _load_ca125_vector(config, sample_ids):
    """CA-125 values aligned with ``sample_ids`` (NaN where unmeasured).

    Mirrors the loader in ``classifier_utils.run_classifiers`` so the fusion
    arms here reproduce Step 15 rather than approximate it. Returns
    ``(values, threshold)``, or ``(None, threshold)`` when CA-125 is disabled,
    absent from the sample Excel, or unavailable for this cohort — which is the
    normal case for every external cohort.
    """
    ml_cfg = (config or {}).get("ml", {})
    thresh = float(ml_cfg.get("ca125_threshold", 35.0))
    if not bool(ml_cfg.get("include_ca125_baseline", True)):
        logger.info("  [leakage] CA-125 arms disabled by [ml] include_ca125_baseline")
        return None, thresh
    try:
        from methylation.sample_grouping import load_ca125_map
        from scripts.modules.classifier_utils import _build_sid_to_case_map, _norm_id

        excel_sheet = (config or {}).get("plasma", {}).get("excel_sheet", "plasma")
        ca125_map = load_ca125_map(category=excel_sheet)
        if not ca125_map:
            logger.info("  [leakage] CA-125 arms skipped — no CA-125 column in sample Excel")
            return None, thresh
        sid_to_case = _build_sid_to_case_map(config)
        vals = np.array(
            [ca125_map.get(_norm_id(sid_to_case.get(s) or s), np.nan) for s in sample_ids],
            dtype=float,
        )
    except Exception as exc:
        logger.warning(f"  [leakage] CA-125 arms skipped (load failure: {exc})")
        return None, thresh

    n_with = int((~np.isnan(vals)).sum())
    if n_with == 0:
        logger.info("  [leakage] CA-125 arms skipped — no sample carries a CA-125 value")
        return None, thresh
    logger.info(f"  [leakage] CA-125: {n_with}/{len(sample_ids)} samples measured "
                f"(clinical cutoff {thresh:.1f} U/mL)")
    return vals, thresh


def _metrics(y_true, y_prob, y_pred):
    """AUC / accuracy / sensitivity / specificity for one arm x classifier."""
    auc = roc_auc_score(y_true, y_prob)
    acc = accuracy_score(y_true, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    return auc, acc, sens, spec


def _auc_ci(y_true, y_prob, n_boot=1000, seed=42):
    """Percentile bootstrap CI for a single AUC, resampling samples."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        vals.append(roc_auc_score(y_true[idx], y_prob[idx]))
    if not vals:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def _delta_ci(y_true, prob_a, prob_b, n_boot=1000, seed=42):
    """Percentile bootstrap CI for the *paired* AUC difference (a - b).

    Both AUCs are recomputed on the same resampled indices in each iteration,
    so the interval describes the difference itself rather than being inferred
    from two marginal intervals — overlapping marginal CIs say nothing about a
    paired difference, which is exactly the misreading this column prevents.
    """
    rng = np.random.default_rng(seed)
    n = len(y_true)
    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        vals.append(roc_auc_score(y_true[idx], prob_a[idx])
                    - roc_auc_score(y_true[idx], prob_b[idx]))
    if not vals:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def _delong_p(y_true, prob_a, prob_b):
    """Two-sided DeLong test for two paired AUCs. Returns (auc_a, auc_b, z, p)."""
    from scipy import stats as _stats
    aucs, cov = _fastdelong_cov(np.asarray(y_true), np.vstack([prob_a, prob_b]))
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return float(aucs[0]), float(aucs[1]), float("nan"), float("nan")
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    p = 2.0 * (1.0 - _stats.norm.cdf(abs(z)))
    return float(aucs[0]), float(aucs[1]), float(z), float(p)


def _run_arm(arm, X, y, region_in_88, config, splitter, classifiers,
             ca125=None, make_fusion=None):
    """Cross-validate every classifier under one leakage protocol.

    Predictions are written back at the held-out samples' own positions, so the
    returned vectors are aligned across arms and can be compared pairwise.

    When ``ca125`` is supplied, each fold additionally fits the
    ``<classifier>+CA125`` fusion models on the cfDNA block with a
    training-fold-standardised CA-125 column appended. ``make_fusion`` must be a
    zero-argument factory returning fresh, unfitted classifier instances.
    """
    select_on_full = ARMS[arm]["select_full"]
    transform_on_full = ARMS[arm]["transform_full"]
    clf_names = list(classifiers.keys())
    fusion_names = ([f"{c}{CA125_SUFFIX}" for c in clf_names]
                    if ca125 is not None and make_fusion is not None else [])
    n = X.shape[0]

    y_prob = {c: np.full(n, np.nan) for c in clf_names + fusion_names}
    y_pred = {c: np.full(n, np.nan) for c in clf_names + fusion_names}
    fold_n_feat = []

    full_sel = None
    if select_on_full:
        # The leak: selection sees every label, including the held-out ones.
        full_sel = _select_features_in_fold(X, y, region_in_88, config)
        if len(full_sel) == 0:
            full_sel = np.arange(X.shape[1])

    for train_idx, test_idx in splitter.split(X, y):
        X_tr_raw, X_te_raw = X[train_idx], X[test_idx]
        y_tr = y[train_idx]

        if select_on_full:
            sel_idx = full_sel
        else:
            sel_idx = _select_features_in_fold(X_tr_raw, y_tr, region_in_88, config)
            if len(sel_idx) == 0:
                sel_idx = np.arange(X.shape[1])
        fold_n_feat.append(len(sel_idx))

        X_tr = X_tr_raw[:, sel_idx]
        X_te = X_te_raw[:, sel_idx]

        if transform_on_full:
            # Second leak: medians and scaling statistics are pooled over the
            # whole cohort, so held-out samples influence their own transform.
            _, meds = _impute_fit(X[:, sel_idx])
            X_tr_imp = _impute_transform(X_tr, meds)
            X_te_imp = _impute_transform(X_te, meds)
            scaler = StandardScaler().fit(_impute_transform(X[:, sel_idx], meds))
        else:
            X_tr_imp, meds = _impute_fit(X_tr)
            X_te_imp = _impute_transform(X_te, meds)
            scaler = StandardScaler().fit(X_tr_imp)

        X_tr_sc = scaler.transform(X_tr_imp)
        X_te_sc = scaler.transform(X_te_imp)

        for clf_name, clf in classifiers.items():
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    clf.fit(X_tr_sc, y_tr)
                    y_prob[clf_name][test_idx] = clf.predict_proba(X_te_sc)[:, 1]
                    y_pred[clf_name][test_idx] = clf.predict(X_te_sc)
            except Exception as exc:
                logger.debug(f"[leakage] {arm}/{clf_name} fold error: {exc}")

        if not fusion_names:
            continue

        # ── Fusion models: cfDNA block + CA-125 as one extra column ──────────
        # The CA-125 median and scaler are fitted on the training fold only, so
        # the fusion models add no leakage of their own beyond whatever the arm
        # already defines. CA-125 is in U/mL, orders of magnitude away from beta
        # values, hence its own scaler rather than a shared one.
        ca_tr, ca_te = ca125[train_idx], ca125[test_idx]
        med = np.nanmedian(ca_tr)
        if med != med:                       # every training value missing
            med = 0.0
        ca_tr_i = np.where(np.isnan(ca_tr), med, ca_tr).reshape(-1, 1)
        ca_te_i = np.where(np.isnan(ca_te), med, ca_te).reshape(-1, 1)
        ca_scaler = StandardScaler().fit(ca_tr_i)
        X_tr_aug = np.hstack([X_tr_sc, ca_scaler.transform(ca_tr_i)])
        X_te_aug = np.hstack([X_te_sc, ca_scaler.transform(ca_te_i)])

        for clf_name, clf in make_fusion().items():
            key = f"{clf_name}{CA125_SUFFIX}"
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    clf.fit(X_tr_aug, y_tr)
                    y_prob[key][test_idx] = clf.predict_proba(X_te_aug)[:, 1]
                    y_pred[key][test_idx] = clf.predict(X_te_aug)
            except Exception as exc:
                logger.debug(f"[leakage] {arm}/{key} fold error: {exc}")

    n_feat_mean = int(round(np.mean(fold_n_feat))) if fold_n_feat else 0
    return y_prob, y_pred, n_feat_mean


def _save_plot(df, result_dir, cohort_label, out_name="leakage_demo.png",
               title=None):
    """Grouped AUC bars, one group per classifier, one bar per arm.

    ``df`` is expected to be pre-filtered to a single family of models (the
    cfDNA-only classifiers or the CA-125 fusion models); mixing the two would
    make the figure too wide to read.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        logger.warning(f"[leakage] plot skipped: {exc}")
        return

    arms = [a for a in ARMS if a in set(df["arm"])]
    clfs = list(dict.fromkeys(df["classifier"]))
    if not arms or not clfs:
        return
    colors = {"in_fold": "#4dac26", "naive_selection": "#e08214",
              "naive_full": "#d6604d", "panel_then_in_fold": "#5e8fc7"}

    x = np.arange(len(clfs))
    width = 0.8 / max(len(arms), 1)
    fig, ax = plt.subplots(figsize=(1.7 * len(clfs) + 3, 4.5))

    for i, arm in enumerate(arms):
        sub = df[df["arm"] == arm].set_index("classifier")
        vals = np.array([sub.loc[c, "auc"] if c in sub.index else np.nan for c in clfs])
        lo = np.array([sub.loc[c, "auc_ci_low"] if c in sub.index else np.nan for c in clfs])
        hi = np.array([sub.loc[c, "auc_ci_high"] if c in sub.index else np.nan for c in clfs])
        # Bootstrap CIs are asymmetric around the point estimate, so pass both legs.
        err = np.vstack([np.clip(vals - lo, 0, None), np.clip(hi - vals, 0, None)])
        ax.bar(x + i * width - 0.4 + width / 2, vals, width,
               yerr=err, capsize=2.5,
               error_kw={"lw": 0.9, "ecolor": "#333333"},
               label=arm, color=colors.get(arm, "#888888"), edgecolor="white")

    # Annotate the reference arm and the largest inflation observed per classifier.
    ref = df[df["arm"] == REFERENCE_ARM].set_index("classifier")
    for j, c in enumerate(clfs):
        if c not in ref.index:
            continue
        base = ref.loc[c, "auc"]
        ax.text(x[j] - 0.4 + width / 2, base + 0.015, f"{base:.3f}",
                ha="center", fontsize=8)
        others = df[(df["arm"] != REFERENCE_ARM) & (df["classifier"] == c)]
        if not others.empty:
            top = others.loc[others["auc"].idxmax()]
            delta = top["auc"] - base
            # Clear the tallest CI whisker in this group so the label never
            # lands on top of an error bar.
            ceiling = max(others["auc_ci_high"].max(), top["auc"])
            ax.text(x[j] + width, ceiling + 0.025, f"{delta:+.3f}",
                    ha="center", fontsize=8,
                    color="#b2182b" if delta > 0.02 else "#666666")

    ax.set_xticks(x)
    ax.set_xticklabels(clfs)
    ax.set_ylabel("AUC (out-of-fold, 95% bootstrap CI)")
    ax.set_ylim(0.4, 1.09)
    ax.axhline(0.5, color="#999999", lw=0.8, ls="--")
    ax.set_title(title or
                 f"Effect of feature-selection protocol on apparent AUC — {cohort_label}")
    ax.legend(frameon=False, ncol=len(arms), loc="upper center",
              bbox_to_anchor=(0.5, -0.10))
    ax.spines[["top", "right"]].set_visible(False)

    out = os.path.join(result_dir, out_name)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    logger.info(f"  [leakage] Saved: {out}")


def _write_ca125_gain(result_dir, cohort_label, y, probs, ca125, n_boot):
    """Quantify what each half of a fusion model contributes, within each arm.

    Two paired comparisons per (arm, classifier), both against the same
    ``<clf>+CA125`` model and both inside a single arm, so neither is affected
    by the protocol question §4 of the GPB write-up settles:

    ``vs_cfdna``   — fusion vs its own cfDNA-only classifier: what CA-125 adds.
    ``vs_ca125``   — fusion vs CA-125 alone: what the cfDNA panel adds over the
                     existing standard of care, which is the clinically
                     load-bearing claim.

    ``delta_pct`` is the relative increase over the baseline AUC.
    ``headroom_pct`` is the share of the baseline's remaining distance to a
    perfect AUC that the fusion closes — the fairer of the two, because AUC is
    floored at 0.5 and a relative percentage flatters a weak baseline.

    Writes ``leakage_demo_ca125_gain.csv``; a no-op for cohorts without CA-125.
    """
    if ca125 is None or not probs:
        return pd.DataFrame()

    rows = []
    for arm in ARMS:
        if arm not in probs:
            continue
        for key, p_fus in probs[arm].items():
            if not key.endswith(CA125_SUFFIX):
                continue
            base_clf = key[: -len(CA125_SUFFIX)]
            baselines = [("vs_cfdna", base_clf, probs[arm].get(base_clf)),
                         ("vs_ca125", "CA-125", np.asarray(ca125, dtype=float))]
            for comparison, base_name, p_base in baselines:
                if p_base is None:
                    continue
                ok = ~np.isnan(p_fus) & ~np.isnan(p_base)
                if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
                    continue
                auc_f, auc_b, z, p = _delong_p(y[ok], p_fus[ok], p_base[ok])
                lo, hi = _delta_ci(y[ok], p_fus[ok], p_base[ok], n_boot=n_boot)
                delta = auc_f - auc_b
                head = 1.0 - auc_b
                rows.append({
                    "cohort": cohort_label, "arm": arm, "comparison": comparison,
                    "model": key, "baseline": base_name, "n_samples": int(ok.sum()),
                    "auc_model": round(auc_f, 4), "auc_baseline": round(auc_b, 4),
                    "delta_auc": round(delta, 4),
                    "delta_ci_low": round(lo, 4), "delta_ci_high": round(hi, 4),
                    "delta_pct": round(100.0 * delta / auc_b, 2) if auc_b > 0 else np.nan,
                    "headroom_pct": round(100.0 * delta / head, 2) if head > 1e-9 else np.nan,
                    "delong_z": round(z, 4) if z == z else np.nan,
                    "delong_p": p,
                })

    gdf = pd.DataFrame(rows)
    if gdf.empty:
        return gdf
    out = os.path.join(result_dir, "leakage_demo_ca125_gain.csv")
    gdf.to_csv(out, index=False, encoding="utf-8")
    logger.info(f"  [leakage] Saved: {out}")

    ref = gdf[gdf["arm"] == REFERENCE_ARM]
    for comparison, label in (("vs_cfdna", "CA-125 added to cfDNA"),
                              ("vs_ca125", "cfDNA added to CA-125")):
        sub = ref[ref["comparison"] == comparison]
        if sub.empty:
            continue
        logger.info(f"  [leakage] {label} ({REFERENCE_ARM} arm):")
        for _, r in sub.iterrows():
            sig = "*" if r["delong_p"] < 0.05 else " "
            logger.info(f"      {r['model']:<12} {r['auc_baseline']:.3f} -> "
                        f"{r['auc_model']:.3f}  {r['delta_auc']:+.3f} "
                        f"({r['delta_pct']:+.1f}%, {r['headroom_pct']:+.0f}% of headroom) "
                        f"CI [{r['delta_ci_low']:+.3f}, {r['delta_ci_high']:+.3f}] "
                        f"p={r['delong_p']:.2e}{sig}")
    return gdf


def run_leakage_demo(workdir, result_dir, config, benign_samples=None,
                     malignant_samples=None, cohort_label=None, n_boot=1000):
    """Run the three leakage arms and write the comparison tables + figure.

    Returns (metrics_df, delong_df).
    """
    logger.info("━━ Leakage demonstration (naive vs in-fold feature selection) ━━")

    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples or [], malignant_samples or [], config
    )
    if X.shape[1] == 0:
        logger.error("[leakage] Feature matrix empty — aborting")
        return pd.DataFrame(), pd.DataFrame()

    ml_cfg = (config or {}).get("ml", {})
    exclude_set = set(ml_cfg.get("exclude_samples", []))
    if exclude_set:
        keep = np.array([sid not in exclude_set for sid in sample_ids])
        X, y = X[keep], y[keep]
        sample_ids = [s for s, k in zip(sample_ids, keep) if k]

    n_neg, n_pos = int((y == 0).sum()), int((y == 1).sum())
    spw = n_neg / n_pos if n_pos > 0 else 1.0
    ca125, ca125_thresh = _load_ca125_vector(config, sample_ids)
    cohort_label = cohort_label or os.path.basename(os.path.normpath(result_dir))
    logger.info(f"  cohort={cohort_label}  n={X.shape[0]} "
                f"({n_neg} benign / {n_pos} malignant)  features={X.shape[1]}")

    cv_folds = ml_cfg.get("cv_folds", 5)
    splitter = (LeaveOneOut() if cv_folds <= 0
                else StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=42))
    logger.info(f"  CV: {'LOO' if cv_folds <= 0 else f'stratified {cv_folds}-fold'} "
                f"(identical partitions across all arms)")

    # The panel arm reproduces Step 15, which restricts the matrix to the
    # cohort-defined panel before cross-validating. Absent a panel file (any
    # external cohort) the filter is a no-op, so the arm is dropped rather than
    # silently duplicating in_fold.
    X_panel, _, r88_panel, panel_mode = _apply_panel_filter(
        X, feature_names, region_in_88, result_dir)
    has_panel = not panel_mode.startswith("no panel restriction")
    if has_panel:
        logger.info(f"  panel arm: {panel_mode} — {X_panel.shape[1]} features")
    else:
        logger.info("  panel arm skipped — no panel_cpg file in this result dir")

    rows, probs = [], {}

    # CA-125 alone needs no folds: it is a lookup, so its out-of-fold prediction
    # is the measurement itself and is identical under every protocol. Emit it
    # once, outside the arm loop, under its own pseudo-arm.
    if ca125 is not None:
        ok = ~np.isnan(ca125)
        if ok.sum() > 0 and len(np.unique(y[ok])) >= 2:
            hard = (ca125[ok] >= ca125_thresh).astype(float)
            auc, acc, sens, spec = _metrics(y[ok], ca125[ok], hard)
            lo, hi = _auc_ci(y[ok], ca125[ok], n_boot=n_boot)
            rows.append({
                "cohort": cohort_label, "arm": CA125_ARM, "classifier": "CA-125",
                "n_samples": int(ok.sum()), "n_features_mean": 1,
                "auc": round(auc, 4), "auc_ci_low": round(lo, 4),
                "auc_ci_high": round(hi, 4), "accuracy": round(acc, 4),
                "sensitivity": round(sens, 4), "specificity": round(spec, 4),
            })
            logger.info(f"  [leakage] arm={CA125_ARM}  CA-125={auc:.3f} "
                        f"(protocol-invariant — measured, not learned)")

    for arm in ARMS:
        if ARMS[arm]["panel"] and not has_panel:
            continue
        logger.info(f"  [leakage] arm={arm} ...")
        # Rebuild classifiers per arm so no fitted state carries over.
        classifiers = _build_classifiers(ml_cfg, spw)
        X_arm, r88_arm = (X_panel, r88_panel) if ARMS[arm]["panel"] else (X, region_in_88)
        y_prob, y_pred, n_feat = _run_arm(
            arm, X_arm, y, r88_arm, config, splitter, classifiers,
            ca125=ca125, make_fusion=lambda: _build_classifiers(ml_cfg, spw))
        probs[arm] = y_prob

        for clf_name in y_prob:
            p, q = y_prob[clf_name], y_pred[clf_name]
            ok = ~np.isnan(p)
            if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
                continue
            auc, acc, sens, spec = _metrics(y[ok], p[ok], q[ok])
            lo, hi = _auc_ci(y[ok], p[ok], n_boot=n_boot)
            rows.append({
                "cohort": cohort_label, "arm": arm, "classifier": clf_name,
                "n_samples": int(ok.sum()), "n_features_mean": n_feat,
                "auc": round(auc, 4), "auc_ci_low": round(lo, 4),
                "auc_ci_high": round(hi, 4), "accuracy": round(acc, 4),
                "sensitivity": round(sens, 4), "specificity": round(spec, 4),
            })
        aucs = {r["classifier"]: r["auc"] for r in rows if r["arm"] == arm}
        logger.info("    " + "  ".join(f"{c}={v:.3f}" for c, v in aucs.items()))

    df = pd.DataFrame(rows)
    out_csv = os.path.join(result_dir, "leakage_demo.csv")
    df.to_csv(out_csv, index=False, encoding="utf-8")
    logger.info(f"  [leakage] Saved: {out_csv}")

    # Paired DeLong: each naive arm against the in-fold reference.
    d_rows = []
    for arm in ARMS:
        if arm == REFERENCE_ARM or arm not in probs:
            continue
        for clf_name in probs[REFERENCE_ARM]:
            pa, pb = probs[arm].get(clf_name), probs[REFERENCE_ARM].get(clf_name)
            if pa is None or pb is None:
                continue
            ok = ~np.isnan(pa) & ~np.isnan(pb)
            if ok.sum() == 0 or len(np.unique(y[ok])) < 2:
                continue
            auc_naive, auc_ref, z, p = _delong_p(y[ok], pa[ok], pb[ok])
            # A non-significant p says the data cannot resolve the difference;
            # only the interval says how large a difference is still compatible
            # with it. Both are needed to claim an AUC is protocol-robust.
            d_lo, d_hi = _delta_ci(y[ok], pa[ok], pb[ok], n_boot=n_boot)
            d_rows.append({
                "cohort": cohort_label, "classifier": clf_name,
                "arm": arm, "auc_arm": round(auc_naive, 4),
                "arm_reference": REFERENCE_ARM, "auc_reference": round(auc_ref, 4),
                "delta_auc": round(auc_naive - auc_ref, 4),
                "delta_ci_low": round(d_lo, 4), "delta_ci_high": round(d_hi, 4),
                "delong_z": round(z, 4) if z == z else np.nan,
                "delong_p": p,
            })
    ddf = pd.DataFrame(d_rows)
    d_csv = os.path.join(result_dir, "leakage_demo_delong.csv")
    ddf.to_csv(d_csv, index=False, encoding="utf-8")
    logger.info(f"  [leakage] Saved: {d_csv}")

    _write_ca125_gain(result_dir, cohort_label, y, probs, ca125, n_boot)

    if not df.empty:
        is_fusion = df["classifier"].str.endswith(CA125_SUFFIX)
        _save_plot(df[~is_fusion & (df["arm"] != CA125_ARM)], result_dir, cohort_label)
        if is_fusion.any():
            _save_plot(
                df[is_fusion], result_dir, cohort_label,
                out_name="leakage_demo_ca125.png",
                title=(f"Same comparison for the CA-125 fusion models — {cohort_label}"),
            )

    if not ddf.empty:
        logger.info("\n  AUC difference vs the leakage-free in-fold reference")
        logger.info("  ┌────────────────────┬──────────┬───────────┬───────────┬───────────┬──────────┐")
        logger.info("  │ Arm                │ Clf      │ arm AUC   │ in-fold   │ ΔAUC      │ DeLong p │")
        logger.info("  ├────────────────────┼──────────┼───────────┼───────────┼───────────┼──────────┤")
        for _, r in ddf.iterrows():
            logger.info(f"  │ {r['arm']:<18} │ {r['classifier']:<8} │"
                        f" {r['auc_arm']:>9.3f} │ {r['auc_reference']:>9.3f} │"
                        f" {r['delta_auc']:>+9.3f} │ {r['delong_p']:>8.2e} │")
        logger.info("  └────────────────────┴──────────┴───────────┴───────────┴───────────┴──────────┘")

    return df, ddf
