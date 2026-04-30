#!/usr/bin/env python3
"""
feature_utils.py - Feature selection module (Step 14)

Builds a feature matrix from per-sample _CpG.bedGraph files and runs
LASSO (L1 logistic) + Random Forest feature selection.

Outputs:
    feature_importance.csv  — LASSO coefficient + RF importance per site
    panel_cpg.csv           — Selected panel CpG sites (after coverage QC)
    panel_cpg_depth.csv     — Per-site depth statistics
    feature_matrix.npz      — Full feature matrix (pre-panel, for classify step)
    wilcoxon_stats.csv      — Per-site Wilcoxon test statistics

Dependencies: scikit-learn, numpy, pandas, scipy
"""

import os
import bisect
import numpy as np
import pandas as pd
from collections import defaultdict

from utils.logger import logger


# ──────────────────────────────────────────────────────
# Region membership check
# ──────────────────────────────────────────────────────

def _check_region_membership(feature_names, candidate_regions_file):
    """
    Check whether each CpG site (format "chr:pos") falls within the 88
    tissue candidate regions.  Uses bisect binary search:
    O(n_features × log(n_regions_per_chrom)).

    Returns:
        np.ndarray(bool), shape (n_features,)
    """
    in_region = np.zeros(len(feature_names), dtype=bool)

    if not candidate_regions_file or not os.path.exists(candidate_regions_file):
        return in_region

    # Build sorted interval lists per chromosome
    regions = defaultdict(list)
    with open(candidate_regions_file, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    regions[parts[0]].append((int(parts[1]), int(parts[2])))
                except ValueError:
                    continue
    for chrom in regions:
        regions[chrom].sort()

    # Extract start arrays for bisect
    starts_map = {chrom: [r[0] for r in rs] for chrom, rs in regions.items()}

    for i, name in enumerate(feature_names):
        if ":" not in name:
            continue
        chrom, pos_str = name.split(":", 1)
        pos = int(pos_str)
        chrom_regions = regions.get(chrom)
        if not chrom_regions:
            continue
        starts = starts_map[chrom]
        idx = bisect.bisect_right(starts, pos) - 1
        if idx >= 0 and pos <= chrom_regions[idx][1]:
            in_region[i] = True

    return in_region


# ──────────────────────────────────────────────────────
# Feature matrix construction
# ──────────────────────────────────────────────────────

def build_feature_matrix(workdir, result_dir, benign_samples, malignant_samples, config=None):
    """
    Build feature matrix and depth matrix from per-sample _CpG.bedGraph files.

    Candidate sites are data-driven: CpGs covered in >= min_feature_coverage
    fraction of all samples.  Each site is also annotated with whether it falls
    inside the 88 tissue candidate regions (used for stratified Wilcoxon filtering).

    Args:
        workdir:           Temp dir containing {sample}_CpG.bedGraph files
        result_dir:        Result dir (diff_methylation.txt read for info only)
        benign_samples:    List of benign sample IDs
        malignant_samples: List of malignant sample IDs
        config:            Config dict

    Returns:
        X:             np.ndarray (n_samples, n_features), methylation 0-1, NaN = no coverage
        y:             np.ndarray (n_samples,), 0=benign 1=malignant
        sample_ids:    list[str]
        feature_names: list[str], format "chr:pos"
        D:             np.ndarray (n_samples, n_features), depth, 0 = no coverage
        region_in_88:  np.ndarray(bool) (n_features,), True if in 88 candidate regions
    """
    ml_cfg = (config or {}).get("ml", {})
    min_coverage = ml_cfg.get("min_feature_coverage", 0.5)

    all_samples = [(sid, 0) for sid in benign_samples] + \
                  [(sid, 1) for sid in malignant_samples]
    n_samples = len(all_samples)

    # ── Pass 1: read all bedGraphs, aggregate site coverage counts ──
    pos_count = {}
    raw_data = {}

    missing_samples = []
    for sample_id, _ in all_samples:
        bedgraph = os.path.join(workdir, f"{sample_id}_CpG.bedGraph")
        if not os.path.exists(bedgraph):
            logger.warning(f"  bedGraph not found, skipping: {bedgraph}")
            missing_samples.append(sample_id)
            raw_data[sample_id] = {}
            continue

        bg = pd.read_csv(bedgraph, sep="\t", header=None,
                         names=["chrom", "start", "end", "meth_pct",
                                "meth_count", "unmeth_count"])
        bg = bg[pd.to_numeric(bg["start"], errors="coerce").notna()].copy()
        bg["start"] = bg["start"].astype(int)
        bg["depth"] = (pd.to_numeric(bg["meth_count"], errors="coerce").fillna(0) +
                       pd.to_numeric(bg["unmeth_count"], errors="coerce").fillna(0)).astype(int)
        bg = bg[bg["depth"] > 0]

        sd = {}
        for row in bg.itertuples(index=False):
            key = (row.chrom, row.start)
            sd[key] = (float(row.meth_pct), int(row.depth))
            pos_count[key] = pos_count.get(key, 0) + 1
        raw_data[sample_id] = sd

    if missing_samples:
        logger.warning(f"  {len(missing_samples)} sample(s) missing bedGraph: {missing_samples}")

    # ── Filter candidate sites by coverage rate ──
    min_count = max(1, int(round(n_samples * min_coverage)))
    candidate_positions = sorted(
        [(chrom, pos) for (chrom, pos), cnt in pos_count.items() if cnt >= min_count],
        key=lambda x: (x[0], x[1])
    )
    n_features = len(candidate_positions)

    # Compare against diff_methylation.txt (informational only)
    diff_file = os.path.join(result_dir, "diff_methylation.txt")
    if os.path.exists(diff_file):
        try:
            _ddf = pd.read_csv(diff_file, sep="\t", header=None,
                               names=["chrom", "pos", "m1", "m2", "diff"])
            _ddf = _ddf[pd.to_numeric(_ddf["pos"], errors="coerce").notna()]
            logger.info(f"  Data-driven candidates: {len(pos_count)} total covered sites, "
                        f">={min_coverage:.0%} sample coverage → {n_features} sites  "
                        f"(diff_methylation.txt DMC: {len(_ddf)}, for reference only)")
        except Exception:
            logger.info(f"  Data-driven candidates: {len(pos_count)} total covered sites, "
                        f">={min_coverage:.0%} sample coverage → {n_features} sites")
    else:
        logger.info(f"  Data-driven candidates: {len(pos_count)} total covered sites, "
                    f">={min_coverage:.0%} sample coverage → {n_features} sites")

    if n_features == 0:
        y = np.array([label for _, label in all_samples])
        sample_ids = [sid for sid, _ in all_samples]
        return (np.zeros((n_samples, 0)), y, sample_ids, [], np.zeros((n_samples, 0), dtype=np.int32),
                np.zeros(0, dtype=bool))

    feature_names = [f"{chrom}:{pos}" for chrom, pos in candidate_positions]
    pos_index = {cp: i for i, cp in enumerate(candidate_positions)}

    # ── Pass 2: fill feature matrix ──
    X = np.full((n_samples, n_features), np.nan)
    D = np.zeros((n_samples, n_features), dtype=np.int32)
    y = np.array([label for _, label in all_samples])
    sample_ids = [sid for sid, _ in all_samples]

    for i, (sample_id, _) in enumerate(all_samples):
        sd = raw_data.get(sample_id, {})
        for (chrom, pos), j in pos_index.items():
            if (chrom, pos) in sd:
                meth_pct, depth = sd[(chrom, pos)]
                X[i, j] = meth_pct / 100.0
                D[i, j] = depth

    # ── Annotate 88-region membership ──
    candidate_regions_file = (config or {}).get("data", {}).get("candidate_regions", "")
    region_in_88 = _check_region_membership(feature_names, candidate_regions_file)
    n_in_88  = int(region_in_88.sum())
    n_out_88 = n_features - n_in_88
    logger.info(f"  In 88 tissue regions: {n_in_88}  Outside: {n_out_88}  "
                f"(total {n_features} candidate sites)")

    logger.info(f"  Feature matrix: {X.shape[0]} samples × {X.shape[1]} features  "
                f"(benign={int(sum(y==0))}, malignant={int(sum(y==1))})")
    logger.info(f"  Missing rate: {np.isnan(X).mean():.2%}")

    return X, y, sample_ids, feature_names, D, region_in_88


def _impute(X):
    """Impute NaN with column median."""
    X_imp = X.copy()
    col_medians = np.nanmedian(X_imp, axis=0)
    for j in range(X_imp.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j] if not np.isnan(col_medians[j]) else 0.0
    return X_imp


# ──────────────────────────────────────────────────────
# LASSO feature selection (L1 regularised logistic)
# ──────────────────────────────────────────────────────

def run_lasso_selection(X, y, config=None):
    """
    Select features using L1-penalised logistic regression (LogisticRegressionCV).

    Imputes NaN → standardises → cross-validates to find optimal C.

    Returns:
        selected_indices: np.ndarray, column indices with non-zero coefficients
        coefs:            np.ndarray, full coefficient vector (0 = not selected)
    """
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.preprocessing import StandardScaler

    ml_cfg = (config or {}).get("ml", {})
    cv_folds = ml_cfg.get("lasso_cv_folds", 5)

    X_imp = _impute(X)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_imp)

    clf = LogisticRegressionCV(
        Cs=20,
        cv=cv_folds,
        penalty="l1",
        solver="liblinear",
        class_weight="balanced",
        max_iter=2000,
        random_state=42,
    )
    clf.fit(X_scaled, y)

    coefs = clf.coef_[0]
    selected_indices = np.where(coefs != 0)[0]
    logger.info(f"  LASSO: {len(selected_indices)}/{len(coefs)} features selected  "
                f"(best C={clf.C_[0]:.4f})")
    return selected_indices, coefs


# ──────────────────────────────────────────────────────
# Random Forest feature importance
# ──────────────────────────────────────────────────────

def run_rf_importance(X, y, config=None):
    """
    Compute feature importances with Random Forest.

    Returns:
        selected_indices: np.ndarray, indices with importance > min_importance
        importances:      np.ndarray, full importance scores
    """
    from sklearn.ensemble import RandomForestClassifier

    ml_cfg = (config or {}).get("ml", {})
    n_estimators = ml_cfg.get("rf_n_estimators", 500)
    min_importance = ml_cfg.get("rf_min_importance", 0.0)

    X_imp = _impute(X)
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )
    rf.fit(X_imp, y)

    importances = rf.feature_importances_
    selected_indices = np.where(importances > min_importance)[0]
    logger.info(f"  RF: {len(selected_indices)}/{len(importances)} features selected  "
                f"(min_importance={min_importance})")
    return selected_indices, importances


# ──────────────────────────────────────────────────────
# Sample outlier detection
# ──────────────────────────────────────────────────────

def detect_sample_outliers(X, y, sample_ids, config=None):
    """
    Identify potentially abnormal samples from multiple angles.
    Does NOT automatically exclude them — results are written to
    outlier_report.csv for the user to review.

    Three detection methods (sample flagged as suspect if >= 2 trigger):
      1. Coverage: NaN fraction > (1 - outlier_coverage_min)
      2. Intra-group PCA distance: z-score > outlier_pca_z
      3. Isolation Forest: bottom 25% anomaly score within group

    Args:
        X:          (n_samples, n_features), NaN = no coverage
        y:          (n_samples,), 0=benign 1=malignant
        sample_ids: list[str]
        config:     config dict

    Returns:
        df:          pd.DataFrame with per-sample metrics, sorted (suspect first)
        suspect_ids: list[str], suspected outlier sample IDs
    """
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        from sklearn.ensemble import IsolationForest
        import pandas as pd
    except ImportError as e:
        logger.warning(f"detect_sample_outliers: missing dependency {e}, skipping")
        return pd.DataFrame(), []

    ml_cfg = (config or {}).get("ml", {})
    coverage_min = ml_cfg.get("outlier_coverage_min", 0.30)
    pca_z_thr    = ml_cfg.get("outlier_pca_z",        2.5)

    n_samples, n_features = X.shape
    rows = []

    nan_frac = np.isnan(X).mean(axis=1)

    # Intra-group PCA + Isolation Forest (computed separately per group)
    pca_z_scores = np.full(n_samples, np.nan)
    iforest_flag = np.zeros(n_samples, dtype=bool)

    for label in (0, 1):
        grp_idx = np.where(y == label)[0]
        if len(grp_idx) < 4:
            continue

        X_grp = X[grp_idx].copy()
        col_meds = np.nanmedian(X_grp, axis=0)
        col_meds = np.where(np.isnan(col_meds), 0.0, col_meds)
        for j in range(X_grp.shape[1]):
            mask = np.isnan(X_grp[:, j])
            X_grp[mask, j] = col_meds[j]

        scaler = StandardScaler()
        X_sc = scaler.fit_transform(X_grp)

        # Method 2: PCA intra-group distance z-score
        n_components = min(2, len(grp_idx) - 1, n_features)
        if n_components >= 1:
            pca = PCA(n_components=n_components, random_state=42)
            coords = pca.fit_transform(X_sc)
            dists  = np.linalg.norm(coords, axis=1)
            mu, sd = dists.mean(), dists.std(ddof=0)
            z = (dists - mu) / (sd + 1e-10)
            pca_z_scores[grp_idx] = z

        # Method 3: Isolation Forest
        if len(grp_idx) >= 8:
            iforest = IsolationForest(
                n_estimators=200, contamination=0.25,
                random_state=42
            )
            pred_if = iforest.fit_predict(X_sc)
            iforest_flag[grp_idx] = (pred_if == -1)

    # Combine flags
    flag_cov = nan_frac > (1.0 - coverage_min)
    flag_pca = pca_z_scores > pca_z_thr
    flag_if  = iforest_flag

    n_flags = flag_cov.astype(int) + flag_pca.astype(int) + flag_if.astype(int)
    suspect  = n_flags >= 2

    label_names = {0: "benign", 1: "malignant"}
    for i, sid in enumerate(sample_ids):
        rows.append({
            "sample_id":       sid,
            "group":           label_names.get(int(y[i]), "unknown"),
            "nan_frac":        round(float(nan_frac[i]), 4),
            "flag_coverage":   bool(flag_cov[i]),
            "pca_z_score":     round(float(pca_z_scores[i]), 4) if not np.isnan(pca_z_scores[i]) else None,
            "flag_pca":        bool(flag_pca[i]),
            "flag_isolation":  bool(flag_if[i]),
            "n_flags":         int(n_flags[i]),
            "suspect":         bool(suspect[i]),
        })

    df = pd.DataFrame(rows).sort_values(["suspect", "n_flags"], ascending=[False, False])
    suspect_ids = [r["sample_id"] for r in rows if r["suspect"]]

    n_suspect = len(suspect_ids)
    logger.info(f"  Sample outlier detection: {n_samples} samples, "
                f"{n_suspect} suspect (>= 2 flags triggered)")
    if suspect_ids:
        logger.info(f"  Suspect samples: {suspect_ids}")
        logger.info("  Add confirmed IDs to config.toml [ml] exclude_samples = [...]")
    else:
        logger.info("  No suspect samples detected")

    return df, suspect_ids


# ──────────────────────────────────────────────────────
# Main entry point: Step 14
# ──────────────────────────────────────────────────────

def select_features(workdir, result_dir, config, benign_samples, malignant_samples,
                    extra_exclude_ids=None):
    """
    Step 14 main function: build feature matrix from bedGraph files,
    run LASSO + RF feature selection, write result files.

    Args:
        extra_exclude_ids: sample IDs to exclude from feature selection (but NOT from
                           the saved feature_matrix.npz, which classifiers need in full).
                           Typically the suspect_ids list from Step 13 outlier_detect.

    Outputs (saved to result_dir):
        feature_importance.csv — LASSO coef + RF importance per site
        panel_cpg.csv          — Selected panel sites
        feature_matrix.npz     — Feature matrix (for classify step; all samples)

    Args:
        workdir:           Temp directory
        result_dir:        Result directory
        config:            Config dict
        benign_samples:    List of benign sample IDs
        malignant_samples: List of malignant sample IDs

    Returns:
        panel_features: list[str], selected feature names
        importance_df:  pd.DataFrame, feature importance table
    """
    logger.info(f"Starting feature selection  "
                f"(benign={len(benign_samples)}, malignant={len(malignant_samples)})")

    # 1. Build feature matrix — load from cache if outlier_detect (Step 13) already ran,
    #    otherwise build from bedGraph files.
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if os.path.exists(fm_path):
        logger.info(f"  Loading cached feature matrix: {fm_path}")
        _fm = np.load(fm_path, allow_pickle=True)
        X            = _fm["X"]
        y            = _fm["y"]
        sample_ids   = list(_fm["sample_ids"])
        feature_names = list(_fm["feature_names"])
        _D_raw       = _fm["D"] if "D" in _fm else None
        if _D_raw is not None and _D_raw.shape[0] == X.shape[0]:
            D        = _D_raw
        else:
            if _D_raw is not None:
                logger.warning(
                    f"  Cached D shape {_D_raw.shape} does not match X shape {X.shape} "
                    f"— depth matrix will be recomputed from bedGraph files."
                )
            D = np.zeros((X.shape[0], X.shape[1]))
        region_in_88 = (_fm["region_in_88"] if "region_in_88" in _fm
                        else np.zeros(X.shape[1], dtype=bool))
        logger.info(f"  Cached: {X.shape[0]} samples × {X.shape[1]} features")
    else:
        X, y, sample_ids, feature_names, D, region_in_88 = build_feature_matrix(
            workdir, result_dir, benign_samples, malignant_samples, config
        )

    if X.shape[1] == 0:
        logger.error("Feature matrix is empty (no valid CpG sites) — skipping feature selection")
        return [], pd.DataFrame()

    # 1.5a. Keep references to the FULL matrix so feature_matrix.npz always covers all samples.
    #        Classifiers (Step 15+) load the full matrix and apply their own exclusions.
    X_save          = X
    y_save          = y
    sample_ids_save = list(sample_ids)
    D_save          = D          # full matrix, before exclusion filtering

    # 1.5b. Exclude outlier samples from the selection computation (Wilcoxon / LASSO / RF).
    #        These samples will still appear in the saved feature_matrix.npz.
    if extra_exclude_ids:
        excl_set = set(extra_exclude_ids)
        keep_mask = np.array([sid not in excl_set for sid in sample_ids])
        n_excl = int((~keep_mask).sum())
        if n_excl:
            excl_list = [sid for sid, k in zip(sample_ids, keep_mask) if not k]
            logger.info(f"  Excluding {n_excl} outlier sample(s) from selection: {excl_list}")
            X          = X[keep_mask]
            y          = y[keep_mask]
            D          = D[keep_mask]
            sample_ids = [sid for sid, k in zip(sample_ids, keep_mask) if k]

    # 1.5. Sample-level Wilcoxon pre-filter
    #   diff_methylation.txt compares merged bedGraphs, so between-group mean differences
    #   are inflated (within-group noise is averaged away).
    #   This step runs Mann-Whitney U per site at the individual-sample level,
    #   keeping only sites that pass both p < wilcox_p AND |Δβ| >= wilcox_delta.
    ml_cfg_pre = (config or {}).get("ml", {})
    wilcox_p     = ml_cfg_pre.get("sample_wilcox_p",     0.05)
    wilcox_delta = ml_cfg_pre.get("sample_delta_beta",   0.10)

    # Preserve pre-Wilcoxon arrays for depth stats (selection subset)
    feature_names_all = list(feature_names)
    region_in_88_all = region_in_88.copy()

    wilcox_p       = ml_cfg_pre.get("sample_wilcox_p",      0.05)   # standard threshold
    wilcox_delta   = ml_cfg_pre.get("sample_delta_beta",    0.10)
    wilcox_p_88    = ml_cfg_pre.get("sample_wilcox_p_88",   0.10)   # relaxed for 88-region sites
    wilcox_delta_88 = ml_cfg_pre.get("sample_delta_beta_88", 0.05)

    try:
        from scipy import stats as _stats
        _HAS_SCIPY = True
    except ImportError:
        _HAS_SCIPY = False
        logger.warning("scipy not installed — skipping sample-level Wilcoxon pre-filter")

    if _HAS_SCIPY:
        X_imp_pre = _impute(X)
        pvals  = np.ones(X.shape[1])
        deltas = np.zeros(X.shape[1])
        b_means = np.zeros(X.shape[1])
        m_means = np.zeros(X.shape[1])

        for j in range(X.shape[1]):
            b_vals = X_imp_pre[y == 0, j]
            m_vals = X_imp_pre[y == 1, j]
            if len(b_vals) >= 3 and len(m_vals) >= 3:
                _, p = _stats.mannwhitneyu(b_vals, m_vals, alternative="two-sided")
                d    = abs(float(np.mean(m_vals)) - float(np.mean(b_vals)))
            else:
                p, d = 1.0, 0.0
            pvals[j]   = p
            deltas[j]  = d
            b_means[j] = float(np.mean(X_imp_pre[y == 0, j]))
            m_means[j] = float(np.mean(X_imp_pre[y == 1, j]))

        # ── Stratified Wilcoxon thresholds ──
        # 88-region sites: relaxed (prior knowledge)
        # Other sites:     standard data-driven threshold
        p_thresh = np.where(region_in_88, wilcox_p_88,    wilcox_p)
        d_thresh = np.where(region_in_88, wilcox_delta_88, wilcox_delta)
        wilcox_mask = (pvals < p_thresh) & (deltas >= d_thresh)

        n_pass    = int(wilcox_mask.sum())
        n_in88    = int(region_in_88.sum())
        n_pass_88 = int((wilcox_mask & region_in_88).sum())
        n_pass_ot = int((wilcox_mask & ~region_in_88).sum())

        logger.info(f"Stratified Wilcoxon results: {n_pass}/{X.shape[1]} sites pass")
        logger.info(f"  In 88 regions  ({n_in88} sites, p<{wilcox_p_88}, |Δβ|≥{wilcox_delta_88}): "
                    f"{n_pass_88} pass")
        logger.info(f"  Outside 88 regions ({X.shape[1]-n_in88} sites, p<{wilcox_p}, |Δβ|≥{wilcox_delta}): "
                    f"{n_pass_ot} pass")
        logger.info(f"  Sample-level |Δβ| — mean={deltas.mean():.3f}  "
                    f"median={float(np.median(deltas)):.3f}  max={deltas.max():.3f}")
        logger.info(f"  p-values — median={float(np.median(pvals)):.3f}  "
                    f"<0.05: {int((pvals<0.05).sum())}  <0.10: {int((pvals<0.10).sum())}  "
                    f"<0.20: {int((pvals<0.20).sum())}")

        # Sites close to threshold (borderline)
        near_mask = (~wilcox_mask) & (
            (deltas >= d_thresh * 0.8) | (pvals < p_thresh * 2.0)
        )
        if near_mask.any():
            bl_names  = [feature_names_all[j] for j in np.where(near_mask)[0]]
            bl_deltas = deltas[near_mask]
            bl_pvals  = pvals[near_mask]
            top_bl = sorted(zip(bl_names, bl_deltas, bl_pvals), key=lambda x: -x[1])[:5]
            logger.info("  Borderline sites (near threshold, did not pass):")
            for nm, dd, pp in top_bl:
                logger.info(f"    {nm}  |Δβ|={dd:.3f}  p={pp:.4f}")

        # Save full Wilcoxon stats table
        wilcox_df = pd.DataFrame({
            "feature":         feature_names_all,
            "in_88_regions":   region_in_88_all,
            "benign_mean":     b_means.round(4),
            "malignant_mean":  m_means.round(4),
            "sample_delta":    deltas.round(4),
            "wilcox_pval":     pvals.round(6),
            "p_threshold":     p_thresh.round(4),
            "delta_threshold": d_thresh.round(4),
            "pass_wilcox":     wilcox_mask,
        }).sort_values("sample_delta", ascending=False)
        wpath = os.path.join(result_dir, "wilcoxon_stats.csv")
        wilcox_df.to_csv(wpath, index=False)
        logger.info(f"  Wilcoxon stats saved: {wpath}")

        if n_pass == 0:
            logger.warning("No sites pass stratified Wilcoxon; retrying with relaxed thresholds "
                           "p<0.20 / |Δβ|>=0.05")
            wilcox_mask = (pvals < 0.20) & (deltas >= 0.05)
            n_pass = int(wilcox_mask.sum())
            logger.info(f"  After relaxation: {n_pass} sites pass")

        if n_pass > 0:
            X             = X[:, wilcox_mask]
            D             = D[:, wilcox_mask]
            region_in_88  = region_in_88_all[wilcox_mask]
            feature_names = [fn for fn, v in zip(feature_names_all, wilcox_mask) if v]
            pvals_kept    = pvals[wilcox_mask]
            deltas_kept   = deltas[wilcox_mask]
        else:
            pvals_kept  = pvals
            deltas_kept = deltas
    else:
        pvals_kept  = np.ones(X.shape[1])
        deltas_kept = np.zeros(X.shape[1])

    if X.shape[1] == 0:
        logger.error("Feature matrix empty after Wilcoxon filter — skipping feature selection")
        return [], pd.DataFrame()

    # 2. LASSO feature selection
    logger.info("Running LASSO (L1 logistic) feature selection...")
    lasso_indices, lasso_coefs = run_lasso_selection(X, y, config)

    # 3. RF feature importance
    logger.info("Running Random Forest feature importance...")
    rf_indices, rf_importances = run_rf_importance(X, y, config)

    # 4. Combine selection results
    ml_cfg = (config or {}).get("ml", {})
    selection_mode = ml_cfg.get("feature_selection_mode", "union")

    if selection_mode == "lasso":
        final_indices = lasso_indices
    elif selection_mode == "rf":
        final_indices = rf_indices
    else:  # "union" (default)
        final_indices = np.unique(np.concatenate([lasso_indices, rf_indices]))

    panel_features = [feature_names[i]      for i in final_indices]
    panel_in_88    = [bool(region_in_88[i]) for i in final_indices]
    logger.info(f"Final panel features: {len(panel_features)}  (mode={selection_mode})")
    logger.info(f"  ★ = in 88 tissue regions   ○ = outside")
    for feat, in88 in sorted(zip(panel_features, panel_in_88), key=lambda x: x[0]):
        marker = "★" if in88 else "○"
        logger.info(f"    {marker} {feat}")

    # 5. Build importance DataFrame
    lasso_selected = np.zeros(len(feature_names), dtype=bool)
    lasso_selected[lasso_indices] = True
    rf_selected = np.zeros(len(feature_names), dtype=bool)
    rf_selected[rf_indices] = True
    panel_selected = np.zeros(len(feature_names), dtype=bool)
    panel_selected[final_indices] = True

    importance_df = pd.DataFrame({
        "feature":           feature_names,
        "in_88_regions":     region_in_88,
        "lasso_coef":        lasso_coefs,
        "rf_importance":     rf_importances,
        "sample_pval":       pvals_kept,
        "sample_delta_beta": deltas_kept.round(4),
        "lasso_selected":    lasso_selected,
        "rf_selected":       rf_selected,
        "panel_selected":    panel_selected,
    }).sort_values("rf_importance", ascending=False)

    # 6. Compute depth stats; apply >=10x coverage QC hard filter
    ml_cfg2     = (config or {}).get("ml", {})
    min_10x_cov = ml_cfg2.get("min_10x_coverage_rate", 0.80)

    D_panel         = D[:, final_indices]
    n_samples_total = D_panel.shape[0]

    site_mean    = D_panel.mean(axis=0)
    site_median  = np.median(D_panel, axis=0)
    site_p10     = np.percentile(D_panel, 10, axis=0)
    site_p90     = np.percentile(D_panel, 90, axis=0)
    site_min     = D_panel.min(axis=0)
    site_max     = D_panel.max(axis=0)
    site_n_cov   = (D_panel > 0).sum(axis=0)
    site_n_10x   = (D_panel >= 10).sum(axis=0)
    cov_rate_10x = site_n_10x / n_samples_total

    _depth_available = D_panel.sum() > 0
    if not _depth_available:
        logger.warning(
            "  Depth matrix is all zeros (stale feature_matrix.npz from a previous run). "
            ">=10x coverage filter skipped — delete feature_matrix.npz and re-run to rebuild."
        )
        pass_cov_qc = np.ones(len(panel_features), dtype=bool)
    else:
        pass_cov_qc = cov_rate_10x >= min_10x_cov

    depth_df = pd.DataFrame({
        "feature":           panel_features,
        "n_covered":         site_n_cov.astype(int),
        "coverage_rate":     (site_n_cov / n_samples_total).round(4),
        "n_10x":             site_n_10x.astype(int),
        "coverage_rate_10x": cov_rate_10x.round(4),
        "pass_coverage_qc":  pass_cov_qc,
        "mean_depth":        site_mean.round(1),
        "median_depth":      site_median.round(1),
        "p10_depth":         site_p10.round(1),
        "p90_depth":         site_p90.round(1),
        "min_depth":         site_min.astype(int),
        "max_depth":         site_max.astype(int),
    })

    n_removed_cov = int((~pass_cov_qc).sum())
    if n_removed_cov > 0:
        logger.warning(
            f"  >=10x coverage filter (threshold {min_10x_cov:.0%}): "
            f"removed {n_removed_cov} low-coverage sites"
        )
        fail_rows = depth_df[~pass_cov_qc][["feature", "coverage_rate_10x"]].head(10)
        for _, row in fail_rows.iterrows():
            logger.debug(f"    removed: {row['feature']}  10x_rate={row['coverage_rate_10x']:.2%}")
        panel_features = [f for f, p in zip(panel_features, pass_cov_qc) if p]
    else:
        logger.info(
            f"  >=10x coverage filter (threshold {min_10x_cov:.0%}): "
            f"all {len(panel_features)} sites pass"
        )

    # Update panel_selected to reflect coverage QC
    panel_set = set(panel_features)
    importance_df["panel_selected"] = importance_df["feature"].isin(panel_set)

    # 7. Save feature importance (panel_selected is now final)
    imp_path = os.path.join(result_dir, "feature_importance.csv")
    importance_df.to_csv(imp_path, index=False)
    logger.info(f"  Feature importance saved: {imp_path}")

    # 8. Save panel CpG site list (coverage-QC-filtered)
    panel_df = pd.DataFrame({"feature": panel_features})
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    panel_df.to_csv(panel_path, index=False)
    logger.info(f"  Panel CpG sites saved: {panel_path} ({len(panel_features)} sites after >=10x filter)")

    # 9. Save full feature matrix (pre-panel, for nested CV in classify).
    #    Use X_save / y_save / sample_ids_save so ALL samples (including any
    #    outliers excluded from the selection computation) are present for classifiers.
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    np.savez(fm_path,
             X=X_save, y=y_save,
             sample_ids=np.array(sample_ids_save, dtype=object),
             feature_names=np.array(feature_names_all, dtype=object),
             D=D_save,
             region_in_88=region_in_88_all)
    logger.info(f"  Feature matrix saved ({X_save.shape[0]} samples × "
                f"{X_save.shape[1]} features, pre-Wilcoxon): {fm_path}")

    # 10. Save depth stats
    depth_path = os.path.join(result_dir, "panel_cpg_depth.csv")
    depth_df.to_csv(depth_path, index=False)
    logger.info(f"  Panel depth stats saved: {depth_path}")

    overall_mean   = float(site_mean.mean())
    overall_median = float(np.median(site_median))
    low5_count  = int((site_mean < 5).sum())
    low10_count = int((site_mean < 10).sum())
    ge20_count  = int((site_mean >= 20).sum())
    n_pass_10x  = int(pass_cov_qc.sum())
    logger.info(f"  Panel depth summary: mean={overall_mean:.1f}x  median={overall_median:.1f}x  "
                f"mean<5x: {low5_count}  <10x: {low10_count}  >=20x: {ge20_count}")
    logger.info(f"  Sites with >=10x coverage in >={min_10x_cov:.0%} of samples: "
                f"{n_pass_10x}/{len(depth_df)}")

    return panel_features, importance_df
