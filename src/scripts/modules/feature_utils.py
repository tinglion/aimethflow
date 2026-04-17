#!/usr/bin/env python3
"""
feature_utils.py - 特征选择模块（Step 14）

从每个样本的 _CpG.bedGraph 文件构建特征矩阵，
使用 LASSO（L1 logistic）+ 随机森林重要性进行特征筛选，
输出 panel_cpg.csv（入选位点列表）和 feature_importance.csv。

依赖：scikit-learn、numpy、pandas
安装：pip install scikit-learn matplotlib seaborn
"""

import os
import bisect
import numpy as np
import pandas as pd
from collections import defaultdict

from utils.logger import logger


# ──────────────────────────────────────────────────────
# 区域归属检查
# ──────────────────────────────────────────────────────

def _check_region_membership(feature_names, candidate_regions_file):
    """
    检查每个 CpG 位点（格式 "chr:pos"）是否落在 88 个组织候选区域内。

    使用 bisect 二分查找，复杂度 O(n_features × log(n_regions_per_chrom))。

    Returns:
        np.ndarray(bool), shape (n_features,)
    """
    in_region = np.zeros(len(feature_names), dtype=bool)

    if not candidate_regions_file or not os.path.exists(candidate_regions_file):
        return in_region

    # 按染色体构建排序区间列表
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

    # 为每条染色体提取起始坐标数组（用于 bisect）
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
# 特征矩阵构建
# ──────────────────────────────────────────────────────

def build_feature_matrix(workdir, result_dir, benign_samples, malignant_samples, config=None):
    """
    从每个样本的 _CpG.bedGraph 构建特征矩阵及深度矩阵。

    候选位点由数据驱动发现：在 ≥ min_feature_coverage 比例的样本中有深度覆盖的 CpG。
    不再依赖 diff_methylation.txt，直接从个体 bedGraph 聚合候选集。
    同时标注每个位点是否属于 88 个组织候选区域（用于后续分层 Wilcoxon 筛选）。

    Args:
        workdir:           临时目录（含 {sample}_CpG.bedGraph）
        result_dir:        结果目录（含 diff_methylation.txt，仅供参考/日志）
        benign_samples:    良性样本 ID 列表
        malignant_samples: 恶性样本 ID 列表
        config:            配置对象

    Returns:
        X:             np.ndarray (n_samples, n_features)，meth 值 0-1，NaN 表示无覆盖
        y:             np.ndarray (n_samples,)，0=良性 1=恶性
        sample_ids:    list[str]
        feature_names: list[str]，格式 "chr:pos"
        D:             np.ndarray (n_samples, n_features)，深度，0 表示无覆盖
        region_in_88:  np.ndarray(bool) (n_features,)，是否属于 88 个组织候选区域
    """
    ml_cfg = (config or {}).get("ml", {})
    min_coverage = ml_cfg.get("min_feature_coverage", 0.5)

    all_samples = [(sid, 0) for sid in benign_samples] + \
                  [(sid, 1) for sid in malignant_samples]
    n_samples = len(all_samples)

    # ── Pass 1: 读取所有样本的 bedGraph，聚合位点覆盖计数 ──
    # pos_count: (chrom, pos) → 覆盖该位点的样本数
    # raw_data: sample_id → {(chrom, pos): (meth_pct, depth)}
    pos_count = {}
    raw_data = {}

    missing_samples = []
    for sample_id, _ in all_samples:
        bedgraph = os.path.join(workdir, f"{sample_id}_CpG.bedGraph")
        if not os.path.exists(bedgraph):
            logger.warning(f"  bedGraph不存在，跳过: {bedgraph}")
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
        logger.warning(f"  {len(missing_samples)} 个样本无bedGraph: {missing_samples}")

    # ── 按覆盖率过滤候选位点 ──
    min_count = max(1, int(round(n_samples * min_coverage)))
    candidate_positions = sorted(
        [(chrom, pos) for (chrom, pos), cnt in pos_count.items() if cnt >= min_count],
        key=lambda x: (x[0], x[1])
    )
    n_features = len(candidate_positions)

    # 与 diff_methylation.txt 对比（仅日志，不影响候选集）
    diff_file = os.path.join(result_dir, "diff_methylation.txt")
    if os.path.exists(diff_file):
        try:
            _ddf = pd.read_csv(diff_file, sep="\t", header=None,
                               names=["chrom", "pos", "m1", "m2", "diff"])
            _ddf = _ddf[pd.to_numeric(_ddf["pos"], errors="coerce").notna()]
            logger.info(f"  数据驱动候选位点: 总覆盖 {len(pos_count)} 个，"
                        f"≥{min_coverage:.0%}样本覆盖 → {n_features} 个  "
                        f"(diff_methylation.txt DMC: {len(_ddf)} 个，仅供参考)")
        except Exception:
            logger.info(f"  数据驱动候选位点: 总覆盖 {len(pos_count)} 个，"
                        f"≥{min_coverage:.0%}样本覆盖 → {n_features} 个")
    else:
        logger.info(f"  数据驱动候选位点: 总覆盖 {len(pos_count)} 个，"
                    f"≥{min_coverage:.0%}样本覆盖 → {n_features} 个")

    if n_features == 0:
        y = np.array([label for _, label in all_samples])
        sample_ids = [sid for sid, _ in all_samples]
        return (np.zeros((n_samples, 0)), y, sample_ids, [], np.zeros((n_samples, 0), dtype=np.int32),
                np.zeros(0, dtype=bool))

    feature_names = [f"{chrom}:{pos}" for chrom, pos in candidate_positions]
    pos_index = {cp: i for i, cp in enumerate(candidate_positions)}

    # ── Pass 2: 填充特征矩阵 ──
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

    # ── 88 区域归属标注 ──
    candidate_regions_file = (config or {}).get("data", {}).get("candidate_regions", "")
    region_in_88 = _check_region_membership(feature_names, candidate_regions_file)
    n_in_88  = int(region_in_88.sum())
    n_out_88 = n_features - n_in_88
    logger.info(f"  88 组织区域内: {n_in_88} 个  区域外: {n_out_88} 个  "
                f"(共 {n_features} 个候选位点)")

    logger.info(f"  特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征  "
                f"(良性={int(sum(y==0))}, 恶性={int(sum(y==1))})")
    logger.info(f"  缺失值比例: {np.isnan(X).mean():.2%}")

    return X, y, sample_ids, feature_names, D, region_in_88


def _impute(X):
    """按列中位数插补 NaN。"""
    X_imp = X.copy()
    col_medians = np.nanmedian(X_imp, axis=0)
    for j in range(X_imp.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j] if not np.isnan(col_medians[j]) else 0.0
    return X_imp


# ──────────────────────────────────────────────────────
# LASSO 特征选择（L1 正则化逻辑回归）
# ──────────────────────────────────────────────────────

def run_lasso_selection(X, y, config=None):
    """
    使用带 L1 惩罚的逻辑回归（LogisticRegressionCV）选择特征。

    先插补 NaN，再标准化，再用交叉验证选择最优正则化强度。

    Returns:
        selected_indices: np.ndarray，非零系数的特征列索引
        coefs:            np.ndarray，全部特征的系数（零=未选中）
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

    coefs = clf.coef_[0]  # shape: (n_features,)
    selected_indices = np.where(coefs != 0)[0]
    logger.info(f"  LASSO 选择: {len(selected_indices)}/{len(coefs)} 个特征  "
                f"(最优C={clf.C_[0]:.4f})")
    return selected_indices, coefs


# ──────────────────────────────────────────────────────
# 随机森林特征重要性
# ──────────────────────────────────────────────────────

def run_rf_importance(X, y, config=None):
    """
    使用随机森林计算特征重要性，返回重要性 > min_importance 的特征索引。

    Returns:
        selected_indices: np.ndarray，重要性达标的特征列索引
        importances:      np.ndarray，全部特征的重要性分数
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
    logger.info(f"  RF 选择: {len(selected_indices)}/{len(importances)} 个特征  "
                f"(min_importance={min_importance})")
    return selected_indices, importances


# ──────────────────────────────────────────────────────
# 主入口：Step 14
# ──────────────────────────────────────────────────────

def select_features(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Step 14 主函数：从 bedGraph + diff_methylation.txt 构建特征矩阵，
    运行 LASSO + RF 特征选择，输出结果文件。

    输出（保存到 result_dir）：
        feature_importance.csv — 每个 DMC 位点的 LASSO 系数 + RF 重要性
        panel_cpg.csv          — 最终入选的特征位点列表
        feature_matrix.npz     — 特征矩阵（供 classify 步骤复用）

    Args:
        workdir:           临时目录
        result_dir:        结果目录
        config:            配置对象
        benign_samples:    良性样本 ID 列表
        malignant_samples: 恶性样本 ID 列表

    Returns:
        panel_features: list[str]，入选的特征名列表
        importance_df:  pd.DataFrame，特征重要性表
    """
    logger.info(f"开始特征选择  (良性={len(benign_samples)}, 恶性={len(malignant_samples)})")

    # 1. 构建特征矩阵（含深度矩阵 D 和 88 区域归属）
    X, y, sample_ids, feature_names, D, region_in_88 = build_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    if X.shape[1] == 0:
        logger.error("特征矩阵为空（无有效 DMC 位点），跳过特征选择")
        return [], pd.DataFrame()

    # 1.5. 样本级 Wilcoxon 预筛选
    #   diff_methylation.txt 是对合并 bedGraph 做的比较，组间均值差异虚高（组内噪声被平均掉）。
    #   这一步在个体样本层面做 Mann-Whitney U 检验，去掉在样本级实际无区分度的位点，
    #   只保留同时满足 p<sample_wilcox_p 且 样本级 |Δβ|>sample_delta_beta 的位点。
    ml_cfg_pre = (config or {}).get("ml", {})
    wilcox_p     = ml_cfg_pre.get("sample_wilcox_p",     0.05)
    wilcox_delta = ml_cfg_pre.get("sample_delta_beta",   0.10)

    # 保存全量 Wilcoxon 统计（所有进入筛选的位点，供报告诊断用）
    feature_names_all = list(feature_names)   # 筛选前的完整列表
    X_all = X.copy()   # 覆盖率过滤后、Wilcoxon 前的完整特征矩阵（供 nested CV 使用）
    D_all = D.copy()
    region_in_88_all = region_in_88.copy()   # 与 X_all/feature_names_all 对齐

    ml_cfg_pre = (config or {}).get("ml", {})
    wilcox_p       = ml_cfg_pre.get("sample_wilcox_p",      0.05)   # 非 88 区域标准阈值
    wilcox_delta   = ml_cfg_pre.get("sample_delta_beta",    0.10)
    wilcox_p_88    = ml_cfg_pre.get("sample_wilcox_p_88",   0.10)   # 88 区域宽松阈值
    wilcox_delta_88 = ml_cfg_pre.get("sample_delta_beta_88", 0.05)

    try:
        from scipy import stats as _stats
        _HAS_SCIPY = True
    except ImportError:
        _HAS_SCIPY = False
        logger.warning("scipy 未安装，跳过样本级 Wilcoxon 预筛选")

    if _HAS_SCIPY:
        X_imp_pre = _impute(X)   # 用列中位数插补，仅用于统计检验
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

        # ── 分层 Wilcoxon 阈值 ──
        # 88 组织区域位点：宽松（先验知识保底）
        # 其余位点：标准阈值（数据驱动发现）
        p_thresh = np.where(region_in_88, wilcox_p_88,    wilcox_p)
        d_thresh = np.where(region_in_88, wilcox_delta_88, wilcox_delta)
        wilcox_mask = (pvals < p_thresh) & (deltas >= d_thresh)

        n_pass    = int(wilcox_mask.sum())
        n_in88    = int(region_in_88.sum())
        n_pass_88 = int((wilcox_mask & region_in_88).sum())
        n_pass_ot = int((wilcox_mask & ~region_in_88).sum())

        logger.info(f"样本级 Wilcoxon 分层筛选结果: {n_pass}/{X.shape[1]} 通过")
        logger.info(f"  88区域内 ({n_in88}个, p<{wilcox_p_88}, |Δβ|≥{wilcox_delta_88}): "
                    f"{n_pass_88} 通过")
        logger.info(f"  88区域外 ({X.shape[1]-n_in88}个, p<{wilcox_p}, |Δβ|≥{wilcox_delta}): "
                    f"{n_pass_ot} 通过")
        logger.info(f"  样本级 |Δβ| — 均值={deltas.mean():.3f}  "
                    f"中位数={float(np.median(deltas)):.3f}  最大={deltas.max():.3f}")
        logger.info(f"  p 值        — 中位数={float(np.median(pvals)):.3f}  "
                    f"<0.05: {int((pvals<0.05).sum())}  <0.10: {int((pvals<0.10).sum())}  "
                    f"<0.20: {int((pvals<0.20).sum())}")

        # 边界位点（接近本档阈值但未通过）
        near_mask = (~wilcox_mask) & (
            (deltas >= d_thresh * 0.8) | (pvals < p_thresh * 2.0)
        )
        if near_mask.any():
            bl_names  = [feature_names_all[j] for j in np.where(near_mask)[0]]
            bl_deltas = deltas[near_mask]
            bl_pvals  = pvals[near_mask]
            top_bl = sorted(zip(bl_names, bl_deltas, bl_pvals), key=lambda x: -x[1])[:5]
            logger.info("  接近阈值（未通过）位点示例:")
            for nm, dd, pp in top_bl:
                logger.info(f"    {nm}  |Δβ|={dd:.3f}  p={pp:.4f}")

        # 保存全量 Wilcoxon 统计表（含区域归属）
        wilcox_df = pd.DataFrame({
            "feature":        feature_names_all,
            "in_88_regions":  region_in_88_all,
            "benign_mean":    b_means.round(4),
            "malignant_mean": m_means.round(4),
            "sample_delta":   deltas.round(4),
            "wilcox_pval":    pvals.round(6),
            "p_threshold":    p_thresh.round(4),
            "delta_threshold": d_thresh.round(4),
            "pass_wilcox":    wilcox_mask,
        }).sort_values("sample_delta", ascending=False)
        wpath = os.path.join(result_dir, "wilcoxon_stats.csv")
        wilcox_df.to_csv(wpath, index=False)
        logger.info(f"  全量 Wilcoxon 统计已保存: {wpath}")

        if n_pass == 0:
            logger.warning("分层 Wilcoxon 后无位点通过；全量放宽至 p<0.20 / |Δβ|≥0.05 重试")
            wilcox_mask = (pvals < 0.20) & (deltas >= 0.05)
            n_pass = int(wilcox_mask.sum())
            logger.info(f"  放宽后: {n_pass} 个位点通过")

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
        logger.error("Wilcoxon 筛选后特征矩阵为空，跳过特征选择")
        return [], pd.DataFrame()

    # 2. LASSO 特征选择
    logger.info("运行 LASSO (L1 logistic) 特征选择...")
    lasso_indices, lasso_coefs = run_lasso_selection(X, y, config)

    # 3. RF 特征重要性
    logger.info("运行随机森林特征重要性分析...")
    rf_indices, rf_importances = run_rf_importance(X, y, config)

    # 4. 合并特征选择结果
    ml_cfg = (config or {}).get("ml", {})
    selection_mode = ml_cfg.get("feature_selection_mode", "union")

    if selection_mode == "lasso":
        final_indices = lasso_indices
    elif selection_mode == "rf":
        final_indices = rf_indices
    else:  # "union"（默认）
        final_indices = np.unique(np.concatenate([lasso_indices, rf_indices]))

    panel_features = [feature_names[i] for i in final_indices]
    logger.info(f"最终入选特征: {len(panel_features)} 个  (mode={selection_mode})")

    # 5. 保存特征重要性
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

    imp_path = os.path.join(result_dir, "feature_importance.csv")
    importance_df.to_csv(imp_path, index=False)
    logger.info(f"  特征重要性已保存: {imp_path}")

    # 6. 保存 panel CpG 位点列表
    panel_df = pd.DataFrame({"feature": panel_features})
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    panel_df.to_csv(panel_path, index=False)
    logger.info(f"  Panel CpG位点已保存: {panel_path}")

    # 7. 保存特征矩阵（供 classify 步骤复用，避免重读 bedGraph）
    #    存入 Wilcoxon 筛选前的全量特征（覆盖率过滤后），让 nested CV 在每折内独立筛选。
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    np.savez(fm_path,
             X=X_all, y=y,
             sample_ids=np.array(sample_ids, dtype=object),
             feature_names=np.array(feature_names_all, dtype=object),
             D=D_all,
             region_in_88=region_in_88_all)
    logger.info(f"  特征矩阵已保存 (全量 {X_all.shape[1]} 特征，Wilcoxon 前，含区域归属): {fm_path}")

    # 8. 计算 panel 位点深度统计并保存
    D_panel = D[:, final_indices]            # shape: (n_samples, n_panel)
    n_samples_total = D_panel.shape[0]

    site_mean   = D_panel.mean(axis=0)
    site_median = np.median(D_panel, axis=0)
    site_p10    = np.percentile(D_panel, 10, axis=0)
    site_p90    = np.percentile(D_panel, 90, axis=0)
    site_min    = D_panel.min(axis=0)
    site_max    = D_panel.max(axis=0)
    site_n_cov  = (D_panel > 0).sum(axis=0)  # samples with depth > 0

    depth_df = pd.DataFrame({
        "feature":       panel_features,
        "n_covered":     site_n_cov.astype(int),
        "coverage_rate": (site_n_cov / n_samples_total).round(4),
        "mean_depth":    site_mean.round(1),
        "median_depth":  site_median.round(1),
        "p10_depth":     site_p10.round(1),
        "p90_depth":     site_p90.round(1),
        "min_depth":     site_min.astype(int),
        "max_depth":     site_max.astype(int),
    })
    depth_path = os.path.join(result_dir, "panel_cpg_depth.csv")
    depth_df.to_csv(depth_path, index=False)
    logger.info(f"  Panel位点深度统计已保存: {depth_path}")

    # 深度摘要日志
    overall_mean   = float(site_mean.mean())
    overall_median = float(np.median(site_median))
    low5_count  = int((site_mean < 5).sum())
    low10_count = int((site_mean < 10).sum())
    ge20_count  = int((site_mean >= 20).sum())
    logger.info(f"  Panel深度摘要: 均值={overall_mean:.1f}x  中位数={overall_median:.1f}x  "
                f"平均深度<5x: {low5_count}个  <10x: {low10_count}个  ≥20x: {ge20_count}个")

    return panel_features, importance_df
