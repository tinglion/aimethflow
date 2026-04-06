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
import numpy as np
import pandas as pd

from utils.logger import logger


# ──────────────────────────────────────────────────────
# 特征矩阵构建
# ──────────────────────────────────────────────────────

def build_feature_matrix(workdir, result_dir, benign_samples, malignant_samples, config=None):
    """
    从每个样本的 _CpG.bedGraph 构建特征矩阵。

    特征 = diff_methylation.txt 中的 DMC 位点（行）。
    每个样本在该位点的甲基化率（0~1）作为特征值。
    样本无覆盖的位点填 NaN（后续按列中位数插补）。

    Args:
        workdir:           临时目录（含 {sample}_CpG.bedGraph）
        result_dir:        结果目录（含 diff_methylation.txt）
        benign_samples:    良性样本 ID 列表
        malignant_samples: 恶性样本 ID 列表
        config:            配置对象

    Returns:
        X:            np.ndarray, shape (n_samples, n_features)，meth 值 0-1，NaN 表示无覆盖
        y:            np.ndarray, shape (n_samples,)，0=良性 1=恶性
        sample_ids:   list[str]
        feature_names: list[str]，格式 "chr:pos"
    """
    ml_cfg = (config or {}).get("ml", {})
    min_coverage = ml_cfg.get("min_feature_coverage", 0.5)

    # 1. 读取 DMC 位点
    diff_file = os.path.join(result_dir, "diff_methylation.txt")
    if not os.path.exists(diff_file):
        raise FileNotFoundError(f"差异甲基化文件不存在: {diff_file}")

    diff_df = pd.read_csv(diff_file, sep="\t", header=None,
                          names=["chrom", "pos", "meth1", "meth2", "diff"])
    diff_df["pos"] = diff_df["pos"].astype(int)
    feature_names = [f"{row.chrom}:{row.pos}" for _, row in diff_df.iterrows()]
    n_features = len(feature_names)
    logger.info(f"  DMC位点数: {n_features}")

    # 2. 构建位置索引：(chrom, pos) → 列索引
    pos_index = {(row.chrom, int(row.pos)): i
                 for i, (_, row) in enumerate(diff_df.iterrows())}

    # 3. 样本列表（良性在前，便于后续与标签对齐）
    all_samples = [(sid, 0) for sid in benign_samples] + \
                  [(sid, 1) for sid in malignant_samples]
    n_samples = len(all_samples)

    X = np.full((n_samples, n_features), np.nan)
    y = np.array([label for _, label in all_samples])
    sample_ids = [sid for sid, _ in all_samples]

    # 4. 逐样本填充特征矩阵
    missing_samples = []
    for i, (sample_id, _) in enumerate(all_samples):
        bedgraph = os.path.join(workdir, f"{sample_id}_CpG.bedGraph")
        if not os.path.exists(bedgraph):
            logger.warning(f"  bedGraph不存在，跳过: {bedgraph}")
            missing_samples.append(sample_id)
            continue

        # 读取 bedGraph（chrom, start, end, meth%, methCount, unmethCount）
        bg = pd.read_csv(bedgraph, sep="\t", header=None,
                         names=["chrom", "start", "end", "meth_pct",
                                "meth_count", "unmeth_count"])

        # 构建 (chrom, start) → meth_pct 字典，只保留 DMC 位置
        bg_dict = {}
        for _, row in bg.iterrows():
            key = (row["chrom"], int(row["start"]))
            if key in pos_index:
                bg_dict[key] = row["meth_pct"]

        # 填入特征矩阵（meth% 0-100 → 0-1）
        for (chrom, pos), j in pos_index.items():
            if (chrom, pos) in bg_dict:
                X[i, j] = bg_dict[(chrom, pos)] / 100.0

    if missing_samples:
        logger.warning(f"  {len(missing_samples)} 个样本无bedGraph: {missing_samples}")

    # 5. 按覆盖率过滤低覆盖特征
    coverage_fraction = np.mean(~np.isnan(X), axis=0)
    good_mask = coverage_fraction >= min_coverage
    n_dropped = int(np.sum(~good_mask))
    if n_dropped > 0:
        logger.info(f"  过滤低覆盖特征 (<{min_coverage:.0%} 覆盖率): 丢弃 {n_dropped} 个，"
                    f"保留 {int(np.sum(good_mask))} 个")
        X = X[:, good_mask]
        feature_names = [fn for fn, g in zip(feature_names, good_mask) if g]

    logger.info(f"  特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征  "
                f"(良性={sum(y==0)}, 恶性={sum(y==1)})")
    nan_rate = np.isnan(X).mean()
    logger.info(f"  缺失值比例: {nan_rate:.2%}")

    return X, y, sample_ids, feature_names


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

    # 1. 构建特征矩阵
    X, y, sample_ids, feature_names = build_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    if X.shape[1] == 0:
        logger.error("特征矩阵为空（无有效 DMC 位点），跳过特征选择")
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
        "feature": feature_names,
        "lasso_coef": lasso_coefs,
        "rf_importance": rf_importances,
        "lasso_selected": lasso_selected,
        "rf_selected": rf_selected,
        "panel_selected": panel_selected,
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
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    np.savez(fm_path,
             X=X, y=y,
             sample_ids=np.array(sample_ids, dtype=object),
             feature_names=np.array(feature_names, dtype=object))
    logger.info(f"  特征矩阵已保存: {fm_path}")

    return panel_features, importance_df
