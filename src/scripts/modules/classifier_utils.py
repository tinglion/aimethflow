#!/usr/bin/env python3
"""
classifier_utils.py - 诊断分类器模块（Step 15）

使用留一法交叉验证（LOO-CV）评估五种分类器：
  LR    — L2 逻辑回归（可解释，适合小样本）
  SVM   — RBF 核支持向量机（非线性边界）
  RF    — 随机森林（集成方法，自带特征重要性）
  XGB   — XGBoost（梯度提升树，scale_pos_weight 处理不平衡）
  LGBM  — LightGBM（高效梯度提升，scale_pos_weight 处理不平衡）

输出：
  roc_curve.png              — 分类器 ROC 曲线对比图
  confusion_matrix.png       — 各分类器混淆矩阵
  classifier_predictions.csv — 每个样本的预测概率和标签
  classifier_summary.csv     — AUC / 准确率 / 灵敏度 / 特异度汇总

依赖：scikit-learn、matplotlib、xgboost、lightgbm
"""

import os
import numpy as np
import pandas as pd

from utils.logger import logger


def _load_feature_matrix(workdir, result_dir, benign_samples, malignant_samples, config):
    """
    加载特征矩阵：优先从 feature_matrix.npz 恢复（feature_select 步骤保存的），
    若不存在则从 bedGraph 重新构建（需要 workdir）。

    Returns:
        X, y, sample_ids, feature_names, region_in_88
    """
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if os.path.exists(fm_path):
        logger.info(f"从缓存加载特征矩阵: {fm_path}")
        data = np.load(fm_path, allow_pickle=True)
        X = data["X"]
        y = data["y"]
        sample_ids    = list(data["sample_ids"])
        feature_names = list(data["feature_names"])
        region_in_88  = data["region_in_88"] if "region_in_88" in data else np.zeros(X.shape[1], dtype=bool)
    else:
        logger.info("未找到 feature_matrix.npz，从 bedGraph 重新构建...")
        from scripts.modules.feature_utils import build_feature_matrix
        X, y, sample_ids, feature_names, _, region_in_88 = build_feature_matrix(
            workdir, result_dir, benign_samples, malignant_samples, config
        )
    return X, y, sample_ids, feature_names, region_in_88


def _apply_panel_filter(X, feature_names, result_dir):
    """
    若存在 panel_cpg.csv，将特征矩阵过滤为 panel 位点；否则使用全部特征。

    Returns:
        X_filtered, filtered_feature_names
    """
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if os.path.exists(panel_path):
        panel_df = pd.read_csv(panel_path)
        panel_set = set(panel_df["feature"])
        indices = [i for i, fn in enumerate(feature_names) if fn in panel_set]
        logger.info(f"应用 panel_cpg.csv：使用 {len(indices)} 个特征（共 {len(feature_names)} 个）")
        return X[:, indices], [feature_names[i] for i in indices]
    else:
        logger.info(f"未找到 panel_cpg.csv，使用全部 {len(feature_names)} 个特征")
        return X, feature_names


def _impute(X):
    """按列中位数插补 NaN。"""
    X_imp = X.copy()
    col_medians = np.nanmedian(X_imp, axis=0)
    for j in range(X_imp.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j] if not np.isnan(col_medians[j]) else 0.0
    return X_imp


def _impute_fit(X):
    """
    用训练集列中位数插补 NaN，返回 (插补后矩阵, 列中位数向量)。
    NaN 中位数（全列缺失）兜底为 0.0。
    """
    col_medians = np.nanmedian(X, axis=0)
    col_medians = np.where(np.isnan(col_medians), 0.0, col_medians)
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp, col_medians


def _impute_transform(X, col_medians):
    """用训练集的列中位数插补测试集 NaN。"""
    X_imp = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_medians[j]
    return X_imp


def _select_features_in_fold(X_train, y_train, region_in_88, config):
    """
    在单个 LOO fold 的训练集上运行特征选择，返回选中的列索引。

    流程：分层 Wilcoxon 预筛选（88 区域宽松 / 其余标准）→ LASSO → RF → 并集
    内部自行插补 NaN（用训练集中位数），严格不使用测试样本数据。

    Args:
        X_train:      (n_train, n_features)
        y_train:      (n_train,)
        region_in_88: (n_features,) bool — 是否属于 88 个组织候选区域
        config:       配置对象

    Returns:
        selected_indices: np.ndarray  (相对于传入 X_train 的列索引)
    """
    try:
        from scipy import stats as _stats
        _HAS_SCIPY = True
    except ImportError:
        _HAS_SCIPY = False

    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler

    ml_cfg = (config or {}).get("ml", {})
    wilcox_p        = ml_cfg.get("sample_wilcox_p",       0.05)
    wilcox_delta    = ml_cfg.get("sample_delta_beta",     0.10)
    wilcox_p_88     = ml_cfg.get("sample_wilcox_p_88",    0.10)
    wilcox_delta_88 = ml_cfg.get("sample_delta_beta_88",  0.05)
    selection_mode  = ml_cfg.get("feature_selection_mode", "union")

    n_features = X_train.shape[1]

    # 1. 插补（训练集中位数）
    X_imp, _ = _impute_fit(X_train)

    # 2. 分层 Wilcoxon 预筛选
    if _HAS_SCIPY:
        pvals  = np.ones(n_features)
        deltas = np.zeros(n_features)
        for j in range(n_features):
            b_vals = X_imp[y_train == 0, j]
            m_vals = X_imp[y_train == 1, j]
            if len(b_vals) >= 2 and len(m_vals) >= 2:
                _, p = _stats.mannwhitneyu(b_vals, m_vals, alternative="two-sided")
                d    = abs(float(np.mean(m_vals)) - float(np.mean(b_vals)))
            else:
                p, d = 1.0, 0.0
            pvals[j]  = p
            deltas[j] = d

        p_thresh = np.where(region_in_88, wilcox_p_88,    wilcox_p)
        d_thresh = np.where(region_in_88, wilcox_delta_88, wilcox_delta)
        wilcox_mask = (pvals < p_thresh) & (deltas >= d_thresh)

        if wilcox_mask.sum() == 0:   # 全量放宽兜底
            wilcox_mask = (pvals < 0.20) & (deltas >= 0.05)
        if wilcox_mask.sum() == 0:
            wilcox_mask = np.ones(n_features, dtype=bool)
    else:
        wilcox_mask = np.ones(n_features, dtype=bool)

    wilcox_indices = np.where(wilcox_mask)[0]
    X_w = X_imp[:, wilcox_indices]

    # 3. LASSO（fold 内加速：Cs=10, cv≤3）
    cv_folds = min(ml_cfg.get("lasso_cv_folds", 5), 3)
    scaler = StandardScaler()
    X_ws = scaler.fit_transform(X_w)
    lasso = LogisticRegressionCV(
        Cs=10, cv=cv_folds, penalty="l1", solver="liblinear",
        class_weight="balanced", max_iter=1000, random_state=42,
    )
    lasso.fit(X_ws, y_train)
    lasso_local = np.where(lasso.coef_[0] != 0)[0]

    # 4. RF（100 棵树，fold 内快速估计）
    rf = RandomForestClassifier(
        n_estimators=100, class_weight="balanced",
        random_state=42, n_jobs=-1,
    )
    rf.fit(X_w, y_train)
    rf_local = np.where(rf.feature_importances_ > 0)[0]

    # 5. 并集映射回原始列索引
    if selection_mode == "lasso":
        local_idx = lasso_local
    elif selection_mode == "rf":
        local_idx = rf_local
    else:
        local_idx = np.unique(np.concatenate([lasso_local, rf_local]))

    return wilcox_indices[local_idx]


def _save_roc_plot(results, result_dir):
    """绘制并保存分类器 ROC 曲线对比图。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve
    matplotlib.rcParams.update({"font.family": "DejaVu Sans"})

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {
        "LR":   "#2196F3",
        "SVM":  "#4CAF50",
        "RF":   "#FF5722",
        "XGB":  "#9C27B0",
        "LGBM": "#FF9800",
    }

    for clf_name, res in results.items():
        fpr, tpr, _ = roc_curve(res["y_true"], res["y_prob"])
        ax.plot(fpr, tpr,
                label=f"{clf_name} (AUC={res['auc']:.3f})",
                color=colors.get(clf_name, "gray"), lw=2)

    ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
    ax.set_xlabel("False Positive Rate (1 − Specificity)")
    ax.set_ylabel("True Positive Rate (Sensitivity)")
    ax.set_title("LOO-CV ROC Curves (Benign vs Malignant)")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])

    roc_path = os.path.join(result_dir, "roc_curve.png")
    fig.savefig(roc_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  ROC曲线已保存: {roc_path}")
    return roc_path


def _save_confusion_matrix_plot(results, result_dir):
    """为每个分类器绘制混淆矩阵，超过 3 个时自动换行。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "DejaVu Sans"})

    n_clf = len(results)
    n_cols = min(n_clf, 3)
    n_rows = (n_clf + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4 * n_cols, 4 * n_rows),
                             squeeze=False)
    # 展平为列表，方便按顺序遍历
    axes_flat = [axes[r][c] for r in range(n_rows) for c in range(n_cols)]

    for ax, (clf_name, res) in zip(axes_flat, results.items()):
        cm = res["confusion_matrix"]
        im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
        ax.set_title(f"{clf_name}\nAUC={res['auc']:.3f}")
        tick_labels = ["Benign", "Malignant"]
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(tick_labels)
        ax.set_yticklabels(tick_labels)
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, str(cm[i, j]),
                        ha="center", va="center",
                        color="white" if cm[i, j] > cm.max() / 2 else "black",
                        fontsize=14)

    # 隐藏多余的子图
    for ax in axes_flat[n_clf:]:
        ax.set_visible(False)

    fig.tight_layout()
    cm_path = os.path.join(result_dir, "confusion_matrix.png")
    fig.savefig(cm_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  混淆矩阵已保存: {cm_path}")
    return cm_path


# ──────────────────────────────────────────────────────
# 主入口：Step 15
# ──────────────────────────────────────────────────────

def run_classifiers(workdir, result_dir, config, benign_samples, malignant_samples):
    """
    Step 15 主函数：Nested LOO-CV 评估 LR / SVM / RF / XGB / LGBM 分类性能。

    特征选择（Wilcoxon + LASSO + RF）在每个 LOO fold 的训练集内独立完成，
    测试样本的数据严格不参与任何特征选择或归一化 fit，杜绝数据泄漏。

    Args:
        workdir:           临时目录（含 bedGraph；仅在无缓存时使用）
        result_dir:        结果目录（读取 feature_matrix.npz，写入输出）
        config:            配置对象
        benign_samples:    良性样本 ID 列表
        malignant_samples: 恶性样本 ID 列表

    Returns:
        results:    dict {clf_name: {y_true, y_prob, y_pred, auc, accuracy,
                                    sensitivity, specificity, confusion_matrix}}
        summary_df: pd.DataFrame，各分类器汇总指标
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import (roc_auc_score, confusion_matrix,
                                 accuracy_score)
    from xgboost import XGBClassifier
    from lightgbm import LGBMClassifier

    logger.info(f"开始分类器评估（Nested LOO-CV）  "
                f"(良性={len(benign_samples)}, 恶性={len(malignant_samples)})")

    # 1. 加载全量特征矩阵（覆盖率过滤后、Wilcoxon 前）
    X, y, sample_ids, feature_names, region_in_88 = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    if X.shape[1] == 0:
        logger.error("特征矩阵无特征，跳过分类器评估")
        return {}, pd.DataFrame()

    n_in_88 = int(region_in_88.sum())
    logger.info(f"  全量特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征  "
                f"(88区域内: {n_in_88}  区域外: {X.shape[1]-n_in_88}，"
                f"特征选择将在每折训练集内独立运行)")

    # 2. 分类器定义
    ml_cfg = (config or {}).get("ml", {})

    n_neg = int(sum(y == 0))
    n_pos = int(sum(y == 1))
    spw   = n_neg / n_pos if n_pos > 0 else 1.0

    classifiers = {
        "LR": LogisticRegression(
            C=ml_cfg.get("lr_C", 1.0),
            class_weight="balanced",
            max_iter=2000,
            random_state=42,
        ),
        "SVM": SVC(
            kernel="rbf",
            C=ml_cfg.get("svm_C", 1.0),
            probability=True,
            class_weight="balanced",
            random_state=42,
        ),
        "RF": RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "XGB": XGBClassifier(
            n_estimators=ml_cfg.get("xgb_n_estimators", 200),
            max_depth=ml_cfg.get("xgb_max_depth", 4),
            learning_rate=ml_cfg.get("xgb_learning_rate", 0.1),
            scale_pos_weight=spw,
            use_label_encoder=False,
            eval_metric="logloss",
            random_state=42,
            n_jobs=-1,
            verbosity=0,
        ),
        "LGBM": LGBMClassifier(
            n_estimators=ml_cfg.get("lgbm_n_estimators", 200),
            max_depth=ml_cfg.get("lgbm_max_depth", 4),
            learning_rate=ml_cfg.get("lgbm_learning_rate", 0.1),
            scale_pos_weight=spw,
            random_state=42,
            n_jobs=-1,
            verbosity=-1,
        ),
    }

    # 3. Nested LOO-CV
    loo = LeaveOneOut()
    clf_names = list(classifiers.keys())

    # 每分类器独立收集预测
    clf_y_true = {n: [] for n in clf_names}
    clf_y_prob = {n: [] for n in clf_names}
    clf_y_pred = {n: [] for n in clf_names}

    fold_n_features = []   # 每折选中的特征数

    n_folds = X.shape[0]
    log_interval = max(1, n_folds // 10)

    for fold_idx, (train_idx, test_idx) in enumerate(loo.split(X)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test         = y[train_idx], y[test_idx]

        # 3a. 在训练集上独立选特征（分层 Wilcoxon + LASSO + RF）
        sel_idx = _select_features_in_fold(X_train_raw, y_train, region_in_88, config)
        if len(sel_idx) == 0:
            sel_idx = np.arange(X.shape[1])   # 兜底：全特征
        fold_n_features.append(len(sel_idx))

        X_tr = X_train_raw[:, sel_idx]
        X_te = X_test_raw[:, sel_idx]

        # 3b. 插补（训练集中位数）
        X_tr_imp, col_meds = _impute_fit(X_tr)
        X_te_imp = _impute_transform(X_te, col_meds)

        # 3c. 标准化（在训练集上 fit）
        scaler = StandardScaler()
        X_tr_sc = scaler.fit_transform(X_tr_imp)
        X_te_sc = scaler.transform(X_te_imp)

        # 3d. 所有分类器在同一折特征集上训练 + 预测
        for clf_name, clf in classifiers.items():
            clf.fit(X_tr_sc, y_train)
            prob = clf.predict_proba(X_te_sc)[0, 1]
            pred = clf.predict(X_te_sc)[0]
            clf_y_true[clf_name].append(int(y_test[0]))
            clf_y_prob[clf_name].append(float(prob))
            clf_y_pred[clf_name].append(int(pred))

        if (fold_idx + 1) % log_interval == 0 or fold_idx == n_folds - 1:
            logger.info(f"    fold {fold_idx+1}/{n_folds}  "
                        f"特征数={len(sel_idx)}")

    # 每折特征数统计
    fold_arr = np.array(fold_n_features)
    logger.info(f"  Nested CV 每折特征数: 均值={fold_arr.mean():.1f}  "
                f"中位数={float(np.median(fold_arr)):.1f}  "
                f"最少={fold_arr.min()}  最多={fold_arr.max()}")

    # 4. 汇总指标
    results = {}
    for clf_name in clf_names:
        y_true = np.array(clf_y_true[clf_name])
        y_prob = np.array(clf_y_prob[clf_name])
        y_pred = np.array(clf_y_pred[clf_name])

        auc = roc_auc_score(y_true, y_prob)
        acc = accuracy_score(y_true, y_pred)
        cm  = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

        results[clf_name] = {
            "y_true": y_true,
            "y_prob": y_prob,
            "y_pred": y_pred,
            "auc": auc,
            "accuracy": acc,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "confusion_matrix": cm,
        }
        logger.info(f"    {clf_name}: AUC={auc:.3f}  Acc={acc:.3f}  "
                    f"Sens={sensitivity:.3f}  Spec={specificity:.3f}")

    # 5. 保存 ROC 曲线
    _save_roc_plot(results, result_dir)

    # 6. 保存混淆矩阵图
    _save_confusion_matrix_plot(results, result_dir)

    # 7. 保存逐样本预测结果
    rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(sample_ids):
            rows.append({
                "sample_id": sid,
                "true_label": "malignant" if res["y_true"][i] == 1 else "benign",
                "classifier": clf_name,
                "prob_malignant": round(float(res["y_prob"][i]), 4),
                "predicted": "malignant" if res["y_pred"][i] == 1 else "benign",
                "correct": bool(res["y_true"][i] == res["y_pred"][i]),
            })
    pred_df = pd.DataFrame(rows)
    pred_path = os.path.join(result_dir, "classifier_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"  逐样本预测结果已保存: {pred_path}")

    # 8. 保存汇总指标
    summary_rows = []
    for clf_name, res in results.items():
        cm = res["confusion_matrix"]
        tn, fp, fn, tp = cm.ravel()
        summary_rows.append({
            "classifier": clf_name,
            "auc": round(res["auc"], 4),
            "accuracy": round(res["accuracy"], 4),
            "sensitivity": round(res["sensitivity"], 4),
            "specificity": round(res["specificity"], 4),
            "TP": int(tp), "FP": int(fp),
            "TN": int(tn), "FN": int(fn),
            "n_benign": int(sum(y == 0)),
            "n_malignant": int(sum(y == 1)),
            "n_features_mean": round(float(fold_arr.mean()), 1),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "classifier_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  分类器汇总已保存: {summary_path}")

    # 打印汇总
    logger.info("\n  ┌──────────────────────────────────────────────────────┐")
    logger.info("  │  分类器 Nested LOO-CV 结果汇总                       │")
    logger.info("  ├────────┬───────┬────────┬──────────┬───────────────┤")
    logger.info("  │ 分类器 │  AUC  │ 准确率 │  灵敏度  │   特异度      │")
    logger.info("  ├────────┼───────┼────────┼──────────┼───────────────┤")
    for row in summary_rows:
        logger.info(f"  │ {row['classifier']:<6} │ {row['auc']:.3f} │  {row['accuracy']:.3f} │  "
                    f"{row['sensitivity']:.3f}   │   {row['specificity']:.3f}         │")
    logger.info("  └────────┴───────┴────────┴──────────┴───────────────┘")

    return results, summary_df
