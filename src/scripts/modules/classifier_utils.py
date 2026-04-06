#!/usr/bin/env python3
"""
classifier_utils.py - 诊断分类器模块（Step 15）

使用留一法交叉验证（LOO-CV）评估三种分类器：
  LR  — L2 逻辑回归（可解释，适合小样本）
  SVM — RBF 核支持向量机（非线性边界）
  RF  — 随机森林（集成方法，自带特征重要性）

输出：
  roc_curve.png              — 三分类器 ROC 曲线对比图
  classifier_predictions.csv — 每个样本的预测概率和标签
  classifier_summary.csv     — AUC / 准确率 / 灵敏度 / 特异度汇总

依赖：scikit-learn、matplotlib
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
        X, y, sample_ids, feature_names（全量特征，未做 panel 过滤）
    """
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if os.path.exists(fm_path):
        logger.info(f"从缓存加载特征矩阵: {fm_path}")
        data = np.load(fm_path, allow_pickle=True)
        X = data["X"]
        y = data["y"]
        sample_ids = list(data["sample_ids"])
        feature_names = list(data["feature_names"])
    else:
        logger.info("未找到 feature_matrix.npz，从 bedGraph 重新构建...")
        from scripts.modules.feature_utils import build_feature_matrix
        X, y, sample_ids, feature_names = build_feature_matrix(
            workdir, result_dir, benign_samples, malignant_samples, config
        )
    return X, y, sample_ids, feature_names


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


def _save_roc_plot(results, result_dir):
    """绘制并保存三分类器 ROC 曲线对比图。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = {"LR": "#2196F3", "SVM": "#4CAF50", "RF": "#FF5722"}

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
    """为每个分类器绘制混淆矩阵。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_clf = len(results)
    fig, axes = plt.subplots(1, n_clf, figsize=(4 * n_clf, 4))
    if n_clf == 1:
        axes = [axes]

    for ax, (clf_name, res) in zip(axes, results.items()):
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
    Step 15 主函数：LOO-CV 评估 LR / SVM / RF 分类性能。

    Args:
        workdir:           临时目录（含 bedGraph；仅在无缓存时使用）
        result_dir:        结果目录（读取 feature_matrix.npz / panel_cpg.csv，写入输出）
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

    logger.info(f"开始分类器评估  (良性={len(benign_samples)}, 恶性={len(malignant_samples)})")

    # 1. 加载特征矩阵
    X, y, sample_ids, feature_names = _load_feature_matrix(
        workdir, result_dir, benign_samples, malignant_samples, config
    )

    # 2. Panel 过滤
    X, feature_names = _apply_panel_filter(X, feature_names, result_dir)

    if X.shape[1] == 0:
        logger.error("特征矩阵无特征，跳过分类器评估")
        return {}, pd.DataFrame()

    # 3. 插补 NaN（全局中位数，与 LOO 无关）
    X = _impute(X)
    logger.info(f"  特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征")

    # 4. 分类器定义
    ml_cfg = (config or {}).get("ml", {})
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
    }

    # 5. LOO-CV
    loo = LeaveOneOut()
    results = {}

    for clf_name, clf in classifiers.items():
        logger.info(f"  LOO-CV: {clf_name} ...")
        y_true_list, y_prob_list, y_pred_list = [], [], []

        for train_idx, test_idx in loo.split(X):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # 标准化（仅在训练集上 fit，防止数据泄漏）
            scaler = StandardScaler()
            X_train_sc = scaler.fit_transform(X_train)
            X_test_sc = scaler.transform(X_test)

            clf.fit(X_train_sc, y_train)
            prob = clf.predict_proba(X_test_sc)[0, 1]  # P(恶性)
            pred = clf.predict(X_test_sc)[0]

            y_true_list.append(int(y_test[0]))
            y_prob_list.append(float(prob))
            y_pred_list.append(int(pred))

        y_true = np.array(y_true_list)
        y_prob = np.array(y_prob_list)
        y_pred = np.array(y_pred_list)

        auc = roc_auc_score(y_true, y_prob)
        acc = accuracy_score(y_true, y_pred)
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
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

    # 6. 保存 ROC 曲线
    _save_roc_plot(results, result_dir)

    # 7. 保存混淆矩阵图
    _save_confusion_matrix_plot(results, result_dir)

    # 8. 保存逐样本预测结果
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

    # 9. 保存汇总指标
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
            "n_features": X.shape[1],
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "classifier_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  分类器汇总已保存: {summary_path}")

    # 打印汇总到控制台
    logger.info("\n  ┌─────────────────────────────────────────────────────┐")
    logger.info("  │  分类器 LOO-CV 结果汇总                            │")
    logger.info("  ├──────┬───────┬────────┬──────────┬────────────────┤")
    logger.info("  │ 分类器│  AUC  │ 准确率 │  灵敏度  │   特异度       │")
    logger.info("  ├──────┼───────┼────────┼──────────┼────────────────┤")
    for row in summary_rows:
        logger.info(f"  │ {row['classifier']:<4} │ {row['auc']:.3f} │  {row['accuracy']:.3f} │  "
                    f"{row['sensitivity']:.3f}   │   {row['specificity']:.3f}          │")
    logger.info("  └──────┴───────┴────────┴──────────┴────────────────┘")

    return results, summary_df
