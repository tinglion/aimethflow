#!/usr/bin/env python3
"""
subtype_utils.py - 亚型分类器模块（Step 16）

对分类为"恶性"或"交界性"的样本，按"亚型"列进行多分类。
支持的亚型：高级别浆液性癌、子宫内膜样癌、透明细胞癌、低级别浆液性癌、黏液性

特征矩阵：优先复用 feature_matrix.npz（Step 14 缓存），不重读 bedGraph。

分类器（均使用 LOO-CV，小样本标准做法）：
  RF  — 随机森林（n_estimators=200）
        原生支持多分类；class_weight="balanced" 应对类别不均衡；
        特征重要性可解释。
  SVM — RBF 核支持向量机，One-vs-Rest（OvR）策略
        高维小样本表现稳健；probability=True（Platt scaling）输出概率；
        class_weight="balanced" 应对类别不均衡。

输出（保存至 result_dir）：
    subtype_predictions.csv       — 每个样本的预测亚型和各类概率
    subtype_summary.csv           — 每个分类器的准确率汇总
    subtype_confusion_{RF|SVM}.png — 混淆矩阵图
    subtype_roc.png               — OvR ROC 曲线（各亚型 AUC）
"""

import os
import numpy as np
import pandas as pd

from utils.logger import logger

# 亚型标签：Excel 中"亚型"列的规范写法
SUBTYPE_CLASSES = [
    "高级别浆液性癌",
    "子宫内膜样癌",
    "透明细胞癌",
    "低级别浆液性癌",
    "黏液性",
]


# ──────────────────────────────────────────────────────
# 标签读取：从 Excel 构建 sample_id → 亚型 映射
# ──────────────────────────────────────────────────────

def _build_subtype_map(config):
    """
    从 config.toml sample_info_excel 读取所有 sheet，
    返回 {sample_id: subtype_str} 字典。
    仅保留分类为"恶性"或"交界性"且亚型非空的行。
    """
    excel_file = config.get("data", {}).get("sample_info_excel", "")
    if not excel_file or not os.path.exists(excel_file):
        logger.error(f"sample_info_excel 不存在: {excel_file}")
        return {}

    xl = pd.ExcelFile(excel_file)
    subtype_map = {}
    for sheet in xl.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet)
        logger.info(f"  sheet '{sheet}' 列名: {list(df.columns)}")
        if "fq_R1" not in df.columns:
            logger.info(f"  → 无 fq_R1 列，跳过")
            continue

        if "分类" in df.columns:
            df["分类"] = df["分类"].astype(str).str.strip()
            cats = df["分类"].value_counts().to_dict()
            logger.info(f"  → 分类分布: {cats}")
        else:
            logger.warning(f"  → 无 '分类' 列")

        if "亚型" in df.columns:
            df["亚型"] = df["亚型"].astype(str).str.strip()
            subtypes = df["亚型"].value_counts().to_dict()
            logger.info(f"  → 亚型分布: {subtypes}")
        else:
            logger.warning(f"  → 无 '亚型' 列，跳过")
            continue

        for _, row in df.iterrows():
            category = str(row.get("分类", "")).strip() if "分类" in df.columns else "恶性"
            if category not in ("恶性", "交界性"):
                continue
            subtype = str(row.get("亚型", "")).strip()
            if not subtype or subtype in ("nan", ""):
                continue
            r1 = str(row.get("fq_R1", "")).strip()
            if not r1 or r1 == "nan":
                continue
            pid = r1.split(".")[0]
            subtype_map[pid] = subtype

    logger.info(f"从Excel读取亚型标签: {len(subtype_map)} 个样本")
    _dist = {}
    for v in subtype_map.values():
        _dist[v] = _dist.get(v, 0) + 1
    for k, cnt in sorted(_dist.items(), key=lambda x: -x[1]):
        logger.info(f"  {k}: {cnt}")
    return subtype_map


# ──────────────────────────────────────────────────────
# 特征矩阵加载（复用 Step 14 缓存）
# ──────────────────────────────────────────────────────

def _load_feature_matrix(result_dir):
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    if not os.path.exists(fm_path):
        raise FileNotFoundError(
            f"未找到 feature_matrix.npz，请先运行 feature_select 步骤: {fm_path}"
        )
    data = np.load(fm_path, allow_pickle=True)
    X = data["X"]
    sample_ids = list(data["sample_ids"])
    feature_names = list(data["feature_names"])
    logger.info(f"加载特征矩阵: {X.shape[0]} 样本 × {X.shape[1]} 特征")
    return X, sample_ids, feature_names


def _apply_panel_filter(X, feature_names, result_dir):
    panel_path = os.path.join(result_dir, "panel_cpg.csv")
    if os.path.exists(panel_path):
        panel_set = set(pd.read_csv(panel_path)["feature"])
        indices = [i for i, fn in enumerate(feature_names) if fn in panel_set]
        logger.info(f"应用 panel_cpg.csv：使用 {len(indices)} 个特征")
        return X[:, indices], [feature_names[i] for i in indices]
    return X, feature_names


def _impute(X):
    X_imp = X.copy()
    col_med = np.nanmedian(X_imp, axis=0)
    for j in range(X_imp.shape[1]):
        mask = np.isnan(X_imp[:, j])
        X_imp[mask, j] = col_med[j] if not np.isnan(col_med[j]) else 0.0
    return X_imp


# ──────────────────────────────────────────────────────
# 可视化
# ──────────────────────────────────────────────────────

def _save_confusion_plot(cm, classes, clf_name, result_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(6, len(classes)), max(5, len(classes))))
    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    ax.set_title(f"{clf_name} — 亚型混淆矩阵")
    tick_labels = [c[:4] for c in classes]  # 缩短中文标签
    ax.set_xticks(range(len(classes)))
    ax.set_yticks(range(len(classes)))
    ax.set_xticklabels(tick_labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(tick_labels, fontsize=9)
    ax.set_xlabel("预测亚型")
    ax.set_ylabel("真实亚型")
    for i in range(len(classes)):
        for j in range(len(classes)):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black", fontsize=11)
    fig.tight_layout()
    path = os.path.join(result_dir, f"subtype_confusion_{clf_name}.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  混淆矩阵已保存: {path}")


def _save_roc_plot(results, classes, result_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, auc
    from sklearn.preprocessing import label_binarize

    fig, axes = plt.subplots(1, len(results), figsize=(6 * len(results), 5))
    if len(results) == 1:
        axes = [axes]

    colors = plt.cm.tab10.colors

    for ax, (clf_name, res) in zip(axes, results.items()):
        y_true_bin = label_binarize(res["y_true"], classes=list(range(len(classes))))
        y_prob = res["y_prob"]  # (n_samples, n_classes)
        for i, cls in enumerate(classes):
            if y_true_bin.shape[1] <= i:
                continue
            fpr, tpr, _ = roc_curve(y_true_bin[:, i], y_prob[:, i])
            roc_auc = auc(fpr, tpr)
            ax.plot(fpr, tpr, color=colors[i % len(colors)], lw=2,
                    label=f"{cls[:4]} (AUC={roc_auc:.2f})")
        ax.plot([0, 1], [0, 1], "k--", lw=1, alpha=0.5)
        ax.set_title(f"{clf_name} OvR ROC")
        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.legend(loc="lower right", fontsize=8)
        ax.grid(alpha=0.3)

    fig.tight_layout()
    path = os.path.join(result_dir, "subtype_roc.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  ROC曲线已保存: {path}")


# ──────────────────────────────────────────────────────
# 主入口：Step 16
# ──────────────────────────────────────────────────────

def run_subtype_classifier(workdir, result_dir, config):
    """
    Step 16 主函数：对恶性+交界性样本按亚型做多分类 LOO-CV。

    Args:
        workdir:    临时目录（feature_matrix.npz 路径的 result_dir 优先）
        result_dir: 结果目录（读取特征矩阵，写入输出）
        config:     配置对象

    Returns:
        results:    dict {clf_name: {y_true, y_prob, y_pred, accuracy, confusion_matrix}}
        summary_df: pd.DataFrame
    """
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.model_selection import LeaveOneOut
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

    # 1. 构建亚型标签映射
    subtype_map = _build_subtype_map(config)
    if not subtype_map:
        logger.error("无法读取亚型标签，跳过亚型分类")
        return {}, pd.DataFrame()

    # 2. 加载特征矩阵
    try:
        X_all, all_sample_ids, feature_names = _load_feature_matrix(result_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        return {}, pd.DataFrame()

    X_all, feature_names = _apply_panel_filter(X_all, feature_names, result_dir)

    import re  # noqa: F401 (re 已在模块顶层可用，此处保留以防直接调用)

    # 3. 筛选有亚型标签的样本（精确匹配 sample_id）
    indices, labels, sample_ids = [], [], []
    missing = []
    for i, sid in enumerate(all_sample_ids):
        subtype = subtype_map.get(sid)
        if subtype is not None:
            if subtype in SUBTYPE_CLASSES:
                indices.append(i)
                labels.append(subtype)
                sample_ids.append(sid)
            else:
                logger.warning(f"  样本 {sid} 亚型 '{subtype}' 不在已知列表，跳过")
        else:
            missing.append(sid)

    if missing:
        logger.info(f"  {len(missing)} 个样本无亚型标签（可能为良性），已跳过")
    # 诊断：打印前5个 sample_id 与 subtype_map 前5个 key 对比
    logger.info(f"  feature_matrix sample_ids 前5个: {all_sample_ids[:5]}")
    logger.info(f"  subtype_map keys 前5个: {list(subtype_map.keys())[:5]}")

    if len(indices) < 5:
        logger.error(f"有效亚型样本不足（{len(indices)} 个），无法训练分类器")
        return {}, pd.DataFrame()

    X = _impute(X_all[indices])
    le = LabelEncoder()
    le.fit(SUBTYPE_CLASSES)
    y = le.transform(labels)
    classes = le.classes_

    logger.info(f"亚型分类样本: {len(sample_ids)} 个，亚型数: {len(set(labels))}")
    for cls in classes:
        cnt = labels.count(cls)
        if cnt > 0:
            logger.info(f"  {cls}: {cnt}")

    # 4. 分类器（多分类，LOO-CV）
    # RF : 随机森林，原生支持多分类，class_weight="balanced" 处理类别不均衡
    # SVM: RBF 核支持向量机，decision_function_shape="ovr"（一对多）策略，probability=True 输出概率
    ml_cfg = (config or {}).get("ml", {})
    classifiers = {
        "RF": RandomForestClassifier(
            n_estimators=ml_cfg.get("rf_n_estimators", 200),
            class_weight="balanced",
            random_state=42,
            n_jobs=-1,
        ),
        "SVM": SVC(
            kernel="rbf",
            C=ml_cfg.get("svm_C", 1.0),
            probability=True,          # Platt scaling 输出概率
            class_weight="balanced",
            decision_function_shape="ovr",   # One-vs-Rest 多分类
            random_state=42,
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

            scaler = StandardScaler()
            X_train_sc = scaler.fit_transform(X_train)
            X_test_sc = scaler.transform(X_test)

            clf.fit(X_train_sc, y_train)
            prob = clf.predict_proba(X_test_sc)[0]   # shape: (n_classes,)
            pred = clf.predict(X_test_sc)[0]

            y_true_list.append(int(y_test[0]))
            y_prob_list.append(prob)
            y_pred_list.append(int(pred))

        y_true = np.array(y_true_list)
        y_prob = np.array(y_prob_list)   # (n_samples, n_classes)
        y_pred = np.array(y_pred_list)

        acc = accuracy_score(y_true, y_pred)
        cm = confusion_matrix(y_true, y_pred, labels=list(range(len(classes))))

        results[clf_name] = {
            "y_true": y_true,
            "y_prob": y_prob,
            "y_pred": y_pred,
            "accuracy": acc,
            "confusion_matrix": cm,
        }
        logger.info(f"    {clf_name}: 准确率={acc:.3f}")
        logger.info("\n" + classification_report(
            y_true, y_pred,
            labels=list(range(len(classes))),
            target_names=classes,
            zero_division=0,
        ))

    # 6. 可视化
    for clf_name, res in results.items():
        _save_confusion_plot(res["confusion_matrix"], classes, clf_name, result_dir)
    _save_roc_plot(results, classes, result_dir)

    # 7. 逐样本预测结果
    rows = []
    for clf_name, res in results.items():
        for i, sid in enumerate(sample_ids):
            row = {
                "sample_id": sid,
                "true_subtype": le.inverse_transform([res["y_true"][i]])[0],
                "classifier": clf_name,
                "predicted_subtype": le.inverse_transform([res["y_pred"][i]])[0],
                "correct": bool(res["y_true"][i] == res["y_pred"][i]),
            }
            for j, cls in enumerate(classes):
                row[f"prob_{cls[:4]}"] = round(float(res["y_prob"][i, j]), 4)
            rows.append(row)
    pred_df = pd.DataFrame(rows)
    pred_path = os.path.join(result_dir, "subtype_predictions.csv")
    pred_df.to_csv(pred_path, index=False)
    logger.info(f"  逐样本预测已保存: {pred_path}")

    # 8. 汇总
    summary_rows = []
    for clf_name, res in results.items():
        summary_rows.append({
            "classifier": clf_name,
            "accuracy": round(res["accuracy"], 4),
            "n_samples": len(sample_ids),
            "n_subtypes": len(set(labels)),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(result_dir, "subtype_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"  亚型分类汇总已保存: {summary_path}")

    return results, summary_df
