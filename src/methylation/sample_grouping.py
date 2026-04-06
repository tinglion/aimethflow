#!/usr/bin/env python3
"""
sample_grouping.py - 样本分组模块

该模块用于处理样本分组逻辑，包括基于索引的分组和基于Excel文件的分组。
"""

import pandas as pd
import os
import re
from utils.logger import logger
from utils.config import load_config

# 加载配置文件
config = load_config()


def _to_posix_path(path):
    """
    将 Windows 路径（如 F:\\meth\\OC_plasma）转换为 WSL 可用的 POSIX 路径
    （/mnt/f/meth/OC_plasma）。非 Windows 路径原样返回。
    """
    if not path:
        return path
    # 匹配 X:\ 或 X:/ 开头的 Windows 绝对路径
    m = re.match(r"^([A-Za-z]):[/\\](.*)", path)
    if m:
        drive = m.group(1).lower()
        rest  = m.group(2).replace("\\", "/")
        return f"/mnt/{drive}/{rest}"
    # 已是 POSIX 路径，只做反斜杠替换
    return path.replace("\\", "/")


def group_methylation_files_v1(methylation_files):
    """
    将甲基化文件分为两组（样本1和样本2），使用mod 2均匀分组（v1版本）

    Args:
        methylation_files: 甲基化文件列表

    Returns:
        tuple: (sample1_methylation_files, sample2_methylation_files)
    """
    sample1_methylation_files = [f for i, f in enumerate(methylation_files) if i % 2 == 0]
    sample2_methylation_files = [f for i, f in enumerate(methylation_files) if i % 2 == 1]
    return sample1_methylation_files, sample2_methylation_files


def group_methylation_files_v2(methylation_files, fastq_files):
    """
    将甲基化文件分为两组（样本1和样本2），根据Excel文件中的良/恶性进行分组（v2版本）

    Args:
        methylation_files: 甲基化文件列表
        fastq_files: FastQ文件列表

    Returns:
        tuple: (benign_methylation_files, malignant_methylation_files)
    """
    # 从配置文件中读取Excel文件路径
    excel_file = config.get("data", {}).get(
        "sample_info_excel", "data/按分期处理_卵巢癌fq编号与样本信息对应表.xlsx"
    )
    if not os.path.exists(excel_file):
        logger.error(f"错误: 文件 {excel_file} 不存在")
        # 如果Excel文件不存在，回退到v1版本
        return group_methylation_files_v1(methylation_files)

    try:
        # 使用parse_sample_info函数解析Excel文件
        sample_info = parse_sample_info(excel_file)
        if not sample_info:
            logger.error("解析Excel文件失败，回退到v1版本")
            return group_methylation_files_v1(methylation_files)

        # 提取良性和恶性样本的FQ文件名
        benign_fq = set(sample_info["benign"]["R1"])
        malignant_fq = set(sample_info["malignant"]["R1"])
        logger.info(f"根据Excel解析到 {len(benign_fq)} 个良性样本, {len(malignant_fq)} 个恶性样本")

        # 根据FastQ文件名匹配甲基化文件
        benign_methylation_files = []
        malignant_methylation_files = []

        logger.info(f"开始根据Excel文件分组 {len(methylation_files)} 个甲基化文件 ")
        for i, methylation_file in enumerate(methylation_files):
            # 提取样本ID（去掉路径、后缀及 _CpG / _CHG 等 MethylDackel 后缀）
            sample_id = os.path.basename(methylation_file)
            for ext in [".txt", ".bedGraph", ".cov", ".gz"]:
                sample_id = sample_id.replace(ext, "")
            # 先按 . 取第一段，再按 _ 取第一段，得到纯样本ID如 B10002
            sample_id = sample_id.split(".")[0].split("_")[0]

            def _in_set(sample_id, fq_set):
                """检查样本ID是否对应fq_set中任何一个文件（兼容 B10002.R1.fq.gz 等命名）。"""
                for r1 in fq_set:
                    r1_base = os.path.splitext(os.path.splitext(r1)[0])[0]  # 去除双后缀
                    r1_sample = r1_base.split(".")[0].split("_")[0]
                    if r1_sample == sample_id:
                        return True
                return False

            if _in_set(sample_id, benign_fq):
                benign_methylation_files.append(methylation_file)
            elif _in_set(sample_id, malignant_fq):
                malignant_methylation_files.append(methylation_file)

        logger.info(
            f"根据Excel文件分组结果: 良性 {len(benign_methylation_files)} 个文件, 恶性 {len(malignant_methylation_files)} 个文件"
        )
        return benign_methylation_files, malignant_methylation_files

    except Exception as e:
        logger.error(f"根据Excel文件分组时出错: {str(e)}")
        # 如果出错，回退到v1版本
        return group_methylation_files_v1(methylation_files)


def parse_sample_info(excel_file=None):
    """
    解析Excel文件中的样本信息，支持多个sheet

    Args:
        excel_file: Excel文件路径，如果为None则从配置文件中读取

    Returns:
        dict: 包含良性和恶性样本分组的字典
    """
    try:
        # 如果未提供Excel文件路径，从配置文件中读取
        if excel_file is None:
            excel_file = config.get("data", {}).get(
                "sample_info_excel", "data/按分期处理_卵巢癌fq编号与样本信息对应表.xlsx"
            )

        # 获取Excel文件的所有sheet名称
        xl = pd.ExcelFile(excel_file)
        sheet_names = xl.sheet_names
        logger.info(f"Excel文件包含以下sheet: {sheet_names}")

        # 存储所有sheet的数据
        all_benign_R1 = []
        all_benign_R2 = []
        all_malignant_R1 = []
        all_malignant_R2 = []

        # 遍历所有sheet
        for sheet_name in sheet_names:
            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            logger.info(f"处理sheet: {sheet_name}, 包含 {len(df)} 行数据")

            if "分类" not in df.columns:
                logger.warning(f"sheet {sheet_name} 中没有 '分类' 列，跳过")
                continue

            # 清理字符串列中的空白字符
            df["分类"] = df["分类"].astype(str).str.strip()
            if "亚型" in df.columns:
                df["亚型"] = df["亚型"].astype(str).str.strip()

            # 明确排除"交界性"（仅保留"良性"和"恶性"用于分类器）
            borderline = df[df["分类"] == "交界性"]
            if not borderline.empty:
                logger.info(f"  排除交界性样本 {len(borderline)} 个（不参与分类器训练）")

            benign_samples    = df[df["分类"] == "良性"]
            malignant_samples = df[df["分类"] == "恶性"]

            logger.info(f"  良性: {len(benign_samples)}  恶性: {len(malignant_samples)}  "
                        f"交界性: {len(borderline)}")

            if "fq_R1" in df.columns:
                all_benign_R1.extend(benign_samples["fq_R1"].dropna().tolist())
                all_malignant_R1.extend(malignant_samples["fq_R1"].dropna().tolist())
            if "fq_R2" in df.columns:
                all_benign_R2.extend(benign_samples["fq_R2"].dropna().tolist())
                all_malignant_R2.extend(malignant_samples["fq_R2"].dropna().tolist())

        # 去重，避免重复的样本
        all_benign_R1 = list(set(all_benign_R1))
        all_benign_R2 = list(set(all_benign_R2))
        all_malignant_R1 = list(set(all_malignant_R1))
        all_malignant_R2 = list(set(all_malignant_R2))

        # 构建返回结果
        benign_fq = {"R1": all_benign_R1, "R2": all_benign_R2}
        malignant_fq = {"R1": all_malignant_R1, "R2": all_malignant_R2}

        logger.info(
            f"解析完成，共获取 {len(all_benign_R1)} 个良性样本, {len(all_malignant_R1)} 个恶性样本"
        )
        return {"benign": benign_fq, "malignant": malignant_fq}

    except Exception as e:
        logger.error(f"解析Excel文件时出错: {str(e)}")
        return None


def resolve_patient_files(csv_path, excel_file=None, data_dir=None):
    """
    从患者ID CSV 解析对应的 R1/R2 FastQ 文件绝对路径。

    Args:
        csv_path:   患者ID列表文件（每行一个ID，无表头），如 data/patient3.csv
        excel_file: Excel文件路径（含 fq_R1/fq_R2 列），None 时从配置读取
        data_dir:   FastQ 文件所在目录，None 时从配置读取（plasma_dir）

    Returns:
        list: FastQ 文件路径列表（R1、R2 交替排列）
    """
    # 读取患者 ID 列表（去空行/空白）
    with open(csv_path, encoding="utf-8") as f:
        patient_ids = [line.strip() for line in f if line.strip()]
    logger.info(f"从 {csv_path} 读取到 {len(patient_ids)} 个患者ID")

    # Excel 文件路径
    if excel_file is None:
        excel_file = config.get("data", {}).get(
            "sample_info_excel", "data/按分期处理_卵巢癌fq编号与样本信息对应表.xlsx"
        )

    # FastQ 数据目录（自动将 Windows 路径转换为 WSL POSIX 路径）
    if data_dir is None:
        data_dir = config.get("data", {}).get("plasma_dir", "")
    data_dir = _to_posix_path(data_dir)
    # 从 Excel 所有 sheet 构建 patient_id → {R1, R2} 映射
    # patient_id = fq_R1 文件名第一个 "." 之前的部分（如 B10186.R1.fq.gz → B10186）
    xl = pd.ExcelFile(excel_file)
    pid_to_files = {}
    for sheet_name in xl.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        if "fq_R1" not in df.columns:
            continue
        if "分类" in df.columns:
            df["分类"] = df["分类"].astype(str).str.strip()
        if "亚型" in df.columns:
            df["亚型"] = df["亚型"].astype(str).str.strip()
        for _, row in df.iterrows():
            r1 = str(row["fq_R1"]).strip() if pd.notna(row.get("fq_R1")) else ""
            r2 = str(row.get("fq_R2", "")).strip() if pd.notna(row.get("fq_R2")) else ""
            if not r1 or r1 == "nan":
                continue
            pid = r1.split(".")[0]
            category = str(row.get("分类", "")).strip() if "分类" in df.columns else ""
            subtype  = str(row.get("亚型",  "")).strip() if "亚型"  in df.columns else ""
            if pid not in pid_to_files:
                pid_to_files[pid] = {"R1": r1, "R2": r2,
                                     "分类": category, "亚型": subtype}

    logger.info(f"从Excel解析到 {len(pid_to_files)} 个样本映射")

    # 按患者ID列表顺序解析文件路径
    fastq_files = []
    missing = []
    skipped_borderline = []
    for pid in patient_ids:
        if pid not in pid_to_files:
            logger.warning(f"  患者ID {pid} 在Excel中未找到")
            missing.append(pid)
            continue
        info = pid_to_files[pid]
        # 跳过交界性样本（不参与比对及后续分析）
        if info.get("分类") == "交界性":
            skipped_borderline.append(pid)
            continue
        r1_name = info["R1"]
        r2_name = info["R2"]
        r1_path = os.path.join(data_dir, r1_name)
        r2_path = os.path.join(data_dir, r2_name) if r2_name else ""
        if not os.path.exists(r1_path):
            logger.warning(f"  R1 文件不存在: {r1_path}")
            missing.append(pid)
            continue
        if r2_name and not os.path.exists(r2_path):
            logger.warning(f"  R2 文件不存在: {r2_path}")
            missing.append(pid)
            continue
        fastq_files.append(r1_path)
        if r2_name:
            fastq_files.append(r2_path)

    found = len(patient_ids) - len(missing) - len(skipped_borderline)
    logger.info(f"解析完成: {found}/{len(patient_ids)} 个患者ID找到对应文件")
    if skipped_borderline:
        logger.info(f"  跳过交界性样本 ({len(skipped_borderline)} 个): {skipped_borderline}")
    if missing:
        logger.warning(f"  缺失/未找到的患者ID: {missing}")

    return fastq_files


if __name__ == "__main__":
    # 从配置文件中读取Excel文件路径
    excel_file = config.get("data", {}).get(
        "sample_info_excel", "data/按分期处理_卵巢癌fq编号与样本信息对应表.xlsx"
    )
    print(parse_sample_info(excel_file))