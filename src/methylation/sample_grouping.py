#!/usr/bin/env python3
"""
sample_grouping.py - 样本分组模块

该模块用于处理样本分组逻辑，包括基于索引的分组和基于Excel文件的分组。
"""

import pandas as pd
import os
from utils.logger import logger
from utils.config import load_config

# 加载配置文件
config = load_config()


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
        # logger.info(f"{methylation_files}\n{benign_fq}\n{malignant_fq}")
        for i, methylation_file in enumerate(methylation_files):
            # 提取FQ文件名（不包含路径）
            fastq_basename = os.path.basename(methylation_file).replace(".txt", "") + ".fq.gz"
            # 检查是否在良性或恶性样本中
            if fastq_basename in benign_fq:
                benign_methylation_files.append(methylation_file)
            elif fastq_basename in malignant_fq:
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
            # 读取当前sheet
            df = pd.read_excel(excel_file, sheet_name=sheet_name)
            logger.info(f"处理sheet: {sheet_name}, 包含 {len(df)} 行数据")

            # 清理分类列中的空白字符
            if "分类" in df.columns:
                df["分类"] = df["分类"].str.strip()

                # 根据分类进行分组
                benign_samples = df[df["分类"] == "良性"]
                malignant_samples = df[df["分类"] == "恶性"]

                # 提取FQ文件名并添加到对应列表
                if "fq_R1" in df.columns:
                    all_benign_R1.extend(benign_samples["fq_R1"].tolist())
                    all_malignant_R1.extend(malignant_samples["fq_R1"].tolist())
                if "fq_R2" in df.columns:
                    all_benign_R2.extend(benign_samples["fq_R2"].tolist())
                    all_malignant_R2.extend(malignant_samples["fq_R2"].tolist())
            else:
                logger.warning(f"sheet {sheet_name} 中没有 '分类' 列")

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


if __name__ == "__main__":
    # 从配置文件中读取Excel文件路径
    excel_file = config.get("data", {}).get(
        "sample_info_excel", "data/按分期处理_卵巢癌fq编号与样本信息对应表.xlsx"
    )
    print(parse_sample_info(excel_file))
