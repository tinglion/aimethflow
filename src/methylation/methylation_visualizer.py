#!/usr/bin/env python3
"""
methylation_visualizer.py - 甲基化数据可视化脚本

该脚本用于可视化甲基化数据，生成甲基化水平分布图。

使用方法:
    python methylation_visualizer.py <input_file> <output_file> [--type <plot_type>]

参数:
    input_file: 甲基化数据文件路径
    output_file: 输出图片文件路径（支持.png, .pdf, .svg格式）
    --type: 图表类型，可选值：density, scatter, heatmap，默认为density
"""

import sys
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import platform

# 设置中文字体
system = platform.system()
if system == 'Windows':
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False
elif system == 'Darwin':  # macOS
    plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'PingFang SC', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
else:  # Linux
    plt.rcParams['font.sans-serif'] = ['WenQuanYi Micro Hei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False


def load_methylation_data(file_path):
    """
    加载甲基化数据

    Args:
        file_path: 甲基化数据文件路径

    Returns:
        list: 甲基化水平列表
        list: 位置列表
    """
    meth_levels = []
    positions = []

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                pos = int(parts[1])
                meth_level = float(parts[3])
                positions.append(pos)
                meth_levels.append(meth_level)

    return meth_levels, positions


def plot_density(meth_levels, output_file):
    """
    绘制甲基化水平密度图

    Args:
        meth_levels: 甲基化水平列表
        output_file: 输出文件路径
    """
    plt.figure(figsize=(10, 6))
    plt.hist(meth_levels, bins=50, density=True, alpha=0.7, color="blue")
    plt.xlabel("甲基化水平")
    plt.ylabel("密度")
    plt.title("甲基化水平分布")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"密度图已保存到: {output_file}")


def plot_scatter(positions, meth_levels, output_file):
    """
    绘制甲基化水平散点图

    Args:
        positions: 位置列表
        meth_levels: 甲基化水平列表
        output_file: 输出文件路径
    """
    plt.figure(figsize=(12, 6))
    plt.scatter(positions, meth_levels, s=5, alpha=0.5, color="green")
    plt.xlabel("基因组位置")
    plt.ylabel("甲基化水平")
    plt.title("甲基化水平分布（散点图）")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"散点图已保存到: {output_file}")


def plot_heatmap(meth_levels, output_file):
    """
    绘制甲基化水平热图

    Args:
        meth_levels: 甲基化水平列表
        output_file: 输出文件路径
    """
    # 将数据整理成矩阵形式
    n = int(np.ceil(np.sqrt(len(meth_levels))))
    # 填充数据到矩阵
    matrix = np.zeros((n, n))
    for i, val in enumerate(meth_levels):
        if i < n * n:
            matrix[i // n, i % n] = val

    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap="viridis", aspect="auto")
    plt.colorbar(label="甲基化水平")
    plt.title("甲基化水平热图")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"热图已保存到: {output_file}")


def visualize_methylation(input_file, output_file, plot_type="density"):
    """
    可视化甲基化数据

    Args:
        input_file: 甲基化数据文件路径
        output_file: 输出图片文件路径
        plot_type: 图表类型
    """
    try:
        # 检查输入文件是否存在
        if not os.path.exists(input_file):
            print(f"错误: 输入文件 '{input_file}' 不存在")
            sys.exit(1)

        # 检查输出文件格式
        ext = os.path.splitext(output_file)[1].lower()
        if ext not in [".png", ".pdf", ".svg"]:
            print("错误: 输出文件必须是.png, .pdf或.svg格式")
            sys.exit(1)

        print(f"加载甲基化数据: {input_file}")
        meth_levels, positions = load_methylation_data(input_file)

        if not meth_levels:
            print("错误: 未加载到甲基化数据")
            sys.exit(1)

        print(f"共加载 {len(meth_levels)} 个甲基化位点")

        # 根据图表类型绘制
        if plot_type == "density":
            plot_density(meth_levels, output_file)
        elif plot_type == "scatter":
            plot_scatter(positions, meth_levels, output_file)
        elif plot_type == "heatmap":
            plot_heatmap(meth_levels, output_file)
        else:
            print(f"错误: 不支持的图表类型 '{plot_type}'")
            sys.exit(1)

    except Exception as e:
        print(f"错误: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="甲基化数据可视化")
    parser.add_argument("input_file", help="甲基化数据文件路径")
    parser.add_argument("output_file", help="输出图片文件路径（支持.png, .pdf, .svg格式）")
    parser.add_argument(
        "--type",
        choices=["density", "scatter", "heatmap"],
        default="density",
        help="图表类型，默认为density",
    )

    args = parser.parse_args()

    # 执行可视化
    visualize_methylation(args.input_file, args.output_file, args.type)
