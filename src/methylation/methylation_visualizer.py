#!/usr/bin/env python3
"""
methylation_visualizer.py - Methylation data visualization script

This script is used to visualize methylation data and generate methylation level distribution plots.

Usage:
    python methylation_visualizer.py <input_file> <output_file> [--type <plot_type>]

Parameters:
    input_file: Methylation data file path
    output_file: Output image file path (supports .png, .pdf, .svg formats)
    --type: Plot type, optional values: density, scatter, heatmap, default is density
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
    Load methylation data

    Args:
        file_path: Methylation data file path

    Returns:
        list: Methylation levels list
        list: Positions list
    """
    meth_levels = []
    positions = []

    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                pos = int(parts[1])
                # 处理不同格式的文件
                if len(parts) == 4:
                    # 4列格式：chrom, start, end, methylation_percentage
                    meth_level = float(parts[3]) / 100  # 转换为0-1范围
                elif len(parts) == 6:
                    # 6列格式：chrom, start, end, methylation_percentage, meth_count, unmeth_count
                    try:
                        meth_count = int(parts[4])
                        unmeth_count = int(parts[5])
                        total_count = meth_count + unmeth_count
                        if total_count > 0:
                            meth_level = meth_count / total_count
                        else:
                            meth_level = 0.0
                    except ValueError:
                        # 如果解析失败，尝试使用第4列作为甲基化水平
                        meth_level = float(parts[3]) / 100
                else:
                    # 其他格式，使用第4列作为甲基化水平
                    meth_level = float(parts[3])
                positions.append(pos)
                meth_levels.append(meth_level)

    return meth_levels, positions


def plot_density(meth_levels, output_file):
    """
    Plot methylation level density plot

    Args:
        meth_levels: Methylation levels list
        output_file: Output file path
    """
    plt.figure(figsize=(10, 6))
    plt.hist(meth_levels, bins=50, density=True, alpha=0.7, color="blue")
    plt.xlabel("Methylation Level")
    plt.ylabel("Density")
    plt.title("Methylation Level Distribution")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Density plot saved to: {output_file}")


def plot_scatter(positions, meth_levels, output_file):
    """
    Plot methylation level scatter plot

    Args:
        positions: Positions list
        meth_levels: Methylation levels list
        output_file: Output file path
    """
    plt.figure(figsize=(12, 6))
    plt.scatter(positions, meth_levels, s=5, alpha=0.5, color="green")
    plt.xlabel("Genomic Position")
    plt.ylabel("Methylation Level")
    plt.title("Methylation Level Distribution (Scatter Plot)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Scatter plot saved to: {output_file}")


def plot_heatmap(meth_levels, output_file):
    """
    Plot methylation level heatmap

    Args:
        meth_levels: Methylation levels list
        output_file: Output file path
    """
    # Organize data into matrix form
    n = int(np.ceil(np.sqrt(len(meth_levels))))
    # Fill data into matrix
    matrix = np.zeros((n, n))
    for i, val in enumerate(meth_levels):
        if i < n * n:
            matrix[i // n, i % n] = val

    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap="viridis", aspect="auto")
    plt.colorbar(label="Methylation Level")
    plt.title("Methylation Level Heatmap")
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Heatmap saved to: {output_file}")


def visualize_methylation(input_file, output_file, plot_type="density"):
    """
    Visualize methylation data

    Args:
        input_file: Methylation data file path
        output_file: Output image file path
        plot_type: Plot type
    """
    try:
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Error: Input file '{input_file}' does not exist")
            sys.exit(1)

        # Check output file format
        ext = os.path.splitext(output_file)[1].lower()
        if ext not in [".png", ".pdf", ".svg"]:
            print("Error: Output file must be in .png, .pdf, or .svg format")
            sys.exit(1)

        print(f"Loading methylation data: {input_file}")
        meth_levels, positions = load_methylation_data(input_file)

        if not meth_levels:
            print("Error: No methylation data loaded")
            sys.exit(1)

        print(f"Loaded {len(meth_levels)} methylation sites")

        # Plot according to plot type
        if plot_type == "density":
            plot_density(meth_levels, output_file)
        elif plot_type == "scatter":
            plot_scatter(positions, meth_levels, output_file)
        elif plot_type == "heatmap":
            plot_heatmap(meth_levels, output_file)
        else:
            print(f"Error: Unsupported plot type '{plot_type}'")
            sys.exit(1)

    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Methylation data visualization")
    parser.add_argument("input_file", help="Methylation data file path")
    parser.add_argument("output_file", help="Output image file path (supports .png, .pdf, .svg formats)")
    parser.add_argument(
        "--type",
        choices=["density", "scatter", "heatmap"],
        default="density",
        help="Plot type, default is density",
    )

    args = parser.parse_args()

    # Execute visualization
    visualize_methylation(args.input_file, args.output_file, args.type)
