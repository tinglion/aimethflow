#!/usr/bin/env python3
"""
methylation_utils.py - Methylation processing module

Functions for methylation extraction and processing.
"""

import os
import logging

logger = logging.getLogger(__name__)


def extract_methylation_info(temp_dir, fastq_files, group_method="v2", config=None, bam_files=None, category="plasma"):
    """
    提取甲基化信息并按样本分组
    
    Args:
        temp_dir (str): 工作目录，用于存储临时文件和输出结果
        fastq_files (list): FastQ 文件列表（当未提供 bam_files 时用于生成 BAM 文件）
        group_method (str, optional): 分组方法，可选 "v1"（基于索引）或 "v2"（基于 Excel）。默认为 "v2"
        config (dict, optional): 配置字典，包含参考基因组和线程数等参数。默认为 None
        bam_files (list, optional): 预计算的 BAM 文件列表（血浆模式：由比对→去重→靶向捕获步骤生成）。
                                     当提供时，跳过内部 BAM 生成（避免对 BS-seq 数据运行普通 bwa mem）。默认为 None
        category (str, optional): 分组类别，例如 "plasma" 或 "tissue"。默认为 "plasma"
    
    Returns:
        tuple: 包含三个元素的元组 (sample1_methylation, sample2_methylation, borderline_methylation)，
               分别为样本组1、样本组2和边界样本的甲基化文件路径
    """
    """
    Extract methylation information.

    Args:
        temp_dir:     working directory
        fastq_files:  FastQ file list (used to generate BAM when bam_files is None)
        group_method: grouping method — "v1" (index-based) or "v2" (Excel-based)
        config:       config dict
        bam_files:    pre-computed BAM list (plasma mode: produced by align → dedup →
                      on-target steps). When provided, internal BAM generation is skipped
                      (avoids running plain bwa mem on BS-seq data).
        category:     category for grouping (e.g., "plasma" or "tissue")

    Returns:
        tuple: (sample1_methylation, sample2_methylation) file paths
    """
    logger.info("Starting methylation extraction")

    reference_genome = config.get("extraction", {}).get("reference_genome")
    threads = config.get("extraction", {}).get("threads", 4)

    all_bam_files = []
    if bam_files:
        # use externally supplied BAMs (plasma mode: already bismark-aligned, UMI-deduped, on-target)
        logger.info(f"Using {len(bam_files)} pre-computed BAM file(s) — skipping internal alignment")
        all_bam_files = bam_files
    elif fastq_files:
        # fallback: generate BAM from FastQ (only for non-BS-seq data or testing)
        logger.warning("No pre-computed BAMs provided — generating BAM with bwa mem (NOTE: not suitable for BS-seq data!)")
        from methylation.bam_generator import generate_bam
        all_bam_files = generate_bam(fastq_files, reference_genome, threads=threads)

    all_methylation_files = []
    if all_bam_files:
        logger.info(f"Extracting methylation from {len(all_bam_files)} BAM file(s)")
        from methylation.methylation_extractor import extract_methylation
        for i, bam_file in enumerate(all_bam_files):
            individual_output = bam_file.replace(".bam", ".txt")
            logger.info(f"extracting methylation: {bam_file}")
            extract_methylation(bam_file, individual_output, reference_genome=reference_genome)
            validate_methylation_output(individual_output)
            all_methylation_files.append(individual_output)

    sample1_methylation = os.path.join(temp_dir, "sample1_methylation.txt")
    sample2_methylation = os.path.join(temp_dir, "sample2_methylation.txt")
    borderline_methylation = os.path.join(temp_dir, "borderline_methylation.txt")

    if group_method == "v2":
        logger.info("using Excel-based grouping (v2)")
        from methylation.sample_grouping import group_methylation_files_v2
        sample1_methylation_files, sample2_methylation_files, borderline_methylation_files = group_methylation_files_v2(
            all_methylation_files, category=category
        )
    else:
        logger.info("using index-based grouping (v1)")
        from methylation.sample_grouping import group_methylation_files_v1
        sample1_methylation_files, sample2_methylation_files = group_methylation_files_v1(
            all_methylation_files
        )
        borderline_methylation_files = []
    logger.info(f"Grouping done\ngroup1: {sample1_methylation_files}\ngroup2: {sample2_methylation_files}")

    if sample1_methylation_files:
        logger.info(f"merging {len(sample1_methylation_files)} methylation file(s) for group 1")
        merge_methylation_files(sample1_methylation_files, sample1_methylation)
    logger.info(f"group 1 methylation ready: {sample1_methylation}")

    if sample2_methylation_files:
        logger.info(f"merging {len(sample2_methylation_files)} methylation file(s) for group 2")
        merge_methylation_files(sample2_methylation_files, sample2_methylation)
    logger.info(f"group 2 methylation ready: {sample2_methylation}")

    if borderline_methylation_files:
        logger.info(f"merging {len(borderline_methylation_files)} methylation file(s) for borderline samples")
        merge_methylation_files(borderline_methylation_files, borderline_methylation)
        logger.info(f"borderline samples methylation ready: {borderline_methylation}")

    return sample1_methylation, sample2_methylation, borderline_methylation


def validate_methylation_output(output_file):
    """
    Validate methylation extraction output.

    Args:
        output_file: path to the methylation extraction output file
    """
    logger.info(f"Validating methylation output: {output_file}")

    try:
        if not os.path.exists(output_file):
            logger.error(f"output file not found: {output_file}")
            return False

        file_size = os.path.getsize(output_file)
        if file_size == 0:
            logger.error(f"output file is empty: {output_file}")
            return False

        with open(output_file, "r", encoding="utf-8") as f:
            header = f.readline()
            if not header.startswith("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数"):
                logger.warning(f"unexpected header in {output_file}")

            line_count = 0
            for line in f:
                if not line.startswith("#"):
                    parts = line.strip().split()
                    if len(parts) != 6:
                        logger.warning(f"malformed line {line_count+2} in {output_file}")
                    line_count += 1

        logger.info(f"Validation done: {output_file} has {line_count} methylation sites")
        return True
    except Exception as e:
        logger.error(f"validation error: {str(e)}")
        return False


def generate_mock_methylation_data(output_file, min_level, max_level):
    """
    Generate mock methylation data for testing.

    Args:
        output_file: output file path
        min_level:   minimum methylation level
        max_level:   maximum methylation level
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")
        num_sites = 1000
        range_level = max_level - min_level
        for i in range(1, num_sites + 1):
            chrom = "chr1"
            pos = 100000 + i * 100
            strand = "+" if i % 2 == 0 else "-"
            meth_level = round(min_level + (i % 50) / 50 * range_level, 2)
            meth_count = int(meth_level * 20)
            unmeth_count = 20 - meth_count
            f.write(f"{chrom}\t{pos}\t{strand}\t{meth_level}\t{meth_count}\t{unmeth_count}\n")


def merge_methylation_files(input_files, output_file):
    """
    Merge multiple methylation files.

    Args:
        input_files: list of input methylation file paths
        output_file: output file path
    """
    logger.info(f"Merging {len(input_files)} methylation file(s) → {output_file}")

    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write("# 染色体\t位置\t链\t甲基化水平\t甲基化计数\t非甲基化计数\n")

        for input_file in input_files:
            logger.info(f"  merging: {input_file}")
            with open(input_file, "r", encoding="utf-8") as infile:
                # skip header and track lines
                header_skipped = False
                for line in infile:
                    if line.startswith("#") or line.startswith("track"):
                        if not header_skipped and line.startswith("#"):
                            header_skipped = True
                        continue
                    outfile.write(line)

    logger.info(f"Merge done: {output_file}")


def filter_methylation_data(methylation_file, temp_dir, config):
    """
    Filter methylation data.

    Args:
        methylation_file: methylation file path
        temp_dir:         working directory
        config:           config dict

    Returns:
        str: path to the filtered methylation file
    """
    logger.info(f"Filtering methylation data: {methylation_file}")

    # load filtering config
    repeat_bed = config.get("filtering", {}).get("repeat_bed")
    promoter_bed = config.get("filtering", {}).get("promoter_bed")
    depth_min = config.get("filtering", {}).get("min_depth", 5)
    depth_max = config.get("filtering", {}).get("max_depth", 50)

    try:
        base_name = os.path.basename(methylation_file)
        filtered_file = os.path.join(temp_dir, f"filtered_{base_name}")
        if os.path.exists(methylation_file):
            from methylation.filter_methylation import filter_methylation
            filter_methylation(
                input_file=methylation_file,
                output_file=filtered_file,
                output_dir=temp_dir,
                repeat_bed=repeat_bed,
                promoter_bed=promoter_bed,
                depth_min=depth_min,
                depth_max=depth_max,
            )
            logger.info(f"Filtering done: {filtered_file}")
            return filtered_file
        else:
            logger.warning(f"input file not found, skipping filter: {methylation_file}")
            return methylation_file
    except Exception as e:
        error_msg = f"error filtering methylation data: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def filter_to_candidate_regions_wrapper(methylation_file, temp_dir, config):
    """
    Filter methylation data to candidate regions.

    Args:
        methylation_file: methylation file path
        temp_dir:         working directory
        config:           config dict

    Returns:
        str: path to the candidate-region-filtered methylation file
    """
    logger.info(f"Filtering to candidate regions: {methylation_file}")

    candidate_regions = config.get("data", {}).get("candidate_regions")

    try:
        base_name = os.path.basename(methylation_file)
        candidate_output = os.path.join(temp_dir, f"candidate_{base_name}")

        if os.path.exists(methylation_file) and os.path.exists(candidate_regions):
            from methylation.candidate_regions import filter_to_candidate_regions
            if filter_to_candidate_regions(
                input_file=methylation_file,
                output_file=candidate_output,
                candidate_regions=candidate_regions,
                output_dir=temp_dir,
            ):
                logger.info(f"Candidate-region filter done: {candidate_output}")
                return candidate_output
            else:
                logger.warning("Candidate-region filter failed, using original file")
                return methylation_file
        else:
            if not os.path.exists(candidate_regions):
                logger.warning(f"Candidate regions file not found, skipping: {candidate_regions}")
            return methylation_file
    except Exception as e:
        error_msg = f"error filtering to candidate regions: {str(e)}"
        logger.error(error_msg)
        return methylation_file


def perform_diff_analysis(sample1_methylation, sample2_methylation, result_dir):
    """
    Run differential methylation analysis.

    Args:
        sample1_methylation: group 1 methylation file path
        sample2_methylation: group 2 methylation file path
        result_dir:          results directory
    """
    logger.info("Starting differential methylation analysis")
    diff_result = os.path.join(result_dir, "diff_methylation.txt")

    try:
        from methylation.diff_methylation import run_diff_analysis
        run_diff_analysis(sample1_methylation, sample2_methylation, diff_result, threshold=0.2)
        logger.info("Differential methylation analysis done")
    except Exception as e:
        error_msg = f"differential methylation analysis error: {str(e)}"
        logger.error(error_msg)


def plot_sample_clustered_heatmap(result_dir, output_file):
    """
    Build a samples × CpG-sites heatmap, grouped by benign/malignant class,
    with hierarchical clustering applied within each class so similar samples
    are placed adjacent. Reads `feature_matrix.npz` from result_dir (written
    by Step 14 feature_analysis); falls back to the old single-sample reshape
    heatmap when the matrix is unavailable.

    Returns
    -------
    True  — sample-clustered heatmap was generated.
    False — feature_matrix.npz missing; caller should use the fallback path.
    """
    import os as _os
    npz_path = _os.path.join(result_dir, "feature_matrix.npz")
    if not _os.path.exists(npz_path):
        return False

    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist
    except ImportError as exc:
        logger.warning(f"  Heatmap dependencies missing ({exc}); falling back.")
        return False

    data = np.load(npz_path, allow_pickle=True)
    X            = data["X"].astype(float)
    y            = np.asarray(data["y"])
    sample_ids   = list(data["sample_ids"])
    feature_names = list(data["feature_names"])

    if X.size == 0 or len(y) == 0:
        logger.warning("  feature_matrix.npz is empty; skipping clustered heatmap.")
        return False

    # Restrict to the panel CpG sites when available — smaller, cleaner heatmap
    # focused on the diagnostically informative subset. Without the panel, fall
    # back to all features (still valid, just larger).
    panel_path = _os.path.join(result_dir, "panel_cpg.csv")
    if _os.path.exists(panel_path):
        try:
            import pandas as pd
            panel_set = set(pd.read_csv(panel_path)["feature"])
            keep_cols = np.array([fn in panel_set for fn in feature_names])
            if keep_cols.any():
                X = X[:, keep_cols]
                feature_names = [f for f, k in zip(feature_names, keep_cols) if k]
                logger.info(f"  Heatmap restricted to {len(feature_names)} panel CpG sites")
        except Exception as exc:
            logger.debug(f"  Panel restriction failed ({exc}); using full feature matrix.")

    # Median-impute NaNs column-wise so clustering metrics are well-defined.
    col_medians = np.nanmedian(X, axis=0)
    # Columns that are entirely NaN → impute with 0.5 (uninformative midpoint).
    col_medians = np.where(np.isnan(col_medians), 0.5, col_medians)
    nan_mask = np.isnan(X)
    if nan_mask.any():
        X = X.copy()
        X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

    def _cluster_rows(rows):
        """Return a row permutation that places similar rows adjacent.
        Uses correlation distance + average linkage; gracefully returns the
        identity permutation when the data is too small or pdist fails."""
        if len(rows) < 2 or rows.shape[1] < 2:
            return np.arange(len(rows))
        try:
            D = pdist(rows, metric="correlation")
            # pdist can produce NaN distances when rows have zero variance;
            # replace with the max finite distance so linkage still works.
            if np.any(np.isnan(D)):
                D = np.nan_to_num(D, nan=np.nanmax(D[~np.isnan(D)]) if np.any(~np.isnan(D)) else 1.0)
            Z = linkage(D, method="average")
            return np.array(leaves_list(Z))
        except Exception as exc:
            logger.debug(f"    Row clustering failed ({exc}); using natural order.")
            return np.arange(len(rows))

    benign_idx = np.where(y == 0)[0]
    malig_idx  = np.where(y == 1)[0]
    order_b = benign_idx[_cluster_rows(X[benign_idx])] if len(benign_idx) else np.array([], dtype=int)
    order_m = malig_idx [_cluster_rows(X[malig_idx ])] if len(malig_idx ) else np.array([], dtype=int)
    row_order = np.concatenate([order_b, order_m])

    # Cluster columns too (similar CpG sites adjacent) — purely cosmetic but
    # makes block-structured methylation patterns visible to reviewers.
    if X.shape[1] >= 2:
        try:
            D_c = pdist(X.T, metric="correlation")
            if np.any(np.isnan(D_c)):
                D_c = np.nan_to_num(D_c, nan=np.nanmax(D_c[~np.isnan(D_c)]) if np.any(~np.isnan(D_c)) else 1.0)
            Z_c = linkage(D_c, method="average")
            col_order = np.array(leaves_list(Z_c))
        except Exception:
            col_order = np.arange(X.shape[1])
    else:
        col_order = np.arange(X.shape[1])

    X_plot = X[row_order][:, col_order]
    n_samp, n_feat = X_plot.shape

    # Figure geometry — wide enough that adjacent CpGs are distinguishable,
    # tall enough that ~120 sample rows are individually readable.
    fig_w = max(9.0, min(20.0, n_feat / 35.0))
    fig_h = max(4.5, min(14.0, n_samp / 9.0 + 1.5))

    fig, (ax_anno, ax_heat) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={"width_ratios": [1, 60], "wspace": 0.015},
    )

    # ── Class annotation strip (left of the main heatmap) ──
    cls_strip = np.array([[0]] * len(order_b) + [[1]] * len(order_m), dtype=float)
    anno_cmap = LinearSegmentedColormap.from_list("anno", ["#4878CF", "#D65F5F"])
    ax_anno.imshow(cls_strip, aspect="auto", cmap=anno_cmap, vmin=0, vmax=1, interpolation="nearest")
    ax_anno.set_xticks([])
    ax_anno.set_yticks([])
    ax_anno.set_title("Class", fontsize=9, pad=4)
    # Boundary line between classes
    if len(order_b) and len(order_m):
        ax_anno.axhline(len(order_b) - 0.5, color="black", lw=0.9)

    # ── Main heatmap ──
    im = ax_heat.imshow(
        X_plot, aspect="auto", cmap="RdYlBu_r",
        vmin=0, vmax=1, interpolation="nearest",
    )
    ax_heat.set_yticks([])
    ax_heat.set_xticks([])
    if len(order_b) and len(order_m):
        ax_heat.axhline(len(order_b) - 0.5, color="black", lw=1.2)

    ax_heat.set_xlabel(
        f"{n_feat} CpG sites (correlation-clustered; similar adjacent)",
        fontsize=10,
    )
    ax_heat.set_title(
        f"Methylation Heatmap — {len(order_b)} benign (blue strip) / "
        f"{len(order_m)} malignant (red strip); "
        f"rows hierarchically clustered within each class",
        fontsize=11, pad=8,
    )

    # Colorbar
    cbar = fig.colorbar(im, ax=ax_heat, shrink=0.55, pad=0.015, fraction=0.025)
    cbar.set_label("Methylation β (0–1)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    # Legend for the class strip
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="#4878CF", label=f"Benign (n={len(order_b)})"),
        Patch(facecolor="#D65F5F", label=f"Malignant (n={len(order_m)})"),
    ]
    ax_heat.legend(
        handles=legend_handles, loc="upper right",
        bbox_to_anchor=(1.0, -0.04), ncol=2, frameon=False, fontsize=9,
    )

    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  Sample-clustered heatmap saved: {output_file} "
                f"({n_samp} samples × {n_feat} CpG sites)")
    return True


def visualize_results(sample1_methylation, result_dir):
    """
    Visualize analysis results.

    Args:
        sample1_methylation: group 1 methylation file path
        result_dir:          results directory
    """
    logger.info("Starting result visualization")
    density_plot = os.path.join(result_dir, "density_plot.png")
    scatter_plot = os.path.join(result_dir, "scatter_plot.png")
    heatmap_plot = os.path.join(result_dir, "heatmap_plot.png")

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, density_plot, plot_type="density")
        logger.info("Density plot saved")
    except Exception as e:
        error_msg = f"error generating density plot: {str(e)}"
        logger.error(error_msg)

    try:
        from methylation.methylation_visualizer import visualize_methylation
        visualize_methylation(sample1_methylation, scatter_plot, plot_type="scatter")
        logger.info("Scatter plot saved")
    except Exception as e:
        error_msg = f"error generating scatter plot: {str(e)}"
        logger.error(error_msg)

    # Prefer the sample-clustered heatmap (samples × CpG sites, grouped by class)
    # when Step 14's feature_matrix.npz is available. Fall back to the legacy
    # single-sample reshape only when the matrix doesn't exist yet — e.g. when
    # Step 12 is run standalone before Step 14.
    try:
        if plot_sample_clustered_heatmap(result_dir, heatmap_plot):
            pass  # already logged inside
        else:
            from methylation.methylation_visualizer import visualize_methylation
            visualize_methylation(sample1_methylation, heatmap_plot, plot_type="heatmap")
            logger.info("Heatmap saved (legacy single-sample reshape — "
                        "run --step feature_analysis first for sample-clustered version)")
    except Exception as e:
        error_msg = f"error generating heatmap: {str(e)}"
        logger.error(error_msg)
