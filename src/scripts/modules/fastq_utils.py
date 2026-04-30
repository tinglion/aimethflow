#!/usr/bin/env python3
"""
fastq_utils.py - FastQ file processing module

Quality analysis and preprocessing functions for FastQ files.
"""

import os
import subprocess
import random
import gzip
from Bio import SeqIO
from utils.logger import logger


def analyze_quality(fastq_files):
    """
    Analyze FastQ file quality.

    Args:
        fastq_files: list of FastQ file paths
    """
    logger.info("Starting FastQ quality analysis")
    from fq_inspect import inspect_fastq

    for fastq_file in fastq_files:
        logger.info(f"analyzing: {fastq_file}")
        try:
            inspect_fastq(fastq_file)
            logger.info(f"done: {fastq_file}")
        except Exception as e:
            error_msg = f"error analyzing {fastq_file}: {str(e)}"
            logger.error(error_msg)


def get_fastq_stats(fastq_file, sample_size=None, random_seed=None, config=None):
    """
    Get basic statistics from a FastQ file (via random sampling).

    Args:
        fastq_file:  path to the FastQ file
        sample_size: number of records to sample (default: from config or 1000)
        random_seed: random seed (default: from config or 42)
        config:      config dict

    Returns:
        dict with keys: total_records, total_bases, avg_length, avg_quality
    """
    # use config values or defaults
    if sample_size is None and config:
        sample_size = config.get("fastq", {}).get("sample_size", 1000)
    else:
        sample_size = 1000

    if random_seed is None and config:
        random_seed = config.get("fastq", {}).get("random_seed", 42)
    else:
        random_seed = 42

    stats = {"total_records": 0, "total_bases": 0, "avg_length": 0, "avg_quality": 0}

    try:
        random.seed(random_seed)

        # quickly estimate total record count by reading a prefix of the file
        total_sequences = 0
        try:
            file_size = os.path.getsize(fastq_file)
            with gzip.open(fastq_file, "rt") as handle:
                sample_data = handle.read(100000)
                sample_sequences = sample_data.count("\n@")
                if sample_sequences > 0:
                    total_sequences = int(
                        (file_size / len(sample_data.encode("utf-8"))) * sample_sequences
                    )
                else:
                    total_sequences = 1000000
        except Exception as e:
            logger.warning(f"could not estimate record count, using default: {str(e)}")
            total_sequences = 1000000

        # compute sampling probability
        if total_sequences <= sample_size:
            sample_prob = 1.0
        else:
            sample_prob = sample_size / total_sequences

        logger.info(f"sample {fastq_file} p={sample_prob:.8f}={sample_size}/{total_sequences}")
        with gzip.open(fastq_file, "rt") as handle:
            lengths = []
            quality_scores = []

            for record in SeqIO.parse(handle, "fastq"):
                if random.random() <= sample_prob:
                    stats["total_records"] += 1
                    length = len(record)
                    stats["total_bases"] += length
                    lengths.append(length)

                    if record.letter_annotations.get("phred_quality"):
                        avg_quality = sum(record.letter_annotations["phred_quality"]) / length
                        quality_scores.append(avg_quality)

                if stats["total_records"] >= sample_size:
                    break

            if lengths:
                stats["avg_length"] = stats["total_bases"] / stats["total_records"]

            if quality_scores:
                stats["avg_quality"] = sum(quality_scores) / len(quality_scores)

        stats["original_total_records"] = total_sequences

    except Exception as e:
        logger.error(f"error getting FastQ stats: {str(e)}")

    return stats


def find_r2_for_r1(r1_file):
    """
    Find the R2 file matching an R1 path, compatible with .R1. and _R1 naming styles.

    Returns:
        str or None: R2 path, or None if not found
    """
    for r1_pat, r2_pat in [(".R1.", ".R2."), ("_R1.", "_R2."), ("_R1_", "_R2_"),
                            (".r1.", ".r2."), ("_r1.", "_r2.")]:
        if r1_pat in r1_file:
            r2_file = r1_file.replace(r1_pat, r2_pat)
            if os.path.exists(r2_file):
                return r2_file
    return None


def is_r1_file(fastq_file):
    """Return True if this file is an R1 read (to avoid processing R2 files twice)."""
    for pat in [".R1.", "_R1.", "_R1_", ".r1.", "_r1."]:
        if pat in fastq_file:
            return True
    return False


def preprocess_fastq(fastq_files, temp_dir, plasma=False, config=None):
    """
    Data preprocessing: quality control and adapter trimming.
    Uses fastp instead of TrimGalore; supports UMI extraction in plasma mode.
    Idempotent: skips a sample if the .clean.R1.fq.gz output already exists.

    Args:
        fastq_files: list of FastQ file paths
        temp_dir:    working directory for intermediate files
        plasma:      True to enable plasma mode (UMI extraction + 5' low-complexity
                     trimming + 3' poly-C trimming)
        config:      config dict

    Returns:
        list: preprocessed FastQ paths (R1 then R2 for each pair)
    """
    logger.info(f"Starting preprocessing ({'plasma' if plasma else 'standard'} mode)")
    processed_files = []

    # plasma config
    plasma_cfg = (config or {}).get("plasma", {})
    trim_front1 = plasma_cfg.get("trim_front1", 12)
    trim_front2 = plasma_cfg.get("trim_front2", 12)
    trim_tail2  = plasma_cfg.get("trim_tail2", 5)
    min_length  = plasma_cfg.get("min_length", 35)
    min_quality = plasma_cfg.get("min_quality", 20)
    umi_length  = plasma_cfg.get("umi_length", 6)

    # thread count from config [threads]
    n_threads = str((config or {}).get("threads", {}).get("fastp", 8))

    # find all R1 files (to avoid processing R2 files twice)
    r1_files = [f for f in fastq_files if is_r1_file(f)]
    # collect unpaired files (single-end)
    unpaired_files = [f for f in fastq_files if not is_r1_file(f)
                      and not any(p in f for p in [".R2.", "_R2.", "_R2_", ".r2.", "_r2."])]

    for r1_file in r1_files:
        r2_file = find_r2_for_r1(r1_file)
        sample_name = os.path.basename(r1_file)
        # strip all common suffixes to get sample name
        for pat in [".R1.fq.gz", ".R1.fastq.gz", "_R1.fq.gz", "_R1.fastq.gz",
                    ".r1.fq.gz", ".r1.fastq.gz"]:
            if sample_name.endswith(pat):
                sample_name = sample_name[:-len(pat)]
                break
        else:
            sample_name = sample_name.split(".")[0]

        out_r1 = os.path.join(temp_dir, f"{sample_name}.clean.R1.fq.gz")
        out_r2 = os.path.join(temp_dir, f"{sample_name}.clean.R2.fq.gz")
        json_report = os.path.join(temp_dir, f"{sample_name}.fastp.json")
        html_report = os.path.join(temp_dir, f"{sample_name}.fastp.html")

        # idempotency check
        if r2_file and os.path.exists(out_r1) and os.path.exists(out_r2):
            logger.info(f"fastp output exists, skipping: {sample_name}")
            processed_files.extend([out_r1, out_r2])
            continue
        if not r2_file and os.path.exists(out_r1):
            logger.info(f"fastp output exists, skipping: {sample_name}")
            processed_files.append(out_r1)
            continue

        if r2_file:
            cmd = [
                "fastp",
                "-i", r1_file, "-I", r2_file,
                "-o", out_r1, "-O", out_r2,
                "-q", str(min_quality),
                "--length_required", str(min_length),
                "--json", json_report,
                "--html", html_report,
                "--thread", n_threads,
            ]
            logger.info(f"processing paired sample ({n_threads} threads): {sample_name} (R1: {r1_file}, R2: {r2_file})")
            if plasma:
                cmd += [
                    "--umi", "--umi_loc", "read2",
                    "--umi_len", str(umi_length),
                    "--umi_prefix", "UMI",
                    "--trim_front1", str(trim_front1),
                    "--trim_front2", str(trim_front2),
                    "--trim_tail2", str(trim_tail2),
                ]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"fastp paired done: {sample_name}")
                if os.path.exists(out_r1) and os.path.exists(out_r2):
                    processed_files.extend([out_r1, out_r2])
                else:
                    logger.error(f"fastp did not produce expected output: {sample_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"fastp failed ({sample_name}): {e.stderr}")
        else:
            # single-end
            cmd = [
                "fastp",
                "-i", r1_file,
                "-o", out_r1,
                "-q", str(min_quality),
                "--length_required", str(min_length),
                "--json", json_report,
                "--html", html_report,
                "--thread", n_threads,
            ]
            if plasma:
                cmd += ["--trim_front1", str(trim_front1)]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                logger.info(f"fastp single-end done: {sample_name}")
                if os.path.exists(out_r1):
                    processed_files.append(out_r1)
                else:
                    logger.error(f"fastp did not produce expected output: {sample_name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"fastp failed ({sample_name}): {e.stderr}")

    for fastq_file in unpaired_files:
        sample_name = os.path.basename(fastq_file).split(".")[0]
        out_r1 = os.path.join(temp_dir, f"{sample_name}.clean.fq.gz")
        json_report = os.path.join(temp_dir, f"{sample_name}.fastp.json")
        html_report = os.path.join(temp_dir, f"{sample_name}.fastp.html")
        if os.path.exists(out_r1):
            logger.info(f"fastp output exists, skipping: {sample_name}")
            processed_files.append(out_r1)
            continue
        cmd = [
            "fastp", "-i", fastq_file, "-o", out_r1,
            "-q", str(min_quality), "--length_required", str(min_length),
            "--json", json_report, "--html", html_report, "--thread", n_threads,
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if os.path.exists(out_r1):
                processed_files.append(out_r1)
        except subprocess.CalledProcessError as e:
            logger.error(f"fastp failed ({fastq_file}): {e.stderr}")

    logger.info(f"Preprocessing done — {len(processed_files)} file(s) produced")
    return processed_files
