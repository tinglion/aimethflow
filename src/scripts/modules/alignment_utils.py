#!/usr/bin/env python3
"""
alignment_utils.py - Sequence alignment module

Supports two alignment tools (select via config.toml [plasma] aligner):
  bismark  — accurate, slow (~3-4h / 30M pairs). Use for initial validation
             or when OT/OB directionality reports are needed.
  bwameth  — 3-5x faster (~40-60min / 30M pairs). Recommended for batches
             of 100+ samples.

Methylation extraction:
  bismark  → bismark_methylation_extractor (.bismark.cov)
  bwameth  → MethylDackel extract (_CpG.bedGraph)
"""

import json
import os
import re
import shutil
import subprocess
import traceback

from utils.logger import logger


def _rm(path):
    """Safely remove a single file; log a warning on failure instead of raising."""
    try:
        if os.path.exists(path):
            os.remove(path)
            logger.debug(f"  [cleanup] removing intermediate file: {path}")
    except OSError as e:
        logger.warning(f"  [cleanup] failed to remove: {path} — {e}")


def _resolve_tool(name, config):
    """
    Resolve the absolute path to a tool executable.
    Priority order:
      1. config [tools] conda_env_bin/{name} (explicit absolute path, most reliable)
      2. shutil.which(name) (tool is on PATH)
      3. return name as-is (let subprocess report the error)
    """
    conda_bin = (config or {}).get("tools", {}).get("conda_env_bin", "")
    if conda_bin:
        candidate = os.path.join(conda_bin, name)
        if os.path.exists(candidate):
            return candidate
    found = shutil.which(name)
    if found:
        return found
    return name


# ──────────────────────────────────────────────
# Internal helpers: filename parsing
# ──────────────────────────────────────────────

def _find_paired_fastq(fastq_files):
    """
    Find R1/R2 pairs from a file list, compatible with multiple naming conventions.
    Prefers .clean.R1/.clean.R2 (fastp output); also handles _val_1/_val_2
    (TrimGalore output).

    Returns:
        list of (r1, r2_or_None)
    """
    pairs = []
    r2_seen = set()
    for f in fastq_files:
        if f in r2_seen:
            continue
        r2 = None
        for r1_pat, r2_pat in [
            (".clean.R1.", ".clean.R2."),
            (".R1.", ".R2."),
            ("_R1.", "_R2."),
            ("_R1_", "_R2_"),
            ("_val_1", "_val_2"),
        ]:
            if r1_pat in f:
                candidate = f.replace(r1_pat, r2_pat)
                if candidate in fastq_files:
                    r2 = candidate
                    r2_seen.add(candidate)
                break
        pairs.append((f, r2))
    return pairs


def _sample_base_from_r1(r1_file):
    """Extract the sample base name from an R1 filename (strips all known suffixes)."""
    base = os.path.basename(r1_file)
    for pat in [".clean.R1.fq.gz", ".R1.fq.gz", "_R1.fq.gz",
                ".clean.R1.fastq.gz", ".R1.fastq.gz", "_val_1.fq.gz", "_trimmed.fq.gz"]:
        if base.endswith(pat):
            return base[:-len(pat)]
    return base.split(".")[0]


# ──────────────────────────────────────────────
# bismark alignment (precise mode)
# ──────────────────────────────────────────────

def _align_bismark_pair(r1_file, r2_file, sample_base, temp_dir, config, non_directional=False):
    """
    bismark --bowtie2 paired-end alignment.
    Speed: ~3-4h / 30M pairs (8 cores).
    Advantage: outputs OT/OB/CTOT/CTOB strand-ratio report for library directionality validation.
    """
    bismark_index = os.path.dirname(config.get("extraction", {}).get("reference_genome", ""))
    cmd = [
        _resolve_tool("bismark", config), "--bowtie2", "--bam", "--multicore", "4",
        bismark_index, "-1", r1_file, "-2", r2_file,
        "--output_dir", temp_dir,
    ]
    if non_directional:
        cmd.append("--non_directional")

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        for fname in sorted(os.listdir(temp_dir)):
            if sample_base in fname and fname.endswith(".bam") and "bismark_bt2" in fname and "sorted" not in fname and "temp" not in fname:
                bam_path = os.path.join(temp_dir, fname)
                sorted_bam = bam_path.replace(".bam", ".sorted.bam")
                subprocess.run([_resolve_tool("samtools", config), "sort", "-o", sorted_bam, bam_path],
                                check=True, capture_output=True)
                subprocess.run([_resolve_tool("samtools", config), "index", sorted_bam],
                                check=True, capture_output=True)
                logger.info(f"bismark alignment + sort done: {sorted_bam}")
                _rm(bam_path)  # remove unsorted intermediate BAM
                for fname2 in os.listdir(temp_dir):
                    if sample_base in fname2 and "PE_report" in fname2:
                        logger.info(f"  ↳ directionality report (check OT/OB ratios): {os.path.join(temp_dir, fname2)}")
                return sorted_bam
        logger.error(f"bismark did not produce expected BAM: {sample_base}")
        return None
    except subprocess.CalledProcessError as e:
        logger.error(f"bismark alignment failed ({r1_file}): {e.stderr[-500:]}")
        return None


# ──────────────────────────────────────────────
# bwameth alignment (fast mode, recommended for batches)
# ──────────────────────────────────────────────

def _align_bwameth_pair(r1_file, r2_file, sample_base, temp_dir, config, non_directional=False):  # noqa: ARG001
    """
    bwameth paired-end alignment (BWA-MEM + C→T conversion).
    Speed: ~40-60min / 30M pairs (8 cores), 3-5x faster than bismark.
    Note: does not output OT/OB strand-ratio reports; validate directionality
    with bismark first before switching to bwameth for batch runs.

    Index requirement — run once before first use:
        bwameth.py index /data/meth/bwa/hg38.fa
    """
    reference_genome = config.get("extraction", {}).get("reference_genome", "")
    sorted_bam = os.path.join(temp_dir, f"{sample_base}.bwameth_pe.sorted.bam")

    # thread counts from config [threads]
    threads_cfg = config.get("threads", {})
    bwameth_t   = str(threads_cfg.get("bwameth", 8))
    sort_t      = str(threads_cfg.get("samtools_sort", 4))
    sort_mem    = str(threads_cfg.get("samtools_sort_mem", "8G"))

    # idempotency check: skip if already present
    if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
        logger.info(f"bwameth BAM exists, skipping: {sorted_bam}")
        return sorted_bam

    cmd_align   = ["bwameth.py", "--reference", reference_genome, "-t", bwameth_t,
                    r1_file, r2_file]
    unsorted_bam = sorted_bam.replace(".sorted.bam", ".unsorted.bam")

    align_log = sorted_bam.replace(".sorted.bam", ".align.log")
    try:
        # step 1: bwameth → samtools view (streaming SAM→BAM)
        # write stderr to file to prevent bwameth blocking stdout when PIPE buffer fills
        cmd_view   = ["samtools", "view", "-bS", "-o", unsorted_bam, "-"]
        logger.info(f"  [1/2] bwameth aligning (output to unsorted BAM)")
        with open(align_log, "w") as f_log:
            logger.info(f"    bwameth cmd: {' '.join(cmd_align)}")
            proc_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=f_log)
            logger.info(f"    samtools view cmd: {' '.join(cmd_view)}")
            proc_view  = subprocess.Popen(cmd_view, stdin=proc_align.stdout,
                                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            proc_align.stdout.close()
            proc_view.wait()
            proc_align.wait()

        if proc_view.returncode != 0:
            raise RuntimeError(f"samtools view failed, exit code: {proc_view.returncode}")
        if proc_align.returncode != 0:
            raise RuntimeError(f"bwameth failed, exit code: {proc_align.returncode}, log: {align_log}")

        # step 2: sort from file (no pipe; bwameth has exited — deadlock impossible)
        logger.info(f"  [2/2] samtools sort running...")
        subprocess.run(
            ["samtools", "sort", "-@", sort_t, "-m", sort_mem, "-o", sorted_bam, unsorted_bam],
            check=True, capture_output=True
        )
        _rm(unsorted_bam)

        subprocess.run(["samtools", "index", sorted_bam], check=True, capture_output=True)
        _rm(align_log)
        logger.info(f"bwameth alignment + sort done: {sorted_bam}")
        return sorted_bam
    except Exception as e:
        logger.error(f"bwameth alignment failed ({r1_file}): {e}")
        traceback.print_exc()
        _rm(unsorted_bam)
        return None


# ──────────────────────────────────────────────
# Public interface: align_methylation
# ──────────────────────────────────────────────

def align_methylation(fastq_files, temp_dir, config, plasma=False):
    """
    Methylation sequence alignment.

    Alignment tool is determined by config.toml [plasma] aligner:
      aligner = "bismark"  — slow but provides directionality report
      aligner = "bwameth"  — 3-5x faster, recommended for batches
                             (requires a pre-built c2t index)

    Args:
        fastq_files: preprocessed FastQ file list
        temp_dir:    working directory for intermediate files
        config:      config dict
        plasma:      True for plasma mode

    Returns:
        list: coordinate-sorted, indexed BAM file paths
    """
    aligner = config.get("plasma", {}).get("aligner", "bismark") if plasma else "bismark"
    non_directional = config.get("plasma", {}).get("non_directional", False)
    logger.info(f"Starting methylation alignment [tool: {aligner}, "
                f"{'plasma' if plasma else 'standard'} mode"
                f"{', non-directional' if plasma and non_directional else ''}]")

    bam_files = []
    paired_files = _find_paired_fastq(fastq_files)

    for r1_file, r2_file in paired_files:
        logger.info(f"processing: {r1_file}")
        if r2_file is None:
            logger.warning(f"no R2 found, skipping unpaired file: {r1_file}")
            continue
        sample_base = _sample_base_from_r1(r1_file)
        logger.info(f"  aligning: {sample_base}  ({r1_file})")

        if aligner == "bwameth":
            bam = _align_bwameth_pair(r1_file, r2_file, sample_base,
                                      temp_dir, config, non_directional)
        else:
            bam = _align_bismark_pair(r1_file, r2_file, sample_base,
                                      temp_dir, config, non_directional)
        if bam:
            bam_files.append(bam)
        logger.info(f"finished {r1_file}")

    logger.info(f"Alignment done — {len(bam_files)} BAM file(s) produced")
    return bam_files


# ──────────────────────────────────────────────
# Deduplication
# ──────────────────────────────────────────────

def deduplicate_bam(bam_files, temp_dir, plasma=False, config=None):
    """
    PCR duplicate removal.
    Plasma mode: umi_tools dedup (UMI-aware, avoids discarding genuine cfDNA fragments).
    Standard mode: deduplicate_bismark (coordinate-based deduplication).
    """
    logger.info(f"Starting PCR deduplication ({'umi_tools plasma mode' if plasma else 'bismark standard mode'})")
    delete_prev = (config or {}).get("cleanup", {}).get("delete_after_dedup", False)
    deduplicated_files = []

    for bam_file in bam_files:
        logger.info(f"processing: {bam_file}")
        if plasma:
            dedup_bam = bam_file.replace(".sorted.bam", ".dedup.bam")
            if dedup_bam == bam_file:
                dedup_bam = bam_file.replace(".bam", ".dedup.bam")
            if os.path.exists(dedup_bam) and os.path.exists(dedup_bam + ".bai"):
                logger.info(f"dedup BAM exists, skipping: {dedup_bam}")
                deduplicated_files.append(dedup_bam)
                continue
            stats_file = dedup_bam.replace(".bam", ".stats")
            cmd = ["umi_tools", "dedup", "--paired",
                   "--stdin", bam_file, "--stdout", dedup_bam, "--log", stats_file]
            try:
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                subprocess.run(["samtools", "index", dedup_bam],
                                check=True, capture_output=True)
                deduplicated_files.append(dedup_bam)
                logger.info(f"umi_tools dedup done: {dedup_bam}")
                if delete_prev:
                    _rm(bam_file)
                    _rm(bam_file + ".bai")
            except subprocess.CalledProcessError as e:
                logger.error(f"umi_tools dedup failed: {e.stderr[-300:]}")
                logger.warning("falling back to deduplicate_bismark (not recommended for plasma data)")
                deduplicated_files.extend(_deduplicate_bismark([bam_file], temp_dir))
        else:
            deduplicated_files.extend(_deduplicate_bismark([bam_file], temp_dir))

    logger.info(f"Deduplication done — {len(deduplicated_files)} deduplicated BAM file(s)")
    return deduplicated_files


def _deduplicate_bismark(bam_files, temp_dir):
    """Coordinate-based deduplication via deduplicate_bismark (tissue mode)."""
    result = []
    for bam_file in bam_files:
        dedup_file = bam_file.replace(".bam", ".deduplicated.bam")
        if os.path.exists(dedup_file):
            logger.info(f"deduplicate_bismark BAM exists, skipping: {dedup_file}")
            result.append(dedup_file)
            continue
        cmd = ["deduplicate_bismark", "--bam", bam_file, "--output_dir", temp_dir]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if os.path.exists(dedup_file):
                result.append(dedup_file)
                logger.info(f"deduplicate_bismark done: {dedup_file}")
            else:
                logger.error(f"deduplicate_bismark did not produce expected output: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"deduplicate_bismark failed: {e.stderr[-300:]}")
    return result


# ──────────────────────────────────────────────
# On-target filtering
# ──────────────────────────────────────────────

def extract_target_regions(bam_files, temp_dir, config):
    """bedtools intersect: retain only reads overlapping the panel BED."""
    logger.info("Starting on-target extraction")
    target_files = []
    delete_prev = (config or {}).get("cleanup", {}).get("delete_after_target", False)

    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if not probe_panel or not os.path.exists(probe_panel):
        logger.warning("Panel BED not found, skipping on-target extraction")
        return bam_files

    for bam_file in bam_files:
        target_file_sorted = os.path.join(
            temp_dir, os.path.basename(bam_file).replace(".bam", ".targeted.sorted.bam")
        )
        target_file_sorted_bai = os.path.join(
            temp_dir, os.path.basename(bam_file).replace(".bam", ".targeted.sorted.bam.bai")
        )
        if os.path.exists(target_file_sorted) and os.path.exists(target_file_sorted_bai):
            logger.info(f"targeted BAM exists, skipping: {target_file_sorted}")
            target_files.append(target_file_sorted)
            continue

        logger.info(f"processing: {bam_file}")
        output_file = os.path.join(
            temp_dir, os.path.basename(bam_file).replace(".bam", ".targeted.bam")
        )
        cmd = ["bedtools", "intersect", "-abam", bam_file,
               "-b", probe_panel, "-wa", "-header"]
        try:
            with open(output_file, "wb") as fout:
                subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.PIPE)
            sorted_file = output_file.replace(".bam", ".sorted.bam")
            subprocess.run(["samtools", "sort", output_file, "-o", sorted_file],
                            check=True, capture_output=True)
            subprocess.run(["samtools", "index", sorted_file],
                            check=True, capture_output=True)
            if os.path.exists(sorted_file + ".bai"):
                target_files.append(sorted_file)
                logger.info(f"on-target extraction done: {sorted_file}")
                _rm(output_file)       # remove unsorted intermediate targeted.bam
                if delete_prev:
                    _rm(bam_file)
                    _rm(bam_file + ".bai")
            else:
                logger.error(f"on-target extraction failed: {bam_file}")
        except subprocess.CalledProcessError as e:
            logger.error(f"bedtools intersect failed: {e.stderr}")

    logger.info(f"On-target extraction done — {len(target_files)} targeted BAM file(s)")
    return target_files


# ──────────────────────────────────────────────
# M-bias helpers
# ──────────────────────────────────────────────

def _get_read_length_from_bam(bam_file: str, samtools: str = "samtools") -> int | None:
    """Return read length by reading the SEQ field of the first mapped alignment."""
    try:
        proc = subprocess.Popen(
            [samtools, "view", "-F", "4", bam_file],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        )
        line = proc.stdout.readline()
        proc.kill()
        proc.wait()
        if line:
            fields = line.split("\t")
            if len(fields) > 9:
                return len(fields[9])
    except Exception as exc:
        logger.debug(f"  Could not determine read length from BAM: {exc}")
    return None


def _parse_mbias_suggestion(stderr: str, read_len: int = None) -> dict:
    """Parse MethylDackel's suggested trim values from its stderr output.

    MethylDackel writes a single line to stderr in *inclusion-position* format:
        Suggested inclusion options: --OT a,b,c,d --OB a,b,c,d

    Meaning of the four 0-based values (a, b, c, d):
        a  R1 5' start position  → skip_5'_R1 = a
        b  R1 3' end   position  → skip_3'_R1 = 0 if b==0 else read_len - b - 1
        c  R2 5' start position  → skip_5'_R2 = c
        d  R2 3' end   position  → skip_3'_R2 = 0 if d==0 else read_len - d - 1

    Example: --OT 0,0,0,110 with read_len=150 → --nOT 0,0,0,39
             --OB 0,0,17,0               → --nOB 0,0,17,0  (no 3' conversion needed)

    Returns {nOT: "R1_5,R1_3,R2_5,R2_3", nOB: "..."} as --nOT/--nOB skip-count
    strings ready for MethylDackel extract.  If read_len is None and a 3' position
    is non-zero, stores the raw position value with a warning (incorrect for extract).
    """
    result = {}
    for md_flag, cfg_key in (("OT", "nOT"), ("OB", "nOB")):
        m = re.search(rf"--{md_flag}\s+([\d,]+)", stderr, re.IGNORECASE)
        if not m:
            continue
        try:
            vals = [int(x) for x in m.group(1).split(",") if x.strip()]
        except ValueError:
            continue
        if len(vals) != 4:
            continue
        r1_5, r1_3_pos, r2_5, r2_3_pos = vals
        if read_len:
            r1_3_skip = 0 if r1_3_pos == 0 else read_len - r1_3_pos - 1
            r2_3_skip = 0 if r2_3_pos == 0 else read_len - r2_3_pos - 1
        else:
            if r1_3_pos > 0 or r2_3_pos > 0:
                logger.warning(
                    f"  [{cfg_key}] 3' positions (R1={r1_3_pos}, R2={r2_3_pos}) require "
                    "read_len for correct conversion to --nOT skip counts — stored as raw "
                    "positions (INCORRECT for extract). Fix: skip_3' = read_len - pos - 1"
                )
            r1_3_skip = r1_3_pos
            r2_3_skip = r2_3_pos
        result[cfg_key] = f"{r1_5},{r1_3_skip},{r2_5},{r2_3_skip}"
    return result


def _aggregate_trim(
    trim_list: list,
    default: str = "0,0,0,0",
    percentile: float = 0.95,
) -> tuple[str, list[dict]]:
    """Per-position p-th percentile of trim values across samples.
    Returns (trim_str, per_position_stats).

    Using a percentile instead of max so that a handful of low-quality BAMs with
    extreme end bias do not inflate the result for all samples.  The 95th percentile
    keeps the value conservative (≥95% of samples are fully covered) while excluding
    the worst outliers.  Configurable via config.toml [plasma] mbias_trim_percentile.
    """
    _LABELS = ("R1_5'", "R1_3'", "R2_5'", "R2_3'")
    if not trim_list:
        return default, []
    try:
        arrays = [[int(x) for x in s.split(",") if x.strip()] for s in trim_list]
        if not arrays or any(len(a) != 4 for a in arrays):
            return default, []
        n = len(arrays)
        p_pct = int(round(percentile * 100))
        result, stats = [], []
        for col in range(4):
            vals    = sorted(a[col] for a in arrays)
            idx_p   = min(int(percentile * n + 0.9999), n - 1)  # ceiling → conservative
            idx_med = n // 2
            result.append(vals[idx_p])
            stats.append({
                "label":    _LABELS[col],
                "min":      vals[0],
                "median":   vals[idx_med],
                f"p{p_pct}": vals[idx_p],
                "max":      vals[-1],
                "n_outliers_excluded": n - 1 - idx_p,
            })
        return ",".join(str(v) for v in result), stats
    except (ValueError, IndexError):
        return default, []


def _resolve_auto_trim(value: str, key: str, workdir: str) -> str:
    """Return value as-is unless it equals 'auto', in which case read from
    mbias_auto_trim.json written by check_mbias().  Falls back to '0,0,0,0'
    with a warning if the file is missing or malformed."""
    if value.strip().lower() != "auto":
        return value
    auto_path = os.path.join(workdir, "mbias_auto_trim.json")
    if not os.path.exists(auto_path):
        logger.warning(
            f"  {key}='auto' but {auto_path} not found — defaulting to '0,0,0,0'. "
            "Run --step mbias first to generate the auto-trim file."
        )
        return "0,0,0,0"
    try:
        with open(auto_path) as f:
            data = json.load(f)
        resolved = data.get(key, "0,0,0,0")
        logger.info(f"  {key} resolved from auto-detection: '{resolved}'  (source: {auto_path})")
        return resolved
    except Exception as e:
        logger.warning(f"  Failed to read {auto_path}: {e} — defaulting to '0,0,0,0'")
        return "0,0,0,0"


# ──────────────────────────────────────────────
# M-bias check
# ──────────────────────────────────────────────

def check_mbias(bam_files, temp_dir, reference_genome, config=None):
    """
    MethylDackel mbias: check for position-level methylation bias across reads.

    Generates SVG plots per sample and parses MethylDackel's built-in suggested
    trim values from stdout.  After all BAMs are processed the per-sample
    suggestions are aggregated (element-wise maximum = most conservative) and
    saved to {temp_dir}/mbias_auto_trim.json.

    Usage with auto-detection:
        Set nOT = "auto" and/or nOB = "auto" in config.toml [plasma].
        The extract step will then load the detected values from
        mbias_auto_trim.json automatically.  Explicit numeric values
        (e.g. nOT = "0,0,0,20") are always used as-is regardless of this file.
    """
    logger.info("Starting M-bias check")
    n_threads = str((config or {}).get("threads", {}).get("methyldackel", 4))
    samtools  = _resolve_tool("samtools", config or {})
    report_files = []
    per_sample: list[dict] = []
    detected_read_lens: list[int] = []

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split(".")[0]
        output_prefix = os.path.join(temp_dir, f"{sample_name}.mbias")
        cmd = ["MethylDackel", "mbias", "--nOT", "0,0,0,0", "-@", n_threads,
               reference_genome, bam_file, output_prefix]
        try:
            read_len = _get_read_length_from_bam(bam_file, samtools)
            if read_len:
                detected_read_lens.append(read_len)
                logger.info(f"  Read length ({sample_name}): {read_len}bp")

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            svg_file = f"{output_prefix}_OT.svg"
            if os.path.exists(svg_file):
                report_files.append(svg_file)

            suggestion = _parse_mbias_suggestion(result.stderr, read_len=read_len)
            if suggestion:
                per_sample.append(suggestion)
                logger.info(
                    f"  M-bias suggestion ({sample_name}): "
                    f"--nOT {suggestion.get('nOT', 'n/a')}  "
                    f"--nOB {suggestion.get('nOB', 'n/a')}"
                )
            else:
                logger.warning(f"  No suggested trim found in MethylDackel stderr for {sample_name}")

            logger.info(
                f"M-bias check done: {sample_name}\n"
                f"  → SVG files: {output_prefix}_OT.svg / _OB.svg\n"
                f"  → Inspect SVGs; set nOT/nOB = \"auto\" in config.toml [plasma] to apply "
                f"auto-detected values, or fill in explicit values"
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"MethylDackel mbias failed ({sample_name}): {e.stderr[-300:]}")

    # Aggregate per-sample suggestions and save auto-trim file
    if per_sample:
        p     = float((config or {}).get("plasma", {}).get("mbias_trim_percentile", 95)) / 100
        p_pct = int(round(p * 100))
        agg_nOT, stats_nOT = _aggregate_trim(
            [s["nOT"] for s in per_sample if "nOT" in s], percentile=p)
        agg_nOB, stats_nOB = _aggregate_trim(
            [s["nOB"] for s in per_sample if "nOB" in s], percentile=p)

        for strand, _, stats in (("OT", agg_nOT, stats_nOT), ("OB", agg_nOB, stats_nOB)):
            logger.info(
                f"  M-bias {strand} trim distribution "
                f"(n={len(per_sample)} samples, p{p_pct} chosen):"
            )
            for st in stats:
                logger.info(
                    f"    {st['label']}: "
                    f"min={st['min']:4d}  median={st['median']:4d}  "
                    f"p{p_pct}={st[f'p{p_pct}']:4d}  max={st['max']:4d}  "
                    f"(outliers excluded: {st['n_outliers_excluded']})"
                )

        auto_path = os.path.join(temp_dir, "mbias_auto_trim.json")
        read_len_used = max(detected_read_lens) if detected_read_lens else None
        with open(auto_path, "w") as f:
            json.dump({
                "nOT": agg_nOT,
                "nOB": agg_nOB,
                "read_len": read_len_used,
                "n_samples": len(per_sample),
                "percentile_used": p_pct,
                "aggregation": (
                    f"element-wise {p_pct}th percentile across all samples "
                    f"(outliers above the {p_pct}th percentile excluded per position)"
                ),
                "nOT_distribution": stats_nOT,
                "nOB_distribution": stats_nOB,
                "note": (
                    f"--OT/--OB inclusion positions converted to --nOT/--nOB skip counts "
                    f"using read_len={read_len_used}. "
                    "Formula: skip_3' = read_len - pos - 1  (pos=0 → skip=0)."
                ),
                "usage": (
                    "Set nOT = \"auto\" and/or nOB = \"auto\" in config.toml [plasma] "
                    "to use these values. Adjust mbias_trim_percentile (default 95) to "
                    "control how aggressively outlier samples are excluded."
                ),
            }, f, indent=2)
        logger.info(
            f"  M-bias auto-trim saved: {auto_path}\n"
            f"    {p_pct}th percentile across {len(per_sample)} sample(s): "
            f"--nOT {agg_nOT}  --nOB {agg_nOB}\n"
            f"  → Set nOT = \"auto\" / nOB = \"auto\" in config.toml [plasma] to apply automatically"
        )
    else:
        logger.warning(
            "  M-bias: no suggested trim values found in MethylDackel output — "
            "mbias_auto_trim.json not written"
        )

    return report_files


# ──────────────────────────────────────────────
# cfDNA fragment size distribution
# ──────────────────────────────────────────────

def plot_fragment_size(bam_files, result_dir, config=None):
    """
    Plot insert-size distribution from paired BAM files (TLEN field).
    Samples up to 200 000 reads per BAM to stay fast.
    Saves fragment_size_distribution.png and fragment_size_stats.csv to result_dir.
    """
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        logger.warning(f"[fragment_size] Cannot import plotting libraries: {exc}")
        return None

    samtools  = _resolve_tool("samtools", config or {})
    MAX_READS = 200_000
    SIZE_MIN, SIZE_MAX = 0, 600

    all_sizes: list   = []
    sample_names: list = []
    stats_rows: list  = []

    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split(".")[0]
        sample_names.append(sample_name)
        sizes: list = []

        try:
            # Proper pairs only (-f 0x2); exclude unmapped / secondary / supplementary (-F 0x904)
            cmd  = [samtools, "view", "-f", "0x2", "-F", "0x904", bam_file]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, text=True)
            for line in proc.stdout:
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                try:
                    tlen = int(parts[8])
                except ValueError:
                    continue
                if tlen > 0:          # keep one end per pair (positive TLEN)
                    sizes.append(tlen)
                    if len(sizes) >= MAX_READS:
                        proc.stdout.close()
                        proc.wait()
                        break
            else:
                proc.wait()
        except Exception as exc:
            logger.warning(f"[fragment_size] samtools view failed for {sample_name}: {exc}")
            all_sizes.append([])
            stats_rows.append({"sample": sample_name, "n_reads": 0,
                                "median": "NA", "mean": "NA", "pct_100_220": "NA"})
            continue

        if sizes:
            arr = np.array(sizes, dtype=np.int32)
            arr = arr[(arr >= SIZE_MIN) & (arr <= SIZE_MAX)]
            in_peak = np.sum((arr >= 100) & (arr <= 220)) / max(len(arr), 1) * 100
            stats_rows.append({
                "sample":      sample_name,
                "n_reads":     len(arr),
                "median":      int(np.median(arr)),
                "mean":        round(float(np.mean(arr)), 1),
                "pct_100_220": round(in_peak, 1),
            })
            all_sizes.append(arr.tolist())
            logger.info(
                f"[fragment_size] {sample_name}: n={len(arr):,d}  "
                f"median={int(np.median(arr))}bp  "
                f"mean={np.mean(arr):.1f}bp  "
                f"pct_100-220bp={in_peak:.1f}%"
            )
        else:
            all_sizes.append([])
            stats_rows.append({"sample": sample_name, "n_reads": 0,
                                "median": "NA", "mean": "NA", "pct_100_220": "NA"})

    # ── plot ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5))
    bins     = np.arange(SIZE_MIN, SIZE_MAX + 2, 2)
    centres  = (bins[:-1] + bins[1:]) / 2
    agg_all: list = []

    n_valid = sum(1 for s in all_sizes if s)
    first_done = False
    for sizes_i in all_sizes:
        if not sizes_i:
            continue
        agg_all.extend(sizes_i)
        hist, _ = np.histogram(np.array(sizes_i), bins=bins, density=True)
        # Single combined legend entry for all per-sample curves; subsequent
        # lines suppressed from the legend so the box stays compact.
        lbl = f"Per-sample (n={n_valid})" if not first_done else None
        ax.plot(centres, hist, color="#3a86ff", alpha=0.3, linewidth=0.8, label=lbl)
        first_done = True

    if agg_all:
        hist_agg, _ = np.histogram(np.array(agg_all), bins=bins, density=True)
        ax.plot(centres, hist_agg, color="#e63946", linewidth=2.0, label="Aggregate")

    ax.axvspan(100, 220, color="#f4a261", alpha=0.15, label="cfDNA zone 100–220 bp")
    ax.axvline(167, color="#f4a261", linestyle="--", linewidth=1.0,
               label="167 bp mono-nucleosome")

    ax.set_xlabel("Fragment size (bp)", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    # ax.set_title("cfDNA Fragment Size Distribution", fontsize=13)
    ax.set_xlim(SIZE_MIN, SIZE_MAX)

    # Deduplicate legend: keep only the first occurrence of each label
    handles, labels = ax.get_legend_handles_labels()
    seen: dict = {}
    unique_h, unique_l = [], []
    for h, lbl in zip(handles, labels):
        if lbl not in seen:
            seen[lbl] = True
            unique_h.append(h)
            unique_l.append(lbl)
    ax.legend(unique_h, unique_l, fontsize=9, loc="upper right")

    os.makedirs(result_dir, exist_ok=True)
    png_path = os.path.join(result_dir, "fragment_size_distribution.png")
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    logger.info(f"[fragment_size] Plot saved: {png_path}")

    # ── stats CSV ──────────────────────────────────────────────────────
    if stats_rows:
        import csv as _csv
        csv_path = os.path.join(result_dir, "fragment_size_stats.csv")
        fieldnames = ["sample", "n_reads", "median", "mean", "pct_100_220"]
        with open(csv_path, "w", newline="") as f:
            writer = _csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(stats_rows)
        logger.info(f"[fragment_size] Stats saved: {csv_path}")

    return png_path


# ──────────────────────────────────────────────
# Methylation extraction: MethylDackel (for bwameth alignments)
# ──────────────────────────────────────────────

def extract_methylation_methyldackel(bam_files, temp_dir, config):
    """
    Extract CpG methylation with MethylDackel extract.
    Used with bwameth alignments; outputs bedGraph format.

    Output format (one CpG site per line):
        chrom  start  end  methylation%  meth_reads  total_reads

    Args:
        bam_files: on-target BAM files (sorted and indexed)
        temp_dir:  working directory for intermediate files
        config:    config dict

    Returns:
        list: bedGraph file paths ({sample}_CpG.bedGraph)
    """
    logger.info("Starting methylation extraction [MethylDackel extract]")
    n_threads = str(config.get("threads", {}).get("methyldackel", 4))
    plasma_cfg = config.get("plasma", {})
    nOT_raw = plasma_cfg.get("nOT", "0,0,0,0")
    nOB_raw = plasma_cfg.get("nOB", "0,0,0,0")
    nOT = _resolve_auto_trim(nOT_raw, "nOT", temp_dir)
    nOB = _resolve_auto_trim(nOB_raw, "nOB", temp_dir)
    nOT_src = "[auto-detected]" if nOT_raw.strip().lower() == "auto" else "[config]"
    nOB_src = "[auto-detected]" if nOB_raw.strip().lower() == "auto" else "[config]"
    logger.info(f"  M-bias trim: --nOT {nOT} {nOT_src}  --nOB {nOB} {nOB_src}")

    # First pass: collect already-extracted files and queue missing ones.
    # Reference genome is only needed for the missing subset.
    result_files = []
    pending = []
    for bam_file in bam_files:
        sample_name = os.path.basename(bam_file).split(".")[0]
        output_file = os.path.join(temp_dir, f"{sample_name}_CpG.bedGraph")
        if os.path.exists(output_file):
            logger.info(f"MethylDackel bedGraph exists, skipping: {output_file}")
            result_files.append(output_file)
        else:
            pending.append(bam_file)

    if not pending:
        logger.info(f"Methylation extraction done — {len(result_files)} bedGraph file(s) (all pre-existing)")
        return result_files

    reference_genome = config.get("extraction", {}).get("reference_genome", "")
    if not reference_genome or not os.path.exists(reference_genome):
        logger.error(f"reference genome not found: {reference_genome}")
        if result_files:
            logger.warning(f"  Returning {len(result_files)} pre-existing bedGraph(s) — {len(pending)} BAM(s) skipped")
            return result_files
        return []

    for bam_file in pending:
        sample_name = os.path.basename(bam_file).split(".")[0]
        output_prefix = os.path.join(temp_dir, sample_name)
        output_file = f"{output_prefix}_CpG.bedGraph"
        cmd = [
            _resolve_tool("MethylDackel", config), "extract",
            "--minDepth", "1",
            "--mergeContext",
            "--nOT", nOT,
            "--nOB", nOB,
            "-@", n_threads,
            reference_genome, bam_file,
            "-o", output_prefix,
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            if os.path.exists(output_file):
                result_files.append(output_file)
                logger.info(f"MethylDackel extraction done: {output_file}")
            else:
                logger.error(f"MethylDackel did not produce expected output: {output_prefix}_CpG.bedGraph")
        except subprocess.CalledProcessError as e:
            logger.error(f"MethylDackel extract failed ({sample_name}): {e.stderr[-300:]}")

    logger.info(f"Methylation extraction done — {len(result_files)} bedGraph file(s)")
    return result_files
