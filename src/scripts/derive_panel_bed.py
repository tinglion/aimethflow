#!/usr/bin/env python3
"""
derive_panel_bed.py - Derive a targeted-panel BED file from FASTQ alignments

Given a set of representative FASTQ files (or pre-aligned BAMs), this script:
  1. Aligns each R1/R2 pair with bwameth (skipped when --bam is used)
  2. Merges all BAMs into one
  3. Computes per-base coverage with bedtools genomecov
  4. Retains regions with depth >= min_depth (default 10×)
  5. Merges nearby regions with bedtools merge -d merge_dist (default 100 bp)
  6. Filters to standard chromosomes (chr1–22, chrX, chrY)
  7. Writes the resulting intervals as a 3-column BED file

Typical use (first-time panel derivation, 3–5 representative samples):
    python src/scripts/derive_panel_bed.py \\
        --fq /data/plasma/B10001.R1.fq.gz /data/plasma/B10002.R1.fq.gz \\
        --out data/panel_B10308.bed

Or from BAM files directly:
    python src/scripts/derive_panel_bed.py \\
        --bam temp/run1/*.bwameth_pe.sorted.bam \\
        --out data/panel_B10308.bed
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

_SRC = os.path.join(os.path.dirname(__file__), "..")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from utils.logger import logger


# Standard chromosomes to retain
_KEEP_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def _resolve_tool(name, config):
    conda_bin = (config or {}).get("tools", {}).get("conda_env_bin", "")
    if conda_bin:
        candidate = os.path.join(conda_bin, name)
        if os.path.exists(candidate):
            return candidate
    found = shutil.which(name)
    return found or name


def _pair_fastq(fastq_files):
    """Pair R1/R2 FASTQ files. Returns list of (r1, r2) tuples."""
    r1 = sorted(f for f in fastq_files if
                ".R1." in f or "_R1." in f or "_val_1." in f or ".1.fq" in f)
    r2 = sorted(f for f in fastq_files if
                ".R2." in f or "_R2." in f or "_val_2." in f or ".2.fq" in f)
    if not r1:
        # fall back: interleaved — treat all as R1, no R2
        return [(f, None) for f in sorted(fastq_files)]
    if len(r1) != len(r2):
        logger.warning(f"R1 count ({len(r1)}) ≠ R2 count ({len(r2)}) — using R1-only pairs")
        return [(f, None) for f in r1]
    return list(zip(r1, r2))


def _align_pair(r1, r2, sample_name, tmpdir, ref, threads, config):
    """Align one R1/R2 pair with bwameth; return sorted BAM path or None."""
    out_bam = os.path.join(tmpdir, f"{sample_name}.sorted.bam")
    if os.path.exists(out_bam) and os.path.exists(out_bam + ".bai"):
        logger.info(f"  BAM exists, reusing: {out_bam}")
        return out_bam

    bwameth = _resolve_tool("bwameth.py", config)
    samtools = _resolve_tool("samtools", config)

    unsorted = out_bam.replace(".sorted.bam", ".unsorted.bam")
    log_path = out_bam.replace(".sorted.bam", ".align.log")

    sort_t   = str((config or {}).get("threads", {}).get("samtools_sort", 4))
    sort_mem = str((config or {}).get("threads", {}).get("samtools_sort_mem", "4G"))

    cmd_align = [bwameth, "--reference", ref, "-t", str(threads), r1]
    if r2:
        cmd_align.append(r2)
    cmd_view = [samtools, "view", "-bS", "-o", unsorted, "-"]

    logger.info(f"  Aligning {sample_name}: {' '.join(cmd_align)}")
    try:
        with open(log_path, "w") as flog:
            p_align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE, stderr=flog)
            p_view  = subprocess.Popen(cmd_view, stdin=p_align.stdout,
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            p_align.stdout.close()
            p_view.wait()
            p_align.wait()

        if p_align.returncode != 0:
            logger.error(f"  bwameth failed (exit {p_align.returncode}), see {log_path}")
            return None
        if p_view.returncode != 0:
            logger.error(f"  samtools view failed (exit {p_view.returncode})")
            return None

        subprocess.run(
            [samtools, "sort", "-@", sort_t, "-m", sort_mem, "-o", out_bam, unsorted],
            check=True, capture_output=True,
        )
        subprocess.run([samtools, "index", out_bam], check=True, capture_output=True)
        os.remove(unsorted)
        return out_bam
    except Exception as exc:
        logger.error(f"  Alignment error for {sample_name}: {exc}")
        return None


def _coverage_to_bed(bam_list, out_bed, min_depth, merge_dist, config):
    """
    Compute per-sample coverage, union covered intervals, sort and merge.

    No merged BAM is written — each BAM is streamed through bedtools genomecov
    independently, covered intervals are collected in memory (typically a few MB
    for targeted panels), sorted, then piped to a single bedtools merge.
    """
    import threading

    bedtools = _resolve_tool("bedtools", config)
    logger.info(
        f"  genomecov ×{len(bam_list)} BAM(s), depth≥{min_depth}, merge -d {merge_dist}"
    )

    intervals = []
    for bam in bam_list:
        p = subprocess.Popen(
            [bedtools, "genomecov", "-ibam", bam, "-bg"],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        for raw in p.stdout:
            parts = raw.split(b"\t")
            if len(parts) >= 4:
                try:
                    if int(parts[3]) >= min_depth:
                        intervals.append((parts[0].decode(), int(parts[1]), int(parts[2])))
                except ValueError:
                    pass
        p.wait()

    if not intervals:
        logger.error("No intervals above depth threshold in any BAM")
        return False

    intervals.sort(key=lambda x: (x[0], x[1]))
    logger.info(f"  {len(intervals):,} raw intervals collected, merging...")

    cmd_merge = [bedtools, "merge", "-d", str(merge_dist)]
    p_mrg = subprocess.Popen(
        cmd_merge, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )

    # Feed intervals in a thread to avoid deadlock while reading stdout
    def _feed():
        for chrom, start, end in intervals:
            p_mrg.stdin.write(f"{chrom}\t{start}\t{end}\n".encode())
        p_mrg.stdin.close()

    t = threading.Thread(target=_feed, daemon=True)
    t.start()

    kept = 0
    os.makedirs(os.path.dirname(os.path.abspath(out_bed)), exist_ok=True)
    with open(out_bed, "w") as fout:
        for raw in p_mrg.stdout:
            line = raw.decode()
            if line.split("\t")[0] in _KEEP_CHROMS:
                fout.write(line)
                kept += 1

    t.join()
    p_mrg.wait()
    logger.info(f"  → {kept} BED intervals written to {out_bed}")
    return kept > 0


def derive_panel_bed(fastq_files, output_bed, config=None, workdir=None,
                     min_depth=10, merge_dist=100, max_samples=None):
    """
    Derive a panel BED from FASTQ files.

    Args:
        fastq_files:  List of FASTQ paths (R1 + R2 interleaved, or all R1s — pairing is auto-detected)
        output_bed:   Output BED file path
        config:       Loaded config dict (uses load_config() if None)
        workdir:      Working directory for intermediate files (auto-created temp dir if None)
        min_depth:    Minimum coverage depth threshold (default 10)
        merge_dist:   Gap to merge adjacent intervals (bp, default 100)
        max_samples:  Use only the first N pairs (for quick panel derivation)

    Returns:
        True on success, False on failure
    """
    if config is None:
        from utils.config import load_config
        config = load_config()

    ref = config.get("extraction", {}).get("reference_genome", "")
    if not ref or not os.path.exists(ref):
        logger.error(f"Reference genome not found: {ref!r} — set [extraction] reference_genome in config")
        return False

    threads = config.get("threads", {}).get("bwameth", 8)

    own_tmpdir = workdir is None
    if own_tmpdir:
        tmpdir = tempfile.mkdtemp(prefix="panel_derive_")
    else:
        tmpdir = os.path.join(workdir, "panel_derive")
        os.makedirs(tmpdir, exist_ok=True)

    try:
        pairs = _pair_fastq(fastq_files)
        if max_samples:
            pairs = pairs[:max_samples]
        logger.info(f"Panel derivation: {len(pairs)} sample pair(s), ref={ref}")

        bam_list = []
        for r1, r2 in pairs:
            name = os.path.basename(r1).split(".")[0]
            bam = _align_pair(r1, r2, name, tmpdir, ref, threads, config)
            if bam:
                bam_list.append(bam)

        if not bam_list:
            logger.error("No BAMs produced — panel derivation failed")
            return False

        ok = _coverage_to_bed(bam_list, output_bed, min_depth, merge_dist, config)
        if ok:
            logger.info(f"✓ Panel BED written: {output_bed}")
        return ok

    finally:
        if own_tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)

def derive_panel_bed_from_bams(bam_files, output_bed, config=None,
                               min_depth=10, merge_dist=100, max_samples=5):
    """Derive panel BED directly from pre-aligned BAM files (no temp files written)."""
    if config is None:
        from utils.config import load_config
        config = load_config()

    bam_list = list(bam_files)
    if max_samples:
        bam_list = bam_list[:max_samples]
    logger.info(f"Panel derivation bam files: {bam_list}")
    ok = _coverage_to_bed(bam_list, output_bed, min_depth, merge_dist, config)
    if ok:
        logger.info(f"✓ Panel BED written: {output_bed}")
    return ok


def _bams_from_state(workdir):
    """
    Read pipeline_state.json and return the best available BAM list.

    Priority (highest first — most processed):
      deduplicated_files  (after step 6 dedup)
      bam_files           (after step 5 align)
      target_files        (after step 7 target)
    """
    import json
    state_path = os.path.join(workdir, "pipeline_state.json")
    if not os.path.exists(state_path):
        logger.error(f"pipeline_state.json not found in {workdir}")
        return None, None

    with open(state_path, encoding="utf-8") as f:
        state = json.load(f)

    for key in ("deduplicated_files", "bam_files", "target_files"):
        bams = state.get(key)
        if bams:
            existing = [b for b in bams if os.path.exists(b)]
            if existing:
                logger.info(f"Using {key} from pipeline state ({len(existing)} BAM(s))")
                return existing, key

    logger.error(
        f"No usable BAMs in {state_path}. "
        "Run at least '--step align' first, or pass --fq / --bam directly."
    )
    return None, None


def main():
    _SRC2 = os.path.join(os.path.dirname(__file__), "..")
    if _SRC2 not in sys.path:
        sys.path.insert(0, _SRC2)
    from utils.config import load_config

    parser = argparse.ArgumentParser(
        description="Derive a targeted-panel BED from FASTQ or BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    src_group = parser.add_mutually_exclusive_group(required=True)
    src_group.add_argument(
        "--fq", nargs="+", metavar="FILE",
        help="FASTQ files (R1 and R2 interleaved). Pairing is auto-detected from filenames.",
    )
    src_group.add_argument(
        "--bam", nargs="+", metavar="FILE",
        help="Pre-aligned BAM files — skips alignment entirely.",
    )
    src_group.add_argument(
        "--workdir", metavar="DIR",
        help="Existing pipeline workdir. BAMs are read from pipeline_state.json "
             "(priority: deduplicated > aligned > on-target). "
             "Example: --workdir temp/20260426_120000",
    )
    parser.add_argument(
        "--out", default="data/panel_B10308.bed", metavar="BED",
        help="Output BED file (default: data/panel_B10308.bed)",
    )
    parser.add_argument(
        "--ref", metavar="FASTA",
        help="Reference genome FASTA (overrides config [extraction] reference_genome)",
    )
    parser.add_argument(
        "--min-depth", type=int, default=10, metavar="INT",
        help="Minimum per-base coverage to include a site (default: 10)",
    )
    parser.add_argument(
        "--merge-dist", type=int, default=100, metavar="INT",
        help="Gap (bp) within which adjacent intervals are merged (default: 100)",
    )
    parser.add_argument(
        "--threads", type=int, default=None, metavar="INT",
        help="bwameth alignment threads — only used with --fq (default: from config)",
    )
    parser.add_argument(
        "--samples", type=int, default=None, metavar="N",
        help="Use only the first N FASTQ pairs — only used with --fq (default: all)",
    )
    parser.add_argument(
        "--tmpdir", metavar="DIR",
        help="Directory for intermediate files (default: system temp)",
    )
    parser.add_argument(
        "--config", metavar="TOML",
        help="Config file (default: auto-detect config.toml)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config) if args.config else load_config()
    if args.ref:
        cfg.setdefault("extraction", {})["reference_genome"] = args.ref
    if args.threads:
        cfg.setdefault("threads", {})["bwameth"] = args.threads

    if args.workdir:
        bam_list, source = _bams_from_state(args.workdir)
        if not bam_list:
            sys.exit(1)
        logger.info(f"Source: {source} from {args.workdir}")
        ok = derive_panel_bed_from_bams(
            bam_list, args.out, config=cfg,
            min_depth=args.min_depth, merge_dist=args.merge_dist,
            max_samples=args.samples,
        )
    elif args.bam:
        ok = derive_panel_bed_from_bams(
            args.bam, args.out, config=cfg,
            min_depth=args.min_depth, merge_dist=args.merge_dist,
            max_samples=args.samples,
        )
    else:
        ok = derive_panel_bed(
            args.fq, args.out, config=cfg, workdir=args.tmpdir,
            min_depth=args.min_depth, merge_dist=args.merge_dist,
            max_samples=args.samples,
        )

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
