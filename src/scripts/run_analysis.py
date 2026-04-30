#!/usr/bin/env python3
"""
run_analysis.py - Full methylation analysis pipeline

Two execution modes:
  [Full pipeline]  Run all steps in sequence (default)
  [Step mode]      Run a single named step; state is persisted between steps

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Basic usage:
  python src/scripts/run_analysis.py --plasma <file/folder>

Step-by-step (plasma pipeline):
  # Step 1: fastp + UMI extraction → creates workdir, prints path
  python src/scripts/run_analysis.py --plasma --step preprocess <folder>

  # Subsequent steps: pass the workdir created above
  python src/scripts/run_analysis.py --plasma --step align     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step dedup     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step target    --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step mbias     --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step extract   --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step filter    --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step diff      --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step visualize       --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step report          --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step feature_select  --workdir temp/20240101_120000
  python src/scripts/run_analysis.py --plasma --step classify        --workdir temp/20240101_120000

  # List step definitions and completion status
  python src/scripts/run_analysis.py --list-steps
  python src/scripts/run_analysis.py --list-steps --workdir temp/20240101_120000

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Steps:
  preprocess         4    fastp QC + UMI extraction (plasma: --umi --trim_front)
  align              5    bwameth paired-end alignment
  dedup              6    UMI-aware deduplication (umi_tools dedup --paired)
  target             7    On-target filter (bedtools intersect -abam panel.bed)
  mbias              8    M-bias check (MethylDackel mbias); confirm nOT/nOB trim params
  extract            9    Methylation extraction (MethylDackel → CpG bedGraph)
  filter             10   Depth/repeat filter + candidate-region filter
  diff               11   Differential methylation analysis
  visualize          12   Plots: density / scatter / heatmap
  outlier_detect     13   Sample outlier detection (coverage / PCA z-score / Isolation Forest)
  feature_select     14   Feature selection (LASSO + RF union; outliers excluded; ≥10× coverage QC)
  classify           15   Diagnostic classifiers LOO-CV: LR / SVM / RF / XGBoost / LightGBM
  subtype_classify   16   Subtype multi-class LOO-CV: RF + SVM OvR
  temporal_validate  17   Temporal split validation (plasma only): train <2024-07-01, test >=2024-07-01
  annotate           18   Site annotation (Ensembl + EWAS Catalog+Atlas + CpG island + PubMed)
  external_validate  19   Cross-validation: GSE40279 blood baseline + TCGA-OV tumour signal
  report             20   Generate HTML + JSON report (all optional sections auto-embedded)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import os
import time
import argparse
import logging
import shutil


# Inject conda env bin into PATH before importing other modules
def _inject_conda_path():
    import toml
    _cfg_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "config.toml")
    try:
        _cfg = toml.load(_cfg_path)
        _bin = _cfg.get("tools", {}).get("conda_env_bin", "")
        if _bin and os.path.isdir(_bin) and _bin not in os.environ.get("PATH", ""):
            os.environ["PATH"] = _bin + ":" + os.environ.get("PATH", "")
    except Exception:
        pass

_inject_conda_path()

from utils.logger import setup_logger
from utils.config import load_config, get_output_dirs, create_output_dirs
from scripts.modules import (
    find_files,
    preprocess_fastq,
    align_methylation,
    deduplicate_bam,
    extract_target_regions,
    check_mbias,
    extract_methylation_methyldackel,
    extract_methylation_info,
    merge_methylation_files,
    filter_methylation_data,
    filter_to_candidate_regions_wrapper,
    perform_diff_analysis,
    visualize_results,
    generate_report,
)
from methylation.sample_grouping import resolve_patient_files

config = load_config()

log_dir = config["directories"]["log_dir"]
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/analysis_{time.strftime('%Y%m%d_%H%M%S')}.log"
logger = setup_logger(log_file=log_file, console_level=logging.INFO, file_level=logging.DEBUG)


def check_presets():
    """Verify reference files and required tools before starting the pipeline."""
    logger.info("Checking pipeline prerequisites")

    # 1. Probe panel BED
    probe_panel = config.get("data", {}).get("probe_panel_bed")
    if probe_panel:
        if not os.path.exists(probe_panel):
            logger.error(f"Probe panel BED not found: {probe_panel}")
            logger.warning("On-target filter step will be skipped")
        else:
            logger.info(f"Probe panel: {probe_panel}")
    else:
        logger.warning("No probe_panel_bed configured — on-target filter will be skipped")

    # 2. Reference genome
    reference_genome = config.get("extraction", {}).get("reference_genome")
    if not reference_genome or not os.path.exists(reference_genome):
        logger.error(f"Reference genome not found: {reference_genome}")
        return False
    logger.info(f"Reference genome: {reference_genome}")

    # 3. Tool checks (aligner-dependent)
    aligner = config.get("plasma", {}).get("aligner", "bismark")
    if aligner == "bwameth":
        required_tools = ["fastp", "bwameth", "samtools", "bedtools", "MethylDackel"]
    else:
        bismark_index = os.path.dirname(reference_genome)
        bismark_index_file = os.path.join(
            bismark_index, "Bisulfite_Genome", "CT_conversion", "BS_CT.1.bt2"
        )
        if not os.path.exists(bismark_index_file):
            logger.warning(f"Bismark index not found: {bismark_index_file}")
            logger.warning("Run first: bismark_genome_preparation --bowtie2 " + bismark_index)
        else:
            logger.info(f"Bismark index: {bismark_index}/Bisulfite_Genome/")
        required_tools = ["fastp", "bismark", "deduplicate_bismark", "samtools", "bedtools", "MethylDackel"]

    # umi_tools may live outside the conda env — detect separately
    umi_tools_bin = shutil.which("umi_tools") or config.get("tools", {}).get("umi_tools_fallback", "")
    if umi_tools_bin and os.path.isfile(umi_tools_bin):
        logger.info(f"umi_tools: {umi_tools_bin}")
        umi_dir = os.path.dirname(umi_tools_bin)
        if umi_dir and umi_dir not in os.environ.get("PATH", ""):
            os.environ["PATH"] = umi_dir + ":" + os.environ.get("PATH", "")
    else:
        logger.warning("umi_tools not found — plasma dedup will fall back to deduplicate_bismark (not recommended)")

    missing_tools = []
    for tool in required_tools:
        if shutil.which(tool):
            logger.info(f"Tool {tool}: {shutil.which(tool)}")
        else:
            logger.warning(f"Tool {tool} not installed or not in PATH")
            missing_tools.append(tool)

    if missing_tools:
        logger.warning(f"Missing tools: {', '.join(missing_tools)}")

    logger.info("Prerequisite check complete")
    return True


# ──────────────────────────────────────────────
# Step definitions and state management
# ──────────────────────────────────────────────

# (step_key, display_num, description, depends_on)
PIPELINE_STEPS = [
    ("data_check",         "1",   "Validate config: sample Excel, ref genome, panel BED, FASTQs",  "—"),
    ("ref_index",          "2",   "Build/verify bwameth hg38 c2t index (~17 GB, one-time setup)",  "—"),
    ("panel_derive",       "3",   "Derive panel BED (reuse dedup BAMs via --workdir, or --fq)",   "ref_index/dedup"),
    ("preprocess",         "4",   "fastp QC + UMI extraction (R2 5' 6N)",                         "—"),
    ("align",              "5",   "bwameth alignment → sorted BAM (hg38)",                        "preprocess"),
    ("dedup",              "6",   "UMI-aware dedup (umi_tools dedup --paired)",                   "align"),
    ("target",             "7",   "On-target filter (bedtools intersect panel.bed)",               "dedup"),
    ("mbias",              "8",   "M-bias check (MethylDackel mbias); review nOT/nOB SVGs",       "target"),
    ("extract",            "9",   "Methylation extraction (MethylDackel → CpG bedGraph)",         "mbias"),
    ("filter",             "10",  "Depth+repeat filter; candidate-region filter (88 regions)",     "extract"),
    ("diff",               "11",  "Differential methylation: merged bedGraphs, |Δβ| ≥ 0.2",      "filter"),
    ("visualize",          "12",  "Plots: density / scatter / heatmap (candidate regions)",        "diff"),
    ("outlier_detect",     "13",  "Sample QC: coverage fraction / PCA z / Isolation Forest",       "diff"),
    ("feature_select",     "14",  "Feature selection: LASSO(L1) ∪ RF; outliers excluded; ≥10× QC", "outlier_detect"),
    ("classify",           "15",  "Classifiers LOO-CV: LR / SVM / RF / XGBoost / LightGBM",       "feature_select"),
    ("subtype_classify",   "16",  "Subtype LOO-CV (malignant+borderline): RF + SVM OvR",          "classify"),
    ("temporal_validate",  "17",  "Temporal split validation (plasma only): train/test by date",  "classify"),
    ("annotate",           "18",  "Site annotation: gene / CpG island / EWAS Catalog+Atlas",      "feature_select"),
    ("external_validate",  "19",  "Cross-validation: GSE40279 blood baseline + TCGA-OV signal",   "feature_select"),
    ("report",             "20",  "HTML + JSON report (annotation + ext-validate auto-embedded)",  "diff"),
]
STEP_KEYS = [s[0] for s in PIPELINE_STEPS]
_STATE_FILE = "pipeline_state.json"


def _state_path(workdir):
    return os.path.join(workdir, _STATE_FILE)


def _load_state(workdir):
    """Load step outputs from workdir/pipeline_state.json."""
    import json
    path = _state_path(workdir)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def _save_state(workdir, state):
    """Merge state into workdir/pipeline_state.json (non-destructive update)."""
    import json
    path = _state_path(workdir)
    existing = _load_state(workdir)
    existing.update(state)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(existing, f, ensure_ascii=False, indent=2)


def _run_ref_index(cfg):
    """Check bwameth c2t index; build it if missing."""
    ref = cfg.get("extraction", {}).get("reference_genome", "")
    if not ref:
        logger.error("reference_genome not set in config — cannot check/build index")
        return False
    c2t = ref + ".bwameth.c2t"
    if os.path.exists(c2t):
        logger.info(f"✓ bwameth index found: {c2t}")
        return True
    logger.warning(f"bwameth c2t index not found: {c2t}")
    logger.info(f"Building index: bwameth --index {ref}  (this may take 30–60 min)")
    import subprocess
    ret = subprocess.run(["bwameth", "--index", ref])
    if ret.returncode == 0:
        logger.info("✓ bwameth index built successfully")
        return True
    logger.error(f"bwameth --index failed (exit {ret.returncode})")
    return False


def _print_step_list(workdir=None):
    """Print all pipeline steps with completion status, description, and dependencies."""
    state     = _load_state(workdir) if workdir else {}
    completed = state.get("completed_steps", [])

    W_STEP, W_NUM, W_DESC, W_REQ = 22, 5, 60, 18

    top = f"┌{'─'*(W_STEP+2)}┬{'─'*(W_NUM+2)}┬{'─'*(W_DESC+2)}┬{'─'*(W_REQ+2)}┐"
    mid = f"├{'─'*(W_STEP+2)}┼{'─'*(W_NUM+2)}┼{'─'*(W_DESC+2)}┼{'─'*(W_REQ+2)}┤"
    bot = f"└{'─'*(W_STEP+2)}┴{'─'*(W_NUM+2)}┴{'─'*(W_DESC+2)}┴{'─'*(W_REQ+2)}┘"

    def _cell(text, width, align="<"):
        text = str(text)
        if len(text) > width:
            text = text[:width - 1] + "…"
        return f"{text:{align}{width}}"

    def _row(key, num, desc, req):
        mark = "✓" if key in completed else " "
        return (
            f"│ {_cell(f'[{mark}] {key}', W_STEP)} "
            f"│ {_cell(num, W_NUM, '^')} "
            f"│ {_cell(desc, W_DESC)} "
            f"│ {_cell(req, W_REQ)} │"
        )

    print()
    print(top)
    print(
        f"│ {_cell('Step  (--step arg)', W_STEP)} "
        f"│ {_cell('#', W_NUM, '^')} "
        f"│ {_cell('Description', W_DESC)} "
        f"│ {_cell('Requires', W_REQ)} │"
    )
    print(mid)
    for entry in PIPELINE_STEPS:
        print(_row(*entry))
    print(bot)
    if workdir:
        n_done = sum(1 for e in PIPELINE_STEPS if e[0] in completed)
        print(f"\n  Workdir  : {workdir}")
        print(f"  Progress : {n_done} / {len(PIPELINE_STEPS)} steps completed")
        print(f"  State    : {_state_path(workdir)}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="AIMethFlow — methylation analysis pipeline (full or step-by-step)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "files", nargs="*",
        help="FastQ files or folders (may be omitted in step mode when workdir already exists)",
    )

    # Step mode
    parser.add_argument(
        "--step", metavar="STEP",
        choices=STEP_KEYS,
        help=f"Run a single step (choices: {', '.join(STEP_KEYS)})",
    )
    parser.add_argument(
        "--workdir", metavar="DIR",
        help="Resume step-by-step execution from an existing workdir (reads pipeline_state.json)",
    )
    parser.add_argument(
        "--list-steps", action="store_true",
        help="Print step table with completion status and exit",
    )

    # Input / data options
    parser.add_argument(
        "--patient-csv", metavar="CSV",
        help="CSV file with one patient ID per line; R1/R2 paths are resolved from the sample info Excel",
    )
    parser.add_argument(
        "--data-dir", metavar="DIR",
        help="FastQ data directory (used with --patient-csv; defaults to config [data] plasma_dir)",
    )
    parser.add_argument(
        "--temp-dir", metavar="DIR",
        help="Override temp directory from config (e.g. /mnt/g/temp)",
    )
    parser.add_argument(
        "--mode", choices=["plasma", "tissue"], default="plasma",
        help="Pipeline mode: 'plasma' (cfDNA) or 'tissue' (biopsy). Default: plasma",
    )
    parser.add_argument(
        "--group-method", choices=["v1", "v2"], default="v2",
        help="Sample grouping method: v1 (index-based) or v2 (Excel-based)",
    )

    # Skip flags for full-pipeline mode (step mode ignores these)
    parser.add_argument("--skip-quality",           action="store_true", default=True)
    parser.add_argument("--skip-preprocessing",     action="store_true", default=False)
    parser.add_argument("--skip-alignment",         action="store_true", default=False)
    parser.add_argument("--skip-deduplication",     action="store_true", default=False)
    parser.add_argument("--skip-target-extraction", action="store_true", default=False)
    parser.add_argument("--skip-extraction",        action="store_true", default=False)
    parser.add_argument("--skip-filtering",         action="store_true", default=False)
    parser.add_argument("--skip-candidate-regions", action="store_true", default=False)
    parser.add_argument("--skip-diff",              action="store_true", default=False)
    parser.add_argument("--skip-visualization",     action="store_true", default=False)
    parser.add_argument("--skip-report",            action="store_true", default=False)
    parser.add_argument("--clean-temp",             action="store_true", default=False)
    parser.add_argument("--config",  type=str, help="Path to a custom config file")
    parser.add_argument("--verbose", action="store_true", default=False)

    args = parser.parse_args()
    mode = args.mode  # "plasma" | "tissue"

    if args.list_steps:
        _print_step_list(args.workdir)
        return

    global config
    if args.config and os.path.exists(args.config):
        from utils.config import load_config as _load_cfg
        config = _load_cfg(args.config)

    if args.temp_dir:
        config["directories"]["temp_dir"] = args.temp_dir
        os.makedirs(args.temp_dir, exist_ok=True)
        logger.info(f"Temp dir overridden: {args.temp_dir}")

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        for h in logger.handlers:
            if getattr(h, "name", None) == "console":
                h.setLevel(logging.DEBUG)

    if not check_presets():
        logger.error("Prerequisite check failed — aborting")
        return

    # Early-exit setup steps (no workdir / FASTQ needed)
    if args.step == "data_check":
        logger.info("✓  Step [data_check] complete — prerequisites verified above")
        return
    if args.step == "ref_index":
        _run_ref_index(config)
        return

    # Resolve working directory
    if args.workdir:
        workdir = args.workdir
        os.makedirs(workdir, exist_ok=True)
        timestamp  = os.path.basename(os.path.normpath(workdir))
        output_dirs = get_output_dirs(config, timestamp=timestamp)
        result_dir  = output_dirs["result"]
        os.makedirs(result_dir, exist_ok=True)
    else:
        output_dirs = get_output_dirs(config)
        create_output_dirs(output_dirs)
        workdir    = output_dirs["temp"]
        result_dir = output_dirs["result"]

    logger.info(f"Working directory: {workdir}")

    state = _load_state(workdir)

    # Resolve input files
    # Only preprocess / align / panel_derive (or full-pipeline start) strictly need FastQ files
    _STEPS_NEED_FASTQ = {"preprocess", "align", "panel_derive"}

    if getattr(args, "patient_csv", None):
        data_dir = getattr(args, "data_dir", None)
        fastq_files = resolve_patient_files(args.patient_csv, data_dir=data_dir)
        bam_files_input = []
        methylation_files_input = []
        if not fastq_files:
            logger.error(f"--patient-csv resolved no files — check {args.patient_csv} and data dir")
            return
        _save_state(workdir, {"fastq_files": fastq_files})
    elif args.files:
        found_files = find_files(args.files)
        fastq_files = found_files["fastq"]
        bam_files_input = found_files["bam"]
        methylation_files_input = found_files["methylation"]
        if fastq_files:
            state["fastq_files"] = fastq_files
        _save_state(workdir, {"fastq_files": fastq_files})
    else:
        fastq_files = state.get("fastq_files", [])
        bam_files_input = []
        methylation_files_input = []
        if not fastq_files and (not single_step or single_step in _STEPS_NEED_FASTQ):
            logger.error("No input files provided and no fastq_files recorded in workdir state")
            return

    if fastq_files:
        logger.info(f"Input FastQ files ({len(fastq_files)}): {fastq_files}")

    # Restore per-step outputs from state
    processed_files       = state.get("processed_files", fastq_files)
    bam_files             = state.get("bam_files", bam_files_input)
    deduplicated_files    = state.get("deduplicated_files", bam_files)
    target_files          = state.get("target_files", deduplicated_files)
    sample1_methylation   = state.get("sample1_methylation")
    sample2_methylation   = state.get("sample2_methylation")
    filtered_sample1      = state.get("filtered_sample1")
    filtered_sample2      = state.get("filtered_sample2")
    candidate_sample1     = state.get("candidate_sample1")
    candidate_sample2     = state.get("candidate_sample2")
    benign_samples        = state.get("benign_samples", [])
    malignant_samples     = state.get("malignant_samples", [])

    single_step = args.step  # None → full pipeline

    def _should_run(step_key, skip_flag):
        if single_step:
            return step_key == single_step
        return not skip_flag

    # =========================================================
    # Step 3: panel_derive  (setup step — skipped in full pipeline)
    # =========================================================
    if _should_run("panel_derive", True) and fastq_files:
        logger.info("━━ Step 3: Derive panel BED from FASTQ alignments ━━")
        from scripts.derive_panel_bed import derive_panel_bed as _derive_panel
        out_bed = config.get("data", {}).get("probe_panel_bed") or "data/panel_B10308.bed"
        ok = _derive_panel(fastq_files, out_bed, config=config, workdir=workdir)
        if ok:
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["panel_derive"],
                                  "panel_bed": out_bed})
            logger.info(f"✓  Step [panel_derive] complete — {out_bed}")
        else:
            logger.error("panel_derive failed")
        if single_step:
            return

    # =========================================================
    # Step 4: preprocess
    # =========================================================
    if _should_run("preprocess", args.skip_preprocessing) and fastq_files:
        logger.info("━━ Step 4: fastp QC + UMI extraction ━━")
        processed_files = preprocess_fastq(fastq_files, workdir, plasma=(mode == "plasma"), config=config)
        if processed_files:
            _save_state(workdir, {"processed_files": processed_files,
                                  "completed_steps": state.get("completed_steps", []) + ["preprocess"]})
            state["processed_files"] = processed_files
            logger.info(f"  → Preprocessed: {processed_files}")
        else:
            logger.error("Preprocessing failed — processed_files is empty")
            if single_step:
                return
    elif not single_step:
        logger.info("4. Skipping preprocessing")

    # =========================================================
    # Step 5: align
    # =========================================================
    if _should_run("align", args.skip_alignment) and processed_files and not bam_files:
        logger.info("━━ Step 5: bwameth paired-end alignment ━━")
        bam_files = align_methylation(processed_files, workdir, config, plasma=(mode == "plasma"))
        if bam_files:
            _save_state(workdir, {"bam_files": bam_files,
                                  "completed_steps": state.get("completed_steps", []) + ["align"]})
            state["bam_files"] = bam_files
            logger.info(f"  → BAM files: {bam_files}")
        else:
            logger.error("Alignment failed — bam_files is empty")
            if single_step:
                return
    elif single_step and bam_files:
        logger.info(f"✓ align already done — {len(bam_files)} BAM(s) in state, skipping")
    elif not single_step:
        logger.info("5. Skipping alignment" if not bam_files else "5. Using existing BAM files")

    # =========================================================
    # Step 6: dedup
    # =========================================================
    if _should_run("dedup", args.skip_deduplication) and bam_files:
        logger.info("━━ Step 6: UMI-aware deduplication ━━")
        deduplicated_files = deduplicate_bam(bam_files, workdir, plasma=(mode == "plasma"), config=config)
        if deduplicated_files:
            _save_state(workdir, {"deduplicated_files": deduplicated_files,
                                  "completed_steps": state.get("completed_steps", []) + ["dedup"]})
            state["deduplicated_files"] = deduplicated_files
            logger.info(f"  → Deduplicated: {deduplicated_files}")
        else:
            logger.error("Deduplication failed — deduplicated_files is empty")
            if single_step:
                return
    elif not single_step:
        logger.info("6. Skipping deduplication")

    # =========================================================
    # Step 7: target
    # =========================================================
    if _should_run("target", args.skip_target_extraction) and deduplicated_files:
        logger.info("━━ Step 7: On-target filter ━━")
        target_files = extract_target_regions(deduplicated_files, workdir, config)
        if target_files:
            _save_state(workdir, {"target_files": target_files,
                                  "completed_steps": state.get("completed_steps", []) + ["target"]})
            state["target_files"] = target_files
            logger.info(f"  → On-target BAMs: {target_files}")
        else:
            logger.warning("On-target filter produced no output — subsequent steps will use deduplicated BAMs")
    elif not single_step:
        logger.info("7. Skipping on-target filter")

    # =========================================================
    # Step 8: mbias
    # =========================================================
    if _should_run("mbias", False) and mode == "plasma":
        mbias_input = target_files or deduplicated_files
        if mbias_input:
            reference_genome = config.get("extraction", {}).get("reference_genome")
            if reference_genome and os.path.exists(reference_genome):
                logger.info("━━ Step 8: M-bias check ━━")
                check_mbias(mbias_input, workdir, reference_genome, config=config)
                _save_state(workdir, {"mbias_done": True,
                                      "completed_steps": state.get("completed_steps", []) + ["mbias"]})
                logger.info("  → M-bias SVGs saved to workdir — review R2 3' bias before extract")
            else:
                logger.warning("7.5. Skipping M-bias check (reference genome not found)")
        else:
            logger.warning("7.5. Skipping M-bias check (no BAM available)")

    # =========================================================
    # Step 9: extract
    # =========================================================
    if _should_run("extract", args.skip_extraction) and not methylation_files_input:
        bam_input = target_files or deduplicated_files
        if bam_input:
            logger.info("━━ Step 9: Methylation extraction ━━")
            if mode == "plasma":
                aligner = config.get("plasma", {}).get("aligner", "bismark")
            else:
                aligner = config.get("tissue", {}).get("aligner", "bismark")
            logger.info(f"  [Mode: {mode}  Aligner: {aligner}]")

            if aligner == "bwameth":
                logger.info("  [Extractor: MethylDackel (aligner=bwameth)]")
                meth_files = extract_methylation_methyldackel(bam_input, workdir, config)
                if meth_files:
                    from methylation.sample_grouping import group_methylation_files_v2
                    _mode_cfg = config.get(mode, {})
                    _excel_sheet = _mode_cfg.get("excel_sheet", mode)
                    logger.info(f"  [Excel sheet: {_excel_sheet}]")
                    benign_files, malignant_files = group_methylation_files_v2(
                        meth_files, category=_excel_sheet
                    )
                    logger.info(f"  → Grouped: {len(benign_files)} benign, {len(malignant_files)} malignant")

                    s1_path = os.path.join(workdir, "sample1_methylation.bedGraph")
                    s2_path = os.path.join(workdir, "sample2_methylation.bedGraph")

                    def _sample_id(path):
                        return os.path.basename(path).split("_CpG")[0].split(".")[0].split("_")[0]

                    if benign_files:
                        merge_methylation_files(benign_files, s1_path)
                        sample1_methylation = s1_path
                        benign_samples = [_sample_id(f) for f in benign_files]
                    else:
                        sample1_methylation = meth_files[0]
                        benign_samples = [_sample_id(meth_files[0])]
                        logger.warning("  → No benign files found — using first file as sample1")

                    if malignant_files:
                        merge_methylation_files(malignant_files, s2_path)
                        sample2_methylation = s2_path
                        malignant_samples = [_sample_id(f) for f in malignant_files]
                    else:
                        sample2_methylation = meth_files[1] if len(meth_files) > 1 else meth_files[0]
                        malignant_samples = [_sample_id(sample2_methylation)]
                        logger.warning("  → No malignant files found — using second file as sample2")
            else:
                logger.info("  [Extractor: bismark_methylation_extractor (aligner=bismark)]")
                precomputed_bams = bam_input if mode == "plasma" else None
                sample1_methylation, sample2_methylation = extract_methylation_info(
                    workdir, fastq_files,
                    group_method=args.group_method,
                    config=config,
                    bam_files=precomputed_bams,
                    category=mode,
                )

            _save_state(workdir, {
                "sample1_methylation": sample1_methylation,
                "sample2_methylation": sample2_methylation,
                "benign_samples": benign_samples,
                "malignant_samples": malignant_samples,
                "completed_steps": state.get("completed_steps", []) + ["extract"],
            })
            state["sample1_methylation"] = sample1_methylation
            state["sample2_methylation"] = sample2_methylation
            state["benign_samples"] = benign_samples
            state["malignant_samples"] = malignant_samples
            logger.info(f"  → Methylation files: {sample1_methylation}, {sample2_methylation}")
        else:
            logger.error("No BAM files available for methylation extraction")
            if single_step:
                return
    elif methylation_files_input and not single_step:
        logger.info("8. Using provided methylation files")
        sample1_methylation = methylation_files_input[0] if methylation_files_input else None
        sample2_methylation = methylation_files_input[1] if len(methylation_files_input) > 1 else sample1_methylation
    elif not single_step:
        logger.info("8. Skipping methylation extraction")

    # =========================================================
    # Step 10: filter
    # =========================================================
    if _should_run("filter", args.skip_filtering) and (sample1_methylation or sample2_methylation):
        logger.info("━━ Step 10: Methylation filtering + candidate regions ━━")
        if sample1_methylation:
            filtered_sample1 = filter_methylation_data(sample1_methylation, workdir, config)
        if sample2_methylation:
            filtered_sample2 = filter_methylation_data(sample2_methylation, workdir, config)
        if not args.skip_candidate_regions:
            if filtered_sample1:
                candidate_sample1 = filter_to_candidate_regions_wrapper(filtered_sample1, workdir, config)
            if filtered_sample2:
                candidate_sample2 = filter_to_candidate_regions_wrapper(filtered_sample2, workdir, config)
        _save_state(workdir, {
            "filtered_sample1": filtered_sample1, "filtered_sample2": filtered_sample2,
            "candidate_sample1": candidate_sample1, "candidate_sample2": candidate_sample2,
            "completed_steps": state.get("completed_steps", []) + ["filter"],
        })
        state.update({
            "filtered_sample1": filtered_sample1, "filtered_sample2": filtered_sample2,
            "candidate_sample1": candidate_sample1, "candidate_sample2": candidate_sample2,
        })
        logger.info(f"  → Candidate region files: {candidate_sample1}, {candidate_sample2}")
    elif not single_step:
        logger.info("9. Skipping methylation filtering")

    # =========================================================
    # Step 11: diff
    # =========================================================
    if _should_run("diff", args.skip_diff) and candidate_sample1 and candidate_sample2:
        logger.info("━━ Step 11: Differential methylation analysis ━━")
        perform_diff_analysis(candidate_sample1, candidate_sample2, result_dir)
        _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["diff"]})
    elif not single_step:
        logger.info("10. Skipping differential analysis")

    # =========================================================
    # Step 12: visualize
    # =========================================================
    if _should_run("visualize", args.skip_visualization):
        vis_input = candidate_sample1 or sample1_methylation
        if vis_input and os.path.exists(vis_input):
            logger.info("━━ Step 12: Visualization ━━")
            visualize_results(vis_input, result_dir)
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["visualize"]})
        else:
            logger.warning("11. Skipping visualization (no methylation file available)")
    elif not single_step:
        logger.info("11. Skipping visualization")

    # =========================================================
    # Step 13: outlier_detect  (runs BEFORE feature_select so that
    #           outlier samples are excluded from LASSO/RF selection)
    # =========================================================
    if _should_run("outlier_detect", False):
        if benign_samples and malignant_samples:
            logger.info("━━ Step 13: Sample outlier detection ━━")
            import numpy as np
            from scripts.modules.feature_utils import build_feature_matrix, detect_sample_outliers
            fm_path = os.path.join(result_dir, "feature_matrix.npz")
            if os.path.exists(fm_path):
                data    = np.load(fm_path, allow_pickle=True)
                X_od    = data["X"]
                y_od    = data["y"]
                sids_od = list(data["sample_ids"])
                logger.info(f"  Loaded cached feature matrix: {X_od.shape[0]} × {X_od.shape[1]}")
            else:
                logger.info("  Building feature matrix from bedGraphs...")
                X_od, y_od, sids_od, feat_od, D_od, ri88_od = build_feature_matrix(
                    workdir, result_dir, benign_samples, malignant_samples, config
                )
                # Save so feature_select (Step 14) can load from cache
                np.savez(fm_path,
                         X=X_od, y=y_od,
                         sample_ids=np.array(sids_od, dtype=object),
                         feature_names=np.array(feat_od, dtype=object),
                         D=D_od,
                         region_in_88=ri88_od)
                logger.info(f"  Feature matrix cached: {fm_path} "
                            f"({X_od.shape[0]} samples × {X_od.shape[1]} features)")
            outlier_df, suspect_ids = detect_sample_outliers(X_od, y_od, sids_od, config)
            if not outlier_df.empty:
                outlier_path = os.path.join(result_dir, "outlier_report.csv")
                outlier_df.to_csv(outlier_path, index=False)
                logger.info(f"  → Outlier report saved: {outlier_path}")
            if suspect_ids:
                logger.warning(f"  ⚠ Suspected outliers (≥2 flags): {suspect_ids}")
                logger.warning(
                    "  ⚠ These samples will be EXCLUDED from feature selection (Step 14)."
                )
                logger.warning(
                    "  ⚠ Review outlier_report.csv; add confirmed IDs to "
                    "config.toml [ml] exclude_samples to also exclude from classifiers."
                )
            else:
                logger.info("  No suspect samples detected")
            _save_state(workdir, {
                "suspect_ids": suspect_ids,
                "completed_steps": state.get("completed_steps", []) + ["outlier_detect"],
            })
            state["suspect_ids"] = suspect_ids
        else:
            logger.error("outlier_detect requires benign_samples and malignant_samples — "
                         "complete the extract step first")
            if single_step:
                return

    # =========================================================
    # Step 14: feature_select  (reads suspect_ids from Step 13 state
    #           to exclude outlier samples from LASSO/RF selection)
    # =========================================================
    if _should_run("feature_select", False):
        if benign_samples and malignant_samples:
            logger.info("━━ Step 14: Feature selection (LASSO + RF) ━━")
            from scripts.modules.feature_utils import select_features
            suspect_ids = state.get("suspect_ids") or []
            if suspect_ids:
                logger.info(
                    f"  Excluding {len(suspect_ids)} outlier sample(s) from feature selection: "
                    f"{suspect_ids}"
                )
            panel_features, _ = select_features(
                workdir, result_dir, config, benign_samples, malignant_samples,
                extra_exclude_ids=suspect_ids or None,
            )
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["feature_select"]})
            logger.info(f"  → {len(panel_features)} features selected — saved to {result_dir}")
        else:
            logger.error("feature_select requires benign_samples and malignant_samples — "
                         "complete the extract step first (or ensure pipeline_state.json has grouping info)")
            if single_step:
                return

    # =========================================================
    # Step 15: classify
    # =========================================================
    if _should_run("classify", False):
        if benign_samples and malignant_samples:
            logger.info("━━ Step 15: Diagnostic classifiers LOO-CV ━━")
            from scripts.modules.classifier_utils import run_classifiers
            _, _ = run_classifiers(
                workdir, result_dir, config, benign_samples, malignant_samples
            )
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["classify"]})
            logger.info(f"  → Classifier results saved to {result_dir}")
        else:
            logger.error("classify requires benign_samples and malignant_samples — "
                         "complete the extract step first")
            if single_step:
                return

    # =========================================================
    # Step 16: subtype_classify
    # =========================================================
    if _should_run("subtype_classify", False):
        logger.info("━━ Step 16: Subtype classifier LOO-CV ━━")
        from scripts.modules.subtype_utils import run_subtype_classifier
        sub_results, _ = run_subtype_classifier(workdir, result_dir, config)
        _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["subtype_classify"]})
        if sub_results:
            logger.info(f"  → Subtype classification results saved to {result_dir}")

    # =========================================================
    # Step 17: temporal_validate (plasma only)
    # =========================================================
    if _should_run("temporal_validate", False):
        logger.info("━━ Step 17: Temporal split validation ━━")
        if mode != "plasma":
            logger.warning("temporal_validate is only applicable in plasma mode — skipping")
        elif not (benign_samples and malignant_samples):
            logger.error("temporal_validate requires benign_samples and malignant_samples — "
                         "complete the classify step first")
            if single_step:
                return
        else:
            from scripts.modules.classifier_utils import run_temporal_validation
            from scripts.modules.subtype_utils import run_subtype_temporal_validation
            run_temporal_validation(workdir, result_dir, config, benign_samples, malignant_samples)
            # run_subtype_temporal_validation(result_dir, config)
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["temporal_validate"]})
            logger.info(f"  → Temporal validation results saved to {result_dir}")

    # =========================================================
    # Step 18: annotate
    # =========================================================
    if _should_run("annotate", False):
        logger.info("━━ Step 18: Site annotation (Ensembl + EWAS Catalog) ━━")
        from scripts.modules.annotation_utils import run_annotation
        ann_df, ann_summary = run_annotation(result_dir, config)
        if not ann_df.empty:
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["annotate"]})
            logger.info(f"  → Annotation saved to {result_dir}/panel_annotation.csv")
            logger.info(
                f"  → EWAS hits: {ann_summary.get('n_with_ewas', 0)}/{ann_summary.get('n_sites', 0)} sites "
                f"(window ±{ann_summary.get('ewas_window_bp', 500)} bp)"
            )
        else:
            logger.error("  Annotation produced no results — check that panel_cpg.csv exists")

    # =========================================================
    # Step 19: external_validate
    # =========================================================
    if _should_run("external_validate", False):
        logger.info("━━ Step 19: External dataset cross-validation ━━")
        from scripts.external_validate import run_external_validation
        data_dir_ext = config.get("directories", {}).get("data_dir", "")
        external_dir = os.path.join(data_dir_ext, "external")
        if not os.path.isdir(external_dir):
            logger.error(f"  external_dir not found: {external_dir} — "
                         "place reference data under data/external/ and retry")
            if single_step:
                return
        else:
            run_external_validation(result_dir, data_dir_ext, external_dir)
            _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["external_validate"]})
            logger.info(f"  → External validation saved: {os.path.join(result_dir, 'external_validation.csv')}")
            logger.info("  → Run --step report to embed Section 13 into the HTML report")

    # =========================================================
    # Step 20: report
    # =========================================================
    if _should_run("report", args.skip_report):
        logger.info("━━ Step 20: Report generation ━━")
        generate_report(workdir, result_dir,
                        sample1_methylation, sample2_methylation,
                        candidate_sample1, candidate_sample2,
                        benign_samples=benign_samples,
                        malignant_samples=malignant_samples)
        _save_state(workdir, {"completed_steps": state.get("completed_steps", []) + ["report"]})

    # Completion message
    if single_step:
        print(f"\n✓  Step [{single_step}] complete")
        print(f"   Workdir (pass to next step via --workdir): {workdir}")
        _print_step_list(workdir)
    else:
        print(f"\n✓  Full pipeline complete — results in: {result_dir}")

    if args.clean_temp and not single_step:
        try:
            shutil.rmtree(workdir)
            logger.info(f"Temp directory cleaned: {workdir}")
        except Exception as e:
            logger.error(f"Failed to clean temp directory: {e}")


if __name__ == "__main__":
    main()
