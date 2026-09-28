# Quick Start

This guide assumes you have already completed [installation.md](installation.md) (or have the Docker container running — see [../DOCKER.md](../DOCKER.md)). It walks through the Web UI, the one-line full-pipeline command, the step-by-step execution flow, and a complete shell-script template for batch jobs.

## Web UI (recommended)

Visual interface for configuring, running individual steps, and browsing results and reports.

**Docker mode** (auto-starts after `docker compose up -d`):

```text
http://<server-ip>:9999
```

**Local development mode**:

```bash
# Run from the project root; METH_WORKDIR tells the Web UI where to find config.toml / temp / result
export METH_WORKDIR=$(pwd)              # Linux / macOS / WSL
# $env:METH_WORKDIR = (Get-Location).Path  # PowerShell

conda activate meth
python src/scripts/webui.py
# Open http://localhost:9999 in your browser
```

> The Web UI reads and writes `$METH_WORKDIR/config.toml` — in local dev mode this is the project root `config.toml`.
> Inside the Docker container it uses `config.docker.toml`; the two paths are independent.

## Full pipeline (one command)

```bash
run-analysis --mode plasma /mnt/f/meth/OC_plasma/
```

## Step-by-step execution

### Phase 1 — one-time setup (per environment / per panel design)

```bash
# Step 1: validate config prerequisites (ref genome, panel BED, sample Excel)
run-analysis --step data_check

# Step 2: build bwameth c2t index for hg38 (~17 GB, 30–60 min, run once)
run-analysis --step ref_index

# Step 3: derive panel BED
# Recommended: run 2–3 pilot samples through align + dedup first, then reuse their BAMs.
# This avoids redundant alignment — the dedup BAMs are already computed.

PILOT_WORKDIR=temp/pilot
run-analysis --mode plasma --step preprocess \
    /mnt/g/OC_tissue/P10001.R1.fq.gz \
    /mnt/g/OC_tissue/P10002.R1.fq.gz \
    --workdir $PILOT_WORKDIR
run-analysis --mode plasma --step align --workdir $PILOT_WORKDIR
run-analysis --mode plasma --step dedup --workdir $PILOT_WORKDIR

python src/scripts/derive_panel_bed.py \
    --workdir $PILOT_WORKDIR \
    --out data/panel_B10308.bed

# Alternative: derive directly from FASTQ (re-aligns, slower)
# python src/scripts/derive_panel_bed.py \
#     --fq /mnt/g/OC_tissue/P10001.R1.fq.gz /mnt/g/OC_tissue/P10001.R2.fq.gz \
#     --out data/panel_B10308.bed
```

### Phase 2 — per-run analysis (all samples)

```bash
WORKDIR=temp/20260406

# Step 4: preprocessing (creates workdir; subsequent steps reuse it)
run-analysis --mode plasma --step preprocess /mnt/f/meth/OC_plasma/  --workdir $WORKDIR
# WORKDIR=$(ls -dt temp/20[0-9]* | head -1)   # capture generated workdir

# Steps 5–12: alignment / dedup / extraction / analysis
run-analysis --mode plasma --step align          --workdir $WORKDIR
run-analysis --mode plasma --step dedup          --workdir $WORKDIR
run-analysis --mode plasma --step target         --workdir $WORKDIR
run-analysis --mode plasma --step mbias          --workdir $WORKDIR
# ── mbias writes mbias_auto_trim.json; review SVGs if you want to override ──
# ── set nOT/nOB = "auto" in config.toml to apply auto-detected values      ──
run-analysis --mode plasma --step extract        --workdir $WORKDIR
run-analysis --mode plasma --step filter         --workdir $WORKDIR
run-analysis --mode plasma --step diff           --workdir $WORKDIR
run-analysis --mode plasma --step visualize      --workdir $WORKDIR

# Step 13: outlier detection → review result/outlier_report.csv
run-analysis --mode plasma --step outlier_detect --workdir $WORKDIR
#   Add confirmed outlier IDs to config.toml [ml] exclude_samples = [...]

# Step 14: feature analysis (LASSO + RF selection; outliers auto-excluded)
run-analysis --mode plasma --step feature_analysis --workdir $WORKDIR

# Step 14a (optional): drop GSE40279 blood-confounded sites before classifying
run-analysis --mode plasma --step blood_filter     --workdir $WORKDIR

# Steps 15–16: classification
run-analysis --mode plasma --step classify         --workdir $WORKDIR
# Step 14b (optional): feature selection ablation — wilcoxon / lasso / rf / union AUC comparison
run-analysis --mode plasma --step fs_benchmark     --workdir $WORKDIR
# Step 15b (optional): leakage demo — naive vs in-fold selection ΔAUC, paired DeLong
run-analysis --mode plasma --step leakage_demo     --workdir $WORKDIR
# Step 15c (optional): train production models on the full cohort → result/$WORKDIR/model/
run-analysis --mode plasma --step final_model      --workdir $WORKDIR
run-analysis --mode plasma --step subtype_classify --workdir $WORKDIR

# Step 17: temporal split validation (plasma only; requires data/plasma_exp.xlsx)
run-analysis --mode plasma --step temporal_validate --workdir $WORKDIR

# Step 18: annotate selected panel sites (network access required on first run; cached after)
run-analysis --mode plasma --step annotate          --workdir $WORKDIR
# Step 19 (optional): cross-validate against public datasets (requires data/external/)
run-analysis --mode plasma --step external_validate --workdir $WORKDIR
# Step 19b (optional): re-run Step 15 on independent cohorts listed in [external_cohort] dirs
run-analysis --mode plasma --step external_cohort   --workdir $WORKDIR
# Step 20 (last): generate report — auto-embeds the Step 18 / 19 / 15b / 19b sections
run-analysis --mode plasma --step report            --workdir $WORKDIR

# Show all steps and current progress
run-analysis --list-steps --workdir $WORKDIR
zip -r 20260406_report.zip result/20260406/20260406_report/
```

## Shell script (background execution)

Suitable for batch jobs run unattended with logging.

```bash
#!/bin/bash
# Usage: bash run_oc_plasma.sh 2>&1 | tee logs/oc_plasma_run.log
set -e
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
PY=/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin/python

# ── Phase 1: one-time setup (comment out after first run) ────────────────────
$PY src/scripts/run_analysis.py --step data_check
$PY src/scripts/run_analysis.py --step ref_index

# Run 2–3 pilot samples through align + dedup to get BAMs for panel derivation
PILOT=temp/pilot
$PY src/scripts/run_analysis.py --mode plasma --step preprocess \
    /mnt/g/OC_tissue/P10001.R1.fq.gz \
    /mnt/g/OC_tissue/P10001.R2.fq.gz \
    --workdir $PILOT
$PY src/scripts/run_analysis.py --mode plasma --step align --workdir $PILOT
$PY src/scripts/run_analysis.py --mode plasma --step dedup --workdir $PILOT

# Step 3: derive panel BED from those dedup BAMs (no re-alignment needed)
$PY src/scripts/derive_panel_bed.py --samples 5 --workdir $PILOT \
  --out temp/panel_P10001.bed
# ─────────────────────────────────────────────────────────────────────────────

# ── Phase 2: full analysis (all samples) ─────────────────────────────────────
$PY src/scripts/run_analysis.py --mode plasma --step preprocess /mnt/f/meth/OC_plasma/
WORKDIR=$(ls -dt temp/20[0-9]* | head -1)
echo "workdir: $WORKDIR"

$PY src/scripts/run_analysis.py --mode plasma --step align          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step dedup          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step target         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step mbias          --workdir $WORKDIR
# mbias writes mbias_auto_trim.json with suggested trim values (see workdir/).
# Review SVGs if you want to override. nOT/nOB = "auto" in config.toml applies them.

$PY src/scripts/run_analysis.py --mode plasma --step extract        --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step filter         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step diff           --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step visualize      --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step outlier_detect --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step feature_analysis --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step blood_filter     --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step classify         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step fs_benchmark     --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step leakage_demo     --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step final_model      --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step subtype_classify --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step temporal_validate --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step annotate        --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step external_validate --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step external_cohort   --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step report          --workdir $WORKDIR

echo "=== Done! Report: result/analysis_report.html ==="
```

## Common options

| Option | Description |
| --- | --- |
| `--mode plasma\|tissue` | Pipeline mode: plasma (bwameth + UMI dedup) or tissue (default: plasma) |
| `--step <name>` | Run a single step only |
| `--workdir <path>` | Existing workdir (required for single-step mode) |
| `--skip-report` | Skip report generation |
| `--clean-temp` | Delete temp files after completion |
| `--list-steps` | Print step list with completion status |
| `--config <path>` | Use a custom config file |

For the full step inventory (descriptions + tools per step) see [pipeline.md](pipeline.md).
