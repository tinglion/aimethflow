# AIMethFlow — Ovarian Cancer Plasma cfDNA Targeted Methylation Analysis

Early diagnosis project for ovarian cancer — targeted methylation sequencing analysis pipeline based on plasma cfDNA, covering the full workflow from raw FastQ to machine-learning diagnostic classification (18 steps).

## Analysis Goals

- Detect differential methylation at CpG sites associated with ovarian cancer in plasma cell-free DNA (cfDNA)
- Data-driven discovery of candidate sites (individual-sample bedGraph direct aggregation; 88 candidate regions used as soft prior weighting)
- Train benign/malignant binary classification models (LR / SVM / RF / XGB / LGBM, Nested LOO-CV evaluation, feature selection performed within each fold's training set)
- Subtype multi-classification of malignant + borderline samples (high-grade serous carcinoma / endometrioid carcinoma / clear cell carcinoma / low-grade serous carcinoma / mucinous)

---

## Library Preparation Prerequisites

Before configuring the pipeline, confirm the following library parameters — they determine several critical step settings.

| Parameter | This Lab | Affected Steps |
| --- | --- | --- |
| UMI | **6N at R2 5' end** | Step 4 extraction, Step 6 deduplication method |
| Library type | **Single-strand, highly directional (OB ≈ 1.2%)** | Step 5 `non_directional = false` |
| End repair | **None; P7 adapter introduces poly-C (unmethylated)** | Step 4 `trim_tail2`, Step 8 M-bias check |
| 5' adapters | **12 nt low-complexity oligonucleotides on both R1 and R2** | Step 4 `trim_front1/2 = 12` |
| Reference genome | hg38 (`/data/meth/bwa/hg38.fa`) | All steps |
| bwameth c2t index | `/data/meth/bwa/hg38.fa.bwameth.c2t` (built, ~17 GB) | Step 5 |
| Panel BED | `data/panel_B10308.bed` | Step 7 on-target filtering |
| Candidate regions | `data/OC_tissue_88_regions.hg38.txt` (88 regions) | Step 9 filtering |

---

## Methodological Improvements

Key differences vs. the prior (Youxun) approach:

| Issue | Prior Approach | This Pipeline |
| --- | --- | --- |
| **Deduplication** | Picard coordinate dedup — removes genuine cfDNA fragments (147–167 bp) | **umi_tools dedup --paired**, UMI-aware |
| **UMI** | UMI built but ignored | fastp `--umi` extraction + umi_tools dedup |
| **On-target filter** | None | **bedtools intersect panel.bed** |
| **Poly-C tail** | R2 3' not trimmed | `trim_tail2=5`, confirmed by M-bias SVG |
| **Methylation threshold** | 80% (stricter than tissue) | **60%**, same as tissue |
| **Differential filter** | p<0.05 only | **p<0.05 AND \|Δβ\|>0.1** |

---

## Environment Requirements

> **Linux / WSL2 only** — bioinformatics tools (bismark, bwameth, MethylDackel, etc.) do not support native Windows environments.

### Option A — Docker (recommended, out-of-the-box)

**Prerequisites**: Docker ≥ 24, Docker Compose Plugin (`docker compose` command available)

```bash
# 1. Copy and fill in paths
cp .env.example .env
# Edit .env: set the host machine's data / reference genome / results directories

# 2. Build the image (first time ~10–20 minutes; cached on subsequent runs)
docker compose --progress=plain build

# 3. Start the persistent container
docker compose up -d
# Submit an analysis job
docker compose exec analysis \
    python /app/src/scripts/run_analysis.py --steps 1-19
# Enter the container shell
docker compose exec analysis bash

# One-shot run (removed after completion)
docker compose --profile run run --rm analysis-run \
    python /app/src/scripts/run_analysis.py --steps 1-19

# 4. Start JupyterLab to browse results step by step
docker compose --profile jupyter up jupyter
#    Open http://localhost:8888 in your browser
```

For details (mounted volumes, tool versions, config.docker.toml path mappings) → [DOCKER.md](DOCKER.md)

---

### Option B — Conda (local development)

#### 1. Create the conda environment

```bash
conda create -n meth python=3.10.14 -y
conda activate meth

# Bioinformatics tools (bioconda)
conda install -y -c bioconda -c conda-forge \
    fastp=0.23.4 \
    bwameth \
    samtools \
    bedtools \
    methyldackel \
    bismark \
    bowtie2 \
    umi_tools \
    pysam=0.22.0

# Python dependencies (install from pyproject.toml, single source of truth)
pip install -e ".[ml]" -i https://pypi.tuna.tsinghua.edu.cn/simple
# Also includes Jupyter: pip install -e ".[ml,jupyter]"
# Dev tools:            pip install -e ".[ml,jupyter,dev]"
```

`[ml]` = xgboost + lightgbm; `[jupyter]` = jupyterlab + ipykernel; `[dev]` = testing tools

#### 2. Reference genome index (first run only)

```bash
# bwameth index (plasma pipeline, recommended)
bwameth.py index /data/meth/bwa/hg38.fa

# bismark index (tissue pipeline)
bismark_genome_preparation /data/meth/bwa/
```

bwameth c2t index files (all required, total ~17 GB):

| File | Size | Purpose |
| --- | --- | --- |
| `hg38.fa.bwameth.c2t` | 6.1 GB | C→T converted reference |
| `hg38.fa.bwameth.c2t.bwt` | 6.0 GB | BWT index (alignment core) |
| `hg38.fa.bwameth.c2t.sa` | 3.0 GB | Suffix array (seed anchoring) |
| `hg38.fa.bwameth.c2t.pac` | 1.5 GB | 2-bit compressed sequence (CIGAR) |
| `hg38.fa.bwameth.c2t.amb` | 42 KB | Ambiguous base table |
| `hg38.fa.bwameth.c2t.ann` | 44 KB | Chromosome annotations (SAM header) |

#### 3. Auxiliary files

| File | Description |
| --- | --- |
| `data/OC_tissue_88_regions.hg38.txt` | Candidate methylation regions (88 tissue candidate regions, used as soft prior) |
| `data/ovarian_fq_series_samples_v4.xlsx` | Sample grouping Excel (contains classification / subtype columns) |
| `data/plasma_exp.xlsx` | Plasma sample collection dates (columns: `case_id`, `datetime` in `YYYYMMDD` format); required for temporal split validation in Step 17 (plasma mode only) |
| `/mnt/f/meth/hg38.repeat.bed` | hg38 repeat element BED |
| `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` | Promoter region BED |

---

## Project Structure

```text
├── config.toml                    # Global configuration
├── requirements.txt               # pip dependencies
├── pyproject.toml                 # Package configuration
├── data/                          # Reference data (BED, Excel)
├── src/
│   ├── methylation/               # Core analysis modules
│   │   ├── diff_methylation.py    # Differential methylation analysis
│   │   ├── sample_grouping.py     # Excel sample grouping
│   │   └── methylation_visualizer.py
│   ├── scripts/
│   │   ├── run_analysis.py        # Main entry point (CLI)
│   │   ├── webui.py               # Web UI (FastAPI, port 9999)
│   │   └── modules/
│   │       ├── alignment_utils.py     # Alignment / dedup / on-target filtering / extraction
│   │       ├── fastq_utils.py         # fastp preprocessing
│   │       ├── feature_utils.py       # Feature selection (Step 14) + outlier detection (Step 13)
│   │       ├── classifier_utils.py    # Diagnostic classification (Step 15)
│   │       ├── subtype_utils.py       # Subtype classification (Step 16) + temporal validation (Step 17)
│   │       ├── annotation_utils.py    # Site annotation (Step 18, Ensembl + EWAS + CpG island)
│   │       └── report_utils.py        # Report generation (Step 20)
│   ├── external_validate.py           # Public dataset cross-validation (Step 19, optional)
│   ├── test/
│   │   └── test_parse_sample_info.py
│   └── utils/
│       ├── config.py
│       └── logger.py
└── docs/
```

---

## Analysis Pipeline

### Plasma mode (`--mode plasma`, default)

| Step | Key | Description | Tool |
| --- | --- | --- | --- |
| 1 | `data_check` | Validate config: sample Excel, ref genome, panel BED, FASTQ paths | Built-in |
| 2 | `ref_index` | Build / verify bwameth hg38 c2t index (~17 GB, one-time setup) | bwameth |
| 3 | `panel_derive` | Derive panel BED from representative FASTQ alignments | bwameth + bedtools |
| 4 | `preprocess` | fastp QC + UMI extraction | fastp |
| 5 | `align` | Paired-end alignment (two-step, avoids pipe deadlock) | bwameth (default) / bismark |
| 6 | `dedup` | UMI-aware deduplication | umi_tools dedup --paired |
| 7 | `target` | On-target filtering | bedtools intersect |
| 8 | `mbias` | M-bias check; auto-detects suggested trim values → `mbias_auto_trim.json` | MethylDackel mbias |
| 9 | `extract` | CpG methylation extraction | MethylDackel extract --mergeContext |
| 10 | `filter` | Filter: repeat regions / depth / candidate regions | Custom + bedtools |
| 11 | `diff` | Differential methylation (Wilcoxon, p<0.05 and \|Δβ\|>0.1) | Custom |
| 12 | `visualize` | Density / scatter / heatmap plots | matplotlib |
| 13 | `outlier_detect` | Sample outlier detection (coverage / PCA z-score / Isolation Forest); suspect samples excluded from feature selection | scikit-learn |
| 14 | `feature_select` | Feature selection (LASSO + RF union); outlier samples excluded from computation; hard-filter sites with <80% samples at ≥10x depth | scikit-learn |
| 15 | `classify` | Diagnostic classifiers, nested CV. Auto-evaluates **11 models** in the same CV loop: 5 cfDNA-only (LR/SVM/RF/XGB/LGBM), CA-125 alone (clinical baseline, ≥35 U/mL cutoff), and 5 combined cfDNA+CA-125 fusion models (`LR+CA125` … `LGBM+CA125`). Disable the CA-125 paths with `ml.include_ca125_baseline = false`. | LR / SVM / RF / XGB / LGBM |
| 15b | `fs_benchmark` | **Optional.** Feature selection ablation: compare Wilcoxon-only / LASSO-only / RF-only / LASSO∪RF AUC on same CV | scikit-learn |
| 15c | `ca125_compare` | **Optional verification.** CA-125 AUC computed two ways (in-loop pseudo-classifier vs Excel-direct lookup) for cross-checking. Mostly redundant after Step 15 — useful as an audit. | scikit-learn |
| 16 | `subtype_classify` | Subtype multi-class LOO-CV (malignant + borderline) | RF + SVM OvR |
| 17 | `temporal_validate` | **Plasma only.** Temporal split validation: train on samples collected before 2024-07-01, evaluate on samples collected on/after that date. Same 11-model lineup as Step 15 (5 cfDNA-only + CA-125 baseline + 5 cfDNA+CA-125 fusion) so head-to-head comparison carries through to the prospective evaluation. Covers binary + subtype classification. Requires `data/plasma_exp.xlsx` (columns: `case_id`, `datetime`). | scikit-learn / XGB / LGBM |
| 18 | `annotate` | Panel CpG annotation: gene context / CpG island / EWAS Catalog / EWAS Atlas / PubMed plasma hits | Ensembl REST + UCSC + EWAS Catalog + CNCB NGDC + NCBI |
| 19 | `external_validate` *(optional)* | Cross-validate panel sites against public 450K datasets (blood baseline + TCGA-OV) | numpy / scipy |
| 20 | `report` | HTML + JSON summary report (annotation + external-validation sections auto-embedded) | Custom |

### Tissue mode (`--mode tissue`)

Steps 4–12 available (bismark + deduplicate_bismark); steps 13–19 require the sample info Excel with benign/malignant groups.

---

## Quick Start

### Web UI (recommended)

Visual interface for configuring, running steps, and browsing results and reports.

**Docker mode** (auto-starts after `docker compose up -d`):

```text
http://<server-ip>:9999
```

**Local development mode**:

```bash
# Run from the project root; METH_WORKDIR tells the Web UI where to find config.toml / temp / result
export METH_WORKDIR=$(pwd)          # Linux / macOS / WSL
# $env:METH_WORKDIR = (Get-Location).Path  # PowerShell

conda activate meth
python src/scripts/webui.py
# Open http://localhost:9999 in your browser
```

> The Web UI reads and writes `$METH_WORKDIR/config.toml` — in local dev mode this is the project root `config.toml`.
> Inside the Docker container it uses `config.docker.toml`; the two paths are independent.

---

### Full pipeline (one command)

```bash
run-analysis --mode plasma /mnt/f/meth/OC_plasma/
```

### Step-by-step execution

#### Phase 1 — one-time setup (per environment / per panel design)

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

#### Phase 2 — per-run analysis (all samples)

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

# Step 14: feature selection (outlier samples automatically excluded from LASSO/RF)
run-analysis --mode plasma --step feature_select --workdir $WORKDIR

# Steps 15–16: classification
run-analysis --mode plasma --step classify         --workdir $WORKDIR
# Step 15b (optional): feature selection ablation — wilcoxon / lasso / rf / union AUC comparison
run-analysis --mode plasma --step fs_benchmark     --workdir $WORKDIR
# Step 15c (optional): CA-125 vs classifier AUC — requires 'ca125' column in sample Excel
run-analysis --mode plasma --step ca125_compare    --workdir $WORKDIR
run-analysis --mode plasma --step subtype_classify --workdir $WORKDIR

# Step 17: temporal split validation (plasma only; requires data/plasma_exp.xlsx)
run-analysis --mode plasma --step temporal_validate --workdir $WORKDIR

# Step 18: annotate selected panel sites (network access required on first run; cached after)
run-analysis --mode plasma --step annotate          --workdir $WORKDIR
# Step 19 (optional): cross-validate against public datasets (requires data/external/)
run-analysis --mode plasma --step external_validate --workdir $WORKDIR
# Step 20 (last): generate report — auto-embeds annotation + external-validation sections
run-analysis --mode plasma --step report            --workdir $WORKDIR

# Show all steps and current progress
run-analysis --list-steps --workdir $WORKDIR
zip -r 20260406_report.zip result/20260406/20260406_report/
```

### Shell Script (background execution)

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
$PY src/scripts/run_analysis.py --mode plasma --step feature_select --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step classify         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step fs_benchmark     --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step ca125_compare    --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step subtype_classify --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step temporal_validate --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step annotate        --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step external_validate --workdir $WORKDIR
$PY src/scripts/run_analysis.py --mode plasma --step report          --workdir $WORKDIR

echo "=== Done! Report: result/analysis_report.html ==="
```

### Common options

| Option | Description |
| --- | --- |
| `--mode plasma\|tissue` | Pipeline mode: plasma (bwameth + UMI dedup) or tissue (default: plasma) |
| `--step <name>` | Run a single step only |
| `--workdir <path>` | Existing workdir (required for single-step mode) |
| `--skip-report` | Skip report generation |
| `--clean-temp` | Delete temp files after completion |
| `--list-steps` | Print step list with completion status |
| `--config <path>` | Use a custom config file |

---

## Config File `config.toml`

### Key paths

```toml
[data]
plasma_dir         = "/mnt/f/meth/OC_plasma"
sample_info_excel  = "data/ovarian_fq_series_samples_v3.xlsx"
probe_panel_bed    = "data/panel_B10308.bed"
candidate_regions  = "data/OC_tissue_88_regions.hg38.txt"

[extraction]
reference_genome   = "/data/meth/bwa/hg38.fa"

[tools]
conda_env_bin      = "/path/to/miniconda3/envs/meth/bin"
```

### Plasma parameters `[plasma]`

```toml
aligner       = "bwameth"   # or "bismark"
umi_length    = 6
trim_front1   = 12          # R1 5' low-complexity oligonucleotides
trim_front2   = 12          # R2 5' oligonucleotides after UMI
trim_tail2    = 5           # R2 3' poly-C trim — adjust after reviewing M-bias SVG
nOT           = "auto"      # "auto" → mbias_auto_trim.json; or explicit skip counts e.g. "0,0,0,39"
nOB           = "auto"      # "auto" → mbias_auto_trim.json; or explicit skip counts e.g. "0,0,17,0"
mbias_trim_percentile = 95  # percentile for aggregating trim values across samples (50=median, 100=max)
non_directional = false     # confirmed highly directional (OB ≈ 1.2%)
min_mapq      = 40
methylated_read_threshold = 0.60   # same threshold as tissue analysis
diff_delta_beta = 0.10      # |Δβ| threshold for differential analysis
```

> **M-bias auto-detection**: run `--step mbias` first — it writes `mbias_auto_trim.json` with skip counts ready for `--nOT`/`--nOB` in the `extract` step. MethylDackel outputs *inclusion positions* (`--OT a,b,c,d`), not skip counts; the pipeline converts them automatically: 5' skip = position directly; 3' skip = `read_len − pos − 1` (e.g. `--OT 0,0,0,110` + `read_len=150` → `nOT "0,0,0,39"`). Set `nOT = "auto"` / `nOB = "auto"` in `config.toml [plasma]` to apply; use explicit skip-count strings (e.g. `"0,0,0,39"`) to override.

### M-bias Re-run Scope

| Scenario | Adjust | Steps to re-run |
| --- | --- | --- |
| `trim_tail2` too small | Increase `trim_tail2` in `[plasma]` | `preprocess` → `dedup` → `target` → `extract` |
| Positional bias visible | Adjust `nOT` / `nOB` in `[plasma]` | `extract` onwards (no re-alignment) |

### Intermediate file cleanup `[cleanup]`

```toml
delete_after_dedup  = false   # delete sorted.bam after dedup
delete_after_target = false   # delete dedup.bam after on-target filter
```

### Thread configuration `[threads]`

```toml
fastp           = 8
bwameth         = 8
samtools_sort   = 4
methyldackel    = 4
```

### Machine learning `[ml]`

```toml
# Feature selection
feature_selection_mode = "union"   # "union" | "lasso" | "rf" | "wilcoxon"
min_feature_coverage   = 0.5       # keep CpG sites covered in ≥50% of samples

# Two-tier Wilcoxon pre-filter (relaxed thresholds for 88 candidate regions)
sample_wilcox_p        = 0.10      # p threshold (non-88 sites)
sample_delta_beta      = 0.10      # |Δβ| threshold (non-88 sites)
sample_wilcox_p_88     = 0.10      # p threshold (88 candidate region sites, relaxed)
sample_delta_beta_88   = 0.05      # |Δβ| threshold (88 candidate region sites, relaxed)

# Overfitting control
max_features_clf  = 15             # max features entering classifiers (≈ n_samples/10)
lasso_cv_folds    = 5

# Coverage QC hard filter (applied at end of feature_select)
min_10x_coverage_rate = 0.80       # remove panel sites with <80% of samples at ≥10x depth

# Sample outlier detection (Step 13)
outlier_coverage_min = 0.30        # flag if <30% of CpG sites are covered
outlier_pca_z        = 2.5         # flag if within-group PCA distance z-score exceeds this
exclude_samples      = []          # fill after reviewing outlier_report.csv

# CA-125 baseline + cfDNA+CA-125 combined models in Step 15 (auto-skipped if no `ca125`
# column in the sample info Excel). Set to false to evaluate cfDNA classifiers only.
include_ca125_baseline = true
ca125_threshold        = 35.0      # clinical cutoff (U/mL) for CA-125 hard predictions

# Classifier hyperparameters
rf_n_estimators   = 500
lr_C              = 1.0
svm_C             = 1.0
xgb_n_estimators  = 200
xgb_max_depth     = 4
xgb_learning_rate = 0.1
lgbm_n_estimators = 200
lgbm_max_depth    = 4
lgbm_learning_rate = 0.1
```

### Feature Annotation `[annotation]`

```toml
ewas_window      = 500    # EWAS Catalog search window ±bp around each CpG site
gene_window      = 2000   # Ensembl gene overlap search window ±bp
api_delay        = 0.25   # seconds between API calls (respects rate limits)
api_timeout      = 15     # HTTP request timeout (seconds)
cpgi_bed         = ""     # path to local hg38 CpG island BED; empty = auto-download from UCSC
ewas_atlas       = true   # query EWAS Atlas (CNCB NGDC) for positional associations
pubmed_search    = true   # search PubMed: "{gene} plasma methylation cancer"
ncbi_api_key     = ""     # NCBI API key (optional; raises limit 3→10 req/s)
pubmed_max_hits  = 20     # max PubMed results per gene
```

API responses are cached in `result/annotation_cache.json` — re-runs are instant.

---

## Output Files

| Dir | File | Description |
| --- | --- | --- |
| `result/<ts>/<ts>_report/` | `analysis_report.html` | Full pipeline HTML report (self-contained, ready to share) |
| `result/<ts>/<ts>_report/` | `analysis_report.json` | Structured report data |
| `result/<ts>/<ts>_report/imgs/` | `*.png`, `*.svg` | All images embedded in the report |
| `result/` | `diff_methylation.txt` | Differentially methylated CpG list |
| `result/` | `wilcoxon_stats.csv` | Sample-level Wilcoxon pre-filter results |
| `result/` | `feature_importance.csv` | LASSO + RF feature scores |
| `result/` | `panel_cpg.csv` | Final selected CpG panel |
| `result/` | `panel_cpg_depth.csv` | Per-panel-site depth stats: `coverage_rate`, `n_10x`, `coverage_rate_10x`, `pass_coverage_qc` |
| `result/` | `feature_matrix.npz` | Feature matrix cache (reused by Step 15 classify) |
| `result/` | `feature_selection_freq.csv` | LOO-CV per-fold feature selection frequency |
| `result/` | `outlier_report.csv` | Outlier report (coverage / PCA z / Isolation Forest) |
| `result/` | `classifier_summary.csv` | AUC / Acc / Sens / Spec for **up to 11 models** when CA-125 data present: 5 cfDNA-only + CA-125 baseline + 5 `*+CA125` combined models. Includes `n_used` column showing how many samples contributed to each row (full cohort for cfDNA-only, CA-125-available subset for CA-125 / combined). |
| `result/` | `feature_selection_benchmark.csv` | Feature selection ablation: AUC / Acc / Sens / Spec for wilcoxon / lasso / rf / union × 5 classifiers |
| `result/` | `ca125_comparison.csv` | (Step 15c verification) CA-125 vs classifier AUC computed via Excel-direct lookup; cross-checks the in-loop CA-125 path |
| `result/` | `ca125_values.csv` | (Step 15c) Per-sample CA-125 value + true label for the matched subset |
| `result/` | `ca125_auc_inputs.csv` | (Step 15c debug cache) Full join trace for the in-loop CA-125 evaluation: sample_id → case_id → ca125 + final y/score arrays passed to `roc_auc_score` |
| `result/` | `classifier_predictions.csv` | Per-sample predictions; one row per (sample_id × model). For CA-125 / `*+CA125`, only samples with a CA-125 measurement are included |
| `result/` | `roc_curve.png` | All-models ROC overlay (cfDNA solid, CA-125 dashed grey, `*+CA125` dotted same colour as base) |
| `result/` | `confusion_matrix.png` | Confusion matrices (auto-wraps to 4 rows × 3 cols when 11 models are evaluated) |
| `result/` | `temporal_validation_summary.csv` | Temporal split validation metrics (plasma only): AUC / Acc / Sens / Spec per classifier on samples collected on/after 2024-07-01. Includes 5 cfDNA-only + CA-125 baseline + 5 `*+CA125` combined models when CA-125 data is available; `n_validation` reflects what each row was evaluated on (full validation cohort for cfDNA / combined; matched subset for CA-125 alone). |
| `result/` | `temporal_validation_predictions.csv` | Per-sample predictions for temporal validation set (plasma only); one row per (sample × model). CA-125 / `*+CA125` rows omit samples without CA-125 measurements. |
| `result/` | `subtype_summary.csv` | Subtype classifier RF + SVM summary |
| `result/` | `subtype_predictions.csv` | Per-sample subtype predictions (with confidence) |
| `result/` | `subtype_roc.png` | Subtype OvR ROC curves |
| `result/` | `subtype_temporal_summary.csv` | Temporal split subtype validation metrics (plasma only): accuracy per classifier |
| `result/` | `subtype_temporal_predictions.csv` | Per-sample subtype predictions for temporal validation set (plasma only) |
| `result/` | `panel_annotation.csv` | Per-site annotation: gene / CpG island / EWAS Catalog / EWAS Atlas / `has_plasma_evidence` / PubMed hits |
| `result/` | `annotation_cache.json` | Cached API responses (speeds up re-runs) |
| `result/` | `external_validation.csv` | P3 cross-validation: per-site blood baseline + TCGA-OV β stats (`blood_mean_beta`, `blood_pass`, `tcga_mean_beta`, `tcga_n_samples`) |
| `temp/` | `*.bwameth_pe.sorted.bam` | Sorted alignment BAM |
| `temp/` | `*.dedup.bam` | UMI-deduplicated BAM |
| `temp/` | `*_CpG.bedGraph` | Per-sample methylation bedGraph |
| `temp/` | `*.mbias_OT.svg` | M-bias plot (OT strand) |
| `temp/` | `mbias_auto_trim.json` | Auto-detected trim values (nOT/nOB) aggregated across all samples |

---

## Sample Excel Format

`data/ovarian_fq_series_samples_v4.xlsx` must contain the following columns (multiple sheets supported):

| Column | Content |
| --- | --- |
| `fq_R1` | R1 filename, e.g. `B10001.R1.fq.gz` |
| `fq_R2` | R2 filename |
| `category` | `benign` / `malignant` / `borderline` |
| `subtype` | Subtype, e.g. `high-grade serous carcinoma` (fill in for malignant / borderline samples) |

- Borderline samples are excluded from binary classification training but are included in subtype classification
- Duplicate entries for the same sample across multiple sheets are automatically deduplicated

---

## Resource Requirements

Reference values measured on 8-core / 32 GB machine (N = 16 M read pairs):

| Step | Tool | Wall time | Peak memory | Output size |
| --- | --- | --- | --- | --- |
| preprocess | fastp | ~2 min | 0.5 GB | ×0.90 of raw FASTQ |
| align | bwameth + samtools | **~47 min** | **13–14 GB** | ≈ clean FASTQ × 0.85 |
| dedup | umi_tools | **~10 min** | 2–4 GB | ×0.38–0.50 of sorted BAM |
| target | bedtools | ~4 min | 1–2 GB | ×0.25–0.36 of dedup BAM |
| mbias | MethylDackel | ~5 min | 2–4 GB | ~2 MB (SVG) |
| extract | MethylDackel | ~10 min | 2–4 GB | ~20–50 MB |
| filter + diff + visualize | Python | <5 min | 1 GB | <5 MB |
| feature_select | scikit-learn | 5–30 min | 2–4 GB | <10 MB |
| classify | sklearn / XGB / LGBM | 10–60 min | 2–8 GB | <5 MB |
| subtype_classify | scikit-learn | 5–30 min | 2–4 GB | <5 MB |
| report | Python | <2 min | 1 GB | ~5–15 MB |
| annotate | urllib (REST APIs) | 2–30 min* | <1 GB | <1 MB |

\* `annotate` time depends on panel size and API latency; cached on re-runs (instant).

**Memory bottleneck**: bwameth alignment (BWT 6 GB + SA 3 GB + working buffer ≈ 13–14 GB).

**Disk (single sample, N=16 M pairs)**:

- Peak (all intermediate files): ~5.8 GB
- Final retained (targeted BAM + results): ~2.3 GB

**Data flow (B10002, N=16 M pairs, measured)**:

```text
Raw FASTQ       1.88 GB  (R1+R2)
  ↓ fastp ×0.90
Clean FASTQ     1.80 GB
  ↓ bwameth ×0.79
Sorted BAM      1.50 GB
  ↓ umi_tools ×0.38
Dedup BAM       0.56 GB
  ↓ bedtools ×0.36
On-target BAM   0.20 GB  ← primary analysis data (11% of raw FASTQ)
  ↓ MethylDackel
CpG.bedGraph    ~1 MB
  ↓ 88-region filter
Candidate BG    <1 MB    ← final analysis data (<0.1%)
```

---

## Helper tools

```bash
# Quick FastQ statistics
fq-inspect /mnt/f/meth/OC_plasma/B10001.R1.fq.gz

# Derive panel BED from FASTQ files (Step 3, preferred)
python src/scripts/derive_panel_bed.py \
  --fq /mnt/f/meth/OC_plasma/*.R1.fq.gz \
  --out data/panel_B10308.bed --samples 5

# Derive panel BED from BAM (manual one-liner, if you already have aligned BAMs)
bedtools genomecov -ibam sample.dedup.bam -bg \
  | awk '$4>10' | bedtools merge -d 100 \
  | awk 'length($1)<6' > panel.bed

# Clean up large intermediate files after analysis
rm $WORKDIR/*.bwameth_pe.sorted.bam $WORKDIR/*.bwameth_pe.sorted.bam.bai
rm $WORKDIR/*.bwameth_pe.dedup.bam  $WORKDIR/*.bwameth_pe.dedup.bam.bai

# Unit tests
pytest src/test/ -v
```

### Step 19 — Public Dataset Cross-Validation (optional)

#### Validation Logic

Panel sites are mapped to public 450K datasets via the following four steps:

```text
        │
        ▼  ① liftOver hg38 → hg19 (UCSC chain file, per-chromosome block mapping)
Panel sites (hg19)
        │
        ▼  ② Match Illumina 450K cg probes within ±2 bp (GPL13534 probe coordinates)
cg probe ID (e.g. cg07881041)
        │
        ├─► ③ Look up GSE40279 whole-blood β matrix (mean of 656 samples) → blood_mean_beta
        └─► ④ Look up TCGA-OV β matrix (mean of ~590 samples) → tcga_mean_beta
```

**Blood baseline (GSE40279)**: cfDNA originates from plasma; blood cell cfDNA is the primary background noise. If a panel site is already highly methylated in normal whole blood (`blood_mean_beta ≥ 0.10`), the methylation signal detected in plasma may be primarily from blood cell contamination rather than tumor — the tumor specificity of that site is questionable and it should be considered for removal or down-weighting from the panel (`blood_pass = False` rows highlighted in yellow).

**Tumor signal (TCGA-OV)**: Confirms that panel sites have genuine methylation signal in ovarian cancer tissue (`tcga_mean_beta` is elevated), providing mutual corroboration with tissue-derived differential methylation results.

**Why does "Probe match: 0%" occur?**  
Coordinate system mismatch: panel sites are in hg38 while 450K array coordinates are in hg19, requiring liftOver conversion before matching. If liftOver succeeds (Liftover success shows ~96%) but probe matching is 0%, it means **`GPL13534_platform.txt` has not been placed in `data/external/`** — when the probe coordinate file is missing the probe_map is empty and no sites can be matched. Download it as described below and re-run.

Run after the main pipeline (Step 20 report) completes to validate panel sites against:

- **Blood baseline** (normal whole-blood cfDNA background methylation < 0.10) — based on GSE40279 (656 whole-blood samples)
- **Tumor signal** (methylation signal present in ovarian cancer tissue) — based on TCGA-OV HumanMethylation450 (~590 TCGA samples)

**Required files** (place in `data/external/`):

| File | Source | Description |
| --- | --- | --- |
| `GSE40279_series_matrix.txt.gz` | [GEO GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) → Download full table | Whole-blood Illumina 450K β matrix, 1.2 GB |
| `HumanMethylation450.gz` | [Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Ovarian%20Cancer%20(OV)) → DNA methylation → HumanMethylation450 → Download | TCGA-OV 450K β matrix, full version ~300 MB; click Download rather than preview |
| `hg38ToHg19.over.chain.gz` | [UCSC downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz) | Coordinate conversion chain file |
| `GPL13534_platform.txt` | [NCBI GEO GPL13534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534) → Data table → Download full table… → save as `GPL13534_platform.txt` | Illumina 450K probe hg19 coordinates (recommended); the script also accepts `HM450.hg19.manifest.tsv.gz` (Zhou lab, ~10 MB, auto-download attempted); **if the file is missing, the probe match column will be all False but the pipeline will not abort** |

**Run**:

```bash
# result-dir should be the timestamp directory for this analysis run
run-analysis --mode plasma --step external_validate --workdir result/20260406
# Equivalent (calling the script directly):
# python src/scripts/external_validate.py --result-dir result/20260406 --data-dir data
# Output: result/20260406/external_validation.csv, and automatically regenerates the HTML report (embeds Section 13)
```

**Output column descriptions**:

| Column | Meaning |
| --- | --- |
| `cg_probe` | Matched Illumina 450K probe ID (liftOver hg38→hg19, ±2 bp) |
| `probe_matched` | Whether a cg probe was successfully matched |
| `blood_mean_beta` | Mean β value across 656 whole-blood samples |
| `blood_pass` | `blood_mean_beta < 0.10` (True = meets low blood methylation criterion) |
| `tcga_mean_beta` | TCGA-OV mean β |
| `tcga_n_samples` | Number of TCGA-OV samples with coverage |

**Interpretation**: Sites with `blood_pass = False` may be subject to blood cell interference in plasma — consider removing or down-weighting them from the panel.

---

## Notes

- `pysam` is Linux/WSL only; cannot be installed on native Windows
- XGBoost / LightGBM are optional dependencies; if not installed they are automatically skipped and the remaining classifiers run normally
- Inter-step state is passed via `temp/<workdir>/pipeline_state.json` — no need to re-specify input files when re-running a single step
- Both the dedup and target steps support idempotent re-runs: if output files already exist they are automatically skipped
- **M-bias auto-detection**: `--step mbias` writes `mbias_auto_trim.json` with converted skip counts. MethylDackel outputs *inclusion positions* (`--OT a,b,c,d`); the pipeline converts them: 5' skip = position directly; 3' skip = `read_len − pos − 1` (e.g. `--OT 0,0,0,110` + `read_len=150` → `nOT "0,0,0,39"`). Set `nOT = "auto"` / `nOB = "auto"` to use in `extract`; explicit skip-count strings override auto. Changing `trim_tail2` requires re-running from `preprocess`; changing `nOT`/`nOB` only requires re-running from `extract`.
- **Outlier exclusion**: Step 13 generates the report and saves suspect sample IDs to pipeline state. Step 14 (`feature_select`) automatically excludes these samples from LASSO/RF computation — the saved `feature_matrix.npz` still contains all samples so classifiers can apply their own exclusion via `exclude_samples`. Add confirmed outlier IDs to `exclude_samples` and re-run `classify`.
- **CA-125 head-to-head (Step 15)**: with a `ca125` column (numeric, U/mL) in the sample info Excel, the `classify` step automatically evaluates 11 models in the same CV folds — 5 cfDNA-only, **CA-125 alone** (clinical baseline; raw value as score, ≥35 U/mL cutoff for accuracy/sens/spec), and **5 combined cfDNA + CA-125 fusion models** (`LR+CA125`, `SVM+CA125`, `RF+CA125`, `XGB+CA125`, `LGBM+CA125`). Combined models concatenate CA-125 (training-fold-median imputed, separately standardized) as an extra column to the cfDNA feature matrix — leakage-free. All 11 ROC curves render together in `roc_curve.png` (cfDNA solid, CA-125 dashed grey, combined dotted same colour as base). Disable via `ml.include_ca125_baseline = false`; tune the binary cutoff via `ml.ca125_threshold = 35.0`. The optional Step 15c `ca125_compare` is now a verification step that cross-checks the in-loop CA-125 AUC against an Excel-direct computation — keep it as an audit but it's no longer the primary CA-125 path.

- **Standalone CA-125 AUC** (`src/scripts/compute_ca125_auc.py`): for quick CA-125-only AUC without running any pipeline step. Reads `category` + `ca125` columns directly from the Excel, applies `ml.exclude_samples` (override with `--exclude` flag), prints AUC + 95% bootstrap CI. Useful as a sanity check against Step 15's in-loop CA-125 AUC; the two should agree to 3 decimal places.

- **Stratified k-fold CV correctness**: prior to v1.0.0 the per-sample `classifier_predictions.csv` had a label-shuffling bug under `cv_folds > 0` (StratifiedKFold with shuffle): `true_label` and `prob_malignant` were attached to the wrong `sample_id` because predictions were appended in fold-iteration order while the CSV writer iterated `sample_ids` in original order. **Aggregate AUC/Acc/Sens/Spec metrics were unaffected**, but downstream per-sample analyses (CA-125 join, confusion-matrix breakdowns, temporal-validation warm-start) saw scrambled labels. Fixed by pre-allocating per-classifier arrays of length `n_samples` and writing back at `test_idx` positions. LOO-CV (`cv_folds=0`) was unaffected because LeaveOneOut iterates 0→n-1. **If you have a workdir from before this fix, re-run `--step classify`** to regenerate `classifier_predictions.csv` with correctly aligned rows.
- **Feature selection ablation (Step 15b, optional)**: run `--step fs_benchmark` after `classify` to compare four feature selection strategies (wilcoxon-only / lasso-only / rf-only / lasso∪rf) using the same CV, classifiers, and cohort. Results in `feature_selection_benchmark.csv`; automatically embedded in Step 15 of the HTML report. In the 116-sample ovarian cancer cohort, ML-based methods (best AUC 0.804–0.809) substantially outperformed wilcoxon-only top-Δβ ranking (best AUC 0.759). `feature_selection_mode = "wilcoxon"` is also available in `config.toml` for single-mode runs.
- **Temporal split validation (Step 17, plasma only)**: run `--step temporal_validate` after Steps 15 and 16 to evaluate the pre-selected panel on prospectively collected samples. Requires `data/plasma_exp.xlsx` with columns `case_id` and `datetime` (format `YYYYMMDD`). Samples with `datetime < 2024-07-01` form the train set; samples on/after that date form the validation set. Mirrors Step 15's **11-model lineup** when a `ca125` column is present: 5 cfDNA-only classifiers + CA-125 baseline (≥35 U/mL cutoff, evaluated only on samples with a measurement) + 5 `*+CA125` combined models (CA-125 imputed on the train set with the training median, separately standardized, then concatenated as an extra column). Combined models give a prospective-deployment view of the marginal value of fusing CA-125 with cfDNA. Toggle via `ml.include_ca125_baseline`; threshold via `ml.ca125_threshold`. If the date Excel is absent the step logs a warning and exits cleanly. Outputs: `temporal_validation_summary.csv`, `temporal_validation_predictions.csv`, `subtype_temporal_summary.csv`, `subtype_temporal_predictions.csv`.
- The 88 candidate regions (`OC_tissue_88_regions.hg38.txt`) are now a **soft prior**: sites within these regions use more relaxed Wilcoxon thresholds (p<0.10, |Δβ|≥0.05) and are no longer a hard filter condition
- **Annotation (Step 18)**: runs before `report` so the annotation table is embedded in the HTML. Requires internet access on first run; all API responses (Ensembl, EWAS Catalog, EWAS Atlas, PubMed) are cached in `annotation_cache.json` — re-runs are instant. If `annotate` is skipped, the report section is simply omitted.
- **Annotation — NCBI API key**: set `[annotation] ncbi_api_key` to increase PubMed rate limit from 3 to 10 req/s. Without a key the search still works but is slower. Register free at [NCBI account portal](https://www.ncbi.nlm.nih.gov/account/).
- **Coverage QC (Step 14)**: at the end of `feature_select`, sites where fewer than `min_10x_coverage_rate` (default 80%) of samples have ≥10x depth are hard-removed from `panel_cpg.csv`. Details recorded in `panel_cpg_depth.csv` (`n_10x`, `coverage_rate_10x`, `pass_coverage_qc` columns). Re-run `classify` after adjusting the threshold — no need to re-run feature selection.
- **Public dataset cross-validation (Step 19)**: run `--step external_validate` after `annotate`, then run `--step report` last to embed Section 13 (blood baseline + TCGA-OV table) in the HTML. Sites with `blood_pass = False` (mean β ≥ 0.10 in whole blood) are highlighted yellow — consider removing or down-weighting them from the panel.
- **Relation to tissue analysis**: the 88 candidate regions derived from tissue DMR analysis serve as a soft prior — they cannot be directly applied to plasma. Plasma cfDNA has lower signal-to-noise; statistical validation in plasma data is always required before modelling
