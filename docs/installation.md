# Installation (Conda / Local Development)

> **Linux / WSL2 only** — bioinformatics tools (bismark, bwameth, MethylDackel, etc.) do not support native Windows environments.

For Docker-based deployment, see [../DOCKER.md](../DOCKER.md) instead.

## 1. Create the conda environment

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

`[ml]` = xgboost + lightgbm; `[jupyter]` = jupyterlab + ipykernel; `[dev]` = testing tools.

### 1a. Protected core (optional, required for Steps 14 / 15 / 19)

Steps 14 (feature selection), 15 (diagnostic classifiers), and 19 (external cross-validation) depend on the `aimethflow-core` package, temporarily held in a private companion GitHub repository during the pre-publication phase.

**For local conda development (this guide):** journal-assigned reviewers receive collaborator access to the private repository from the corresponding author and install editable:

```bash
git clone https://github.com/<user>/aimethflow-core.git
cd aimethflow-core
pip install -e .[ml]
```

After installation the public stubs transparently substitute the real implementations at import time — no further configuration needed.

**For Docker deployments:** the wheel-mount path is documented in [../DOCKER.md](../DOCKER.md) → *Protected core via mounted wheel*. Docker users do NOT need GitHub collaborator access.

**Without `aimethflow-core`**: Steps 1-13, 16-18, and 20 still work end-to-end; attempting Step 14 / 15 / 19 raises an `ImportError` pointing to the corresponding author for access. See the *Status — pre-publication* section in [../README.md](../README.md) for details.

## 2. Reference genome index (first run only)

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

## 3. Auxiliary files

| File | Description |
| --- | --- |
| `data/OC_tissue_88_regions.hg38.txt` | Candidate methylation regions (88 tissue candidate regions, used as soft prior) |
| `data/ovarian_fq_series_samples_v4.xlsx` | Sample grouping Excel (contains classification / subtype columns) |
| `data/plasma_exp.xlsx` | Plasma sample collection dates (columns: `case_id`, `datetime` in `YYYYMMDD` format); required for temporal split validation in Step 17 (plasma mode only) |
| `/mnt/f/meth/hg38.repeat.bed` | hg38 repeat element BED |
| `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` | Promoter region BED |

### 3a. External reference files (optional, Step 19 only)

Step 19 (`external_validate`) cross-references panel CpGs against the GSE40279 whole-blood baseline and the TCGA-OV HumanMethylation450 tumour cohort. **4 reference files (~1.4 GB total) must live under `data/external/`**. Fetch them once with:

```bash
bash download_external_data.sh                # writes into data/external/
```

Manual sources and troubleshooting notes are in [pipeline.md → Required files](pipeline.md#required-files-place-in-dataexternal). All other pipeline steps run without these files — you can defer this setup until you actually want to run `--step external_validate`.

## Verification

Once installed, validate the configuration prerequisites:

```bash
run-analysis --step data_check
```

This confirms the sample Excel, reference genome, panel BED, and FASTQ paths declared in `config.toml` are reachable.

Next: head to [quickstart.md](quickstart.md) for the first analysis run.
