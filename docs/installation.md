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

## Verification

Once installed, validate the configuration prerequisites:

```bash
run-analysis --step data_check
```

This confirms the sample Excel, reference genome, panel BED, and FASTQ paths declared in `config.toml` are reachable.

Next: head to [quickstart.md](quickstart.md) for the first analysis run.
