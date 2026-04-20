# AIMethFlow — 卵巢癌血浆 cfDNA 靶向甲基化分析

卵巢癌早期诊断项目——基于血浆 cfDNA 的靶向甲基化测序分析流程，涵盖从原始 FastQ 到机器学习诊断分类的全流程（18 个步骤）。

## 分析目标

- 检测血浆游离 DNA（cfDNA）中与卵巢癌相关的 CpG 位点差异甲基化
- 数据驱动发现候选位点（个体样本 bedGraph 直接聚合，88 个候选区域作为软先验加权）
- 训练良恶性二分类模型（LR / SVM / RF / XGB / LGBM，Nested LOO-CV 评估，特征选择在每折训练集内完成）
- 对恶性 + 交界性样本进行亚型多分类（高级别浆液性癌 / 子宫内膜样癌 / 透明细胞癌 / 低级别浆液性癌 / 黏液性）

---

## Library Preparation Prerequisites

Before configuring the pipeline, confirm the following library parameters — they determine several critical step settings.

| Parameter | This Lab | Affected Steps |
| --- | --- | --- |
| UMI | **6N at R2 5' end** | Step 4 extraction, Step 6 deduplication method |
| Library type | **Single-strand, highly directional (OB ≈ 1.2%)** | Step 5 `non_directional = false` |
| End repair | **None; P7 adapter introduces poly-C (unmethylated)** | Step 4 `trim_tail2`, Step 7.5 M-bias check |
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

## 环境要求

> **仅支持 Linux / WSL2**，生信工具（bismark、bwameth、MethylDackel 等）不支持 Windows 原生环境。

### 1. 创建 conda 环境

```bash
conda create -n meth python=3.10.14 -y
conda activate meth

# 生信工具（bioconda）
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

# Python 依赖
pip install -e ".[ml]" -i https://pypi.tuna.tsinghua.edu.cn/simple
```

`[ml]` 包含 `xgboost` 和 `lightgbm`；不需要 ML 步骤时可省略：`pip install -e .`

### 2. 参考基因组索引（首次运行）

```bash
# bwameth 索引（血浆流程，推荐）
bwameth.py index /data/meth/bwa/hg38.fa

# bismark 索引（组织流程）
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

### 3. 辅助文件

| 文件 | 说明 |
| --- | --- |
| `data/OC_tissue_88_regions.hg38.txt` | 候选甲基化区域（88 个组织候选区域，作为软先验） |
| `data/卵巢癌fq编号与样本信息对应表v3.xlsx` | 样本分组 Excel（含分类/亚型列） |
| `/mnt/f/meth/hg38.repeat.bed` | hg38 重复序列 BED |
| `/mnt/f/meth/hg38_upstream_2k_st.bed.gz` | 启动子区域 BED |

---

## 项目结构

```text
├── config.toml                    # 全局配置
├── requirements.txt               # pip 依赖
├── pyproject.toml                 # 包配置
├── data/                          # 参考数据（BED、Excel）
├── src/
│   ├── methylation/               # 核心分析模块
│   │   ├── diff_methylation.py    # 差异甲基化分析
│   │   ├── sample_grouping.py     # Excel 样本分组
│   │   └── methylation_visualizer.py
│   ├── scripts/
│   │   ├── run_analysis.py        # 主入口（CLI）
│   │   └── modules/
│   │       ├── alignment_utils.py     # 比对 / 去重 / 靶区过滤 / 提取
│   │       ├── fastq_utils.py         # fastp 预处理
│   │       ├── feature_utils.py       # 特征选择（Step 14）+ 异常检测（Step 14.5）
│   │       ├── classifier_utils.py    # 诊断分类（Step 15）
│   │       ├── subtype_utils.py       # 亚型分类（Step 16）
│   │       ├── annotation_utils.py    # 位点注释（Step 17，Ensembl + EWAS + CpG island）
│   │       └── report_utils.py        # 报告生成（Step 19）
│   ├── external_validate.py           # 公共数据集交叉验证（Step 18，可选）
│   ├── test/
│   │   └── test_parse_sample_info.py
│   └── utils/
│       ├── config.py
│       └── logger.py
└── docs/
```

---

## 分析流程

### 血浆模式（`--plasma`，推荐）

| Step | Key | Description | Tool |
| --- | --- | --- | --- |
| 4 | `preprocess` | fastp QC + UMI extraction | fastp |
| 5 | `align` | Paired-end alignment (two-step, avoids pipe deadlock) | bwameth (default) / bismark |
| 6 | `dedup` | UMI-aware deduplication | umi_tools dedup --paired |
| 7 | `target` | On-target filtering | bedtools intersect |
| 7.5 | `mbias` | M-bias check (optional but recommended) | MethylDackel mbias |
| 8 | `extract` | CpG methylation extraction | MethylDackel extract --mergeContext |
| 9 | `filter` | Filter: repeat regions / depth / candidate regions | Custom + bedtools |
| 10 | `diff` | Differential methylation (Wilcoxon, p<0.05 and \|Δβ\|>0.1) | Custom |
| 11 | `visualize` | Density / scatter / heatmap plots | matplotlib |
| 14 | `feature_select` | Feature selection (LASSO + RF union); hard-filter sites with <80% samples at ≥10x depth | scikit-learn |
| 14.5 | `outlier_detect` | Sample outlier detection (coverage / PCA z-score / Isolation Forest) | scikit-learn |
| 15 | `classify` | Diagnostic classifiers, Nested LOO-CV | LR / SVM / RF / XGB / LGBM |
| 16 | `subtype_classify` | Subtype multi-class LOO-CV (malignant + borderline) | RF + SVM OvR |
| 17 | `annotate` | Panel CpG annotation: gene context / CpG island / EWAS Catalog / EWAS Atlas / PubMed plasma hits | Ensembl REST + UCSC + EWAS Catalog + CNCB NGDC + NCBI |
| 18 | `external_validate` *(optional)* | Cross-validate panel sites against public 450K datasets (blood baseline + TCGA-OV) | numpy / scipy |
| 19 | `report` | HTML + JSON summary report (annotation + external-validation sections auto-embedded) | Custom |

### 组织模式（无 `--plasma`）

步骤 4–11 可用（bismark + deduplicate_bismark）；步骤 14–18 需要 Excel 含良/恶性分组。

---

## 快速开始

### 全流程一键运行

```bash
run-analysis --plasma /mnt/f/meth/OC_plasma/
```

### 分步执行

```bash
# Step 4: preprocessing (creates workdir; subsequent steps reuse it)
run-analysis --plasma --step preprocess /mnt/f/meth/OC_plasma/
WORKDIR=$(ls -dt temp/20[0-9]* | head -1)   # capture generated workdir

# Steps 5–11: alignment / dedup / extraction / analysis
run-analysis --plasma --step align          --workdir $WORKDIR
run-analysis --plasma --step dedup          --workdir $WORKDIR
run-analysis --plasma --step target         --workdir $WORKDIR
run-analysis --plasma --step mbias          --workdir $WORKDIR
# ── review M-bias SVGs ($WORKDIR/*.mbias_OT.svg) before continuing ──
run-analysis --plasma --step extract        --workdir $WORKDIR
run-analysis --plasma --step filter         --workdir $WORKDIR
run-analysis --plasma --step diff           --workdir $WORKDIR
run-analysis --plasma --step visualize      --workdir $WORKDIR

# Step 14: feature selection
run-analysis --plasma --step feature_select --workdir $WORKDIR

# Step 14.5: outlier detection → review result/outlier_report.csv
run-analysis --plasma --step outlier_detect --workdir $WORKDIR
#   Add confirmed outlier IDs to config.toml [ml] exclude_samples = [...]

# Steps 15–17: classifiers + annotation
run-analysis --plasma --step classify         --workdir $WORKDIR
run-analysis --plasma --step subtype_classify --workdir $WORKDIR
# Step 17: annotate selected panel sites (network access required on first run; cached after)
run-analysis --plasma --step annotate         --workdir $WORKDIR
# Step 18 (optional): cross-validate against public datasets (requires data/external/)
run-analysis --plasma --step external_validate --workdir $WORKDIR
# Step 19 (last): generate report — auto-embeds annotation + external-validation sections
run-analysis --plasma --step report           --workdir $WORKDIR

# Show all steps and current progress
run-analysis --list-steps --workdir $WORKDIR
# 
zip -r 20260406_report.zip .\result\20260406\20260406_report\
```

### Shell Script (background execution)

```bash
#!/bin/bash
# Usage: bash run_oc_plasma.sh 2>&1 | tee logs/oc_plasma_run.log
set -e
cd /mnt/e/projects/hq/meth
export PYTHONPATH=src
PY=/mnt/e/programs/ubuntu/miniconda3/envs/meth/bin/python

$PY src/scripts/run_analysis.py --plasma --step preprocess /mnt/f/meth/OC_plasma/
WORKDIR=$(ls -dt temp/20[0-9]* | head -1)
echo "workdir: $WORKDIR"

$PY src/scripts/run_analysis.py --plasma --step align          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step dedup          --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step target         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step mbias          --workdir $WORKDIR

echo "=== Review M-bias SVGs ($WORKDIR/*.mbias_OT.svg). Press Enter to continue, Ctrl+C to adjust. ==="
read

$PY src/scripts/run_analysis.py --plasma --step extract        --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step filter         --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step diff           --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step visualize      --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step feature_select --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step classify       --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step subtype_classify --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step annotate       --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step external_validate --workdir $WORKDIR
$PY src/scripts/run_analysis.py --plasma --step report         --workdir $WORKDIR

# Step 19 (optional): cross-validate panel sites against public 450K datasets
# Requires data/external/ files — see "辅助工具 → Step 19" section
# Runs external_validate.py then auto-regenerates report with Section 13 embedded
$PY src/scripts/run_analysis.py --plasma --step external_validate --workdir $WORKDIR

echo "=== Done! Report: result/analysis_report.html ==="
```

### 常用选项

| Option | Description |
| --- | --- |
| `--plasma` | Plasma mode (bwameth + UMI dedup) |
| `--step <name>` | Run a single step only |
| `--workdir <path>` | Existing workdir (required for single-step mode) |
| `--skip-report` | Skip report generation |
| `--clean-temp` | Delete temp files after completion |
| `--list-steps` | Print step list with completion status |
| `--config <path>` | Use a custom config file |

---

## 配置文件 `config.toml`

### 关键路径

```toml
[data]
plasma_dir         = "/mnt/f/meth/OC_plasma"
sample_info_excel  = "data/卵巢癌fq编号与样本信息对应表v3.xlsx"
probe_panel_bed    = "data/panel_B10308.bed"
candidate_regions  = "data/OC_tissue_88_regions.hg38.txt"

[extraction]
reference_genome   = "/data/meth/bwa/hg38.fa"

[tools]
conda_env_bin      = "/path/to/miniconda3/envs/meth/bin"
```

### 血浆参数 `[plasma]`

```toml
aligner       = "bwameth"   # or "bismark"
umi_length    = 6
trim_front1   = 12          # R1 5' low-complexity oligonucleotides
trim_front2   = 12          # R2 5' oligonucleotides after UMI
trim_tail2    = 5           # R2 3' poly-C trim — adjust after reviewing M-bias SVG
nOT           = "0,0,0,20"  # MethylDackel OT-strand positional bias correction
nOB           = "0,0,17,0"  # MethylDackel OB-strand positional bias correction
non_directional = false     # confirmed highly directional (OB ≈ 1.2%)
min_mapq      = 40
methylated_read_threshold = 0.60   # same threshold as tissue analysis
diff_delta_beta = 0.10      # |Δβ| threshold for differential analysis
```

> **M-bias adjustment**: run `--step mbias` first and inspect the SVG. If the R2 3' methylation rate drops sharply, increase `trim_tail2` (e.g. 8 or 10) and re-run from `extract` (no re-alignment needed). If `nOT`/`nOB` need tuning, only `extract` onwards must be re-run.

### M-bias Re-run Scope

| Scenario | Adjust | Steps to re-run |
| --- | --- | --- |
| `trim_tail2` too small | Increase `trim_tail2` in `[plasma]` | `preprocess` → `dedup` → `target` → `extract` |
| Positional bias visible | Adjust `nOT` / `nOB` in `[plasma]` | `extract` onwards (no re-alignment) |

### 中间文件清理 `[cleanup]`

```toml
delete_after_dedup  = false   # delete sorted.bam after dedup
delete_after_target = false   # delete dedup.bam after on-target filter
```

### 线程配置 `[threads]`

```toml
fastp           = 8
bwameth         = 8
samtools_sort   = 4
methyldackel    = 4
```

### 机器学习 `[ml]`

```toml
# Feature selection
feature_selection_mode = "union"   # "union" | "lasso" | "rf"
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

# Sample outlier detection (Step 14.5)
outlier_coverage_min = 0.30        # flag if <30% of CpG sites are covered
outlier_pca_z        = 2.5         # flag if within-group PCA distance z-score exceeds this
exclude_samples      = []          # fill after reviewing outlier_report.csv

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

## 输出文件

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
| `result/` | `feature_matrix.npz` | Feature matrix cache (reused by Step 15) |
| `result/` | `feature_selection_freq.csv` | LOO-CV per-fold feature selection frequency |
| `result/` | `outlier_report.csv` | Outlier report (coverage / PCA z / Isolation Forest) |
| `result/` | `classifier_summary.csv` | 5 classifiers AUC / Acc / Sens / Spec |
| `result/` | `classifier_predictions.csv` | Per-sample predictions |
| `result/` | `roc_curve.png` | ROC curves comparison |
| `result/` | `confusion_matrix.png` | Confusion matrices (2 rows × 3 cols) |
| `result/` | `subtype_summary.csv` | Subtype classifier RF + SVM summary |
| `result/` | `subtype_predictions.csv` | Per-sample subtype predictions (with confidence) |
| `result/` | `subtype_roc.png` | Subtype OvR ROC curves |
| `result/` | `panel_annotation.csv` | Per-site annotation: gene / CpG island / EWAS Catalog / EWAS Atlas / `has_plasma_evidence` / PubMed hits |
| `result/` | `annotation_cache.json` | Cached API responses (speeds up re-runs) |
| `result/` | `external_validation.csv` | P3 cross-validation: per-site blood baseline + TCGA-OV β stats (`blood_mean_beta`, `blood_pass`, `tcga_mean_beta`, `tcga_n_samples`) |
| `temp/` | `*.bwameth_pe.sorted.bam` | Sorted alignment BAM |
| `temp/` | `*.dedup.bam` | UMI-deduplicated BAM |
| `temp/` | `*_CpG.bedGraph` | Per-sample methylation bedGraph |
| `temp/` | `*.mbias_OT.svg` | M-bias plot (OT strand) |

---

## 样本 Excel 格式

`data/卵巢癌fq编号与样本信息对应表v3.xlsx` 须包含以下列（支持多 sheet）：

| 列名 | 内容 |
| --- | --- |
| `fq_R1` | R1 文件名，如 `B10001.R1.fq.gz` |
| `fq_R2` | R2 文件名 |
| `分类` | `良性` / `恶性` / `交界性` |
| `亚型` | `高级别浆液性癌` 等（恶性/交界性填写） |

- `交界性` 不参与二分类训练，但参与亚型分类
- 同一样本在多个 sheet 出现时自动去重

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

## 辅助工具

```bash
# 快速查看 FastQ 统计
fq-inspect /mnt/f/meth/OC_plasma/B10001.R1.fq.gz

# 从 BAM 文件生成 panel.bed（若无捕获 BED）
bedtools genomecov -ibam sample.dedup.bam -bg \
  | awk '$4>10' | bedtools merge -d 100 \
  | awk 'length($1)<6' > panel.bed

# Clean up large intermediate files after analysis
rm $WORKDIR/*.bwameth_pe.sorted.bam $WORKDIR/*.bwameth_pe.sorted.bam.bai
rm $WORKDIR/*.bwameth_pe.dedup.bam  $WORKDIR/*.bwameth_pe.dedup.bam.bai

# 单元测试
pytest src/test/ -v
```

### Step 18 — 公共数据集交叉验证（可选）

#### 验证逻辑

panel 位点经以下四步映射到公共 450K 数据集：

```text
        │
        ▼  ① liftOver hg38 → hg19（UCSC chain 文件，每条染色体独立块映射）
Panel 位点 (hg19)
        │
        ▼  ② 在 ±2 bp 范围内匹配 Illumina 450K cg 探针（GPL13534 探针坐标）
cg 探针 ID（如 cg07881041）
        │
        ├─► ③ 查 GSE40279 全血 β 矩阵（656 样本均值）→ blood_mean_beta
        └─► ④ 查 TCGA-OV β 矩阵（~590 样本均值）→ tcga_mean_beta
```

**血液基线（GSE40279）**：cfDNA 来自血浆，血细胞 cfDNA 是主要背景噪声。若某 panel 位点在正常全血中已高甲基化（`blood_mean_beta ≥ 0.10`），血浆中检测到的甲基化信号可能主要源于血细胞污染而非肿瘤，该位点的肿瘤特异性存疑——应考虑从 panel 移除或降权（`blood_pass = False` 行黄色高亮）。

**肿瘤信号（TCGA-OV）**：确认 panel 位点在卵巢癌组织中确实存在甲基化信号（`tcga_mean_beta` 较高），与组织来源的差异甲基化结果形成互证。

**为什么会出现"Probe match: 0%"？**  
坐标系不同：panel 位点为 hg38，450K 数组坐标为 hg19，需 liftOver 转换后才能匹配。liftOver 成功（Liftover success 显示 ~96%）但探针匹配为 0，说明 **`GPL13534_platform.txt` 未放入 `data/external/`**——探针坐标文件缺失时 probe_map 为空，所有位点无从比对。按下方说明下载后重跑即可。

在主流程（Step 19 report）完成后运行，检验 panel 位点的：

- **血液基线**（正常全血 cfDNA 背景甲基化 < 0.10）——基于 GSE40279（656 名全血样本）
- **肿瘤信号**（卵巢癌组织中存在甲基化信号）——基于 TCGA-OV HumanMethylation450（~590 TCGA 样本）

**前置文件**（放入 `data/external/`）：

| 文件 | 来源 | 说明 |
| --- | --- | --- |
| `GSE40279_series_matrix.txt.gz` | [GEO GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) → Download full table | 全血 Illumina 450K β 矩阵，1.2 GB |
| `HumanMethylation450.gz` | [Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Ovarian%20Cancer%20(OV)) → DNA methylation → HumanMethylation450 → 下载 | TCGA-OV 450K β 矩阵，完整版 ~300 MB，注意要点击 Download 而非预览 |
| `hg38ToHg19.over.chain.gz` | [UCSC downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz) | 坐标转换链文件 |
| `GPL13534_platform.txt` | [NCBI GEO GPL13534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534) → Data table → Download full table… → 另存为 `GPL13534_platform.txt` | Illumina 450K 探针 hg19 坐标（推荐）；脚本也接受 `HM450.hg19.manifest.tsv.gz`（Zhou lab，~10 MB，自动尝试下载）；**若文件缺失，探针匹配列全为 False，不中止流程** |

**运行**：

```bash
# result-dir 填写本次分析的时间戳目录
run-analysis --plasma --step external_validate --workdir result/20260406
# 等价于（直接调用脚本）:
# python src/scripts/external_validate.py --result-dir result/20260406 --data-dir data
# 输出: result/20260406/external_validation.csv，并自动重新生成 HTML 报告（嵌入 Section 13）
```

**结果列说明**：

| 列 | 含义 |
| --- | --- |
| `cg_probe` | 匹配到的 Illumina 450K probe ID（liftOver hg38→hg19，±2 bp） |
| `probe_matched` | 是否成功匹配到 cg probe |
| `blood_mean_beta` | 全血 656 样本平均 β 值 |
| `blood_pass` | `blood_mean_beta < 0.10`（True = 满足血液低甲基化条件） |
| `tcga_mean_beta` | TCGA-OV 均值 β |
| `tcga_n_samples` | TCGA-OV 有覆盖的样本数 |

**解读**：`blood_pass = False` 的位点在血浆中可能受血细胞干扰，建议从 panel 中移除或降权。

---

## 注意事项

- `pysam` 仅支持 Linux/WSL，Windows 本机无法安装
- XGBoost / LightGBM 为可选依赖，未安装时自动跳过，其余分类器正常运行
- 步骤间状态通过 `temp/<workdir>/pipeline_state.json` 传递，重跑单步时不需要重新指定输入文件
- dedup / target 步骤均支持幂等重跑：输出文件已存在时自动跳过
- **M-bias 调参**：修改 `trim_tail2` 需从 `preprocess` 重跑；修改 `nOT`/`nOB` 只需从 `extract` 重跑
- **Outlier exclusion**: Step 14.5 only generates the report — no automatic exclusion. Add confirmed IDs to `exclude_samples` and re-run `classify` (no need to re-run feature selection)
- 88 个候选区域（`OC_tissue_88_regions.hg38.txt`）现为**软先验**：该区域内位点使用更宽松的 Wilcoxon 阈值（p<0.10, |Δβ|≥0.05），不再作为硬过滤条件
- **Annotation (Step 17)**: runs before `report` so the annotation table is embedded in the HTML. Requires internet access on first run; all API responses (Ensembl, EWAS Catalog, EWAS Atlas, PubMed) are cached in `annotation_cache.json` — re-runs are instant. If `annotate` is skipped, the report section is simply omitted.
- **Annotation — NCBI API key**: set `[annotation] ncbi_api_key` to increase PubMed rate limit from 3 to 10 req/s. Without a key the search still works but is slower. Register free at [NCBI account portal](https://www.ncbi.nlm.nih.gov/account/).
- **Coverage QC (Step 14)**: at the end of `feature_select`, sites where fewer than `min_10x_coverage_rate` (default 80%) of samples have ≥10x depth are hard-removed from `panel_cpg.csv`. Details recorded in `panel_cpg_depth.csv` (`n_10x`, `coverage_rate_10x`, `pass_coverage_qc` columns). Re-run `classify` after adjusting the threshold — no need to re-run feature selection.
- **Public dataset cross-validation (Step 18)**: run `--step external_validate` after `annotate`, then run `--step report` last to embed Section 13 (blood baseline + TCGA-OV table) in the HTML. Sites with `blood_pass = False` (mean β ≥ 0.10 in whole blood) are highlighted yellow — consider removing or down-weighting them from the panel.
- **Relation to tissue analysis**: the 88 candidate regions derived from tissue DMR analysis serve as a soft prior — they cannot be directly applied to plasma. Plasma cfDNA has lower signal-to-noise; statistical validation in plasma data is always required before modelling
