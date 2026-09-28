# Analysis Pipeline

The pipeline has two modes (plasma is the default). Each step is invocable via `--step <name>` for incremental runs, or chained automatically by the one-line `run-analysis --mode plasma <input>` entry point. See [quickstart.md](quickstart.md) for invocation examples.

## Plasma mode (`--mode plasma`, default)

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
| 13 | `outlier_detect` | Sample outlier detection (coverage / PCA z-score / Isolation Forest); suspect samples excluded from feature analysis | scikit-learn |
| 14 | `feature_analysis` | Feature analysis: LASSO + RF union selection; outlier samples excluded from computation; hard-filter sites with <80% samples at ≥10x depth | scikit-learn |
| 14a | `blood_filter` | **Optional.** Blood-baseline filter: drop panel sites with GSE40279 whole-blood mean β ≥ 0.10 (leukocyte-cfDNA confounders). Writes `panel_cpg_blood_filtered.csv`, which Steps 15/16/17 auto-prefer. Requires `data/external/{GSE40279_series_matrix.txt.gz, GPL13534_platform.txt, hg38ToHg19.over.chain.gz}`. | numpy / scipy |
| 14b | `fs_benchmark` | **Optional.** Feature selection ablation: compare Wilcoxon-only / LASSO-only / RF-only / LASSO∪RF AUC on same CV | scikit-learn |
| 15 | `classify` | Diagnostic classifiers, nested CV. Auto-evaluates **11 models** in the same CV loop: 5 cfDNA-only (LR/SVM/RF/XGB/LGBM), CA-125 alone (clinical baseline, ≥35 U/mL cutoff), and 5 combined cfDNA+CA-125 fusion models (`LR+CA125` … `LGBM+CA125`). Disable the CA-125 paths with `ml.include_ca125_baseline = false`. Auto-prefers `panel_cpg_blood_filtered.csv` (Step 14a output) when available. | LR / SVM / RF / XGB / LGBM |
| 15b | `leakage_demo` | **Optional, analysis-only.** Quantifies the feature-selection leakage this pipeline is built to avoid: re-runs the same classifiers over the **same CV partitions** under four protocols that differ only in what selection / imputation / scaling may see (`in_fold`, `naive_selection`, `naive_full`, `panel_then_in_fold`), then compares each against the leakage-free reference with the paired DeLong test. When the cohort has CA-125, every arm also runs the `<classifier>+CA125` fusion models and CA-125 alone is reported once as the protocol-invariant `ca125_only` row. Every ΔAUC carries a paired bootstrap CI beside its p-value, so a non-significant result reads as a bound on the effect rather than as its absence. With CA-125 present it also writes a gain table quantifying what each half of a fusion model contributes — against its own cfDNA classifier and against CA-125 alone. Writes `leakage_demo.csv`, `leakage_demo_delong.csv`, `leakage_demo.png` (+ `leakage_demo_ca125_gain.csv`, `leakage_demo_ca125.png`); rendered as its own HTML report section. Bootstrap resamples via `ml.leakage_n_boot`. | scikit-learn / scipy |
| 15c | `final_model` | **Optional.** Train production models on the full cohort (no CV holdout); saves artefacts to `result/<wd>/model/` for inference and the WebUI diagnosis page. | scikit-learn |
| 16 | `subtype_classify` | Subtype multi-class LOO-CV (malignant + borderline) | RF + SVM OvR |
| 17 | `temporal_validate` | **Plasma only.** Temporal split validation: train on samples collected before 2024-01-01, evaluate on samples collected on/after that date. Same 11-model lineup as Step 15 (5 cfDNA-only + CA-125 baseline + 5 cfDNA+CA-125 fusion) so head-to-head comparison carries through to the prospective evaluation. Covers binary + subtype classification. Requires `data/plasma_exp.xlsx` (columns: `case_id`, `datetime`). | scikit-learn / XGB / LGBM |
| 18 | `annotate` | Panel CpG annotation: gene context / CpG island / EWAS Catalog / EWAS Atlas / PubMed plasma hits | Ensembl REST + UCSC + EWAS Catalog + CNCB NGDC + NCBI |
| 19 | `external_validate` *(optional)* | Cross-validate panel sites against public 450K datasets (blood baseline + TCGA-OV) — details below | numpy / scipy |
| 19b | `external_cohort` | **Optional.** Independent external cohort re-run: executes the production Step 15 **unmodified** in each result directory listed in `config.toml [external_cohort] dirs` (each prepared by `src/scripts/prepare_external_matrix.py`), optionally deriving a panel de novo per cohort, plus any sensitivity arms declared as `[[external_cohort.arms]]` (a second build of the same deposit under a different sample selection, e.g. `hgsoc_matched` / `with_healthy_donors` for OvaPrint), and writes `external_cohort_summary.csv` back into the main result directory for a side-by-side internal/external HTML table. Validates the **framework**, not the panel — no public cohort shares this assay's probe set. Empty `dirs` ⇒ no-op. | scikit-learn |
| 20 | `report` | HTML + JSON summary report. Sections for Steps 18, 19, 15b and 19b are auto-embedded when their outputs exist, and show a "not yet run" placeholder otherwise | Custom |

## Tissue mode (`--mode tissue`)

Steps 4–12 available (bismark + deduplicate_bismark); steps 13–19 require the sample info Excel with benign/malignant groups.

## Step 19 — Public Dataset Cross-Validation (optional)

### Validation logic

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

### Required files (place in `data/external/`)

| File | Source | Description |
| --- | --- | --- |
| `GSE40279_series_matrix.txt.gz` | [GEO GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) → Download full table | Whole-blood Illumina 450K β matrix, 1.2 GB |
| `HumanMethylation450.gz` | [Xena Browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Ovarian%20Cancer%20(OV)) → DNA methylation → HumanMethylation450 → Download | TCGA-OV 450K β matrix, full version ~300 MB; click Download rather than preview |
| `hg38ToHg19.over.chain.gz` | [UCSC downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz) | Coordinate conversion chain file |
| `GPL13534_platform.txt` | [NCBI GEO GPL13534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534) → Data table → Download full table… → save as `GPL13534_platform.txt` | Illumina 450K probe hg19 coordinates (recommended); the script also accepts `HM450.hg19.manifest.tsv.gz` (Zhou lab, ~10 MB, auto-download attempted); **if the file is missing, the probe match column will be all False but the pipeline will not abort** |

### Run

```bash
# result-dir should be the timestamp directory for this analysis run
run-analysis --mode plasma --step external_validate --workdir result/20260406
# Equivalent (calling the script directly):
# python src/scripts/external_validate.py --result-dir result/20260406 --data-dir data
# Output: result/20260406/external_validation.csv, and automatically regenerates the HTML report (embeds Section 13)
```

### Output column descriptions

| Column | Meaning |
| --- | --- |
| `cg_probe` | Matched Illumina 450K probe ID (liftOver hg38→hg19, ±2 bp) |
| `probe_matched` | Whether a cg probe was successfully matched |
| `blood_mean_beta` | Mean β value across 656 whole-blood samples |
| `blood_pass` | `blood_mean_beta < 0.10` (True = meets low blood methylation criterion) |
| `tcga_mean_beta` | TCGA-OV mean β |
| `tcga_n_samples` | Number of TCGA-OV samples with coverage |

**Interpretation**: Sites with `blood_pass = False` may be subject to blood cell interference in plasma — consider removing or down-weighting them from the panel.
