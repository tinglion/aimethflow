# Configuration Reference (`config.toml`)

The pipeline reads its settings from a TOML file at the project root by default (`config.toml`). The Docker image instead uses `config.docker.toml` with container-internal paths.

Override the location with `--config <path>`.

## Key paths

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

## Plasma parameters `[plasma]`

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

### M-bias re-run scope

| Scenario | Adjust | Steps to re-run |
| --- | --- | --- |
| `trim_tail2` too small | Increase `trim_tail2` in `[plasma]` | `preprocess` → `dedup` → `target` → `extract` |
| Positional bias visible | Adjust `nOT` / `nOB` in `[plasma]` | `extract` onwards (no re-alignment) |

## Intermediate file cleanup `[cleanup]`

```toml
delete_after_dedup  = false   # delete sorted.bam after dedup
delete_after_target = false   # delete dedup.bam after on-target filter
```

## Thread configuration `[threads]`

```toml
fastp           = 8
bwameth         = 8
samtools_sort   = 4
methyldackel    = 4
```

## Machine learning `[ml]`

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

# Coverage QC hard filter (applied at end of feature_analysis)
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

## Feature Annotation `[annotation]`

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
