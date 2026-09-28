# Configuration Reference (`config.toml`)

The pipeline reads its settings from a TOML file at the project root by default (`config.toml`). The Docker image instead uses `config.docker.toml` with container-internal paths.

Override the location with `--config <path>`.

## Key paths

```toml
[data]
plasma_dir         = "/mnt/f/meth/OC_plasma"
sample_info_excel  = "ovarian_fq_series_samples_v4.xlsx"   # resolved under data_dir
plasma_exp_excel   = "data/plasma_exp.xlsx"                # Step 17: case_id + datetime
probe_panel_bed    = "panel_B10308.bed"
candidate_regions  = "OC_tissue_88_regions.hg38.txt"

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
sample_wilcox_p        = 0.05      # p threshold (non-88 sites)
sample_delta_beta      = 0.10      # |Δβ| threshold (non-88 sites)
sample_wilcox_p_88     = 0.05      # p threshold (88 candidate region sites, relaxed)
sample_delta_beta_88   = 0.10      # |Δβ| threshold (88 candidate region sites, relaxed)

# Overfitting control
max_features_clf  = 20             # max features entering classifiers (≈ n_samples/10)
cv_folds          = 5              # Step 15 outer loop; 0 = leave-one-out CV
lasso_cv_folds    = 5              # inner loop within each outer training fold

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

# Step 15b (leakage_demo): bootstrap resamples per arm for the AUC confidence interval.
# Lower to ~200 for a quick look; 1000 is what the published tables use.
leakage_n_boot    = 1000

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

## Independent external cohorts `[external_cohort]` (Step 19b)

```toml
[external_cohort]
# Result dirs prepared by src/scripts/prepare_external_matrix.py; each must contain
# feature_matrix.npz. An empty list makes Step 19b a clean no-op.
dirs      = []      # e.g. ["result/external/ovaprint"]
labels    = []      # display names, positional; falls back to the directory basename
discovery = false   # also derive a panel de novo per cohort (Step 14) and report the
                    # paired DeLong comparison against that cohort's full-matrix run

# Optional sensitivity arms: another build of the SAME deposit under a different
# sample selection, run through Step 15 and reported beside the primary arm.
# `cohort` must match the cohort's display label.
[[external_cohort.arms]]
cohort = "OvaPrint"
arm    = "hgsoc_matched"          # prepare_external_matrix.py --subset hgsoc_matched
dir    = "result/external/ovaprint_hgsoc"
```

Two arms are defined for the OvaPrint deposit. `hgsoc_matched` restricts it to the comparison
the original authors reported (HGSOC versus benign masses), so the reproduction can be judged
against their published figures; `with_healthy_donors` (built with `--include-normal`) folds
the deposit's healthy donors into the benign class. Both are sensitivity arms — the primary
arm is the all-comers question, and it is the only one comparable to the internal cohort.

Step 19b runs the **unmodified** production Step 15 inside each listed directory and writes
`external_cohort_summary.csv` back into the main result directory, so the HTML report renders
an internal-vs-external table side by side. Three internal-cohort bindings are neutralised
automatically: `include_ca125_baseline` (CA-125 is keyed by internal sample IDs, so no
external sample resolves to a value), `exclude_samples`
(those IDs are internal), and the panel filter (no `panel_cpg*.csv` exists externally, so
per-fold selection runs over the full matrix).

## Temporal validation `[temporal_validate]` (Step 17)

```toml
cutoff_date          = "2024-01-01"  # collected before → train; on/after → validation
reselect_features    = true          # re-run selection on training-period samples only
                                     # (false = reuse feature_selection_freq.csv — leaks)
top_n_features       = 20            # used only when reselect_features = false
borderline_as_benign = false         # add validation-period borderline samples as y=0
low_ca125_threshold  = 105.0         # extra evaluation on the CA-125 < 105 subset;
                                     # 0 or negative disables it
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
