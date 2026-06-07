# Output Files and Sample Excel Format

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
| `result/` | `panel_cpg_blood_filtered.csv` | Step 14a output: `panel_cpg.csv` minus sites with GSE40279 mean β ≥ 0.10. Auto-preferred by Steps 15/16/17. |
| `result/` | `panel_cpg_blood_summary.csv` | Step 14a per-site blood baseline: `cg_probe`, `blood_mean_beta`, `blood_pass` |
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
