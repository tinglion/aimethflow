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
| `result/` | `feature_selection_freq.csv` | Per-fold frequency of the **leakage-free in-fold** feature selection (Wilcoxon + LASSO ∪ RF), run inside Step 15's CV over the **full candidate matrix**, not over the committed panel. Folds follow `[ml] cv_folds` (5-fold stratified by default; LOO when `cv_folds <= 0`) and each fold contributes `max_features_clf` rows, so `selection_count` maxes out at the fold count. It measures the **stability of the selection procedure**, not the stability of `panel_cpg.csv` |
| `result/` | `outlier_report.csv` | Outlier report (coverage / PCA z / Isolation Forest) |
| `result/` | `classifier_summary.csv` | AUC / Acc / Sens / Spec for **up to 11 models** when CA-125 data present: 5 cfDNA-only + CA-125 baseline + 5 `*+CA125` combined models. Includes `n_used` column showing how many samples contributed to each row (full cohort for cfDNA-only, CA-125-available subset for CA-125 / combined). |
| `result/` | `feature_selection_benchmark.csv` | Feature selection ablation: AUC / Acc / Sens / Spec for wilcoxon / lasso / rf / union × 5 classifiers |
| `result/` | `ca125_comparison.csv` | (Step 15c verification) CA-125 vs classifier AUC computed via Excel-direct lookup; cross-checks the in-loop CA-125 path |
| `result/` | `ca125_values.csv` | (Step 15c) Per-sample CA-125 value + true label for the matched subset |
| `result/` | `ca125_auc_inputs.csv` | (Step 15c debug cache) Full join trace for the in-loop CA-125 evaluation: sample_id → case_id → ca125 + final y/score arrays passed to `roc_auc_score` |
| `result/` | `classifier_predictions.csv` | Per-sample predictions; one row per (sample_id × model). For CA-125 / `*+CA125`, only samples with a CA-125 measurement are included |
| `result/` | `leakage_demo.csv` | (Step 15b) AUC / CI / Acc / Sens / Spec / features-per-fold for each protocol arm × classifier. `arm` ∈ `in_fold` \| `naive_selection` \| `naive_full` \| `panel_then_in_fold` \| `ca125_only`. With CA-125 present each protocol arm also carries the five `<classifier>+CA125` fusion models; `ca125_only` holds the single CA-125-alone row, which is protocol-invariant |
| `result/` | `leakage_demo_delong.csv` | (Step 15b) Paired DeLong test of every arm against the leakage-free `in_fold` reference: `delta_auc`, `delta_ci_low` / `delta_ci_high`, `delong_z`, `delong_p`. The CI is a paired bootstrap on the difference itself — both AUCs are recomputed on each resample — so it is not obtainable by comparing the two arms' separate CIs, which overlap even when a paired difference is large. Read it with the p-value: p says whether the difference is resolvable, the interval says how large a difference the data still permit, which is what turns "not significant" into "not large". `ca125_only` is excluded — it could only be compared against itself |
| `result/` | `leakage_demo_ca125_gain.csv` | (Step 15b, when CA-125 is available) Within each arm, what each half of a fusion model contributes. One row per arm × `<classifier>+CA125` × baseline, with `comparison` ∈ `vs_cfdna` (against its own cfDNA-only classifier — what CA-125 adds) \| `vs_ca125` (against CA-125 alone — what the cfDNA panel adds over the existing standard of care). Both comparisons sit inside a single arm, so neither is affected by the panel-provenance question. `delta_pct` is the relative increase over the baseline AUC; `headroom_pct` is the share of the baseline's remaining distance to AUC 1.0 that the fusion closes — the fairer figure, since AUC is floored at 0.5 and a relative percentage flatters a weak baseline. Same paired bootstrap CI and DeLong p as above |
| `result/` | `leakage_demo.png` | (Step 15b) Per-arm AUC comparison figure, cfDNA-only classifiers; embedded in the report's Step 15b section |
| `result/` | `leakage_demo_ca125.png` | (Step 15b, when CA-125 is available) The same comparison for the `+CA125` fusion models, as a separate figure so neither is too wide to read |
| `result/` | `roc_curve.png` | All-models ROC overlay (cfDNA solid, CA-125 dashed grey, `*+CA125` dotted same colour as base) |
| `result/` | `confusion_matrix.png` | Confusion matrices (auto-wraps to 4 rows × 3 cols when 11 models are evaluated) |
| `result/` | `temporal_validation_summary.csv` | Temporal split validation metrics (plasma only): AUC / Acc / Sens / Spec per classifier on samples collected on/after 2024-01-01. Includes 5 cfDNA-only + CA-125 baseline + 5 `*+CA125` combined models when CA-125 data is available; `n_validation` reflects what each row was evaluated on (full validation cohort for cfDNA / combined; matched subset for CA-125 alone). |
| `result/` | `temporal_validation_predictions.csv` | Per-sample predictions for temporal validation set (plasma only); one row per (sample × model). CA-125 / `*+CA125` rows omit samples without CA-125 measurements. |
| `result/` | `subtype_summary.csv` | Subtype classifier RF + SVM summary |
| `result/` | `subtype_predictions.csv` | Per-sample subtype predictions (with confidence) |
| `result/` | `subtype_roc.png` | Subtype OvR ROC curves |
| `result/` | `subtype_temporal_summary.csv` | Temporal split subtype validation metrics (plasma only): accuracy per classifier |
| `result/` | `subtype_temporal_predictions.csv` | Per-sample subtype predictions for temporal validation set (plasma only) |
| `result/` | `panel_annotation.csv` | Per-site annotation: gene / CpG island / EWAS Catalog / EWAS Atlas / `has_plasma_evidence` / PubMed hits |
| `result/` | `annotation_cache.json` | Cached API responses (speeds up re-runs) |
| `result/` | `external_validation.csv` | P3 cross-validation: per-site blood baseline + TCGA-OV β stats (`blood_mean_beta`, `blood_pass`, `tcga_mean_beta`, `tcga_n_samples`) |
| `result/` | `external_cohort_summary.csv` | (Step 19b) One row per cohort × arm × classifier, written into the **main** result dir. `arm` ∈ `panel_restricted` (the internal cohort's own rows, prepended) \| `full_matrix` (the primary external result) \| `denovo_panel` \| any sensitivity arm named in `[[external_cohort.arms]]`, e.g. `hgsoc_matched`, `with_healthy_donors`. Carries `delta_auc_vs_full` + `delong_p` for the de-novo arm |
| `result/external/<cohort>/` | `feature_matrix.npz`, `classifier_summary.csv`, `classifier_predictions.csv`, … | (Step 19b) Each external cohort's own result dir, prepared by `prepare_external_matrix.py` and then populated by the unmodified Step 15 |
| `result/external/<cohort>/` | `panel_discovery/` | (Step 19b, `discovery = true`) The de-novo arm: a copy of the cohort matrix with Step 14 selection applied, plus that arm's own Step 15 outputs |
| `result/<wd>/panel_concordance/` | `panel_concordance_sites.csv` | (`run_panel_concordance.py`, standalone) One row per committed panel site: its `region` cluster id, hg19 liftover position, matched EPIC probe and `match_distance_bp`, the plasma group means and signed `plasma_delta`, the external cohort's `tissue_delta` and p, and whether the two directions agree. `probe_w<N>` columns record the match at each window so a reader can re-cut the analysis at their own threshold |
| `result/<wd>/panel_concordance/` | `panel_concordance_summary.csv` | (`run_panel_concordance.py`) One row per matching window × counting `level` ∈ `site` (pseudo-replicated) \| `probe` (distinct EPIC probes) \| `region` (independent clusters — **the level to quote**). Carries `n_concordant` / `pct_concordant`, `pct_expected_null` measured from random EPIC probes rather than assumed to be 50%, and both `binom_p_vs_chance` and `binom_p_vs_null` |
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
