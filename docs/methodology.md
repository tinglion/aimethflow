# Methodology

## Library preparation prerequisites

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

## Methodological improvements

Key differences vs. the prior (Youxun) approach:

| Issue | Prior Approach | This Pipeline |
| --- | --- | --- |
| **Deduplication** | Picard coordinate dedup — removes genuine cfDNA fragments (147–167 bp) | **umi_tools dedup --paired**, UMI-aware |
| **UMI** | UMI built but ignored | fastp `--umi` extraction + umi_tools dedup |
| **On-target filter** | None | **bedtools intersect panel.bed** |
| **Poly-C tail** | R2 3' not trimmed | `trim_tail2=5`, confirmed by M-bias SVG |
| **Methylation threshold** | 80% (stricter than tissue) | **60%**, same as tissue |
| **Differential filter** | p<0.05 only | **p<0.05 AND \|Δβ\|>0.1** |

## Why UMI-aware deduplication matters for cfDNA

**The problem.** Plasma cfDNA is not randomly sheared — circulating fragments are protected by histone octamers and concentrate at the ~147 bp mono-nucleosome and ~167 bp mono-nucleosome + linker sizes. Many genuinely independent cfDNA molecules therefore share identical start/end coordinates because they were shaped by the same nucleosome footprint. Coordinate-only deduplication (Picard `MarkDuplicates`, `deduplicate_bismark`) cannot tell these apart from true PCR duplicates and collapses them all into a single read, silently discarding real molecular complexity. (Tissue DNA, fragmented by sonication or restriction enzymes, carries no such bias, so coordinate dedup is safe there.)

**How this pipeline handles it (plasma mode):**

| Step | Tool & flags | What it does |
| --- | --- | --- |
| 4 — preprocess ([fastq_utils.py:220-228](../src/scripts/modules/fastq_utils.py#L220-L228)) | `fastp --umi --umi_loc read2 --umi_len 6 --umi_prefix UMI` | Extracts the 6 N inline UMI from the 5' end of R2 and appends it to the read name (`@HEADER_UMI:NNNNNN`). Also applies `trim_front1=12 / trim_front2=12 / trim_tail2=5 / Q20 / min_length=35`. |
| 5 — align ([alignment_utils.py:147](../src/scripts/modules/alignment_utils.py#L147)) | `bwameth` paired-end | UMI rides through alignment as part of the read name into the BAM. |
| 6 — dedup ([alignment_utils.py:293-294](../src/scripts/modules/alignment_utils.py#L293-L294)) | `umi_tools dedup --paired` | Deduplication key is the **(start, end, UMI) triple** — fragments with identical coordinates but distinct UMIs are kept as separate molecules; only entries identical on all three are collapsed as PCR duplicates. |

**Tissue mode (`plasma=False`)** bypasses UMIs and uses `deduplicate_bismark` for coordinate-only dedup ([alignment_utils.py:315](../src/scripts/modules/alignment_utils.py#L315)). If `umi_tools` is missing or fails in plasma mode the code falls back to `deduplicate_bismark` and emits a `not recommended for plasma data` warning ([alignment_utils.py:305-306](../src/scripts/modules/alignment_utils.py#L305-L306)).

**Sanity-check telemetry.** The report surfaces the mean UMI deduplication rate ([report/main.py:330](../src/scripts/modules/report/main.py#L330)) and a fragment-size distribution figure that highlights the 100–220 bp cfDNA zone with the 167 bp mono-nucleosome line marked ([alignment_utils.py:754-756](../src/scripts/modules/alignment_utils.py#L754-L756)) — both serve as cross-validation that the input is genuinely cfDNA and that the dedup step has not flattened the nucleosome peak.

## Quantitative validation — two complementary audits

The pipeline ships two stand-alone scripts that quantify how much real signal UMI-aware dedup recovers (vs. simulated coordinate-only dedup) and how much of that recovery reaches the final diagnostic panel:

| Script | Granularity | What it measures | Reads from |
|---|---|---|---|
| [`quantify_umi_rescue.py`](../src/scripts/quantify_umi_rescue.py) | BAM-wide, per-sample | Parses `umi_tools dedup --output-stats` logs; computes per-library rescue% = `(UMI_kept − coord_only_would_keep) / UMI_kept`. Outputs `docs/figures/umi_rescue_summary.csv` + plot. | `temp/<wd>/*.dedup.stats` log files |
| [`quantify_panel_umi_rescue_bam.py`](../src/scripts/quantify_panel_umi_rescue_bam.py) | Per-(sample, panel CpG) | Walks each UMI-deduped BAM with `pysam`, groups reads at each panel CpG position by `(R1_start, R2_start)`, computes *local* per-CpG rescue rate, and checks whether each panel site would still pass the ≥10× ≥80%-sample threshold under simulated coord-only dedup. | `<bam-dir>/*.bwameth_pe.dedup.targeted.sorted.bam` + `panel_cpg_blood_filtered.csv` |

See [tools.md](tools.md#umi-rescue-audits) for the run commands.

**Reference result on the 118-sample ovarian-cancer cohort (`result/20260406`):**

| Quantity | Value |
|---|---|
| BAM-wide rescue rate (146 sequencing libraries, 3,744 M input read pairs) | median 20.83%, range 0.17–33.32%, up to 2,551 unique UMIs / single genomic position |
| Median *local* rescue rate across the 60 final panel CpGs | **9.20%** |
| Maximum local rescue at any (sample, CpG) cell | **50.0%** |
| Cancer-related panel genes with the highest median local rescue | **SOX1 (16.4%), C2CD4D (16.2%)** |
| Panel CpGs that would have failed the ≥10× ≥80% threshold under coord-only dedup | **0 of 60** |

**Interpretation.** Panel CpGs sit in regions with substantial depth headroom (median depth ~73× under UMI-aware dedup), so the ≥10× ≥80% hard filter is not the binding constraint of UMI rescue at panel positions — even sites with 50% local rescue retain ≥10× in ≥80% of samples after the simulated coord-only collapse. The rescued reads contribute primarily to β-estimate precision and the statistical stability of the LASSO / RF feature ranking, rather than to gating panel inclusion through the hard depth threshold.

## Operational notes

- `pysam` is Linux/WSL only; cannot be installed on native Windows
- XGBoost / LightGBM are optional dependencies; if not installed they are automatically skipped and the remaining classifiers run normally
- Inter-step state is passed via `temp/<workdir>/pipeline_state.json` — no need to re-specify input files when re-running a single step
- Both the dedup and target steps support idempotent re-runs: if output files already exist they are automatically skipped
- **M-bias auto-detection**: `--step mbias` writes `mbias_auto_trim.json` with converted skip counts. MethylDackel outputs *inclusion positions* (`--OT a,b,c,d`); the pipeline converts them: 5' skip = position directly; 3' skip = `read_len − pos − 1` (e.g. `--OT 0,0,0,110` + `read_len=150` → `nOT "0,0,0,39"`). Set `nOT = "auto"` / `nOB = "auto"` to use in `extract`; explicit skip-count strings override auto. Changing `trim_tail2` requires re-running from `preprocess`; changing `nOT`/`nOB` only requires re-running from `extract`.
- **Outlier exclusion**: Step 13 generates the report and saves suspect sample IDs to pipeline state. Step 14 (`feature_analysis`) automatically excludes these samples from LASSO/RF computation — the saved `feature_matrix.npz` still contains all samples so classifiers can apply their own exclusion via `exclude_samples`. Add confirmed outlier IDs to `exclude_samples` and re-run `classify`.
- **CA-125 head-to-head (Step 15)**: with a `ca125` column (numeric, U/mL) in the sample info Excel, the `classify` step automatically evaluates 11 models in the same CV folds — 5 cfDNA-only, **CA-125 alone** (clinical baseline; raw value as score, ≥35 U/mL cutoff for accuracy/sens/spec), and **5 combined cfDNA + CA-125 fusion models** (`LR+CA125`, `SVM+CA125`, `RF+CA125`, `XGB+CA125`, `LGBM+CA125`). Combined models concatenate CA-125 (training-fold-median imputed, separately standardized) as an extra column to the cfDNA feature matrix — leakage-free. All 11 ROC curves render together in `roc_curve.png` (cfDNA solid, CA-125 dashed grey, combined dotted same colour as base). Disable via `ml.include_ca125_baseline = false`; tune the binary cutoff via `ml.ca125_threshold = 35.0`. The optional Step 15c `ca125_compare` is now a verification step that cross-checks the in-loop CA-125 AUC against an Excel-direct computation — keep it as an audit but it's no longer the primary CA-125 path.
- **Standalone CA-125 AUC** (`src/scripts/compute_ca125_auc.py`): for quick CA-125-only AUC without running any pipeline step. Reads `category` + `ca125` columns directly from the Excel, applies `ml.exclude_samples` (override with `--exclude` flag), prints AUC + 95% bootstrap CI. Useful as a sanity check against Step 15's in-loop CA-125 AUC; the two should agree to 3 decimal places.
- **Stratified k-fold CV correctness**: prior to v1.0.0 the per-sample `classifier_predictions.csv` had a label-shuffling bug under `cv_folds > 0` (StratifiedKFold with shuffle): `true_label` and `prob_malignant` were attached to the wrong `sample_id` because predictions were appended in fold-iteration order while the CSV writer iterated `sample_ids` in original order. **Aggregate AUC/Acc/Sens/Spec metrics were unaffected**, but downstream per-sample analyses (CA-125 join, confusion-matrix breakdowns, temporal-validation warm-start) saw scrambled labels. Fixed by pre-allocating per-classifier arrays of length `n_samples` and writing back at `test_idx` positions. LOO-CV (`cv_folds=0`) was unaffected because LeaveOneOut iterates 0→n-1. **If you have a workdir from before this fix, re-run `--step classify`** to regenerate `classifier_predictions.csv` with correctly aligned rows.
- **Blood-baseline filter (Step 14a, optional)**: run `--step blood_filter` after `feature_analysis` to drop panel sites with GSE40279 (n = 656 whole-blood samples) mean β ≥ 0.10 — these are leukocyte-cfDNA confounders that cannot distinguish tumor-derived signal from blood-cell background. Writes `panel_cpg_blood_filtered.csv` and `panel_cpg_blood_summary.csv` to `result/<wd>/`. Steps 15 (`classify`), 16 (`subtype_classify`), and 17 (`temporal_validate`) automatically prefer the blood-filtered panel when present and silently fall back to `panel_cpg.csv` if absent. Requires `data/external/{GSE40279_series_matrix.txt.gz, GPL13534_platform.txt, hg38ToHg19.over.chain.gz}` — see [pipeline.md § Step 19](pipeline.md#step-19--public-dataset-cross-validation-optional) for download instructions. The same filter is also re-applied internally by Step 19 (`external_validate`) so manuscripts that don't formally run Step 14a still produce a filtered panel after Step 19.
- **Feature selection ablation (Step 14b, optional)**: run `--step fs_benchmark` after `classify` to compare four feature selection strategies (wilcoxon-only / lasso-only / rf-only / lasso∪rf) using the same CV, classifiers, and cohort. Results in `feature_selection_benchmark.csv`; automatically embedded in Step 15 of the HTML report. In the 116-sample ovarian cancer cohort, ML-based methods (best AUC 0.804–0.809) substantially outperformed wilcoxon-only top-Δβ ranking (best AUC 0.759). `feature_selection_mode = "wilcoxon"` is also available in `config.toml` for single-mode runs.
- **Temporal split validation (Step 17, plasma only)**: run `--step temporal_validate` after Steps 15 and 16 to evaluate the pre-selected panel on prospectively collected samples. Requires `data/plasma_exp.xlsx` with columns `case_id` and `datetime` (format `YYYYMMDD`). Samples with `datetime < 2024-07-01` form the train set; samples on/after that date form the validation set. Mirrors Step 15's **11-model lineup** when a `ca125` column is present: 5 cfDNA-only classifiers + CA-125 baseline (≥35 U/mL cutoff, evaluated only on samples with a measurement) + 5 `*+CA125` combined models (CA-125 imputed on the train set with the training median, separately standardized, then concatenated as an extra column). Combined models give a prospective-deployment view of the marginal value of fusing CA-125 with cfDNA. Toggle via `ml.include_ca125_baseline`; threshold via `ml.ca125_threshold`. If the date Excel is absent the step logs a warning and exits cleanly. Outputs: `temporal_validation_summary.csv`, `temporal_validation_predictions.csv`, `subtype_temporal_summary.csv`, `subtype_temporal_predictions.csv`.
- The 88 candidate regions (`OC_tissue_88_regions.hg38.txt`) are now a **soft prior**: sites within these regions use more relaxed Wilcoxon thresholds (p<0.10, |Δβ|≥0.05) and are no longer a hard filter condition
- **Annotation (Step 18)**: runs before `report` so the annotation table is embedded in the HTML. Requires internet access on first run; all API responses (Ensembl, EWAS Catalog, EWAS Atlas, PubMed) are cached in `annotation_cache.json` — re-runs are instant. If `annotate` is skipped, the report section is simply omitted.
- **Annotation — NCBI API key**: set `[annotation] ncbi_api_key` to increase PubMed rate limit from 3 to 10 req/s. Without a key the search still works but is slower. Register free at [NCBI account portal](https://www.ncbi.nlm.nih.gov/account/).
- **Coverage QC (Step 14)**: at the end of `feature_analysis`, sites where fewer than `min_10x_coverage_rate` (default 80%) of samples have ≥10x depth are hard-removed from `panel_cpg.csv`. Details recorded in `panel_cpg_depth.csv` (`n_10x`, `coverage_rate_10x`, `pass_coverage_qc` columns). Re-run `classify` after adjusting the threshold — no need to re-run feature analysis.
- **Public dataset cross-validation (Step 19)**: run `--step external_validate` after `annotate`, then run `--step report` last to embed Section 13 (blood baseline + TCGA-OV table) in the HTML. Sites with `blood_pass = False` (mean β ≥ 0.10 in whole blood) are highlighted yellow — consider removing or down-weighting them from the panel.
- **Relation to tissue analysis**: the 88 candidate regions derived from tissue DMR analysis serve as a soft prior — they cannot be directly applied to plasma. Plasma cfDNA has lower signal-to-noise; statistical validation in plasma data is always required before modelling
