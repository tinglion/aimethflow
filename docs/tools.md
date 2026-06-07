# Helper Tools

Auxiliary scripts and one-liners that complement the main `run-analysis` pipeline.

## FastQ inspection

```bash
# Quick FastQ statistics
fq-inspect /mnt/f/meth/OC_plasma/B10001.R1.fq.gz
```

## Panel BED derivation

```bash
# Derive panel BED from FASTQ files (Step 3, preferred)
python src/scripts/derive_panel_bed.py \
    --fq /mnt/f/meth/OC_plasma/*.R1.fq.gz \
    --out data/panel_B10308.bed --samples 5

# Derive panel BED from BAM (manual one-liner, if you already have aligned BAMs)
bedtools genomecov -ibam sample.dedup.bam -bg \
    | awk '$4>10' | bedtools merge -d 100 \
    | awk 'length($1)<6' > panel.bed
```

## Intermediate file cleanup

```bash
# Clean up large intermediate files after analysis
rm $WORKDIR/*.bwameth_pe.sorted.bam $WORKDIR/*.bwameth_pe.sorted.bam.bai
rm $WORKDIR/*.bwameth_pe.dedup.bam  $WORKDIR/*.bwameth_pe.dedup.bam.bai
```

## UMI rescue audits

See [methodology.md § Quantitative validation](methodology.md#quantitative-validation--two-complementary-audits) for the methodology behind these audits.

### BAM-wide per-sample rescue rate

Parses `umi_tools dedup` stats logs, writes `docs/figures/umi_rescue_summary.csv` AND a two-panel diagnostic PNG in one shot:

```bash
python src/scripts/quantify_umi_rescue.py temp/<workdir>            # extract + plot
python src/scripts/quantify_umi_rescue.py temp/<workdir> --no-plot   # CSV only
python src/scripts/quantify_umi_rescue.py --plot-only                # re-plot from existing CSV
# Optional overrides: --out-csv <path>  --out-png <path>  --csv <existing-csv>
```

The diagnostic PNG has two panels:

| Panel | Content |
|---|---|
| (a) | Per-sample rescue % vs sequencing depth (log x-axis), coloured by mean unique UMIs per genomic position. Horizontal dashed line marks the cohort-wide median rescue; corner stats box reports n, median, and range. |
| (b) | Per-sample scatter of UMI-aware vs coord-only duplicate rate, with the y = x identity line. Distance below the diagonal = rescue magnitude in percentage points. Scales gracefully to any cohort size (no per-sample bars to crowd at 100+ samples). |

### Per-(sample, panel CpG) local rescue rate

Walks each UMI-deduped BAM with `pysam` to compute the *local* rescue rate at every final panel CpG position, then checks whether each site would still pass the ≥10× ≥80%-sample threshold under simulated coord-only dedup.

```bash
python src/scripts/quantify_panel_umi_rescue_bam.py \
    --result-dir result/<workdir> \
    --bam-dir    /path/to/<workdir>          \
    --config     config.toml \
    --threads    16
```

Runtime: ~1–3 min for a typical 118-sample cohort with 16 threads.

Outputs:

- `<result-dir>/umi_rescue_bam_attribution/per_sample_per_cpg_rescue.csv` — raw per-cell rescue table
- `<result-dir>/umi_rescue_bam_attribution/per_cpg_summary.csv` — aggregated per-CpG summary
- `<result-dir>/umi_rescue_bam_attribution/SUMMARY.md` — human-readable summary with rescue-dependent CpGs called out

## Unit tests

```bash
pytest src/test/ -v
```
