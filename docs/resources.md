# Resource Requirements

Reference values measured on an 8-core / 32 GB machine (N = 16 M read pairs):

| Step | Tool | Wall time | Peak memory | Output size |
| --- | --- | --- | --- | --- |
| preprocess | fastp | ~2 min | 0.5 GB | ×0.90 of raw FASTQ |
| align | bwameth + samtools | **~47 min** | **13–14 GB** | ≈ clean FASTQ × 0.85 |
| dedup | umi_tools | **~10 min** | 2–4 GB | ×0.38–0.50 of sorted BAM |
| target | bedtools | ~4 min | 1–2 GB | ×0.25–0.36 of dedup BAM |
| mbias | MethylDackel | ~5 min | 2–4 GB | ~2 MB (SVG) |
| extract | MethylDackel | ~10 min | 2–4 GB | ~20–50 MB |
| filter + diff + visualize | Python | <5 min | 1 GB | <5 MB |
| feature_analysis | scikit-learn | 5–30 min | 2–4 GB | <10 MB |
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
