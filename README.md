# AIMethFlow — Ovarian Cancer Plasma cfDNA Targeted Methylation Analysis

End-to-end targeted bisulfite sequencing pipeline for ovarian-cancer early-diagnosis research, covering raw FASTQ → on-target methylation extraction → leakage-free machine-learning classification → temporal prospective validation → multi-source biological annotation, with a Docker-containerized deployment and a built-in Web UI.

## Status — pre-publication

The 3 modules `feature_utils` (5-stage feature selection), `classifier_utils` (11-model nested CV with cfDNA + CA-125 fusion), and `external_validate` (GSE40279 blood-baseline filter) are protected pending publication of the AIMethFlow paper in *iMeta* (2026) and the corresponding Chinese patent application. They are distributed as a private companion package `aimethflow-core`.

- **Reviewers** receive access to `aimethflow-core` via the corresponding author during peer review (see `LICENSE`).
- **All other pipeline steps (1-13, 16-18, 20)** run with the public code only.
- **Steps 14, 15, and 19** require `pip install aimethflow-core` (private). Without it, importing those modules raises an `ImportError` with the access instructions.
- On paper acceptance, `aimethflow-core` will be merged back into this repository under MIT.

## Analysis goals

- Detect differential methylation at CpG sites associated with ovarian cancer in plasma cell-free DNA (cfDNA)
- Data-driven discovery of candidate sites (individual-sample bedGraph direct aggregation; 88 candidate regions used as soft prior weighting)
- Train benign/malignant binary classification models (LR / SVM / RF / XGB / LGBM, Nested LOO-CV evaluation, feature selection performed within each fold's training set)
- Subtype multi-classification of malignant + borderline samples (high-grade serous carcinoma / endometrioid carcinoma / clear cell carcinoma / low-grade serous carcinoma / mucinous)

## Key features

- **Plasma-aware**: 6 nt UMI extraction (fastp) + UMI-aware deduplication (umi_tools), preserving genuine ~167 bp nucleosome-protected cfDNA fragments that coordinate-only dedup would discard
- **Leakage-free**: feature selection, imputation, and standardization confined strictly inside each CV training fold
- **Multimarker fusion**: 11-model benchmark in a single CV loop — 5 cfDNA classifiers + CA-125 baseline + 5 cfDNA + CA-125 fusion models
- **Two-stage signal preservation**: upstream UMI rescue + downstream GSE40279 whole-blood baseline filter remove leukocyte-cfDNA confounders before classification
- **Temporal prospective validation**: chronological held-out split confirms internal-CV gains translate to chronologically independent cohorts
- **One-command Dockerized deployment** with embedded Web UI on port 9999

## Quick start (Docker)

```bash
cp .env.example .env                                       # set host paths
docker compose --progress=plain build                       # ~10-20 min, first time
docker compose up -d                                        # persistent container
docker compose exec analysis \
    python /app/src/scripts/run_analysis.py --steps 1-19
# Web UI: http://<server-ip>:9999
```

See [DOCKER.md](DOCKER.md) for the full Docker reference, or [docs/installation.md](docs/installation.md) for Conda / local-development setup.

## Documentation

| Topic | Where |
|---|---|
| Web UI walkthrough + step-by-step pipeline + shell-script template | [docs/quickstart.md](docs/quickstart.md) |
| 20-step pipeline reference (plasma + tissue) | [docs/pipeline.md](docs/pipeline.md) |
| `config.toml` field-by-field reference | [docs/configuration.md](docs/configuration.md) |
| Output file catalogue + sample Excel format | [docs/outputs.md](docs/outputs.md) |
| Library-preparation prerequisites + methodological improvements + UMI-rescue audits | [docs/methodology.md](docs/methodology.md) |
| Helper scripts (FastQ inspection, panel-BED derivation, UMI audits, cleanup) | [docs/tools.md](docs/tools.md) |
| Hardware sizing, peak memory, disk-usage measurements | [docs/resources.md](docs/resources.md) |
| Documentation index | [docs/README.md](docs/README.md) |

## Project layout

```text
├── config.toml                     # Default configuration
├── config.docker.toml              # Container-internal-path variant
├── pyproject.toml                  # Package + dependency manifest
├── DOCKER.md                       # Docker deployment reference
├── README.md                       # (this file)
├── data/                           # Reference data (BED, Excel, externals)
├── docs/                           # User documentation (see table above)
└── src/
    ├── methylation/                # Differential methylation analysis modules
    ├── scripts/
    │   ├── run_analysis.py         # Main CLI entry point
    │   ├── webui.py                # FastAPI Web UI (port 9999)
    │   ├── external_validate.py    # Step 19 public-dataset cross-validation
    │   ├── quantify_umi_rescue.py            # UMI-rescue BAM-wide audit + plot
    │   ├── quantify_panel_umi_rescue_bam.py  # Per-CpG local UMI-rescue audit
    │   └── modules/
    │       ├── alignment_utils.py           # Alignment / dedup / on-target / extraction
    │       ├── fastq_utils.py               # fastp preprocessing
    │       ├── feature_utils.py             # Feature analysis (Step 14) + outlier (Step 13)
    │       ├── classifier_utils.py          # Diagnostic classification (Step 15)
    │       ├── subtype_utils.py             # Subtype + temporal validation (Steps 16/17)
    │       ├── annotation_utils.py          # Site annotation (Step 18)
    │       └── report_utils.py              # Report generation (Step 20)
    ├── test/
    └── utils/                      # Config + logger helpers
```

## Environment

> **Linux / WSL2 only** — bioinformatics tools (bismark, bwameth, MethylDackel, etc.) do not support native Windows environments.

## License

See `LICENSE` for full terms.
