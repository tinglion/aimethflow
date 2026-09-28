#!/usr/bin/env bash
# download_external_data.sh — fetch external reference data.
#
#   Part A (~1.4 GB) — the 4 files required by Step 19 (external_validate).
#   Part B (~1.1 GB) — external validation cohorts (OvaPrint ~4.6 MB, GSE19711
#                      ~58 MB, GSE133556 ~755 MB + its EPIC manifest ~272 MB).
#                      Optional: nothing in the main pipeline needs them, but the
#                      GPB external-cohort work does.
#                      See docs/GPB/02_external_validation_plan.md.
#
# Usage:   bash download_external_data.sh [DEST_DIR]
# Default: DEST_DIR = data/external (relative to the current working directory)
#
# Behaviour:
#   - Creates DEST_DIR (and any per-cohort subdirectory) if absent.
#   - Skips any file already present with non-zero size.
#   - Resumes a partial download with curl -C -.
#   - Prints a per-file OK / FAIL status. Exits non-zero only if a *required*
#     file failed; a failed optional cohort is reported but does not fail the run.
#
# See docs/pipeline.md → "Required files (place in data/external/)" for the
# full source table and troubleshooting notes.

set -uo pipefail

DEST="${1:-data/external}"
mkdir -p "$DEST"

echo "Downloading external reference files to $DEST/ (~1.4 GB total)"
echo "  (existing non-zero files are skipped)"
echo

FAILED=0
OPT_FAILED=0
# Set by _check_or_get: which counter a failure increments.
OPTIONAL=0

_check_or_get() {
    local name="$1" url="$2" note="${3:-}"
    local target="$DEST/$name"

    if [[ -s "$target" ]]; then
        printf "  [skip]  %-40s  (%s already present)\n" "$name" "$(du -h "$target" | cut -f1)"
        return 0
    fi

    # $name may carry a subdirectory (e.g. ovaprint/metadata.csv.gz).
    mkdir -p "$(dirname "$target")"

    printf "  [get ]  %s\n          from  %s\n" "$name" "$url"
    if curl -L --fail --retry 3 --continue-at - -o "$target.part" "$url"; then
        mv "$target.part" "$target"
        printf "  [ok  ]  %s  (%s)\n\n" "$name" "$(du -h "$target" | cut -f1)"
        return 0
    else
        rm -f "$target.part"
        printf "  [FAIL]  %s\n" "$name"
        [[ -n "$note" ]] && printf "          %s\n" "$note"
        echo
        if [[ $OPTIONAL -eq 1 ]]; then
            OPT_FAILED=$((OPT_FAILED + 1))
        else
            FAILED=$((FAILED + 1))
        fi
        return 1
    fi
}

# ── 1. UCSC hg38 → hg19 liftOver chain (~1.2 MB) ──────────────────────────────
_check_or_get "hg38ToHg19.over.chain.gz" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"

# ── 2. GSE40279 whole-blood series matrix (~1.2 GB) ───────────────────────────
_check_or_get "GSE40279_series_matrix.txt.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz"

# ── 3. GPL13534 Illumina 450K probe hg19 coordinates (~213 MB) ─────────────────
_check_or_get "GPL13534_platform.txt" \
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&view=data&form=text&acc=GPL13534" \
    "GEO's text-download endpoint sometimes rate-limits or returns HTML. If the file above is empty or looks wrong, download manually: (1) open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534 in a browser, (2) at the bottom under 'Data table', click 'Download full table…', (3) save as GPL13534_platform.txt in $DEST/."

# ── 4. TCGA-OV HumanMethylation450 β matrix (~300 MB in full version) ─────────
_check_or_get "HumanMethylation450.gz" \
    "https://tcga.xenahubs.net/download/TCGA.OV.sampleMap/HumanMethylation450.gz" \
    "If the Xena URL fails, browse to https://xenabrowser.net/datapages/?cohort=TCGA%20Ovarian%20Cancer%20(OV) → DNA methylation → HumanMethylation450 → Download, save the .gz file to $DEST/."

# ══════════════════════════════════════════════════════════════════════════════
# Part B — external validation cohorts (optional, ~4.6 MB)
# ══════════════════════════════════════════════════════════════════════════════
echo "──────────────────────────────────────────────────────────────────────"
echo "Optional: external validation cohorts (GPB external-cohort work)"
echo
OPTIONAL=1

# ── OvaPrint plasma cfDNA (Buckley et al., Clin Cancer Res 2023;29(24):5196) ──
#    Zenodo 10.5281/zenodo.8395817, CC BY 4.0. 372 samples x 1677 features;
#    313 plasma cfDNA. Feature IDs are anonymised (FEATURE_n) — no coordinates.
_ZEN="https://zenodo.org/records/8395817/files"
_check_or_get "ovaprint/beta_value_matrix.csv.gz" \
    "$_ZEN/beta_value_matrix.csv.gz?download=1" \
    "Browse to https://doi.org/10.5281/zenodo.8395817 and download manually into $DEST/ovaprint/."
_check_or_get "ovaprint/metadata.csv.gz" \
    "$_ZEN/metadata.csv.gz?download=1" \
    "Browse to https://doi.org/10.5281/zenodo.8395817 and download manually into $DEST/ovaprint/."

if [[ -s "$DEST/ovaprint/beta_value_matrix.csv.gz" && -s "$DEST/ovaprint/metadata.csv.gz" ]]; then
    echo "  OvaPrint ready. Build the pipeline entry point with:"
    echo "     python src/scripts/prepare_external_matrix.py --dataset ovaprint \\"
    echo "         --input-dir $DEST/ovaprint --result-dir result/external/ovaprint"
    echo
fi

# ── GSE19711 / UKOPS whole blood (Teschendorff et al., Genome Res 2010) ───────
#    Illumina 27K, 540 samples: 274 cancer-free postmenopausal controls and 266
#    ovarian cancer patients (131 pre-treatment, 135 post-treatment). Sample
#    metadata is embedded in the series matrix's own !Sample_* header block, so
#    this single file is the whole cohort.
#
#    NOTE — whole blood, not plasma cfDNA, and no probe overlap with this
#    assay's panel. It is a framework-replication cohort only; see
#    docs/GPB/07_gse19711_results.md.
_check_or_get "GSE19711/GSE19711_series_matrix.txt.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19711/matrix/GSE19711_series_matrix.txt.gz"

if [[ -s "$DEST/GSE19711/GSE19711_series_matrix.txt.gz" ]]; then
    echo "  GSE19711 ready. Build the primary arm and its two sensitivity arms with:"
    echo "     python src/scripts/prepare_external_matrix.py --dataset gse19711 \\"
    echo "         --input-dir $DEST/GSE19711 --result-dir result/external/gse19711"
    echo "     python src/scripts/prepare_external_matrix.py --dataset gse19711 \\"
    echo "         --subset all_cases --input-dir $DEST/GSE19711 \\"
    echo "         --result-dir result/external/gse19711_all_cases"
    echo "     python src/scripts/prepare_external_matrix.py --dataset gse19711 \\"
    echo "         --subset balanced_batches --input-dir $DEST/GSE19711 \\"
    echo "         --result-dir result/external/gse19711_balanced"
    echo
fi

# ── GSE133556 / EPIC ovarian tissue (Buckley et al., Clin Cancer Res 2023) ───
#    Illumina EPIC (GPL21145), 111 samples: 99 primary high-grade serous ovarian
#    carcinomas and 12 normal fallopian tube controls. Sample labels live in the
#    series matrix's own !Sample_characteristics_ch1 "tissue:" line.
#
#    Used ONLY for region-level panel concordance — the one analysis in the study
#    that tests the panel rather than the framework. EPIC probes carry genomic
#    coordinates, so the committed panel's hg38 sites can be lifted to hg19 and
#    located here; no classifier is ever run on this cohort. See
#    docs/GPB/08_gse133556_panel_concordance.md.
#
#    Large: the series matrix is ~755 MB because betas are deposited at full
#    precision, and the EPIC manifest is ~272 MB. Both are read once and streamed.
_check_or_get "GSE133556/GSE133556_series_matrix.txt.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133556/matrix/GSE133556_series_matrix.txt.gz"

# GPL21145 EPIC probe hg19 coordinates. Same GEO text endpoint as GPL13534 above;
# the file carries Genome_Build = 37 for all 864,935 cg probes.
_check_or_get "GPL21145_platform.txt" \
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&view=data&form=text&acc=GPL21145" \
    "GEO's text-download endpoint sometimes rate-limits or returns HTML. If the file above is empty or looks wrong, download manually: (1) open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21145 in a browser, (2) at the bottom under 'Data table', click 'Download full table…', (3) save as GPL21145_platform.txt in $DEST/."

if [[ -s "$DEST/GSE133556/GSE133556_series_matrix.txt.gz" && -s "$DEST/GPL21145_platform.txt" ]]; then
    echo "  GSE133556 ready. Run the region-level panel concordance analysis with:"
    echo "     python src/scripts/run_panel_concordance.py \\"
    echo "         --result-dir result/20260406 --external-dir $DEST \\"
    echo "         --series-matrix $DEST/GSE133556/GSE133556_series_matrix.txt.gz"
    echo
fi

OPTIONAL=0

echo "──────────────────────────────────────────────────────────────────────"
if [[ $OPT_FAILED -ne 0 ]]; then
    echo "  $OPT_FAILED optional cohort file(s) failed — the main pipeline is unaffected."
fi
if [[ $FAILED -eq 0 ]]; then
    echo "All 4 required files present in $DEST/:"
    ls -lh "$DEST"
    exit 0
else
    echo "  $FAILED required file(s) failed to download automatically."
    echo "  Follow the manual instructions printed above, or see:"
    echo "     docs/pipeline.md → Required files (place in data/external/)"
    exit 1
fi
