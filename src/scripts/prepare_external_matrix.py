#!/usr/bin/env python3
"""
Build a pipeline-compatible ``feature_matrix.npz`` from a public methylation
beta matrix, so that an external cohort can enter the AIMethFlow analysis at
Step 13 without passing through the FASTQ-based Steps 1-12.

Why this exists
---------------
Every AIMethFlow entry point below Step 13 starts from raw reads. Public
methylation cohorts are distributed as processed beta matrices, so validating
the *framework* on an independent cohort requires a matrix-level ingress.
``feature_matrix.npz`` is already the interchange format the pipeline uses
internally: ``feature_utils.select_features`` (Step 14) and
``classifier_utils._load_feature_matrix`` (Step 15 / 14b) both prefer it and
only rebuild from bedGraph when it is absent. Writing a valid ``.npz`` is
therefore sufficient to run Steps 13-15 unchanged on external data.

What this does NOT do
---------------------
This is *framework* validation, not *panel* validation. Public cohorts do not
share the 60-CpG targeted panel, and the OvaPrint deposit in particular ships
anonymised feature identifiers (``FEATURE_1`` ...) with no genomic
coordinates. Consequently:

- the resulting matrix carries ``region_in_88 = False`` for every feature, so
  in-fold selection applies its standard (non-relaxed) thresholds uniformly;
- the depth matrix ``D`` is all zeros, which the >=10x coverage filter in
  ``feature_utils`` detects and skips with a warning — correct behaviour for
  already-processed beta values;
- coordinate-dependent steps (Step 14a blood-baseline filter, Step 18/19
  annotation and cross-reference) cannot be applied and must be skipped.

Datasets
--------
``ovaprint``
    Buckley et al., *Clin Cancer Res* 2023;29(24):5196 (PMC10722131).
    Zenodo 10.5281/zenodo.8395817, CC BY 4.0. 372 samples x 1677 features;
    313 are plasma cfDNA assayed by hybridisation probe capture bisulfite
    sequencing — the same modality and the same benign-versus-malignant
    clinical question as the AIMethFlow cohort.

``gse19711``
    Teschendorff et al., *PLoS ONE* 2009;4(12):e8274 — the UKOPS cohort.
    GEO GSE19711, Illumina 27K, 540 whole-blood samples x 27,578 probes:
    266 ovarian-cancer cases (131 pre-treatment, 135 post-treatment) and 274
    cancer-free postmenopausal controls.

    **This is whole-blood leukocyte gDNA, not plasma cfDNA, and the negative
    class is healthy women rather than benign adnexal masses.** It therefore
    tests the analytical framework on an independent, well-powered ovarian
    case/control cohort; it says nothing about the assay, and nothing about
    the benign-versus-malignant triage question the internal cohort answers.
    Any sentence that blurs those is worse than not running the experiment.

    Its value is that it sits at a third features-per-sample ratio: ~68
    features per sample against the internal cohort's ~160 and OvaPrint's
    ~6.3, which turns the two-point dimensionality claim in the leakage
    demonstration into a three-point one.

Label convention
----------------
Mirrors the AIMethFlow cohort definition: ``0`` benign, ``1`` malignant,
``2`` borderline (retained in the matrix but dropped from binary analyses by
``_load_feature_matrix``). Diagnoses with no counterpart in the AIMethFlow
cohort are excluded by default — for OvaPrint that means the 26 healthy
``Normal`` plasma donors, since every AIMethFlow participant presented with an
adnexal mass. ``--include-normal`` folds them into the benign class as a
sensitivity analysis.

Sub-cohorts
-----------
``--subset`` builds a named row subset instead of the full set. For OvaPrint,
``--subset hgsoc_matched`` reproduces the original publication's own model
definition — 82 HGSOC versus 100 benign masses, with stage IV, borderline,
LGSOC and the non-serous histotypes excluded. The default build keeps every
malignant histotype, which is the harder all-comers question an adnexal-mass
triage clinic actually faces; the two are answers to different questions and
their AUCs must not be compared to each other, only each to its own reference.

For GSE19711, ``--subset all_cases`` adds the 135 post-treatment bloods back
in, reaching the deposit's full n = 540. The default build excludes them: a
blood sample drawn *after* treatment asks whether a treated patient can be
recognised, not whether an untreated one can be detected, and only the second
is the question this pipeline is for. Quote n = 540 for the deposit if useful,
but report the primary result on the 405-sample pre-treatment build.

Outputs
-------
- ``<result-dir>/feature_matrix.npz``            — X, y, sample_ids, feature_names, D, region_in_88
- ``<result-dir>/external_cohort_manifest.csv``  — per-sample id, label, source diagnosis, batch

Example
-------
::

    python src/scripts/prepare_external_matrix.py \\
        --dataset ovaprint \\
        --input-dir data/external/ovaprint \\
        --result-dir result/external/ovaprint

    # the original publication's own comparison, for a like-for-like check
    python src/scripts/prepare_external_matrix.py \\
        --dataset ovaprint --subset hgsoc_matched \\
        --input-dir data/external/ovaprint \\
        --result-dir result/external/ovaprint_hgsoc

    # a second cohort, from a GEO series matrix rather than a CSV pair
    python src/scripts/prepare_external_matrix.py \\
        --dataset gse19711 \\
        --input-dir data/external/GSE19711 \\
        --result-dir result/external/gse19711
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

try:
    from utils.logger import logger
except ImportError:  # standalone use outside the package layout
    import logging
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    logger = logging.getLogger("prepare_external_matrix")


# ─────────────────────────────────────────────────────────────────────────────
# Dataset presets
# ─────────────────────────────────────────────────────────────────────────────

def _gse19711_status(meta):
    """Collapse GSE19711's case/control flag and treatment timing into one label.

    The deposit records these separately: ``sample type`` is free text ending in
    "from Case" / "from Control", and ``pre-treatment sample`` is Yes/No (controls
    are No, since they were never treated). Merging them into a single column lets
    the ordinary ``benign`` / ``malignant`` label lists express the pre-treatment
    default and the ``all_cases`` subset without any extra filter machinery.
    """
    grp = meta["sample type"].str.extract(r"from (Case|Control)\s*$")[0]
    pre = meta["pre-treatment sample"].astype(str).str.strip().eq("Yes")
    return np.where(grp.eq("Control"), "Control",
                    np.where(pre, "Case (pre-treatment)", "Case (post-treatment)"))


DATASETS = {
    "ovaprint": {
        "format":        "csv_pair",
        "beta_file":     "beta_value_matrix.csv.gz",
        "meta_file":     "metadata.csv.gz",
        "id_col":        "ID",
        "label_col":     "Diagnosis",
        "batch_col":     "Source",
        # keep only plasma cfDNA; the 59 tissue samples were used by the
        # original authors for feature discovery, not model evaluation
        "subset_col":    "Sample.type",
        "subset_value":  "PLASMA",
        # beta matrix is features x samples and needs transposing
        "features_are_rows": True,
        "benign":     ["Benign"],
        "malignant":  ["HGSOC", "HGSOC (STAGE IV)", "LGSOC", "Endometrioid",
                       "Clear Cell", "Mucinous", "Epithelial Mixed"],
        "borderline": ["Borderline Serous"],
        # healthy donors — no counterpart in the AIMethFlow cohort
        "optional_benign": ["Normal"],
        # Named sub-cohorts, selectable with --subset. Each one replaces the
        # label lists above; the source files, the plasma subsetting and the
        # matrix orientation are unchanged, so the resulting matrix is a strict
        # row subset of the default build.
        "subsets": {
            "hgsoc_matched": {
                "description":
                    "the original publication's own model set — 82 HGSOC vs "
                    "100 benign masses, stage IV / borderline / LGSOC / "
                    "non-serous histotypes all excluded, exactly as Buckley "
                    "et al. defined it. Use this to compare against their "
                    "reported figures; the default build answers the harder "
                    "all-comers question and is not comparable to them.",
                "benign":          ["Benign"],
                "malignant":       ["HGSOC"],
                "borderline":      [],
                "optional_benign": [],
            },
        },
    },
    "gse19711": {
        "format":        "geo_series_matrix",
        "beta_file":     "GSE19711_series_matrix.txt.gz",
        # Metadata is embedded in the same file's !Sample_* header block, so
        # there is no separate meta_file to read.
        "id_col":        "sample_id",
        "label_col":     "status",
        "batch_col":     "batch",
        "derive_label":  _gse19711_status,
        "features_are_rows": True,
        "benign":     ["Control"],
        # Post-treatment bloods are deliberately absent here and land in the
        # "not in this build's label map" log line; --subset all_cases adds them.
        "malignant":  ["Case (pre-treatment)"],
        "borderline": [],
        "optional_benign": [],
        "subsets": {
            "all_cases": {
                "description":
                    "every case regardless of treatment timing — the deposit's "
                    "full 540 samples, 266 cases vs 274 controls. The 135 "
                    "post-treatment bloods it adds answer a different question "
                    "from detection, so this is a sensitivity arm and the "
                    "405-sample pre-treatment build stays the primary result.",
                "benign":          ["Control"],
                "malignant":       ["Case (pre-treatment)", "Case (post-treatment)"],
                "borderline":      [],
                "optional_benign": [],
            },
            "balanced_batches": {
                "description":
                    "the eight hybridisation batches that contain both classes "
                    "(1-8), dropping batches 9-12. Batches 10 and 11 are 100% "
                    "control (47 samples each) and 9 and 12 are 85-90% control, "
                    "so in the full build a probe that tracks the beadchip "
                    "tracks the label for free. This arm removes that route "
                    "entirely: 144 controls vs 126 pre-treatment cases, every "
                    "batch between 41% and 60% case. If the discrimination "
                    "survives here it is not a batch artefact.",
                "subset_col":   "batch",
                "subset_value": ["1", "2", "3", "4", "5", "6", "7", "8"],
            },
        },
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# Core
# ─────────────────────────────────────────────────────────────────────────────

_TABLE_BEGIN = "!series_matrix_table_begin"


def _read_geo_series_matrix(path):
    """Split a GEO series-matrix file into (beta_df, meta_df).

    The format carries both in one file: a header of ``!Sample_<field>`` lines,
    each holding one value per sample in column order, then the beta table
    between the ``!series_matrix_table_begin`` / ``_end`` markers.

    Sample-level annotation is spread across repeated
    ``!Sample_characteristics_ch1`` lines of the form ``"<field>: <value>"``.
    Each such line is unpacked into its own column when it names a single field
    consistently; lines that mix field names across samples are skipped, since
    there is no unambiguous column to put them in.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"missing input: {path}")

    header, n_skip = {}, 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            n_skip += 1
            if line.startswith(_TABLE_BEGIN):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][len("!Sample_"):]
            header.setdefault(key, []).append([v.strip('"') for v in parts[1:]])
        else:
            raise ValueError(
                f"{path}: no {_TABLE_BEGIN} marker — not a series matrix?")

    if "geo_accession" not in header:
        raise ValueError(f"{path}: no !Sample_geo_accession line")
    meta = {"sample_id": header["geo_accession"][0]}
    for block in header.get("characteristics_ch1", []):
        fields = {v.split(":", 1)[0].strip() for v in block if ":" in v}
        if len(fields) != 1:
            logger.warning(f"  skipping a characteristics line spanning "
                           f"{len(fields)} field names: {sorted(fields)[:4]}")
            continue
        meta[fields.pop()] = [v.split(":", 1)[1].strip() if ":" in v else ""
                              for v in block]
    # Plain single-value fields worth keeping alongside the characteristics.
    for key in ("title", "source_name_ch1"):
        if key in header:
            meta.setdefault(key, header[key][0])
    meta = pd.DataFrame(meta)

    beta = pd.read_csv(path, sep="\t", skiprows=n_skip, index_col=0,
                       low_memory=False)
    # The closing "!series_matrix_table_end" marker parses as a final data row.
    beta = beta[~beta.index.astype(str).str.startswith("!")]
    beta = beta.apply(pd.to_numeric, errors="coerce")
    return beta, meta


def _load_dataset(spec, input_dir):
    """Read the beta matrix and metadata; return (beta_df, meta_df).

    beta_df is returned oriented as samples x features.
    """
    beta_path = os.path.join(input_dir, spec["beta_file"])

    if spec.get("format") == "geo_series_matrix":
        beta, meta = _read_geo_series_matrix(beta_path)
    else:
        meta_path = os.path.join(input_dir, spec["meta_file"])
        for p in (beta_path, meta_path):
            if not os.path.exists(p):
                raise FileNotFoundError(f"missing input: {p}")
        beta = pd.read_csv(beta_path, index_col=0)
        meta = pd.read_csv(meta_path)

    if spec["features_are_rows"]:
        beta = beta.T
    logger.info(f"  beta matrix : {beta.shape[0]} samples x {beta.shape[1]} features")
    logger.info(f"  metadata    : {meta.shape[0]} rows, columns {list(meta.columns)}")
    return beta, meta


def _resolve_subset(spec, subset):
    """Overlay a named sub-cohort's label lists onto *spec*.

    Returns ``(spec, description)``; ``description`` is None when no subset was
    requested, in which case *spec* is returned unchanged.
    """
    if not subset:
        return spec, None
    subsets = spec.get("subsets", {}) or {}
    if subset not in subsets:
        raise KeyError(f"unknown subset {subset!r}; known: {sorted(subsets)}")
    overlay = dict(subsets[subset])
    description = overlay.pop("description", subset)
    return {**spec, **overlay}, description


def _assign_labels(spec, meta, include_normal):
    """Map source diagnoses onto the AIMethFlow 0/1/2 convention.

    Returns a DataFrame with columns [sample_id, y, diagnosis, batch],
    restricted to the samples that carry a usable label.
    """
    benign = list(spec["benign"])
    if include_normal:
        benign += list(spec.get("optional_benign", []))

    mapping = {}
    mapping.update({d: 0 for d in benign})
    mapping.update({d: 1 for d in spec["malignant"]})
    mapping.update({d: 2 for d in spec["borderline"]})

    df = meta.copy()
    if spec.get("subset_col"):
        # A list keeps the several-values case (a set of batches) in the same
        # declarative form as the single-value case (one sample type).
        wanted = spec["subset_value"]
        wanted = list(wanted) if isinstance(wanted, (list, tuple, set)) else [wanted]
        before = len(df)
        df = df[df[spec["subset_col"]].astype(str).isin([str(v) for v in wanted])]
        shown = wanted[0] if len(wanted) == 1 else f"one of {wanted}"
        logger.info(f"  subset {spec['subset_col']}=={shown}: "
                    f"{len(df)}/{before} samples")

    # Some deposits spread the clinical grouping across several raw columns.
    # A dataset may supply a function that folds them into one label column,
    # which the label lists below then address like any other diagnosis field.
    if spec.get("derive_label"):
        df[spec["label_col"]] = spec["derive_label"](df)
        counts = df[spec["label_col"]].value_counts(dropna=False)
        logger.info(f"  derived {spec['label_col']}: "
                    + ", ".join(f"{k} = {v}" for k, v in counts.items()))

    unknown = sorted(set(df[spec["label_col"]]) - set(mapping))
    if unknown:
        n_drop = int(df[spec["label_col"]].isin(unknown).sum())
        logger.info(f"  excluded diagnoses (not in this build's label map): "
                    f"{unknown} — {n_drop} samples")

    df = df[df[spec["label_col"]].isin(mapping)].copy()
    out = pd.DataFrame({
        "sample_id": df[spec["id_col"]].astype(str).values,
        "y":         df[spec["label_col"]].map(mapping).astype(int).values,
        "diagnosis": df[spec["label_col"]].values,
        "batch":     (df[spec["batch_col"]].values
                      if spec.get("batch_col") else ""),
    })
    return out


def build_external_matrix(dataset, input_dir, result_dir, include_normal=False,
                          subset=None):
    """Write ``feature_matrix.npz`` + a manifest for one external cohort."""
    if dataset not in DATASETS:
        raise KeyError(f"unknown dataset {dataset!r}; known: {sorted(DATASETS)}")
    spec, subset_desc = _resolve_subset(DATASETS[dataset], subset)

    logger.info(f"━━ Preparing external cohort: {dataset} ━━")
    if subset_desc:
        logger.info(f"  subset      : {subset} — {subset_desc}")
        if include_normal:
            logger.warning("  --include-normal has no effect on this subset "
                           "(it defines its own benign class)")
    beta, meta = _load_dataset(spec, input_dir)
    labels = _assign_labels(spec, meta, include_normal)

    missing = sorted(set(labels["sample_id"]) - set(beta.index.astype(str)))
    if missing:
        logger.warning(f"  {len(missing)} labelled samples absent from the beta "
                       f"matrix and dropped, e.g. {missing[:5]}")
        labels = labels[~labels["sample_id"].isin(missing)]

    beta.index = beta.index.astype(str)
    X = beta.loc[labels["sample_id"].tolist()].to_numpy(dtype=np.float64)
    y = labels["y"].to_numpy(dtype=np.int64)
    sample_ids    = labels["sample_id"].tolist()
    feature_names = [str(c) for c in beta.columns]

    # D (per-site read depth) is meaningless for pre-processed beta values.
    # All-zero is the documented sentinel: feature_utils detects it and skips
    # the >=10x coverage filter rather than discarding every site.
    D = np.zeros_like(X, dtype=np.float64)
    # No 88-candidate-region prior exists for an external panel, so in-fold
    # selection applies its standard thresholds to every feature.
    region_in_88 = np.zeros(X.shape[1], dtype=bool)

    os.makedirs(result_dir, exist_ok=True)
    fm_path = os.path.join(result_dir, "feature_matrix.npz")
    np.savez(fm_path,
             X=X, y=y,
             sample_ids=np.array(sample_ids, dtype=object),
             feature_names=np.array(feature_names, dtype=object),
             D=D,
             region_in_88=region_in_88)

    man_path = os.path.join(result_dir, "external_cohort_manifest.csv")
    labels.to_csv(man_path, index=False, encoding="utf-8")

    n0, n1, n2 = int((y == 0).sum()), int((y == 1).sum()), int((y == 2).sum())
    logger.info(f"  → {fm_path}")
    logger.info(f"    X = {X.shape[0]} samples x {X.shape[1]} features, "
                f"{np.isnan(X).mean() * 100:.2f}% missing")
    logger.info(f"    y = {n0} benign / {n1} malignant / {n2} borderline "
                f"(borderline excluded from binary analyses)")
    logger.info(f"  → {man_path}")
    return X, y, sample_ids, feature_names


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build feature_matrix.npz from a public methylation beta matrix.")
    ap.add_argument("--dataset", default="ovaprint", choices=sorted(DATASETS),
                    help="dataset preset (default: ovaprint)")
    ap.add_argument("--input-dir", required=True,
                    help="directory holding the downloaded beta matrix + metadata")
    ap.add_argument("--result-dir", required=True,
                    help="output directory; feature_matrix.npz is written here")
    ap.add_argument("--include-normal", action="store_true",
                    help="fold healthy-donor samples into the benign class "
                         "(sensitivity analysis; excluded by default because the "
                         "AIMethFlow cohort contains no healthy controls)")
    ap.add_argument("--subset", default=None,
                    help="named sub-cohort to build instead of the full set. "
                         "Known: " + "; ".join(
                             f"{ds} → {', '.join(sorted(spec['subsets']))}"
                             for ds, spec in sorted(DATASETS.items())
                             if spec.get("subsets")))
    args = ap.parse_args(argv)

    build_external_matrix(args.dataset, args.input_dir, args.result_dir,
                          include_normal=args.include_normal,
                          subset=args.subset)
    return 0


if __name__ == "__main__":
    sys.exit(main())
