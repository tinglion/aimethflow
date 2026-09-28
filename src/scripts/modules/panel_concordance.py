#!/usr/bin/env python3
"""panel_concordance — region-level Δβ direction concordance for the committed panel.

The one analysis in this study that tests the **panel** rather than the framework.

Every other external cohort enters at the matrix level and validates the pipeline:
OvaPrint ships anonymised ``FEATURE_n`` identifiers with no coordinates, and
GSE19711's 27K array overlaps the panel at a single site. Neither can say anything
about the 60-CpG panel itself. An EPIC array can, because EPIC probes carry genomic
coordinates: the panel's hg38 positions can be lifted to hg19 and matched against
the manifest, so the *same loci* can be measured on independent data.

What this establishes is **direction concordance, not panel transfer**. The panel is
still not scoreable on another platform — probe-level β from a single EPIC probe is
not the region-level methylation the targeted bisulfite assay reports, and the
matching window used here is deliberately region-scale rather than base-exact. What
can be asked is narrower and still worth asking: at the loci the panel commits to,
does methylation move in the *same direction* in independent tumour tissue as it does
in plasma?

Design notes:

* **Two analytes, one direction.** Plasma Δβ is malignant minus benign among women
  presenting with an adnexal mass. Tissue Δβ is HGSOC tumour minus normal fallopian
  tube. These are different comparisons on different material; only the sign is
  comparable, which is exactly why the output is a sign table and not a correlation
  of effect sizes.
* **Precedent.** OvaPrint's own paper (Buckley et al., *Clin Cancer Res*
  2023;29(24):5196) used this same deposit in this same way as its external
  consistency check.
* **Pseudo-replication is the main statistical hazard here, and it is handled
  explicitly.** The panel is a *targeted* design: its 65 sites sit in a couple of
  dozen small genomic clusters, and once matched at region scale several panel sites
  routinely resolve to the same EPIC probe. Counting agreements per *site* would
  therefore count the same measurement several times and report a p-value far smaller
  than the data support. Concordance is counted at three nested levels — site, distinct
  probe, and genomic region — and the **region level is the headline**, because
  neighbouring CpGs inside one differentially methylated region are not independent
  observations of anything.
* **The core package is read-only.** Coordinate machinery — chain parsing, liftover,
  manifest parsing, series-matrix streaming — already exists in
  ``aimethflow_core.external_validate`` and is imported here rather than
  reimplemented or modified.

Outputs (written to ``--out``):

``panel_concordance_sites.csv``
    One row per panel site: liftover result, matched EPIC probe and its distance,
    plasma Δβ and p, tissue Δβ and p, and whether the signs agree.
``panel_concordance_summary.csv``
    One row per matching window × counting level (site / probe / region): how many
    units match, how many are testable, the concordant fraction and its binomial
    p-value against chance.
"""

import os
import logging

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, binomtest

# The core package ships the coordinate machinery this analysis needs. It is
# installed as a wheel and is not modified by this repository; these helpers are
# private by naming convention but are the established reuse path here.
from aimethflow_core.external_validate import (
    _parse_chain,
    _liftover,
    _build_probe_map,
    _parse_series_matrix,
)

logger = logging.getLogger(__name__)

# Matching windows (bp) reported as a ladder rather than a single choice. ±2 bp is
# the base-exact criterion Step 19 already uses against the 450K array, and it is
# carried here so the EPIC result can be read against that published 10.8% figure.
# The wider windows are region-scale: a panel site and an EPIC probe a few hundred
# bases apart sit inside the same differentially methylated region, which is the
# resolution at which the claim is being made.
MATCH_WINDOWS = [2, 50, 100, 250, 500]

# The window the headline result is quoted at.
PRIMARY_WINDOW = 500

# Panel sites closer together than this are treated as one region — one independent
# observation — when concordance is counted at region level. 1 kb is the scale of the
# amplicons the targeted assay reads, so CpGs within it are measuring the same
# methylation event rather than independent ones.
REGION_GAP_BP = 1000


def _assign_regions(features, gap=REGION_GAP_BP):
    """Group panel sites into genomic clusters. Returns {feature: region_id}."""
    parsed = []
    for feat in features:
        chrom, pos = feat.split(":", 1)
        parsed.append((chrom, int(pos), feat))
    parsed.sort()

    regions, rid, last = {}, -1, None
    for chrom, pos, feat in parsed:
        if last is None or chrom != last[0] or pos - last[1] > gap:
            rid += 1
        regions[feat] = rid
        last = (chrom, pos)
    return regions


def _binom_greater(n_conc, n_test, p_null=0.5):
    """One-sided binomial p for concordance exceeding ``p_null``. NaN when n_test == 0."""
    if not n_test or not (0.0 < p_null < 1.0):
        return float("nan")
    return float(binomtest(n_conc, n_test, p_null, alternative="greater").pvalue)


# Random EPIC probes sampled to establish what "agreement by chance" actually means
# in this tumour comparison. See _expected_concordance below.
N_BACKGROUND_PROBES = 5000
BACKGROUND_SEED = 0


def _expected_concordance(plasma_signs, hyper_rate):
    """
    The concordance a *directionally biased* panel would reach against a
    *directionally biased* tumour by chance alone.

    Testing agreement against 0.5 assumes both sides are coin flips. Neither is.
    This panel is overwhelmingly hypermethylation-driven, and HGSOC tissue is itself
    globally hypermethylated at the kind of CpG-dense loci a methylation panel
    selects — so a panel that had never seen this disease would still agree with it
    most of the time. ``hyper_rate`` is the observed fraction of random EPIC probes
    hypermethylated in this comparison, and the expected concordance is the chance
    that a randomly directed tumour effect happens to match each panel unit's sign:

        E = f_hyper * hyper_rate + f_hypo * (1 - hyper_rate)

    Testing against E rather than 0.5 is the difference between showing the panel's
    loci behave like ovarian cancer loci and merely showing that hypermethylation is
    common in tumours.
    """
    signs = [s for s in plasma_signs if s != 0]
    if not signs or hyper_rate != hyper_rate:
        return float("nan")
    f_hyper = sum(1 for s in signs if s > 0) / len(signs)
    return f_hyper * hyper_rate + (1.0 - f_hyper) * (1.0 - hyper_rate)

_EPIC_MANIFEST_CANDIDATES = [
    "GPL21145_platform.txt",
    "GPL21145_platform.txt.gz",
    "EPIC.hg19.manifest.tsv.gz",
]


def _find_epic_manifest(data_dir):
    """Locate the EPIC probe coordinate file. Returns a path, or "" if absent."""
    for fname in _EPIC_MANIFEST_CANDIDATES:
        path = os.path.join(data_dir, fname)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            logger.info(f"  EPIC manifest: {path}")
            return path
    return ""


def _parse_series_labels(matrix_gz, label_field="tissue"):
    """
    Read a GEO series matrix header and return {gsm_id: characteristic_value}.

    Only the ``!Sample_geo_accession`` line and the ``!Sample_characteristics_ch1``
    line carrying ``label_field`` are needed, so the (very large) table body is
    never touched — the loop stops at ``!series_matrix_table_begin``.
    """
    import gzip

    opener = gzip.open if matrix_gz.endswith(".gz") else open
    gsms, labels = None, None

    with opener(matrix_gz, "rt", encoding="latin-1", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key, values = parts[0], parts[1:]
            if key == "!Sample_geo_accession":
                gsms = values
            elif key == "!Sample_characteristics_ch1" and values:
                # A deposit emits one such line per characteristic; keep the one
                # whose values carry the field we were asked for.
                if values[0].lower().startswith(f"{label_field}:"):
                    labels = [v.split(":", 1)[1].strip() if ":" in v else ""
                              for v in values]

    if gsms is None or labels is None:
        raise ValueError(f"{os.path.basename(matrix_gz)}: could not read "
                         f"sample accessions and '{label_field}' characteristics")
    if len(gsms) != len(labels):
        raise ValueError(f"{os.path.basename(matrix_gz)}: {len(gsms)} accessions "
                         f"but {len(labels)} '{label_field}' values")
    return dict(zip(gsms, labels))


def _match_probes_windowed(panel_features, chain_blocks, probe_map, windows):
    """
    Liftover each hg38 panel site to hg19 and find the nearest EPIC probe within
    each window in ``windows``.

    A thin generalisation of the core's ``_match_probes``, which resolves a single
    fixed window. Returns ``{feature: {"hg19": (chrom, pos) or None,
    "matches": {window: (cg_id, distance) or (None, None)}}}``.
    """
    from collections import defaultdict
    import bisect

    probes_by_chrom = defaultdict(list)
    for (chrom, pos), cg_id in probe_map.items():
        probes_by_chrom[chrom].append((pos, cg_id))
    for chrom in probes_by_chrom:
        probes_by_chrom[chrom].sort()

    max_window = max(windows)
    out = {}

    for feat in panel_features:
        chrom38, pos_str = feat.split(":", 1)
        pos38 = int(pos_str)
        rec = {"hg19": None, "matches": {w: (None, None) for w in windows}}

        hg19 = _liftover(chrom38, pos38, chain_blocks)
        if hg19 is None:
            out[feat] = rec
            continue
        rec["hg19"] = hg19
        chrom19, pos19 = hg19

        probe_list = probes_by_chrom.get(chrom19, [])
        if not probe_list:
            out[feat] = rec
            continue

        # Scan once at the widest window, then assign the nearest hit to every
        # window it also satisfies — the nearest probe within ±w is by definition
        # the nearest probe within ±max_window whenever the distance allows.
        positions = [p[0] for p in probe_list]
        idx = bisect.bisect_left(positions, pos19 - max_window)
        best_cg, best_dist = None, max_window + 1
        while idx < len(probe_list) and probe_list[idx][0] <= pos19 + max_window:
            dist = abs(probe_list[idx][0] - pos19)
            if dist < best_dist:
                best_dist, best_cg = dist, probe_list[idx][1]
            idx += 1

        if best_cg is not None:
            for w in windows:
                if best_dist <= w:
                    rec["matches"][w] = (best_cg, best_dist)

        out[feat] = rec

    return out


def run_panel_concordance(panel_csv, wilcoxon_csv, series_matrix, chain_gz,
                          epic_manifest, out_dir,
                          case_label="High grade serous ovarian cancer",
                          control_label="Normal Fallopian tube control",
                          windows=None, n_background=N_BACKGROUND_PROBES):
    """
    Region-level Δβ direction concordance between the committed panel's plasma
    effects and an independent EPIC tissue cohort.

    Returns (sites_df, summary_df).
    """
    windows = list(windows or MATCH_WINDOWS)
    os.makedirs(out_dir, exist_ok=True)

    # ── 1. The panel, and its plasma-side direction ──────────────────────────
    panel = pd.read_csv(panel_csv)
    features = panel["feature"].astype(str).tolist()
    regions = _assign_regions(features)
    logger.info(f"Panel: {len(features)} committed CpG sites (hg38) in "
                f"{len(set(regions.values()))} genomic clusters "
                f"(>{REGION_GAP_BP} bp apart)")

    wilcox = pd.read_csv(wilcoxon_csv).set_index("feature")
    # sample_delta in wilcoxon_stats.csv is unsigned; the direction that matters
    # here is recovered from the two group means.
    plasma = {}
    for feat in features:
        if feat not in wilcox.index:
            continue
        row = wilcox.loc[feat]
        plasma[feat] = {
            "plasma_benign_mean": float(row["benign_mean"]),
            "plasma_malignant_mean": float(row["malignant_mean"]),
            "plasma_delta": float(row["malignant_mean"]) - float(row["benign_mean"]),
            "plasma_p": float(row["wilcox_pval"]),
            "in_88_regions": bool(row["in_88_regions"]),
        }
    logger.info(f"  Plasma direction resolved for {len(plasma)}/{len(features)} sites")

    # ── 2. Liftover hg38 → hg19, then match to EPIC probes ──────────────────
    logger.info("Parsing liftOver chain (hg38 → hg19) ...")
    chain_blocks = _parse_chain(chain_gz)
    if not chain_blocks:
        raise RuntimeError(f"Could not parse chain file: {chain_gz}")

    logger.info("Parsing EPIC manifest ...")
    probe_map = _build_probe_map(epic_manifest)
    if not probe_map:
        raise RuntimeError(f"Could not parse EPIC manifest: {epic_manifest}")
    logger.info(f"  EPIC probes with hg19 coordinates: {len(probe_map):,}")

    matched = _match_probes_windowed(features, chain_blocks, probe_map, windows)
    n_lifted = sum(1 for r in matched.values() if r["hg19"])
    logger.info(f"  Liftover succeeded for {n_lifted}/{len(features)} sites")
    for w in windows:
        n = sum(1 for r in matched.values() if r["matches"][w][0])
        logger.info(f"  Matched within ±{w:>4} bp: {n:>3}/{len(features)} "
                    f"({100.0 * n / len(features):.1f}%)")

    # ── 3. Tissue β for every probe any window reached ───────────────────────
    probes_needed = {cg for r in matched.values()
                     for cg, _ in r["matches"].values() if cg}
    if not probes_needed:
        raise RuntimeError("No panel site matched any EPIC probe at any window — "
                           "nothing to compare.")

    # Random probes for the empirical null, read in the same pass over the file —
    # the series matrix is ~800 MB, so it is scanned once, not twice.
    rng = np.random.default_rng(BACKGROUND_SEED)
    all_probes = sorted(set(probe_map.values()) - probes_needed)
    n_bg = min(n_background, len(all_probes))
    background = set(rng.choice(np.array(all_probes, dtype=object),
                                size=n_bg, replace=False).tolist()) if n_bg else set()
    logger.info(f"  Background: {len(background)} random EPIC probes "
                f"for the empirical null")

    logger.info(f"Reading {os.path.basename(series_matrix)} "
                f"for {len(probes_needed)} panel + {len(background)} background "
                f"probes ...")
    labels = _parse_series_labels(series_matrix)
    betas = _parse_series_matrix(series_matrix, probes_needed | background)
    if betas.empty:
        raise RuntimeError("Series matrix returned no rows for the matched probes.")

    case_gsms = [g for g in betas.columns if labels.get(g) == case_label]
    ctrl_gsms = [g for g in betas.columns if labels.get(g) == control_label]
    logger.info(f"  Tissue cohort: {len(case_gsms)} '{case_label}' vs "
                f"{len(ctrl_gsms)} '{control_label}'")
    if not case_gsms or not ctrl_gsms:
        raise RuntimeError("Case or control group is empty — check the labels.")

    tissue = {}
    for cg in betas.index:
        a = pd.to_numeric(betas.loc[cg, case_gsms], errors="coerce").dropna()
        b = pd.to_numeric(betas.loc[cg, ctrl_gsms], errors="coerce").dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        try:
            p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
        except ValueError:
            p = float("nan")
        tissue[cg] = {
            "tissue_case_mean": float(a.mean()),
            "tissue_control_mean": float(b.mean()),
            "tissue_delta": float(a.mean() - b.mean()),
            "tissue_p": p,
            "tissue_n_case": len(a),
            "tissue_n_control": len(b),
        }
    n_panel_probes = sum(1 for cg in tissue if cg in probes_needed)
    logger.info(f"  Tissue Δβ computed for {len(tissue)} probes "
                f"({n_panel_probes} panel, {len(tissue) - n_panel_probes} background)")

    bg_deltas = [tissue[cg]["tissue_delta"] for cg in tissue if cg in background]
    bg_deltas = [d for d in bg_deltas if d == d and np.sign(d) != 0]
    hyper_rate = (sum(1 for d in bg_deltas if d > 0) / len(bg_deltas)
                  if bg_deltas else float("nan"))
    logger.info(f"  Genome-wide null: {len(bg_deltas)} background probes, "
                f"{100.0 * hyper_rate:.1f}% hypermethylated in "
                f"{case_label.lower()}")

    # ── 4. Per-site table ────────────────────────────────────────────────────
    rows = []
    for feat in features:
        rec = matched[feat]
        hg19 = rec["hg19"]
        pl = plasma.get(feat, {})
        cg_primary, dist_primary = rec["matches"].get(PRIMARY_WINDOW, (None, None))
        ti = tissue.get(cg_primary, {}) if cg_primary else {}

        pd_, td_ = pl.get("plasma_delta"), ti.get("tissue_delta")
        concordant = ""
        if pd_ is not None and td_ is not None and pd_ == pd_ and td_ == td_:
            concordant = bool(np.sign(pd_) == np.sign(td_) and np.sign(pd_) != 0)

        row = {
            "feature": feat,
            "region": regions[feat],
            "in_88_regions": pl.get("in_88_regions", ""),
            "hg19": f"{hg19[0]}:{hg19[1]}" if hg19 else "",
            "epic_probe": cg_primary or "",
            "match_distance_bp": dist_primary if dist_primary is not None else "",
            "plasma_benign_mean": pl.get("plasma_benign_mean", ""),
            "plasma_malignant_mean": pl.get("plasma_malignant_mean", ""),
            "plasma_delta": pd_ if pd_ is not None else "",
            "plasma_p": pl.get("plasma_p", ""),
            "tissue_control_mean": ti.get("tissue_control_mean", ""),
            "tissue_case_mean": ti.get("tissue_case_mean", ""),
            "tissue_delta": td_ if td_ is not None else "",
            "tissue_p": ti.get("tissue_p", ""),
            "concordant": concordant,
        }
        for w in windows:
            row[f"probe_w{w}"] = rec["matches"][w][0] or ""
        rows.append(row)

    sites = pd.DataFrame(rows)
    sites_path = os.path.join(out_dir, "panel_concordance_sites.csv")
    sites.to_csv(sites_path, index=False)
    logger.info(f"  Wrote {sites_path}")

    # ── 5. Window ladder, counted at three nested levels ─────────────────────
    # site   — every matched panel site, pseudo-replicated and reported only for
    #          completeness;
    # probe  — distinct EPIC probes, removing sites that resolve to the same
    #          measurement;
    # region — distinct genomic clusters, the only level at which the observations
    #          are plausibly independent, and the one the headline is quoted from.
    summary_rows = []
    for w in windows:
        n_match = sum(1 for f in features if matched[f]["matches"][w][0])

        # Collect testable (feature → probe) pairs once, then aggregate three ways.
        pairs = []
        for feat in features:
            cg, _ = matched[feat]["matches"][w]
            if not cg:
                continue
            pl, ti = plasma.get(feat), tissue.get(cg)
            if not pl or not ti:
                continue
            a, b = pl["plasma_delta"], ti["tissue_delta"]
            if a != a or b != b or np.sign(a) == 0 or np.sign(b) == 0:
                continue
            pairs.append((feat, cg, np.sign(a), np.sign(b)))

        # Each level stores (concordant, plasma_sign) per unit — the sign is needed
        # to work out what concordance that unit would reach by chance.
        levels = {}
        levels["site"] = [(sa == sb, sa) for _f, _c, sa, sb in pairs]

        # One vote per distinct EPIC probe. Where several panel sites share a probe
        # the tissue sign is identical by construction, so the plasma sign is
        # averaged over those sites to give the probe a single direction.
        by_probe = {}
        for feat, cg, sa, sb in pairs:
            by_probe.setdefault(cg, {"plasma": [], "tissue": sb})
            by_probe[cg]["plasma"].append(plasma[feat]["plasma_delta"])
        levels["probe"] = []
        for cg, rec in by_probe.items():
            sa = np.sign(float(np.mean(rec["plasma"])))
            if sa == 0:
                continue
            levels["probe"].append((sa == rec["tissue"], sa))

        # One vote per genomic region, averaging both sides within the region.
        by_region = {}
        for feat, cg, _sa, _sb in pairs:
            r = regions[feat]
            by_region.setdefault(r, {"plasma": [], "tissue": {}})
            by_region[r]["plasma"].append(plasma[feat]["plasma_delta"])
            by_region[r]["tissue"][cg] = tissue[cg]["tissue_delta"]
        levels["region"] = []
        for r, rec in by_region.items():
            sa = np.sign(float(np.mean(rec["plasma"])))
            sb = np.sign(float(np.mean(list(rec["tissue"].values()))))
            if sa == 0 or sb == 0:
                continue
            levels["region"].append((sa == sb, sa))

        for level, flags in levels.items():
            n_test = len(flags)
            n_conc = sum(1 for c, _s in flags if c)
            expected = _expected_concordance([s for _c, s in flags], hyper_rate)
            summary_rows.append({
                "window_bp": w,
                "level": level,
                "n_panel_sites": len(features),
                "n_matched_sites": n_match,
                "pct_matched_sites": round(100.0 * n_match / len(features), 1),
                "n_testable": n_test,
                "n_concordant": n_conc,
                "pct_concordant": (round(100.0 * n_conc / n_test, 1)
                                   if n_test else float("nan")),
                "pct_expected_null": (round(100.0 * expected, 1)
                                      if expected == expected else float("nan")),
                "binom_p_vs_chance": _binom_greater(n_conc, n_test, 0.5),
                "binom_p_vs_null": _binom_greater(n_conc, n_test, expected),
            })

    summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(out_dir, "panel_concordance_summary.csv")
    summary.to_csv(summary_path, index=False)
    logger.info(f"  Wrote {summary_path}")

    logger.info("")
    logger.info("Direction concordance (plasma malignant-vs-benign "
                "vs tissue HGSOC-vs-normal FT)")
    logger.info("  levels: site = pseudo-replicated | probe = distinct EPIC probes | "
                "region = independent clusters (headline)")
    for r in summary_rows:
        pct = ("    —" if r["pct_concordant"] != r["pct_concordant"]
               else f"{r['pct_concordant']:5.1f}%")
        exp = ("  —" if r["pct_expected_null"] != r["pct_expected_null"]
               else f"{r['pct_expected_null']:4.1f}%")
        logger.info(f"  ±{r['window_bp']:>4} bp  {r['level']:<6}: "
                    f"testable {r['n_testable']:>3}, "
                    f"concordant {r['n_concordant']:>3} ({pct}) "
                    f"vs null {exp}, p = {r['binom_p_vs_null']:.3g}")
    logger.info("")
    logger.info("This is direction concordance, not panel transfer: the panel "
                "cannot be scored on another platform.")

    return sites, summary
