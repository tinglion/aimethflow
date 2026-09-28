#!/usr/bin/env python3
"""report.cohort — Step 19b independent external cohort re-run section.

Renders ``external_cohort_summary.csv`` (produced by
``scripts.modules.external_cohort``): the production Step 15 re-run unmodified
on one or more independent cohorts, shown side by side with the internal
cohort, plus the optional de-novo panel comparison.

Distinct from ``report.external`` (Step 19), which cross-references the panel
against public *reference* datasets rather than running an independent cohort.
"""

import os
import logging

from scripts.modules.report.helpers import _read_csv_safe

logger = logging.getLogger(__name__)

_CLF_ORDER = ["LR", "SVM", "RF", "XGB", "LGBM"]

ARM_LABEL = {
    "full_matrix": ("full matrix, in-fold selection",
                    "No panel prior: selection runs over the whole feature pool, confined to "
                    "each training fold. This is the primary external result."),
    "denovo_panel": ("de-novo panel",
                     "Step 14 run on this cohort to derive a panel from scratch, then Step 15 "
                     "restricted to it — an independent replication of the panel-definition "
                     "question raised on the internal cohort."),
    "hgsoc_matched": ("source publication's own definition",
                      "Sensitivity arm: the cohort restricted to the comparison the "
                      "original authors reported — high-grade serous versus benign "
                      "masses, with stage IV, borderline and the non-serous histotypes "
                      "excluded. The primary arm keeps every malignant histotype, which "
                      "is the harder all-comers question; this arm exists so the "
                      "reproduction can be judged against the published figures on the "
                      "definition those figures were computed under."),
    "all_cases": ("post-treatment cases included",
                  "Sensitivity arm: cases drawn after treatment added back to the "
                  "pre-treatment ones, reaching the deposit's full sample count. "
                  "Excluded from the primary arm because a blood sample drawn after "
                  "treatment asks whether a treated patient can be recognised, not "
                  "whether an untreated one can be detected — and only the second is "
                  "the question this pipeline is built for."),
    "balanced_batches": ("batch-balanced subset",
                         "Sensitivity arm: the hybridisation batches that contain both "
                         "classes, dropping those that are wholly or almost wholly control. "
                         "In the full deposit a probe that tracks the beadchip tracks the "
                         "label for free; this arm removes that route entirely. Read the "
                         "lower AUC with care &mdash; the discriminating probes keep their "
                         "univariate effect here, but their &Delta;&beta; falls below the "
                         "<code>sample_delta_beta</code> gate, so the pre-filter relaxes and "
                         "hands the classifiers a larger, weaker feature set. That is a "
                         "threshold-portability effect, not a batch artefact."),
    "with_healthy_donors": ("healthy donors included",
                            "Sensitivity arm: the deposit's healthy donors folded into the "
                            "benign class. Excluded from the primary arm because every "
                            "AIMethFlow participant presented with an adnexal mass — adding "
                            "easy negatives raises AUC for a reason unrelated to the method, "
                            "so this number is not comparable to the internal cohort."),
    "panel_restricted": ("committed panel",
                         "The internal cohort's reported configuration: Step 15 restricted to "
                         "the panel Step 14 committed, which was defined from this cohort's "
                         "own labels."),
}


def _num(val, fmt="{:.3f}"):
    try:
        f = float(val)
    except (TypeError, ValueError):
        return "—"
    if f != f:
        return "—"
    return fmt.format(f)


def _fmt_p(val):
    try:
        f = float(val)
    except (TypeError, ValueError):
        return "—"
    if f != f:
        return "—"
    if f < 1e-3:
        mant, exp = f"{f:.1e}".split("e")
        return f"{mant} &times; 10<sup>{int(exp)}</sup>"
    return f"{f:.3f}"


def _is_sig(val, alpha=0.05):
    """True only when the p-value parses to a real number below alpha."""
    try:
        f = float(val)
    except (TypeError, ValueError):
        return False
    return f == f and f < alpha


def _build_external_cohort_section(result_dir):
    """
    Read external_cohort_summary.csv and build the Step 19b HTML section.
    Returns (html_str, data_dict). data_dict is None when the step has not run.
    """
    records = _read_csv_safe(os.path.join(result_dir, "external_cohort_summary.csv"))
    if not records:
        return (
            '<div class="section"><h2>Step 19b. Independent External Cohort Re-run</h2>'
            '<p style="color:#aaa;font-size:.9em">Step 19b not yet run — prepare a cohort with '
            '<code>src/scripts/prepare_external_matrix.py</code>, list its result directory in '
            '<code>config.toml [external_cohort] dirs</code>, then execute '
            '<code>--step external_cohort</code> to populate this section.</p></div>',
            None,
        )

    # Group by (cohort, arm), preserving first-seen order — the writer emits the
    # internal cohort first, then each external cohort in configured order.
    groups, order = {}, []
    for r in records:
        key = (r.get("cohort", ""), r.get("arm", ""))
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(r)
    for key in groups:
        groups[key].sort(key=lambda r: (_CLF_ORDER.index(r.get("classifier"))
                                        if r.get("classifier") in _CLF_ORDER
                                        else len(_CLF_ORDER)))

    has_discovery = any(arm == "denovo_panel" for _, arm in order)

    # ── AUC table: one row per (cohort, arm), one column per classifier ──────
    # Arms go down rather than across: the count of arms grows with every cohort
    # and every sensitivity analysis, while the classifier lineup is fixed at
    # five, so this is the axis that stays readable as the study grows.
    header_cells = "".join(
        f"<th style='text-align:right;border-left:1px solid #e5e7eb'>{clf}</th>"
        for clf in _CLF_ORDER)

    rows_per_cohort = {}
    for cohort, _arm in order:
        rows_per_cohort[cohort] = rows_per_cohort.get(cohort, 0) + 1

    body, seen_cohorts = "", set()
    for cohort, arm in order:
        label, _ = ARM_LABEL.get(arm, (arm, ""))
        n = groups[(cohort, arm)][0].get("n_samples", "?")

        cohort_cell = ""
        if cohort not in seen_cohorts:
            seen_cohorts.add(cohort)
            cohort_cell = (f"<td rowspan='{rows_per_cohort[cohort]}' "
                           f"style='vertical-align:top;border-top:2px solid #cbd5e1'>"
                           f"<strong>{cohort}</strong></td>")
        top = "border-top:2px solid #cbd5e1" if cohort_cell else ""

        cells = ""
        for clf in _CLF_ORDER:
            row = next((r for r in groups[(cohort, arm)]
                        if r.get("classifier") == clf), None)
            if row is None:
                cells += (f"<td style='text-align:right;"
                          f"border-left:1px solid #e5e7eb;{top}'>—</td>")
                continue
            ci = f"{_num(row.get('auc_ci_low'))}–{_num(row.get('auc_ci_high'))}"
            cells += (
                f"<td style='text-align:right;border-left:1px solid #e5e7eb;{top}'>"
                f"<strong>{_num(row.get('auc'), '{:.4f}')}</strong>"
                f"<br/><span style='color:#666;font-size:.82em'>{ci}</span></td>"
            )
        body += (
            f"<tr>{cohort_cell}"
            f"<td style='{top}'>{label}"
            f"<br/><span style='color:#666;font-size:.82em'><code>{arm}</code></span></td>"
            f"<td style='text-align:right;{top}'>{_num(n, '{:.0f}')}</td>"
            f"{cells}</tr>"
        )

    # ── Optional: the de-novo panel effect, per external cohort ─────────────
    disc_html = ""
    if has_discovery:
        disc_rows = ""
        for cohort, arm in order:
            if arm != "denovo_panel":
                continue
            full = {r.get("classifier"): r for r in groups.get((cohort, "full_matrix"), [])}
            for r in groups[(cohort, arm)]:
                clf = r.get("classifier")
                base = full.get(clf, {})
                try:
                    d = float(r.get("delta_auc_vs_full"))
                    p = float(r.get("delong_p"))
                except (TypeError, ValueError):
                    d, p = float("nan"), float("nan")
                if p == p and p < 0.05:
                    colour, weight = ("#b91c1c" if d > 0 else "#1d4ed8"), "bold"
                else:
                    colour, weight = "#6b7280", "normal"
                delta_txt = f"{d:+.4f}" if d == d else "—"
                disc_rows += (
                    f"<tr><td>{cohort}</td><td><strong>{clf}</strong></td>"
                    f"<td style='text-align:right'>{_num(base.get('auc'), '{:.4f}')}</td>"
                    f"<td style='text-align:right'>{_num(r.get('auc'), '{:.4f}')}</td>"
                    f"<td style='text-align:right;color:{colour};font-weight:{weight}'>"
                    f"{delta_txt}</td>"
                    f"<td style='text-align:right'>{_fmt_p(r.get('delong_p'))}</td></tr>"
                )
        if disc_rows:
            # State the outcome rather than leaving the reader to scan p-values.
            sig = [r for r in records if r.get("arm") == "denovo_panel"
                   and _is_sig(r.get("delong_p"))]
            n_tested = sum(1 for r in records if r.get("arm") == "denovo_panel")
            if sig:
                verdict = (
                    f"<strong>{len(sig)} of {n_tested}</strong> classifiers show a significant "
                    f"difference (p &lt; 0.05): on this cohort too, a panel defined from the whole "
                    f"cohort's labels carries information into the held-out folds.")
            else:
                verdict = (
                    f"<strong>None of the {n_tested}</strong> classifiers show a significant "
                    f"difference (all p &gt; 0.05): on this cohort the de-novo panel neither helps "
                    f"nor inflates, so the panel-definition step is not detectably leaking here. "
                    f"That is consistent with its lower dimensionality burden — the external "
                    f"matrix selects fewer features per fold from a larger sample set, leaving a "
                    f"supervised filter far less room to fit the label vector by chance.")
            disc_html = f"""
<h3 style="margin-bottom:4px">Effect of a cohort-defined panel, tested externally</h3>
<p style="color:#555;font-size:.9em;line-height:1.7">
  Deriving a panel de novo on the external cohort and restricting the classifiers to it
  replays, on independent data, the structure Step 15b probes internally: a supervised
  selection stage sitting upstream of the cross-validation. Both runs score the same samples
  over the same fold partition, so the AUCs are paired and the difference is testable —
  the test decides whether the effect is present here, and a null result is as informative
  as a positive one, since the size of this leak tracks the cohort's features-per-sample
  burden rather than being a fixed property of the protocol.
</p>
<table>
<thead><tr><th>Cohort</th><th>Classifier</th><th style="text-align:right">Full matrix</th>
<th style="text-align:right">De-novo panel</th><th style="text-align:right">&Delta;AUC</th>
<th style="text-align:right">DeLong p</th></tr></thead>
<tbody>{disc_rows}</tbody>
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Red = significantly higher with the panel (p &lt; 0.05); blue = significantly lower;
  grey = not significant.
</p>
<p style="color:#555;font-size:.9em;line-height:1.7;background:#f8fafc;
          border-left:3px solid #94a3b8;padding:8px 12px;margin:8px 0">{verdict}</p>"""

    n_external = len({c for c, a in order if a != "panel_restricted"})

    html = f"""<div class="section">
<h2>Step 19b. Independent External Cohort Re-run</h2>
<p style="color:#555;line-height:1.7">
  This step re-runs the <strong>production Step 15 routine, unmodified</strong>, on an
  independent cohort whose feature matrix was prepared by
  <code>prepare_external_matrix.py</code>. Because <code>feature_matrix.npz</code> is the
  pipeline's interchange format, Steps 1&ndash;12 are bypassed entirely and the external data
  enters at the matrix level &mdash; no aligner, reference genome or panel BED is involved.
</p>
<p style="color:#7c2d12;background:#fff7ed;border-left:3px solid #fb923c;padding:8px 12px;
          line-height:1.7;font-size:.92em;margin:10px 0">
  <strong>This validates the framework, not the panel.</strong> No public cohort shares this
  assay's probe set, so the committed panel cannot be mapped onto external data and nothing
  here speaks to the assay's transferability. What is demonstrated is that the same code, the
  same hyperparameters and the same leakage-free protocol reproduce comparable discrimination
  on a different cohort, a different assay chemistry and a different feature space.
</p>
<p style="color:#555;font-size:.9em;line-height:1.7">
  Exactly three internal-cohort bindings are overridden per external run, all of them cohort
  bindings rather than method choices: <code>[ml] include_ca125_baseline</code> is forced off
  (CA-125 is read from the internal sample Excel and keyed by internal sample IDs, so no
  external sample can resolve to a value &mdash; even for a deposit that ships CA-125 &mdash; and the
  benchmark reduces to the five cfDNA-only classifiers); <code>[ml] exclude_samples</code> is
  cleared (those IDs were resolved against the internal cohort); and no <code>panel_cpg</code>
  file exists externally, so selection runs over the full feature pool inside each training
  fold. Everything else is the unmodified pipeline reading <code>config.toml</code>.
</p>
<p style="color:#666;font-size:.87em;margin:0 0 10px">
  External cohorts: <strong>{n_external}</strong> &middot;
  comparison restricted to the five cfDNA-only classifiers, the models both the internal and
  external cohorts have.
</p>

<h3 style="margin-bottom:4px">Internal vs external, side by side</h3>
<table>
<thead>
<tr><th>Cohort</th><th>Arm</th><th style="text-align:right">n</th>{header_cells}</tr>
</thead>
<tbody>{body}</tbody>
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Arms differ in feature space, so read across a row before reading down a column:
  {" &middot; ".join(f"<code>{a}</code> = {ARM_LABEL.get(a, (a, ''))[0]}"
                     for a in dict.fromkeys(a for _, a in order))}.
  Full data: <code>external_cohort_summary.csv</code>
</p>
{disc_html}
</div>"""

    return html, {"rows": records, "cohorts": [{"cohort": c, "arm": a} for c, a in order]}
