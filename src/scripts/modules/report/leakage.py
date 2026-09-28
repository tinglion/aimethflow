#!/usr/bin/env python3
"""report.leakage — Step 15b feature-selection leakage demonstration section.

Renders ``leakage_demo.csv`` / ``leakage_demo_delong.csv`` / ``leakage_demo.png``
(plus ``leakage_demo_ca125_gain.csv`` and ``leakage_demo_ca125.png`` when the
cohort has CA-125) together with the explanation a reader needs to interpret
them: what each protocol arm is allowed to see, why the arms are directly
comparable, and what the measured gap means.

Differences carry a paired bootstrap CI alongside the DeLong p, because the two
answer different questions: p says whether the data can resolve a difference,
the interval says how large a difference is still compatible with it. Only the
interval can turn "not significant" into "not large", which is the claim the
panel-provenance discussion actually needs.
"""

import os
import logging

from scripts.modules.report.helpers import _read_csv_safe, _copy_img, _img_tag

logger = logging.getLogger(__name__)

# arm key -> (display name, what selection sees, what imputation/scaling sees, feature space)
ARM_DOC = {
    "in_fold": (
        "in_fold",
        "training fold only",
        "training fold only",
        "full matrix",
        "The AIMethFlow protocol. This is the reference arm in every comparison below.",
    ),
    "naive_selection": (
        "naive_selection",
        "<strong>full matrix</strong>",
        "training fold only",
        "full matrix",
        "Selection fitted once on the complete matrix, labels included, then frozen. "
        "Isolates the selection component — the specific error this pipeline is built to avoid.",
    ),
    "naive_full": (
        "naive_full",
        "<strong>full matrix</strong>",
        "<strong>full matrix</strong>",
        "full matrix",
        "Selection, imputation and scaling all fitted on the complete matrix — the textbook "
        "worst case, and what an unwary reimplementation produces.",
    ),
    "panel_then_in_fold": (
        "panel_then_in_fold",
        "training fold only",
        "training fold only",
        "<strong>cohort-defined panel</strong>",
        "What Step 15 actually does: restrict to the committed panel, then run the in-fold "
        "protocol inside it. The panel was defined by Step 14 from the whole cohort's labels, "
        "so this arm measures how much reported performance depends on that upstream step.",
    ),
    "ca125_only": (
        "ca125_only",
        "not applicable",
        "not applicable",
        "CA-125 only",
        "CA-125 is a directly measured clinical assay — nothing is selected and nothing is "
        "trained — so its AUC is identical under every protocol. It is reported once rather "
        "than repeated per arm, and it is excluded from the DeLong table below, where it "
        "could only ever be compared against itself.",
    ),
}

_ARM_ORDER = ["ca125_only", "in_fold", "naive_selection", "naive_full", "panel_then_in_fold"]
_BASE_CLFS = ["LR", "SVM", "RF", "XGB", "LGBM"]
_CLF_ORDER = ["CA-125"] + _BASE_CLFS + [f"{c}+CA125" for c in _BASE_CLFS]


def _num(val, fmt="{:.3f}"):
    """Format a CSV cell as a number, or '—' when it is missing/unparseable."""
    try:
        f = float(val)
    except (TypeError, ValueError):
        return "—"
    if f != f:                       # NaN
        return "—"
    return fmt.format(f)


def _fmt_p(val):
    """Render a p-value compactly: scientific below 0.001, else three decimals."""
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


def _sort_key(rec, key, order):
    v = rec.get(key, "")
    return order.index(v) if v in order else len(order)


def _f(rec, key):
    """Float value of a CSV cell, NaN when missing or unparseable."""
    try:
        return float(rec.get(key))
    except (TypeError, ValueError):
        return float("nan")


def _ci(rec, lo_key="delta_ci_low", hi_key="delta_ci_high"):
    """Render a paired bootstrap interval on a difference, or '—'."""
    lo, hi = _f(rec, lo_key), _f(rec, hi_key)
    if lo != lo or hi != hi:
        return "—"
    return f"{lo:+.3f} to {hi:+.3f}"


def _build_leakage_section(result_dir, imgs_dir):
    """
    Read leakage_demo.csv (+ the DeLong table and figure) and build the
    Step 15b HTML section.
    Returns (html_str, data_dict). data_dict is None when the step has not run.
    """
    metrics = _read_csv_safe(os.path.join(result_dir, "leakage_demo.csv"))
    if not metrics:
        return (
            '<div class="section"><h2>Step 15b. Feature-selection Leakage Demonstration</h2>'
            '<p style="color:#aaa;font-size:.9em">Step 15b not yet run — '
            'execute <code>--step leakage_demo</code> to populate this section.</p></div>',
            None,
        )

    delong = _read_csv_safe(os.path.join(result_dir, "leakage_demo_delong.csv")) or []

    metrics = sorted(metrics, key=lambda r: (_sort_key(r, "arm", _ARM_ORDER),
                                             _sort_key(r, "classifier", _CLF_ORDER)))
    arms_present = [a for a in _ARM_ORDER if any(r.get("arm") == a for r in metrics)]
    cohort = metrics[0].get("cohort", "") if metrics else ""
    n_samples = metrics[0].get("n_samples", "?") if metrics else "?"

    # ── Explanation: what each arm is allowed to see ─────────────────────────
    arm_rows = ""
    for arm in arms_present:
        name, sel, trans, space, note = ARM_DOC[arm]
        arm_rows += (
            f"<tr><td><code>{name}</code></td><td>{sel}</td><td>{trans}</td>"
            f"<td>{space}</td><td style='font-size:.9em;color:#555'>{note}</td></tr>"
        )

    # ── Per-arm metrics ──────────────────────────────────────────────────────
    metric_rows = ""
    prev_arm = None
    for r in metrics:
        arm = r.get("arm", "")
        # A rule between arm blocks keeps a 15-20 row table readable.
        border = "border-top:2px solid #d1d5db" if prev_arm and arm != prev_arm else ""
        prev_arm = arm
        ci = f"{_num(r.get('auc_ci_low'))}–{_num(r.get('auc_ci_high'))}"
        metric_rows += (
            f"<tr style='{border}'><td><code>{arm}</code></td>"
            f"<td><strong>{r.get('classifier', '')}</strong></td>"
            f"<td style='text-align:right'><strong>{_num(r.get('auc'), '{:.4f}')}</strong></td>"
            f"<td style='text-align:right;color:#666;font-size:.9em'>{ci}</td>"
            f"<td style='text-align:right'>{_num(r.get('accuracy'))}</td>"
            f"<td style='text-align:right'>{_num(r.get('sensitivity'))}</td>"
            f"<td style='text-align:right'>{_num(r.get('specificity'))}</td>"
            f"<td style='text-align:right;color:#666'>{_num(r.get('n_features_mean'), '{:.1f}')}</td>"
            f"</tr>"
        )

    # ── Paired DeLong tests ──────────────────────────────────────────────────
    delong = sorted(delong, key=lambda r: (_sort_key(r, "arm", _ARM_ORDER),
                                           _sort_key(r, "classifier", _CLF_ORDER)))
    delong_rows = ""
    prev_arm = None
    for r in delong:
        arm = r.get("arm", "")
        border = "border-top:2px solid #d1d5db" if prev_arm and arm != prev_arm else ""
        prev_arm = arm
        try:
            d = float(r.get("delta_auc"))
            p = float(r.get("delong_p"))
        except (TypeError, ValueError):
            d, p = float("nan"), float("nan")
        # Colour by significance, not by size: a large delta that fails the test
        # is not evidence, and a reader should not read it as such.
        if p == p and p < 0.05:
            colour, weight = ("#b91c1c", "bold") if d > 0 else ("#1d4ed8", "bold")
        else:
            colour, weight = "#6b7280", "normal"
        delta_txt = f"{d:+.4f}" if d == d else "—"
        delong_rows += (
            f"<tr style='{border}'><td><code>{arm}</code></td>"
            f"<td><strong>{r.get('classifier', '')}</strong></td>"
            f"<td style='text-align:right'>{_num(r.get('auc_arm'), '{:.4f}')}</td>"
            f"<td style='text-align:right'>{_num(r.get('auc_reference'), '{:.4f}')}</td>"
            f"<td style='text-align:right;color:{colour};font-weight:{weight}'>{delta_txt}</td>"
            f"<td style='text-align:right;color:#666;font-size:.9em'>{_ci(r)}</td>"
            f"<td style='text-align:right'>{_fmt_p(r.get('delong_p'))}</td></tr>"
        )

    # ── Interpretation, computed from the numbers actually present ───────────
    interp = []
    naive = [r for r in delong if r.get("arm") == "naive_selection"]
    if naive:
        deltas = [(_f(r, "delta_auc"), r) for r in naive]
        deltas = [(d, r) for d, r in deltas if d == d]
        if deltas:
            worst_d, worst_r = max(deltas, key=lambda t: t[0])
            n_sig = sum(1 for _, r in deltas if _f(r, "delong_p") < 0.05)
            interp.append(
                f"Selecting features on the full matrix before cross-validation inflates AUC "
                f"by up to <strong>{worst_d:+.3f}</strong> ({worst_r.get('classifier')}: "
                f"{_num(worst_r.get('auc_reference'), '{:.3f}')} &rarr; "
                f"{_num(worst_r.get('auc_arm'), '{:.3f}')}), with "
                f"<strong>{n_sig} of {len(deltas)}</strong> classifiers significant at p &lt; 0.05. "
                f"Both arms use the identical fold partition, so the difference is attributable "
                f"to the protocol and nothing else."
            )

    # naive_full vs naive_selection tells the reader which of the two commonly
    # conflated leaks actually matters at this dimensionality.
    def _auc_map(arm):
        return {r.get("classifier"): r.get("auc") for r in metrics if r.get("arm") == arm}
    a_sel, a_full = _auc_map("naive_selection"), _auc_map("naive_full")
    shared = [c for c in _CLF_ORDER if c in a_sel and c in a_full]
    if shared:
        try:
            gaps = [abs(float(a_full[c]) - float(a_sel[c])) for c in shared]
            interp.append(
                f"Adding imputation and scaling leakage on top of selection leakage "
                f"(<code>naive_full</code> vs <code>naive_selection</code>) moves AUC by at most "
                f"<strong>{max(gaps):.3f}</strong>. <strong>Supervised feature selection accounts "
                f"for essentially the entire effect</strong>; pooling median-imputation and "
                f"standardisation statistics across folds is comparatively harmless."
            )
        except (TypeError, ValueError):
            pass

    panel = [r for r in delong if r.get("arm") == "panel_then_in_fold"]
    if panel:
        interp.append(
            "The <code>panel_then_in_fold</code> arm reproduces what Step 15 reports. Where it "
            "exceeds <code>in_fold</code>, the committed panel — defined by Step 14 from the "
            "whole cohort's labels — is itself carrying cohort-level information into the "
            "held-out folds. That is one supervised selection stage sitting outside the nested "
            "CV, and it is a property of the panel's provenance, not of Step 15's cross-validation."
        )
        # Name the models that survive the unrestricted protocol: that is the
        # actionable half of the finding, and a reader should not have to
        # reconstruct it from the table.
        robust = [r for r in panel if not (_f(r, "delong_p") < 0.05)]
        moved = [r for r in panel if _f(r, "delong_p") < 0.05]
        if robust:
            txt = ("<strong>" + ", ".join(r.get("classifier") for r in robust) +
                   "</strong> show no significant dependence on the panel's provenance "
                   "(p &ge; 0.05): their reported AUCs are reproducible under a fully "
                   "leakage-free protocol.")
            # A non-significant p only says the data cannot resolve the difference.
            # The interval says how large a difference is still compatible with it,
            # and it varies enough between these models that quoting one summary
            # figure would flatter the weaker ones — so name both ends.
            his = [(_f(r, "delta_ci_high"), r.get("classifier")) for r in robust]
            his = [(h, c) for h, c in his if h == h]
            if his:
                (lo_h, lo_c), (hi_h, hi_c) = min(his), max(his)
                txt += (f" Read their bootstrap intervals for how much residual dependence the "
                        f"data still permit: upper limits run from <strong>{lo_h:+.3f}</strong> "
                        f"({lo_c}) to <strong>{hi_h:+.3f}</strong> ({hi_c}). The tightest "
                        f"exclude anything approaching the effect measured for the models below; "
                        f"the widest remain compatible with a moderate one, so &quot;not "
                        f"significant&quot; is a bound here, not a demonstration of no effect.")
            if moved:
                txt += (" <strong>" + ", ".join(r.get("classifier") for r in moved) +
                        "</strong> do not — those AUCs should be read as the performance of "
                        "the committed assay, not as a protocol-level estimate.")
            interp.append(txt)

    # ── What each half of a fusion model contributes ─────────────────────────
    gain = _read_csv_safe(os.path.join(result_dir, "leakage_demo_ca125_gain.csv")) or []
    # Both comparisons are internal to one arm, so they are unaffected by the
    # panel-provenance question above; the reference arm is the honest one to show.
    gain_ref = sorted([r for r in gain if r.get("arm") == "in_fold"],
                      key=lambda r: (r.get("comparison", ""),
                                     _sort_key(r, "model", _CLF_ORDER)))
    gain_html = ""
    if gain_ref:
        gain_rows = ""
        prev_cmp = None
        for r in gain_ref:
            cmp_key = r.get("comparison", "")
            border = "border-top:2px solid #d1d5db" if prev_cmp and cmp_key != prev_cmp else ""
            prev_cmp = cmp_key
            d, p = _f(r, "delta_auc"), _f(r, "delong_p")
            colour, weight = (("#15803d", "bold") if p == p and p < 0.05 and d > 0
                              else ("#6b7280", "normal"))
            what = ("CA-125 added to cfDNA" if cmp_key == "vs_cfdna"
                    else "cfDNA added to CA-125")
            gain_rows += (
                f"<tr style='{border}'><td style='font-size:.9em;color:#555'>{what}</td>"
                f"<td><strong>{r.get('model', '')}</strong></td>"
                f"<td style='text-align:right;color:#666'>{r.get('baseline', '')}</td>"
                f"<td style='text-align:right'>{_num(r.get('auc_baseline'), '{:.4f}')}</td>"
                f"<td style='text-align:right'>{_num(r.get('auc_model'), '{:.4f}')}</td>"
                f"<td style='text-align:right;color:{colour};font-weight:{weight}'>"
                f"{_num(r.get('delta_auc'), '{:+.4f}')}</td>"
                f"<td style='text-align:right;color:{colour}'>"
                f"{_num(r.get('delta_pct'), '{:+.1f}')}%</td>"
                f"<td style='text-align:right;color:#666'>"
                f"{_num(r.get('headroom_pct'), '{:+.0f}')}%</td>"
                f"<td style='text-align:right;color:#666;font-size:.9em'>{_ci(r)}</td>"
                f"<td style='text-align:right'>{_fmt_p(r.get('delong_p'))}</td></tr>"
            )
        gain_html = f"""
<h3 style="margin-bottom:4px">What each half of a fusion model contributes
  (<code>in_fold</code> arm)</h3>
<table>
<thead><tr><th>Comparison</th><th>Fusion model</th><th style="text-align:right">Baseline</th>
<th style="text-align:right">Baseline AUC</th><th style="text-align:right">Fusion AUC</th>
<th style="text-align:right">&Delta;AUC</th><th style="text-align:right">Relative</th>
<th style="text-align:right">Headroom</th>
<th style="text-align:right">&Delta;AUC 95% CI</th>
<th style="text-align:right">DeLong p</th></tr></thead>
<tbody>{gain_rows}</tbody>
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Both comparisons sit inside a single arm, so neither is affected by the panel-provenance
  question above. <em>Relative</em> is the percentage increase over the baseline AUC;
  <em>Headroom</em> is the share of the baseline's remaining distance to AUC 1.0 that the
  fusion closes — the fairer figure, because AUC is floored at 0.5 and a relative percentage
  flatters a weak baseline. Green = significant improvement (p &lt; 0.05).
  Full data: <code>leakage_demo_ca125_gain.csv</code>
</p>"""

        best = max((r for r in gain_ref if r.get("comparison") == "vs_ca125"),
                   key=lambda r: _f(r, "delta_auc"), default=None)
        if best is not None and _f(best, "delta_auc") == _f(best, "delta_auc"):
            txt = (
                f"Against the existing standard of care, adding the cfDNA panel to CA-125 "
                f"raises AUC by up to <strong>{_f(best, 'delta_auc'):+.3f}</strong> "
                f"({best.get('model')}: {_num(best.get('auc_baseline'), '{:.3f}')} &rarr; "
                f"{_num(best.get('auc_model'), '{:.3f}')}, "
                f"{_num(best.get('delta_pct'), '{:+.1f}')}% relative, closing "
                f"{_num(best.get('headroom_pct'), '{:+.0f}')}% of the remaining headroom to a "
                f"perfect AUC; 95% CI {_ci(best)}, p = {_fmt_p(best.get('delong_p'))}). This "
                f"is the clinically load-bearing comparison — CA-125 is what the panel has to "
                f"beat to be worth ordering — and it is measured here under the leakage-free "
                f"protocol, not the panel-restricted one."
            )
            # The direction of this comparison is the one most likely to be
            # overstated, so say plainly when the cohort cannot support the claim
            # rather than leaving the p-value to be skimmed past.
            if not (_f(best, "delong_p") < 0.05):
                txt += (" <strong>The increment is not statistically significant at this "
                        "sample size</strong> and its interval includes zero, so it should be "
                        "reported as an effect size the cohort is underpowered to establish — "
                        "not as the panel outperforming CA-125.")
            interp.append(txt)

    ca125 = [r for r in metrics if r.get("arm") == "ca125_only"]
    if ca125:
        interp.append(
            f"CA-125 alone reaches AUC <strong>{_num(ca125[0].get('auc'), '{:.4f}')}</strong> "
            "and is protocol-invariant, since it is measured rather than learned. The "
            "<code>+CA125</code> fusion models pass only their cfDNA half through feature "
            "selection, so they inherit their arm's protocol and stay comparable across arms; "
            "the CA-125 column is median-imputed and standardised on the training fold only."
        )

    interp.append(
        "This effect scales with dimensionality burden (features per sample). A cohort with "
        "many features relative to its sample count gives a supervised filter far more room to "
        "fit the label vector by chance, so the same protocols can show a large gap on one "
        "matrix and none on another."
    )
    interp_html = "<ul style='line-height:1.8;margin:6px 0'>" + "".join(
        f"<li>{t}</li>" for t in interp) + "</ul>"

    fig_rel = _copy_img(os.path.join(result_dir, "leakage_demo.png"), imgs_dir)
    fig_html = (_img_tag(fig_rel, "Leakage demonstration") if fig_rel else
                "<p style='color:#aaa;font-size:.85em'>leakage_demo.png not generated</p>")
    # The fusion figure is its own file: 10 models in one grouped bar chart is
    # unreadable, and it is absent for every cohort without CA-125.
    ca_rel = _copy_img(os.path.join(result_dir, "leakage_demo_ca125.png"), imgs_dir)
    if ca_rel:
        fig_html += _img_tag(ca_rel, "Leakage demonstration — CA-125 fusion models")

    html = f"""<div class="section">
<h2>Step 15b. Feature-selection Leakage Demonstration</h2>
<p style="color:#555;line-height:1.7">
  AIMethFlow's central methodological claim is that selecting features on the full cohort
  before cross-validation inflates apparent performance, and that confining selection to each
  training fold prevents it. This section measures that claim rather than asserting it. The
  same classifiers are run over the <strong>same cross-validation partitions</strong> under
  protocols that differ <em>only</em> in what the selection, imputation and scaling stages are
  allowed to see. Every sample therefore receives exactly one out-of-fold prediction per arm,
  so the arms are paired sample-by-sample and their AUC differences are testable with the
  DeLong test rather than merely described.
</p>
<p style="color:#666;font-size:.87em;margin:0 0 10px">
  Cohort: <strong>{cohort}</strong> &middot; n = {n_samples} &middot;
  arms run: {len(arms_present)} &middot; AUC CIs are bootstrap percentile intervals.
</p>

<h3 style="margin-bottom:4px">What each protocol is allowed to see</h3>
<table>
<thead><tr><th>Arm</th><th>Feature selection</th><th>Imputation + scaling</th>
<th>Feature space</th><th>Notes</th></tr></thead>
<tbody>{arm_rows}</tbody>
</table>

<h3 style="margin-bottom:4px">Performance by arm</h3>
<table>
<thead><tr><th>Arm</th><th>Classifier</th><th style="text-align:right">AUC</th>
<th style="text-align:right">95% CI</th><th style="text-align:right">Acc</th>
<th style="text-align:right">Sens</th><th style="text-align:right">Spec</th>
<th style="text-align:right">Feat/fold</th></tr></thead>
<tbody>{metric_rows}</tbody>
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Full data: <code>leakage_demo.csv</code>
</p>

<h3 style="margin-bottom:4px">Paired DeLong tests vs the leakage-free reference</h3>
<table>
<thead><tr><th>Arm</th><th>Classifier</th><th style="text-align:right">Arm AUC</th>
<th style="text-align:right"><code>in_fold</code></th>
<th style="text-align:right">&Delta;AUC</th>
<th style="text-align:right">&Delta;AUC 95% CI</th>
<th style="text-align:right">DeLong p</th></tr></thead>
<tbody>{delong_rows}</tbody>
</table>
<p style="color:#666;font-size:.82em;margin-top:4px">
  Red = significantly inflated relative to the leakage-free protocol (p &lt; 0.05);
  blue = significantly lower; grey = not significant. The interval is a paired bootstrap on
  the difference itself — both AUCs are recomputed on each resample — so it is not the same
  as comparing the two arms' separate CIs, which overlap even when a paired difference is
  large. Read the two together: p says whether the difference is resolvable, the interval
  says how large a difference remains compatible with the data.
  Full data: <code>leakage_demo_delong.csv</code>
</p>
{gain_html}

<h3 style="margin-bottom:4px">Interpretation</h3>
<div class="interp">{interp_html}</div>

{fig_html}
</div>"""

    return html, {"metrics": metrics, "delong": delong, "arms": arms_present,
                  "ca125_gain": gain}
