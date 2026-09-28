"""build_gpb_supplementary_tables.py — generate the GPB submission's missing
Supplementary Tables (S3 Part B, S4, S5, S6) from existing result CSVs and
append them as new sheets to docs/GPB/resources/supplementary_tables.xlsx.

Run once from the project root:
    python src/scripts/build_gpb_supplementary_tables.py
"""

import os

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULT_DIR = os.path.join(PROJECT_ROOT, "result", "20260406")
EXTERNAL_DIR = os.path.join(PROJECT_ROOT, "result", "external")
CONCORDANCE_DIR = os.path.join(RESULT_DIR, "panel_concordance")
XLSX_PATH = os.path.join(PROJECT_ROOT, "docs", "GPB", "resources", "supplementary_tables.xlsx")

HEADER_FONT = Font(bold=True, color="FFFFFF", size=11)
HEADER_FILL = PatternFill(start_color="34495E", end_color="34495E", fill_type="solid")
HEADER_ALIGN = Alignment(horizontal="center", vertical="center", wrap_text=True)
ZEBRA_FILL = PatternFill(start_color="F4F6F8", end_color="F4F6F8", fill_type="solid")


def write_sheet(wb, title: str, df: pd.DataFrame) -> None:
    """Write ``df`` as a new (or replaced) formatted sheet named ``title``."""
    if title in wb.sheetnames:
        del wb[title]
    ws = wb.create_sheet(title)

    headers = list(df.columns)
    ws.append(headers)
    for c in range(1, len(headers) + 1):
        cell = ws.cell(row=1, column=c)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = HEADER_ALIGN

    for r_idx, row in enumerate(df.itertuples(index=False), start=2):
        for c_idx, value in enumerate(row, start=1):
            if isinstance(value, float) and np.isnan(value):
                value = None
            cell = ws.cell(row=r_idx, column=c_idx, value=value)
            cell.alignment = Alignment(horizontal="center", vertical="center")
            if r_idx % 2 == 0:
                cell.fill = ZEBRA_FILL

    ws.freeze_panes = "A2"
    ws.auto_filter.ref = ws.dimensions
    for c_idx, header in enumerate(headers, start=1):
        max_len = len(str(header))
        for r in range(2, min(ws.max_row + 1, 200)):
            v = ws.cell(r, c_idx).value
            if v is not None:
                max_len = max(max_len, len(str(v)))
        ws.column_dimensions[get_column_letter(c_idx)].width = min(max_len + 2, 40)


# ─────────────────────────────────────────────────────────────────────────
# Table S3 Part B — McNemar operating-point comparisons
# ─────────────────────────────────────────────────────────────────────────

def build_s3_partb() -> pd.DataFrame:
    path = os.path.join(RESULT_DIR, "pub_mcnemar_operating.csv")
    if not os.path.exists(path):
        path = os.path.join(RESULT_DIR, "imeta_mcnemar_operating.csv")
    df = pd.read_csv(path)
    df = df.rename(columns={
        "reference": "Reference",
        "comparator": "Comparator",
        "metric": "Metric",
        "stratum": "Stratum",
        "n": "N",
        "rate_reference": "Rate (Reference)",
        "rate_comparator": "Rate (Comparator)",
        "delta": "Delta",
        "discordant_ref_only": "Discordant (Ref only)",
        "discordant_cmp_only": "Discordant (Cmp only)",
        "p_mcnemar_exact": "p (exact McNemar)",
        "significant": "Significant",
    })
    return df


# ─────────────────────────────────────────────────────────────────────────
# Table S4 — complete four-arm leakage benchmark (internal cohort)
# ─────────────────────────────────────────────────────────────────────────

def build_s4() -> pd.DataFrame:
    perf = pd.read_csv(os.path.join(RESULT_DIR, "leakage_demo.csv"))
    delong = pd.read_csv(os.path.join(RESULT_DIR, "leakage_demo_delong.csv"))

    perf = perf.rename(columns={
        "arm": "Arm",
        "classifier": "Classifier",
        "n_samples": "N",
        "n_features_mean": "Mean Features Selected",
        "auc": "AUC",
        "auc_ci_low": "AUC CI Low",
        "auc_ci_high": "AUC CI High",
        "accuracy": "Accuracy",
        "sensitivity": "Sensitivity",
        "specificity": "Specificity",
    })
    perf = perf.drop(columns=["cohort"], errors="ignore")

    delong = delong.rename(columns={
        "classifier": "Classifier",
        "arm": "Arm",
        "auc_arm": "AUC (Arm)",
        "arm_reference": "Reference Arm",
        "auc_reference": "AUC (Reference)",
        "delta_auc": "Delta AUC",
        "delta_ci_low": "Delta AUC CI Low",
        "delta_ci_high": "Delta AUC CI High",
        "delong_z": "DeLong Z",
        "delong_p": "DeLong p",
    })
    delong = delong.drop(columns=["cohort"], errors="ignore")

    merged = perf.merge(delong, on=["Classifier", "Arm"], how="left")
    order = ["in_fold", "naive_selection", "naive_full", "panel_then_in_fold", "ca125_only"]
    merged["_order"] = merged["Arm"].apply(lambda a: order.index(a) if a in order else len(order))
    merged = merged.sort_values(["Classifier", "_order"]).drop(columns=["_order"]).reset_index(drop=True)
    return merged


# ─────────────────────────────────────────────────────────────────────────
# Table S5 — complete external cohort results (+ per-histology breakdown)
# ─────────────────────────────────────────────────────────────────────────

def build_s5_arms() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(RESULT_DIR, "external_cohort_summary.csv"))
    df = df.rename(columns={
        "cohort": "Cohort",
        "arm": "Arm",
        "classifier": "Classifier",
        "n_samples": "N",
        "n_features_mean": "Mean Features",
        "auc": "AUC",
        "auc_ci_low": "AUC CI Low",
        "auc_ci_high": "AUC CI High",
        "accuracy": "Accuracy",
        "sensitivity": "Sensitivity",
        "specificity": "Specificity",
        "delta_auc_vs_full": "Delta AUC vs Full Matrix",
        "delong_p": "DeLong p",
    })
    return df


def _per_histology_auc(cohort_dir: str) -> pd.DataFrame:
    """One-vs-benign AUC per histology label, per classifier, for one external run."""
    from sklearn.metrics import roc_auc_score

    pred_path = os.path.join(cohort_dir, "classifier_predictions.csv")
    manifest_path = os.path.join(cohort_dir, "external_cohort_manifest.csv")
    if not (os.path.exists(pred_path) and os.path.exists(manifest_path)):
        return pd.DataFrame()

    pred = pd.read_csv(pred_path)
    manifest = pd.read_csv(manifest_path)
    merged = pred.merge(manifest[["sample_id", "diagnosis"]], on="sample_id", how="left")

    benign_labels = {"Benign", "benign"}
    rows = []
    for classifier, grp in merged.groupby("classifier"):
        benign_mask = grp["diagnosis"].isin(benign_labels)
        benign_probs = grp.loc[benign_mask, "prob_malignant"]
        for histology, hgrp in grp.groupby("diagnosis"):
            if histology in benign_labels:
                continue
            if len(hgrp) < 3:
                continue
            y = np.concatenate([np.zeros(len(benign_probs)), np.ones(len(hgrp))])
            scores = np.concatenate([benign_probs.values, hgrp["prob_malignant"].values])
            if len(np.unique(y)) < 2:
                continue
            try:
                auc = roc_auc_score(y, scores)
            except ValueError:
                continue
            rows.append({
                "Histology": histology,
                "Classifier": classifier,
                "N (histology)": len(hgrp),
                "N (benign)": len(benign_probs),
                "AUC vs benign": round(auc, 4),
            })
    return pd.DataFrame(rows)


def build_s5_histology() -> pd.DataFrame:
    df = _per_histology_auc(os.path.join(EXTERNAL_DIR, "ovaprint"))
    if df.empty:
        return df
    df.insert(0, "Cohort", "OvaPrint")
    df = df.sort_values(["Histology", "Classifier"]).reset_index(drop=True)
    return df


# ─────────────────────────────────────────────────────────────────────────
# Table S6 — panel direction concordance against GSE133556
# ─────────────────────────────────────────────────────────────────────────

def build_s6_sites() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(CONCORDANCE_DIR, "panel_concordance_sites.csv"))
    keep = [
        "feature", "hg19", "epic_probe", "match_distance_bp",
        "plasma_benign_mean", "plasma_malignant_mean", "plasma_delta", "plasma_p",
        "tissue_control_mean", "tissue_case_mean", "tissue_delta", "tissue_p",
        "concordant",
    ]
    df = df[[c for c in keep if c in df.columns]].rename(columns={
        "feature": "Site (hg38)",
        "hg19": "Site (hg19, liftOver)",
        "epic_probe": "Matched EPIC Probe",
        "match_distance_bp": "Match Distance (bp)",
        "plasma_benign_mean": "Plasma Mean (Benign)",
        "plasma_malignant_mean": "Plasma Mean (Malignant)",
        "plasma_delta": "Plasma Delta",
        "plasma_p": "Plasma p",
        "tissue_control_mean": "Tissue Mean (Control)",
        "tissue_case_mean": "Tissue Mean (Carcinoma)",
        "tissue_delta": "Tissue Delta",
        "tissue_p": "Tissue p",
        "concordant": "Direction Concordant",
    })
    return df


def build_s6_summary() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(CONCORDANCE_DIR, "panel_concordance_summary.csv"))
    df = df.rename(columns={
        "window_bp": "Matching Window (bp)",
        "level": "Level",
        "n_panel_sites": "N Panel Sites",
        "n_matched_sites": "N Matched",
        "pct_matched_sites": "% Matched",
        "n_testable": "N Testable",
        "n_concordant": "N Concordant",
        "pct_concordant": "% Concordant",
        "pct_expected_null": "% Expected (Null)",
        "binom_p_vs_chance": "Binomial p (vs 50%)",
        "binom_p_vs_null": "Binomial p (vs measured null)",
    })
    return df


def main() -> None:
    wb = load_workbook(XLSX_PATH)

    s3b = build_s3_partb()
    write_sheet(wb, "Table_S3_mcnemar", s3b)
    print(f"Table S3 Part B (McNemar): {len(s3b)} rows")

    s4 = build_s4()
    write_sheet(wb, "Table_S4_leakage_benchmark", s4)
    print(f"Table S4 (leakage benchmark): {len(s4)} rows")

    s5_arms = build_s5_arms()
    write_sheet(wb, "Table_S5_external_arms", s5_arms)
    print(f"Table S5 (external cohort arms): {len(s5_arms)} rows")

    s5_hist = build_s5_histology()
    if not s5_hist.empty:
        write_sheet(wb, "Table_S5_histology", s5_hist)
        print(f"Table S5 (per-histology breakdown): {len(s5_hist)} rows")
    else:
        print("Table S5 (per-histology breakdown): SKIPPED — source data not found")

    s6_sites = build_s6_sites()
    write_sheet(wb, "Table_S6_concordance_sites", s6_sites)
    print(f"Table S6 (concordance, per-site): {len(s6_sites)} rows")

    s6_summary = build_s6_summary()
    write_sheet(wb, "Table_S6_concordance_summary", s6_summary)
    print(f"Table S6 (concordance, per-window summary): {len(s6_summary)} rows")

    wb.save(XLSX_PATH)
    print(f"\nSaved: {XLSX_PATH}")
    print(f"Sheets: {wb.sheetnames}")


if __name__ == "__main__":
    main()
