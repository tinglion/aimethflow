"""split_supplementary_tables.py — split the bundled supplementary_tables.xlsx
into one standalone file per supplementary table (GPB requires "1 figure/table
in 1 file"). Sheets are grouped by table number, not one-sheet-per-file, since
the manuscript's own captions describe some tables as having multiple parts
sharing one table number (S3: Part A/B; S5: primary arms + per-histology; S6:
per-site + per-window summary) — those parts become tabs within that table's
one file, not separate table numbers.

Run once from the project root:
    python src/scripts/split_supplementary_tables.py
"""

import os
import shutil

from openpyxl import load_workbook

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SOURCE_XLSX = os.path.join(PROJECT_ROOT, "docs", "GPB", "resources", "supplementary_tables.xlsx")
OUT_DIR = os.path.join(PROJECT_ROOT, "docs", "GPB", "resources", "supplementary_tables_split")

# table_number -> (output filename, [(source_sheet_name, new_tab_name), ...])
TABLE_GROUPS = {
    "S1": ("Table_S1.xlsx", [("Table_S1_steps", "Steps")]),
    "S2": ("Table_S2.xlsx", [("Table_S2_annotation", "Annotation")]),
    "S3": ("Table_S3.xlsx", [("Table_S3_delong", "Part A - DeLong"),
                              ("Table_S3_mcnemar", "Part B - McNemar")]),
    "S4": ("Table_S4.xlsx", [("Table_S4_leakage_benchmark", "Leakage benchmark")]),
    "S5": ("Table_S5.xlsx", [("Table_S5_external_arms", "External arms"),
                              ("Table_S5_histology", "Per-histology")]),
    "S6": ("Table_S6.xlsx", [("Table_S6_concordance_sites", "Per-site"),
                              ("Table_S6_concordance_summary", "Per-window summary")]),
}


def main():
    src_wb = load_workbook(SOURCE_XLSX)
    os.makedirs(OUT_DIR, exist_ok=True)

    for table_num, (out_name, sheet_map) in TABLE_GROUPS.items():
        missing = [s for s, _ in sheet_map if s not in src_wb.sheetnames]
        if missing:
            print(f"{table_num}: SKIPPED — missing source sheet(s) {missing}")
            continue

        out_path = os.path.join(OUT_DIR, out_name)
        shutil.copy(SOURCE_XLSX, out_path)
        dst_wb = load_workbook(out_path)

        # Rename the tabs we want, in the requested order.
        for src_name, new_name in sheet_map:
            ws = dst_wb[src_name]
            ws.title = new_name

        # Drop every sheet from the copy that doesn't belong to this table.
        keep = {new_name for _, new_name in sheet_map}
        for ws in list(dst_wb.worksheets):
            if ws.title not in keep:
                del dst_wb[ws.title]

        # openpyxl requires at least one visible sheet and a valid active index.
        dst_wb.active = 0
        dst_wb.save(out_path)
        print(f"{table_num}: saved {out_path} with tabs {[n for _, n in sheet_map]}")

    print(f"\nDone. Files are in: {OUT_DIR}")
    print("Original bundled workbook left untouched at:", SOURCE_XLSX)


if __name__ == "__main__":
    main()
