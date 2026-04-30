#!/usr/bin/env python3
"""
test_parse_sample_info.py - unit + integration tests for parse_sample_info

Run (from project root):
    pytest src/test/test_parse_sample_info.py -v
    python src/test/test_parse_sample_info.py
"""

import gc
import os
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock

import pandas as pd

_SRC = os.path.join(os.path.dirname(__file__), "..")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

_mock_config = MagicMock(return_value={})
with patch("utils.config.load_config", _mock_config):
    from methylation.sample_grouping import parse_sample_info


# ── helpers ───────────────────────────────────────────────────────────────────

def _make_excel(sheets: dict, suffix=".xlsx") -> str:
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.close()
    with pd.ExcelWriter(tmp.name, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    return tmp.name


def _unlink_safe(path):
    gc.collect()
    try:
        os.unlink(path)
    except PermissionError:
        pass


def _basic_df(rows):
    """Build a DataFrame with columns: fq_R1, fq_R2, category, subtype."""
    return pd.DataFrame(rows, columns=["fq_R1", "fq_R2", "category", "subtype"])


# ── unit tests ────────────────────────────────────────────────────────────────

class TestParseSampleInfoUnit(unittest.TestCase):

    def test_file_not_found_returns_none(self):
        result = parse_sample_info("/nonexistent/path/sample.xlsx")
        self.assertIsNone(result)

    def test_sheet_without_category_skipped(self):
        df = pd.DataFrame({"fq_R1": ["A.R1.fq.gz"], "fq_R2": ["A.R2.fq.gz"]})
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            self.assertEqual(result["benign"]["R1"], [])
            self.assertEqual(result["malignant"]["R1"], [])
        finally:
            _unlink_safe(path)

    def test_benign_malignant_split(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "benign",    ""),
            ("B002.R1.fq.gz", "B002.R2.fq.gz", "benign",    ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "malignant", "high-grade serous carcinoma"),
        ])
        path = _make_excel({"plasma": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            self.assertEqual(len(result["benign"]["R1"]), 2)
            self.assertEqual(len(result["malignant"]["R1"]), 1)
            self.assertIn("B001.R1.fq.gz", result["benign"]["R1"])
            self.assertIn("M001.R1.fq.gz", result["malignant"]["R1"])
        finally:
            _unlink_safe(path)

    def test_borderline_excluded(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "benign",     ""),
            ("J001.R1.fq.gz", "J001.R2.fq.gz", "borderline", "mucinous"),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "malignant",  "clear cell carcinoma"),
        ])
        path = _make_excel({"tissue": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            all_r1 = result["benign"]["R1"] + result["malignant"]["R1"]
            self.assertNotIn("J001.R1.fq.gz", all_r1)
            self.assertEqual(len(result["benign"]["R1"]), 1)
            self.assertEqual(len(result["malignant"]["R1"]), 1)
        finally:
            _unlink_safe(path)

    def test_deduplication(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "benign",    ""),
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "benign",    ""),   # duplicate
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "malignant", "endometrioid carcinoma"),
        ])
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(result["benign"]["R1"].count("B001.R1.fq.gz"), 1)
        finally:
            _unlink_safe(path)

    def test_multiple_sheets_merged(self):
        df1 = _basic_df([("B001.R1.fq.gz", "B001.R2.fq.gz", "benign", "")])
        df2 = _basic_df([
            ("B002.R1.fq.gz", "B002.R2.fq.gz", "benign",    ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "malignant", "low-grade serous carcinoma"),
        ])
        empty = pd.DataFrame()
        path = _make_excel({"sheet1": df1, "sheet2": df2, "Sheet3": empty})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 2)
            self.assertEqual(len(result["malignant"]["R1"]), 1)
            self.assertIn("B001.R1.fq.gz", result["benign"]["R1"])
            self.assertIn("B002.R1.fq.gz", result["benign"]["R1"])
        finally:
            _unlink_safe(path)

    def test_missing_fq_r1_column(self):
        df = pd.DataFrame({"category": ["benign", "malignant"], "subtype": ["", ""]})
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            self.assertEqual(result["benign"]["R1"], [])
            self.assertEqual(result["malignant"]["R1"], [])
        finally:
            _unlink_safe(path)

    def test_category_with_whitespace(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", " benign ",   ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "malignant",  ""),
        ])
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 1,
                             "Leading/trailing whitespace in 'benign' should be stripped")
        finally:
            _unlink_safe(path)

    def test_nan_fq_r1_excluded(self):
        df = pd.DataFrame({
            "fq_R1":    ["B001.R1.fq.gz", None],
            "fq_R2":    ["B001.R2.fq.gz", None],
            "category": ["benign",         "benign"],
            "subtype":  ["",               ""],
        })
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 1)
        finally:
            _unlink_safe(path)

    def test_return_structure(self):
        df = _basic_df([("B001.R1.fq.gz", "B001.R2.fq.gz", "benign", "")])
        path = _make_excel({"s": df})
        try:
            result = parse_sample_info(path)
            self.assertIn("benign", result)
            self.assertIn("malignant", result)
            self.assertIn("R1", result["benign"])
            self.assertIn("R2", result["benign"])
            self.assertIn("R1", result["malignant"])
            self.assertIn("R2", result["malignant"])
        finally:
            _unlink_safe(path)


# ── integration tests (require real Excel file) ───────────────────────────────

_REAL_EXCEL = os.path.join(
    os.path.dirname(__file__), "..", "..", "data",
    "ovarian_fq_series_samples_v3.xlsx"
)


@unittest.skipUnless(os.path.exists(_REAL_EXCEL), "Real Excel file not found, skipping integration tests")
class TestParseSampleInfoIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = parse_sample_info(_REAL_EXCEL)

    def test_result_not_none(self):
        self.assertIsNotNone(self.result)

    def test_has_benign_samples(self):
        self.assertGreater(len(self.result["benign"]["R1"]), 0)

    def test_has_malignant_samples(self):
        self.assertGreater(len(self.result["malignant"]["R1"]), 0)

    def test_no_borderline_in_groups(self):
        xl = pd.ExcelFile(_REAL_EXCEL)
        borderline_r1 = set()
        for sheet in xl.sheet_names:
            df = pd.read_excel(_REAL_EXCEL, sheet_name=sheet)
            if "category" not in df.columns or "fq_R1" not in df.columns:
                continue
            df["category"] = df["category"].astype(str).str.strip()
            borderline_r1.update(
                df.loc[df["category"] == "borderline", "fq_R1"].dropna().tolist()
            )
        all_r1 = set(self.result["benign"]["R1"]) | set(self.result["malignant"]["R1"])
        leaked = borderline_r1 & all_r1
        self.assertEqual(len(leaked), 0, f"Borderline samples leaked into groups: {leaked}")

    def test_no_duplicates_in_benign(self):
        r1_list = self.result["benign"]["R1"]
        self.assertEqual(len(r1_list), len(set(r1_list)))

    def test_no_duplicates_in_malignant(self):
        r1_list = self.result["malignant"]["R1"]
        self.assertEqual(len(r1_list), len(set(r1_list)))

    def test_r1_r2_same_length(self):
        b_r1, b_r2 = self.result["benign"]["R1"], self.result["benign"]["R2"]
        m_r1, m_r2 = self.result["malignant"]["R1"], self.result["malignant"]["R2"]
        self.assertEqual(len(b_r1), len(b_r2))
        self.assertEqual(len(m_r1), len(m_r2))

    def test_known_counts(self):
        # tissue: benign=49, malignant=74, borderline=12
        # plasma: benign=58, malignant=62, borderline=14
        self.assertGreaterEqual(len(self.result["benign"]["R1"]), 49)
        self.assertGreaterEqual(len(self.result["malignant"]["R1"]), 62)


if __name__ == "__main__":
    unittest.main(verbosity=2)
