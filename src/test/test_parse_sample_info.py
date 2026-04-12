#!/usr/bin/env python3
"""
test_parse_sample_info.py - parse_sample_info 单元 + 集成测试

运行方式（在项目根目录）:
    pytest src/test/test_parse_sample_info.py -v
    python src/test/test_parse_sample_info.py       # 无 pytest 时直接运行
"""

import io
import os
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock

import pandas as pd

# ── 路径修正：让 `from methylation.sample_grouping import ...` 可用 ──
_SRC = os.path.join(os.path.dirname(__file__), "..")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

# 在 import sample_grouping 之前 mock load_config，避免依赖 config.toml
_mock_config = MagicMock(return_value={})
with patch("utils.config.load_config", _mock_config):
    from methylation.sample_grouping import parse_sample_info


# ─────────────────────────────────────────────────────────────────────────────
#  辅助函数：构建临时 Excel 文件
# ─────────────────────────────────────────────────────────────────────────────

def _make_excel(sheets: dict, suffix=".xlsx") -> str:
    """
    将 sheets = {"sheet_name": dataframe, ...} 写入临时 Excel 文件，
    返回文件路径（调用方负责删除）。
    """
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.close()
    with pd.ExcelWriter(tmp.name, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    return tmp.name


def _unlink_safe(path):
    """
    安全删除临时文件。
    Windows 上 pandas ExcelFile 有时持有文件句柄，
    强制 GC 后再尝试删除；失败则静默跳过（进程退出时 OS 自动清理）。
    """
    import gc
    gc.collect()
    try:
        os.unlink(path)
    except PermissionError:
        pass


def _basic_df(rows):
    """快速构建含 fq_R1 / fq_R2 / 分类 / 亚型 的 DataFrame。"""
    return pd.DataFrame(rows, columns=["fq_R1", "fq_R2", "分类", "亚型"])


# ─────────────────────────────────────────────────────────────────────────────
#  单元测试
# ─────────────────────────────────────────────────────────────────────────────

class TestParseSampleInfoUnit(unittest.TestCase):

    # ── 1. 文件不存在 → 返回 None ──────────────────────────────────────────
    def test_file_not_found_returns_none(self):
        result = parse_sample_info("/nonexistent/path/sample.xlsx")
        self.assertIsNone(result)

    # ── 2. sheet 无 '分类' 列 → 跳过，返回空组 ────────────────────────────
    def test_sheet_without_fen_lei_skipped(self):
        df = pd.DataFrame({"fq_R1": ["A.R1.fq.gz"], "fq_R2": ["A.R2.fq.gz"]})
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            self.assertEqual(result["benign"]["R1"], [])
            self.assertEqual(result["malignant"]["R1"], [])
        finally:
            _unlink_safe(path)

    # ── 3. 良性 / 恶性 正确分组 ─────────────────────────────────────────────
    def test_benign_malignant_split(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", ""),
            ("B002.R1.fq.gz", "B002.R2.fq.gz", "良性", ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "恶性", "高级别浆液性癌"),
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

    # ── 4. 交界性样本不出现在任何分组 ──────────────────────────────────────
    def test_borderline_excluded(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", ""),
            ("J001.R1.fq.gz", "J001.R2.fq.gz", "交界性", "黏液性"),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "恶性", "透明细胞癌"),
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

    # ── 5. 重复行去重 ────────────────────────────────────────────────────────
    def test_deduplication(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", ""),
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", ""),   # 重复
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "恶性", "子宫内膜样癌"),
        ])
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(result["benign"]["R1"].count("B001.R1.fq.gz"), 1)
        finally:
            _unlink_safe(path)

    # ── 6. 多 sheet 数据合并 ─────────────────────────────────────────────────
    def test_multiple_sheets_merged(self):
        df1 = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", ""),
        ])
        df2 = _basic_df([
            ("B002.R1.fq.gz", "B002.R2.fq.gz", "良性", ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "恶性", "低级别浆液性癌"),
        ])
        empty = pd.DataFrame()   # 空 sheet（无列）
        path = _make_excel({"sheet1": df1, "sheet2": df2, "Sheet3": empty})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 2)
            self.assertEqual(len(result["malignant"]["R1"]), 1)
            self.assertIn("B001.R1.fq.gz", result["benign"]["R1"])
            self.assertIn("B002.R1.fq.gz", result["benign"]["R1"])
        finally:
            _unlink_safe(path)

    # ── 7. fq_R1 列缺失时不报错 ──────────────────────────────────────────────
    def test_missing_fq_r1_column(self):
        df = pd.DataFrame({"分类": ["良性", "恶性"], "亚型": ["", ""]})
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertIsNotNone(result)
            # 无 fq_R1 列时应安全返回空列表
            self.assertEqual(result["benign"]["R1"], [])
            self.assertEqual(result["malignant"]["R1"], [])
        finally:
            _unlink_safe(path)

    # ── 8. 分类列含前后空白应正确匹配 ──────────────────────────────────────
    def test_classification_with_whitespace(self):
        df = _basic_df([
            ("B001.R1.fq.gz", "B001.R2.fq.gz", " 良性 ", ""),
            ("M001.R1.fq.gz", "M001.R2.fq.gz", "恶性",   ""),
        ])
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 1,
                             "前后空白的 '良性' 应被正确识别")
        finally:
            _unlink_safe(path)

    # ── 9. NaN 行不进入分组 ──────────────────────────────────────────────────
    def test_nan_fq_r1_excluded(self):
        df = pd.DataFrame({
            "fq_R1":  ["B001.R1.fq.gz", None],
            "fq_R2":  ["B001.R2.fq.gz", None],
            "分类":   ["良性",           "良性"],
            "亚型":   ["",               ""],
        })
        path = _make_excel({"sheet1": df})
        try:
            result = parse_sample_info(path)
            self.assertEqual(len(result["benign"]["R1"]), 1)
        finally:
            _unlink_safe(path)

    # ── 10. 返回值结构正确 ───────────────────────────────────────────────────
    def test_return_structure(self):
        df = _basic_df([("B001.R1.fq.gz", "B001.R2.fq.gz", "良性", "")])
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


# ─────────────────────────────────────────────────────────────────────────────
#  集成测试（依赖真实 Excel 文件）
# ─────────────────────────────────────────────────────────────────────────────

_REAL_EXCEL = os.path.join(
    os.path.dirname(__file__), "..", "..", "data",
    "卵巢癌fq编号与样本信息对应表v3.xlsx"
)


@unittest.skipUnless(os.path.exists(_REAL_EXCEL), "真实 Excel 文件不存在，跳过集成测试")
class TestParseSampleInfoIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.result = parse_sample_info(_REAL_EXCEL)

    def test_result_not_none(self):
        self.assertIsNotNone(self.result)

    def test_has_benign_samples(self):
        self.assertGreater(len(self.result["benign"]["R1"]), 0,
                           "至少应有一个良性样本")

    def test_has_malignant_samples(self):
        self.assertGreater(len(self.result["malignant"]["R1"]), 0,
                           "至少应有一个恶性样本")

    def test_no_borderline_in_groups(self):
        """交界性样本（fq_R1 含关键字）不应出现在分组结果中。"""
        xl = pd.ExcelFile(_REAL_EXCEL)
        borderline_r1 = set()
        for sheet in xl.sheet_names:
            df = pd.read_excel(_REAL_EXCEL, sheet_name=sheet)
            if "分类" not in df.columns or "fq_R1" not in df.columns:
                continue
            df["分类"] = df["分类"].astype(str).str.strip()
            borderline_r1.update(
                df.loc[df["分类"] == "交界性", "fq_R1"].dropna().tolist()
            )

        all_r1 = set(self.result["benign"]["R1"]) | set(self.result["malignant"]["R1"])
        leaked = borderline_r1 & all_r1
        self.assertEqual(len(leaked), 0,
                         f"以下交界性样本泄漏进了分组：{leaked}")

    def test_no_duplicates_in_benign(self):
        r1_list = self.result["benign"]["R1"]
        self.assertEqual(len(r1_list), len(set(r1_list)), "良性 R1 列表中有重复")

    def test_no_duplicates_in_malignant(self):
        r1_list = self.result["malignant"]["R1"]
        self.assertEqual(len(r1_list), len(set(r1_list)), "恶性 R1 列表中有重复")

    def test_r1_r2_same_length(self):
        b_r1, b_r2 = self.result["benign"]["R1"], self.result["benign"]["R2"]
        m_r1, m_r2 = self.result["malignant"]["R1"], self.result["malignant"]["R2"]
        self.assertEqual(len(b_r1), len(b_r2), "良性 R1/R2 数量不一致")
        self.assertEqual(len(m_r1), len(m_r2), "恶性 R1/R2 数量不一致")

    def test_known_counts(self):
        """
        Excel 中 tissue sheet: 良性 49 / 恶性 74 / 交界性 12
                 plasma sheet: 良性 58 / 恶性 62 / 交界性 14
        合计（去重后）良性 ≥ 49，恶性 ≥ 62
        """
        self.assertGreaterEqual(len(self.result["benign"]["R1"]), 49)
        self.assertGreaterEqual(len(self.result["malignant"]["R1"]), 62)


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    unittest.main(verbosity=2)
