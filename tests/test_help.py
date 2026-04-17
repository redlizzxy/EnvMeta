"""S8-ux 新手落地包 — 数据字典完整性测试。

验证三个数据模块（NAVIGATOR / INTERPRETATIONS / FILE_TO_ANALYSIS）的交叉引用
一致性，避免新增分析时漏填字典被用户从 UI 端踩到。
"""
from __future__ import annotations

from envmeta.file_manager.detector import FileType
from envmeta.help.file_analysis_map import (
    ANALYSIS_INPUTS, FILE_TO_ANALYSIS, analyses_ready,
)
from envmeta.help.interpretations import INTERPRETATIONS
from envmeta.help.research_navigator import NAVIGATOR, all_analysis_ids


# 预期全覆盖的 14 个分析 ID
EXPECTED_IDS = {
    "stackplot", "pcoa", "gene_heatmap", "alpha_boxplot", "log2fc", "rda", "lefse",
    "mag_quality", "mag_heatmap", "pathway", "gene_profile", "network",
    "cycle_diagram", "hypothesis_score",
}


# ══════════════════════════════════════════════════════════════
# ANALYSIS_INPUTS — 主数据源完整性
# ══════════════════════════════════════════════════════════════

def test_analysis_inputs_covers_14_analyses():
    assert set(ANALYSIS_INPUTS.keys()) == EXPECTED_IDS


def test_analysis_inputs_schema_valid():
    for aid, spec in ANALYSIS_INPUTS.items():
        assert "name" in spec and isinstance(spec["name"], str)
        assert "page" in spec and isinstance(spec["page"], str)
        assert "required" in spec and isinstance(spec["required"], list)
        assert "optional" in spec and isinstance(spec["optional"], list)
        for ft in spec["required"] + spec["optional"]:
            assert isinstance(ft, FileType), f"{aid}: {ft} 不是 FileType 枚举"


def test_analysis_inputs_no_duplicate_filetype():
    """同一分析不能同时把一个 FileType 放在 required 和 optional。"""
    for aid, spec in ANALYSIS_INPUTS.items():
        both = set(spec["required"]) & set(spec["optional"])
        assert not both, f"{aid}: {both} 同时出现在 required 和 optional"


# ══════════════════════════════════════════════════════════════
# INTERPRETATIONS — 14 覆盖 + schema
# ══════════════════════════════════════════════════════════════

def test_interpretations_covers_all_analyses():
    assert set(INTERPRETATIONS.keys()) == EXPECTED_IDS


def test_interpretations_schema_valid():
    required_keys = {"title", "what_it_shows", "how_to_read", "good_signal",
                     "warning", "caveats"}
    for aid, content in INTERPRETATIONS.items():
        missing = required_keys - set(content.keys())
        assert not missing, f"{aid}: 缺 {missing}"
        assert isinstance(content["how_to_read"], list) and content["how_to_read"]


def test_interpretations_content_length_reasonable():
    """每字段至少 10 字符，避免空串漏填。"""
    for aid, content in INTERPRETATIONS.items():
        for key in ("what_it_shows", "good_signal", "warning", "caveats"):
            assert len(content[key]) >= 10, f"{aid}.{key} 过短"


# ══════════════════════════════════════════════════════════════
# NAVIGATOR — 所有 analysis_id 合法
# ══════════════════════════════════════════════════════════════

def test_navigator_references_valid_analyses():
    nav_ids = all_analysis_ids()
    invalid = nav_ids - set(ANALYSIS_INPUTS.keys())
    assert not invalid, f"NAVIGATOR 引用了未注册的 analysis_id: {invalid}"


def test_navigator_schema_valid():
    for cat in NAVIGATOR:
        assert {"category", "icon", "description", "subquestions"} <= set(cat.keys())
        assert cat["subquestions"], f"{cat['category']} subquestions 为空"
        for sq in cat["subquestions"]:
            assert "q" in sq and "recommended" in sq and "required_files" in sq
            assert sq["recommended"], f"{cat['category']}/{sq['q']} 无推荐"
            for rec in sq["recommended"]:
                assert {"analysis_id", "reason", "priority"} <= set(rec.keys())
            # required_files 必须是合法 FileType.value
            valid_values = {ft.value for ft in FileType}
            for rf in sq["required_files"]:
                assert rf in valid_values, f"{cat['category']}/{sq['q']}: {rf} 非法"


def test_navigator_covers_most_analyses():
    """所有 Reads + MAG + 循环图分析都应该有至少一个 NAVIGATOR 入口。"""
    covered = all_analysis_ids()
    # 14 个分析里至少 12 个应该被 navigator 覆盖
    assert len(covered & EXPECTED_IDS) >= 12, \
        f"NAVIGATOR 只覆盖 {len(covered)} 个分析（应 ≥ 12）"


# ══════════════════════════════════════════════════════════════
# FILE_TO_ANALYSIS — 反向索引正确
# ══════════════════════════════════════════════════════════════

def test_file_to_analysis_all_filetypes_indexed():
    """FILE_TO_ANALYSIS 应为每个 FileType 生成 key（值可以是空列表）。"""
    for ft in FileType:
        assert ft in FILE_TO_ANALYSIS


def test_file_to_analysis_entries_consistent():
    """反向索引的内容和 ANALYSIS_INPUTS 的 required/optional 完全对应。"""
    for ft, entries in FILE_TO_ANALYSIS.items():
        for aid, name, role in entries:
            assert aid in ANALYSIS_INPUTS, f"反向索引里的 {aid} 不在 ANALYSIS_INPUTS"
            assert role in ("required", "optional")
            spec = ANALYSIS_INPUTS[aid]
            if role == "required":
                assert ft in spec["required"], f"{aid}-{ft}: 索引标 required 但实际不是"
            else:
                assert ft in spec["optional"]


def test_analyses_ready_detects_complete_sets():
    """给出所有必需 FileType，analyses_ready 应返回全部 14 个 analysis。"""
    all_ft = set(FileType)
    all_ft.discard(FileType.UNKNOWN)
    ready = set(analyses_ready(all_ft))
    assert ready == EXPECTED_IDS


def test_analyses_ready_partial_set():
    """只给 metadata + abundance_wide，应能跑 stackplot / lefse，但不能跑 pcoa（缺 distance）。"""
    partial = {FileType.METADATA, FileType.ABUNDANCE_WIDE}
    ready = set(analyses_ready(partial))
    assert "stackplot" in ready
    assert "lefse" in ready
    assert "pcoa" not in ready
    assert "rda" not in ready  # 缺 env_factors
    assert "mag_heatmap" in ready  # 只需 abundance_wide


def test_analyses_ready_empty_returns_nothing_useful():
    ready = analyses_ready(set())
    # 空文件集下只有不需要任何 required 的分析会返回（预期 0 个）
    assert len(ready) == 0
