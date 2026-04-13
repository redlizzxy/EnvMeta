"""堆叠图分析器端到端 smoke 测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")  # 非交互后端，测试用

from envmeta.analysis import stackplot
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def phylum_and_meta():
    ab = pd.read_csv(SAMPLE / "Phylum.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ab, md


def test_stackplot_sample_style(phylum_and_meta):
    ab, md = phylum_and_meta
    r = stackplot.analyze(ab, md, {"style": "sample", "top_n": 10})
    assert r.figure is not None
    assert r.stats is not None
    # 列数 = 样本数（10）
    assert r.stats.shape[1] == 10
    # 每列百分比之和约为 100
    assert (r.stats.sum(axis=0).round(1) == 100.0).all()


def test_stackplot_group_style(phylum_and_meta):
    ab, md = phylum_and_meta
    r = stackplot.analyze(ab, md, {"style": "group", "top_n": 10})
    # 列 = 3 个分组 (CK/A/B)
    assert list(r.stats.columns) == ["CK", "A", "B"]
    assert r.stats.shape[1] == 3


def test_stackplot_topn_plus_others(phylum_and_meta):
    """用 Genus（88 属）确认会产生 Others 行。"""
    ab = pd.read_csv(SAMPLE / "Genus.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    r = stackplot.analyze(ab, md, {"style": "group", "top_n": 10})
    # Top-10 + Others = 11 行
    assert r.stats.shape[0] == 11
    assert "Others" in r.stats.index


def test_export_to_bytes(phylum_and_meta):
    ab, md = phylum_and_meta
    r = stackplot.analyze(ab, md, {"style": "sample", "top_n": 10})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
