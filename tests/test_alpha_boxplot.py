"""α 多样性箱线图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import alpha_boxplot
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def alpha_and_meta():
    a = pd.read_csv(SAMPLE / "alpha.txt", sep="\t")
    m = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return a, m


def test_alpha_smoke(alpha_and_meta):
    a, m = alpha_and_meta
    r = alpha_boxplot.analyze(a, m)
    assert r.figure is not None
    assert r.stats is not None
    assert {"metric", "group_a", "group_b", "p", "padj", "significance"}.issubset(r.stats.columns)


def test_alpha_pair_count(alpha_and_meta):
    """5 指数 × 3 组对 = 15 行；padj 在 [0, 1]。"""
    a, m = alpha_and_meta
    r = alpha_boxplot.analyze(a, m)
    assert len(r.stats) == 5 * 3
    assert r.stats["padj"].between(0, 1).all()


def test_alpha_metrics_subset(alpha_and_meta):
    """显式指定 metrics 时只画和统计选中列。"""
    a, m = alpha_and_meta
    r = alpha_boxplot.analyze(a, m, {"metrics": ["shannon", "simpson"]})
    assert set(r.stats["metric"]) == {"shannon", "simpson"}


def test_alpha_export(alpha_and_meta):
    a, m = alpha_and_meta
    r = alpha_boxplot.analyze(a, m)
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
