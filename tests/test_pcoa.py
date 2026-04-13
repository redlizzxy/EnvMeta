"""PCoA + PERMANOVA 测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import pcoa
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def bray_and_meta():
    dist = pd.read_csv(SAMPLE / "beta_bray.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return dist, md


def test_pcoa_smoke(bray_and_meta):
    dist, md = bray_and_meta
    r = pcoa.analyze(dist, md)
    assert r.figure is not None
    assert r.stats is not None
    assert not r.stats.empty


def test_pcoa_variance_explained_reasonable(bray_and_meta):
    """前两轴应至少解释 30% 总方差（Bray 距离 10 样本应该足够）。"""
    dist, md = bray_and_meta
    r = pcoa.analyze(dist, md)
    tables = r.params["_stats_tables"]
    coords = tables["coords"]
    explained = coords.shape[0]  # 轴数 ≤ 样本数
    pc1 = r.stats.loc[r.stats["metric"] == "PC1_explained", "value"].iloc[0]
    pc2 = r.stats.loc[r.stats["metric"] == "PC2_explained", "value"].iloc[0]
    assert (pc1 + pc2) > 0.3


def test_pcoa_permanova_significant(bray_and_meta):
    """三组微生物群落应该有显著差异（P < 0.1）。"""
    dist, md = bray_and_meta
    r = pcoa.analyze(dist, md, {"n_permutations": 499})
    p_val = r.stats.loc[r.stats["metric"] == "PERMANOVA_p_value", "value"].iloc[0]
    assert p_val < 0.1


def test_pcoa_export_png(bray_and_meta):
    dist, md = bray_and_meta
    r = pcoa.analyze(dist, md)
    png = export_to_bytes(r.figure, "png")
    assert png.startswith(b"\x89PNG")
