"""元素循环基因热图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import gene_heatmap
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def ko_and_meta():
    ko = pd.read_csv(SAMPLE / "ko_tpm.spf", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ko, md


def test_heatmap_smoke(ko_and_meta):
    ko, md = ko_and_meta
    r = gene_heatmap.analyze(ko, md)
    assert r.figure is not None
    assert r.stats is not None
    # 3 组（CK / A / B）
    assert r.stats.shape[1] == 3


def test_heatmap_present_kos(ko_and_meta):
    """ko_tpm.spf 有 51 个 KO，热图应覆盖 ≤ 57 且 > 30。"""
    ko, md = ko_and_meta
    r = gene_heatmap.analyze(ko, md)
    assert 30 < r.stats.shape[0] <= 57


def test_heatmap_element_filter(ko_and_meta):
    """只看 arsenic 时 KO 数应该大幅减少。"""
    ko, md = ko_and_meta
    r_all = gene_heatmap.analyze(ko, md)
    r_as = gene_heatmap.analyze(ko, md, {"element_filter": ["arsenic"]})
    assert r_as.stats.shape[0] < r_all.stats.shape[0]
    # 所有行对应的元素都应是 arsenic
    row_meta = r_as.params["_row_meta"]
    assert all(el == "arsenic" for el, _, _ in row_meta)


def test_heatmap_export_png(ko_and_meta):
    ko, md = ko_and_meta
    r = gene_heatmap.analyze(ko, md)
    png = export_to_bytes(r.figure, "png")
    assert png.startswith(b"\x89PNG")
