"""log2FC 差异柱图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import log2fc
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def ko_and_meta():
    ab = pd.read_csv(SAMPLE / "ko_tpm.spf", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ab, md


def test_log2fc_smoke(ko_and_meta):
    ab, md = ko_and_meta
    r = log2fc.analyze(ab, md, {"group_a": "A", "group_b": "CK"})
    assert r.figure is not None
    assert not r.stats.empty
    expected = {"ko", "gene", "pathway", "element",
                "mean_a", "mean_b", "log2fc", "t", "p", "padj",
                "significant", "significance"}
    assert expected.issubset(r.stats.columns)


def test_log2fc_padj_range(ko_and_meta):
    ab, md = ko_and_meta
    r = log2fc.analyze(ab, md, {"group_a": "B", "group_b": "CK"})
    assert r.stats["padj"].between(0, 1).all()
    # log2fc 可以为负/正，但应为有限数
    import numpy as np
    assert np.isfinite(r.stats["log2fc"]).all()


def test_log2fc_element_filter(ko_and_meta):
    ab, md = ko_and_meta
    r = log2fc.analyze(ab, md, {"group_a": "A", "group_b": "CK",
                                 "element_filter": ["arsenic", "sulfur"]})
    assert set(r.stats["element"]) == {"arsenic", "sulfur"}


def test_log2fc_export(ko_and_meta):
    ab, md = ko_and_meta
    r = log2fc.analyze(ab, md, {"group_a": "A", "group_b": "CK"})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")


def test_log2fc_same_group_rejected(ko_and_meta):
    ab, md = ko_and_meta
    with pytest.raises(ValueError):
        log2fc.analyze(ab, md, {"group_a": "A", "group_b": "A"})
