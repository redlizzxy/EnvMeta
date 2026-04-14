"""LEfSe 差异特征测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import lefse
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def lefse_inputs():
    ab = pd.read_csv(SAMPLE / "Genus.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ab, md


def test_lefse_smoke(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {"alpha_kw": 0.1, "lda_threshold": 2.0})
    assert r.figure is not None
    assert not r.stats.empty
    assert {"feature", "group", "lda", "kw_p"}.issubset(r.stats.columns)


def test_lefse_lda_above_threshold(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {"alpha_kw": 0.1, "lda_threshold": 2.5})
    assert (r.stats["lda"] >= 2.5).all()


def test_lefse_kw_below_alpha(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {"alpha_kw": 0.1, "lda_threshold": 2.0})
    assert (r.stats["kw_p"] < 0.1).all()


def test_lefse_group_assignment(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {"alpha_kw": 0.1, "lda_threshold": 2.0})
    # enriched 组应是 metadata 里的组之一
    valid = set(md["Group"].astype(str))
    assert set(r.stats["group"].astype(str)).issubset(valid)


def test_lefse_tax_level_filter(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {
        "alpha_kw": 0.1, "lda_threshold": 2.0,
        "tax_levels": ["Genus"],
    })
    assert (r.stats["tax_level"] == "Genus").all()


def test_lefse_export(lefse_inputs):
    ab, md = lefse_inputs
    r = lefse.analyze(ab, md, {"alpha_kw": 0.1, "lda_threshold": 2.0})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
