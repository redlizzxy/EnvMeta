"""MAG 质量散点图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import mag_quality
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def mq_inputs():
    quality = pd.read_csv(SAMPLE / "quality_report.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    return quality, tax, ks


def test_mag_quality_smoke(mq_inputs):
    q, t, k = mq_inputs
    r = mag_quality.analyze(q, t, k)
    assert r.figure is not None
    assert not r.stats.empty
    assert {"Quality", "Count"}.issubset(r.stats.columns)


def test_mag_quality_summary_has_three_rows(mq_inputs):
    q, t, k = mq_inputs
    r = mag_quality.analyze(q, t, k)
    summary = r.stats[r.stats["type"] == "summary"]
    assert set(summary["Quality"]) == {"High", "Medium", "Low"}
    assert summary["Count"].sum() == len(q)


def test_mag_quality_no_taxonomy(mq_inputs):
    q, _, _ = mq_inputs
    r = mag_quality.analyze(q)
    # 无 taxonomy 时所有点归 Other
    detail = r.stats[r.stats["type"] == "detail"]
    assert (detail["Phylum"] == "Unknown").all()


def test_mag_quality_keystone_marked(mq_inputs):
    q, t, k = mq_inputs
    r = mag_quality.analyze(q, t, k)
    detail = r.stats[r.stats["type"] == "detail"]
    assert detail["is_keystone"].any()


def test_mag_quality_export(mq_inputs):
    q, t, k = mq_inputs
    r = mag_quality.analyze(q, t, k)
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
