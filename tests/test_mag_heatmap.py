"""MAG 丰度热图测试（S6）。"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import mag_heatmap
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def mh_inputs():
    ab = pd.read_csv(SAMPLE / "abundance.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ab, tax, ks, md


# ── 基础 ────────────────────────────────────────────────────

def test_mag_heatmap_smoke(mh_inputs):
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md)
    assert r.figure is not None
    assert not r.stats.empty
    assert {"MAG", "Phylum", "is_keystone", "row_order",
            "selection_score"}.issubset(r.stats.columns)


def test_top_n_filter_correct(mh_inputs):
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md, params={"top_n": 20})
    assert len(r.stats) == 20


def test_selection_by_mean_is_descending(mh_inputs):
    """selection_score 按选择标准降序（top_n 截断前）。"""
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"top_n": 30, "selection_by": "mean",
                                    "cluster_rows": False})
    # 关闭聚类时，row_order == 选择顺序 == score 降序
    s = r.stats["selection_score"].to_numpy()
    assert (s[:-1] >= s[1:] - 1e-9).all(), "selection_score 必须降序"


def test_cluster_order_deterministic(mh_inputs):
    """同输入 → 相同聚类顺序（凸显可重现性）。"""
    ab, tax, ks, md = mh_inputs
    r1 = mag_heatmap.analyze(ab, tax, ks, md, params={"top_n": 30})
    r2 = mag_heatmap.analyze(ab, tax, ks, md, params={"top_n": 30})
    assert r1.stats["MAG"].tolist() == r2.stats["MAG"].tolist()


# ── 可选输入 ────────────────────────────────────────────────

def test_works_without_taxonomy(mh_inputs):
    ab, _, _, _ = mh_inputs
    r = mag_heatmap.analyze(ab, params={"top_n": 15})
    assert (r.stats["Phylum"] == "Unknown").all()
    assert not r.stats["is_keystone"].any()


def test_works_without_metadata(mh_inputs):
    ab, tax, ks, _ = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, params={"top_n": 15})
    assert r.figure is not None
    assert len(r.stats) == 15


def test_phylum_bar_toggle_does_not_crash(mh_inputs):
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"top_n": 12, "show_phylum_bar": False,
                                    "show_group_bar": False})
    assert r.figure is not None


# ── 配色 ────────────────────────────────────────────────────

def test_three_stage_colormap_breakpoints_valid(mh_inputs):
    """自定义 breakpoints 不崩。"""
    ab, tax, _, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, params={"top_n": 10,
                                             "color_breakpoints": (0.1, 0.4)})
    assert r.figure is not None


def test_invalid_breakpoints_still_runs(mh_inputs):
    """极端 breakpoints（两值相等）应 fallback 正常运行。"""
    ab, _, _, _ = mh_inputs
    r = mag_heatmap.analyze(ab, params={"top_n": 5,
                                         "color_breakpoints": (0.0, 0.5)})
    assert r.figure is not None


# ── 错误输入 ────────────────────────────────────────────────

def test_empty_abundance_raises():
    empty = pd.DataFrame({"Genome": [], "S1": []})
    with pytest.raises(ValueError, match="无有效 MAG"):
        mag_heatmap.analyze(empty)


def test_missing_sample_cols_raises():
    only_mag = pd.DataFrame({"Genome": ["m1", "m2"]})
    with pytest.raises(ValueError, match="样本列"):
        mag_heatmap.analyze(only_mag)


# ── 导出 ────────────────────────────────────────────────────

def test_export_png_pdf(mh_inputs):
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md, params={"top_n": 10})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
