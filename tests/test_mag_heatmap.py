"""MAG 丰度热图测试（S6）。"""
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt  # noqa: F401 — used by legend assertion
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
    # filter_mode="top_n" → 精确 Top-N（不与 keystone 取并集）
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 20})
    assert len(r.stats) == 20


def test_selection_by_mean_is_descending(mh_inputs):
    """metric_desc 行排序下 selection_score 应降序。"""
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 30,
                                    "top_n_by": "mean",
                                    "row_order": "metric_desc"})
    s = r.stats["selection_score"].to_numpy()
    assert (s[:-1] >= s[1:] - 1e-9).all(), "selection_score 必须降序"


def test_cluster_order_deterministic(mh_inputs):
    """同输入 → 相同聚类顺序（凸显可重现性）。"""
    ab, tax, ks, md = mh_inputs
    r1 = mag_heatmap.analyze(ab, tax, ks, md,
                             params={"filter_mode": "top_n", "top_n_count": 30})
    r2 = mag_heatmap.analyze(ab, tax, ks, md,
                             params={"filter_mode": "top_n", "top_n_count": 30})
    assert r1.stats["MAG"].tolist() == r2.stats["MAG"].tolist()


# ── 可选输入 ────────────────────────────────────────────────

def test_works_without_taxonomy(mh_inputs):
    ab, _, _, _ = mh_inputs
    r = mag_heatmap.analyze(ab, params={"filter_mode": "top_n",
                                         "top_n_count": 15})
    assert (r.stats["Phylum"] == "Unknown").all()
    assert not r.stats["is_keystone"].any()


def test_works_without_metadata(mh_inputs):
    ab, tax, ks, _ = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks,
                            params={"filter_mode": "top_n", "top_n_count": 15})
    assert r.figure is not None
    assert len(r.stats) == 15


def test_phylum_bar_toggle_does_not_crash(mh_inputs):
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 12,
                                    "show_phylum_bar": False,
                                    "show_group_bar": False})
    assert r.figure is not None


# ── 配色 ────────────────────────────────────────────────────

def test_three_stage_colormap_breakpoints_valid(mh_inputs):
    ab, tax, _, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax,
                            params={"filter_mode": "top_n", "top_n_count": 10,
                                    "color_breakpoints": (0.1, 0.4)})
    assert r.figure is not None


def test_invalid_breakpoints_still_runs(mh_inputs):
    ab, _, _, _ = mh_inputs
    r = mag_heatmap.analyze(ab,
                            params={"filter_mode": "top_n", "top_n_count": 5,
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
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 10})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")


# ── S6-fix: Genus 标签 + 图例 ────────────────────────────────

def test_mag_label_uses_genus_when_tax_provided(mh_inputs):
    """上传 taxonomy 后 stats 里 label 列应包含 Genus（不是纯 Mx_All_XX）。"""
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 30})
    assert "label" in r.stats.columns
    not_equal = (r.stats["label"] != r.stats["MAG"]).sum()
    assert not_equal >= 1, "上传 taxonomy 后应至少有一条 MAG 拿到 Genus 标签"


def test_mag_label_falls_back_without_tax(mh_inputs):
    """没 taxonomy 时 label == MAG id。"""
    ab, _, _, _ = mh_inputs
    r = mag_heatmap.analyze(ab,
                            params={"filter_mode": "top_n", "top_n_count": 10})
    assert (r.stats["label"] == r.stats["MAG"]).all()


def test_phylum_legend_rendered(mh_inputs):
    """右侧图例区应至少包含一个 'Phylum' 标题的 legend。"""
    ab, tax, ks, md = mh_inputs
    r = mag_heatmap.analyze(ab, tax, ks, md,
                            params={"filter_mode": "top_n", "top_n_count": 30})
    axes = r.figure.axes
    # 第三个 Axes 是图例区，查找其 legend
    legends = []
    for a in axes:
        for child in a.get_children():
            if isinstance(child, plt.matplotlib.legend.Legend):
                legends.append(child)
    titles = [lg.get_title().get_text() for lg in legends]
    assert any("Phylum" in t for t in titles), \
        f"应有 Phylum 图例，实际 legend 标题：{titles}"
