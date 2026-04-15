"""级联细胞渲染 primitive 单元测试（S2.5-2b）。"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from envmeta.geocycle.cell_renderer import (
    draw_cascade_cell, genes_to_steps, mini_heatmap,
)


def _make_ax():
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 6)
    ax.axis("off")
    return fig, ax


def test_draw_cascade_cell_basic_single_step():
    fig, ax = _make_ax()
    r = draw_cascade_cell(
        ax,
        cy=3.0, cell_x0=5.0, cell_w=4.0, cell_h=1.3,
        title="Test single", mag_label="Genus_A",
        steps=[{"gene": "narG", "substrate": "NO3-",
                "product": "NO2-", "color": "#3498DB",
                "corr": [0.2, 0.5, -0.1, 0.6]}],
        element_color="#3498DB",
    )
    assert r["substrate_pos"] is not None
    assert r["product_pos"] is not None
    assert r["cell_box"] == (5.0, 9.0, 3.0)
    assert len(r["gene_positions"]) == 1
    plt.close(fig)


def test_draw_cascade_cell_multi_step_chain():
    fig, ax = _make_ax()
    r = draw_cascade_cell(
        ax,
        cy=3.0, cell_x0=3.0, cell_w=10.0, cell_h=1.3,
        title="Denitrifier", mag_label="Gallionella",
        steps=[
            {"gene": "narG", "substrate": "NO3-", "product": "NO2-"},
            {"gene": "nirK", "substrate": "NO2-", "product": "NO"},
            {"gene": "norB", "substrate": "NO",   "product": "N2O"},
            {"gene": "nosZ", "substrate": "N2O",  "product": "N2"},
        ],
        element_color="#3498DB",
        show_heatmap=False,
    )
    assert len(r["gene_positions"]) == 4
    # 起始底物 / 最终产物应有
    assert r["substrate_pos"] is not None
    assert r["product_pos"] is not None
    plt.close(fig)


def test_empty_steps_fallback():
    fig, ax = _make_ax()
    r = draw_cascade_cell(
        ax,
        cy=3.0, cell_x0=5.0, cell_w=4.0, cell_h=1.3,
        title="empty", mag_label="", steps=[],
        element_color="#888",
    )
    assert r["substrate_pos"] is None
    assert r["product_pos"] is None
    assert r["gene_positions"] == []
    plt.close(fig)


def test_regulator_step_has_no_outside_chems():
    """substrate/product=None 时（如 arsR 调控），不画外部化学物节点。"""
    fig, ax = _make_ax()
    r = draw_cascade_cell(
        ax,
        cy=3.0, cell_x0=5.0, cell_w=4.0, cell_h=1.3,
        title="regulator", mag_label="Foo",
        steps=[{"gene": "arsR", "substrate": None, "product": None}],
        element_color="#E74C3C",
    )
    # 外部节点应为 None
    assert r["substrate_pos"] is None
    assert r["product_pos"] is None
    # 基因仍然画出
    assert len(r["gene_positions"]) == 1
    plt.close(fig)


def test_show_outside_chems_false_skips_external_nodes():
    fig, ax = _make_ax()
    r = draw_cascade_cell(
        ax,
        cy=3.0, cell_x0=5.0, cell_w=4.0, cell_h=1.3,
        title="t", mag_label="",
        steps=[{"gene": "x", "substrate": "A", "product": "B"}],
        show_outside_chems=False,
    )
    assert r["substrate_pos"] is None
    assert r["product_pos"] is None
    plt.close(fig)


def test_genes_to_steps_conversion():
    gene_list = [
        {"ko": "K00370", "name": "narG", "substrate": "NO3-", "product": "NO2-"},
        {"ko": "K00368", "name": "nirK", "substrate": "NO2-", "product": "NO"},
    ]
    env_corr = {"K00370": [0.1, 0.2, 0.3, 0.4]}
    steps = genes_to_steps(gene_list, default_color="#3498DB", env_corr=env_corr)
    assert len(steps) == 2
    assert steps[0]["gene"] == "narG"
    assert steps[0]["substrate"] == "NO3-"
    assert steps[0]["corr"] == [0.1, 0.2, 0.3, 0.4]
    assert steps[1]["corr"] is None   # env_corr 没给 K00368


def test_mini_heatmap_handles_nan():
    """NaN 值不应崩；应被画成灰色。"""
    fig, ax = _make_ax()
    mini_heatmap(ax, 1.0, 1.0, [0.5, float("nan"), -0.3, None])
    plt.close(fig)
