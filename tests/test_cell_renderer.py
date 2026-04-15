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


# ── S2.5-6 化学式 mathtext + 重复中间产物折叠 ──────────────

def test_pretty_formula_overrides():
    from envmeta.geocycle.cell_renderer import _pretty_formula
    assert "SO_4^{2-}" in _pretty_formula("SO4-2")
    assert "N_2" in _pretty_formula("N2")
    assert "N_2O" in _pretty_formula("N2O")
    assert "As_2S_3" in _pretty_formula("As2S3")
    assert _pretty_formula("As(III)").startswith("$\\mathbf{")
    # 兜底：未知字符串保持原样
    assert _pretty_formula("cysteine") == "cysteine"
    # 空值
    assert _pretty_formula(None) == ""
    assert _pretty_formula("") == ""


def test_pretty_formula_regex_charge():
    from envmeta.geocycle.cell_renderer import _pretty_formula
    # NO3- 自动识别
    s = _pretty_formula("NO3-")
    assert "NO_3" in s and "-" in s


def test_collapse_identical_intermediates():
    """Fe transport 式级联（所有产物 = Fe_internal）触发 parallel_complex 合并为单 bundle。"""
    fig, ax = _make_ax()
    steps = [
        {"gene": f"G{i}", "substrate": "Fe(III)", "product": "Fe_internal"}
        for i in range(5)
    ]
    r = draw_cascade_cell(
        ax, cy=3.0, cell_x0=3.0, cell_w=10.0, cell_h=1.3,
        title="Fe transport", mag_label="Foo",
        steps=steps, element_color="#E67E22",
    )
    # S2.5-10 post: parallel_complex → 1 bundle 椭圆（而非 5 个小椭圆）
    assert len(r["gene_positions"]) == 1
    plt.close(fig)


def test_mixed_intermediates_not_collapsed():
    """不同中间产物的 denitrification 级联应保留中间产物显示。"""
    fig, ax = _make_ax()
    steps = [
        {"gene": "narG", "substrate": "NO3-", "product": "NO2-"},
        {"gene": "nirK", "substrate": "NO2-", "product": "NO"},
        {"gene": "norB", "substrate": "NO",   "product": "N2O"},
    ]
    r = draw_cascade_cell(
        ax, cy=3.0, cell_x0=3.0, cell_w=10.0, cell_h=1.3,
        title="denitrifier", mag_label="x", steps=steps,
    )
    assert len(r["gene_positions"]) == 3
    plt.close(fig)


# ── S2.5-7 排版 & parallel complex ──────────────────────

def test_pretty_formula_is_bold():
    """_pretty_formula 应输出 \\mathbf 而非 \\mathrm。"""
    from envmeta.geocycle.cell_renderer import _pretty_formula
    s = _pretty_formula("SO4-2")
    assert "\\mathbf" in s
    assert "\\mathrm" not in s


def test_long_gene_name_autosize():
    """TC.FEV.OM 这类长名应有更大椭圆 + 更小字号。"""
    from envmeta.geocycle.cell_renderer import _gene_size
    w_short, _, fs_short = _gene_size("narG")
    w_long, _, fs_long = _gene_size("TC.FEV.OM")
    assert w_long > w_short
    assert fs_long <= fs_short


def test_parallel_complex_no_internal_arrows():
    """6 个 KO 全部 substrate=NO3- product=NO2- → 合并为单 bundle 椭圆，无内部箭头。"""
    from matplotlib.patches import FancyArrowPatch
    fig, ax = _make_ax()
    steps = [
        {"gene": g, "substrate": "NO3-", "product": "NO2-"}
        for g in ("narG", "narH", "narI", "napA", "napB", "narB")
    ]
    r = draw_cascade_cell(
        ax, cy=3.0, cell_x0=3.0, cell_w=12.0, cell_h=1.3,
        title="Nitrate reduction", mag_label="test", steps=steps,
    )
    # S2.5-10 post: 6 个并联 gene 合并为 1 个 bundle ellipse
    assert len(r["gene_positions"]) == 1
    # 外部 substrate→cell 1 + cell→product 1 = 2
    n_arrows = sum(1 for p in ax.patches if isinstance(p, FancyArrowPatch))
    assert n_arrows <= 3, f"parallel complex 应无内部箭头，实际 {n_arrows}"
    plt.close(fig)


def test_denitrification_has_internal_arrows():
    """denitrification 3 基因串联应有内部箭头。"""
    from matplotlib.patches import FancyArrowPatch
    fig, ax = _make_ax()
    steps = [
        {"gene": "narG", "substrate": "NO3-", "product": "NO2-"},
        {"gene": "nirK", "substrate": "NO2-", "product": "NO"},
        {"gene": "norB", "substrate": "NO",   "product": "N2O"},
    ]
    draw_cascade_cell(
        ax, cy=3.0, cell_x0=3.0, cell_w=10.0, cell_h=1.3,
        title="denitrifier", mag_label="x", steps=steps,
    )
    # 2 内部箭头 + 2 外部箭头 = 4
    n_arrows = sum(1 for p in ax.patches if isinstance(p, FancyArrowPatch))
    assert n_arrows >= 4
    plt.close(fig)


# ── S2.5-10d 段内复合体分段 ───────────────────────────────

def test_segment_by_complex_basic():
    from envmeta.geocycle.cell_renderer import _segment_by_complex
    steps = [
        {"gene": "sqr", "complex": None},
        {"gene": "soxA", "complex": "M00595"},
        {"gene": "soxB", "complex": "M00595"},
        {"gene": "foo", "complex": None},
        {"gene": "bar", "complex": "M99999"},
    ]
    segs = _segment_by_complex(steps)
    assert segs == [[0], [1, 2], [3], [4]]


def test_segment_two_different_complexes_not_merged():
    from envmeta.geocycle.cell_renderer import _segment_by_complex
    steps = [
        {"gene": "a", "complex": "M001"},
        {"gene": "b", "complex": "M001"},
        {"gene": "c", "complex": "M002"},
        {"gene": "d", "complex": "M002"},
    ]
    assert _segment_by_complex(steps) == [[0, 1], [2, 3]]


def test_segment_none_always_single():
    from envmeta.geocycle.cell_renderer import _segment_by_complex
    steps = [{"complex": None}, {"complex": None}]
    # 两个连续 None 不合并
    assert _segment_by_complex(steps) == [[0], [1]]


def test_sulfide_oxidation_sqr_plus_sox_renders_without_crash():
    """真实场景：sqr (complex=None) + Sox 5 亚基 (complex=M00595)。
    S2.5-10 post: 2 segments → 1 sqr gene + 1 Sox bundle = 2 gene_positions
    """
    fig, ax = _make_ax()
    steps = [
        {"gene": "sqr", "substrate": "S-2", "product": "S0", "complex": None},
        {"gene": "soxA", "substrate": "S2O3-2", "product": "SO4-2", "complex": "M00595"},
        {"gene": "soxX", "substrate": "S2O3-2", "product": "SO4-2", "complex": "M00595"},
        {"gene": "soxB", "substrate": "S2O3-2", "product": "SO4-2", "complex": "M00595"},
        {"gene": "soxY", "substrate": "S2O3-2", "product": "SO4-2", "complex": "M00595"},
        {"gene": "soxZ", "substrate": "S2O3-2", "product": "SO4-2", "complex": "M00595"},
    ]
    r = draw_cascade_cell(
        ax, cy=5, cell_x0=1, cell_w=15, cell_h=2,
        title="Sulfide oxidation", mag_label="Thiobacillus",
        steps=steps,
    )
    # 2 段：单个 sqr + 5 Sox bundle
    assert len(r["gene_positions"]) == 2
    plt.close(fig)


def test_bundle_label_common_prefix():
    from envmeta.geocycle.cell_renderer import _bundle_label
    # narG/narH/narI → narG/H/I（共同前缀缩写）
    steps = [{"gene": n} for n in ("narG", "narH", "narI")]
    label = _bundle_label(steps)
    assert label == "narG/H/I"


def test_bundle_label_no_truncation_mixed():
    from envmeta.geocycle.cell_renderer import _bundle_label
    # narG/napA 等混合前缀不可缩 → 全拼
    steps = [{"gene": n} for n in ("narG", "narH", "narI", "napA", "napB", "narB")]
    label = _bundle_label(steps)
    # 不得截断 ("..." 或 "/N)")
    assert "..." not in label
    # 全部 6 个基因名都出现
    for g in ("narG", "narH", "narI", "napA", "napB", "narB"):
        assert g in label


def test_bundle_label_no_common_prefix():
    from envmeta.geocycle.cell_renderer import _bundle_label
    steps = [{"gene": n} for n in ("aoxA", "aoxB")]
    label = _bundle_label(steps)
    assert "aox" in label
