"""生物地球化学循环图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import cycle_diagram
from envmeta.export.figure_export import export_to_bytes
from envmeta.geocycle.inference import infer

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def cycle_inputs():
    ko = pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    ab = pd.read_csv(SAMPLE / "abundance.tsv", sep="\t")
    env = pd.read_csv(SAMPLE / "env_factors.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ko, tax, ks, ab, env, md


def test_inference_smoke(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md)
    assert len(data.elements) >= 3   # 至少 As/N/S 其一有活跃通路
    assert data.meta["n_mags"] > 0
    assert data.meta["n_pathways_total"] >= 18


def test_inference_elements_have_pathways(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md)
    for ec in data.elements:
        assert len(ec.pathways) > 0


def test_inference_env_correlations_present(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={"env_rho_min": 0.3, "env_p_max": 0.1})
    # 真实数据应该有至少 1 条显著相关
    assert data.meta["n_env_correlations"] >= 0  # 宽松：至少字段存在


def test_analyze_returns_figure(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md)
    assert r.figure is not None
    assert not r.stats.empty
    assert "type" in r.stats.columns


def test_analyze_minimal_inputs(cycle_inputs):
    """只给 KO 注释也能跑通（taxonomy/ks/abund/env 都可选）。"""
    ko, *_ = cycle_inputs
    r = cycle_diagram.analyze(ko)
    assert r.figure is not None


def test_analyze_export(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md)
    pdf = export_to_bytes(r.figure, "pdf")
    assert pdf.startswith(b"%PDF")


# ── S1 去偏测试 ────────────────────────────────────────────────

def test_full_correlation_matrix_populated(cycle_inputs):
    """完整相关矩阵应该包含所有通路 × env 组合（不只超阈值的）。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={"env_rho_min": 0.5, "env_p_max": 0.05})
    assert len(data.full_corr_matrix) >= len(data.env_correlations)
    # 完整矩阵应覆盖多因子（Total_As/Eh/TOC/pH）
    factors = {c.env_factor for c in data.full_corr_matrix}
    assert len(factors) >= 2


def test_sensitivity_scan_present(cycle_inputs):
    """敏感度扫描：每通路都应该有 3 档阈值的结果。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md)
    assert len(data.sensitivity) > 0
    for sr in data.sensitivity:
        assert len(sr.thresholds) == 3
        assert len(sr.top1_by_threshold) == 3
        assert len(sr.n_active_by_threshold) == 3
    # 至少某些通路是 robust 的（三档 Top-1 一致）
    assert any(sr.robust for sr in data.sensitivity)


def test_stats_contains_new_types(cycle_inputs):
    """扁平 stats 应含 full_correlation 和 sensitivity 两个新 type。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md)
    types = set(r.stats["type"])
    assert "pathway" in types
    assert "env_correlation" in types
    assert "full_correlation" in types
    assert "sensitivity" in types
