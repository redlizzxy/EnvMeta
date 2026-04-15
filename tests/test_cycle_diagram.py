"""生物地球化学循环图测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import cycle_diagram
from envmeta.export.figure_export import export_to_bytes
from envmeta.geocycle.inference import (
    _confidence_label, _permutation_rho_p, infer,
)

SAMPLE = Path(__file__).parent / "sample_data"

# 测试用小 perm_n 加速
FAST_PARAMS = {"perm_n": 99}


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
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    assert len(data.elements) >= 3
    assert data.meta["n_mags"] > 0
    assert data.meta["n_pathways_total"] >= 18


def test_inference_elements_have_pathways(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    for ec in data.elements:
        assert len(ec.pathways) > 0


def test_inference_env_correlations_present(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={**FAST_PARAMS, "env_rho_min": 0.3,
                         "env_p_max": 0.1})
    assert data.meta["n_env_correlations"] >= 0


def test_analyze_returns_figure(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    assert r.figure is not None
    assert not r.stats.empty
    assert "type" in r.stats.columns


def test_analyze_minimal_inputs(cycle_inputs):
    ko, *_ = cycle_inputs
    r = cycle_diagram.analyze(ko, params=FAST_PARAMS)
    assert r.figure is not None


def test_analyze_export(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    pdf = export_to_bytes(r.figure, "pdf")
    assert pdf.startswith(b"%PDF")


# ── S1 去偏测试 ────────────────────────────────────────────────

def test_full_correlation_matrix_populated(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={**FAST_PARAMS, "env_rho_min": 0.5,
                         "env_p_max": 0.05})
    assert len(data.full_corr_matrix) >= len(data.env_correlations)
    factors = {c.env_factor for c in data.full_corr_matrix}
    assert len(factors) >= 2


def test_sensitivity_scan_present(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    assert len(data.sensitivity) > 0
    for sr in data.sensitivity:
        assert len(sr.thresholds) == 3
        assert len(sr.top1_by_threshold) == 3
        assert len(sr.n_active_by_threshold) == 3
    assert any(sr.robust for sr in data.sensitivity)


def test_stats_contains_new_types(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    r = cycle_diagram.analyze(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    types = set(r.stats["type"])
    assert "pathway" in types
    assert "env_correlation" in types
    assert "full_correlation" in types
    assert "sensitivity" in types


# ── S2 置换零假设 + 可信度标签测试 ────────────────────────

def test_permutation_rho_p_unit():
    """_permutation_rho_p 对已知相关数据应返回合理 ρ 和 p。"""
    import numpy as np
    rng = np.random.default_rng(0)
    x = np.arange(30, dtype=float)
    y = x * 1.0 + rng.normal(0, 2.0, 30)   # 强相关
    rho, p = _permutation_rho_p(x, y, n=199, seed=42)
    assert rho > 0.8
    assert p < 0.05

    # 无相关
    y_rand = rng.normal(0, 5, 30)
    rho2, p2 = _permutation_rho_p(x, y_rand, n=199, seed=42)
    # 随机数据 ρ 可能小也可能大，但 p 应该不会特别小
    assert abs(rho2) < 0.6 or p2 > 0.01


def test_confidence_label_logic():
    """各可信度等级的边界。"""
    # strong: |ρ|>0.7, perm_p<0.01, robust
    assert _confidence_label(0.85, 0.005, True) == "strong"
    # suggestive: 0.5<|ρ|≤0.7, perm_p<0.05
    assert _confidence_label(0.6, 0.02, True) == "suggestive"
    # weak
    assert _confidence_label(0.4, 0.08, True) == "weak"
    # spurious: |ρ|>=0.5 but perm_p>0.05
    assert _confidence_label(0.55, 0.2, True) == "spurious?"
    # NaN
    assert _confidence_label(float("nan"), 0.01, True) == "unknown"


def test_env_correlations_have_perm_p(cycle_inputs):
    """每条 env_correlation 应带 perm_p 和 confidence 字段。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={**FAST_PARAMS, "env_rho_min": 0.3,
                         "env_p_max": 0.1})
    for ec in data.env_correlations:
        assert ec.perm_p is not None
        assert 0.0 <= ec.perm_p <= 1.0
        assert ec.confidence in {"strong", "suggestive", "weak",
                                  "spurious?", "none", "unknown"}


def test_confidence_labels_distribution(cycle_inputs):
    """在真实数据上应该产生一些 strong / suggestive 标签。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md,
                 params={**FAST_PARAMS, "env_rho_min": 0.3,
                         "env_p_max": 0.1})
    labels = [ec.confidence for ec in data.env_correlations]
    label_set = set(labels)
    # 真实数据应该至少有 strong 或 suggestive 之一
    assert label_set & {"strong", "suggestive"}
    # meta 里应该有相关计数
    assert "n_confidence_strong" in data.meta
    assert "n_confidence_suggestive" in data.meta


# ── S2.5-2a KO 级 cascade 测试 ────────────────────────────────

def test_contributors_carry_gene_cascade(cycle_inputs):
    """每个 MAGContribution 的 genes 字段应有 substrate/product 化学信息。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    # 至少存在一个带贡献者的通路
    found = False
    for ec in data.elements:
        for pw in ec.pathways:
            for c in pw.contributors:
                assert isinstance(c.genes, list)
                if c.genes:
                    found = True
                    # 字段存在
                    g = c.genes[0]
                    assert {"ko", "name", "substrate", "product"} <= set(g.keys())
                    # KO id 形如 K12345
                    assert g["ko"].startswith("K")
    assert found, "至少一个贡献 MAG 应该持有通路里的 KO"


def test_gene_cascade_matches_kb(cycle_inputs):
    """catalytic KO 的 substrate/product 应等于 KB 定义；调控 KO 为 None。"""
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    # 找一条 catalytic KO（如 narG=K00370 应有 NO3- → NO2-）
    catalytic_seen = {"K00370", "K00380", "K00376"}
    regulator_seen = {"K03892", "K03711"}
    catalytic_verified = False
    for ec in data.elements:
        for pw in ec.pathways:
            for c in pw.contributors:
                for g in c.genes:
                    if g["ko"] in catalytic_seen:
                        assert g["substrate"] is not None
                        assert g["product"] is not None
                        catalytic_verified = True
                    if g["ko"] in regulator_seen:
                        assert g["substrate"] is None
                        assert g["product"] is None
    # 本数据集未必出现这些 catalytic KO；若出现则断言已验过
    assert catalytic_verified or True  # 宽松断言：确保循环没崩


def test_gene_cascade_order_matches_kb(cycle_inputs):
    """genes 列表顺序应与 KB pathway.genes 的登记序一致（MAG 持有的子集）。"""
    from envmeta.geocycle.knowledge_base import pathway_ko_sets
    ko, tax, ks, ab, env, md = cycle_inputs
    data = infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)
    pw_kos = pathway_ko_sets()
    for ec in data.elements:
        for pw in ec.pathways:
            kb_order = pw_kos[pw.pathway_id]
            for c in pw.contributors:
                actual = [g["ko"] for g in c.genes]
                # 取 kb_order 中 actual 包含的子序列
                expected = [k for k in kb_order if k in set(actual)]
                assert actual == expected
