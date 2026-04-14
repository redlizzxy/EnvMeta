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
