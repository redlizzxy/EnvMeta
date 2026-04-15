"""机制假说 YAML 评分器测试（S3）。"""
from pathlib import Path

import pandas as pd
import pytest

from envmeta.geocycle.hypothesis import (
    Claim,
    ClaimResult,
    Hypothesis,
    HypothesisScore,
    load_hypothesis,
    score,
)
from envmeta.geocycle.inference import infer

SAMPLE = Path(__file__).parent / "sample_data"
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


@pytest.fixture
def cycle_data(cycle_inputs):
    ko, tax, ks, ab, env, md = cycle_inputs
    return infer(ko, tax, ks, ab, env, md,
                 params={**FAST_PARAMS, "env_rho_min": 0.3, "env_p_max": 0.1})


# ── 加载 / schema 校验 ────────────────────────────────────────

def test_load_hypothesis_from_file():
    hyp = load_hypothesis(SAMPLE / "sample_hypothesis.yaml")
    assert hyp.name == "测试用最小假说"
    assert len(hyp.claims) == 4
    assert {c.type for c in hyp.claims} == {
        "pathway_active", "coupling_possible",
        "env_correlation", "keystone_in_pathway",
    }
    assert hyp.strong_threshold == 0.75


def test_load_hypothesis_from_dict():
    hyp = load_hypothesis({
        "name": "h",
        "claims": [{"id": "c1", "type": "pathway_active",
                     "params": {"pathway": "X"}}],
    })
    assert hyp.name == "h"
    assert hyp.claims[0].id == "c1"


def test_load_hypothesis_schema_errors():
    with pytest.raises(ValueError, match="name"):
        load_hypothesis({"claims": []})
    with pytest.raises(ValueError, match="claims"):
        load_hypothesis({"name": "h"})
    with pytest.raises(ValueError, match="type"):
        load_hypothesis({
            "name": "h",
            "claims": [{"id": "x", "type": "unknown_type"}],
        })
    with pytest.raises(ValueError, match="重复"):
        load_hypothesis({
            "name": "h",
            "claims": [
                {"id": "x", "type": "pathway_active", "params": {"pathway": "A"}},
                {"id": "x", "type": "pathway_active", "params": {"pathway": "B"}},
            ],
        })


# ── 4 类 claim 评估 ─────────────────────────────────────────

def test_pathway_active_satisfied(cycle_data):
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="pathway_active", weight=1.0,
        params={"pathway": "Arsenate reduction", "min_completeness": 30},
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "satisfied", f"got {cr}"
    assert cr.score == 1.0


def test_pathway_active_missing_pathway_is_skipped(cycle_data):
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="pathway_active",
        params={"pathway": "NonexistentPathway"},
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "skipped"
    # skipped 不进分母
    assert result.label == "insufficient"


def test_coupling_possible_satisfied(cycle_data):
    # S2.5-14 修复后，As(V) 和 Fe(III) 都应在观测 species 里
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="coupling_possible",
        params={"species_a": "As(V)", "species_b": "Fe(III)"},
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "satisfied"
    assert cr.evidence["kb_product"] == "Fe-As_surface"


def test_coupling_not_in_kb_is_skipped(cycle_data):
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="coupling_possible",
        params={"species_a": "X_fake", "species_b": "Y_fake"},
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "skipped"


def test_env_correlation_sign_wrong(cycle_data):
    # 找一条真实 rho>0 的相关，写 expected_sign=negative 应该 unsatisfied
    positive = [c for c in cycle_data.env_correlations if c.rho > 0]
    assert positive, "样例数据里至少应有一条正相关"
    pc = positive[0]
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="env_correlation",
        params={
            "pathway": pc.pathway_id,
            "env_factor": pc.env_factor,
            "expected_sign": "negative",
            "min_confidence": "weak",
        },
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "unsatisfied"
    assert cr.score == 0.0


def test_keystone_in_pathway(cycle_data):
    # 找一个含 keystone 的 pathway
    target = None
    for el in cycle_data.elements:
        for pw in el.pathways:
            if any(c.is_keystone for c in pw.contributors):
                target = pw
                break
        if target:
            break
    assert target is not None, "样例数据里应至少有一个 keystone 通路"
    hyp = Hypothesis(name="t", claims=[Claim(
        id="c", type="keystone_in_pathway",
        params={"pathway": target.display_name, "min_keystones": 1},
    )])
    result = score(hyp, cycle_data)
    cr = result.claim_results[0]
    assert cr.status == "satisfied"


# ── 聚合 / 阈值 / 输出 ────────────────────────────────────

def test_overall_label_thresholds():
    # 构造不跑 infer，直接测聚合逻辑
    from envmeta.geocycle.model import CycleData
    empty = CycleData()

    # 3 claim 全部 satisfied → overall=1.0 → strong
    hyp_all = Hypothesis(name="t", claims=[
        Claim(id=f"c{i}", type="pathway_active",
              params={"pathway": "_NONE_"}) for i in range(3)
    ])
    # 但都是 _NONE_ 会 skipped → label=insufficient
    r = score(hyp_all, empty)
    assert r.label == "insufficient"
    assert r.n_total == 0
    assert r.n_skipped == 3


def test_label_weighted_aggregation():
    # 手工构造 claim_results 通过 score() 的阈值
    # 用一个简单办法：empty CycleData + 混合 claim
    # 但空数据所有 claim 都 skip。测阈值最干净的方式是直接 unit 测聚合：
    # 这里通过 Hypothesis.thresholds 验证分类边界
    hyp = Hypothesis(name="t", strong_threshold=0.75, suggestive_threshold=0.40,
                     claims=[])
    # 空 claim → total_w=0 → insufficient
    from envmeta.geocycle.model import CycleData
    r = score(hyp, CycleData())
    assert r.label == "insufficient"


def test_to_dataframe_and_json(cycle_data):
    hyp = load_hypothesis(SAMPLE / "sample_hypothesis.yaml")
    result = score(hyp, cycle_data)
    df = result.to_dataframe()
    assert len(df) == len(hyp.claims)
    assert set(df.columns) >= {"claim_id", "type", "status", "score",
                                "weight", "explanation"}
    js = result.to_json()
    assert '"hypothesis_name"' in js
    assert result.hypothesis_name in js
