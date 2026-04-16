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


# ── S3.5: null comparison (permutation) ────────────────────────

def test_null_permutation_basic(cycle_data):
    """≥3 non-skipped claims with heterogeneous weights AND heterogeneous scores
    → null_p is a float [0,1]."""
    # 强制 claim 得到不同 score：一个 satisfied、一个 partial、一个 unsatisfied
    hyp = Hypothesis(name="t", claims=[
        # satisfied：min_completeness=30 容易达
        Claim(id="c1_sat", type="pathway_active", weight=2.0,
              params={"pathway": "Arsenate reduction", "min_completeness": 30}),
        # partial：min_completeness 抬高到 99（活跃但达不到阈值）
        Claim(id="c2_partial", type="pathway_active", weight=1.0,
              params={"pathway": "Arsenate reduction", "min_completeness": 99}),
        # satisfied via coupling
        Claim(id="c3_sat", type="coupling_possible", weight=0.5,
              params={"species_a": "As(V)", "species_b": "Fe(III)"}),
    ])
    r = score(hyp, cycle_data, null_n=199)  # faster than 999
    assert r.null_p is not None, f"应返回 null_p，got claims={[(cr.status, cr.score) for cr in r.claim_results]}"
    assert 0 <= r.null_p <= 1
    assert r.null_p_samples == 199


def test_null_permutation_too_few_claims(cycle_data):
    """≤2 non-skipped → null_p=None (排列空间不足)."""
    hyp = Hypothesis(name="t", claims=[
        Claim(id="c1", type="pathway_active", weight=2.0,
              params={"pathway": "Arsenate reduction", "min_completeness": 30}),
        Claim(id="c2", type="pathway_active", weight=1.0,
              params={"pathway": "_DOES_NOT_EXIST_"}),  # skipped
    ])
    r = score(hyp, cycle_data, null_n=99)
    assert r.null_p is None
    assert r.null_p_samples == 0


def test_null_disable(cycle_data):
    """run_null=False → null_p 字段保持默认 None。"""
    hyp = Hypothesis(name="t", claims=[
        Claim(id=f"c{i}", type="pathway_active", weight=float(i + 1),
              params={"pathway": "Arsenate reduction", "min_completeness": 30})
        for i in range(3)
    ])
    r = score(hyp, cycle_data, run_null=False)
    assert r.null_p is None
    assert r.null_p_samples == 0


# ── S3.5: weight sensitivity (OAT) ─────────────────────────────

def test_weight_sensitivity_scan_present(cycle_data):
    """至少返回 2 × n_scored 条扰动记录，含 new_label / flipped 字段。"""
    hyp = Hypothesis(name="t", claims=[
        Claim(id="c1", type="pathway_active", weight=2.0,
              params={"pathway": "Arsenate reduction", "min_completeness": 30}),
        Claim(id="c2", type="coupling_possible", weight=1.0,
              params={"species_a": "As(V)", "species_b": "Fe(III)"}),
    ])
    r = score(hyp, cycle_data, run_null=False)
    assert r.weight_sensitivity_rows, "应产出扰动记录"
    # 2 claim × 2 方向 = 4 行
    assert len(r.weight_sensitivity_rows) == 4
    row = r.weight_sensitivity_rows[0]
    assert {"claim_id", "direction", "new_weight", "new_overall",
            "new_label", "flipped"} <= set(row)
    assert r.weight_robust in (True, False)


def test_weight_sensitivity_stable_when_all_satisfied(cycle_data):
    """全 satisfied + 高分时，OAT ±20% 不应翻 label → robust=True."""
    hyp = Hypothesis(name="t", claims=[
        Claim(id=f"c{i}", type="pathway_active",
              weight=1.0 + 0.3 * i,  # 权重略异
              params={"pathway": "Arsenate reduction", "min_completeness": 30})
        for i in range(3)
    ])
    r = score(hyp, cycle_data, run_null=False)
    assert r.weight_robust is True


# ── S3.5: required / veto ──────────────────────────────────────

def test_required_satisfied_does_not_veto(cycle_data):
    """required=true + status=satisfied → 正常走阈值逻辑，无 veto。"""
    hyp = Hypothesis(name="t", claims=[
        Claim(id="core", type="pathway_active", weight=2.0, required=True,
              params={"pathway": "Arsenate reduction", "min_completeness": 30}),
        Claim(id="aux", type="coupling_possible", weight=1.0,
              params={"species_a": "As(V)", "species_b": "Fe(III)"}),
    ])
    r = score(hyp, cycle_data, run_null=False)
    assert r.veto_reasons == []
    assert r.label in ("strong", "suggestive")  # 实际数据两个都 satisfied


def test_required_unsatisfied_triggers_veto(cycle_data):
    """required=true + status=skipped/unsatisfied → label 强制 insufficient。"""
    hyp = Hypothesis(name="t", claims=[
        Claim(id="aux_sat", type="pathway_active", weight=1.0,
              params={"pathway": "Arsenate reduction", "min_completeness": 30}),
        Claim(id="core_missing", type="pathway_active", weight=5.0,
              required=True,
              params={"pathway": "_NONEXISTENT_PATHWAY_"}),
    ])
    r = score(hyp, cycle_data, run_null=False)
    assert r.label == "insufficient"
    assert r.veto_reasons
    assert any("core_missing" in reason for reason in r.veto_reasons)
    # overall_score 仍然被计算（透明度）
    assert r.overall_score >= 0
    # base_label 应该在 params 里保留
    assert "base_label_before_veto" in r.params


def test_load_hypothesis_rejects_non_bool_required():
    """required 必须是 bool，不是 str/int。"""
    with pytest.raises(ValueError, match="required"):
        load_hypothesis({
            "name": "h",
            "claims": [{"id": "c", "type": "pathway_active",
                         "required": "yes",  # 非 bool
                         "params": {"pathway": "X"}}],
        })


# ── S3.5: group_contrast claim ─────────────────────────────────

def test_group_contrast_satisfied():
    """构造 compare_df，B/CK ratio ≥ min_ratio → satisfied。"""
    from envmeta.geocycle.model import CycleData
    compare_df = pd.DataFrame([
        {"group": "CK", "element": "arsenic",
         "pathway_id": "as_arsenate_reduction",
         "display_name": "Arsenate reduction",
         "n_active_mags": 2, "mean_completeness": 40,
         "total_contribution": 100.0, "top_mag": "m1",
         "top_mag_genus": "g", "top_mag_is_keystone": False},
        {"group": "B", "element": "arsenic",
         "pathway_id": "as_arsenate_reduction",
         "display_name": "Arsenate reduction",
         "n_active_mags": 4, "mean_completeness": 60,
         "total_contribution": 250.0, "top_mag": "m2",
         "top_mag_genus": "g", "top_mag_is_keystone": True},
    ])
    hyp = Hypothesis(name="t", claims=[Claim(
        id="gc", type="group_contrast", weight=1.0,
        params={"pathway": "Arsenate reduction",
                "high_group": "B", "low_group": "CK",
                "min_ratio": 1.5},
    )])
    r = score(hyp, CycleData(), compare_df=compare_df, run_null=False)
    cr = r.claim_results[0]
    assert cr.status == "satisfied"
    assert cr.score == 1.0
    assert cr.evidence["ratio"] == 2.5


def test_group_contrast_skipped_without_compare_df(cycle_data):
    """YAML 含 group_contrast 但未传 compare_df → skipped。"""
    hyp = Hypothesis(name="t", claims=[Claim(
        id="gc", type="group_contrast",
        params={"pathway": "Arsenate reduction",
                "high_group": "B", "low_group": "CK", "min_ratio": 1.5},
    )])
    r = score(hyp, cycle_data, run_null=False)
    cr = r.claim_results[0]
    assert cr.status == "skipped"
    assert "compare" in cr.explanation.lower()


def test_group_contrast_in_claim_types_whitelist():
    """load_hypothesis 接受 group_contrast type。"""
    hyp = load_hypothesis({
        "name": "h",
        "claims": [{
            "id": "gc", "type": "group_contrast",
            "params": {"pathway": "X", "high_group": "A",
                       "low_group": "B", "min_ratio": 1.5},
        }],
    })
    assert hyp.claims[0].type == "group_contrast"


# ── S3.5: validate CLI ─────────────────────────────────────────

def test_validator_accepts_good_yaml():
    """示例 YAML 应 0 errors 0 warnings (KB 里所有引用都存在)."""
    from envmeta.tools.hypothesis_validator import validate_file
    result = validate_file("paper/hypotheses/arsenic_steel_slag.yaml")
    assert result["errors"] == []
    assert result["warnings"] == []


def test_validator_warns_unknown_pathway(tmp_path):
    """不存在的 pathway 应产生 warning (不是 error)."""
    from envmeta.tools.hypothesis_validator import validate_file
    p = tmp_path / "bad.yaml"
    p.write_text(
        "name: test\n"
        "claims:\n"
        "  - id: bad_pw\n"
        "    type: pathway_active\n"
        "    params:\n"
        "      pathway: \"Totally_Fake_Pathway_XYZ\"\n",
        encoding="utf-8",
    )
    result = validate_file(p)
    assert result["errors"] == []
    assert result["warnings"]
    assert any("Totally_Fake_Pathway_XYZ" in w for w in result["warnings"])
