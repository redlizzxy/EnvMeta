"""机制假说 YAML 评分器（S3 / L2 层）。

用户写一份假说 YAML，列出若干 claim（期望在数据里观察到的现象），
本模块对着 CycleData 打分，给 strong / suggestive / weak / insufficient 标签。

**设计原则**：
- 假说评估与推断**架构性解耦**：inference.py 不依赖本模块；本模块只读 CycleData
- 描述性：评分反映"数据支持度"，不下"因果结论"
- skipped 不扣分：YAML claim 指向 KB 没有 / data 里没跑到的对象 → 分母剔除，
  不让用户因"写多了无关 claim"被压低分数

支持 4 类 claim（v1）：
1. `pathway_active`        —— 某通路活跃
2. `coupling_possible`     —— 某两物种 KB 有耦合 AND 数据里两端 species 都观测到
3. `env_correlation`       —— (通路, 环境因子) 相关性符号 + confidence 达阈
4. `keystone_in_pathway`   —— 某通路包含 ≥N 个 keystone contributor
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

from envmeta.geocycle.knowledge_base import couplings as kb_couplings
from envmeta.geocycle.model import CycleData, EnvCorrelation, PathwayActivity


# =============================================================================
# 数据模型
# =============================================================================

CLAIM_TYPES = (
    "pathway_active",
    "coupling_possible",
    "env_correlation",
    "keystone_in_pathway",
)

CONFIDENCE_RANK = {
    "spurious?": -1,
    "none": 0,
    "unknown": 0,
    "weak": 1,
    "suggestive": 2,
    "strong": 3,
}


@dataclass
class Claim:
    id: str
    type: str
    weight: float = 1.0
    description: str = ""
    params: dict = field(default_factory=dict)


@dataclass
class Hypothesis:
    name: str
    description: str = ""
    author: str = ""
    version: str = "1.0"
    claims: list[Claim] = field(default_factory=list)
    strong_threshold: float = 0.75
    suggestive_threshold: float = 0.40


@dataclass
class ClaimResult:
    claim_id: str
    claim_type: str
    status: str           # satisfied | partial | unsatisfied | skipped
    score: float          # [0, 1]
    weight: float
    evidence: dict = field(default_factory=dict)
    explanation: str = ""


@dataclass
class HypothesisScore:
    hypothesis_name: str
    overall_score: float
    label: str            # strong | suggestive | weak | insufficient
    n_satisfied: int
    n_total: int
    n_skipped: int
    claim_results: list[ClaimResult] = field(default_factory=list)
    params: dict = field(default_factory=dict)

    def to_dataframe(self) -> pd.DataFrame:
        rows = []
        for cr in self.claim_results:
            rows.append({
                "claim_id": cr.claim_id,
                "type": cr.claim_type,
                "status": cr.status,
                "score": round(cr.score, 3),
                "weight": cr.weight,
                "explanation": cr.explanation,
            })
        return pd.DataFrame(rows)

    def to_json(self) -> str:
        return json.dumps(asdict(self), ensure_ascii=False, indent=2)


# =============================================================================
# 加载 / schema 校验
# =============================================================================

def load_hypothesis(src: str | Path | dict) -> Hypothesis:
    """从 YAML 文件路径 / YAML 字符串 / 已解析 dict 加载假说。

    做最小 schema 校验：
    - 必须有 `name` (非空字符串)
    - 必须有 `claims` 列表（≥1 条）
    - 每个 claim 必须有 `id`（非空） + `type`（在 CLAIM_TYPES 中）
    出错时 raise `ValueError`，消息包含具体字段。
    """
    import yaml

    if isinstance(src, dict):
        raw = src
    elif isinstance(src, (str, Path)):
        p = Path(src) if not isinstance(src, Path) else src
        if p.exists() and p.is_file():
            with open(p, "r", encoding="utf-8") as f:
                raw = yaml.safe_load(f)
        else:
            # 当作 YAML 字符串
            raw = yaml.safe_load(str(src))
    else:
        raise TypeError(f"load_hypothesis: 不支持的输入类型 {type(src)}")

    if not isinstance(raw, dict):
        raise ValueError("YAML 根必须是映射（dict），请检查格式")

    name = raw.get("name")
    if not name or not isinstance(name, str):
        raise ValueError("YAML 缺失必填字段 `name`")

    claims_raw = raw.get("claims")
    if not claims_raw or not isinstance(claims_raw, list):
        raise ValueError("YAML 必须有 `claims` 列表，且至少 1 条")

    claims: list[Claim] = []
    seen_ids: set[str] = set()
    for i, c in enumerate(claims_raw):
        if not isinstance(c, dict):
            raise ValueError(f"claim[{i}] 必须是 dict，实际 {type(c).__name__}")
        cid = c.get("id")
        if not cid or not isinstance(cid, str):
            raise ValueError(f"claim[{i}] 缺失 `id` 或 id 非字符串")
        if cid in seen_ids:
            raise ValueError(f"claim id 重复: {cid!r}")
        seen_ids.add(cid)
        ctype = c.get("type")
        if ctype not in CLAIM_TYPES:
            raise ValueError(
                f"claim[{cid}] 的 type={ctype!r} 非法；合法值: {CLAIM_TYPES}"
            )
        claims.append(Claim(
            id=cid,
            type=ctype,
            weight=float(c.get("weight", 1.0)),
            description=str(c.get("description", "")),
            params=dict(c.get("params") or {}),
        ))

    return Hypothesis(
        name=str(name),
        description=str(raw.get("description", "")),
        author=str(raw.get("author", "")),
        version=str(raw.get("version", "1.0")),
        claims=claims,
        strong_threshold=float(raw.get("strong_threshold", 0.75)),
        suggestive_threshold=float(raw.get("suggestive_threshold", 0.40)),
    )


# =============================================================================
# 辅助 —— pathway 名模糊匹配
# =============================================================================

def _norm(s: str | None) -> str:
    if s is None:
        return ""
    return str(s).strip().lower().replace("_", " ").replace("-", " ")


def _find_pathway(data: CycleData, name: str) -> tuple[PathwayActivity, str] | None:
    """在 CycleData 里按 display_name / pathway_id 模糊匹配 pathway。

    返回 (pathway, element_id) 或 None。
    """
    target = _norm(name)
    if not target:
        return None
    for el in data.elements:
        for pw in el.pathways:
            if _norm(pw.display_name) == target or _norm(pw.pathway_id) == target:
                return pw, el.element_id
    # 退化：子串匹配（避免用户写 "arsenate reduction" vs "Arsenate reduction (ars)"）
    for el in data.elements:
        for pw in el.pathways:
            if target in _norm(pw.display_name) or target in _norm(pw.pathway_id):
                return pw, el.element_id
    return None


def _collect_observed_species(data: CycleData) -> set[str]:
    """扫描所有 contributors 的 genes，收集被观测到的 substrate/product species 集合。"""
    out: set[str] = set()
    for el in data.elements:
        for pw in el.pathways:
            for c in pw.contributors:
                for g in c.genes or []:
                    s = g.get("substrate")
                    p = g.get("product")
                    if s:
                        out.add(str(s))
                    if p:
                        out.add(str(p))
    return out


# =============================================================================
# Claim 评估器
# =============================================================================

def _eval_pathway_active(claim: Claim, data: CycleData) -> ClaimResult:
    p = claim.params
    pw_name = p.get("pathway")
    if not pw_name:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"reason": "missing params.pathway"},
            explanation="未指定 pathway 名",
        )
    match = _find_pathway(data, pw_name)
    if match is None:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"pathway_query": pw_name},
            explanation=f"数据里找不到通路 {pw_name!r}（可能被 element_filter 过滤）",
        )
    pw, element = match
    min_comp = float(p.get("min_completeness", 50))
    min_contrib = float(p.get("min_contribution", 0))
    ev = {
        "pathway_id": pw.pathway_id,
        "element": element,
        "n_active_mags": pw.n_active_mags,
        "mean_completeness": round(pw.mean_completeness, 2),
        "total_contribution": round(pw.total_contribution, 2),
    }
    if pw.n_active_mags <= 0:
        return ClaimResult(
            claim.id, claim.type, "unsatisfied", 0.0, claim.weight,
            evidence=ev,
            explanation=f"{pw.display_name}: 无活跃 MAG",
        )
    ok_comp = pw.mean_completeness >= min_comp
    ok_contrib = pw.total_contribution >= min_contrib
    if ok_comp and ok_contrib:
        return ClaimResult(
            claim.id, claim.type, "satisfied", 1.0, claim.weight,
            evidence=ev,
            explanation=(
                f"{pw.display_name}: {pw.n_active_mags} 个 MAG 活跃，"
                f"平均完整度 {pw.mean_completeness:.0f}%，"
                f"总贡献 {pw.total_contribution:.1f}"
            ),
        )
    return ClaimResult(
        claim.id, claim.type, "partial", 0.5, claim.weight,
        evidence=ev,
        explanation=(
            f"{pw.display_name}: 活跃但未达阈值 "
            f"(completeness {pw.mean_completeness:.0f}%/{min_comp:.0f}%, "
            f"contribution {pw.total_contribution:.1f}/{min_contrib:.1f})"
        ),
    )


def _eval_coupling_possible(claim: Claim, data: CycleData) -> ClaimResult:
    p = claim.params
    a = p.get("species_a")
    b = p.get("species_b")
    if not a or not b:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"reason": "missing species_a/species_b"},
            explanation="未指定 species_a / species_b",
        )
    # KB 里是否有这对配对（允许 (a,b) 或 (b,a)）
    kb = kb_couplings()
    matched = None
    for cp in kb:
        sa, sb = cp.get("species_a"), cp.get("species_b")
        if (sa == a and sb == b) or (sa == b and sb == a):
            matched = cp
            break
    if matched is None:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"species_a": a, "species_b": b},
            explanation=f"KB 里没有 {a}↔{b} 的耦合定义（skip，不扣分）",
        )
    observed = _collect_observed_species(data)
    obs_a = a in observed
    obs_b = b in observed
    ev = {
        "kb_product": matched.get("product"),
        "kb_type": matched.get("type"),
        "obs_a": obs_a,
        "obs_b": obs_b,
    }
    if obs_a and obs_b:
        return ClaimResult(
            claim.id, claim.type, "satisfied", 1.0, claim.weight,
            evidence=ev,
            explanation=(
                f"{a} + {b} → {matched.get('product')} "
                f"({matched.get('type')})：两端物种都被观测到"
            ),
        )
    if obs_a or obs_b:
        missing = b if not obs_b else a
        return ClaimResult(
            claim.id, claim.type, "partial", 0.5, claim.weight,
            evidence=ev,
            explanation=f"{a}↔{b} 耦合：数据里缺 {missing}（仅一端观测到）",
        )
    return ClaimResult(
        claim.id, claim.type, "unsatisfied", 0.0, claim.weight,
        evidence=ev,
        explanation=f"{a}↔{b} 耦合：数据里两端物种都未观测到",
    )


def _eval_env_correlation(claim: Claim, data: CycleData) -> ClaimResult:
    p = claim.params
    pw_name = p.get("pathway")
    env_factor = p.get("env_factor")
    expected = (p.get("expected_sign") or "any").lower()
    min_conf = (p.get("min_confidence") or "suggestive").lower()
    if not pw_name or not env_factor:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"reason": "missing pathway / env_factor"},
            explanation="未指定 pathway 或 env_factor",
        )
    match = _find_pathway(data, pw_name)
    target_pid = match[0].pathway_id if match else None
    target_disp = match[0].display_name if match else None

    corr: EnvCorrelation | None = None
    # 先在 filtered 里找；找不到再到 full_corr_matrix 里兜底
    for pool in (data.env_correlations, data.full_corr_matrix):
        for ec in pool:
            if ec.env_factor != env_factor:
                continue
            if (target_pid and ec.pathway_id == target_pid) or (
                target_disp and _norm(ec.pathway_id) == _norm(target_disp)
            ) or _norm(ec.pathway_id) == _norm(pw_name):
                corr = ec
                break
        if corr is not None:
            break

    if corr is None:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"pathway_query": pw_name, "env_factor": env_factor},
            explanation=f"找不到 ({pw_name}, {env_factor}) 的相关性记录",
        )

    ev = {
        "pathway_id": corr.pathway_id,
        "env_factor": corr.env_factor,
        "rho": round(corr.rho, 3),
        "perm_p": (None if corr.perm_p is None else round(corr.perm_p, 4)),
        "confidence": corr.confidence,
    }
    sign_ok = (
        expected == "any"
        or (expected == "positive" and corr.rho > 0)
        or (expected == "negative" and corr.rho < 0)
    )
    conf_rank = CONFIDENCE_RANK.get(corr.confidence, 0)
    min_rank = CONFIDENCE_RANK.get(min_conf, 2)
    if not sign_ok:
        return ClaimResult(
            claim.id, claim.type, "unsatisfied", 0.0, claim.weight,
            evidence=ev,
            explanation=(
                f"({pw_name}, {env_factor}): ρ={corr.rho:.2f} 与期望方向 "
                f"{expected!r} 不符"
            ),
        )
    if conf_rank >= min_rank and corr.confidence != "spurious?":
        return ClaimResult(
            claim.id, claim.type, "satisfied", 1.0, claim.weight,
            evidence=ev,
            explanation=(
                f"({pw_name}, {env_factor}): ρ={corr.rho:+.2f} "
                f"[{corr.confidence}]，方向 {expected} 成立"
            ),
        )
    # 符号对但 confidence 低一档 → partial
    if conf_rank == min_rank - 1 and conf_rank >= 0:
        return ClaimResult(
            claim.id, claim.type, "partial", 0.5, claim.weight,
            evidence=ev,
            explanation=(
                f"({pw_name}, {env_factor}): ρ={corr.rho:+.2f} "
                f"[{corr.confidence}]，置信度低于要求 [{min_conf}]"
            ),
        )
    return ClaimResult(
        claim.id, claim.type, "unsatisfied", 0.0, claim.weight,
        evidence=ev,
        explanation=(
            f"({pw_name}, {env_factor}): ρ={corr.rho:+.2f} "
            f"但 confidence={corr.confidence}，未达 {min_conf} 要求"
        ),
    )


def _eval_keystone_in_pathway(claim: Claim, data: CycleData) -> ClaimResult:
    p = claim.params
    pw_name = p.get("pathway")
    min_n = int(p.get("min_keystones", 1))
    if not pw_name:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"reason": "missing params.pathway"},
            explanation="未指定 pathway",
        )
    match = _find_pathway(data, pw_name)
    if match is None:
        return ClaimResult(
            claim.id, claim.type, "skipped", 0.0, claim.weight,
            evidence={"pathway_query": pw_name},
            explanation=f"数据里找不到通路 {pw_name!r}",
        )
    pw, _ = match
    keystones = [c for c in pw.contributors if c.is_keystone]
    ev = {
        "pathway_id": pw.pathway_id,
        "n_keystones": len(keystones),
        "keystone_labels": [c.label for c in keystones[:5]],
        "min_required": min_n,
    }
    if len(keystones) >= min_n:
        return ClaimResult(
            claim.id, claim.type, "satisfied", 1.0, claim.weight,
            evidence=ev,
            explanation=(
                f"{pw.display_name}: 含 {len(keystones)} 个 keystone "
                f"({', '.join(c.label for c in keystones[:3])})"
            ),
        )
    if len(keystones) >= 1:
        return ClaimResult(
            claim.id, claim.type, "partial", 0.5, claim.weight,
            evidence=ev,
            explanation=(
                f"{pw.display_name}: 仅 {len(keystones)} keystone (要求 ≥{min_n})"
            ),
        )
    return ClaimResult(
        claim.id, claim.type, "unsatisfied", 0.0, claim.weight,
        evidence=ev,
        explanation=f"{pw.display_name}: 无 keystone contributor",
    )


_DISPATCH = {
    "pathway_active": _eval_pathway_active,
    "coupling_possible": _eval_coupling_possible,
    "env_correlation": _eval_env_correlation,
    "keystone_in_pathway": _eval_keystone_in_pathway,
}


# =============================================================================
# 聚合入口
# =============================================================================

def score(hypothesis: Hypothesis, data: CycleData) -> HypothesisScore:
    """评估假说对 CycleData 的支持度。"""
    results: list[ClaimResult] = []
    for claim in hypothesis.claims:
        evaluator = _DISPATCH.get(claim.type)
        if evaluator is None:
            results.append(ClaimResult(
                claim.id, claim.type, "skipped", 0.0, claim.weight,
                evidence={"reason": "unknown claim type"},
                explanation=f"未知 claim type: {claim.type}",
            ))
            continue
        try:
            results.append(evaluator(claim, data))
        except Exception as e:  # noqa: BLE001 - 评估器不应让整个评分崩掉
            results.append(ClaimResult(
                claim.id, claim.type, "skipped", 0.0, claim.weight,
                evidence={"error": str(e)},
                explanation=f"评估出错（已跳过）：{e}",
            ))

    scored = [r for r in results if r.status != "skipped"]
    total_w = sum(r.weight for r in scored)
    earned = sum(r.weight * r.score for r in scored)
    overall = (earned / total_w) if total_w > 0 else 0.0

    if total_w == 0:
        label = "insufficient"
    elif overall >= hypothesis.strong_threshold:
        label = "strong"
    elif overall >= hypothesis.suggestive_threshold:
        label = "suggestive"
    elif overall > 0:
        label = "weak"
    else:
        label = "insufficient"

    n_satisfied = sum(1 for r in scored if r.status == "satisfied")
    n_total = len(scored)
    n_skipped = len(results) - n_total
    return HypothesisScore(
        hypothesis_name=hypothesis.name,
        overall_score=round(overall, 3),
        label=label,
        n_satisfied=n_satisfied,
        n_total=n_total,
        n_skipped=n_skipped,
        claim_results=results,
        params={
            "strong_threshold": hypothesis.strong_threshold,
            "suggestive_threshold": hypothesis.suggestive_threshold,
        },
    )
