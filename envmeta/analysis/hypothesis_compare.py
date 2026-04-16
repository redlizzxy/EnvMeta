"""跨组假说评分对比（S3.5-ui 扩展）。

对每组分别跑 infer + score，汇总成一张对比表。
这是论文 Results 里"B 组独特性"叙事的直接素材。

用法：
    df = score_by_groups(hyp, ko, tax, ks, ab, env, md)
    # df 每行一个组，列含 overall / label / null_p / robust / interpretation
"""
from __future__ import annotations

import pandas as pd

from envmeta.analysis.cycle_compare import compare_groups
from envmeta.geocycle.hypothesis import Hypothesis, HypothesisScore, score
from envmeta.geocycle.inference import infer


def _interpret(hs: HypothesisScore) -> str:
    """一行人话解读（与 app.py._interpret_hyp_score 语义对齐，纯文本版）。"""
    overall = hs.overall_score
    label = hs.label
    p = hs.null_p
    robust = hs.weight_robust
    vetoed = bool(hs.veto_reasons)

    if vetoed:
        return "VETO: 必要前提失败，假说未通过"

    specific = (p is not None and p < 0.05)
    degenerate_ok = (p is None and overall >= 0.95)

    if label == "strong":
        if (specific or degenerate_ok) and robust:
            return "最强支持：label=strong + 权重-数据一致 + robust"
        if robust and p is not None and p >= 0.20:
            return (
                f"strong 但不特异（null_p={p:.2f} ≥ 0.20，通过率运气主导）"
            )
        if not robust:
            return "strong 但权重敏感（±20% 扰动下 label 翻转）"
        if p is not None and 0.05 <= p < 0.20:
            return f"strong 边界（null_p={p:.2f}）"
        return "strong"

    if label == "suggestive":
        if specific or degenerate_ok:
            return "中等支持：suggestive + 权重-数据一致"
        return "中等支持：证据部分吻合"

    if label == "weak":
        return "证据薄弱：overall > 0 但低于 suggestive 阈值"

    return "证据不足：overall ≈ 0 或全 skipped"


def score_by_groups(
    hypothesis: Hypothesis,
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    env_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    groups: list[str] | None = None,
    params: dict | None = None,
    null_n: int = 499,
    null_seed: int | None = 42,
    sensitivity_delta: float = 0.2,
) -> pd.DataFrame:
    """对每组跑 infer + score，返回每组一行的 DataFrame。

    参数
    -----
    hypothesis          已加载的 Hypothesis 对象
    groups              要评分的组；None → 从 metadata_df.Group 自动列出
    null_n              每组 permutation 次数（默认 499，权衡速度）
    其余同 infer() / score()

    返回列：
      group / overall_score / label / null_p / null_p_samples /
      weight_robust / n_satisfied / n_total / n_skipped / n_veto /
      interpretation
    """
    if groups is None:
        if metadata_df is not None and "Group" in metadata_df.columns:
            groups = sorted(set(metadata_df["Group"].astype(str)))
        else:
            groups = ["All"]

    # 若 YAML 含 group_contrast claim，需要 compare_df；每组共用
    needs_compare = any(c.type == "group_contrast" for c in hypothesis.claims)
    compare_df = None
    if needs_compare:
        cmp_params = dict(params or {})
        cmp_params.pop("group_filter", None)
        try:
            compare_df = compare_groups(
                ko_annotation_df, taxonomy_df, keystone_df,
                abundance_df, env_df, metadata_df,
                params=cmp_params,
            )
        except Exception:
            compare_df = None   # 失败 → group_contrast claim 会 skipped

    rows: list[dict] = []
    for g in groups:
        gparams = {**(params or {}), "group_filter": g}
        data = infer(
            ko_annotation_df, taxonomy_df, keystone_df,
            abundance_df, env_df, metadata_df, params=gparams,
        )
        hs = score(
            hypothesis, data,
            compare_df=compare_df,
            null_n=null_n,
            null_seed=null_seed,
            sensitivity_delta=sensitivity_delta,
        )
        rows.append({
            "group": g,
            "overall_score": round(hs.overall_score, 3),
            "label": hs.label,
            "null_p": (None if hs.null_p is None else round(hs.null_p, 3)),
            "null_p_samples": hs.null_p_samples,
            "weight_robust": hs.weight_robust,
            "n_satisfied": hs.n_satisfied,
            "n_total": hs.n_total,
            "n_skipped": hs.n_skipped,
            "n_veto": len(hs.veto_reasons),
            "interpretation": _interpret(hs),
        })
    return pd.DataFrame(rows)
