"""生物地球化学循环图（Phase 3 v1）。

把 geocycle.inference + geocycle.renderer 包装为和其他 analysis 模块一致的
`analyze()` 接口，返回 AnalysisResult。

输入：
    ko_annotation_df: MAG + KEGG_ko
    taxonomy_df:      可选 MAG + classification
    keystone_df:      可选 MAG + Genus
    abundance_df:     可选 MAG × sample
    env_df:           可选 SampleID + Group + 环境因子
    metadata_df:      可选 SampleID + Group

输出：
    figure — 2×2 元素象限 + env 相关性面板
    stats  — 扁平长表（type=pathway / env_correlation）
"""
from __future__ import annotations

import pandas as pd

from envmeta.analysis.base import AnalysisResult
from envmeta.geocycle.inference import DEFAULTS as INF_DEFAULTS, infer
from envmeta.geocycle.renderer import DEFAULTS as REN_DEFAULTS, render


def _compute_most_active_pathways(
    ko_annotation_df, taxonomy_df, keystone_df, abundance_df,
    env_df, metadata_df, inf_params, current_group,
) -> set[str]:
    """对所有组跑推断，找每条通路 total_contribution 最大的组；
    返回当前组"最活"的 pathway_id 集合。

    metadata 无 Group 列 → 返回空 set。
    """
    if (metadata_df is None or metadata_df.empty
            or "Group" not in metadata_df.columns):
        return set()
    from envmeta.analysis.cycle_compare import compare_groups
    cmp = compare_groups(
        ko_annotation_df, taxonomy_df, keystone_df,
        abundance_df, env_df, metadata_df,
        params={**inf_params, "group_filter": None},
    )
    if cmp.empty:
        return set()
    best = cmp.loc[cmp.groupby("pathway_id")["total_contribution"].idxmax()]
    return set(best.loc[best["group"] == str(current_group), "pathway_id"])


def analyze(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    env_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> AnalysisResult:
    p = {**INF_DEFAULTS, **REN_DEFAULTS, **(params or {})}

    # 推断参数 vs 渲染参数
    inf_keys = set(INF_DEFAULTS.keys())
    inf_params = {k: v for k, v in p.items() if k in inf_keys}
    ren_params = {k: v for k, v in p.items() if k not in inf_keys}

    data = infer(
        ko_annotation_df, taxonomy_df, keystone_df, abundance_df,
        env_df, metadata_df,
        params=inf_params,
    )

    # S2.5-8：跨组"最活"标注（仅在单组模式 + 明确开启时做）
    if (p.get("annotate_cross_group")
            and inf_params.get("group_filter")
            and str(inf_params["group_filter"]).lower() not in ("all", "none")):
        try:
            most_active = _compute_most_active_pathways(
                ko_annotation_df, taxonomy_df, keystone_df,
                abundance_df, env_df, metadata_df,
                inf_params, inf_params["group_filter"],
            )
            ren_params["most_active_pathways"] = most_active
        except Exception:
            ren_params["most_active_pathways"] = set()

    fig = render(data, params=ren_params)
    return AnalysisResult(
        figure=fig,
        stats=data.to_flat_stats(),
        params=p,
    )
