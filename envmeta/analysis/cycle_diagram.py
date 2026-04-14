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
    fig = render(data, params=ren_params)
    return AnalysisResult(
        figure=fig,
        stats=data.to_flat_stats(),
        params=p,
    )
