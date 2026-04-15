"""跨组循环图对比（S2.5-7d）。

回答"三张单组循环图看起来一样、如何描述组间差异"：
对每组分别跑推断，把 pathway × group 结果汇成长表。差别藏在三个维度：
  1. total_contribution  — 丰度加权活性（同一通路在哪组更活跃）
  2. top_mag_genus       — 通路承载者身份（是否换了物种）
  3. top_mag_is_keystone — 承载者是否为关键物种

单独的 cycle_diagram.analyze() 接受 group_filter=None 时看全部样本；
此函数是对多个 group_filter 值的批处理 + DataFrame 汇总。
"""
from __future__ import annotations

import pandas as pd

from envmeta.geocycle.inference import infer


def compare_groups(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    env_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    groups: list[str] | None = None,
    params: dict | None = None,
) -> pd.DataFrame:
    """对每组分别跑 infer，返回 pathway × group 的长表。

    列：group, element, pathway_id, display_name, n_active_mags,
         mean_completeness, total_contribution, top_mag, top_mag_genus,
         top_mag_is_keystone
    """
    if groups is None:
        if metadata_df is not None and "Group" in metadata_df.columns:
            groups = sorted(set(metadata_df["Group"].astype(str)))
        else:
            groups = ["All"]

    keystone_mags: set[str] = set()
    if keystone_df is not None and not keystone_df.empty:
        mag_col = next((c for c in keystone_df.columns
                        if c.lower() in ("mag", "bin", "bin_id")),
                       keystone_df.columns[0])
        keystone_mags = set(keystone_df[mag_col].astype(str))

    rows: list[dict] = []
    for g in groups:
        gparams = {**(params or {}), "group_filter": g}
        data = infer(
            ko_annotation_df, taxonomy_df, keystone_df,
            abundance_df, env_df, metadata_df, params=gparams,
        )
        for ec in data.elements:
            for pw in ec.pathways:
                top = pw.contributors[0] if pw.contributors else None
                rows.append({
                    "group": g,
                    "element": ec.element_id,
                    "pathway_id": pw.pathway_id,
                    "display_name": pw.display_name,
                    "n_active_mags": pw.n_active_mags,
                    "mean_completeness": pw.mean_completeness,
                    "total_contribution": pw.total_contribution,
                    "top_mag": top.mag if top else None,
                    "top_mag_genus": top.label if top else None,
                    "top_mag_is_keystone": (
                        bool(top and top.mag in keystone_mags)
                        if top else False
                    ),
                })
    return pd.DataFrame(rows)
