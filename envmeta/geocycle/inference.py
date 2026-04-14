"""生物地球化学循环图推断引擎。

从 MAG × KO 注释（+ 可选分类 / 丰度 / 环境因子）推断哪些元素循环通路在
当前数据集里活跃、由哪些 MAG 承载，以及环境因子与通路活性的相关性。

输入（与 pathway / gene_profile 保持一致）：
    ko_annotation_df: MAG + KEGG_ko 长表
    taxonomy_df:      可选 MAG + classification（GTDB）
    keystone_df:      可选 MAG + Genus / Phylum
    abundance_df:     可选 MAG × sample 丰度
    env_df:           可选 SampleID + Group + 数值环境因子
    metadata_df:      可选 SampleID + Group（env 对齐用）

输出：
    CycleData
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from envmeta.analysis.pathway import (
    _extract_phylum, _mag_col, _parse_ko_annotation,
)
from envmeta.geocycle.knowledge_base import (
    element_colors, element_display,
    flat_ko_map, load_kb, pathway_display, pathway_element_map,
    pathway_ko_sets,
)
from envmeta.geocycle.model import (
    CycleData, ElementCycle, EnvCorrelation, MAGContribution,
    PathwayActivity, SensitivityRow,
)

DEFAULTS = {
    "completeness_threshold": 50.0,   # MAG 完整度 ≥ 此值视为"承载该通路"
    "top_n_contributors": 5,          # 每条通路保留 Top-N 贡献 MAG
    "env_rho_min": 0.5,               # 环境相关性最低 |rho|
    "env_p_max": 0.05,                # 相关性显著性阈值
    "sensitivity_thresholds": [30.0, 50.0, 70.0],  # S1 新增：完整度敏感度扫描
}


def _classify_mags(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None,
    keystone_df: pd.DataFrame | None,
    abundance_df: pd.DataFrame | None,
) -> pd.DataFrame:
    """构造 MAG 基础表：MAG | Phylum | Genus | abundance_mean | is_keystone。"""
    mag_kos = _parse_ko_annotation(ko_annotation_df)
    all_mags = set(mag_kos.keys())

    if taxonomy_df is not None and not taxonomy_df.empty:
        tax = taxonomy_df.copy()
        tax = tax.rename(columns={_mag_col(tax): "MAG"})
        cls_col = next((c for c in tax.columns
                        if "classif" in c.lower() or "taxonomy" in c.lower()),
                       tax.columns[1] if len(tax.columns) > 1 else None)
        if cls_col:
            tax["Phylum"] = tax[cls_col].apply(_extract_phylum)
            tax["Genus"] = tax[cls_col].apply(
                lambda s: next((p[3:] for p in str(s).split(";")
                                if p.strip().startswith("g__")), "") or ""
            )
        else:
            tax["Phylum"] = "Unknown"; tax["Genus"] = ""
        all_mags.update(tax["MAG"].astype(str))
        base = tax[["MAG", "Phylum", "Genus"]].copy()
    else:
        base = pd.DataFrame({"MAG": sorted(all_mags),
                             "Phylum": "Unknown", "Genus": ""})

    base = base[base["MAG"].isin(all_mags)].drop_duplicates("MAG").reset_index(drop=True)

    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={_mag_col(ks): "MAG"})
        base["is_keystone"] = base["MAG"].isin(set(ks["MAG"].astype(str)))
        if "Genus" in ks.columns:
            ks_genus = dict(zip(ks["MAG"].astype(str), ks["Genus"].astype(str)))
            base["Genus"] = base.apply(
                lambda r: ks_genus.get(r["MAG"], r["Genus"]) or r["Genus"], axis=1)
    else:
        base["is_keystone"] = False

    if abundance_df is not None and not abundance_df.empty:
        ab = abundance_df.copy()
        ab = ab.rename(columns={_mag_col(ab): "MAG"})
        scols = [c for c in ab.columns if c != "MAG"]
        ab["abundance_mean"] = ab[scols].apply(
            pd.to_numeric, errors="coerce").fillna(0).mean(axis=1)
        base = base.merge(ab[["MAG", "abundance_mean"]], on="MAG", how="left")
        base["abundance_mean"] = base["abundance_mean"].fillna(0.0)
    else:
        base["abundance_mean"] = 0.0
    return base


def _pathway_activity(
    mag_kos: dict[str, set[str]],
    base: pd.DataFrame,
    pw_kos: dict[str, list[str]],
    threshold: float,
    top_n: int,
) -> dict[str, PathwayActivity]:
    pw_display = pathway_display(lang="en")
    pw_elem = pathway_element_map()
    kb = load_kb()
    reactions = {pw_id: pw.get("reaction")
                 for el in kb["elements"].values()
                 for pw_id, pw in el["pathways"].items()}

    mag_meta = base.set_index("MAG")
    activities: dict[str, PathwayActivity] = {}
    for pw_id, kos in pw_kos.items():
        k_set = set(kos)
        contributors: list[MAGContribution] = []
        for mag, owned in mag_kos.items():
            comp = len(owned & k_set) / len(k_set) * 100 if kos else 0.0
            if comp < threshold:
                continue
            m = mag_meta.loc[mag] if mag in mag_meta.index else None
            phy = str(m["Phylum"]) if m is not None else "Unknown"
            ab = float(m["abundance_mean"]) if m is not None else 0.0
            genus = str(m["Genus"]) if (m is not None and "Genus" in m.index) else ""
            label = genus if genus else mag
            contributors.append(MAGContribution(
                mag=mag, phylum=phy, completeness=comp,
                abundance_mean=ab, label=label,
            ))
        # 按 "completeness × (log1p abundance)" 排序
        contributors.sort(
            key=lambda c: c.completeness * np.log1p(c.abundance_mean),
            reverse=True,
        )
        total = sum(c.completeness * c.abundance_mean for c in contributors)
        mean_comp = float(np.mean([c.completeness for c in contributors])) if contributors else 0.0
        activities[pw_id] = PathwayActivity(
            pathway_id=pw_id,
            element=pw_elem.get(pw_id, "unknown"),
            display_name=pw_display.get(pw_id, pw_id),
            reaction=reactions.get(pw_id),
            n_active_mags=len(contributors),
            mean_completeness=mean_comp,
            total_contribution=total,
            contributors=contributors[:top_n],
        )
    return activities


def _env_correlations(
    mag_kos: dict[str, set[str]],
    abundance_df: pd.DataFrame | None,
    env_df: pd.DataFrame | None,
    metadata_df: pd.DataFrame | None,
    pw_kos: dict[str, list[str]],
    rho_min: float,
    p_max: float,
) -> list[EnvCorrelation]:
    """对每条通路计算"KO 丰度加权总和" vs 每个环境因子的 Spearman。

    样本级通路活性 = Σ_MAG (abundance_MAG × |MAG 的 KO ∩ 通路 KO|)。
    """
    if (env_df is None or env_df.empty
            or abundance_df is None or abundance_df.empty):
        return [], []
    if metadata_df is None:
        return [], []

    ab = abundance_df.copy()
    ab = ab.rename(columns={_mag_col(ab): "MAG"})
    sample_cols = [c for c in ab.columns if c != "MAG"]
    ab_mat = ab.set_index("MAG")[sample_cols].apply(
        pd.to_numeric, errors="coerce").fillna(0.0)

    env_sample_col = "SampleID" if "SampleID" in env_df.columns else env_df.columns[0]
    group_col = "Group" if "Group" in env_df.columns else None
    env_num_cols = [c for c in env_df.columns
                    if c not in (env_sample_col, group_col, "Replicate")
                    and pd.to_numeric(env_df[c], errors="coerce").notna().mean() > 0.9]
    if not env_num_cols:
        return [], []

    # env 对齐到 abundance 列；fallback 按 group 顺序（与 rda 同思路）
    env_ids = env_df[env_sample_col].astype(str).tolist()
    common = [c for c in sample_cols if c in set(env_ids)]
    if len(common) < 4:
        # Group 兜底
        meta_sc = "SampleID" if "SampleID" in metadata_df.columns else metadata_df.columns[0]
        mgmap = dict(zip(metadata_df[meta_sc].astype(str),
                         metadata_df.get("Group", pd.Series(dtype=str)).astype(str)))
        if group_col:
            env_per_g: dict[str, list[int]] = {}
            for i, g in enumerate(env_df[group_col].astype(str)):
                env_per_g.setdefault(g, []).append(i)
            seen: dict[str, int] = {}
            keep_cols, env_idx = [], []
            for c in sample_cols:
                g = mgmap.get(c)
                if g and g in env_per_g and seen.get(g, 0) < len(env_per_g[g]):
                    keep_cols.append(c)
                    env_idx.append(env_per_g[g][seen.get(g, 0)])
                    seen[g] = seen.get(g, 0) + 1
            common = keep_cols
            env_aligned = env_df.iloc[env_idx].reset_index(drop=True)
        else:
            return [], []
    else:
        env_aligned = env_df.set_index(env_sample_col).loc[common].reset_index()

    if len(common) < 4:
        return [], []
    ab_mat = ab_mat[common]

    filtered: list[EnvCorrelation] = []
    full: list[EnvCorrelation] = []
    for pw_id, kos in pw_kos.items():
        # 样本级活性：Σ_MAG abundance_MAG × |MAG 的 KO ∩ 通路 KO|
        per_sample = np.zeros(len(common))
        for mag in ab_mat.index:
            owned = mag_kos.get(str(mag), set())
            overlap = len(owned & set(kos))
            if overlap == 0:
                continue
            per_sample += ab_mat.loc[mag].to_numpy() * overlap
        if per_sample.sum() == 0:
            continue
        for factor in env_num_cols:
            env_vec = pd.to_numeric(env_aligned[factor], errors="coerce").to_numpy()
            try:
                rho, pv = sp_stats.spearmanr(per_sample, env_vec)
            except Exception:
                continue
            if np.isnan(rho):
                continue
            ec = EnvCorrelation(
                pathway_id=pw_id, env_factor=factor,
                rho=float(rho), p_value=float(pv),
                n_samples=len(common),
            )
            full.append(ec)   # 所有非 NaN 都进完整矩阵
            if abs(rho) >= rho_min and pv <= p_max:
                filtered.append(ec)
    return filtered, full


def _sensitivity_scan(
    mag_kos: dict[str, set[str]],
    base: pd.DataFrame,
    pw_kos: dict[str, list[str]],
    thresholds: list[float],
    pw_elem: dict[str, str],
    pw_display: dict[str, str],
) -> list[SensitivityRow]:
    """在多档完整度阈值下扫描每条通路的 Top-1 contributor 稳定性。

    Top-1 = 活跃 MAG 里 completeness × log1p(abundance) 最高者的 label。
    robust = 所有阈值档下 Top-1 一致。
    """
    mag_meta = base.set_index("MAG")
    rows: list[SensitivityRow] = []
    for pw_id, kos in pw_kos.items():
        k_set = set(kos)
        top1_list: list[str | None] = []
        n_active_list: list[int] = []
        for thr in thresholds:
            contributors = []
            for mag, owned in mag_kos.items():
                comp = len(owned & k_set) / len(k_set) * 100 if kos else 0.0
                if comp < thr:
                    continue
                m = mag_meta.loc[mag] if mag in mag_meta.index else None
                ab = float(m["abundance_mean"]) if m is not None else 0.0
                genus = str(m["Genus"]) if (m is not None and "Genus" in m.index) else ""
                label = genus if genus else mag
                contributors.append((comp * np.log1p(ab), label))
            contributors.sort(reverse=True)
            top1_list.append(contributors[0][1] if contributors else None)
            n_active_list.append(len(contributors))
        non_none = [t for t in top1_list if t is not None]
        robust = len(set(non_none)) <= 1 and len(non_none) == len(top1_list)
        rows.append(SensitivityRow(
            pathway_id=pw_id, element=pw_elem.get(pw_id, "unknown"),
            display_name=pw_display.get(pw_id, pw_id),
            thresholds=list(thresholds),
            top1_by_threshold=top1_list,
            n_active_by_threshold=n_active_list,
            robust=robust,
        ))
    return rows


def infer(
    ko_annotation_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    abundance_df: pd.DataFrame | None = None,
    env_df: pd.DataFrame | None = None,
    metadata_df: pd.DataFrame | None = None,
    params: dict | None = None,
) -> CycleData:
    p = {**DEFAULTS, **(params or {})}

    mag_kos = _parse_ko_annotation(ko_annotation_df)
    base = _classify_mags(ko_annotation_df, taxonomy_df, keystone_df, abundance_df)
    pw_kos = pathway_ko_sets()

    activities = _pathway_activity(
        mag_kos, base, pw_kos,
        threshold=p["completeness_threshold"],
        top_n=p["top_n_contributors"],
    )

    # 按元素组织
    elem_disp = element_display(lang="en")
    elem_color = element_colors()
    elements: dict[str, ElementCycle] = {}
    pw_elem = pathway_element_map()
    for pw_id, act in activities.items():
        eid = pw_elem.get(pw_id, "unknown")
        if eid not in elements:
            elements[eid] = ElementCycle(
                element_id=eid,
                display_name=elem_disp.get(eid, eid),
                color=elem_color.get(eid, "#888888"),
            )
        elements[eid].pathways.append(act)

    # 组内按 total_contribution 降序
    for ec in elements.values():
        ec.pathways.sort(key=lambda x: x.total_contribution, reverse=True)

    env_corrs, full_corr = _env_correlations(
        mag_kos, abundance_df, env_df, metadata_df, pw_kos,
        rho_min=p["env_rho_min"], p_max=p["env_p_max"],
    )

    # 敏感度扫描
    sensitivity = _sensitivity_scan(
        mag_kos, base, pw_kos,
        thresholds=p["sensitivity_thresholds"],
        pw_elem=pw_elem, pw_display=pathway_display(lang="en"),
    )

    meta = {
        "n_mags": len(base),
        "n_pathways_total": len(pw_kos),
        "n_pathways_active": sum(
            1 for a in activities.values() if a.n_active_mags > 0),
        "n_env_correlations": len(env_corrs),
        "n_full_corr": len(full_corr),
        "n_sensitivity_robust": sum(1 for s in sensitivity if s.robust),
    }

    return CycleData(
        elements=list(elements.values()),
        env_correlations=env_corrs,
        full_corr_matrix=full_corr,
        sensitivity=sensitivity,
        params=p,
        meta=meta,
    )
