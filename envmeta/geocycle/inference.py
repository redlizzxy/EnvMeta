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
    flat_ko_map, ko_complex_map, ko_substrate_product_map, load_kb,
    pathway_display, pathway_element_map,
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
    "sensitivity_thresholds": [30.0, 50.0, 70.0],  # S1：完整度敏感度扫描
    "perm_n": 999,                    # S2：置换检验次数
    "perm_seed": 123,                 # S2：置换随机种子
    "group_filter": None,             # S2.5-4：None/All 或组名（"CK"/"A"/"B"）
    "contributor_ranking": "abundance",  # S2.5-9：MAG 选择判据
    # 可选："abundance" / "completeness" / "keystone_only" / "keystone_priority"
}


_VALID_RANKINGS = {"abundance", "completeness", "keystone_only", "keystone_priority"}


def _contributor_sort_key(c, ranking: str):
    """按 ranking 给 MAGContribution 打分（越大越排前）。"""
    base = c.completeness * np.log1p(c.abundance_mean)
    if ranking == "completeness":
        return c.completeness
    if ranking == "keystone_priority":
        return base * (10.0 if c.is_keystone else 1.0)
    return base   # abundance / keystone_only 用同一个 base


def _apply_group_filter(
    group: str | None,
    metadata_df: pd.DataFrame | None,
    abundance_df: pd.DataFrame | None,
    env_df: pd.DataFrame | None,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None, pd.DataFrame | None]:
    """按 Group 过滤 metadata/abundance/env。缺 metadata 或 Group 列则不过滤。"""
    if group is None or str(group).lower() in ("all", "", "none"):
        return metadata_df, abundance_df, env_df
    if metadata_df is None or metadata_df.empty or "Group" not in metadata_df.columns:
        return metadata_df, abundance_df, env_df

    sample_col = ("SampleID" if "SampleID" in metadata_df.columns
                  else metadata_df.columns[0])
    keep_samples = set(
        metadata_df.loc[metadata_df["Group"].astype(str) == str(group), sample_col]
        .astype(str)
    )
    if not keep_samples:
        return metadata_df, abundance_df, env_df  # 找不到组名就原样返回

    md_f = metadata_df[metadata_df[sample_col].astype(str).isin(keep_samples)].copy()

    ab_f = abundance_df
    if abundance_df is not None and not abundance_df.empty:
        mag_col = next(
            (c for c in abundance_df.columns if c.lower() in ("mag", "bin_id", "bin")),
            abundance_df.columns[0],
        )
        keep_cols = [mag_col] + [c for c in abundance_df.columns
                                  if c != mag_col and c in keep_samples]
        if len(keep_cols) > 1:
            ab_f = abundance_df[keep_cols].copy()

    env_f = env_df
    if env_df is not None and not env_df.empty:
        env_sc = ("SampleID" if "SampleID" in env_df.columns
                  else env_df.columns[0])
        mask = env_df[env_sc].astype(str).isin(keep_samples)
        if mask.any():
            env_f = env_df[mask].copy()
        # 若 env 用组级别写法（SampleID=CK/A/B），按 Group 值过滤
        elif "Group" in env_df.columns:
            mask2 = env_df["Group"].astype(str) == str(group)
            if mask2.any():
                env_f = env_df[mask2].copy()
    return md_f, ab_f, env_f


def _permutation_rho_p(x: np.ndarray, y: np.ndarray,
                        n: int = 999, seed: int = 123) -> tuple[float, float]:
    """对 (x, y) 做 Spearman 相关 + 置换零假设 empirical p-value。

    置换方法：保持 x 不动，随机打乱 y，重算 Spearman ρ，n 次后
    empirical p = (sum(|perm_rho| >= |obs_rho|) + 1) / (n + 1)
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs_rho, _ = sp_stats.spearmanr(x, y)
    if np.isnan(obs_rho):
        return float("nan"), 1.0
    rng = np.random.default_rng(seed)
    count = 0
    abs_obs = abs(obs_rho)
    for _ in range(n):
        y_perm = rng.permutation(y)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            perm_rho, _ = sp_stats.spearmanr(x, y_perm)
        if np.isnan(perm_rho):
            continue
        if abs(perm_rho) >= abs_obs:
            count += 1
    emp_p = (count + 1) / (n + 1)
    return float(obs_rho), float(emp_p)


def _confidence_label(rho: float, perm_p: float,
                       sensitivity_robust: bool = True) -> str:
    """基于 |ρ|、置换 p、敏感度稳健性打可信度标签。

    strong       : |ρ|>0.7, perm_p<0.01, robust sensitivity
    suggestive   : 0.5<|ρ|≤0.7, perm_p<0.05
    weak         : 0.3<|ρ|≤0.5, perm_p<0.1
    spurious?    : |ρ|>=0.5 但 perm_p>0.05（相关大但置换验证不支持）
    none         : 其余
    """
    if np.isnan(rho) or np.isnan(perm_p):
        return "unknown"
    a = abs(rho)
    if a > 0.7 and perm_p < 0.01 and sensitivity_robust:
        return "strong"
    if a > 0.5 and perm_p < 0.05:
        return "suggestive"
    if a >= 0.5 and perm_p >= 0.05:
        return "spurious?"
    if 0.3 < a <= 0.5 and perm_p < 0.1:
        return "weak"
    return "none"


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
            tax["Species"] = tax[cls_col].apply(
                lambda s: next((p[3:] for p in str(s).split(";")
                                if p.strip().startswith("s__")), "") or ""
            )
        else:
            tax["Phylum"] = "Unknown"; tax["Genus"] = ""; tax["Species"] = ""
        all_mags.update(tax["MAG"].astype(str))
        base = tax[["MAG", "Phylum", "Genus", "Species"]].copy()
    else:
        base = pd.DataFrame({"MAG": sorted(all_mags),
                             "Phylum": "Unknown", "Genus": "",
                             "Species": ""})

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
    ranking: str = "abundance",
) -> dict[str, PathwayActivity]:
    pw_display = pathway_display(lang="en")
    pw_elem = pathway_element_map()
    kb = load_kb()
    reactions = {pw_id: pw.get("reaction")
                 for el in kb["elements"].values()
                 for pw_id, pw in el["pathways"].items()}
    ko_flat = flat_ko_map(kb)                     # {ko: (gene_name, pw_id, elem)}
    ko_sub_prod = ko_substrate_product_map(kb)    # {ko: {substrate, product}}
    ko_complex = ko_complex_map(kb)               # {ko: complex_id|None}

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
            species = str(m["Species"]) if (m is not None and "Species" in m.index) else ""
            # S2.5-13: 标签格式 Genus + species | Genus sp. Mx_XX | MAG_id
            if genus and species:
                label = f"{genus} {species}"
            elif genus:
                mag_tail = str(mag).split("_")[-1] if "_" in str(mag) else str(mag)
                label = f"{genus} sp. Mx_{mag_tail}"
            else:
                label = str(mag)
            is_keystone = bool(m["is_keystone"]) if (
                m is not None and "is_keystone" in m.index) else False
            # KO 级级联：按 KB 登记顺序保留该 MAG 实际持有的通路 KO
            genes_list: list[dict] = []
            for ko in kos:  # kos 保持 KB 里 pathway.genes 的插入序
                if ko in owned:
                    sp = ko_sub_prod.get(ko, {})
                    name = ko_flat[ko][0] if ko in ko_flat else ko
                    genes_list.append({
                        "ko": ko,
                        "name": name,
                        "substrate": sp.get("substrate"),
                        "product": sp.get("product"),
                        "complex": ko_complex.get(ko),
                    })
            contributors.append(MAGContribution(
                mag=mag, phylum=phy, completeness=comp,
                abundance_mean=ab, label=label,
                genes=genes_list,
                is_keystone=is_keystone,
            ))
        # S2.5-9: keystone_only 过滤（只保留 is_keystone 承载者）
        if ranking == "keystone_only":
            contributors = [c for c in contributors if c.is_keystone]
        # 按所选 ranking 排序
        contributors.sort(
            key=lambda c: _contributor_sort_key(c, ranking),
            reverse=True,
        )
        # total_contribution / mean_completeness 始终基于过滤后的 contributors
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
    perm_n: int = 999,
    perm_seed: int = 123,
) -> tuple[list[EnvCorrelation], list[EnvCorrelation]]:
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

            # 置换检验（S2）：打乱 env 重算 ρ 的经验 p
            try:
                _, perm_p = _permutation_rho_p(
                    per_sample, env_vec, n=perm_n, seed=perm_seed,
                )
            except Exception:
                perm_p = float("nan")

            ec = EnvCorrelation(
                pathway_id=pw_id, env_factor=factor,
                rho=float(rho), p_value=float(pv),
                n_samples=len(common),
                perm_p=float(perm_p) if not np.isnan(perm_p) else None,
                confidence="unknown",   # 先占位，外部用 sensitivity 信息打标签
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

    # S2.5-4：按组过滤样本级输入（metadata/abundance/env）
    metadata_df, abundance_df, env_df = _apply_group_filter(
        p.get("group_filter"), metadata_df, abundance_df, env_df,
    )

    mag_kos = _parse_ko_annotation(ko_annotation_df)
    base = _classify_mags(ko_annotation_df, taxonomy_df, keystone_df, abundance_df)
    pw_kos = pathway_ko_sets()

    ranking = p.get("contributor_ranking", "abundance")
    if ranking not in _VALID_RANKINGS:
        raise ValueError(
            f"contributor_ranking={ranking!r} invalid；可选 {sorted(_VALID_RANKINGS)}"
        )
    activities = _pathway_activity(
        mag_kos, base, pw_kos,
        threshold=p["completeness_threshold"],
        top_n=p["top_n_contributors"],
        ranking=ranking,
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
        perm_n=p["perm_n"], perm_seed=p["perm_seed"],
    )

    # 敏感度扫描
    sensitivity = _sensitivity_scan(
        mag_kos, base, pw_kos,
        thresholds=p["sensitivity_thresholds"],
        pw_elem=pw_elem, pw_display=pathway_display(lang="en"),
    )

    # 给每条 env_correlation 打可信度标签（S2），参考其对应通路的敏感度稳健性
    robust_map = {s.pathway_id: s.robust for s in sensitivity}
    for ec in env_corrs + full_corr:
        robust = robust_map.get(ec.pathway_id, True)
        pp = ec.perm_p if ec.perm_p is not None else ec.p_value
        ec.confidence = _confidence_label(ec.rho, pp, sensitivity_robust=robust)

    meta = {
        "n_mags": len(base),
        "contributor_ranking": ranking,
        "group_filter": p.get("group_filter"),
        "n_samples_used": (
            abundance_df.shape[1] - 1 if abundance_df is not None and
            not abundance_df.empty else 0
        ),
        "n_pathways_total": len(pw_kos),
        "n_pathways_active": sum(
            1 for a in activities.values() if a.n_active_mags > 0),
        "n_env_correlations": len(env_corrs),
        "n_full_corr": len(full_corr),
        "n_sensitivity_robust": sum(1 for s in sensitivity if s.robust),
        "n_confidence_strong": sum(1 for e in env_corrs if e.confidence == "strong"),
        "n_confidence_suggestive": sum(
            1 for e in env_corrs if e.confidence == "suggestive"),
        "n_confidence_spurious": sum(
            1 for e in env_corrs if e.confidence == "spurious?"),
    }

    return CycleData(
        elements=list(elements.values()),
        env_correlations=env_corrs,
        full_corr_matrix=full_corr,
        sensitivity=sensitivity,
        params=p,
        meta=meta,
    )
