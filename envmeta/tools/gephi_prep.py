"""Gephi CSV 预处理 + 格式校验。

用户痛点：在 Gephi 里要一个个删非 keystone 的标签。
本模块帮用户把 nodes/edges CSV 预处理成"Gephi 就绪"格式：
- Label 列只保留 keystone 的 Genus 标签（其余清空）
- 补 Genus 标签（4 档 fallback，与 MAG 图一致）
- 校验必需列 / 引用完整性 / 重复边

用法：
    from envmeta.tools.gephi_prep import prepare_gephi_csv, validate_gephi_format
    nodes_out, edges_out = prepare_gephi_csv(nodes, edges,
                                              taxonomy_df=tax,
                                              label_mode="keystone_only")
    warnings = validate_gephi_format(nodes_out, edges_out)
"""
from __future__ import annotations

import pandas as pd

from envmeta.analysis import _mag_common as _mc


def validate_gephi_format(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
) -> list[str]:
    """检查 CSV 是否符合 Gephi 导入要求。返回 warning/error 列表（空 = OK）。"""
    issues: list[str] = []

    # --- Nodes ---
    if nodes_df is None or nodes_df.empty:
        issues.append("[ERROR] nodes 表为空")
        return issues
    id_col = next((c for c in nodes_df.columns
                   if c.lower() in ("id", "mag", "name", "genome")), None)
    if id_col is None:
        issues.append("[ERROR] nodes 表缺少 Id / MAG / Name / Genome 列")
    else:
        n_dup = nodes_df[id_col].duplicated().sum()
        if n_dup:
            issues.append(f"[WARN] nodes 有 {n_dup} 条重复 Id")

    # --- Edges ---
    if edges_df is None or edges_df.empty:
        issues.append("[ERROR] edges 表为空")
        return issues
    src_col = next((c for c in edges_df.columns
                    if c.lower() in ("source", "from", "mag1")), None)
    tgt_col = next((c for c in edges_df.columns
                    if c.lower() in ("target", "to", "mag2")), None)
    if src_col is None or tgt_col is None:
        issues.append("[ERROR] edges 表缺少 Source + Target 列")
    else:
        # 引用完整性
        if id_col:
            node_ids = set(nodes_df[id_col].astype(str))
            edge_ids = set(edges_df[src_col].astype(str)) | set(edges_df[tgt_col].astype(str))
            orphan = edge_ids - node_ids
            if orphan:
                issues.append(
                    f"[WARN] edges 引用了 {len(orphan)} 个不在 nodes 里的 Id："
                    f" {', '.join(sorted(orphan)[:5])}"
                    + ("..." if len(orphan) > 5 else "")
                )
        # 重复边
        if src_col and tgt_col:
            edge_pairs = edges_df[[src_col, tgt_col]].apply(
                lambda r: tuple(sorted([str(r.iloc[0]), str(r.iloc[1])])), axis=1)
            n_dup_e = edge_pairs.duplicated().sum()
            if n_dup_e:
                issues.append(f"[WARN] edges 有 {n_dup_e} 条重复边")

    # Weight 列类型
    wt_col = next((c for c in edges_df.columns if c.lower() == "weight"), None)
    if wt_col:
        non_num = pd.to_numeric(edges_df[wt_col], errors="coerce").isna().sum()
        if non_num:
            issues.append(f"[WARN] edges.Weight 有 {non_num} 条非数值")

    return issues


def prepare_gephi_csv(
    nodes_df: pd.DataFrame,
    edges_df: pd.DataFrame,
    *,
    taxonomy_df: pd.DataFrame | None = None,
    keystone_df: pd.DataFrame | None = None,
    label_mode: str = "keystone_only",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """把原始 nodes/edges CSV → Gephi 导入就绪格式。

    核心操作：
    - 规范 Id / Label / Source / Target 列名（Gephi 要求首字母大写）
    - 补 Genus 标签（4 档 fallback，与其他 MAG 图一致）
    - label_mode="keystone_only" → 非 keystone Label 列清空
    - 确保 is_keystone 列为 boolean

    参数
    -----
    label_mode: "keystone_only" | "all" | "none"
        - keystone_only: 只有 keystone 的 Label 有值（推荐；省去 Gephi 手动删标签）
        - all: 所有节点都加 Label
        - none: Label 列全部清空
    """
    # --- Nodes ---
    nodes = nodes_df.copy()
    id_col = next((c for c in nodes.columns
                   if c.lower() in ("id", "mag", "name", "genome")),
                  nodes.columns[0])
    nodes = nodes.rename(columns={id_col: "Id"})
    nodes["Id"] = nodes["Id"].astype(str)

    # --- 生成 Genus 标签 ---
    # 优先用 nodes CSV 自带的 Genus/Family 列（避免依赖外部 taxonomy_df）
    # 只有 nodes 里没 Genus 时才 fallback 到 taxonomy_df
    nodes_has_genus = any(c.lower() == "genus" for c in nodes.columns)
    nodes_has_family = any(c.lower() == "family" for c in nodes.columns)

    if nodes_has_genus or nodes_has_family:
        # 直接从 nodes 已有列生成 label
        genus_col = next((c for c in nodes.columns if c.lower() == "genus"), None)
        family_col = next((c for c in nodes.columns if c.lower() == "family"), None)
        species_col = next((c for c in nodes.columns if c.lower() == "species"), None)
        phylum_col = next((c for c in nodes.columns if c.lower() == "phylum"), None)
        nodes["_genus"] = nodes[genus_col].fillna("").astype(str) if genus_col else ""
        nodes["_family"] = nodes[family_col].fillna("").astype(str) if family_col else ""
        nodes["_species"] = nodes[species_col].fillna("").astype(str) if species_col else ""
        nodes["_label"] = nodes.apply(
            lambda r: _mc.mag_display_label(
                r["Id"], r["_genus"], r["_species"], r["_family"]),
            axis=1,
        )
    elif taxonomy_df is not None and not taxonomy_df.empty:
        # fallback: 用外部 taxonomy_df 解析
        tmp = nodes[["Id"]].rename(columns={"Id": "MAG"})
        tmp = _mc.annotate_taxonomy(tmp, taxonomy_df)
        nodes["_label"] = tmp["label"].values
    else:
        nodes["_label"] = nodes["Id"]

    # keystone 标注
    if keystone_df is not None and not keystone_df.empty:
        ks = keystone_df.copy()
        ks = ks.rename(columns={_mc.mag_col(ks): "MAG"})
        ks_set = set(ks["MAG"].astype(str))
        nodes["is_keystone"] = nodes["Id"].isin(ks_set)
    elif "is_keystone" in nodes.columns:
        ks_col = nodes["is_keystone"]
        if ks_col.dtype == object:
            nodes["is_keystone"] = ks_col.str.lower().isin(("true", "1", "yes"))

    # Label 列（Gephi 大小写不敏感，只保留一个 Label）
    if label_mode == "none":
        nodes["Label"] = ""
    elif label_mode == "all":
        nodes["Label"] = nodes["_label"]
    else:  # keystone_only
        nodes["Label"] = nodes.apply(
            lambda r: r["_label"] if r.get("is_keystone", False) else "",
            axis=1,
        )

    # 清理：删除中间列 + 避免 Gephi "repeated column" 报错
    # 删除所有临时列和可能与 Label 重复的 label 列
    drop_cols = {"_genus", "_family", "_species", "_label"}
    # 删除原始的小写 label 列（若存在）避免与 Label 冲突
    if "label" in nodes.columns and "Label" in nodes.columns:
        drop_cols.add("label")
    # 删除 _x / _y 后缀列（merge 产物）
    for c in list(nodes.columns):
        if c.endswith("_x") or c.endswith("_y"):
            drop_cols.add(c)
    out_nodes = nodes.drop(columns=[c for c in drop_cols if c in nodes.columns])

    # 确保 is_keystone 可被 Gephi 识别（Boolean → string）
    if "is_keystone" in out_nodes.columns:
        out_nodes["is_keystone"] = out_nodes["is_keystone"].map(
            {True: "True", False: "False"})

    # --- Edges ---
    edges = edges_df.copy()
    # 规范列名
    col_renames = {}
    for c in edges.columns:
        cl = c.lower()
        if cl in ("source", "from", "mag1"):
            col_renames[c] = "Source"
        elif cl in ("target", "to", "mag2"):
            col_renames[c] = "Target"
        elif cl == "weight":
            col_renames[c] = "Weight"
    edges = edges.rename(columns=col_renames)

    # Type（Gephi 导入需要）
    if "Type" not in edges.columns:
        edges["Type"] = "Undirected"

    return out_nodes, edges
