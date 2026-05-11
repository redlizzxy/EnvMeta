# -*- coding: utf-8 -*-
"""
从 tests/sample_data/（168 MAG）抽取 30 MAG 子集到 tests/sample_data_demo/，
用于在线 demo（降低 Streamlit Cloud 并发崩溃风险）。

挑选规则：
1. 全部 14 keystone MAGs（保证 keystone_in_pathway claim 可评分）
2. 从剩余 MAG 中按总丰度降序补足 30 个
3. 同时保证 4 元素（As/N/S/Fe）每元素至少 5 个活跃 MAG

输出：tests/sample_data_demo/ 全套文件
"""
from __future__ import annotations
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "tests" / "sample_data"
DST = ROOT / "tests" / "sample_data_demo"
DST.mkdir(parents=True, exist_ok=True)

TARGET_N = 30


def load_kb_ko_to_element() -> dict[str, str]:
    kb = json.load(open(ROOT / "envmeta" / "geocycle" / "knowledge_base" / "elements.json",
                        encoding="utf-8"))
    out: dict[str, str] = {}
    for elem_id, elem in kb["elements"].items():
        for pw in elem["pathways"].values():
            for ko_id in pw["genes"]:
                out[ko_id] = elem_id
    return out


def pick_mags() -> list[str]:
    abund = pd.read_csv(SRC / "abundance.tsv", sep="\t")
    keystone = pd.read_csv(SRC / "keystone_species.txt", sep="\t")
    ko_long = pd.read_csv(SRC / "kegg_target_only.tsv", sep="\t")

    ko_to_elem = load_kb_ko_to_element()
    ko_long["ko_clean"] = ko_long["KEGG_ko"].str.replace("ko:", "", regex=False)
    ko_long["element"] = ko_long["ko_clean"].map(ko_to_elem)
    mag_elem = (ko_long.dropna(subset=["element"])
                .groupby("MAG")["element"].apply(set))

    # 1. 所有 keystone MAGs（与 abundance / quality 表的交集）
    ks_set = set(keystone["MAG"])
    abund_set = set(abund["Genome"])
    keystone_present = sorted(ks_set & abund_set)

    # 2. 计算非 keystone MAG 的总丰度
    abund_sum = abund.set_index("Genome").sum(axis=1).sort_values(ascending=False)
    non_ks_top = [m for m in abund_sum.index if m not in ks_set]

    picked = list(keystone_present)
    elements_target = {"arsenic", "nitrogen", "sulfur", "iron"}
    elem_count = {e: 0 for e in elements_target}
    for m in picked:
        for e in mag_elem.get(m, set()):
            if e in elem_count:
                elem_count[e] += 1

    # 3. 先补足每个元素 ≥ 5 个活跃 MAG
    for elem in elements_target:
        while elem_count[elem] < 5:
            for m in non_ks_top:
                if m in picked:
                    continue
                if elem in mag_elem.get(m, set()):
                    picked.append(m)
                    for e in mag_elem.get(m, set()):
                        if e in elem_count:
                            elem_count[e] += 1
                    break
            else:
                break  # 没找到，跳出

    # 4. 按总丰度补足到 TARGET_N
    for m in non_ks_top:
        if len(picked) >= TARGET_N:
            break
        if m not in picked:
            picked.append(m)

    print(f"Picked {len(picked)} MAGs:")
    print(f"  Keystone: {len(keystone_present)}")
    print(f"  Element coverage (count of MAGs active in element):")
    final_count = {e: 0 for e in elements_target}
    for m in picked:
        for e in mag_elem.get(m, set()):
            if e in final_count:
                final_count[e] += 1
    for e, c in final_count.items():
        print(f"    {e}: {c}")
    return picked


def subset_files(mags: list[str]) -> None:
    mags_set = set(mags)

    # abundance: filter rows
    abund = pd.read_csv(SRC / "abundance.tsv", sep="\t")
    abund_sub = abund[abund["Genome"].isin(mags_set)].reset_index(drop=True)
    abund_sub.to_csv(DST / "abundance.tsv", sep="\t", index=False)
    print(f"  abundance.tsv: {abund.shape} -> {abund_sub.shape}")

    # quality_report
    qual = pd.read_csv(SRC / "quality_report.tsv", sep="\t")
    qual_sub = qual[qual["Name"].isin(mags_set)].reset_index(drop=True)
    qual_sub.to_csv(DST / "quality_report.tsv", sep="\t", index=False)
    print(f"  quality_report.tsv: {qual.shape} -> {qual_sub.shape}")

    # mag_taxonomy_labels (no header, columns 0=MAG 1=classification)
    tax = pd.read_csv(SRC / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "classification"])
    tax_sub = tax[tax["MAG"].isin(mags_set)].reset_index(drop=True)
    tax_sub.to_csv(DST / "mag_taxonomy_labels.tsv", sep="\t",
                   index=False, header=False)
    print(f"  mag_taxonomy_labels.tsv: {tax.shape} -> {tax_sub.shape}")

    # kegg_target_only (long format, filter MAG column)
    ko_long = pd.read_csv(SRC / "kegg_target_only.tsv", sep="\t")
    ko_sub = ko_long[ko_long["MAG"].isin(mags_set)].reset_index(drop=True)
    ko_sub.to_csv(DST / "kegg_target_only.tsv", sep="\t", index=False)
    print(f"  kegg_target_only.tsv: {ko_long.shape} -> {ko_sub.shape}")

    # keystone_species (filter MAG column)
    keystone = pd.read_csv(SRC / "keystone_species.txt", sep="\t")
    ks_sub = keystone[keystone["MAG"].isin(mags_set)].reset_index(drop=True)
    ks_sub.to_csv(DST / "keystone_species.txt", sep="\t", index=False)
    print(f"  keystone_species.txt: {keystone.shape} -> {ks_sub.shape}")

    # gephi_nodes (Id column)
    if (SRC / "gephi_nodes.csv").exists():
        nodes = pd.read_csv(SRC / "gephi_nodes.csv")
        nodes_sub = nodes[nodes["Id"].isin(mags_set)].reset_index(drop=True)
        nodes_sub.to_csv(DST / "gephi_nodes.csv", index=False)
        print(f"  gephi_nodes.csv: {nodes.shape} -> {nodes_sub.shape}")

    # gephi_edges (Source / Target both must be in mags_set)
    if (SRC / "gephi_edges.csv").exists():
        edges = pd.read_csv(SRC / "gephi_edges.csv")
        src_col = "Source" if "Source" in edges.columns else edges.columns[0]
        tgt_col = "Target" if "Target" in edges.columns else edges.columns[1]
        edges_sub = edges[edges[src_col].isin(mags_set) & edges[tgt_col].isin(mags_set)].reset_index(drop=True)
        edges_sub.to_csv(DST / "gephi_edges.csv", index=False)
        print(f"  gephi_edges.csv: {edges.shape} -> {edges_sub.shape}")

    # Sample-level files (no MAG filtering needed; copy as-is)
    import shutil
    for fname in ["metadata.txt", "env_factors.txt", "alpha.txt",
                  "beta_bray.txt", "ko_tpm.spf", "07_element_KO_heatmap.txt",
                  "Phylum.txt", "Genus.txt", "Species.txt",
                  "sample_hypothesis.yaml"]:
        src_f = SRC / fname
        if src_f.exists():
            shutil.copy2(src_f, DST / fname)
            print(f"  copied: {fname}")


def write_readme(picked: list[str]) -> None:
    text = f"""# Sample data — demo subset (30 MAGs)

> **本目录是在线 demo 用的轻量化测试样本**，由 `tests/sample_data/` 完整 168 MAG
> 数据集抽取 30 MAG 而成。
>
> 抽样规则：(1) 保留全部 keystone MAGs；(2) 保证 4 元素（As/N/S/Fe）每元素 ≥ 5
> 个活跃 MAG；(3) 按总丰度补足到 30 个。

## 用途

- **在线 demo**（Streamlit Cloud）：app.py 默认从这里加载示例数据，减少
  并发崩溃风险（30 MAG 内存占用 ≈ 168 MAG 的 1/5-1/6）。
- **本目录不用于**：单元测试（pytest 跑 `tests/sample_data/`）/ 论文 Arm A 校准
  实验 / 性能 benchmark — 这些场景都用完整 168 MAG 数据集。

## 仅作功能演示，不是科学论断的数据

本目录数据经过 MAG 子采样，**不反映原研究的科学结论**。仅供用户在线点几下
"加载示例数据 → 生成图" 看看 EnvMeta 各功能 UI 是否正常。要复现砷渣修复研究
的结果，请用完整 168 MAG 数据集（`tests/sample_data/`）或下载论文 Fork Bundle。

## 生成方式

运行：

```
python tests/sample_data_demo/_build_demo_subset.py
```

会从 `tests/sample_data/` 重新抽取并覆盖本目录所有文件。

## 文件清单（与 tests/sample_data/ 对应）

- MAG-keyed（已抽样）：abundance.tsv / quality_report.tsv /
  mag_taxonomy_labels.tsv / kegg_target_only.tsv / keystone_species.txt /
  gephi_nodes.csv / gephi_edges.csv
- Sample-keyed（原样拷贝）：metadata.txt / env_factors.txt / alpha.txt /
  beta_bray.txt / ko_tpm.spf / 07_element_KO_heatmap.txt
- Taxon-aggregated（原样拷贝）：Phylum.txt / Genus.txt / Species.txt
- 假说 YAML 模板（原样拷贝）：sample_hypothesis.yaml

## 被抽样的 30 个 MAG

{chr(10).join(f"- {m}" for m in picked)}
"""
    (DST / "README.md").write_text(text, encoding="utf-8")
    print(f"  wrote README.md (subset list)")


if __name__ == "__main__":
    print("=== Picking MAGs ===")
    picked = pick_mags()
    print()
    print("=== Writing subset files ===")
    subset_files(picked)
    print()
    write_readme(picked)
    print()
    print(f"Done. {len(picked)} MAGs written to {DST}")
