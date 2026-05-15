"""
Liu et al. 2023 npj Biofilms (DOI 10.1038/s41522-023-00382-8) Suppl xlsx
→ EnvMeta input files

Source license: CC BY 4.0（Suppl 数据）— 比 Wei (CC BY-NC-ND) 宽松
但仍按对照实验 pre-registration 原则不入 git repo，本地 reshape 再复现。

用法:
    python tools/external_benchmarks/liu2023_reshape.py \\
        --moesm2 D:\\download\\41522_2023_382_MOESM2_ESM.xlsx \\
        --moesm4 D:\\download\\41522_2023_382_MOESM4_ESM.xlsx \\
        --moesm5 D:\\download\\41522_2023_382_MOESM5_ESM.xlsx \\
        --moesm6 D:\\download\\41522_2023_382_MOESM6_ESM.xlsx \\
        --out paper/benchmarks/external/liu_2023_coldseep/input_data_local

输出 6 个文件（与 Wei reshape 完全同 schema）:
    metadata.tsv             # 87 samples × {SampleID, Group="All", Replicate}
    env_factors.tsv          # 占位 env (深海冷泉 MOESM2 没发 pH/化学元素)
    mag_taxonomy_labels.tsv  # 1084 As-cycling MAGs × taxonomy
    quality_report.tsv       # 1084 MAGs × CheckM-style schema
    abundance.tsv            # 1084 MAGs × 87 samples（per-sample，无维度退化）
    kegg_target_only.tsv     # As gene MAG-level presence

注: Liu 数据规模与 Wei 不同：
- 1741 全集 MAG（MOESM4），1084 As-cycling 子集发布丰度（MOESM5）→ 取 1084 对齐
- Group="All" 单组（Liu 论文是 pooled 分析所有 87 metagenomes，无 case/control 设计）
- arxA / arsP / arsI 跳过（KB v1.1 未收录，与 Wei 一致策略）
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Liu 11 个 As 基因 → EnvMeta KB v1.1 KO 映射
GENE_TO_KO = {
    "aioA":     "K08356",
    "arxA":     None,        # KB v1.1 未收录（与 Wei 一致）
    "arrA":     "K28466",
    "arsC1":    "K00537",
    "arsC2":    "K03741",
    "arsM":     "K07755",    # As methylation
    "acr3":     "K03325",    # As transport
    "arsB":     "K03893",    # As transport
    "arsB_935": "K03893",    # arsB variant
    "arsP":     None,        # KB v1.1 未收录（KB 有 arsJ K25223 不同基因）
    "arsI":     None,        # KB v1.1 未收录
    "arsH":     "K11811",    # As transport/detox
}

GENE_DESC = {
    "aioA":     "arsenite oxidase large subunit",
    "arxA":     "anaerobic arsenite oxidase (KB unmapped)",
    "arrA":     "respiratory arsenate reductase",
    "arsC1":    "arsenate reductase (glutaredoxin coupled)",
    "arsC2":    "arsenate reductase (thioredoxin coupled)",
    "arsM":     "arsenite S-adenosylmethionine methyltransferase",
    "acr3":     "arsenite efflux pump (ACR3 family)",
    "arsB":     "arsenical pump membrane protein",
    "arsB_935": "arsenical pump variant 935",
    "arsP":     "putative methylarsenite efflux permease (KB unmapped)",
    "arsI":     "C-As lyase organoarsenical lyase (KB unmapped)",
    "arsH":     "arsenical resistance protein H",
}


def reshape_metadata(moesm5: pd.DataFrame) -> pd.DataFrame:
    samples = [c for c in moesm5.columns if c not in ("Bins", "Taxonomy")]
    return pd.DataFrame({
        "SampleID": samples,
        "Group": ["All"] * len(samples),
        "Replicate": list(range(1, len(samples) + 1)),
    })


def reshape_env(samples: list[str]) -> pd.DataFrame:
    """深海冷泉 MOESM2 metadata 太分级（站点级 vs sample 级），且没发布 pH/化学
    元素等数据（深海传感器测的是 depth/temp/methane/H2S 等，未公开）。
    返回最小 env 表（仅 SampleID + Group + Depth_proxy=0），EnvMeta 仍可跑 cycle_diagram。
    """
    return pd.DataFrame({
        "SampleID": samples,
        "Group": ["All"] * len(samples),
        "Depth_proxy": [0.0] * len(samples),
    })


def reshape_taxonomy(moesm5: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "MAG_ID": moesm5["Bins"],
        "taxonomy": moesm5["Taxonomy"],
    })


def reshape_quality(moesm4: pd.DataFrame, target_mags: set) -> pd.DataFrame:
    sub = moesm4[moesm4["Genomes"].isin(target_mags)].copy()
    n = len(sub)
    return pd.DataFrame({
        "Name": sub["Genomes"].values,
        "Completeness": pd.to_numeric(sub["Completeness"], errors="coerce").values,
        "Contamination": pd.to_numeric(sub["Contamination"], errors="coerce").values,
        "Completeness_Model_Used": ["Liu2023_published"] * n,
        "Translation_Table_Used": [11] * n,
        "Coding_Density": [0.92] * n,
        "Contig_N50": [0] * n,
        "Average_Gene_Length": [0] * n,
        "Genome_Size": [0] * n,
        "GC_Content": [0.0] * n,
        "Total_Coding_Sequences": [0] * n,
        "Total_Contigs": [0] * n,
        "Max_Contig_Length": [0] * n,
        "Additional_Notes": ["Liu2023 MOESM4"] * n,
    })


def reshape_abundance(moesm5: pd.DataFrame) -> pd.DataFrame:
    df = moesm5.drop(columns=["Taxonomy"]).copy()
    df = df.rename(columns={"Bins": "Genome"})
    return df


def reshape_ko_long(moesm6: pd.DataFrame) -> pd.DataFrame:
    rows = []
    skipped: list[tuple[str, str]] = []
    for _, r in moesm6.iterrows():
        gene = r.get("Genes")
        mag = r.get("Genomes")
        if pd.isna(gene) or pd.isna(mag):
            continue
        ko = GENE_TO_KO.get(gene)
        if ko is None:
            skipped.append((mag, gene))
            continue
        rows.append({
            "MAG": mag,
            "Gene_ID": str(r.get("Gene-ID", ""))[:80],
            "KEGG_ko": f"ko:{ko}",
            "Description": GENE_DESC.get(gene, ""),
        })
    if skipped:
        print(f"  [Skipped (KB v1.1 unmapped)] {len(skipped)} records "
              f"across {len({s[0] for s in skipped})} MAGs")
        for gene in sorted({s[1] for s in skipped}):
            n = sum(1 for s in skipped if s[1] == gene)
            print(f"    {gene}: {n} hits")
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--moesm2", required=True)
    p.add_argument("--moesm4", required=True)
    p.add_argument("--moesm5", required=True)
    p.add_argument("--moesm6", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Reading MOESM2/4/5/6...")
    moesm5 = pd.read_excel(args.moesm5, sheet_name=0, header=1)
    moesm4 = pd.read_excel(args.moesm4, sheet_name=0, header=1)
    moesm6 = pd.read_excel(args.moesm6, sheet_name=0, header=1)
    # MOESM2 暂不深度解析（站点级 metadata 太复杂；此次只用 87 sample 列名）

    samples = [c for c in moesm5.columns if c not in ("Bins", "Taxonomy")]
    print(f"  {len(samples)} samples in MOESM5")
    print(f"  {len(moesm5)} As-cycling MAGs in MOESM5")
    print(f"  {len(moesm4)} all MAGs in MOESM4")
    print(f"  {len(moesm6)} As gene hits in MOESM6")

    metadata = reshape_metadata(moesm5)
    metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)
    print(f"\nWriting outputs:")
    print(f"  metadata.tsv: {len(metadata)} samples")

    env = reshape_env(samples)
    env.to_csv(out / "env_factors.tsv", sep="\t", index=False)
    print(f"  env_factors.tsv: {len(env)} samples (placeholder env)")

    tax = reshape_taxonomy(moesm5)
    tax.to_csv(out / "mag_taxonomy_labels.tsv", sep="\t",
               index=False, header=False)
    print(f"  mag_taxonomy_labels.tsv: {len(tax)} MAGs")

    target_mags = set(moesm5["Bins"])
    qual = reshape_quality(moesm4, target_mags)
    qual.to_csv(out / "quality_report.tsv", sep="\t", index=False)
    print(f"  quality_report.tsv: {len(qual)} MAGs (subset of MOESM4 1741 in MOESM5 1084)")

    abund = reshape_abundance(moesm5)
    abund.to_csv(out / "abundance.tsv", sep="\t", index=False)
    print(f"  abundance.tsv: {abund.shape[0]} MAGs × {abund.shape[1] - 1} samples (per-sample, no degeneracy)")

    ko = reshape_ko_long(moesm6)
    ko.to_csv(out / "kegg_target_only.tsv", sep="\t", index=False)
    print(f"  kegg_target_only.tsv: {len(ko)} MAG-KO records")

    print(f"\nDone. Output: {out.absolute()}")


if __name__ == "__main__":
    main()
