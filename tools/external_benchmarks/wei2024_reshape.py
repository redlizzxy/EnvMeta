"""
Wei et al. 2024 Microbiome (DOI 10.1186/s40168-024-01952-4) MOESM2_ESM.xlsx
→ EnvMeta input files

Source license: CC BY-NC-ND 4.0 — 原 supplementary 数据不入 git repo（NoDerivatives
限制），脚本仅本地处理，输出落到 paper/benchmarks/external/wei_2024_paddy/
input_data_local/（已加 .gitignore）。

用法:
    python tools/external_benchmarks/wei2024_reshape.py \\
        --xlsx "D:\\download\\40168_2024_1952_MOESM2_ESM (1).xlsx" \\
        --out paper/benchmarks/external/wei_2024_paddy/input_data_local

输出 6 个文件:
    metadata.tsv             # SampleID / Group / Replicate
    env_factors.tsv          # SampleID / Group / pH / 14 化学元素
    mag_taxonomy_labels.tsv  # MAG_ID \\t taxonomy 字符串 (无 header)
    quality_report.tsv       # CheckM2-style schema
    abundance.tsv            # Genome × 36 sample (group-mean 复制到 per-sample)
    kegg_target_only.tsv     # MAG / Gene_ID / KEGG_ko / Description (long)

注:
    abundance.tsv 维度退化（同组 sample 填同 group mean），PCoA/RDA 不可用，
    但 group-aware 图（heatmap / 循环图 / 假说评分）可用。
    arxA / nrfA 跳过（EnvMeta KB v1.1 未收录）。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Wei 14 个目标基因 → EnvMeta KB v1.1 KO 编号
# None 表示 KB 未收录（脚本跳过该 gene 的 KO 注入，README 里说明）
GENE_TO_KO = {
    "aioA":  "K08356",  # aoxB / aioA large subunit
    "arxA":  None,       # ROCker 自定义模型，KEGG 无标准 KO，KB v1.1 未收录
    "arrA":  "K28466",
    "arsC1": "K00537",
    "arsC2": "K03741",
    "napA":  "K02567",
    "narG":  "K00370",
    "nirK":  "K00368",
    "nirS":  "K15864",
    "norB":  "K04561",
    "nosZ":  "K00376",
    "nrfA":  None,       # DNRA pathway，KB v1.1 未收录（论文 Discussion 需说明扩展计划）
    "pmoA":  "K10944",
    "pmoB":  "K10945",
}

GENE_DESCRIPTIONS = {
    "aioA":  "arsenite oxidase large subunit (aoxB family)",
    "arxA":  "anaerobic arsenite oxidase (ROCker custom model)",
    "arrA":  "respiratory arsenate reductase",
    "arsC1": "arsenate reductase (glutaredoxin coupled)",
    "arsC2": "arsenate reductase (thioredoxin coupled)",
    "napA":  "periplasmic nitrate reductase",
    "narG":  "membrane-bound nitrate reductase alpha",
    "nirK":  "Cu-containing nitrite reductase",
    "nirS":  "cd1-type nitrite reductase",
    "norB":  "nitric oxide reductase B",
    "nosZ":  "nitrous oxide reductase",
    "nrfA":  "DNRA periplasmic cytochrome c nitrite reductase",
    "pmoA":  "particulate methane monooxygenase A",
    "pmoB":  "particulate methane monooxygenase B",
}

ENV_COLS_FLAT = [
    "SampleID", "pH", "WH2O_pct",
    "Cr", "Sb", "As", "Cd", "Pb", "Cu", "NH4_N", "NO3_N", "TotalFe",
    "OM", "TC", "TN", "TP", "TS",
]

S7_COLS = [
    "Number", "MAG_ID", "Completeness", "Contamination",
    "AsContam_abundance", "NoContam_abundance",
    "Domain", "Phylum", "Class", "Order", "Family", "Genus",
]

S8_COLS = [
    "MAG_ID", "Taxa",
    "aioA", "arxA", "arrA", "arsC1", "arsC2",
    "napA", "narG", "nirK", "nirS", "norB", "nosZ", "nrfA",
    "pmoA", "pmoB",
]


def reshape_metadata(xlsx_path: Path) -> pd.DataFrame:
    """Table S1 → metadata: 仅取 36 个 metagenome 样本（filter 掉 metatranscriptome）"""
    df = pd.read_excel(
        xlsx_path, sheet_name="Table S1. Sample Information", header=1
    )
    df = df[df["Sample ID"].notna()].copy()
    df = df[df["Type"] == "soil metagenome"].copy()
    df["Group"] = df["Arsenic level"]
    df["Replicate"] = df.groupby("Location (Province, China)").cumcount() + 1
    return df[["Sample ID", "Group", "Replicate"]].rename(
        columns={"Sample ID": "SampleID"}
    )


def reshape_env(xlsx_path: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    """Table S2 → env_factors: 仅取 row0-35（36 个 individual sample），跳过 site-level mean 块"""
    df = pd.read_excel(
        xlsx_path, sheet_name="Table S2. Soil Properties",
        header=None, skiprows=1,
    )
    extra_pad = max(len(df.columns) - len(ENV_COLS_FLAT), 0)
    df.columns = ENV_COLS_FLAT + [f"_extra_{i}" for i in range(extra_pad)]
    df = df.iloc[2:].reset_index(drop=True)
    valid_ids = set(metadata["SampleID"])
    df = df[df["SampleID"].isin(valid_ids)].copy()
    df = df[ENV_COLS_FLAT]
    for c in ENV_COLS_FLAT[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.merge(metadata[["SampleID", "Group"]], on="SampleID", how="left")
    return df[["SampleID", "Group"] + ENV_COLS_FLAT[1:]]


def reshape_mag_basic(xlsx_path: Path) -> pd.DataFrame:
    df = pd.read_excel(
        xlsx_path, sheet_name="Table S7. MAGs Basic Info",
        header=None, skiprows=1,
    )
    df.columns = S7_COLS
    df = df.iloc[2:].reset_index(drop=True)
    df = df[df["MAG_ID"].notna()].copy()
    for c in [
        "Completeness", "Contamination",
        "AsContam_abundance", "NoContam_abundance",
    ]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def build_taxonomy_labels(mag_basic: pd.DataFrame) -> pd.DataFrame:
    prefixes = ["d__", "p__", "c__", "o__", "f__", "g__"]
    levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus"]

    def fmt(row):
        parts = []
        for level, pre in zip(levels, prefixes):
            v = row[level]
            if pd.isna(v) or v == "\\":
                v = ""
            if v and not v.startswith(pre):
                v = pre + v
            parts.append(v)
        return ";".join(parts) + ";s__"

    return pd.DataFrame({
        "MAG_ID": mag_basic["MAG_ID"],
        "taxonomy": mag_basic.apply(fmt, axis=1),
    })


def build_quality_report(mag_basic: pd.DataFrame) -> pd.DataFrame:
    n = len(mag_basic)
    return pd.DataFrame({
        "Name": mag_basic["MAG_ID"],
        "Completeness": mag_basic["Completeness"],
        "Contamination": mag_basic["Contamination"],
        "Completeness_Model_Used": ["Wei2024_published"] * n,
        "Translation_Table_Used": [11] * n,
        "Coding_Density": [0.92] * n,
        "Contig_N50": [0] * n,
        "Average_Gene_Length": [0] * n,
        "Genome_Size": [0] * n,
        "GC_Content": [0.0] * n,
        "Total_Coding_Sequences": [0] * n,
        "Total_Contigs": [0] * n,
        "Max_Contig_Length": [0] * n,
        "Additional_Notes": ["Wei2024 MOESM2 Table S7"] * n,
    })


def build_abundance(mag_basic: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    samples = metadata["SampleID"].tolist()
    sample_to_group = dict(zip(metadata["SampleID"], metadata["Group"]))
    rows = []
    for _, mag in mag_basic.iterrows():
        row = {"Genome": mag["MAG_ID"]}
        for s in samples:
            g = sample_to_group[s]
            if g == "AsContam":
                row[s] = mag["AsContam_abundance"]
            elif g == "NoContam":
                row[s] = mag["NoContam_abundance"]
            else:
                row[s] = 0.0
        rows.append(row)
    return pd.DataFrame(rows)


def reshape_ko_long(xlsx_path: Path) -> pd.DataFrame:
    df = pd.read_excel(
        xlsx_path, sheet_name="Table S8.MAGs With Target Gene",
        header=None, skiprows=1,
    )
    df.columns = S8_COLS
    df = df.iloc[2:].reset_index(drop=True)
    df = df[df["MAG_ID"].notna()].copy()

    skipped = []
    long_rows = []
    for _, row in df.iterrows():
        mag = row["MAG_ID"]
        for gene, ko in GENE_TO_KO.items():
            cnt = row.get(gene, 0)
            try:
                cnt = int(cnt) if pd.notna(cnt) else 0
            except (ValueError, TypeError):
                cnt = 0
            if ko is None:
                if cnt > 0:
                    skipped.append((mag, gene, cnt))
                continue
            for i in range(cnt):
                long_rows.append({
                    "MAG": mag,
                    "Gene_ID": f"{mag}__{gene}_{i + 1}",
                    "KEGG_ko": f"ko:{ko}",
                    "Description": GENE_DESCRIPTIONS[gene],
                })
    if skipped:
        print(f"  [Skipped (KB v1.1 unmapped)] {len(skipped)} records "
              f"across {len({s[0] for s in skipped})} MAGs:")
        for gene in sorted({s[1] for s in skipped}):
            n = sum(1 for s in skipped if s[1] == gene)
            tot_cnt = sum(s[2] for s in skipped if s[1] == gene)
            print(f"    {gene}: {n} MAGs / {tot_cnt} gene copies")
    return pd.DataFrame(long_rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--xlsx", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    xlsx = Path(args.xlsx)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {xlsx}")

    metadata = reshape_metadata(xlsx)
    metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)
    print(f"  metadata.tsv: {len(metadata)} samples")

    env = reshape_env(xlsx, metadata)
    env.to_csv(out / "env_factors.tsv", sep="\t", index=False)
    print(f"  env_factors.tsv: {len(env)} samples × {len(env.columns)} factors")

    mag_basic = reshape_mag_basic(xlsx)
    print(f"  [Table S7 read] {len(mag_basic)} MAGs")

    tax = build_taxonomy_labels(mag_basic)
    tax.to_csv(out / "mag_taxonomy_labels.tsv", sep="\t",
               index=False, header=False)
    print(f"  mag_taxonomy_labels.tsv: {len(tax)} MAGs")

    qual = build_quality_report(mag_basic)
    qual.to_csv(out / "quality_report.tsv", sep="\t", index=False)
    print(f"  quality_report.tsv: {len(qual)} MAGs")

    abund = build_abundance(mag_basic, metadata)
    abund.to_csv(out / "abundance.tsv", sep="\t", index=False)
    print(f"  abundance.tsv: {len(abund)} MAGs × {len(abund.columns) - 1} samples")

    ko_long = reshape_ko_long(xlsx)
    ko_long.to_csv(out / "kegg_target_only.tsv", sep="\t", index=False)
    print(f"  kegg_target_only.tsv: {len(ko_long)} MAG-KO records")

    print(f"\nDone. Output: {out.absolute()}")


if __name__ == "__main__":
    main()
