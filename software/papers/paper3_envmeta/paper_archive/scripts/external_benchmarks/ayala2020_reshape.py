"""
Ayala-Muñoz et al. 2020 Microorganisms (DOI 10.3390/microorganisms8091350)
GhostKOALA user_ko.txt → EnvMeta 6 输入文件

Source license: CC BY 4.0
Path: GhostKOALA 重跑 (2020 论文未公开 KEGG 注释，用 NCBI assembly 重做 ORF + GhostKOALA)

13 MAGs from BioProject PRJNA646106:
  Archaea (5): A_CRE_07, A_EUR_01, A_EUR_06, A_MIC_10, A_NAN_12
  Bacteria (8): B_ACI_09, B_ACT_02, B_ACT_11, B_CHL_03, B_DOR_08,
                 B_NIT_04, B_PAT_13, B_PRO_05

用法:
    python tools/external_benchmarks/ayala2020_reshape.py \\
        --user_ko paper/benchmarks/external/ayala_2020_pitlake/input_data_local/user_ko.txt \\
        --out paper/benchmarks/external/ayala_2020_pitlake/input_data_local

输出 6 个文件:
    metadata.tsv             # 1 pooled sample "All" (Ayala 论文 deep layer 35m，单组)
    env_factors.tsv          # placeholder (pit lake pH 2.5)
    mag_taxonomy_labels.tsv  # 13 MAGs × phylum-level taxonomy (从文件名推)
    quality_report.tsv       # 13 MAGs × placeholder quality (待 Suppl Table 替换)
    abundance.tsv            # 13 MAGs × 1 sample (placeholder=1)
    kegg_target_only.tsv     # MAG × KO long format (从 GhostKOALA user_ko.txt)

注:
- Taxonomy 从 MAG 文件名 phylum 前缀推断（粗 phylum 级，CRE/EUR/MIC/NAN/
  ACI/ACT/CHL/DOR/NIT/PAT/PRO）。完整 taxonomy 待 Ayala 2020 Suppl Table
  S1 提供后用论文实测 GTDB-Tk 替换。
- Quality placeholder（completeness=70, contamination=5）—— hypothesis
  scoring 不依赖该字段，cycle_diagram 也不严格依赖；mag_quality 图会跑
  出 placeholder 值。
- Abundance placeholder=1（13 MAGs × 1 pooled sample 时唯一可解）。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# MAG 文件名前缀 → GTDB phylum (粗推断)
PREFIX_TO_TAXONOMY = {
    "A_CRE": ("d__Archaea", "p__Crenarchaeota"),
    "A_EUR": ("d__Archaea", "p__Euryarchaeota"),
    "A_MIC": ("d__Archaea", "p__Micrarchaeota"),
    "A_NAN": ("d__Archaea", "p__Nanoarchaeota"),
    "B_ACI": ("d__Bacteria", "p__Acidobacteriota"),
    "B_ACT": ("d__Bacteria", "p__Actinobacteriota"),
    "B_CHL": ("d__Bacteria", "p__Chloroflexota"),
    "B_DOR": ("d__Bacteria", "p__Dormibacterota"),
    "B_NIT": ("d__Bacteria", "p__Nitrospirae"),
    "B_PAT": ("d__Bacteria", "p__Patescibacteria"),
    "B_PRO": ("d__Bacteria", "p__Proteobacteria"),
}

TAX_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__"]


def fmt_taxonomy(domain: str, phylum: str) -> str:
    """13 MAGs 仅有 domain + phylum 信息，fill 其余为空。"""
    out = [domain, phylum]
    for pre in TAX_PREFIXES[2:]:  # c__, o__, f__, g__
        out.append(pre)
    return ";".join(out) + ";s__"


def parse_user_ko(path: Path) -> pd.DataFrame:
    """GhostKOALA user_ko.txt → MAG-KO long table.

    Format:
        protein_id <TAB> KO_id  (KO_id 缺失则该列空)
        protein_id 形如 'B_ACI_09|JACWAK010000186.1|gene_2'

    Returns DataFrame with columns: MAG, KEGG_ko, Gene_ID, Description
    """
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            protein_id = parts[0]
            ko = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
            if ko is None:
                continue  # 跳过未注释的
            mag = protein_id.split("|", 1)[0]
            rows.append({
                "MAG": mag,
                "KEGG_ko": f"ko:{ko}",
                "Gene_ID": protein_id,
                "Description": f"GhostKOALA {ko}",
            })
    return pd.DataFrame(rows).drop_duplicates(subset=["MAG", "KEGG_ko", "Gene_ID"])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--user_ko", required=True, help="Path to GhostKOALA user_ko.txt")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    user_ko_path = Path(args.user_ko)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print(f"Loading GhostKOALA output: {user_ko_path}")
    ko_df = parse_user_ko(user_ko_path)
    mags = sorted(ko_df["MAG"].unique())
    print(f"  {len(ko_df)} unique MAG-KO records across {len(mags)} MAGs")
    print(f"  MAGs: {mags}")
    n = len(mags)

    # ── 1. metadata.tsv ──
    samples = ["All"]
    pd.DataFrame({
        "SampleID": samples, "Group": samples, "Replicate": [1],
    }).to_csv(out / "metadata.tsv", sep="\t", index=False)
    print(f"  metadata.tsv: {len(samples)} sample (pooled)")

    # ── 2. env_factors.tsv ──
    pd.DataFrame({
        "SampleID": samples, "Group": samples, "pH_proxy": [2.5],
    }).to_csv(out / "env_factors.tsv", sep="\t", index=False)
    print("  env_factors.tsv: placeholder (pH=2.5 IPB pit lake)")

    # ── 3. mag_taxonomy_labels.tsv ──
    tax_rows = []
    for mag in mags:
        prefix = "_".join(mag.split("_")[:2])  # B_ACI_09 → B_ACI
        domain, phylum = PREFIX_TO_TAXONOMY.get(prefix, ("d__Bacteria", "p__"))
        tax_rows.append({"MAG_ID": mag, "taxonomy": fmt_taxonomy(domain, phylum)})
    pd.DataFrame(tax_rows).to_csv(
        out / "mag_taxonomy_labels.tsv",
        sep="\t", index=False, header=False,
    )
    print(f"  mag_taxonomy_labels.tsv: {n} MAGs (phylum-level)")

    # ── 4. quality_report.tsv ──
    pd.DataFrame({
        "Name": mags,
        "Completeness": [70.0] * n,
        "Contamination": [5.0] * n,
        "Completeness_Model_Used": ["Ayala2020_placeholder"] * n,
        "Translation_Table_Used": [11] * n,
        "Coding_Density": [0.92] * n,
        "Contig_N50": [0] * n,
        "Average_Gene_Length": [0] * n,
        "Genome_Size": [2_500_000] * n,
        "GC_Content": [0.5] * n,
        "Total_Coding_Sequences": [0] * n,
        "Total_Contigs": [0] * n,
        "Max_Contig_Length": [0] * n,
        "Additional_Notes": [
            "Ayala 2020 placeholder (待 Suppl Table 替换)"
        ] * n,
    }).to_csv(out / "quality_report.tsv", sep="\t", index=False)
    print(f"  quality_report.tsv: {n} MAGs (placeholder quality)")

    # ── 5. abundance.tsv ──
    pd.DataFrame({
        "Genome": mags,
        "All": [1.0] * n,
    }).to_csv(out / "abundance.tsv", sep="\t", index=False)
    print(f"  abundance.tsv: {n} MAGs × 1 placeholder sample")

    # ── 6. kegg_target_only.tsv ──
    ko_df[["MAG", "Gene_ID", "KEGG_ko", "Description"]].to_csv(
        out / "kegg_target_only.tsv", sep="\t", index=False,
    )
    print(f"  kegg_target_only.tsv: {len(ko_df)} unique MAG-KO records")

    print(f"\nDone. Output: {out.absolute()}")


if __name__ == "__main__":
    main()
