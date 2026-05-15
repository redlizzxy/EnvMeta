"""
Grettenberger & Hamilton 2021 AEM (DOI 10.1128/AEM.00772-21)
Table 1 (HTML) + Data Set S1 KEGGModuleStepHit (xlsx) → EnvMeta input files

Source license: CC BY 4.0
Path: plug-and-play (METABOLIC step-level data already published)

用法:
    python tools/external_benchmarks/grettenberger2021_reshape.py \\
        --table1 D:\\download\\grettenberger_table1.html \\
        --dataset_s1 D:\\download\\aem.00772-21-s0002.xlsx \\
        --out paper/benchmarks/external/grettenberger_2021_amd_stream/input_data_local

输出 6 个文件:
    metadata.tsv             # 1 pooled sample (Grettenberger 论文是 5 sites
                              # combined, but Data Set S1 仅给 MAG-level 数据，
                              # 无 sample × MAG abundance, 故用单组 pooled)
    env_factors.tsv          # placeholder (AMD 一般 pH 2-4)
    mag_taxonomy_labels.tsv  # 29 MAGs (Data Set S1 子集) × GTDB-Tk taxonomy
    quality_report.tsv       # 29 MAGs × completeness/contamination/size/GC
    abundance.tsv            # 29 MAGs × 1 sample (placeholder=1)
    kegg_target_only.tsv     # MAG × KO long format (METABOLIC step Present 展开)

注:
- KO 展开是 conservative inflation：1 step (K1,K2,K3) Present → 3 条 KO
  records。EnvMeta pathway_completeness 因此偏高，但对 hypothesis scoring 的
  影响是可能 inflate STRONG。INSUFFICIENT 仍可信，STRONG 需在 README 标注。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

TAX_PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__"]


def fmt_taxonomy(s: str | float) -> str:
    """'Bacteria; Actinobacteriota; ...' → 'd__Bacteria;p__Actinobacteriota;...;s__'"""
    if pd.isna(s):
        return ";".join(TAX_PREFIXES) + ";s__"
    parts = [p.strip() for p in str(s).split(";")]
    out = []
    for i, pre in enumerate(TAX_PREFIXES):
        if i < len(parts) and parts[i]:
            out.append(pre + parts[i])
        else:
            out.append(pre)
    return ";".join(out) + ";s__"


def load_table1(html_path: Path) -> pd.DataFrame:
    t = pd.read_html(html_path)[0]
    t["MAG_num"] = t["MAG"].astype(str).str.extract(r"MAG\s*(\d+)").astype(int)
    return t


def load_dataset_s1(xlsx_path: Path) -> tuple[list, list, pd.DataFrame]:
    """Returns (step_ids, ko_lists, mag_data_df with MAG_num/MAG_full + step columns)"""
    df = pd.read_excel(xlsx_path, sheet_name="KEGGModuleStepHit", header=None)
    step_ids = df.iloc[0, 2:].tolist()
    ko_lists = df.iloc[2, 2:].tolist()

    # Row 4+ = MAG data; row 3 = "Module Category" (skip)
    mag_data = df.iloc[4:, :].copy()
    mag_data.columns = ["MAG_label", "MAG_full"] + list(step_ids)
    mag_data = mag_data[mag_data["MAG_label"].notna()
                       & mag_data["MAG_label"].astype(str).str.startswith("MAG")].copy()
    mag_data["MAG_num"] = (
        mag_data["MAG_label"].astype(str).str.extract(r"MAG\s*(\d+)").astype(int)
    )
    return step_ids, ko_lists, mag_data


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--table1", required=True, help="Path to grettenberger_table1.html")
    p.add_argument("--dataset_s1", required=True, help="Path to aem.00772-21-s0002.xlsx")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading Table 1 HTML...")
    t1 = load_table1(Path(args.table1))
    print(f"  Table 1: {len(t1)} MAGs")

    print("Loading Data Set S1 KEGGModuleStepHit...")
    step_ids, ko_lists, mag_data = load_dataset_s1(Path(args.dataset_s1))
    print(f"  Data Set S1: {len(mag_data)} MAGs × {len(step_ids)} steps")

    # Inner-join on MAG_num
    merged = t1.merge(mag_data, on="MAG_num", how="inner")
    n = len(merged)
    print(f"\nJoined: {n} MAGs (intersection)")

    # ── 1. metadata.tsv ──
    samples = ["All"]
    pd.DataFrame({
        "SampleID": samples, "Group": samples, "Replicate": [1],
    }).to_csv(out / "metadata.tsv", sep="\t", index=False)
    print(f"  metadata.tsv: {len(samples)} sample (pooled)")

    # ── 2. env_factors.tsv ──
    pd.DataFrame({
        "SampleID": samples, "Group": samples, "pH_proxy": [3.0],
    }).to_csv(out / "env_factors.tsv", sep="\t", index=False)
    print(f"  env_factors.tsv: placeholder (pH=3.0 proxy)")

    # ── 3. mag_taxonomy_labels.tsv ──
    tax = pd.DataFrame({
        "MAG_ID": merged["MAG_full"],
        "taxonomy": merged["Taxonomy"].apply(fmt_taxonomy),
    })
    tax.to_csv(out / "mag_taxonomy_labels.tsv", sep="\t", index=False, header=False)
    print(f"  mag_taxonomy_labels.tsv: {len(tax)} MAGs")

    # ── 4. quality_report.tsv ──
    qual = pd.DataFrame({
        "Name": merged["MAG_full"],
        "Completeness": pd.to_numeric(merged["Completeness (%)"], errors="coerce"),
        "Contamination": pd.to_numeric(merged["Contamination (%)"], errors="coerce"),
        "Completeness_Model_Used": ["Grettenberger2021_Table1"] * n,
        "Translation_Table_Used": [11] * n,
        "Coding_Density": [0.92] * n,
        "Contig_N50": [0] * n,
        "Average_Gene_Length": [0] * n,
        "Genome_Size": (pd.to_numeric(merged["Size (Mbp)"], errors="coerce") * 1_000_000).fillna(0),
        "GC_Content": pd.to_numeric(merged["GC content (%)"], errors="coerce").fillna(0) / 100,
        "Total_Coding_Sequences": pd.to_numeric(
            merged["No. of protein coding sequences"], errors="coerce"
        ).fillna(0),
        "Total_Contigs": pd.to_numeric(merged["No. of contigs"], errors="coerce").fillna(0),
        "Max_Contig_Length": [0] * n,
        "Additional_Notes": ["Grettenberger 2021 Table 1 + Data Set S1"] * n,
    })
    qual.to_csv(out / "quality_report.tsv", sep="\t", index=False)
    print(f"  quality_report.tsv: {n} MAGs")

    # ── 5. abundance.tsv ──
    pd.DataFrame({
        "Genome": merged["MAG_full"],
        "All": [1.0] * n,
    }).to_csv(out / "abundance.tsv", sep="\t", index=False)
    print(f"  abundance.tsv: {n} MAGs × 1 placeholder sample")

    # ── 6. kegg_target_only.tsv ──
    # Expand each MAG × step (Present) → all KO ids
    ko_records = []
    n_steps_present = 0
    n_kos_total = 0
    for _, row in merged.iterrows():
        mag = row["MAG_full"]
        for step_id, ko_str in zip(step_ids, ko_lists):
            present = row.get(step_id, "Absent")
            if present != "Present":
                continue
            n_steps_present += 1
            kos = re.findall(r"K\d{5}", str(ko_str)) if pd.notna(ko_str) else []
            for ko in set(kos):
                ko_records.append({
                    "MAG": mag,
                    "Gene_ID": f"{mag}__{step_id}__{ko}",
                    "KEGG_ko": f"ko:{ko}",
                    "Description": f"METABOLIC step {step_id}",
                })
                n_kos_total += 1
    ko_df = pd.DataFrame(ko_records).drop_duplicates(subset=["MAG", "KEGG_ko"])
    ko_df.to_csv(out / "kegg_target_only.tsv", sep="\t", index=False)
    print(f"  kegg_target_only.tsv: {len(ko_df)} unique MAG-KO records "
          f"(from {n_steps_present} steps × {n_kos_total} KO assignments, dedup)")

    print(f"\nDone. Output: {out.absolute()}")
    print(f"\nNote: KO assignment is conservative inflation (each Present step "
          f"expands to all listed KOs). EnvMeta pathway_completeness may "
          f"slightly inflate. INSUFFICIENT label is robust; STRONG label "
          f"should be reported with this caveat.")


if __name__ == "__main__":
    main()
