"""Reshape Liu 2023 cold-seep data to fill gaps for performance benchmark.

Generates:
- alpha.tsv         (Shannon / Simpson / Chao1 / Pielou per sample)
- distance_bray.tsv (Bray-Curtis distance matrix from MAG abundance)
- ko_wide.tsv       (KO x sample wide; MAG-level KO presence x MAG abundance)
- phylum_abundance.tsv (Phylum x sample, aggregated from MAG abundance + taxonomy)

Inputs: paper/benchmarks/external/liu_2023_coldseep/input_data_local/
Outputs: paper/benchmarks/performance/liu_reshaped/
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

LIU_IN = Path("paper/benchmarks/external/liu_2023_coldseep/input_data_local")
OUT = Path("paper/benchmarks/performance/liu_reshaped")
OUT.mkdir(parents=True, exist_ok=True)


def compute_alpha(abundance: pd.DataFrame) -> pd.DataFrame:
    """Per-sample alpha diversity from MAG abundance matrix (Genome x sample)."""
    samples = [c for c in abundance.columns if c != "Genome"]
    mat = abundance[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy()
    rows = []
    for j, s in enumerate(samples):
        col = mat[:, j]
        total = col.sum()
        if total <= 0:
            rows.append({"SampleID": s, "Shannon": 0.0, "Simpson": 0.0,
                         "Chao1": 0.0, "Pielou": 0.0})
            continue
        p = col / total
        nz = p[p > 0]
        shannon = -np.sum(nz * np.log(nz))
        simpson = 1.0 - np.sum(nz ** 2)
        observed = (col > 0).sum()
        f1 = (col == 1).sum()
        f2 = (col == 2).sum()
        chao1 = observed + (f1 * (f1 - 1)) / (2 * (f2 + 1)) if observed > 0 else 0.0
        pielou = shannon / np.log(observed) if observed > 1 else 0.0
        rows.append({"SampleID": s, "Shannon": shannon, "Simpson": simpson,
                     "Chao1": chao1, "Pielou": pielou})
    return pd.DataFrame(rows)


def compute_bray_curtis(abundance: pd.DataFrame) -> pd.DataFrame:
    """Bray-Curtis distance matrix from sample x feature matrix.

    Returns: DataFrame with first column 'SampleID', remaining columns are
    sample IDs forming the square distance matrix.
    """
    samples = [c for c in abundance.columns if c != "Genome"]
    mat = abundance[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy().T  # sample x feature
    n = len(samples)
    d = np.zeros((n, n))
    sums = mat.sum(axis=1)
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(mat[i] - mat[j]).sum()
            den = sums[i] + sums[j]
            d[i, j] = d[j, i] = num / den if den > 0 else 0.0
    df = pd.DataFrame(d, index=samples, columns=samples)
    df.index.name = "SampleID"
    return df.reset_index()


def pivot_ko_wide(ko_long: pd.DataFrame, abundance: pd.DataFrame) -> pd.DataFrame:
    """KO long (MAG x KO presence) -> wide (KO x sample, sum over MAG abundance).

    Per sample, KO abundance = sum_{MAG: MAG has KO} abundance[MAG, sample].
    """
    ko_long = ko_long.copy()
    ko_long["KEGG_ko"] = ko_long["KEGG_ko"].astype(str).str.replace("ko:", "", regex=False).str.strip()
    mag_ko = ko_long.groupby("MAG")["KEGG_ko"].apply(set).to_dict()
    samples = [c for c in abundance.columns if c != "Genome"]
    abundance_idx = abundance.set_index("Genome")[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    all_kos = sorted({k for kos in mag_ko.values() for k in kos if k.startswith("K")})
    out = pd.DataFrame(0.0, index=all_kos, columns=samples)
    for mag, kos in mag_ko.items():
        if mag not in abundance_idx.index:
            continue
        valid_kos = [k for k in kos if k in out.index]
        if not valid_kos:
            continue
        ab = abundance_idx.loc[mag].to_numpy()
        out.loc[valid_kos] += ab
    out.index.name = "KO"
    return out.reset_index()


def aggregate_phylum(abundance: pd.DataFrame, taxonomy: pd.DataFrame) -> pd.DataFrame:
    """MAG abundance + taxonomy -> Phylum x sample matrix (long taxonomy semicolon-separated)."""
    tax = taxonomy.copy()
    tax.columns = ["Genome", "Taxonomy"]
    def parse_phylum(s: str) -> str:
        for part in str(s).split(";"):
            part = part.strip()
            if part.startswith("p__"):
                p = part[3:].strip()
                return p if p else "Unclassified"
        return "Unclassified"
    tax["Phylum"] = tax["Taxonomy"].apply(parse_phylum)
    merged = abundance.merge(tax[["Genome", "Phylum"]], on="Genome", how="left")
    merged["Phylum"] = merged["Phylum"].fillna("Unclassified")
    samples = [c for c in abundance.columns if c != "Genome"]
    agg = merged.groupby("Phylum")[samples].sum().reset_index()
    agg = agg.rename(columns={"Phylum": "Taxonomy"})
    return agg


def main() -> int:
    print(f"Loading Liu 2023 data from {LIU_IN}...")
    abundance = pd.read_csv(LIU_IN / "abundance.tsv", sep="\t")
    ko_long = pd.read_csv(LIU_IN / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(LIU_IN / "mag_taxonomy_labels.tsv", sep="\t", header=None)
    n_mags, n_samples = abundance.shape[0], abundance.shape[1] - 1
    print(f"  abundance: {n_mags} MAGs x {n_samples} samples")
    print(f"  ko_long:   {len(ko_long)} MAG-KO records")
    print(f"  taxonomy:  {len(tax)} MAGs")

    print("[1/4] Computing alpha diversity...")
    alpha = compute_alpha(abundance)
    alpha.to_csv(OUT / "alpha.tsv", sep="\t", index=False)
    print(f"  -> {OUT / 'alpha.tsv'}  ({len(alpha)} samples)")

    print("[2/4] Computing Bray-Curtis distance matrix...")
    dist = compute_bray_curtis(abundance)
    dist.to_csv(OUT / "distance_bray.tsv", sep="\t", index=False)
    print(f"  -> {OUT / 'distance_bray.tsv'}  ({len(dist)} x {len(dist)})")

    print("[3/4] Pivoting KO long -> wide...")
    ko_wide = pivot_ko_wide(ko_long, abundance)
    ko_wide.to_csv(OUT / "ko_wide.tsv", sep="\t", index=False)
    print(f"  -> {OUT / 'ko_wide.tsv'}  ({len(ko_wide)} KOs x {ko_wide.shape[1]-1} samples)")

    print("[4/4] Aggregating MAG abundance by Phylum...")
    phy = aggregate_phylum(abundance, tax)
    phy.to_csv(OUT / "phylum_abundance.tsv", sep="\t", index=False)
    print(f"  -> {OUT / 'phylum_abundance.tsv'}  ({len(phy)} phyla x {phy.shape[1]-1} samples)")

    print("\nDone. Reshaped outputs in:", OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())
