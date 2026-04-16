"""Validate mag_heatmap on real sample data (168 MAG × 10 samples)."""
from pathlib import Path

import pandas as pd

from envmeta.analysis import mag_heatmap
from envmeta.export.figure_export import export_figure

SAMPLE = Path("tests/sample_data")
OUT = Path("paper/benchmarks/validation/mag_heatmap")

ab = pd.read_csv(SAMPLE / "abundance.tsv", sep="\t")
tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                  header=None, names=["MAG", "Taxonomy"])
ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")

r = mag_heatmap.analyze(ab, tax, ks, md, params={"top_n": 30})
export_figure(r.figure, OUT / "top30_EN.pdf", "pdf")
r.stats.to_csv(OUT / "top30_stats.tsv", sep="\t", index=False)

print(f"Top {len(r.stats)} MAG extracted")
print(f"  Keystone in Top30: {int(r.stats['is_keystone'].sum())}")
phy_counts = r.stats['Phylum'].value_counts()
print(f"  Phyla: {len(phy_counts)} distinct")
print(phy_counts.head(8).to_string())
print(f"  Top-3 MAG by selection_score:")
for _, row in r.stats.sort_values('selection_score', ascending=False).head(3).iterrows():
    print(f"    {row['MAG']:<15} {row['Phylum']:<25} score={row['selection_score']:.4f}")

# 附加: variance 模式（突出组间差异 MAG）
r2 = mag_heatmap.analyze(ab, tax, ks, md,
                          params={"top_n": 30, "selection_by": "variance"})
export_figure(r2.figure, OUT / "top30_variance_EN.pdf", "pdf")
r2.stats.to_csv(OUT / "top30_variance_stats.tsv", sep="\t", index=False)
print("\n[variance mode] Top-3:")
for _, row in r2.stats.sort_values('selection_score', ascending=False).head(3).iterrows():
    print(f"    {row['MAG']:<15} {row['Phylum']:<25} score={row['selection_score']:.4f}")
