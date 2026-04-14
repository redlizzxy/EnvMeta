"""MAG 元素循环基因谱测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import gene_profile
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def gp_inputs():
    ko = pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    return ko, tax, ks


def test_gene_profile_smoke(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks)
    assert r.figure is not None
    assert not r.stats.empty
    assert {"MAG", "Phylum", "gene_count"}.issubset(r.stats.columns)


def test_gene_profile_copies_non_negative(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks)
    meta = {"MAG", "Phylum", "is_keystone", "abundance_mean", "gene_count"}
    ko_cols = [c for c in r.stats.columns if c not in meta]
    assert len(ko_cols) > 0
    vals = r.stats[ko_cols].values
    assert (vals >= 0).all()


def test_gene_profile_element_filter(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks, params={"element_filter": ["arsenic"]})
    from envmeta.geocycle.knowledge_base import flat_ko_map
    ko_map = flat_ko_map()
    meta = {"MAG", "Phylum", "is_keystone", "abundance_mean", "gene_count"}
    ko_cols = [c for c in r.stats.columns if c not in meta]
    # active_kos 都来自 arsenic 元素
    assert all(ko_map[c][2] == "arsenic" for c in ko_cols if c in ko_map)


def test_gene_profile_max_mags(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks, params={"max_mags": 20})
    assert len(r.stats) <= 20


def test_gene_profile_keystone_marked(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks)
    assert r.stats["is_keystone"].any()


def test_gene_profile_export(gp_inputs):
    ko, tax, ks = gp_inputs
    r = gene_profile.analyze(ko, tax, ks, params={"max_mags": 30})
    pdf = export_to_bytes(r.figure, "pdf")
    assert pdf.startswith(b"%PDF")
