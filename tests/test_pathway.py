"""MAG 通路完整度测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import pathway
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def pathway_inputs():
    ko = pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    return ko, tax, ks


def test_pathway_smoke(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks)
    assert r.figure is not None
    assert not r.stats.empty
    assert {"MAG", "Phylum", "total_completeness"}.issubset(r.stats.columns)


def test_pathway_completeness_range(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks)
    meta_cols = {"MAG", "label", "Phylum", "Genus", "Species",
                 "is_keystone", "abundance_mean", "total_completeness"}
    pw_cols = [c for c in r.stats.columns if c not in meta_cols]
    assert len(pw_cols) >= 10  # 至少 10 条通路
    vals = r.stats[pw_cols].values
    assert (vals >= 0).all() and (vals <= 100).all()


def test_pathway_max_mags_top(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks, params={"max_mags": 20})
    assert len(r.stats) <= 20


def test_pathway_element_filter(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks, params={"element_filter": ["arsenic"]})
    # 剩余通路列应全部属 arsenic
    from envmeta.geocycle.knowledge_base import pathway_element_map
    pw_elem = pathway_element_map()
    pw_cols = [c for c in r.stats.columns if c in pw_elem]
    assert all(pw_elem[c] == "arsenic" for c in pw_cols)


def test_pathway_bubble_style(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks, params={"style": "bubble", "max_mags": 10})
    assert r.figure is not None


def test_pathway_export(pathway_inputs):
    ko, tax, ks = pathway_inputs
    r = pathway.analyze(ko, tax, ks, params={"max_mags": 30})
    pdf = export_to_bytes(r.figure, "pdf")
    assert pdf.startswith(b"%PDF")
