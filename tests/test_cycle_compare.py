"""跨组循环图对比（S2.5-7d）。"""
from __future__ import annotations

from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis.cycle_compare import compare_groups

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def inputs():
    ko = pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    ab = pd.read_csv(SAMPLE / "abundance.tsv", sep="\t")
    env = pd.read_csv(SAMPLE / "env_factors.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ko, tax, ks, ab, env, md


def test_compare_groups_returns_long_table(inputs):
    ko, tax, ks, ab, env, md = inputs
    df = compare_groups(ko, tax, ks, ab, env, md,
                        params={"perm_n": 99})
    assert not df.empty
    required = {
        "group", "element", "pathway_id", "display_name",
        "n_active_mags", "mean_completeness", "total_contribution",
        "top_mag", "top_mag_genus", "top_mag_is_keystone",
    }
    assert required <= set(df.columns)


def test_compare_groups_covers_all_metadata_groups(inputs):
    ko, tax, ks, ab, env, md = inputs
    df = compare_groups(ko, tax, ks, ab, env, md,
                        params={"perm_n": 99})
    expected = set(md["Group"].astype(str))
    assert set(df["group"]) == expected


def test_compare_groups_keystone_flag_is_bool(inputs):
    ko, tax, ks, ab, env, md = inputs
    df = compare_groups(ko, tax, ks, ab, env, md,
                        params={"perm_n": 99})
    assert df["top_mag_is_keystone"].dtype == bool or set(
        df["top_mag_is_keystone"].unique()) <= {True, False}
