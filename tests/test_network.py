"""共现网络辅助 — Degree vs Betweenness 散点图 + Gephi 预处理 测试。"""
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import network
from envmeta.export.figure_export import export_to_bytes
from envmeta.tools.gephi_prep import prepare_gephi_csv, validate_gephi_format

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def net_inputs():
    nodes = pd.read_csv(SAMPLE / "gephi_nodes.csv")
    edges = pd.read_csv(SAMPLE / "gephi_edges.csv")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    return nodes, edges, tax, ks


# ── analyze: Degree vs Betweenness 散点图 ──────────────────

def test_network_smoke(net_inputs):
    nodes, edges, tax, ks = net_inputs
    r = network.analyze(nodes, edges, tax, ks)
    assert r.figure is not None
    assert not r.stats.empty
    assert {"MAG", "Degree", "Betweenness", "is_keystone"}.issubset(r.stats.columns)


def test_network_node_color_by_phylum(net_inputs):
    nodes, edges, tax, ks = net_inputs
    r = network.analyze(nodes, edges, tax, ks,
                        params={"node_color_by": "phylum"})
    assert r.figure is not None


def test_network_filter_mode_keystone_only(net_inputs):
    nodes, edges, tax, ks = net_inputs
    r = network.analyze(nodes, edges, tax, ks,
                        params={"filter_mode": "keystone_only"})
    assert r.stats["is_keystone"].all()


def test_network_export(net_inputs):
    nodes, edges, tax, ks = net_inputs
    r = network.analyze(nodes, edges, tax, ks)
    pdf = export_to_bytes(r.figure, "pdf")
    assert pdf.startswith(b"%PDF")


# ── gephi_prep: 预处理 + 校验 ──────────────────────────────

def test_validate_format_good(net_inputs):
    nodes, edges, _, _ = net_inputs
    issues = validate_gephi_format(nodes, edges)
    errors = [i for i in issues if i.startswith("[ERROR")]
    assert len(errors) == 0, f"Valid CSV should have no errors: {errors}"


def test_validate_format_missing_col():
    bad_nodes = pd.DataFrame({"X": [1, 2]})
    bad_edges = pd.DataFrame({"Y": [1]})
    issues = validate_gephi_format(bad_nodes, bad_edges)
    assert any("[ERROR]" in i for i in issues)


def test_prepare_gephi_keystone_only_labels(net_inputs):
    nodes, edges, tax, ks = net_inputs
    out_nodes, out_edges = prepare_gephi_csv(
        nodes, edges, taxonomy_df=tax, keystone_df=ks,
        label_mode="keystone_only",
    )
    # 非 keystone 的 Label 应该是空字符串
    non_ks = out_nodes[out_nodes["is_keystone"] == "False"]
    assert (non_ks["Label"] == "").all(), "非 keystone Label 应为空"
    # keystone 的 Label 不为空
    ks_rows = out_nodes[out_nodes["is_keystone"] == "True"]
    assert (ks_rows["Label"] != "").all(), "keystone Label 应有值"


def test_prepare_gephi_all_labels(net_inputs):
    nodes, edges, tax, ks = net_inputs
    out_nodes, _ = prepare_gephi_csv(
        nodes, edges, taxonomy_df=tax, keystone_df=ks,
        label_mode="all",
    )
    assert (out_nodes["Label"] != "").all()
