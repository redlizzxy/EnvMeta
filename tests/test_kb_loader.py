"""KO 知识库加载器测试。"""
import pytest

from envmeta.geocycle.knowledge_base import (
    couplings,
    element_colors,
    element_pathway_ko_order,
    element_schematic,
    flat_ko_map,
    ko_substrate_product_map,
    load_kb,
)


def test_kb_loads_4_elements_57_kos():
    kb = load_kb()
    assert set(kb["elements"].keys()) == {"arsenic", "nitrogen", "sulfur", "iron"}
    assert kb["ko_count"] == 57


def test_flat_ko_map_covers_all():
    m = flat_ko_map()
    assert len(m) == 57
    # 抽查几个
    assert m["K00537"] == ("arsC-grx", "Arsenate reduction", "arsenic")
    assert m["K00376"] == ("nosZ", "N2O reduction", "nitrogen")
    assert m["K11180"] == ("dsrA", "Dissim. sulfate red.", "sulfur")


def test_element_colors_all_hex():
    c = element_colors()
    assert len(c) == 4
    assert all(v.startswith("#") and len(v) == 7 for v in c.values())


def test_pathway_order_consistent():
    order = element_pathway_ko_order()
    assert len(order) == 57
    # 同一元素的通路应连续出现
    elements_seen = [eid for eid, _, _ in order]
    # arsenic 的所有 KO 在 nitrogen 之前
    as_idx = [i for i, e in enumerate(elements_seen) if e == "arsenic"]
    n_idx = [i for i, e in enumerate(elements_seen) if e == "nitrogen"]
    assert max(as_idx) < min(n_idx)


# --- v1.1 新增字段（S2.5-1）--------------------------------------------------


def test_ko_substrate_product_map_covers_all():
    m = ko_substrate_product_map()
    assert len(m) == 57
    # 催化型 KO 有 substrate/product
    assert m["K00537"] == {"substrate": "As(V)", "product": "As(III)"}
    assert m["K00380"] == {"substrate": "SO3-2", "product": "S-2"}
    assert m["K00376"] == {"substrate": "N2O", "product": "N2"}
    # 调控型 KO substrate/product 为 None
    assert m["K03892"] == {"substrate": None, "product": None}  # arsR
    assert m["K03711"] == {"substrate": None, "product": None}  # fur


def test_element_schematic_returns_species_and_positions():
    arsenic = element_schematic("arsenic")
    assert "As(V)" in arsenic["species"]
    assert arsenic["positions"]["As(V)"] == "left"
    assert arsenic["positions"]["As(III)"] == "right"

    sulfur = element_schematic("sulfur")
    assert "SO4-2" in sulfur["species"]
    assert sulfur["positions"]["SO4-2"] == "left"

    with pytest.raises(KeyError):
        element_schematic("nonexistent_element")


def test_couplings_returns_cross_element_reactions():
    cs = couplings()
    assert len(cs) >= 3
    types = {c["type"] for c in cs}
    assert {"precipitation", "adsorption"} <= types
    # 包含 As2S3 沉淀这条
    as2s3 = [c for c in cs if c.get("product") == "As2S3"]
    assert len(as2s3) == 1
    assert {as2s3[0]["species_a"], as2s3[0]["species_b"]} == {"As(III)", "S-2"}
    # 每条必须有颜色
    assert all(c.get("color", "").startswith("#") for c in cs)
