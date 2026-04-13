"""KO 知识库加载器测试。"""
from envmeta.geocycle.knowledge_base import (
    element_colors,
    element_pathway_ko_order,
    flat_ko_map,
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
