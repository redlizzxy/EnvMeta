"""元素循环知识库加载工具。

`elements.json` 是扁平的 JSON 文件，本模块提供展开 / 索引工具供 analysis 和 geocycle
共用。
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

_DEFAULT_PATH = Path(__file__).parent / "elements.json"


@lru_cache(maxsize=4)
def load_kb(path: str | Path | None = None) -> dict:
    """加载 JSON 知识库，返回 dict（lru_cache 缓存）。"""
    p = Path(path) if path else _DEFAULT_PATH
    return json.loads(p.read_text(encoding="utf-8"))


def flat_ko_map(kb: dict | None = None) -> dict[str, tuple[str, str, str]]:
    """展开为 {ko: (gene_name, pathway, element)}。"""
    kb = kb or load_kb()
    out: dict[str, tuple[str, str, str]] = {}
    for element_id, el in kb["elements"].items():
        for pw_id, pw in el["pathways"].items():
            for ko, info in pw["genes"].items():
                out[ko] = (info["name"], pw_id, element_id)
    return out


def element_colors(kb: dict | None = None) -> dict[str, str]:
    """{element_id: hex_color}。"""
    kb = kb or load_kb()
    return {eid: el.get("color", "#888888") for eid, el in kb["elements"].items()}


def element_display(kb: dict | None = None, lang: str = "en") -> dict[str, str]:
    kb = kb or load_kb()
    return {eid: el["display_name"].get(lang, eid) for eid, el in kb["elements"].items()}


def pathway_display(kb: dict | None = None, lang: str = "en") -> dict[str, str]:
    kb = kb or load_kb()
    out: dict[str, str] = {}
    for el in kb["elements"].values():
        for pw_id, pw in el["pathways"].items():
            out[pw_id] = pw.get("display_name", {}).get(lang, pw_id)
    return out


def pathway_ko_sets(kb: dict | None = None) -> dict[str, list[str]]:
    """{pathway_id: [ko, ...]}，用于通路完整度计算。"""
    kb = kb or load_kb()
    out: dict[str, list[str]] = {}
    for el in kb["elements"].values():
        for pw_id, pw in el["pathways"].items():
            out[pw_id] = list(pw["genes"].keys())
    return out


def pathway_element_map(kb: dict | None = None) -> dict[str, str]:
    """{pathway_id: element_id}。"""
    kb = kb or load_kb()
    out: dict[str, str] = {}
    for eid, el in kb["elements"].items():
        for pw_id in el["pathways"].keys():
            out[pw_id] = eid
    return out


def element_pathway_ko_order(kb: dict | None = None) -> list[tuple[str, str, str]]:
    """按 element → pathway → KO 的自然顺序返回 [(element, pathway, ko), ...]。"""
    kb = kb or load_kb()
    out = []
    for eid, el in kb["elements"].items():
        for pw_id, pw in el["pathways"].items():
            for ko in pw["genes"].keys():
                out.append((eid, pw_id, ko))
    return out


def ko_substrate_product_map(kb: dict | None = None) -> dict[str, dict[str, str | None]]:
    """{ko: {"substrate": str|None, "product": str|None}}。

    v1.1 新增，支持细胞内级联渲染。若 KO 未声明字段则返回 {"substrate": None, "product": None}。
    """
    kb = kb or load_kb()
    out: dict[str, dict[str, str | None]] = {}
    for el in kb["elements"].values():
        for pw in el["pathways"].values():
            for ko, info in pw["genes"].items():
                out[ko] = {
                    "substrate": info.get("substrate"),
                    "product": info.get("product"),
                }
    return out


def ko_complex_map(kb: dict | None = None) -> dict[str, str | None]:
    """{ko: complex_id or None}。v2.0 (S2.5-10) 新增。

    complex_id 是 KEGG MODULE id（如 "M00595"），标记该 KO 所属的多亚基
    酶复合体。同一 complex 的连续 KO 在 renderer 里应并联显示、无内部箭头。
    """
    kb = kb or load_kb()
    out: dict[str, str | None] = {}
    for el in kb["elements"].values():
        for pw in el["pathways"].values():
            for ko, info in pw["genes"].items():
                out[ko] = info.get("complex")
    return out


def element_schematic(element_id: str, kb: dict | None = None) -> dict:
    """返回指定元素的 schematic 字段 {"species": [...], "positions": {...}}。

    若元素未声明 schematic 则返回 {"species": [], "positions": {}}。
    """
    kb = kb or load_kb()
    if element_id not in kb["elements"]:
        raise KeyError(f"unknown element_id: {element_id}")
    el = kb["elements"][element_id]
    sch = el.get("schematic") or {}
    return {
        "species": list(sch.get("species", [])),
        "positions": dict(sch.get("positions", {})),
    }


def couplings(kb: dict | None = None) -> list[dict]:
    """返回顶级 couplings 列表；缺失时返回 []。

    每条结构：{species_a, species_b, product, type, color, description?}。
    """
    kb = kb or load_kb()
    return list(kb.get("couplings", []))
