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


def element_pathway_ko_order(kb: dict | None = None) -> list[tuple[str, str, str]]:
    """按 element → pathway → KO 的自然顺序返回 [(element, pathway, ko), ...]。"""
    kb = kb or load_kb()
    out = []
    for eid, el in kb["elements"].items():
        for pw_id, pw in el["pathways"].items():
            for ko in pw["genes"].keys():
                out.append((eid, pw_id, ko))
    return out
