"""KB 骨架构建器（S2.5-10b）。

从 `envmeta/geocycle/kegg_snapshot.json`（S2.5-10a 生成）+
`scripts/seed_ko_list.json`（种子表）生成 `elements.json` 骨架。

两种模式：
- **preserve 模式**（推荐）：传入现有 `elements.json`，保留通路分组 / couplings /
  schematic 等用户定制，只把**KEGG 层字段**（complex / kegg_reaction / 兜底
  substrate+product）注入每个 KO
- **greenfield 模式**：根据 seed + KEGG module 自动分组出骨架，由用户后续补
  couplings / schematic

CLI：
    python -m envmeta kb-build --elements arsenic,nitrogen,sulfur,iron \\
        --preserve-from envmeta/geocycle/knowledge_base/elements.json \\
        --output envmeta/geocycle/knowledge_base/elements.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

_PKG_ROOT = Path(__file__).resolve().parent.parent
_DEFAULT_SNAPSHOT = _PKG_ROOT / "geocycle" / "kegg_snapshot.json"
_DEFAULT_SEED = _PKG_ROOT.parent / "scripts" / "seed_ko_list.json"


# --- 辅助 --------------------------------------------------------------------

def _load_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _normalize_compound(name: str | None) -> str | None:
    """把 KEGG 化合物名归一化为现 KB 风格简写（可选的友好化）。

    仅做最常见映射；未命中保持原样供用户人工决定。
    """
    if not name:
        return None
    map_common = {
        "Sulfate": "SO4-2",
        "Sulfite": "SO3-2",
        "Thiosulfate": "S2O3-2",
        "Hydrogen sulfide": "S-2",
        "Sulfur": "S0",
        "Nitrate": "NO3-",
        "Nitrite": "NO2-",
        "Ammonia": "NH3",
        "Nitric oxide": "NO",
        "Nitrous oxide": "N2O",
        "Arsenate": "As(V)",
        "Arsenite": "As(III)",
        "Fe(III)": "Fe(III)",
        "Fe(II)": "Fe(II)",
    }
    return map_common.get(name.strip(), name.strip())


# --- 核心 --------------------------------------------------------------------

def enrich_ko_fields(
    ko: str, gene_info: dict, snapshot: dict, overwrite: bool = False,
) -> dict:
    """对单个 KO 的 gene_info dict 注入 KEGG 字段。

    overwrite=False（默认）：仅补空字段，保留用户已写的值
    overwrite=True：强制覆盖（用于"KEGG-first"迁移）
    """
    kegg = snapshot.get("kos", {}).get(ko)
    if not kegg:
        return gene_info   # KEGG 里没这个 KO，保持原样

    # name：保留用户的短名（简洁），但记录 KEGG 官方长名到 kegg_name
    if kegg.get("symbol") and (overwrite or not gene_info.get("name")):
        gene_info["name"] = kegg["symbol"]
    if kegg.get("name"):
        gene_info["kegg_name"] = kegg["name"]

    # substrate/product：兜底填充（不覆盖用户手写风格如 "As(V)"）
    kegg_sub = _normalize_compound(kegg.get("substrate"))
    kegg_prod = _normalize_compound(kegg.get("product"))
    if overwrite or gene_info.get("substrate") is None:
        if gene_info.get("substrate") is None and kegg_sub:
            gene_info["substrate"] = kegg_sub
    if overwrite or gene_info.get("product") is None:
        if gene_info.get("product") is None and kegg_prod:
            gene_info["product"] = kegg_prod

    # complex: 取 module 第一个（通常最主要）
    mods = kegg.get("modules") or []
    if mods:
        gene_info["complex"] = mods[0]

    # kegg_reaction: 取 reaction 第一个
    rxs = kegg.get("reactions") or []
    if rxs:
        gene_info["kegg_reaction"] = rxs[0]

    return gene_info


def build_kb(
    elements: list[str],
    snapshot_path: Path | None = None,
    seed_path: Path | None = None,
    preserve_from: Path | None = None,
    overwrite: bool = False,
) -> dict:
    """构建 elements.json 骨架。"""
    snapshot = _load_json(snapshot_path or _DEFAULT_SNAPSHOT)
    seed = _load_json(seed_path or _DEFAULT_SEED)

    # Step 1: 初始化 elements 容器
    if preserve_from:
        old = _load_json(preserve_from)
        base = {
            "version": old.get("version", "2.0"),
            "source": old.get("source", ""),
            "elements": {},
            "couplings": old.get("couplings", []),
        }
        for eid in elements:
            if eid in old.get("elements", {}):
                base["elements"][eid] = json.loads(
                    json.dumps(old["elements"][eid])   # deep copy
                )
            else:
                base["elements"][eid] = _greenfield_element(eid, seed, snapshot)
    else:
        base = {
            "version": "2.0",
            "source": "EnvMeta kb-build (greenfield from KEGG snapshot)",
            "elements": {},
            "couplings": [],
        }
        for eid in elements:
            base["elements"][eid] = _greenfield_element(eid, seed, snapshot)

    # Step 2: 为每个元素里的每个 KO 注入 KEGG 字段
    total_ko = 0
    for eid, elem in base["elements"].items():
        for pw_name, pw in elem.get("pathways", {}).items():
            for ko, gene_info in pw.get("genes", {}).items():
                enrich_ko_fields(ko, gene_info, snapshot, overwrite=overwrite)
                total_ko += 1

    # Step 3: 记录来源
    base["kegg_source"] = {
        "version": snapshot.get("version"),
        "fetched_at": snapshot.get("fetched_at"),
        "license": snapshot.get("license", "academic"),
        "citation": snapshot.get("citation", ""),
    }
    base["ko_count"] = total_ko
    return base


def _greenfield_element(eid: str, seed: dict, snapshot: dict) -> dict:
    """没有 preserve_from 时从 KEGG module 自动分组生成元素结构。"""
    ko_list = seed.get(eid, [])
    if not ko_list:
        return {
            "display_name": {"en": eid.capitalize(), "cn": eid},
            "color": "#888888",
            "pathways": {},
        }
    # 按 module 聚合；无 module 的 KO 单独放入 "unclassified"
    from collections import OrderedDict
    pathway_groups: "OrderedDict[str, dict]" = OrderedDict()
    for ko in ko_list:
        k = snapshot.get("kos", {}).get(ko, {})
        mods = k.get("modules") or []
        pw_key = mods[0] if mods else "unclassified"
        pw_name = (snapshot.get("modules", {}).get(pw_key, {}).get("name")
                   if pw_key != "unclassified" else "Unclassified")
        if pw_key not in pathway_groups:
            pathway_groups[pw_key] = {
                "display_name": {"en": pw_name, "cn": pw_name},
                "genes": {},
            }
        pathway_groups[pw_key]["genes"][ko] = {
            "name": k.get("symbol", ko),
        }
    return {
        "display_name": {"en": eid.capitalize(), "cn": eid},
        "color": "#888888",
        "pathways": dict(pathway_groups),
    }


# --- CLI ---------------------------------------------------------------------

def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="envmeta kb-build",
        description="Build elements.json from KEGG snapshot + seed KO list.",
    )
    ap.add_argument("--elements", required=True,
                    help="Comma-separated element IDs, e.g. arsenic,nitrogen,sulfur,iron")
    ap.add_argument("--snapshot", type=Path, default=None,
                    help=f"KEGG snapshot JSON (default: {_DEFAULT_SNAPSHOT})")
    ap.add_argument("--seed", type=Path, default=None,
                    help=f"Seed KO list JSON (default: {_DEFAULT_SEED})")
    ap.add_argument("--preserve-from", type=Path, default=None,
                    help="Existing elements.json to preserve couplings/schematic/pathway groupings")
    ap.add_argument("--output", type=Path, required=True,
                    help="Output elements.json path")
    ap.add_argument("--overwrite", action="store_true",
                    help="Force KEGG values to overwrite existing substrate/product")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    elements = [x.strip() for x in args.elements.split(",") if x.strip()]
    kb = build_kb(
        elements,
        snapshot_path=args.snapshot,
        seed_path=args.seed,
        preserve_from=args.preserve_from,
        overwrite=args.overwrite,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(kb, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    n_mods = sum(
        1 for e in kb["elements"].values()
        for p in e.get("pathways", {}).values()
        for g in p.get("genes", {}).values()
        if "complex" in g
    )
    print(f"[OK] {args.output} written: {kb['ko_count']} KOs, "
          f"{n_mods} with KEGG complex field")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
