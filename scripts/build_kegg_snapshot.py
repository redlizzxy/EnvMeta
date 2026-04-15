"""KEGG 快照构建器（S2.5-10a）。

从 KEGG REST API 拉取种子 KO 列表的官方数据（name / pathways / modules /
reactions / substrate / product），打包成离线 JSON 快照供 `envmeta kb-build`
使用。

**用法**（EnvMeta 维护者做一次，需联网）：

    python scripts/build_kegg_snapshot.py \\
        --seed scripts/seed_ko_list.json \\
        --output envmeta/geocycle/kegg_snapshot.json \\
        [--kegg-version 2024.10]

**限制**：
- KEGG REST 限速 ~3 req/s（默认每请求 sleep 0.35s）
- 偶发 5xx/连接错误重试 3 次
- 本脚本**不做**增量更新；整体重跑一次约 1000 KO × 3 请求 ≈ 15 min

**KEGG 许可**：学术研究免费；商业须付费订阅。参考
https://www.kegg.jp/kegg/legal.html。输出 JSON 带 `license: academic` 标记。

**Citation**：Kanehisa M, Goto S. (2000) Nucleic Acids Res. 28:27-30.
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

KEGG_BASE = "https://rest.kegg.jp"
DEFAULT_SLEEP = 0.35   # 秒；KEGG 建议 ≤ 3 req/s
MAX_RETRIES = 3
TIMEOUT = 15


def _http_get(url: str, sleep: float = DEFAULT_SLEEP) -> str:
    """GET + 重试；失败则 raise。返回 body 文本。"""
    last_err: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT) as resp:
                body = resp.read().decode("utf-8", errors="replace")
            time.sleep(sleep)
            return body
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            last_err = e
            time.sleep(1.0 * (attempt + 1))
    raise RuntimeError(f"KEGG GET failed after {MAX_RETRIES} retries: {url}: {last_err}")


def _parse_flat(body: str) -> dict[str, list[str]]:
    """把 KEGG flat text 解析成 {section: [lines]}。section 开头大写，续行缩进。"""
    out: dict[str, list[str]] = {}
    current = None
    for line in body.splitlines():
        if not line:
            continue
        if line[0] != " ":
            # 新 section
            parts = line.split(None, 1)
            current = parts[0]
            out.setdefault(current, [])
            if len(parts) > 1:
                out[current].append(parts[1])
        else:
            if current is not None:
                out[current].append(line.strip())
    return out


def fetch_ko(ko: str) -> dict[str, Any]:
    """拉 KO 条目。返回 symbol / name / definition / pathways / modules / reactions。

    KEGG KO flat 顶级 sections：SYMBOL（短名）/ NAME（长名+EC） / DEFINITION /
    PATHWAY / MODULE / REACTION / BRITE / DBLINKS / GENES / REFERENCE 等。
    """
    text = _http_get(f"{KEGG_BASE}/get/ko:{ko}")
    sections = _parse_flat(text)

    def _first(sec: str) -> str:
        return (sections.get(sec) or [""])[0].strip()

    symbol = _first("SYMBOL").split(",")[0].strip()  # sqr / narG / ...
    name_full = _first("NAME")                        # long name + EC
    definition = " ".join(sections.get("DEFINITION", [])).strip()

    pathways: list[str] = []
    for ln in sections.get("PATHWAY", []):
        m = re.match(r"(ko|map)(\d{5})\s", ln)
        if m:
            pathways.append("map" + m.group(2))

    modules: list[str] = []
    for ln in sections.get("MODULE", []):
        m = re.match(r"(M\d{5})\s", ln)
        if m:
            modules.append(m.group(1))

    reactions: list[str] = []
    for ln in sections.get("REACTION", []):
        reactions += re.findall(r"R\d{5}", ln)

    return {
        "symbol": symbol,
        "name": name_full,
        "definition": definition,
        "pathways": pathways,
        "modules": modules,
        "reactions": list(dict.fromkeys(reactions)),   # 去重保序
    }


def fetch_reaction(rn: str) -> dict[str, Any]:
    """拉 REACTION 条目，解析 substrate / product equation。"""
    try:
        text = _http_get(f"{KEGG_BASE}/get/rn:{rn}")
    except Exception:
        return {"equation": "", "substrate": None, "product": None}
    sections = _parse_flat(text)
    eq = " ".join(sections.get("EQUATION", [])).strip()
    substrate, product = _parse_equation(eq)
    return {"equation": eq, "substrate": substrate, "product": product}


def _parse_equation(eq: str) -> tuple[str | None, str | None]:
    """从 KEGG EQUATION 提取主 substrate / product 的名字。

    KEGG eq 形如 `C00283 + 2 C00030 <=> C00087 + 2 C00001`（化合物 ID），
    这里只抓首个 substrate / product 的 C ID（上游再 resolve 成名字或直接用 ID）。
    """
    if not eq:
        return None, None
    if "<=>" in eq:
        left, right = eq.split("<=>", 1)
    elif "=>" in eq:
        left, right = eq.split("=>", 1)
    elif "->" in eq:
        left, right = eq.split("->", 1)
    else:
        return None, None

    def _first_c(s: str) -> str | None:
        m = re.search(r"C\d{5}", s)
        return m.group(0) if m else None

    return _first_c(left), _first_c(right)


def fetch_compound_name(cid: str) -> str:
    """把 C##### 化合物 ID 转名字（第一条 NAME）。"""
    if not cid:
        return ""
    try:
        text = _http_get(f"{KEGG_BASE}/get/cpd:{cid}")
    except Exception:
        return cid
    sections = _parse_flat(text)
    if "NAME" in sections:
        return sections["NAME"][0].rstrip(";").strip()
    return cid


def fetch_module(mid: str) -> dict[str, Any]:
    """拉 MODULE 条目，主要取 NAME + ORTHOLOGY 里的 KO 成员。"""
    try:
        text = _http_get(f"{KEGG_BASE}/get/md:{mid}")
    except Exception:
        return {"name": "", "kos": []}
    sections = _parse_flat(text)
    name = " ".join(sections.get("NAME", [])).strip()
    kos: list[str] = []
    for ln in sections.get("ORTHOLOGY", []):
        kos += re.findall(r"K\d{5}", ln)
    return {"name": name, "kos": sorted(set(kos))}


def build_snapshot(seed_path: Path, kegg_version: str) -> dict[str, Any]:
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    all_kos: set[str] = set()
    for k, vs in seed.items():
        if k.startswith("_"):
            continue
        all_kos.update(vs)
    all_kos = sorted(all_kos)

    kos_data: dict[str, Any] = {}
    modules_seen: set[str] = set()
    compound_cache: dict[str, str] = {}

    print(f"Fetching {len(all_kos)} KOs...", file=sys.stderr)
    for i, ko in enumerate(all_kos, 1):
        print(f"  [{i:3d}/{len(all_kos)}] {ko}", file=sys.stderr)
        info = fetch_ko(ko)

        # 取第一个 reaction 当主反应；substrate/product 解析为化合物名
        sub_name = None
        prod_name = None
        if info["reactions"]:
            rn = info["reactions"][0]
            rinfo = fetch_reaction(rn)
            sub_id = rinfo.get("substrate")
            prod_id = rinfo.get("product")
            if sub_id:
                sub_name = compound_cache.get(sub_id) or fetch_compound_name(sub_id)
                compound_cache[sub_id] = sub_name
            if prod_id:
                prod_name = compound_cache.get(prod_id) or fetch_compound_name(prod_id)
                compound_cache[prod_id] = prod_name

        kos_data[ko] = {
            "symbol": info["symbol"],            # 短名 (sqr / narG / ...)
            "name": info["name"],                # 长描述 + EC
            "definition": info["definition"],
            "pathways": info["pathways"],
            "modules": info["modules"],
            "reactions": info["reactions"],
            "substrate": sub_name,
            "product": prod_name,
        }
        modules_seen.update(info["modules"])

    # 拉 module 详情
    modules_data: dict[str, Any] = {}
    print(f"\nFetching {len(modules_seen)} modules...", file=sys.stderr)
    for i, mid in enumerate(sorted(modules_seen), 1):
        print(f"  [{i:3d}/{len(modules_seen)}] {mid}", file=sys.stderr)
        modules_data[mid] = fetch_module(mid)

    return {
        "version": kegg_version,
        "fetched_at": datetime.date.today().isoformat(),
        "license": "academic",
        "source": "KEGG REST API (https://rest.kegg.jp)",
        "citation": "Kanehisa M, Goto S. (2000) Nucleic Acids Res. 28:27-30.",
        "kos": kos_data,
        "modules": modules_data,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    ap.add_argument("--seed", type=Path, required=True,
                    help="种子 KO 列表 JSON（按元素分组）")
    ap.add_argument("--output", type=Path, required=True,
                    help="输出快照 JSON 路径")
    ap.add_argument("--kegg-version", default=None,
                    help="KEGG 版本号/日期标签（默认：今日）")
    args = ap.parse_args()

    version = args.kegg_version or f"kegg_{datetime.date.today().isoformat()}"
    snap = build_snapshot(args.seed, version)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(snap, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    n_kos = len(snap["kos"])
    n_mods = len(snap["modules"])
    size_kb = args.output.stat().st_size / 1024
    print(f"\n[OK] snapshot saved: {args.output} "
          f"({n_kos} KOs, {n_mods} modules, {size_kb:.1f} KB)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
