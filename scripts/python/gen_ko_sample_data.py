"""从 eggnog.KEGG_ko.TPM.spf 过滤出知识库里的 57 个 KO，生成 tests/sample_data/ko_tpm.spf。

一次性生成脚本。重跑无副作用，会覆盖已有文件。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
KB_PATH = ROOT / "envmeta" / "geocycle" / "knowledge_base" / "elements.json"
SOURCE = ROOT / "data" / "raw" / "eggnog.KEGG_ko.TPM.spf"
DEST = ROOT / "tests" / "sample_data" / "ko_tpm.spf"


def collect_kos(kb: dict) -> set[str]:
    kos = set()
    for _, el in kb["elements"].items():
        for _, pw in el["pathways"].items():
            kos.update(pw["genes"].keys())
    return kos


def main() -> int:
    if not SOURCE.exists():
        print(f"❌ 找不到源文件 {SOURCE}", file=sys.stderr)
        return 1

    kb = json.loads(KB_PATH.read_text(encoding="utf-8"))
    target_kos = collect_kos(kb)
    print(f"知识库 KO 数：{len(target_kos)}")

    df = pd.read_csv(SOURCE, sep="\t", dtype=str, keep_default_na=False)
    # 第 2 列 "KEGG_ko" 是 KO 列
    ko_col = df.columns[1]
    filtered = df[df[ko_col].isin(target_kos)].copy()
    print(f"源表 {len(df)} 行 → 过滤后 {len(filtered)} 行")

    DEST.parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(DEST, sep="\t", index=False)
    print(f"✅ 保存到 {DEST} ({DEST.stat().st_size / 1024:.1f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
