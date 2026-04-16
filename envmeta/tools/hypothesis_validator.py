"""假说 YAML 校验器 CLI（S3.5）。

用法：
    python -m envmeta hypothesis-validate path/to/hypothesis.yaml

校验内容：
1. YAML 结构合法（复用 load_hypothesis）
2. 每条 claim 的 params 引用 KB 时是否存在
   - pathway_active / keystone_in_pathway / group_contrast：pathway 名
   - coupling_possible：species_a / species_b
   - env_correlation：env_factor 无法离线校验（只给出信息性提示）

退出码 = errors 数（0 = 完全通过）
"""
from __future__ import annotations

import argparse
import difflib
from pathlib import Path

from envmeta.geocycle.hypothesis import Claim, load_hypothesis
from envmeta.geocycle.knowledge_base import (
    couplings as kb_couplings,
    load_kb,
    pathway_display,
    pathway_ko_sets,
)


def _norm(s: str) -> str:
    return s.strip().lower().replace("_", " ").replace("-", " ")


def _suggest(name: str, candidates: list[str], n: int = 3) -> list[str]:
    """给 name 返回最多 n 个最相近的候选。"""
    norm_map = {_norm(c): c for c in candidates}
    hits = difflib.get_close_matches(_norm(name), list(norm_map.keys()),
                                      n=n, cutoff=0.4)
    return [norm_map[h] for h in hits]


def validate_claim(claim: Claim) -> tuple[list[str], list[str]]:
    """返回 (errors, warnings)。errors = 字段缺失 / 非法；warnings = 引用不到。"""
    errors: list[str] = []
    warnings: list[str] = []
    p = claim.params or {}
    ctype = claim.type

    # 预加载 KB 引用
    kb = load_kb()
    pw_ko_sets = pathway_ko_sets(kb)            # {pathway_id: [ko, ...]}
    pw_disp = pathway_display(kb, lang="en")    # {pathway_id: display_name}
    # 可查询名字集：display_name + pathway_id
    pathway_names = set(pw_disp.values()) | set(pw_ko_sets.keys())
    pathway_names_norm = {_norm(n): n for n in pathway_names}

    cpls = kb_couplings()
    species_set: set[str] = set()
    for cp in cpls:
        a = cp.get("species_a")
        b = cp.get("species_b")
        if a:
            species_set.add(str(a))
        if b:
            species_set.add(str(b))

    def _check_pathway(name: str) -> None:
        if not name:
            errors.append(f"{claim.id}: params.pathway 缺失")
            return
        if _norm(name) not in pathway_names_norm:
            suggestions = _suggest(name, sorted(pathway_names))
            msg = f"{claim.id}: 通路 {name!r} 在 KB 找不到"
            if suggestions:
                msg += f"（相近候选: {', '.join(repr(s) for s in suggestions)}）"
            warnings.append(msg)

    if ctype == "pathway_active":
        _check_pathway(p.get("pathway", ""))
    elif ctype == "keystone_in_pathway":
        _check_pathway(p.get("pathway", ""))
        if "min_keystones" in p:
            try:
                int(p["min_keystones"])
            except (TypeError, ValueError):
                errors.append(f"{claim.id}: min_keystones 必须是整数")
    elif ctype == "coupling_possible":
        a = p.get("species_a")
        b = p.get("species_b")
        if not a or not b:
            errors.append(f"{claim.id}: species_a / species_b 必填")
        else:
            # KB 里是否存在这对（a,b）或 (b,a) 配对
            found = any(
                (cp.get("species_a") == a and cp.get("species_b") == b)
                or (cp.get("species_a") == b and cp.get("species_b") == a)
                for cp in cpls
            )
            if not found:
                # 分别给两端的相近候选
                msg_parts = [f"{claim.id}: KB 里无 {a}↔{b} 耦合定义"]
                for side, val in (("species_a", a), ("species_b", b)):
                    if val not in species_set:
                        hits = _suggest(val, sorted(species_set))
                        if hits:
                            msg_parts.append(
                                f"  {side}={val!r} 不在 KB；相近: "
                                f"{', '.join(repr(s) for s in hits)}"
                            )
                warnings.append("\n".join(msg_parts))
    elif ctype == "env_correlation":
        if not p.get("pathway"):
            errors.append(f"{claim.id}: params.pathway 必填")
        else:
            _check_pathway(p["pathway"])
        if not p.get("env_factor"):
            errors.append(f"{claim.id}: params.env_factor 必填")
        # env_factor 只有数据上下文才能校验
        exp_sign = (p.get("expected_sign") or "any").lower()
        if exp_sign not in ("positive", "negative", "any"):
            errors.append(
                f"{claim.id}: expected_sign={exp_sign!r} 非法，"
                f"应为 positive/negative/any"
            )
        min_conf = (p.get("min_confidence") or "suggestive").lower()
        if min_conf not in ("strong", "suggestive", "weak"):
            errors.append(
                f"{claim.id}: min_confidence={min_conf!r} 非法，"
                f"应为 strong/suggestive/weak"
            )
    elif ctype == "group_contrast":
        _check_pathway(p.get("pathway", ""))
        if not p.get("high_group"):
            errors.append(f"{claim.id}: params.high_group 必填")
        if not p.get("low_group"):
            errors.append(f"{claim.id}: params.low_group 必填")
        try:
            ratio = float(p.get("min_ratio", 1.5))
            if ratio <= 0:
                errors.append(f"{claim.id}: min_ratio 必须 > 0")
        except (TypeError, ValueError):
            errors.append(f"{claim.id}: min_ratio 必须是数字")
    return errors, warnings


def validate_file(path: str | Path) -> dict:
    """返回 {errors, warnings, n_claims}。"""
    path = Path(path)
    all_errors: list[str] = []
    all_warnings: list[str] = []
    try:
        hyp = load_hypothesis(path)
    except ValueError as e:
        return {
            "errors": [f"YAML 解析失败: {e}"],
            "warnings": [],
            "n_claims": 0,
        }
    for claim in hyp.claims:
        errs, warns = validate_claim(claim)
        all_errors.extend(errs)
        all_warnings.extend(warns)
    return {
        "errors": all_errors,
        "warnings": all_warnings,
        "n_claims": len(hyp.claims),
        "hypothesis_name": hyp.name,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="envmeta hypothesis-validate",
        description="校验假说 YAML 对 KB 的引用正确性（S3.5）",
    )
    parser.add_argument("path", type=str, help="假说 YAML 文件路径")
    args = parser.parse_args(argv)

    result = validate_file(args.path)
    print(f"文件: {args.path}")
    if "hypothesis_name" in result:
        print(f"假说: {result['hypothesis_name']}")
    print(f"Claim 数: {result['n_claims']}")

    if result["errors"]:
        print("\n[X] ERRORS:")
        for e in result["errors"]:
            print(f"  [X] {e}")
    if result["warnings"]:
        print("\n[!] WARNINGS:")
        for w in result["warnings"]:
            for line in w.split("\n"):
                print(f"  [!] {line}")
    if not result["errors"] and not result["warnings"]:
        print("\n[OK] all checks passed")
    print(
        f"\nSummary: {len(result['errors'])} error(s), "
        f"{len(result['warnings'])} warning(s)."
    )
    return len(result["errors"])


if __name__ == "__main__":
    raise SystemExit(main())
