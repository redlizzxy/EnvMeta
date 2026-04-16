"""Fork Bundle — 论文-EnvMeta 绑定发布协议（S4 / L4 层）。

把一篇论文用到的 KB + 假说 YAML + 分析配置打包成单个 .zip，
读者/审稿人一键加载复现。

Bundle 结构：
    bundle.zip
    ├── manifest.yaml
    ├── kb/
    │   ├── elements.json
    │   └── kegg_snapshot.json   (可选)
    ├── hypotheses/
    │   └── *.yaml               (可多个)
    ├── config/
    │   └── cycle_params.yaml
    └── README.md                (可选)

用法（Python API）：
    from envmeta.tools.bundle import create_bundle, load_bundle
    create_bundle("paper.zip", kb_path=..., hypothesis_paths=[...],
                  config={...}, manifest={...})
    bundle = load_bundle("paper.zip")
    bundle.kb, bundle.hypotheses, bundle.config, bundle.manifest

CLI：
    python -m envmeta bundle-create --kb ... --hypothesis ... --output ...
    python -m envmeta bundle-inspect paper.zip
"""
from __future__ import annotations

import datetime as _dt
import io
import json
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from envmeta import __version__ as _ENVMETA_VERSION
from envmeta.geocycle.hypothesis import Hypothesis, load_hypothesis


# =============================================================================
# 数据模型
# =============================================================================

MANIFEST_REQUIRED_FIELDS = ("name", "envmeta_version", "created")


@dataclass
class BundleContents:
    """解压后的 bundle 内容。

    kb: parsed elements.json (dict)
    kegg_snapshot: 可选 parsed kegg_snapshot.json (dict) 或 None
    hypotheses: list[Hypothesis]（已通过 load_hypothesis 加载）
    hypothesis_texts: list[tuple[filename, text]]（原始 YAML 文本，保留注释用）
    config: 分析参数字典
    manifest: 元信息字典
    readme: README.md 文本（若存在）
    """
    kb: dict
    kegg_snapshot: dict | None
    hypotheses: list[Hypothesis]
    hypothesis_texts: list[tuple[str, str]]
    config: dict
    manifest: dict
    readme: str | None = None


# =============================================================================
# 内部辅助
# =============================================================================

def _make_manifest(
    name: str,
    *,
    author: str = "",
    paper_doi: str = "",
    description: str = "",
    kegg_snapshot_date: str = "",
    hypothesis_names: list[str] | None = None,
    kb_files: list[str] | None = None,
    extra: dict | None = None,
) -> dict:
    """按约定格式组装 manifest 字典。"""
    m = {
        "name": name,
        "author": author,
        "paper_doi": paper_doi,
        "description": description,
        "envmeta_version": _ENVMETA_VERSION,
        "kegg_snapshot_date": kegg_snapshot_date,
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
        "files": {
            "kb": kb_files or ["kb/elements.json"],
            "hypotheses": hypothesis_names or [],
            "config": "config/cycle_params.yaml",
        },
    }
    if extra:
        m.update(extra)
    return m


def _validate_manifest(manifest: dict) -> list[str]:
    """返回错误列表；空列表 = 通过。"""
    errors: list[str] = []
    if not isinstance(manifest, dict):
        return ["manifest 必须是 dict / yaml mapping"]
    for f in MANIFEST_REQUIRED_FIELDS:
        if not manifest.get(f):
            errors.append(f"manifest 缺少必填字段 `{f}`")
    if "envmeta_version" in manifest and not isinstance(
            manifest["envmeta_version"], str):
        errors.append("manifest.envmeta_version 必须是字符串")
    return errors


# =============================================================================
# create_bundle
# =============================================================================

def create_bundle(
    output_path: str | Path,
    *,
    kb_path: str | Path,
    hypothesis_paths: list[str | Path] | None = None,
    config: dict | None = None,
    manifest: dict | None = None,
    name: str | None = None,
    author: str = "",
    paper_doi: str = "",
    description: str = "",
    kegg_snapshot_path: str | Path | None = None,
    readme_text: str | None = None,
) -> Path:
    """创建 Fork Bundle zip 文件。

    参数
    -----
    output_path        目标 .zip 路径
    kb_path            elements.json 路径（必填）
    hypothesis_paths   一个或多个 .yaml 文件路径
    config             cycle 分析参数字典（将写成 YAML）
    manifest           完整自定义 manifest（若给出则覆盖 name/author 等参数）
    name/author/...    若 manifest=None，用这些字段组装默认 manifest
    kegg_snapshot_path 可选 KEGG 快照 JSON 路径
    readme_text        可选 README.md 文本

    返回
    -----
    实际写入的 Path
    """
    output_path = Path(output_path)
    kb_path = Path(kb_path)
    if not kb_path.is_file():
        raise FileNotFoundError(f"KB 文件不存在: {kb_path}")

    hypothesis_paths = [Path(p) for p in (hypothesis_paths or [])]
    for p in hypothesis_paths:
        if not p.is_file():
            raise FileNotFoundError(f"假说 YAML 不存在: {p}")

    if kegg_snapshot_path is not None:
        kegg_snapshot_path = Path(kegg_snapshot_path)
        if not kegg_snapshot_path.is_file():
            raise FileNotFoundError(f"KEGG 快照不存在: {kegg_snapshot_path}")

    hyp_names = [f"hypotheses/{p.name}" for p in hypothesis_paths]
    kb_files = ["kb/elements.json"] + (
        ["kb/kegg_snapshot.json"] if kegg_snapshot_path else []
    )

    if manifest is None:
        manifest = _make_manifest(
            name=name or output_path.stem,
            author=author,
            paper_doi=paper_doi,
            description=description,
            hypothesis_names=hyp_names,
            kb_files=kb_files,
        )
    else:
        # 确保必填字段存在；补上版本 + 时间戳
        manifest = dict(manifest)
        manifest.setdefault("envmeta_version", _ENVMETA_VERSION)
        manifest.setdefault("created",
                             _dt.datetime.now().isoformat(timespec="seconds"))
    errors = _validate_manifest(manifest)
    if errors:
        raise ValueError("manifest 校验失败：\n  " + "\n  ".join(errors))

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as z:
        # manifest.yaml
        z.writestr(
            "manifest.yaml",
            yaml.safe_dump(manifest, allow_unicode=True, sort_keys=False),
        )
        # kb
        z.write(kb_path, "kb/elements.json")
        if kegg_snapshot_path is not None:
            z.write(kegg_snapshot_path, "kb/kegg_snapshot.json")
        # hypotheses（保留原始文本含注释）
        for p in hypothesis_paths:
            z.write(p, f"hypotheses/{p.name}")
        # config
        z.writestr(
            "config/cycle_params.yaml",
            yaml.safe_dump(config or {}, allow_unicode=True, sort_keys=False),
        )
        # readme（可选）
        if readme_text:
            z.writestr("README.md", readme_text)

    return output_path


# =============================================================================
# load_bundle
# =============================================================================

def load_bundle(bundle_path: str | Path) -> BundleContents:
    """加载 Fork Bundle zip → BundleContents。"""
    bundle_path = Path(bundle_path)
    if not bundle_path.is_file():
        raise FileNotFoundError(f"Bundle 不存在: {bundle_path}")

    with zipfile.ZipFile(bundle_path, "r") as z:
        namelist = set(z.namelist())

        # manifest 必需
        if "manifest.yaml" not in namelist:
            raise ValueError(f"Bundle 缺少 manifest.yaml: {bundle_path}")
        manifest = yaml.safe_load(z.read("manifest.yaml").decode("utf-8"))
        errors = _validate_manifest(manifest)
        if errors:
            raise ValueError("manifest 校验失败：\n  " + "\n  ".join(errors))

        # KB 必需
        if "kb/elements.json" not in namelist:
            raise ValueError("Bundle 缺少 kb/elements.json")
        kb = json.loads(z.read("kb/elements.json").decode("utf-8"))

        # KEGG snapshot 可选
        kegg_snapshot = None
        if "kb/kegg_snapshot.json" in namelist:
            kegg_snapshot = json.loads(
                z.read("kb/kegg_snapshot.json").decode("utf-8"))

        # Hypothesis YAML(s)
        hypotheses: list[Hypothesis] = []
        hypothesis_texts: list[tuple[str, str]] = []
        for name in sorted(namelist):
            if name.startswith("hypotheses/") and name.endswith((".yaml", ".yml")):
                text = z.read(name).decode("utf-8")
                hypothesis_texts.append((name, text))
                hypotheses.append(load_hypothesis(text))

        # Config
        config: dict = {}
        if "config/cycle_params.yaml" in namelist:
            raw = yaml.safe_load(
                z.read("config/cycle_params.yaml").decode("utf-8"))
            config = raw or {}

        # README
        readme = None
        if "README.md" in namelist:
            readme = z.read("README.md").decode("utf-8")

    return BundleContents(
        kb=kb,
        kegg_snapshot=kegg_snapshot,
        hypotheses=hypotheses,
        hypothesis_texts=hypothesis_texts,
        config=config,
        manifest=manifest,
        readme=readme,
    )


# =============================================================================
# inspect_bundle —— 浅层检查，不完整加载（快速审查）
# =============================================================================

def inspect_bundle(bundle_path: str | Path) -> dict:
    """返回 bundle 的 summary 字典（供 CLI / UI 展示）。

    不会反序列化 KB 或假说体 —— 仅清点文件 + 读 manifest。
    """
    bundle_path = Path(bundle_path)
    with zipfile.ZipFile(bundle_path, "r") as z:
        namelist = z.namelist()
        manifest_raw = None
        if "manifest.yaml" in namelist:
            manifest_raw = yaml.safe_load(
                z.read("manifest.yaml").decode("utf-8"))

    hyp_files = [n for n in namelist
                 if n.startswith("hypotheses/")
                 and n.endswith((".yaml", ".yml"))]
    has_kegg = "kb/kegg_snapshot.json" in namelist
    has_readme = "README.md" in namelist
    version_mismatch = False
    if manifest_raw and manifest_raw.get("envmeta_version") not in (
            None, "", _ENVMETA_VERSION):
        version_mismatch = True

    return {
        "path": str(bundle_path),
        "manifest": manifest_raw or {},
        "n_files": len(namelist),
        "has_kb": "kb/elements.json" in namelist,
        "has_kegg_snapshot": has_kegg,
        "has_readme": has_readme,
        "hypothesis_files": hyp_files,
        "envmeta_version_runtime": _ENVMETA_VERSION,
        "envmeta_version_bundle": (manifest_raw or {}).get("envmeta_version"),
        "version_mismatch": version_mismatch,
    }


# =============================================================================
# CLI: bundle-create / bundle-inspect
# =============================================================================

def _cli_create(argv: list[str]) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="envmeta bundle-create",
        description="创建 Fork Bundle zip",
    )
    parser.add_argument("--output", "-o", required=True,
                        help="目标 .zip 路径")
    parser.add_argument("--kb", required=True,
                        help="elements.json 路径")
    parser.add_argument("--hypothesis", "-H", action="append", default=[],
                        help="假说 YAML 路径（可多次）")
    parser.add_argument("--kegg-snapshot", default=None,
                        help="可选 KEGG 快照 JSON 路径")
    parser.add_argument("--config", default=None,
                        help="可选分析参数 YAML 路径")
    parser.add_argument("--name", default=None, help="bundle 名")
    parser.add_argument("--author", default="", help="作者")
    parser.add_argument("--paper-doi", default="", help="论文 DOI")
    parser.add_argument("--description", default="", help="描述")
    parser.add_argument("--readme", default=None,
                        help="可选 README.md 文件路径")
    args = parser.parse_args(argv)

    config: dict = {}
    if args.config:
        config = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}

    readme_text = None
    if args.readme:
        readme_text = Path(args.readme).read_text(encoding="utf-8")

    path = create_bundle(
        args.output,
        kb_path=args.kb,
        hypothesis_paths=args.hypothesis,
        config=config,
        name=args.name,
        author=args.author,
        paper_doi=args.paper_doi,
        description=args.description,
        kegg_snapshot_path=args.kegg_snapshot,
        readme_text=readme_text,
    )
    print(f"[OK] bundle written to {path}")
    return 0


def _cli_inspect(argv: list[str]) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="envmeta bundle-inspect",
        description="检查 Fork Bundle 内容 + 版本兼容性",
    )
    parser.add_argument("path", help="Bundle .zip 路径")
    args = parser.parse_args(argv)

    info = inspect_bundle(args.path)
    print(f"Bundle: {info['path']}")
    m = info["manifest"] or {}
    print(f"  name:         {m.get('name', '(missing)')}")
    print(f"  author:       {m.get('author', '')}")
    print(f"  paper_doi:    {m.get('paper_doi', '')}")
    print(f"  created:      {m.get('created', '(missing)')}")
    print(f"  envmeta ver:  bundle={info['envmeta_version_bundle']} "
          f"runtime={info['envmeta_version_runtime']}"
          + ("  [!] MISMATCH" if info["version_mismatch"] else ""))
    print(f"  files:        {info['n_files']}")
    print(f"    kb:          {'yes' if info['has_kb'] else 'NO'}")
    print(f"    kegg_snap:   {'yes' if info['has_kegg_snapshot'] else 'no'}")
    print(f"    readme:      {'yes' if info['has_readme'] else 'no'}")
    print(f"    hypotheses:  {len(info['hypothesis_files'])}")
    for h in info["hypothesis_files"]:
        print(f"      - {h}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Not used directly — dispatched from envmeta.__main__ via two subcommands."""
    parser = argparse.ArgumentParser(prog="envmeta bundle")
    parser.add_argument("sub", choices=["create", "inspect"])
    args, rest = parser.parse_known_args(argv)
    if args.sub == "create":
        return _cli_create(rest)
    return _cli_inspect(rest)


# 暴露给 __main__.py 的两个 entry point
cli_create = _cli_create
cli_inspect = _cli_inspect

# 避免 NameError（main() 内部用 argparse 但懒加载）
import argparse  # noqa: E402
