"""Fork Bundle 测试（S4）。"""
from pathlib import Path

import pytest
import yaml

from envmeta.tools.bundle import (
    BundleContents,
    create_bundle,
    inspect_bundle,
    load_bundle,
)

SAMPLE_KB = Path("envmeta/geocycle/knowledge_base/elements.json")
SAMPLE_HYP = Path("paper/hypotheses/arsenic_steel_slag.yaml")


# ── 基础 round-trip ─────────────────────────────────────────

def test_create_and_load_round_trip(tmp_path):
    out = tmp_path / "paper.zip"
    create_bundle(
        out,
        kb_path=SAMPLE_KB,
        hypothesis_paths=[SAMPLE_HYP],
        config={"completeness_threshold": 50, "env_rho_min": 0.5},
        name="Test",
        author="A",
        paper_doi="10.xxx/test",
        description="round trip",
    )
    assert out.is_file()
    b = load_bundle(out)
    assert isinstance(b, BundleContents)
    assert b.manifest["name"] == "Test"
    assert b.manifest["author"] == "A"
    assert b.manifest["paper_doi"] == "10.xxx/test"
    assert "elements" in b.kb
    assert len(b.hypotheses) == 1
    assert b.hypotheses[0].name  # 已通过 load_hypothesis
    assert b.config["completeness_threshold"] == 50
    assert b.config["env_rho_min"] == 0.5


def test_multiple_hypotheses(tmp_path):
    # 制作第二个小 YAML
    hyp2 = tmp_path / "alt.yaml"
    hyp2.write_text(
        "name: Alt Hypothesis\n"
        "claims:\n"
        "  - id: c1\n"
        "    type: pathway_active\n"
        "    params: {pathway: X}\n",
        encoding="utf-8",
    )
    out = tmp_path / "multi.zip"
    create_bundle(
        out,
        kb_path=SAMPLE_KB,
        hypothesis_paths=[SAMPLE_HYP, hyp2],
        name="Multi",
    )
    b = load_bundle(out)
    assert len(b.hypotheses) == 2
    names = {h.name for h in b.hypotheses}
    assert "Alt Hypothesis" in names


def test_hypothesis_text_preserves_comments(tmp_path):
    """bundle 应保留 YAML 原始文本（含注释），不是只 load 后重 dump。"""
    out = tmp_path / "p.zip"
    create_bundle(out, kb_path=SAMPLE_KB, hypothesis_paths=[SAMPLE_HYP],
                  name="X")
    b = load_bundle(out)
    assert b.hypothesis_texts
    first_name, first_text = b.hypothesis_texts[0]
    # 示例 YAML 里有 "§ 1 使用前必读"（理论段注释），round-trip 不丢
    assert "使用前必读" in first_text


# ── manifest 校验 ──────────────────────────────────────────

def test_missing_kb_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="KB"):
        create_bundle(
            tmp_path / "p.zip",
            kb_path=tmp_path / "nope.json",
            hypothesis_paths=[],
        )


def test_missing_hypothesis_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="假说"):
        create_bundle(
            tmp_path / "p.zip",
            kb_path=SAMPLE_KB,
            hypothesis_paths=[tmp_path / "nope.yaml"],
        )


def test_bad_manifest_rejected(tmp_path):
    with pytest.raises(ValueError, match="manifest"):
        create_bundle(
            tmp_path / "p.zip",
            kb_path=SAMPLE_KB,
            hypothesis_paths=[],
            manifest={"author": "no-name"},  # 缺 name
        )


# ── inspect ────────────────────────────────────────────────

def test_inspect_reports_contents(tmp_path):
    out = tmp_path / "p.zip"
    create_bundle(out, kb_path=SAMPLE_KB, hypothesis_paths=[SAMPLE_HYP],
                  name="Inspect Test")
    info = inspect_bundle(out)
    assert info["manifest"]["name"] == "Inspect Test"
    assert info["has_kb"] is True
    assert info["has_kegg_snapshot"] is False
    assert len(info["hypothesis_files"]) == 1
    assert info["version_mismatch"] is False


def test_inspect_flags_version_mismatch(tmp_path):
    """manually-crafted bundle with old envmeta_version 应被标记。"""
    import zipfile
    out = tmp_path / "old.zip"
    with zipfile.ZipFile(out, "w") as z:
        z.writestr("manifest.yaml", yaml.safe_dump({
            "name": "Old Bundle",
            "envmeta_version": "0.0.1-ancient",
            "created": "2020-01-01",
        }))
        z.writestr("kb/elements.json", "{}")
    info = inspect_bundle(out)
    assert info["version_mismatch"] is True


# ── load_bundle 错误路径 ────────────────────────────────────

def test_load_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_bundle(tmp_path / "nope.zip")


def test_empty_manifest_fields_are_omitted(tmp_path):
    """空字符串字段（paper_doi / kegg_snapshot_date / author / description）
    不应写入 manifest（保持 YAML 整洁）。"""
    out = tmp_path / "clean.zip"
    create_bundle(
        out,
        kb_path=SAMPLE_KB,
        hypothesis_paths=[],
        name="Clean",
        # author / paper_doi / description / kegg_snapshot_date 全为空
    )
    b = load_bundle(out)
    m = b.manifest
    # 必填字段仍在
    assert m["name"] == "Clean"
    assert "envmeta_version" in m
    assert "created" in m
    # 空字段应被省略
    assert "author" not in m
    assert "paper_doi" not in m
    assert "description" not in m
    assert "kegg_snapshot_date" not in m


def test_load_bundle_missing_kb(tmp_path):
    import zipfile
    bad = tmp_path / "bad.zip"
    with zipfile.ZipFile(bad, "w") as z:
        z.writestr("manifest.yaml", yaml.safe_dump({
            "name": "x", "envmeta_version": "0.1.0", "created": "now",
        }))
        # no kb/
    with pytest.raises(ValueError, match="elements.json"):
        load_bundle(bad)
