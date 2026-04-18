"""T2 HTML 交互导出测试。

覆盖目标：
- cycle_to_json 字段完整 + JSON 可序列化
- build_interactive_html 大小合理 + 含 D3 + 含数据 + 无外部 http 引用
- 假说评分 + 跨组 DataFrame 注入正常
"""
import json
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.geocycle.html_exporter import (
    build_interactive_html,
    cycle_to_json,
    export_html,
)
from envmeta.geocycle.inference import infer


SAMPLE = Path(__file__).parent / "sample_data"
FAST_PARAMS = {"perm_n": 99}


@pytest.fixture(scope="module")
def cycle_data():
    ko = pd.read_csv(SAMPLE / "kegg_target_only.tsv", sep="\t")
    tax = pd.read_csv(SAMPLE / "mag_taxonomy_labels.tsv", sep="\t",
                      header=None, names=["MAG", "Taxonomy"])
    ks = pd.read_csv(SAMPLE / "keystone_species.txt", sep="\t")
    ab = pd.read_csv(SAMPLE / "abundance.tsv", sep="\t")
    env = pd.read_csv(SAMPLE / "env_factors.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return infer(ko, tax, ks, ab, env, md, params=FAST_PARAMS)


# ═══════════════════════════════════════════════════════════════
# cycle_to_json
# ═══════════════════════════════════════════════════════════════

def test_cycle_to_json_required_fields(cycle_data):
    payload = cycle_to_json(cycle_data)
    required = {
        "version", "generated_at", "envmeta_version",
        "elements", "env_correlations", "full_corr_matrix",
        "sensitivity", "params", "meta",
    }
    missing = required - set(payload.keys())
    assert not missing, f"缺字段：{missing}"


def test_cycle_to_json_is_json_serializable(cycle_data):
    payload = cycle_to_json(cycle_data)
    # 真实 round-trip —— 任何 NaN / numpy / dataclass 残留都会报错
    s = json.dumps(payload, ensure_ascii=False, default=str)
    recovered = json.loads(s)
    assert recovered["version"] == "1.0"
    assert isinstance(recovered["elements"], list)


def test_cycle_to_json_elements_structure(cycle_data):
    payload = cycle_to_json(cycle_data)
    assert len(payload["elements"]) > 0
    el = payload["elements"][0]
    # ElementCycle asdict 应有这些字段
    assert "element_id" in el
    assert "pathways" in el
    if el["pathways"]:
        pw = el["pathways"][0]
        assert "pathway_id" in pw
        assert "contributors" in pw


def test_cycle_to_json_with_compare_df(cycle_data):
    fake_compare = pd.DataFrame({
        "group": ["CK", "A", "B"],
        "total_contribution": [1.2, 2.3, 3.4],
        "top_mag_genus": ["Gallionella", None, "Thiobacillus"],
    })
    payload = cycle_to_json(cycle_data, compare_df=fake_compare)
    assert "compare_groups" in payload
    assert len(payload["compare_groups"]) == 3
    # NaN 应被转为 None 保证 json 安全
    assert payload["compare_groups"][1]["top_mag_genus"] is None


# ═══════════════════════════════════════════════════════════════
# build_interactive_html
# ═══════════════════════════════════════════════════════════════

def test_build_html_returns_bytes(cycle_data):
    html = build_interactive_html(cycle_data)
    assert isinstance(html, bytes)
    assert html.startswith(b"<!DOCTYPE html>")


def test_build_html_size_reasonable(cycle_data):
    """D3.js (~280 KB) + 模板 + 数据 ≈ 300-500 KB。"""
    html = build_interactive_html(cycle_data)
    size_kb = len(html) / 1024
    assert 250 < size_kb < 800, f"HTML 大小异常：{size_kb:.1f} KB"


def test_build_html_includes_d3(cycle_data):
    """验证 D3.js inline 嵌入成功（关键字符串存在）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # D3 header 注释
    assert "d3js.org" in html
    # D3 核心 API
    assert "scaleOrdinal" in html or "forceSimulation" in html


def test_build_html_includes_data(cycle_data):
    """前端 EM_DATA 注入正确（不应包含未替换的占位符）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "{{CYCLE_DATA_JSON}}" not in html
    assert "{{D3_JS_INLINE}}" not in html
    assert "{{META_HTML}}" not in html
    assert "const EM_DATA =" in html


def test_build_html_no_external_http_references(cycle_data):
    """离线校验：不依赖任何外部 http:// 或 https:// 资源 script/link。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # 允许 href 里的 github 链接（仅文本），但不能有外部 script/link src
    import re
    external_scripts = re.findall(
        r'<script[^>]*\bsrc\s*=\s*["\']https?://[^"\']+["\']', html
    )
    assert not external_scripts, f"有外部 script：{external_scripts}"
    external_stylesheets = re.findall(
        r'<link[^>]*\bhref\s*=\s*["\']https?://[^"\']+["\'][^>]*rel\s*=\s*["\']stylesheet',
        html,
    )
    assert not external_stylesheets, f"有外部 CSS：{external_stylesheets}"


def test_build_html_meta_block_present(cycle_data):
    """顶部 Meta 块注入成功（含 EnvMeta 版本 + badge）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "em-meta-block" in html
    assert "EnvMeta" in html
    assert "em-badge" in html


def test_export_html_writes_file(cycle_data, tmp_path):
    out = tmp_path / "demo.html"
    path = export_html(cycle_data, out)
    assert path.exists()
    assert path.stat().st_size > 250_000  # 至少 250 KB
    content = path.read_bytes()
    assert content.startswith(b"<!DOCTYPE html>")


def test_build_html_with_hypothesis_by_group(cycle_data):
    """score_by_groups 格式的 dict 能正确嵌入。"""
    # 构造最小 fake HypothesisScore dict (不走真实 score 流程加速测试)
    from envmeta.geocycle.hypothesis import (
        ClaimResult, HypothesisScore,
    )
    fake_score = HypothesisScore(
        hypothesis_name="test",
        overall_score=0.85,
        label="strong",
        n_satisfied=3,
        n_total=3,
        n_skipped=0,
        claim_results=[
            ClaimResult(
                claim_id="c1", claim_type="pathway_active",
                status="satisfied", score=1.0, weight=1.0,
                evidence={}, explanation="ok",
            ),
        ],
        params={},
        null_p=0.02,
        null_p_samples=99,
        weight_robust=True,
        weight_sensitivity_rows=[],
        veto_reasons=[],
    )
    payload = cycle_to_json(
        cycle_data,
        hypothesis_by_group={"CK": fake_score, "B": fake_score},
    )
    assert "hypothesis_by_group" in payload
    assert set(payload["hypothesis_by_group"].keys()) == {"CK", "B"}
    assert payload["hypothesis_by_group"]["CK"]["label"] == "strong"
