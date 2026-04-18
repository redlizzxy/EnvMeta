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
        "sensitivity", "couplings", "params", "meta",
    }
    missing = required - set(payload.keys())
    assert not missing, f"缺字段：{missing}"


def test_cycle_to_json_couplings_from_kb(cycle_data):
    """T2-β: couplings 应从 KB 读取，至少含 4 条（As/N/S/Fe 对应的预设耦合）。"""
    payload = cycle_to_json(cycle_data)
    couplings = payload.get("couplings", [])
    assert len(couplings) >= 3, f"couplings 太少：{len(couplings)}"
    # 每条必须有 species_a/species_b/product/type/color
    for cp in couplings:
        assert "species_a" in cp and "species_b" in cp
        assert "product" in cp and "type" in cp
        assert cp.get("color", "").startswith("#") or cp.get("color") is None


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


def test_cycle_to_json_hypothesis_by_group_df(cycle_data):
    """T2-γ: score_by_groups 返回 DataFrame 时应存为 hypothesis_by_group_summary。"""
    fake_df = pd.DataFrame({
        "group": ["CK", "A", "B"],
        "overall_score": [0.8, 0.7, 0.9],
        "label": ["strong", "suggestive", "strong"],
        "null_p": [0.04, 0.3, 0.01],
    })
    payload = cycle_to_json(cycle_data, hypothesis_by_group=fake_df)
    assert "hypothesis_by_group_summary" in payload
    assert len(payload["hypothesis_by_group_summary"]) == 3
    assert payload["hypothesis_by_group_summary"][0]["group"] == "CK"


# ═══════════════════════════════════════════════════════════════
# T2-ε 精调：节点筛选标准 / env panel / 耦合精确锚点 / 布局反转
# ═══════════════════════════════════════════════════════════════

def test_meta_block_shows_node_selection_standard(cycle_data):
    """ε.1: meta 块应明示筛选标准（completeness 阈值 / top N / ranking）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "节点筛选标准" in html
    assert "completeness ≥" in html


def test_html_contains_render_env_panel(cycle_data):
    """ε.2: env panel 实际渲染函数存在（非占位符）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "renderEnvPanel" in html
    assert "__em_jumpToEnv" in html  # claim 点击穿透到 env tab 的入口
    # 旧占位符文本消失
    assert "T2-β 实现" not in html


def test_html_coupling_only_products(cycle_data):
    """Q2: 耦合线只连"产物↔产物"（不画底物节点间连线，不画中间产物节点）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # 化学物节点系统
    assert "em-chem-node" in html
    assert "computeCouplings" in html
    # 只连产物
    assert "role === 'product'" in html or 'role: \'product\'' in html
    # 旧中间产物节点已移除
    assert "em-coupling-product" not in html


def test_html_cross_group_overview_first(cycle_data):
    """Q4: 有分组时跨组概览逻辑存在（renderHypothesisBlock + 跨组评分概览）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "renderHypothesisBlock" in html
    assert "跨组评分概览" in html
    assert "组间 label 差异" in html


def test_html_has_per_group_cycle_support(cycle_data):
    """Q3: HTML 支持 per-group 循环图切换（cycles_by_group + group selector）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "cycles_by_group" in html
    assert "installGroupSelector" in html
    assert "em-group-selector" in html


def test_html_panel_try_catch_isolation(cycle_data):
    """Q1: 各 panel 渲染用 try-catch 隔离 + showPanelError 函数。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "showPanelError" in html
    # 3 处 try { renderXxxPanel } catch
    assert html.count("catch (e) { showPanelError(") >= 3


def test_cycle_to_json_per_group_cycles(cycle_data):
    """Q3: cycle_to_json 接受 per_group_cycles 参数并输出 cycles_by_group。"""
    payload = cycle_to_json(
        cycle_data,
        per_group_cycles={"A": cycle_data, "B": cycle_data},
    )
    assert "cycles_by_group" in payload
    assert set(payload["cycles_by_group"].keys()) == {"A", "B"}
    assert "elements" in payload["cycles_by_group"]["A"]


# ═══════════════════════════════════════════════════════════════
# Q1-Q6 轮二精调（耦合 tooltip / 化学物 hover 关联 / 角色重合 /
# env panel per-group 切换）
# ═══════════════════════════════════════════════════════════════

def test_per_group_cycles_include_full_corr_and_sensitivity(cycle_data):
    """Q6: per_group_cycles payload 补齐 full_corr_matrix + sensitivity。"""
    payload = cycle_to_json(
        cycle_data,
        per_group_cycles={"A": cycle_data},
    )
    g = payload["cycles_by_group"]["A"]
    assert "full_corr_matrix" in g
    assert "sensitivity" in g


def test_html_chemical_hover_link_and_coupling_fallback(cycle_data):
    """Q2 + Q4: chem-links 常驻容器 + 耦合 fallback 到 substrate。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # Q2: chem-links 容器（v1.3 改为常驻 linkSel）
    assert 'id="em-chem-links"' in html
    assert "em-chem-link" in html
    # Q4: computeCouplings 不再强制 product-only
    assert "product || hit.substrate" in html


def test_html_env_panel_supports_per_group_selector(cycle_data):
    """Q6: env panel wrapper 含 per-group 切换（groups.length 判定）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "_renderEnvForSource" in html
    # wrapper 里对 cycles_by_group 的判断
    assert "data.cycles_by_group" in html


def test_html_chem_links_default_persistent(cycle_data):
    """v1.3: 化学物 → pathway 常驻连接线（default opacity 0.18）+ toggle。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # 常驻 linkSel 渲染（低透明度）
    assert "chemLinkData" in html
    assert "em-chem-link" in html
    # 顶部 toggle 按钮
    assert 'id="em-toggle-chem-links"' in html
    # 拖拽时 updateChemLinks 同步
    assert "updateChemLinks" in html


def test_html_has_formatChemName_subsuper(cycle_data):
    """v1.3: formatChemName 格式化化学式下标/上标。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "formatChemName" in html
    # Unicode 下标/上标数字常量数组
    assert "₀" in html and "₉" in html
    assert "⁰" in html and "⁹" in html


def test_html_coupling_uses_custom_tooltip(cycle_data):
    """Q1: 耦合线使用 d3 自定义 showTooltip（不再依赖 <title>）。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    # coupSel 的 mouseenter 事件绑定存在
    assert "coupSel" in html
    # 耦合 tooltip 显示"锚点："含 role 描述
    assert "锚点：" in html and "roleA" in html


def test_html_compare_tab_renamed_and_enhanced(cycle_data):
    """ε.5: compare tab 改名 + 承载者切换高亮 + 排序控件。"""
    html = build_interactive_html(cycle_data).decode("utf-8")
    assert "通路 × 组对比" in html  # 改名
    assert "em-cmp-sort" in html    # 排序 select
    assert "em-cmp-elem" in html    # 元素 filter
    assert "em-cmp-changed-only" in html  # 承载者切换过滤 checkbox
    assert "承载者切换" in html or "承载者跨组变化" in html


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
