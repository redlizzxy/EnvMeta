"""循环图 + 假说评分的独立交互 HTML 导出（S4.5 / T2）。

核心产品价值：审稿人/读者下载**单个 HTML 文件**（~400 KB，含内嵌 D3.js），
双击浏览器打开 → 亲手操作循环图 + 假说评分 + 跨组对比，零安装。

与 matplotlib 静态图的关键差异：
- 4 象限力导向布局（可拖拽节点调整）
- 节点 hover tooltip（MAG / genus / completeness / contribution）
- 点击 gene 高亮整条通路
- null_p 直方图 + 权重敏感度交互
- 跨组 CK/A/B toggle（无需生成 3 个文件）
- URL state（审稿人可粘贴视角给编辑）
- SVG 导出（可编辑矢量）

对比 Krona / Anvi'o / MicrobiomeAnalyst / iTOL：
- 元素循环专用语义（As/N/S/Fe）+ MAG × 通路 × 元素三维模型：独有
- 假说评分闭环（claim + null_p + 权重敏感度）：独有
- 独立离线 HTML：Krona 做过但非领域专用

模块布局：
- `cycle_to_json(cycle_data, hypothesis=None, compare_df=None) -> dict`
- `build_interactive_html(cycle_data, ...) -> bytes`
- `export_html(cycle_data, output_path, ...) -> Path`（便捷封装，写入文件）

依赖：
- `envmeta/geocycle/templates/cycle_interactive.html`（HTML/CSS/JS 模板）
- `envmeta/geocycle/templates/d3.v7.min.js`（inline 嵌入的 D3.js v7.9.0）

使用示例：
    from envmeta.geocycle.html_exporter import build_interactive_html
    html_bytes = build_interactive_html(cycle_data, hypothesis=hyp_score)
    Path("demo.html").write_bytes(html_bytes)
"""
from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    import pandas as pd

    from envmeta.geocycle.hypothesis import HypothesisScore
    from envmeta.geocycle.model import CycleData


# ═══════════════════════════════════════════════════════════════
# 模板资源定位
# ═══════════════════════════════════════════════════════════════

_TEMPLATE_DIR = Path(__file__).parent / "templates"
_HTML_TEMPLATE_PATH = _TEMPLATE_DIR / "cycle_interactive.html"
_D3_JS_PATH = _TEMPLATE_DIR / "d3.v7.min.js"

# 占位符（严格 {{ABC}} 格式，不用 Python format 避免 D3 代码里的 { 冲突）
_PLACEHOLDER_D3 = "{{D3_JS_INLINE}}"
_PLACEHOLDER_DATA = "{{CYCLE_DATA_JSON}}"
_PLACEHOLDER_META = "{{META_HTML}}"


# ═══════════════════════════════════════════════════════════════
# 数据转换
# ═══════════════════════════════════════════════════════════════


def cycle_to_json(
    cycle_data: "CycleData",
    *,
    hypothesis: Optional["HypothesisScore"] = None,
    compare_df: Optional["pd.DataFrame"] = None,
    hypothesis_by_group: Optional[dict[str, "HypothesisScore"]] = None,
    per_group_cycles: Optional[dict[str, "CycleData"]] = None,
) -> dict:
    """把 CycleData (+ 可选 HypothesisScore / 跨组对比 DataFrame) 扁平为
    JSON-serializable dict，供前端 D3 消费。

    字段约定（前端依赖此结构，修改前评估兼容性）：
    - `version`: str — 导出格式版本（"1.0"）
    - `generated_at`: str — ISO8601 UTC 时间戳
    - `envmeta_version`: str — EnvMeta 版本
    - `elements`: list[ElementCycle asdict]（嵌套 pathways → contributors → genes）
    - `env_correlations`: list[EnvCorrelation asdict]（过滤后）
    - `full_corr_matrix`: list[EnvCorrelation asdict]（不过滤）
    - `sensitivity`: list[SensitivityRow asdict]（S1 去偏产物）
    - `params`: dict — 推断参数
    - `meta`: dict — 样本数 / MAG 数 / 数据集签名
    - `hypothesis`: dict | None — HypothesisScore asdict（含 null_p / 权重敏感度 / veto）
    - `hypothesis_by_group`: dict[str, dict] | None — 跨组假说评分
    - `compare_groups`: list[dict] | None — DataFrame.to_dict(orient='records')

    参数
    -----
    cycle_data: 必需
    hypothesis: 可选 — 单组假说评分
    hypothesis_by_group: 可选 — 跨组假说评分（score_by_groups 返回）
    compare_df: 可选 — cycle_compare.compare_groups 返回的 DataFrame
    """
    # 延迟导入避免 envmeta.__init__ 慢 startup
    from envmeta import __version__ as _env_version
    from envmeta.geocycle.knowledge_base import couplings as _kb_couplings

    payload: dict = {
        "version": "1.0",
        "generated_at": _dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "envmeta_version": _env_version,
        "elements": [asdict(el) for el in cycle_data.elements],
        "env_correlations": [asdict(ec) for ec in cycle_data.env_correlations],
        "full_corr_matrix": [asdict(ec) for ec in cycle_data.full_corr_matrix],
        "sensitivity": [asdict(sr) for sr in cycle_data.sensitivity],
        "couplings": list(_kb_couplings()),  # T2-β: 跨元素化学物耦合（KB 定义）
        "params": dict(cycle_data.params),
        "meta": dict(cycle_data.meta),
    }

    if hypothesis is not None:
        payload["hypothesis"] = asdict(hypothesis)

    if hypothesis_by_group is not None:
        # 支持两种类型：
        # 1) dict[str, HypothesisScore] — 直接 asdict
        # 2) DataFrame (score_by_groups 返回的 summary) — 转 dict[group, row_dict]
        if isinstance(hypothesis_by_group, dict):
            payload["hypothesis_by_group"] = {
                g: asdict(s) for g, s in hypothesis_by_group.items()
            }
        else:
            # 假定是 DataFrame-like
            try:
                _df = hypothesis_by_group
                records = json.loads(_df.to_json(orient="records", default_handler=str))
                payload["hypothesis_by_group_summary"] = records  # 简化版：表格化
            except Exception:  # noqa: BLE001
                # 非 dict / 非 DataFrame 就丢掉
                pass

    if compare_df is not None:
        # DataFrame → records（用 to_json round-trip 稳定处理 NaN → None）
        # `to_dict(orient='records')` 不会自动把 NaN 转成 None；
        # 走 to_json (pandas 输出 null) → json.loads 一趟是最干净的做法
        payload["compare_groups"] = json.loads(
            compare_df.to_json(orient="records", default_handler=str)
        )

    # Q3: per-group cycle 数据（HTML 可切换查看 CK / A / B 的循环图）
    # 每组是独立的 CycleData（来自 infer(group_filter=g)）
    # Q6: 补齐 full_corr_matrix + sensitivity，供 env panel 按组切换
    if per_group_cycles:
        payload["cycles_by_group"] = {
            g: {
                "elements": [asdict(el) for el in cd.elements],
                "env_correlations": [asdict(ec) for ec in cd.env_correlations],
                "full_corr_matrix": [asdict(ec) for ec in cd.full_corr_matrix],
                "sensitivity": [asdict(sr) for sr in cd.sensitivity],
                "meta": dict(cd.meta),
            }
            for g, cd in per_group_cycles.items()
        }

    return payload


# ═══════════════════════════════════════════════════════════════
# HTML 组装
# ═══════════════════════════════════════════════════════════════


def build_interactive_html(
    cycle_data: "CycleData",
    *,
    hypothesis: Optional["HypothesisScore"] = None,
    compare_df: Optional["pd.DataFrame"] = None,
    hypothesis_by_group: Optional[dict[str, "HypothesisScore"]] = None,
    per_group_cycles: Optional[dict[str, "CycleData"]] = None,
    title: str = "EnvMeta — Interactive Biogeochemical Cycle",
) -> bytes:
    """把 CycleData 渲染成独立可交互 HTML（bytes）。

    返回 bytes 便于 Streamlit st.download_button 直接使用；写文件请用
    `export_html` 便捷封装。
    """
    if not _HTML_TEMPLATE_PATH.exists():
        raise FileNotFoundError(
            f"HTML 模板缺失：{_HTML_TEMPLATE_PATH}。"
            "请确认 envmeta/geocycle/templates/cycle_interactive.html 存在。"
        )
    if not _D3_JS_PATH.exists():
        raise FileNotFoundError(
            f"D3.js 资源缺失：{_D3_JS_PATH}。"
            "请从 https://d3js.org/d3.v7.min.js 下载到该路径。"
        )

    template = _HTML_TEMPLATE_PATH.read_text(encoding="utf-8")
    d3_js = _D3_JS_PATH.read_text(encoding="utf-8")

    payload = cycle_to_json(
        cycle_data,
        hypothesis=hypothesis,
        compare_df=compare_df,
        hypothesis_by_group=hypothesis_by_group,
        per_group_cycles=per_group_cycles,
    )
    # 前端消费的 JSON 字符串。用 ensure_ascii=False 保留中文。
    # 不用 indent — 压缩体积 ~30%
    data_json = json.dumps(payload, ensure_ascii=False, default=str)

    # Meta 区块（HTML 顶部审计信息）
    meta_html = _build_meta_html(payload, title=title)

    html = template.replace(_PLACEHOLDER_D3, d3_js)
    html = html.replace(_PLACEHOLDER_DATA, data_json)
    html = html.replace(_PLACEHOLDER_META, meta_html)

    return html.encode("utf-8")


def export_html(
    cycle_data: "CycleData",
    output_path: str | Path,
    **kwargs,
) -> Path:
    """便捷封装：渲染 + 写文件。返回最终路径。"""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(build_interactive_html(cycle_data, **kwargs))
    return out


# ═══════════════════════════════════════════════════════════════
# 辅助
# ═══════════════════════════════════════════════════════════════


def _build_meta_html(payload: dict, *, title: str) -> str:
    """生成 HTML 顶部的元信息块（可审计 footer）。"""
    meta = payload.get("meta", {})
    params = payload.get("params", {})
    n_elements = len(payload.get("elements", []))
    n_pathways = sum(len(el.get("pathways", [])) for el in payload.get("elements", []))
    has_hyp = payload.get("hypothesis") is not None
    has_groups = payload.get("hypothesis_by_group") is not None
    has_compare = payload.get("compare_groups") is not None

    lines = [
        f'<div class="em-meta-block">',
        f'  <h1 class="em-title">{_escape_html(title)}</h1>',
        f'  <div class="em-badges">',
        f'    <span class="em-badge">EnvMeta v{payload.get("envmeta_version", "?")}</span>',
        f'    <span class="em-badge">{n_elements} elements · {n_pathways} pathways</span>',
    ]
    if has_hyp:
        lines.append('    <span class="em-badge em-badge-hyp">✓ hypothesis</span>')
    if has_groups:
        lines.append('    <span class="em-badge em-badge-hyp">✓ hypothesis × groups</span>')
    if has_compare:
        lines.append('    <span class="em-badge em-badge-cmp">✓ cross-group compare</span>')
    lines.append(f'  </div>')
    lines.append(f'  <div class="em-audit">')
    lines.append(f'    生成时间 {_escape_html(payload.get("generated_at", "?"))}  ·  ')
    n_mags = meta.get("n_mags", "?")
    n_samples = meta.get("n_samples", "?")
    lines.append(
        f'    {n_mags} MAG × {n_samples} samples  ·  '
        f'completeness ≥ {params.get("completeness_threshold", "?")}%'
    )
    lines.append(f'  </div>')
    # T2-ε.1: 节点筛选标准明示（用户一眼看到 HTML 展示的是怎么选出来的）
    comp_thr = params.get("completeness_threshold", "?")
    top_n = params.get("top_contributors_per_pathway", params.get("top_n_contributors", "all"))
    ranking = params.get("contributor_ranking", "abundance × completeness")
    min_abund = params.get("min_abundance_mean", None)
    parts = [
        f"completeness ≥ <b>{comp_thr}%</b>",
        f"每通路 Top <b>{top_n}</b> contributors" if top_n != "all" else "展示所有活跃 MAG",
        f"排序 = <b>{ranking}</b>",
    ]
    if min_abund is not None:
        parts.append(f"abundance ≥ <b>{min_abund}</b>")
    lines.append(
        f'  <div class="em-audit" style="margin-top:4px;'
        f'padding-top:6px;border-top:1px dashed #c8d4e3">'
    )
    lines.append(
        f'    <span style="color:#1a365d;font-weight:600">🔍 节点筛选标准：</span>'
        f' {" · ".join(parts)}'
    )
    lines.append(f'  </div>')
    lines.append(f'</div>')
    return "\n".join(lines)


def _escape_html(s: object) -> str:
    """最小 HTML 转义（dataclass 字段里可能有 <>）。"""
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )
