"""
EnvMeta — 环境微生物宏基因组可视化分析平台
Streamlit 主入口：streamlit run app.py
"""

from pathlib import Path

import pandas as pd
import streamlit as st

from envmeta import __version__
from envmeta.analysis import (
    cycle_diagram as cycle_mod,
    gene_heatmap, gene_profile as gene_profile_mod,
    lefse as lefse_mod,
    mag_heatmap as mag_heatmap_mod,
    mag_quality as mag_quality_mod,
    network as network_mod,
    pathway as pathway_mod,
    pcoa, rda as rda_mod, stackplot,
)
from envmeta.tools.gephi_prep import prepare_gephi_csv, validate_gephi_format
from envmeta.export.code_generator import generate as generate_code
from envmeta.export.figure_export import export_to_bytes
from envmeta.file_manager.detector import FileType, detect, read_table
from envmeta.help import (
    ANALYSIS_INPUTS, FILE_TO_ANALYSIS, INTERPRETATIONS, NAVIGATOR,
)
from envmeta.help.file_analysis_map import analyses_ready
from envmeta.params.common import render_figure_size, render_font_controls


def _vector_downloads(figure, basename: str, key_prefix: str) -> None:
    """渲染「SVG（SCI 投稿）」+「TIFF 600dpi（毕业论文 / 期刊位图）」两键。

    与已有的 PNG / PDF 按钮并列出现，让用户按期刊要求一键拿矢量 / 高分位图。
    """
    cs = st.columns(2)
    with cs[0]:
        st.download_button(
            "⬇️ SVG（SCI 投稿）",
            data=export_to_bytes(figure, "svg"),
            file_name=f"{basename}.svg",
            mime="image/svg+xml",
            key=f"{key_prefix}_svg",
            help="纯矢量图形，无损放大；SCI 投稿多数期刊首选",
        )
    with cs[1]:
        st.download_button(
            "⬇️ TIFF 600dpi（毕业论文）",
            data=export_to_bytes(figure, "tiff"),
            file_name=f"{basename}.tiff",
            mime="image/tiff",
            key=f"{key_prefix}_tiff",
            help="LZW 压缩 600dpi TIFF，符合多数期刊位图最低要求",
        )


def _reproduce_button(analysis_id: str, file_paths: dict[str, str],
                      params: dict, key: str, output_base: str | None = None):
    """渲染「下载 .py 复现脚本」按钮。"""
    try:
        src = generate_code(analysis_id, file_paths, params,
                            output_base=output_base or analysis_id)
    except Exception as e:
        st.caption(f"⚠️ 无法生成复现脚本：{e}")
        return
    st.download_button(
        "⬇️ 复现脚本（.py）",
        data=src.encode("utf-8"),
        file_name=f"{analysis_id}_reproduce.py",
        mime="text/x-python",
        key=key,
        help="独立可运行的 Python 脚本，复现当前参数下的图表",
    )


def _interpret_hyp_score(hyp_last) -> tuple[str, str]:
    """根据 (label, overall, null_p, weight_robust, vetoed) 的组合，
    返回 (一句话解读, alert_type)。alert_type ∈ {error, warning, success, info}。

    10 档判断覆盖"最强支持 / 通过但不特异 / 权重敏感 / 边界 / 中等 / 弱 / 未通过"等
    常见情景；"退化式 N/A"(null_p=None AND overall≥0.95 ≈ 全 satisfied) 视为
    可信 signal（兜底有 weight_robust 把关）。
    """
    overall = hyp_last.overall_score
    label = hyp_last.label
    p = hyp_last.null_p
    robust = hyp_last.weight_robust
    vetoed = bool(hyp_last.veto_reasons)

    if vetoed:
        return (
            "🚫 **假说未通过**：必要前提（required=true 的 claim）未被数据支持；"
            "无论 overall 多高都被硬否决。",
            "error",
        )

    specific = (p is not None and p < 0.05)
    degenerate_ok = (p is None and overall >= 0.95)

    if label == "strong":
        if (specific or degenerate_ok) and robust:
            return (
                "⭐ **最强支持**：label=strong，权重设计与数据支持一致"
                + ("（全 claim 通过，排列退化）" if degenerate_ok else "")
                + "，对 ±20% 权重扰动稳健。",
                "success",
            )
        if robust and p is not None and p >= 0.20:
            return (
                f"⚠️ **strong 但不特异**：label 达 strong，但 null_p={p:.2f} 偏高"
                "，说明通过率够高即可达此分，权重设计未被数据特异支持。"
                "建议关注哪些 claim 是真正的核心机制。",
                "warning",
            )
        if not robust:
            return (
                "⚠️ **strong 但权重敏感**：±20% 权重扰动下 label 翻转，"
                "结论不稳，建议审视权重分配依据。",
                "warning",
            )
        if p is not None and 0.05 <= p < 0.20:
            return (
                f"🟢 **strong 边界**：null_p={p:.2f}（0.05-0.20），"
                "权重设计与数据中度一致，建议加入更多独立证据。",
                "info",
            )
        # fallback（应极少命中）
        return (
            "🟢 **strong**：overall 达阈值；具体指标见下方。",
            "info",
        )

    if label == "suggestive":
        if specific or degenerate_ok:
            return (
                "🟡 **中等支持**：未达 strong 但权重设计与数据一致，"
                "可作为支持性证据叙述（非主论证）。",
                "info",
            )
        return (
            "🟡 **中等支持**：证据部分吻合，未达 strong 阈值。",
            "info",
        )

    if label == "weak":
        return (
            "🟠 **证据薄弱**：overall > 0 但低于 suggestive 阈值；"
            "数据仅零星支持假说。",
            "warning",
        )

    # insufficient（未 veto）
    return (
        "🔴 **证据不足**：overall ≈ 0，或无可评分 claim（全 skipped）；"
        "数据未提供任何支持。",
        "error",
    )


# ── 页面配置 ──────────────────────────────────────────────
st.set_page_config(
    page_title="EnvMeta",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

if "files" not in st.session_state:
    st.session_state.files = {}

# ── 侧边栏导航 ────────────────────────────────────────────
st.sidebar.title("EnvMeta")
st.sidebar.caption(f"v{__version__} · 环境微生物宏基因组可视化分析平台")

# 跳转支持：其他页的按钮改 session_state.page_radio（widget key）即可控制 radio
_PAGES = [
    "首页",
    "数据准备指南",
    "文件管理",
    "图表选择向导",
    "Reads-based 分析",
    "MAG-based 分析",
    "生物地球化学循环图",
    "导出中心",
]
if "page_radio" not in st.session_state:
    st.session_state.page_radio = "首页"

page = st.sidebar.radio("功能模块", _PAGES, key="page_radio")

TYPE_BADGES = {
    FileType.METADATA: "🗂️ metadata",
    FileType.ABUNDANCE_WIDE: "📊 abundance (wide)",
    FileType.DISTANCE_MATRIX: "📐 distance matrix",
    FileType.ALPHA_DIVERSITY: "📈 alpha diversity",
    FileType.CHECKM_QUALITY: "🧬 CheckM quality",
    FileType.ENV_FACTORS: "🌡️ env factors",
    FileType.KO_ABUNDANCE_WIDE: "🧪 KO abundance",
    FileType.KO_ANNOTATION_LONG: "🔬 KO 注释（长表）",
    FileType.KEYSTONE_SPECIES: "⭐ keystone species",
    FileType.MAG_TAXONOMY: "🌳 MAG taxonomy",
    FileType.GEPHI_NODES: "🔗 Gephi nodes",
    FileType.GEPHI_EDGES: "🔗 Gephi edges",
    FileType.UNKNOWN: "❓ unknown",
}
TYPE_OPTIONS = [ft.value for ft in FileType]


def _files_of(*types: FileType) -> dict[str, dict]:
    return {n: info for n, info in st.session_state.files.items() if info["type"] in types}


# ══════════════════════════════════════════════════════════
# S8-ux 新手落地包辅助
# ══════════════════════════════════════════════════════════

def _render_interpretation_expander(analysis_id: str, *, key: str | None = None) -> None:
    """渲染「如何解读」expander。数据源：envmeta/help/interpretations.py。

    用法：在每个分析页分支的**开头**调用一次，让用户生成图之前就看到解读指南。
    """
    content = INTERPRETATIONS.get(analysis_id)
    if not content:
        return
    with st.expander(f"📖 {content['title']}", expanded=False):
        st.markdown(f"**这张图回答什么？**  \n{content['what_it_shows']}")
        st.markdown("**怎么看？**")
        for bullet in content["how_to_read"]:
            st.markdown(f"- {bullet}")
        st.success(f"✅ **正面证据**：{content['good_signal']}")
        st.warning(f"⚠️ **常见误判**：{content['warning']}")
        st.info(f"ℹ️ **方法学局限**：{content['caveats']}")


def _jump_callback(analysis_id: str) -> None:
    """on_click 回调 — 在 rerun 前改 widget session_state，避开
    "cannot be modified after widget instantiated" 限制。"""
    spec = ANALYSIS_INPUTS.get(analysis_id)
    if not spec:
        return
    st.session_state.page_radio = spec["page"]
    atype = spec.get("analysis_type")
    if atype:
        if spec["page"] == "Reads-based 分析":
            st.session_state.reads_analysis_type = atype
        elif spec["page"] == "MAG-based 分析":
            st.session_state.mag_analysis_type = atype


def _goto_page_callback(target: str) -> None:
    """on_click 回调 — 切到指定 sidebar page。"""
    st.session_state.page_radio = target


# ── 导出中心注册表（T1 / S8-ui）─────────────────────────────
# 每个分析生成结果后调用 `_register_export`，导出中心从中列出全部产物
# 批量下载 / .py 复现脚本 / 单独下载 等操作都在这里集中

def _register_export(
    analysis_id: str, *, result, file_paths: dict[str, str],
    params: dict, output_base: str | None = None,
) -> None:
    """把分析结果登记到 session_state._export_registry 供导出中心统一访问。

    参数
    -----
    analysis_id: 与 ANALYSIS_INPUTS 一致的键（stackplot / pcoa / ...）
    result: AnalysisResult（需有 .figure + .stats）
    file_paths: 传给 code_generator 的输入文件映射
    params: 分析参数字典
    output_base: 下载文件基础名（默认 analysis_id）
    """
    reg = st.session_state.setdefault("_export_registry", {})
    reg[analysis_id] = {
        "result": result,
        "file_paths": file_paths,
        "params": params,
        "output_base": output_base or analysis_id,
    }


def _jump_to_analysis(analysis_id: str, *, key: str) -> None:
    """渲染「跳转到分析」按钮。通过 on_click 回调安全修改 widget key。"""
    spec = ANALYSIS_INPUTS.get(analysis_id)
    if not spec:
        return
    st.button(
        f"🚀 去跑「{spec['name']}」",
        key=key,
        on_click=_jump_callback,
        args=(analysis_id,),
    )


def _render_file_reverse_index(file_type: FileType, *, key_prefix: str) -> None:
    """文件管理页 — 每个文件卡片下方的「可跑什么分析」区块。"""
    entries = FILE_TO_ANALYSIS.get(file_type, [])
    if not entries:
        st.caption("这个文件类型暂无对应分析。")
        return
    st.markdown("**🎯 这个文件可以跑什么分析**")
    # 已上传的 FileType 集合，用于高亮"全部 required 已齐全"
    have_types = {info["type"] for info in st.session_state.files.values()}
    ready = set(analyses_ready(have_types))
    for aid, name, role in entries:
        role_badge = "🔴 required" if role == "required" else "🟢 optional"
        ready_badge = " ✅ 已齐全" if aid in ready else ""
        c1, c2 = st.columns([4, 1])
        with c1:
            st.markdown(f"- **{name}** · {role_badge}{ready_badge}")
        with c2:
            _jump_to_analysis(aid, key=f"{key_prefix}_jump_{aid}")


def _load_sample_data() -> tuple[int, int]:
    """加载 tests/sample_data/ 里的全部文件到 session_state。

    返回 (成功加载数, 跳过/失败数)。复用现有 detect() 逻辑。
    """
    sample_dir = Path(__file__).parent / "tests" / "sample_data"
    if not sample_dir.exists():
        st.error(f"样例数据目录不存在：{sample_dir}")
        return (0, 0)
    loaded = 0
    skipped = 0
    for fpath in sorted(sample_dir.iterdir()):
        if fpath.suffix.lower() not in (".txt", ".tsv", ".csv", ".spf"):
            continue
        if fpath.name in st.session_state.files:
            skipped += 1
            continue
        try:
            df, enc, sep = read_table(fpath)
            result = detect(fpath, filename=fpath.name)
            st.session_state.files[fpath.name] = {
                "df": df, "type": result.file_type, "result": result,
            }
            loaded += 1
        except Exception:
            skipped += 1
    return loaded, skipped


# ══════════════════════════════════════════════════════════
# MAG-based 共享 UI（4 张图的侧边栏 Layer 1-2 一致）
# ══════════════════════════════════════════════════════════

def _first_of(*ftypes: FileType) -> str | None:
    for n, info in st.session_state.files.items():
        if info["type"] in ftypes:
            return n
    return None


def _abundance_candidates() -> list[tuple[float, str]]:
    """返回所有 ABUNDANCE_WIDE 文件 [(confidence, name)]，按 detector 置信度降序。"""
    out: list[tuple[float, str]] = []
    for n, info in st.session_state.files.items():
        if info["type"] != FileType.ABUNDANCE_WIDE:
            continue
        conf = getattr(info.get("result"), "confidence", 0.0) or 0.0
        out.append((conf, n))
    out.sort(key=lambda t: -t[0])
    return out


def _first_mag_abundance() -> str | None:
    """优先选 MAG 级丰度表（detector conf=0.95）。

    修复 2026-04-19：Genus.txt / Phylum.txt / Species.txt 也会被识别为
    ABUNDANCE_WIDE，`_first_of` 按插入顺序选到 Genus.txt → MAG 分析里
    abundance_mean 全部为 0 → Top-N 退化为 tiebreaker 顺序（bug）。
    """
    mag_headers = {"genome", "mag", "bin", "mag_id", "genome_id"}
    candidates = _abundance_candidates()
    if not candidates:
        return None
    top_conf, top_name = candidates[0]
    if top_conf >= 0.92:
        return top_name
    for _, n in candidates:
        df = st.session_state.files[n]["df"]
        if df.shape[1] > 0 and df.columns[0].strip().lower() in mag_headers:
            return n
    return candidates[0][1]


def _first_taxon_abundance() -> str | None:
    """优先选 TAXON 级丰度表（detector conf=0.88，首列 = Taxonomy/OTU）。

    修复 2026-04-19：堆叠图 / PCoA / RDA / LEfSe 等"物种组成"语义分析
    需要 TAXON 级（Genus/Phylum/Species），而非 MAG 级（abundance.tsv）。
    跨平台一致：不依赖 session_state.files 插入顺序（Windows/Linux 各异），
    改用 detector 置信度（0.88 = TAXON，0.95 = MAG）判定。
    """
    taxon_headers = {"taxonomy", "#otu id", "#otuid", "taxon"}
    candidates = _abundance_candidates()
    if not candidates:
        return None
    # 优先选 TAXON 级（conf ≈ 0.88，即 conf < 0.92）
    taxon_like = [(c, n) for c, n in candidates if c < 0.92]
    if taxon_like:
        # 同为 TAXON 级时，按文件名字母序稳定选择（Windows/Linux 一致）
        taxon_like.sort(key=lambda t: t[1].lower())
        return taxon_like[0][1]
    # 回退：首列匹配 taxonomy/otu
    for _, n in candidates:
        df = st.session_state.files[n]["df"]
        if df.shape[1] > 0 and df.columns[0].strip().lower() in taxon_headers:
            return n
    return candidates[0][1]


def _idx_or_default(options: list[str], name: str | None) -> int:
    return options.index(name) if (name and name in options) else 0


def render_mag_layer1_filter(prefix: str,
                              default_max_mags: int = 0) -> dict:
    """渲染 MAG 图统一 Layer 1 过滤参数。所有 4 张 MAG 图用同一套 widget。"""
    st.markdown("**子集过滤**")
    filter_mode = st.selectbox(
        "MAG 子集",
        ["top_plus_keystone", "top_n", "keystone_only", "all"],
        index=0, key=f"{prefix}_fm",
        help=(
            "• top_plus_keystone（默认）: Top-N 丰度 ∪ 全部 keystone\n"
            "• top_n: 只留按 Top-N 打分选出的 N 条\n"
            "• keystone_only: 只留 keystone 物种\n"
            "• all: 不过滤，显全部 MAG（再由 max_mags 硬截断）"
        ),
    )
    top_n_by = st.selectbox(
        "Top-N 打分依据",
        ["mean", "sum", "variance"],
        index=0, key=f"{prefix}_tnby",
        help=(
            "• mean: 平均丰度（默认）\n"
            "• sum: 累计丰度\n"
            "• variance: 方差（组间差异大的 MAG）"
        ),
    )
    top_n_count = st.slider("Top-N 取几条", 5, 100, 30,
                             step=5, key=f"{prefix}_tnc")
    max_mags = st.slider("max_mags 最终展示上限（0=不截）", 0, 200,
                          default_max_mags, step=10, key=f"{prefix}_mm",
                          help=(
                              "顺序：filter_mode 过滤 → row_order 排序 → "
                              "max_mags 截断 → 画图。若过滤 + 排序后剩的 MAG "
                              "数 ≤ max_mags，此项不生效；超过才截前 N。"
                              "设 0 = 不截断。"
                          ))
    return {
        "filter_mode": filter_mode,
        "top_n_by": top_n_by,
        "top_n_count": top_n_count,
        "max_mags": max_mags,
    }


def render_mag_layer2_visual(prefix: str,
                              default_width: int, default_height: int,
                              show_phylum_bar: bool = True) -> dict:
    """渲染 MAG 图统一 Layer 2 视觉参数。"""
    st.markdown("**视觉**")
    hi_ks = st.checkbox("★ keystone 标记", True, key=f"{prefix}_hiks")
    show_phy = True
    if show_phylum_bar:
        show_phy = st.checkbox("门彩条（左侧）", True, key=f"{prefix}_phy")
    show_leg = st.checkbox("门图例（右侧）", True, key=f"{prefix}_leg")
    size = render_figure_size({"width_mm": default_width,
                               "height_mm": default_height},
                              prefix=prefix)
    return {
        "highlight_keystones": hi_ks,
        "show_phylum_bar": show_phy,
        "show_phylum_legend": show_leg,
        **size,
    }


def render_mag_layer3_order(prefix: str) -> dict:
    """渲染 MAG 热图统一 Layer 3 行排序参数（mag_quality 散点图不用）。"""
    st.markdown("**行排序**")
    row_order = st.selectbox(
        "排序方式",
        ["phylum_cluster", "metric_desc", "abundance"],
        index=0, key=f"{prefix}_ro",
        help=(
            "• phylum_cluster（默认）: 按门分组 → 门内层次聚类\n"
            "• metric_desc: 按该图主 metric（完整度/基因数/丰度分）降序\n"
            "• abundance: 按丰度均值降序"
        ),
    )
    linkage_method = st.selectbox(
        "聚类方法",
        ["average", "ward", "complete"],
        index=0, key=f"{prefix}_lm",
    )
    return {"row_order": row_order, "linkage_method": linkage_method}


def _mag_file_selectors(prefix: str,
                         *, ko_type: FileType | None = None,
                         require_abundance: bool = False) -> dict[str, str]:
    """渲染 MAG 图统一文件上传下拉（4 种标准输入，按 FileType 自动选默认）。

    返回：{"_ko": name | None, "abundance": name | None, "taxonomy": ...,
           "keystone": ..., "metadata": ...}
    """
    file_names = list(st.session_state.files.keys())
    if not file_names:
        st.warning("请先在「文件管理」上传数据文件。")
        st.stop()

    selected: dict[str, str | None] = {}

    # KO 注释表（仅 pathway / gene_profile 需要）
    if ko_type is not None:
        ko_default = _first_of(ko_type)
        ko_name = st.selectbox(
            "KO 注释表（MAG + KEGG_ko 长表）",
            file_names, key=f"{prefix}_ko",
            index=_idx_or_default(file_names, ko_default),
        )
        selected["_ko"] = ko_name
    used = [selected.get("_ko")] if ko_type else []

    # MAG 丰度表（优先 MAG 级 abundance.tsv，避开 Genus.txt / Phylum.txt）
    ab_default = _first_mag_abundance()
    ab_label = ("MAG 丰度表（必需）" if require_abundance
                else "MAG 丰度表（可选，用于 Top-N 打分）")
    if require_abundance:
        ab_options = [n for n in file_names if n not in used]
        ab_name = st.selectbox(ab_label, ab_options, key=f"{prefix}_ab",
                               index=_idx_or_default(ab_options, ab_default))
    else:
        ab_options = ["（无）"] + [n for n in file_names if n not in used]
        ab_name = st.selectbox(ab_label, ab_options, key=f"{prefix}_ab",
                               index=_idx_or_default(ab_options, ab_default))
    selected["abundance"] = ab_name if ab_name != "（无）" else None
    used.append(ab_name if ab_name != "（无）" else None)

    # Taxonomy
    t_default = _first_of(FileType.MAG_TAXONOMY)
    t_options = ["（无）"] + [n for n in file_names if n not in used]
    t_name = st.selectbox(
        "GTDB 分类表（可选，用于 Genus 标签 + 门彩条）",
        t_options, key=f"{prefix}_tax",
        index=_idx_or_default(t_options, t_default),
    )
    selected["taxonomy"] = t_name if t_name != "（无）" else None
    used.append(t_name if t_name != "（无）" else None)

    # Keystone
    k_default = _first_of(FileType.KEYSTONE_SPECIES)
    k_options = ["（无）"] + [n for n in file_names if n not in used]
    k_name = st.selectbox(
        "Keystone 物种列表（可选）",
        k_options, key=f"{prefix}_ks",
        index=_idx_or_default(k_options, k_default),
    )
    selected["keystone"] = k_name if k_name != "（无）" else None
    used.append(k_name if k_name != "（无）" else None)

    # Metadata
    m_default = _first_of(FileType.METADATA)
    m_options = ["（无）"] + [n for n in file_names if n not in used]
    m_name = st.selectbox(
        "样本分组表（可选，提供组彩条 / 样本排序）",
        m_options, key=f"{prefix}_md",
        index=_idx_or_default(m_options, m_default),
    )
    selected["metadata"] = m_name if m_name != "（无）" else None

    return selected


def _load_mag_df(name: str | None):
    """按 name 取 DataFrame；处理 taxonomy 无 header 的特殊格式。"""
    if name is None:
        return None
    t = st.session_state.files[name]["df"]
    if t.shape[1] == 2:
        cols_lower = [c.lower() for c in t.columns]
        if ("classification" not in cols_lower
                and "taxonomy" not in cols_lower
                and not any("mag" in c for c in cols_lower)):
            t = t.copy()
            t.columns = ["MAG", "Taxonomy"]
    return t


# ══════════════════════════════════════════════════════════
# 首页
# ══════════════════════════════════════════════════════════
if page == "首页":
    st.title("EnvMeta")
    st.subheader("环境微生物宏基因组可视化分析平台 · v0.8")
    st.markdown("---")
    st.markdown(
        """
        **解决环境微生物博士生的核心痛点 —— 测序公司给了一堆表格，
        不知道哪个文件能做什么分析。** 文件一键识别 + 14 张发表级图表 +
        元素循环图自动推断 + 假说评分器 + 独立交互 HTML 导出。
        **全部开源免费、离线可用。**

        #### 🚀 三步上手
        1. 首页点「📦 加载砷渣修复示例数据」（或自己上传文件到「文件管理」）
        2. 走「🧭 图表选择向导」按研究问题找分析，或直接进 Reads/MAG/循环图分析
        3. 调参出图 → 导出中心批量下载 PNG/PDF/SVG/TIFF + `.py` 复现脚本

        #### ✨ 核心卖点
        | 能力 | 说明 |
        |---|---|
        | 🔄 **元素循环图自动推断** | 从 KO 注释 + 环境因子自动推断 4 元素 × 18 通路（业界独有） |
        | 🧪 **假说评分 YAML 评分器** | 5 类 claim + 置换 null_p + 权重敏感度 + 9 档解读（业界独有） |
        | 🌐 **独立交互 HTML 导出** | 400 KB 单文件 D3.js 嵌入，审稿人浏览器直接操作（SI 杀手锏） |
        | 📦 **Fork Bundle** | 打包 KB + YAML + config 为 zip，论文-工具绑定复现 |
        | ⚡ **跨元素化学物耦合** | As↔H₂S→As₂S₃ / Fe↔S 等虚线连接（业界独有） |
        | 💡 **新手落地包** | 数据准备指南 + 图表向导 + 每图「如何解读」+ 一键样例数据 |

        #### 📊 功能矩阵
        | 模块 | 支持内容 |
        |------|---------|
        | 📁 文件识别 | 11 类（abundance 分 MAG/Taxon 两级识别 / KO 宽长 / CheckM / GTDB / env / Gephi ...） |
        | 📊 Reads-based（7 图） | 堆叠图 / PCoA / α 多样性 / RDA / LEfSe / 基因热图 / log2FC |
        | 🧬 MAG-based（5 图） | MAG 质量 / 丰度热图 / 通路完整度 / 基因谱 / 共现网络 Gephi-prep |
        | 🔄 循环图 | 4 元素推断 + 假说评分 + 跨组对比 + Fork Bundle + 独立 HTML |
        | 💾 导出 | PNG / PDF / SVG / TIFF 600dpi / TSV / `.py` 脚本，批量 ZIP |

        > 📖 **在线版**：<https://envmeta-3xjhcu7lv2gkj4pjtk8gsb.streamlit.app/>
        > 📖 **源码**：<https://github.com/redlizzxy/EnvMeta>
        > 📖 **数据准备**：左侧「数据准备指南」— 覆盖 CoverM / HUMAnN3 / eggNOG / DRAM / GTDB-Tk / CheckM2 等 11 种上游工具
        """
    )
    st.markdown("---")

    # ── S8-ux 新手落地包入口 ────────────────────────────────
    c_home1, c_home2 = st.columns([1, 1])
    with c_home1:
        with st.expander("🚀 快速体验（一键加载示例数据）", expanded=True):
            st.markdown(
                "没有数据？点击下方按钮加载砷渣-钢渣修复示例数据"
                "（168 MAG / 10 样本 / 3 组），即可跑通全部 14 个分析。"
            )
            if st.button("📦 加载砷渣修复示例数据", key="home_load_sample"):
                loaded, skipped = _load_sample_data()
                if loaded:
                    st.success(
                        f"已加载 {loaded} 个样例文件（跳过 {skipped} 个）。"
                        "切换到左侧「文件管理」或「图表选择向导」继续。"
                    )
                    st.toast("示例数据就绪 ✅", icon="📦")
                else:
                    st.warning(
                        f"没有加载任何新文件（已跳过 {skipped} 个，可能已加载过）。"
                    )
            if st.session_state.files:
                st.caption(f"📂 当前已加载 {len(st.session_state.files)} 个文件。")
    with c_home2:
        with st.expander("🧭 不知道用哪张图？", expanded=True):
            st.markdown(
                "「**图表选择向导**」按研究问题推荐分析（组间差异 / 多样性 / 网络 / "
                "元素循环 / 环境因子关系 等 8 大类）。"
            )
            st.button(
                "📊 打开图表选择向导", key="home_open_nav",
                on_click=_goto_page_callback, args=("图表选择向导",),
            )
            st.button(
                "📚 查看数据准备指南", key="home_open_guide",
                on_click=_goto_page_callback, args=("数据准备指南",),
            )

# ══════════════════════════════════════════════════════════
# 数据准备指南（S8-ux）—— 内嵌 docs/data_preparation_zh.md
# ══════════════════════════════════════════════════════════
elif page == "数据准备指南":
    st.title("📚 数据准备指南")
    st.caption("上游工具 → EnvMeta 输入格式映射 · 覆盖 11 种常见流程")
    guide_path = Path(__file__).parent / "docs" / "data_preparation_zh.md"
    if guide_path.exists():
        # 下载按钮
        guide_bytes = guide_path.read_bytes()
        st.download_button(
            "⬇️ 下载指南（Markdown）",
            data=guide_bytes,
            file_name="data_preparation_zh.md",
            mime="text/markdown",
            key="guide_dl",
        )
        st.markdown("---")
        # 内嵌渲染
        st.markdown(guide_path.read_text(encoding="utf-8"), unsafe_allow_html=False)
    else:
        st.error(f"指南文件未找到：{guide_path}")
        st.info("请确认项目根目录下存在 `docs/data_preparation_zh.md`。")

# ══════════════════════════════════════════════════════════
# 文件管理（模块 A）
# ══════════════════════════════════════════════════════════
elif page == "文件管理":
    st.title("文件管理")

    uploaded = st.file_uploader(
        "上传数据文件（可多选）", accept_multiple_files=True, type=None,
        help="支持 .txt / .tsv / .csv / .spf。系统自动识别类型。",
    )

    for up in uploaded or []:
        if up.name not in st.session_state.files:
            try:
                df, enc, sep = read_table(up)
                up.seek(0)
                result = detect(up, filename=up.name)
                st.session_state.files[up.name] = {
                    "df": df, "type": result.file_type, "result": result,
                }
            except Exception as e:
                st.error(f"读取 {up.name} 失败：{e}")

    if not st.session_state.files:
        st.info("请上传文件。或在 PowerShell 里拖 `tests/sample_data/*` 里的文件做测试。")
    else:
        st.subheader(f"已识别文件（{len(st.session_state.files)} 个）")
        for fname, info in list(st.session_state.files.items()):
            with st.expander(f"📄 {fname} — {TYPE_BADGES[info['type']]}", expanded=False):
                col1, col2, col3 = st.columns([2, 2, 1])
                with col1:
                    reasons = info["result"].reasons
                    st.caption(
                        f"置信度 **{info['result'].confidence:.2f}** · "
                        f"编码 {info['result'].encoding} · "
                        f"分隔符 {'TAB' if info['result'].separator == chr(9) else info['result'].separator!r}"
                    )
                    if reasons:
                        st.caption("识别依据：" + "；".join(reasons))
                with col2:
                    new_type = st.selectbox(
                        "手动修正类型", TYPE_OPTIONS,
                        index=TYPE_OPTIONS.index(info["type"].value),
                        key=f"type_{fname}",
                    )
                    if new_type != info["type"].value:
                        st.session_state.files[fname]["type"] = FileType(new_type)
                with col3:
                    if st.button("移除", key=f"rm_{fname}"):
                        del st.session_state.files[fname]
                        st.rerun()
                st.dataframe(info["result"].preview_df, use_container_width=True, height=180)
                # S8-ux 反向索引 —— 这个文件可以跑什么分析
                st.markdown("---")
                _render_file_reverse_index(info["type"], key_prefix=f"rev_{fname}")

        # 汇总：当前文件组合可直接跑的分析
        have_types = {info["type"] for info in st.session_state.files.values()}
        ready_ids = analyses_ready(have_types)
        if ready_ids:
            st.success(
                f"✅ 当前文件组合可**直接运行**的分析："
                + "、".join(ANALYSIS_INPUTS[aid]["name"] for aid in ready_ids)
            )

# ══════════════════════════════════════════════════════════
# 图表选择向导（S8-ux）
# ══════════════════════════════════════════════════════════
elif page == "图表选择向导":
    st.title("📊 图表选择向导")
    st.caption("按研究问题推荐分析类型。两步 radio → 推荐卡片 → 一键跳转。")

    # 步骤 1 — 选大类
    categories = [f"{cat['icon']} {cat['category']}" for cat in NAVIGATOR]
    cat_idx = st.radio(
        "**第一步 · 你想回答什么研究问题？**",
        list(range(len(NAVIGATOR))),
        format_func=lambda i: categories[i],
        key="nav_cat",
    )
    cat = NAVIGATOR[cat_idx]
    st.caption(f"*{cat['description']}*")
    st.markdown("---")

    # 步骤 2 — 选具体子问题
    sub_opts = [sq["q"] for sq in cat["subquestions"]]
    sub_idx = st.radio(
        "**第二步 · 更具体一点…**",
        list(range(len(sub_opts))),
        format_func=lambda i: sub_opts[i],
        key=f"nav_sub_{cat_idx}",
    )
    sq = cat["subquestions"][sub_idx]
    st.markdown("---")

    # 步骤 3 — 推荐卡片
    st.subheader("🎯 推荐分析")
    have_types = {info["type"] for info in st.session_state.files.values()}
    ready = set(analyses_ready(have_types))
    for i, rec in enumerate(sorted(sq["recommended"], key=lambda r: r["priority"])):
        aid = rec["analysis_id"]
        spec = ANALYSIS_INPUTS.get(aid, {"name": aid})
        priority_label = {1: "🥇 首选", 2: "🥈 次选", 3: "🥉 补充"}.get(rec["priority"], "")
        ready_badge = " ✅ 文件已齐全" if aid in ready else ""
        with st.container(border=True):
            c1, c2 = st.columns([4, 1])
            with c1:
                st.markdown(f"### {priority_label} {spec['name']}{ready_badge}")
                st.markdown(f"**推荐理由**：{rec['reason']}")
                # 需要的文件类型
                req_badges = [TYPE_BADGES.get(ft, ft.value)
                              for ft in ANALYSIS_INPUTS.get(aid, {}).get("required", [])]
                opt_badges = [TYPE_BADGES.get(ft, ft.value)
                              for ft in ANALYSIS_INPUTS.get(aid, {}).get("optional", [])]
                if req_badges:
                    st.caption("**必需文件**：" + " · ".join(req_badges))
                if opt_badges:
                    st.caption("**可选文件**：" + " · ".join(opt_badges))
            with c2:
                _jump_to_analysis(aid, key=f"nav_jump_{cat_idx}_{sub_idx}_{i}")

    if sq.get("tip"):
        st.info(f"💡 **提示**：{sq['tip']}")

# ══════════════════════════════════════════════════════════
# Reads-based 分析（模块 B1）
# ══════════════════════════════════════════════════════════
elif page == "Reads-based 分析":
    st.title("Reads-based 群落分析")
    analysis_type = st.selectbox(
        "选择分析类型",
        [
            "物种组成堆叠图",
            "β多样性 PCoA",
            "功能基因热图",
            "α多样性",
            "RDA/CCA 排序",
            "LEfSe 差异分析",
            "基因差异分析 (log2FC)",
        ],
        key="reads_analysis_type",
    )

    # ── 堆叠图 ───────────────────────────────────────────
    if analysis_type == "物种组成堆叠图":
        _render_interpretation_expander("stackplot", key="interp_stackplot")
        abundance_files = _files_of(FileType.ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not abundance_files or not metadata_files:
            st.warning(
                f"需要 1 个丰度表 + 1 个 metadata。当前：丰度表 {len(abundance_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        ab_options = list(abundance_files.keys())
        ab_default = _first_taxon_abundance()
        md_options = list(metadata_files.keys())
        with c1:
            ab_name = st.selectbox("丰度表", ab_options,
                                   index=_idx_or_default(ab_options, ab_default),
                                   key="stack_ab")
        with c2:
            md_name = st.selectbox("Metadata", md_options, key="stack_md")

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        style = st.sidebar.radio("横轴", ["sample", "group", "combined"],
                                 horizontal=True, key="stack_style")
        top_n = st.sidebar.slider("Top-N", 5, 20, 10, key="stack_top_n")
        sort_by = st.sidebar.selectbox("排序依据", ["mean", "max", "median"], key="stack_sort_by")
        reverse_stack = st.sidebar.checkbox("高丰度在柱顶（反转）", value=False, key="stack_rev")
        drop_unc = st.sidebar.checkbox("过滤 unclassified 行", value=True, key="stack_drop_unc")
        size = render_figure_size({"width_mm": 160, "height_mm": 100}, prefix="stack")

        params = {
            "style": style, "top_n": top_n, "sort_by": sort_by,
            "reverse_stack": reverse_stack, "drop_unclassified": drop_unc, **size,
        }

        if st.button("🎨 生成图表", type="primary", key="stack_go"):
            try:
                result = stackplot.analyze(
                    abundance_files[ab_name]["df"], metadata_files[md_name]["df"], params)
                st.session_state["last_stackplot"] = result
                _register_export(
                    "stackplot", result=result,
                    file_paths={"abundance": ab_name, "metadata": md_name},
                    params=result.params,
                    output_base=f"stackplot_{result.params['style']}",
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_stackplot")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)
            c_d1, c_d2, c_d3 = st.columns(3)
            with c_d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name=f"stackplot_{last.params['style']}.png",
                                   mime="image/png", key="stack_png")
            with c_d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name=f"stackplot_{last.params['style']}.pdf",
                                   mime="application/pdf", key="stack_pdf")
            with c_d3:
                st.download_button("⬇️ 百分比表（TSV）",
                                   data=last.stats.to_csv(sep="\t").encode("utf-8"),
                                   file_name=f"stackplot_{last.params['style']}_percentage.tsv",
                                   mime="text/tab-separated-values", key="stack_tsv")
            _vector_downloads(last.figure, f"stackplot_{last.params['style']}", "stack")
            _reproduce_button("stackplot",
                              {"abundance": ab_name, "metadata": md_name},
                              last.params, key="stack_code",
                              output_base=f"stackplot_{last.params['style']}")
            with st.expander("查看百分比表"):
                st.dataframe(last.stats.round(2), use_container_width=True)

    # ── PCoA ────────────────────────────────────────────
    elif analysis_type == "β多样性 PCoA":
        _render_interpretation_expander("pcoa", key="interp_pcoa")
        dist_files = _files_of(FileType.DISTANCE_MATRIX)
        metadata_files = _files_of(FileType.METADATA)
        if not dist_files or not metadata_files:
            st.warning(
                f"需要 1 个距离矩阵（如 Bray-Curtis）+ 1 个 metadata。"
                f"当前：距离矩阵 {len(dist_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        with c1:
            dist_name = st.selectbox("距离矩阵", list(dist_files.keys()))
        with c2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="pcoa_md")

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        show_labels = st.sidebar.checkbox("显示样本标签", value=True, key="pcoa_show_labels")
        show_perm = st.sidebar.checkbox("显示 PERMANOVA 结果", value=True, key="pcoa_show_perm")
        n_perm = st.sidebar.select_slider("置换次数", options=[99, 499, 999, 4999],
                                          value=999, key="pcoa_n_perm")
        marker = st.sidebar.slider("点大小", 30, 200, 90, key="pcoa_marker")
        size = render_figure_size({"width_mm": 130, "height_mm": 110}, prefix="pcoa")
        fonts = render_font_controls({"label_size": 9, "tick_size": 8}, prefix="pcoa")

        params = {
            "show_labels": show_labels, "show_permanova": show_perm,
            "n_permutations": n_perm, "marker_size": marker, **size, **fonts,
        }

        if st.button("🎨 生成图表", type="primary", key="pcoa_go"):
            try:
                result = pcoa.analyze(
                    dist_files[dist_name]["df"], metadata_files[md_name]["df"], params)
                st.session_state["last_pcoa"] = result
                _register_export(
                    "pcoa", result=result,
                    file_paths={"distance": dist_name, "metadata": md_name},
                    params=result.params,
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_pcoa")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)
            c_d1, c_d2, c_d3 = st.columns(3)
            with c_d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name="pcoa.png", mime="image/png", key="pcoa_png")
            with c_d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="pcoa.pdf", mime="application/pdf", key="pcoa_pdf")
            with c_d3:
                st.download_button("⬇️ 统计汇总（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="pcoa_stats.tsv", mime="text/tab-separated-values",
                                   key="pcoa_tsv")
            _vector_downloads(last.figure, "pcoa", "pcoa")
            _reproduce_button("pcoa",
                              {"distance": dist_name, "metadata": md_name},
                              last.params, key="pcoa_code")
            with st.expander("PERMANOVA 两两对比"):
                st.dataframe(last.params["_stats_tables"]["pairwise_permanova"],
                             use_container_width=True)

    # ── 元素循环基因热图 ───────────────────────────────
    elif analysis_type == "功能基因热图":
        _render_interpretation_expander("gene_heatmap", key="interp_gene_heatmap")
        ko_files = _files_of(FileType.KO_ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not ko_files or not metadata_files:
            st.warning(
                f"需要 1 个 KO 丰度表（含 KEGG_ko 列）+ 1 个 metadata。"
                f"当前：KO 表 {len(ko_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        with c1:
            ko_name = st.selectbox("KO 丰度表", list(ko_files.keys()))
        with c2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="heat_md")

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        elements = st.sidebar.multiselect(
            "显示元素", ["arsenic", "nitrogen", "sulfur", "iron"],
            default=["arsenic", "nitrogen", "sulfur", "iron"], key="heat_elements",
        )
        zscore = st.sidebar.checkbox("行 Z-score 归一", value=True, key="heat_zscore")
        show_genes = st.sidebar.checkbox("显示基因名", value=True, key="heat_show_genes")
        cmap = st.sidebar.selectbox("色阶", ["RdBu_r", "coolwarm", "RdYlBu_r", "viridis"],
                                    key="heat_cmap")
        size = render_figure_size({"width_mm": 180, "height_mm": 220}, prefix="heat")

        params = {
            "element_filter": elements if elements else None,
            "zscore": zscore, "show_gene_names": show_genes, "cmap": cmap, **size,
        }

        if st.button("🎨 生成图表", type="primary", key="heat_go"):
            try:
                result = gene_heatmap.analyze(
                    ko_files[ko_name]["df"], metadata_files[md_name]["df"], params)
                st.session_state["last_heat"] = result
                _register_export(
                    "gene_heatmap", result=result,
                    file_paths={"ko_abundance": ko_name, "metadata": md_name},
                    params=result.params,
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_heat")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            c_d1, c_d2, c_d3 = st.columns(3)
            with c_d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name="gene_heatmap.png", mime="image/png", key="heat_png")
            with c_d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="gene_heatmap.pdf", mime="application/pdf",
                                   key="heat_pdf")
            with c_d3:
                st.download_button("⬇️ Z-score 矩阵（TSV）",
                                   data=last.stats.to_csv(sep="\t").encode("utf-8"),
                                   file_name="gene_heatmap_zscore.tsv",
                                   mime="text/tab-separated-values", key="heat_tsv")
            _vector_downloads(last.figure, "gene_heatmap", "heat")
            _reproduce_button("gene_heatmap",
                              {"ko_abundance": ko_name, "metadata": md_name},
                              last.params, key="heat_code")
            with st.expander("原始分组均值矩阵"):
                st.dataframe(last.params["_raw_means"].round(2),
                             use_container_width=True)

    # ── α 多样性箱线图 ─────────────────────────────────
    elif analysis_type == "α多样性":
        _render_interpretation_expander("alpha_boxplot", key="interp_alpha_boxplot")
        alpha_files = _files_of(FileType.ALPHA_DIVERSITY)
        metadata_files = _files_of(FileType.METADATA)
        if not alpha_files or not metadata_files:
            st.warning(
                f"需要 1 个 α 多样性表 + 1 个 metadata。当前：α 表 {len(alpha_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        with c1:
            alpha_name = st.selectbox("α 多样性表", list(alpha_files.keys()), key="alpha_file")
        with c2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="alpha_md")

        alpha_df = alpha_files[alpha_name]["df"]
        sample_col = next((c for c in alpha_df.columns
                           if c.lower().replace("_", "") == "sampleid"), alpha_df.columns[0])
        candidate_metrics = [c for c in alpha_df.columns
                             if c != sample_col
                             and pd.to_numeric(alpha_df[c], errors="coerce").notna().mean() > 0.95]

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        metrics = st.sidebar.multiselect(
            "展示指数", candidate_metrics, default=candidate_metrics, key="alpha_metrics"
        )
        n_cols = st.sidebar.slider("子图列数", 1, 5, min(3, len(metrics) or 1), key="alpha_ncols")
        show_pvalues = st.sidebar.checkbox("显示组间 p 值", value=True, key="alpha_showp")
        alpha_cut = st.sidebar.slider("p 阈值", 0.01, 0.10, 0.05, step=0.01, key="alpha_cut")
        size = render_figure_size({"width_mm": 200, "height_mm": 120}, prefix="alpha")

        params = {
            "metrics": metrics or None, "n_cols": n_cols,
            "show_pvalues": show_pvalues, "alpha": alpha_cut, **size,
        }

        if st.button("🎨 生成图表", type="primary", key="alpha_go"):
            try:
                from envmeta.analysis import alpha_boxplot
                result = alpha_boxplot.analyze(alpha_df, metadata_files[md_name]["df"], params)
                st.session_state["last_alpha"] = result
                _register_export(
                    "alpha_boxplot", result=result,
                    file_paths={"alpha": alpha_name, "metadata": md_name},
                    params=result.params,
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_alpha")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name="alpha_boxplot.png", mime="image/png", key="alpha_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="alpha_boxplot.pdf", mime="application/pdf", key="alpha_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="alpha_stats.tsv",
                                   mime="text/tab-separated-values", key="alpha_tsv")
            _vector_downloads(last.figure, "alpha_boxplot", "alpha")
            _reproduce_button("alpha_boxplot",
                              {"alpha": alpha_name, "metadata": md_name},
                              last.params, key="alpha_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # ── 基因差异分析 log2FC ─────────────────────────────
    elif analysis_type == "基因差异分析 (log2FC)":
        _render_interpretation_expander("log2fc", key="interp_log2fc")
        ko_files = _files_of(FileType.KO_ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not ko_files or not metadata_files:
            st.warning(
                f"需要 1 个 KO 丰度表 + 1 个 metadata。当前：KO 表 {len(ko_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        with c1:
            ko_name = st.selectbox("KO 丰度表", list(ko_files.keys()), key="log2fc_ko")
        with c2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="log2fc_md")

        md_df = metadata_files[md_name]["df"]
        group_opts = list(dict.fromkeys(md_df["Group"].astype(str).tolist())) if "Group" in md_df.columns else []

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        if len(group_opts) < 2:
            st.error("metadata 的 Group 列至少需要 2 个不同值")
            st.stop()
        group_a = st.sidebar.selectbox("组 A（分子）", group_opts, index=0, key="log2fc_a")
        group_b = st.sidebar.selectbox("组 B（分母）", group_opts,
                                       index=1 if group_opts[0] == group_a else 0,
                                       key="log2fc_b")
        element_filter = st.sidebar.multiselect(
            "元素", ["arsenic", "nitrogen", "sulfur", "iron"],
            default=["arsenic", "nitrogen", "sulfur", "iron"], key="log2fc_el",
        )
        alpha_cut = st.sidebar.slider("padj 阈值", 0.01, 0.10, 0.05, step=0.01, key="log2fc_alpha")
        lfc_thresh = st.sidebar.slider("|log2FC| 阈值", 0.0, 2.0, 1.0, step=0.1, key="log2fc_thresh")
        size = render_figure_size({"width_mm": 220, "height_mm": 200}, prefix="log2fc")

        params = {
            "group_a": group_a, "group_b": group_b,
            "element_filter": element_filter or None,
            "alpha": alpha_cut, "log2fc_threshold": lfc_thresh, **size,
        }

        if st.button("🎨 生成图表", type="primary", key="log2fc_go"):
            try:
                from envmeta.analysis import log2fc
                result = log2fc.analyze(ko_files[ko_name]["df"],
                                         md_df, params)
                st.session_state["last_log2fc"] = result
                _register_export(
                    "log2fc", result=result,
                    file_paths={"ko_abundance": ko_name, "metadata": md_name},
                    params=result.params,
                    output_base=f"log2fc_{group_a}_vs_{group_b}",
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_log2fc")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)
            d1, d2, d3 = st.columns(3)
            suffix = f"{last.params['group_a']}_vs_{last.params['group_b']}"
            with d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name=f"log2fc_{suffix}.png",
                                   mime="image/png", key="log2fc_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name=f"log2fc_{suffix}.pdf",
                                   mime="application/pdf", key="log2fc_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name=f"log2fc_{suffix}_stats.tsv",
                                   mime="text/tab-separated-values", key="log2fc_tsv")
            _vector_downloads(last.figure, f"log2fc_{suffix}", "log2fc")
            _reproduce_button("log2fc",
                              {"ko_abundance": ko_name, "metadata": md_name},
                              last.params, key="log2fc_code",
                              output_base=f"log2fc_{suffix}")
            with st.expander("查看统计表（按 padj 排序）"):
                st.dataframe(last.stats.sort_values("padj").round(4),
                             use_container_width=True)

    # ── RDA 排序图 ─────────────────────────────────
    elif analysis_type == "RDA/CCA 排序":
        _render_interpretation_expander("rda", key="interp_rda")
        abundance_files = _files_of(FileType.ABUNDANCE_WIDE)
        env_files = _files_of(FileType.ENV_FACTORS)
        metadata_files = _files_of(FileType.METADATA)
        if not abundance_files or not env_files or not metadata_files:
            st.warning(
                f"需要 1 丰度表 + 1 环境因子表 + 1 metadata。当前："
                f"丰度 {len(abundance_files)}，env {len(env_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2, c3 = st.columns(3)
        ab_options = list(abundance_files.keys())
        ab_default = _first_taxon_abundance()
        with c1:
            ab_name = st.selectbox("丰度表", ab_options,
                                   index=_idx_or_default(ab_options, ab_default),
                                   key="rda_ab")
        with c2:
            env_name = st.selectbox("环境因子表", list(env_files.keys()), key="rda_env")
        with c3:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="rda_md")

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        n_perm = st.sidebar.select_slider("Mantel 置换次数",
                                          options=[99, 499, 999], value=999, key="rda_n_perm")
        show_arrows = st.sidebar.checkbox("显示环境因子箭头", value=True, key="rda_arrows")
        show_labels = st.sidebar.checkbox("显示样本标签", value=True, key="rda_labels")
        arrow_scale = st.sidebar.slider("箭头缩放", 0.5, 3.0, 1.0, step=0.1, key="rda_scale")
        drop_unc = st.sidebar.checkbox("过滤 unclassified", value=True, key="rda_drop_unc")
        size = render_figure_size({"width_mm": 180, "height_mm": 140}, prefix="rda")

        params = {
            "n_permutations": n_perm, "show_env_arrows": show_arrows,
            "show_sample_labels": show_labels, "arrow_scale": arrow_scale,
            "drop_unclassified": drop_unc, **size,
        }

        if st.button("🎨 生成图表", type="primary", key="rda_go"):
            try:
                result = rda_mod.analyze(
                    abundance_files[ab_name]["df"],
                    env_files[env_name]["df"],
                    metadata_files[md_name]["df"],
                    params,
                )
                st.session_state["last_rda"] = result
                _register_export(
                    "rda", result=result,
                    file_paths={"abundance": ab_name, "env_factors": env_name, "metadata": md_name},
                    params=result.params,
                )
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_rda")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG（300 DPI）", data=export_to_bytes(last.figure, "png"),
                                   file_name="rda.png", mime="image/png", key="rda_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="rda.pdf", mime="application/pdf", key="rda_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="rda_stats.tsv",
                                   mime="text/tab-separated-values", key="rda_tsv")
            _vector_downloads(last.figure, "rda", "rda")
            _reproduce_button("rda",
                              {"abundance": ab_name, "env_factors": env_name, "metadata": md_name},
                              last.params, key="rda_code")
            with st.expander("查看 RDA 统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # ── LEfSe 差异分析 ──────────────────────────────────
    elif analysis_type == "LEfSe 差异分析":
        _render_interpretation_expander("lefse", key="interp_lefse")
        abundance_files = _files_of(FileType.ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not abundance_files or not metadata_files:
            st.warning(
                f"需要 1 个丰度表（Species/Genus/Phylum）+ 1 个 metadata。"
                f"当前：丰度表 {len(abundance_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()
        ab_options = list(abundance_files.keys())
        ab_default = _first_taxon_abundance()
        ab_name = st.selectbox("丰度表", ab_options,
                               index=_idx_or_default(ab_options, ab_default),
                               key="lefse_ab")
        md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="lefse_md")

        with st.sidebar:
            st.subheader("LEfSe 参数")
            alpha_kw = st.slider("Kruskal-Wallis α", 0.01, 0.2, 0.05, step=0.01,
                                 key="lefse_alpha")
            lda_thresh = st.slider("LDA 阈值（log10）", 1.0, 5.0, 2.0, step=0.1,
                                   key="lefse_lda")
            tax_all = ["Kingdom", "Phylum", "Class", "Order", "Family",
                       "Genus", "Species", "Strain"]
            tax_levels = st.multiselect(
                "限制分类层级（空 = 全部）", tax_all, default=[],
                key="lefse_tax",
            )
            max_features = st.slider("最多展示特征数（0 = 全部）", 0, 100, 40,
                                     step=5, key="lefse_max")
            size = render_figure_size({"width_mm": 180, "height_mm": 200}, prefix="lefse")

        if st.button("生成 LEfSe 图", type="primary", key="lefse_go"):
            ab = abundance_files[ab_name]["df"]
            md = metadata_files[md_name]["df"]
            params = {
                "alpha_kw": alpha_kw,
                "lda_threshold": lda_thresh,
                "tax_levels": tax_levels or None,
                "max_features": max_features or None,
                **size,
            }
            try:
                result = lefse_mod.analyze(ab, md, params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_lefse_last"] = result
            _register_export(
                "lefse", result=result,
                file_paths={"abundance": ab_name, "metadata": md_name},
                params=result.params,
            )

        last = st.session_state.get("_lefse_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="lefse.png", mime="image/png", key="lefse_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="lefse.pdf", mime="application/pdf", key="lefse_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="lefse_stats.tsv",
                                   mime="text/tab-separated-values", key="lefse_tsv")
            _vector_downloads(last.figure, "lefse", "lefse")
            _reproduce_button("lefse",
                              {"abundance": ab_name, "metadata": md_name},
                              last.params, key="lefse_code")
            with st.expander("查看 LEfSe 统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # 未实现的子页面
    else:
        st.info(f"「{analysis_type}」将在 Phase 2 实现。")

# ══════════════════════════════════════════════════════════
# 其他占位页面
# ══════════════════════════════════════════════════════════
elif page == "MAG-based 分析":
    st.title("MAG-based 基因组分析")
    analysis_type = st.selectbox(
        "选择分析类型",
        ["MAG 质量评估", "MAG 丰度热图", "代谢通路完整度", "MAG 元素循环基因谱", "共现网络图"],
        key="mag_analysis_type",
    )

    if analysis_type == "MAG 质量评估":
        _render_interpretation_expander("mag_quality", key="interp_mag_quality")
        checkm_files = _files_of(FileType.CHECKM_QUALITY)
        if not checkm_files:
            st.warning("需要 1 个 CheckM / CheckM2 质量表（含 Completeness/Contamination/Genome_Size 列）。")
            st.stop()
        q_name = st.selectbox("CheckM 质量表（必需）", list(checkm_files.keys()),
                              key="mq_checkm")
        t_default = _first_of(FileType.MAG_TAXONOMY)
        t_options = ["（无）"] + [n for n in st.session_state.files
                                 if n != q_name]
        t_name = st.selectbox(
            "GTDB 分类表（可选，用于门分色）", t_options, key="mq_tax",
            index=_idx_or_default(t_options, t_default),
        )
        k_default = _first_of(FileType.KEYSTONE_SPECIES)
        k_options = ["（无）"] + [n for n in st.session_state.files
                                 if n not in (q_name, t_name)]
        k_name = st.selectbox(
            "Keystone 物种列表（可选）", k_options, key="mq_ks",
            index=_idx_or_default(k_options, k_default),
        )

        with st.sidebar:
            st.subheader("MAG 质量参数")
            l1 = render_mag_layer1_filter("mq", default_max_mags=0)
            # 默认 all —— 散点图核心价值 = 看全部 MAG 质量分布
            l1["filter_mode"] = st.selectbox(
                "（覆盖）MAG 子集",
                ["all", "top_n", "keystone_only", "top_plus_keystone"],
                index=0, key="mq_fm_override",
                help="散点图默认 all（展示全部 MAG 质量分布）；选 keystone_only 可单看 keystone。",
            )
            l2 = render_mag_layer2_visual("mq", 260, 140, show_phylum_bar=False)
            st.markdown("**质量阈值**")
            hc = st.slider("高质量 Completeness ≥", 50, 100, 90, key="mq_hc")
            hcon = st.slider("高质量 Contamination <", 1, 10, 5, key="mq_hcon")
            mc = st.slider("中质量 Completeness ≥", 30, 80, 50, key="mq_mc")
            mcon = st.slider("中质量 Contamination <", 5, 20, 10, key="mq_mcon")

        if st.button("生成 MAG 质量图", type="primary", key="mq_go"):
            q = checkm_files[q_name]["df"]
            t = _load_mag_df(t_name if t_name != "（无）" else None)
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            params = {
                **l1, **l2,
                "high_completeness": float(hc), "high_contamination": float(hcon),
                "med_completeness": float(mc), "med_contamination": float(mcon),
            }
            try:
                result = mag_quality_mod.analyze(q, t, k, params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_mq_last"] = result

        last = st.session_state.get("_mq_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="mag_quality.png", mime="image/png", key="mq_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="mag_quality.pdf", mime="application/pdf", key="mq_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="mag_quality_stats.tsv",
                                   mime="text/tab-separated-values", key="mq_tsv")
            files_map = {"quality": q_name}
            if t_name != "（无）": files_map["taxonomy"] = t_name
            if k_name != "（无）": files_map["keystone"] = k_name
            _register_export("mag_quality", result=last, file_paths=files_map, params=last.params)
            _vector_downloads(last.figure, "mag_quality", "mq")
            _reproduce_button("mag_quality", files_map, last.params, key="mq_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "MAG 丰度热图":
        _render_interpretation_expander("mag_heatmap", key="interp_mag_heatmap")
        sel = _mag_file_selectors("mh", require_abundance=True)
        ab_name = sel["abundance"]

        with st.sidebar:
            st.subheader("MAG 丰度热图参数")
            l1 = render_mag_layer1_filter("mh", default_max_mags=0)
            l2 = render_mag_layer2_visual("mh", 240, 220)
            l3 = render_mag_layer3_order("mh")
            st.markdown("**配色 / 聚类（丰度热图特有）**")
            bp_lo = st.slider("配色低段边界 (%)", 0.05, 0.50, 0.20,
                              step=0.05, key="mh_bplo")
            bp_hi = st.slider("配色高段边界 (%)", 0.30, 1.00, 0.50,
                              step=0.05, key="mh_bphi")
            log_tr = st.checkbox("聚类前 log1p 变换（长尾分布）",
                                 True, key="mh_log")
            cluster_cols = st.checkbox("列聚类（样本）", False, key="mh_cc")
            show_grp = st.checkbox("组彩条（顶部）", True, key="mh_grp")

        if st.button("生成 MAG 丰度热图", type="primary", key="mh_go"):
            ab = st.session_state.files[ab_name]["df"]
            t = _load_mag_df(sel["taxonomy"])
            k = _load_mag_df(sel["keystone"])
            m = _load_mag_df(sel["metadata"])
            params = {
                **l1, **l2, **l3,
                "color_breakpoints": (bp_lo, max(bp_hi, bp_lo + 0.05)),
                "log_transform": log_tr,
                "cluster_cols": cluster_cols,
                "show_group_bar": show_grp,
            }
            try:
                result = mag_heatmap_mod.analyze(ab, t, k, m, params=params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_mh_last"] = result

        last = st.session_state.get("_mh_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="mag_heatmap.png", mime="image/png",
                                   key="mh_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="mag_heatmap.pdf",
                                   mime="application/pdf", key="mh_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="mag_heatmap_stats.tsv",
                                   mime="text/tab-separated-values", key="mh_tsv")
            files_map = {"abundance": ab_name}
            for k_, v_ in [("taxonomy", sel["taxonomy"]),
                           ("keystone", sel["keystone"]),
                           ("metadata", sel["metadata"])]:
                if v_:
                    files_map[k_] = v_
            _register_export("mag_heatmap", result=last, file_paths=files_map, params=last.params)
            _vector_downloads(last.figure, "mag_heatmap", "mh")
            _reproduce_button("mag_heatmap", files_map, last.params, key="mh_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "代谢通路完整度":
        _render_interpretation_expander("pathway", key="interp_pathway")
        sel = _mag_file_selectors("pw", ko_type=FileType.KO_ANNOTATION_LONG)
        ko_name = sel["_ko"]

        with st.sidebar:
            st.subheader("通路完整度参数")
            l1 = render_mag_layer1_filter("pw", default_max_mags=50)
            l2 = render_mag_layer2_visual("pw", 320, 220)
            l3 = render_mag_layer3_order("pw")
            st.markdown("**通路特有**")
            style = st.radio("图样式", ["heatmap", "bubble"], key="pw_style")
            elem_opts = ["arsenic", "sulfur", "iron", "nitrogen"]
            elem_sel = st.multiselect("元素过滤（空=全部）", elem_opts,
                                      default=[], key="pw_elem")
            min_comp = st.slider("最低总完整度", 0, 200, 0, step=10,
                                  key="pw_mc",
                                  help="MAG 总完整度 < 此值时不显示")
            bubble_scale = st.slider("bubble 缩放（仅 bubble）", 1.0, 10.0,
                                      4.0, step=0.5, key="pw_bs")

        if st.button("生成通路完整度图", type="primary", key="pw_go"):
            ko = st.session_state.files[ko_name]["df"]
            t = _load_mag_df(sel["taxonomy"])
            k = _load_mag_df(sel["keystone"])
            a = _load_mag_df(sel["abundance"])
            params = {
                **l1, **l2, **l3,
                "style": style,
                "element_filter": elem_sel or None,
                "min_completeness": float(min_comp),
                "bubble_scale": float(bubble_scale),
            }
            try:
                result = pathway_mod.analyze(ko, t, k, a, params=params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_pw_last"] = result

        last = st.session_state.get("_pw_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="pathway.png", mime="image/png",
                                   key="pw_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="pathway.pdf",
                                   mime="application/pdf", key="pw_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="pathway_stats.tsv",
                                   mime="text/tab-separated-values", key="pw_tsv")
            files_map = {"ko_annotation": ko_name}
            for k_, v_ in [("taxonomy", sel["taxonomy"]),
                           ("keystone", sel["keystone"]),
                           ("abundance", sel["abundance"])]:
                if v_:
                    files_map[k_] = v_
            _register_export("pathway", result=last, file_paths=files_map, params=last.params)
            _vector_downloads(last.figure, "pathway", "pw")
            _reproduce_button("pathway", files_map, last.params, key="pw_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "MAG 元素循环基因谱":
        _render_interpretation_expander("gene_profile", key="interp_gene_profile")
        sel = _mag_file_selectors("gp", ko_type=FileType.KO_ANNOTATION_LONG)
        ko_name = sel["_ko"]

        with st.sidebar:
            st.subheader("基因谱参数")
            l1 = render_mag_layer1_filter("gp", default_max_mags=60)
            l2 = render_mag_layer2_visual("gp", 340, 220)
            l3 = render_mag_layer3_order("gp")
            st.markdown("**基因谱特有**")
            cmap_opts = ["viridis", "plasma", "YlGnBu", "YlOrBr"]
            cmap_name = st.selectbox("配色", cmap_opts, index=0, key="gp_cmap",
                help="viridis（默认）感知均匀、色域宽；YlGnBu 双色冷暖过渡")
            blank_zeros = st.checkbox(
                "0 值留白（缺失 vs 低表达一眼分）",
                True, key="gp_bz",
                help=(
                    "推荐开启。开启后 0 拷贝数显示为白色，"
                    "彩色 = 有基因存在；关闭则 0 显 cmap 最低色（深色）"
                    "可能占据大部分视觉。"
                ),
            )
            sort_ko_cov = st.checkbox(
                "KO 列按覆盖率排序（左密右稀）",
                False, key="gp_skc",
                help=(
                    "勾选后 KO 列按「有多少 MAG 持有该 KO」降序；"
                    "高覆盖 KO 集中在左侧，稀疏 KO 在右。"
                    "适合看子集 MAG 的共有功能。"
                ),
            )
            elem_opts = ["arsenic", "sulfur", "iron", "nitrogen"]
            elem_sel = st.multiselect("元素过滤（空=全部）", elem_opts,
                                      default=[], key="gp_elem")
            show_names = st.checkbox("列标签显示基因名", True, key="gp_names")
            drop_zero = st.checkbox("过滤全 0 KO 列", True, key="gp_dz")

        if st.button("生成基因谱图", type="primary", key="gp_go"):
            ko = st.session_state.files[ko_name]["df"]
            t = _load_mag_df(sel["taxonomy"])
            k = _load_mag_df(sel["keystone"])
            a = _load_mag_df(sel["abundance"])
            params = {
                **l1, **l2, **l3,
                "cmap_name": cmap_name,
                "blank_zeros": blank_zeros,
                "sort_ko_by_coverage": sort_ko_cov,
                "element_filter": elem_sel or None,
                "show_gene_names": show_names,
                "drop_zero_kos": drop_zero,
            }
            try:
                result = gene_profile_mod.analyze(ko, t, k, a, params=params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_gp_last"] = result

        last = st.session_state.get("_gp_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="gene_profile.png", mime="image/png",
                                   key="gp_png")
            with d2:
                st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                                   file_name="gene_profile.pdf",
                                   mime="application/pdf", key="gp_pdf")
            with d3:
                st.download_button("⬇️ 统计表（TSV）",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="gene_profile_stats.tsv",
                                   mime="text/tab-separated-values", key="gp_tsv")
            files_map = {"ko_annotation": ko_name}
            for k_, v_ in [("taxonomy", sel["taxonomy"]),
                           ("keystone", sel["keystone"]),
                           ("abundance", sel["abundance"])]:
                if v_:
                    files_map[k_] = v_
            _register_export("gene_profile", result=last, file_paths=files_map, params=last.params)
            _vector_downloads(last.figure, "gene_profile", "gp")
            _reproduce_button("gene_profile", files_map, last.params, key="gp_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "共现网络图":
        _render_interpretation_expander("network", key="interp_network")
        st.caption("Gephi 辅助工具 — 散点图 + CSV 预处理 + 参数指南（不画网络图本体，交给 Gephi）")
        file_names = list(st.session_state.files.keys())
        if not file_names:
            st.warning("请先在「文件管理」上传 Gephi nodes CSV + edges CSV。")
            st.stop()
        n_default = _first_of(FileType.GEPHI_NODES)
        n_name = st.selectbox("Gephi Nodes CSV（必需：Id + Degree + Betweenness）",
                              file_names, key="nw_nodes",
                              index=_idx_or_default(file_names, n_default))
        e_default = _first_of(FileType.GEPHI_EDGES)
        e_options = [n for n in file_names if n != n_name]
        e_name = st.selectbox("Gephi Edges CSV（必需：Source + Target + Weight）",
                              e_options, key="nw_edges",
                              index=_idx_or_default(e_options, e_default))
        t_default = _first_of(FileType.MAG_TAXONOMY)
        t_options = ["（无）"] + [n for n in file_names if n not in (n_name, e_name)]
        t_name = st.selectbox("GTDB 分类表（可选，补 Genus 标签）",
                              t_options, key="nw_tax",
                              index=_idx_or_default(t_options, t_default))
        k_default = _first_of(FileType.KEYSTONE_SPECIES)
        k_options = ["（无）"] + [n for n in file_names if n not in (n_name, e_name, t_name)]
        k_name = st.selectbox("Keystone 物种列表（可选）",
                              k_options, key="nw_ks",
                              index=_idx_or_default(k_options, k_default))

        with st.sidebar:
            st.subheader("网络辅助参数")
            l1 = render_mag_layer1_filter("nw", default_max_mags=0)
            # 网络默认 all
            l1["filter_mode"] = st.selectbox(
                "（覆盖）MAG 子集",
                ["all", "top_n", "keystone_only", "top_plus_keystone"],
                index=0, key="nw_fm_override",
                help="网络图默认 all（保持拓扑完整性）；选 keystone_only 只看 keystone 散点。",
            )
            l2 = render_mag_layer2_visual("nw", 260, 200, show_phylum_bar=False)
            st.markdown("**网络特有**")
            deg_thr = st.slider("Degree 阈值线", 1, 30, 10, key="nw_deg")
            bet_thr = st.slider("Betweenness 阈值线", 50, 1000, 200,
                                step=50, key="nw_bet")
            color_by = st.selectbox("节点着色",
                                    ["keystone", "phylum"], key="nw_clr")

        # ── 散点图 ──────────────────────────────────────────
        if st.button("生成 Degree vs Betweenness 散点图",
                     type="primary", key="nw_go"):
            nodes = st.session_state.files[n_name]["df"]
            edges = st.session_state.files[e_name]["df"]
            t = _load_mag_df(t_name if t_name != "（无）" else None)
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            params = {
                **l1, **l2,
                "degree_threshold": deg_thr,
                "betweenness_threshold": bet_thr,
                "node_color_by": color_by,
            }
            try:
                result = network_mod.analyze(nodes, edges, t, k, params=params)
            except ValueError as e:
                st.error(str(e))
                st.stop()
            st.session_state["_nw_last"] = result

        last = st.session_state.get("_nw_last")
        if last is not None:
            st.pyplot(last.figure, use_container_width=False)
            d1, d2, d3 = st.columns(3)
            with d1:
                st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                                   file_name="network_scatter.png",
                                   mime="image/png", key="nw_png")
            with d2:
                st.download_button("⬇️ PDF",
                                   data=export_to_bytes(last.figure, "pdf"),
                                   file_name="network_scatter.pdf",
                                   mime="application/pdf", key="nw_pdf")
            with d3:
                st.download_button("⬇️ 节点表 TSV",
                                   data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                   file_name="network_node_stats.tsv",
                                   mime="text/tab-separated-values", key="nw_tsv")
            files_map = {"nodes": n_name, "edges": e_name}
            if t_name != "（无）": files_map["taxonomy"] = t_name
            if k_name != "（无）": files_map["keystone"] = k_name
            _register_export("network", result=last, file_paths=files_map, params=last.params)
            _vector_downloads(last.figure, "network", "nw")
            _reproduce_button("network", files_map, last.params, key="nw_code")
            with st.expander("查看节点统计表"):
                st.dataframe(last.stats, use_container_width=True)

        # ── Gephi 预处理导出 ────────────────────────────────
        st.markdown("---")
        st.subheader("Gephi 预处理导出")
        st.caption("帮你清理 CSV 标签 → 导入 Gephi 后 keystone 直接有标签、其余无标签")
        label_mode = st.selectbox(
            "标签模式",
            ["keystone_only", "all", "none"],
            index=0, key="nw_labmode",
            help=(
                "• keystone_only（推荐）: 只有 keystone 的 Label 有值\n"
                "  → 导入 Gephi 后只显 keystone 标签，不用手动一个个删\n"
                "• all: 所有节点加 Genus 标签\n"
                "• none: Label 列全空"
            ),
        )
        if st.button("校验 + 导出 Gephi 就绪 CSV", key="nw_gephi"):
            nodes = st.session_state.files[n_name]["df"]
            edges = st.session_state.files[e_name]["df"]
            t = _load_mag_df(t_name if t_name != "（无）" else None)
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            # 校验
            issues = validate_gephi_format(nodes, edges)
            if issues:
                for iss in issues:
                    if "[ERROR]" in iss:
                        st.error(iss)
                    else:
                        st.warning(iss)
                if any("[ERROR]" in i for i in issues):
                    st.stop()
            else:
                st.success("格式校验通过 — 可安全导入 Gephi")
            # 预处理
            out_n, out_e = prepare_gephi_csv(
                nodes, edges, taxonomy_df=t, keystone_df=k,
                label_mode=label_mode,
            )
            c1, c2 = st.columns(2)
            with c1:
                st.download_button(
                    "⬇️ nodes_gephi.csv",
                    data=out_n.to_csv(index=False).encode("utf-8"),
                    file_name="nodes_gephi.csv", mime="text/csv",
                    key="nw_gn",
                )
            with c2:
                st.download_button(
                    "⬇️ edges_gephi.csv",
                    data=out_e.to_csv(index=False).encode("utf-8"),
                    file_name="edges_gephi.csv", mime="text/csv",
                    key="nw_ge",
                )
            with st.expander("预览 nodes_gephi.csv"):
                st.dataframe(out_n.head(20), use_container_width=True)

        # ── Gephi 推荐参数指南 ──────────────────────────────
        with st.expander("Gephi 操作指南（推荐参数）"):
            st.markdown("""
### 1. 导入 CSV

1. 打开 Gephi → **File → Import spreadsheet**
2. 先导入 `nodes_gephi.csv`（as **Nodes Table**）
3. 再导入 `edges_gephi.csv`（as **Edges Table**, Undirected）

### 2. 布局

| 布局算法 | 参数 | 推荐值 | 说明 |
|---|---|---|---|
| **Fruchterman Reingold** | 区 | 10000 | 控制节点间距 |
|  | 重力 | 10.0 | 中心吸引力 |
|  | 速度 | 1.0 | 收敛速率 |
| ForceAtlas2（备选） | Scaling | 100 | 全局缩放 |
|  | Gravity | 1.0 | |
|  | Edge Weight Influence | 1.0 | |

运行布局直到节点位置稳定（~30s），点「停止」固定。

### 3. 外观

| 设置 | 位置 | 推荐值 |
|---|---|---|
| 节点大小 | 外观 → 节点 → 排名 → Degree | 最小=1, 最大=4 |
| 节点颜色 | 外观 → 节点 → 分区 → is_keystone | false=`#B0D4E8` true=`#1B3A5C` |
| 边颜色 | 外观 → 边 → 统一 | `#CCCCCC`（浅灰） |
| 标签字体 | 底部工具栏 → 字体 | Arial Italic, 32 |
| 标签显示 | 底部工具栏 → T 按钮 | 开启（Label 列已预处理，只有 keystone 有值） |

### 4. 导出

- **File → Export → SVG/PDF/PNG**
- 推荐 **SVG**（矢量，论文投稿后可在 Illustrator 编辑）
- PNG 分辨率 ≥ 300 dpi

### 5. 网络统计参数参考

本项目实测数据：
- 节点 130 / 边 311 / 阈值 |Spearman r| > 0.9, p < 0.05
- Modularity class 8 组
- Keystone 14 个（Degree ≥ 10 或 Betweenness ≥ 200）
""")
    else:
        st.info(f"「{analysis_type}」模块开发中 — Phase 2 继续。")

elif page == "生物地球化学循环图":
    st.title("生物地球化学循环图生成器")
    st.caption("Phase 3 v1 — 从 MAG × KO + 环境因子自动推断活跃元素循环通路")
    _render_interpretation_expander("cycle_diagram", key="interp_cycle_diagram")

    file_names = list(st.session_state.files.keys())
    if not file_names:
        st.warning("请先在「文件管理」上传 KO 注释表（MAG + KEGG_ko 长表）。")
        st.stop()

    # 基于文件识别类型预选最合适的文件（用户可再手动改）
    def _first_of(*ftypes: FileType) -> str | None:
        for n, info in st.session_state.files.items():
            if info["type"] in ftypes:
                return n
        return None

    def _idx_or_default(options: list[str], name: str | None) -> int:
        if name and name in options:
            return options.index(name)
        return 0

    ko_default = _first_of(FileType.KO_ANNOTATION_LONG)
    ko_options = file_names
    ko_name = st.selectbox(
        "KO 注释表（必需）", ko_options, key="cy_ko",
        index=_idx_or_default(ko_options, ko_default),
    )

    t_default = _first_of(FileType.MAG_TAXONOMY)
    t_options = ["（无）"] + [n for n in file_names if n != ko_name]
    t_name = st.selectbox(
        "GTDB 分类表（可选，用于 Phylum/Genus 标签）",
        t_options, key="cy_tax",
        index=_idx_or_default(t_options, t_default),
    )

    k_default = _first_of(FileType.KEYSTONE_SPECIES)
    k_options = ["（无）"] + [n for n in file_names if n not in (ko_name, t_name)]
    k_name = st.selectbox(
        "Keystone 物种列表（可选）", k_options, key="cy_ks",
        index=_idx_or_default(k_options, k_default),
    )

    # 循环图丰度加权是 MAG 级语义 —— 避开 Genus.txt/Phylum.txt/Species.txt
    # （detector 同样标为 ABUNDANCE_WIDE 但 conf=0.88，MAG 级 conf=0.95）
    a_default = _first_mag_abundance()
    a_options = ["（无）"] + [n for n in file_names
                             if n not in (ko_name, t_name, k_name)]
    a_name = st.selectbox(
        "MAG 丰度表（可选，用于通路贡献加权）",
        a_options, key="cy_ab",
        index=_idx_or_default(a_options, a_default),
    )

    e_default = _first_of(FileType.ENV_FACTORS)
    e_options = ["（无）"] + [n for n in file_names
                             if n not in (ko_name, t_name, k_name, a_name)]
    e_name = st.selectbox(
        "环境因子表（可选，用于 env-pathway 相关性）",
        e_options, key="cy_env",
        index=_idx_or_default(e_options, e_default),
    )

    m_default = _first_of(FileType.METADATA)
    m_options = ["（无）"] + [n for n in file_names
                             if n not in (ko_name, t_name, k_name, a_name, e_name)]
    m_name = st.selectbox(
        "Metadata（env 表需要搭配）",
        m_options, key="cy_md",
        index=_idx_or_default(m_options, m_default),
    )

    # 组选择下拉（S2.5-4）— 读取 metadata 的 Group 列可用值
    group_options = ["All"]
    if m_name != "（无）":
        md_df = st.session_state.files[m_name]["df"]
        if "Group" in md_df.columns:
            group_options = ["All"] + sorted(
                set(str(g) for g in md_df["Group"].dropna().unique()))
    group_filter = st.selectbox(
        "组选择（选单组将只用该组样本做通路贡献和 env 相关）",
        group_options, key="cy_group",
        help="选 'All' 用全部样本；选 'CK' / 'A' / 'B' 则仅用该组。"
             "用于生成单组循环图，三张单组图可手动拼接做组对比。",
    )

    with st.sidebar:
        st.subheader("循环图参数")
        comp_thresh = st.slider("通路完整度阈值（%）", 0, 100, 50, step=5,
                                key="cy_comp")
        top_n = st.slider("每通路显示 Top MAG 数", 1, 10, 5, key="cy_topn")
        rho_min = st.slider("env-pathway |ρ| 阈值", 0.0, 1.0, 0.5, step=0.05,
                            key="cy_rho")
        p_max = st.slider("env-pathway p 阈值", 0.001, 0.2, 0.05, step=0.005,
                          key="cy_p")
        max_cells = st.slider("每元素最多细胞数", 1, 6, 3, key="cy_cells")
        show_env = st.checkbox("显示环境耦合面板", True, key="cy_env_panel")
        show_couplings = st.checkbox("显示化学物耦合线", True, key="cy_coup")
        hide_regulator = st.checkbox(
            "🧪 隐藏纯调控型 cell（fur/tonB 等无底物产物）",
            False, key="cy_hide_reg",
            help=(
                "调控型基因（如 fur / tonB / arsR 转录因子）不催化底物→产物反应，"
                "KEGG 里 substrate/product 为空。勾选后隐藏所有基因都为调控型"
                "的 cell，使画面聚焦于催化型通路。"
            ),
        )
        annotate_cg = st.checkbox(
            "🏆 标注跨组最活 ★ + keystone ⭐（仅单组模式有效）",
            False, key="cy_annot_cg",
            help="选 CK/A/B 单组后勾选；会额外跑 3 次推断，耗时稍增。",
        )
        ranking = st.radio(
            "MAG 选择判据",
            options=["abundance", "completeness",
                     "keystone_only", "keystone_priority"],
            index=0, key="cy_ranking",
            help=(
                "• abundance（默认）: 完整度 × log(丰度)，偏重占主导者\n"
                "• completeness: 纯通路覆盖深度，忽略丰度\n"
                "• keystone_only: 只看关键物种；某通路无 keystone 则省略 cell\n"
                "• keystone_priority: keystone 10× 加权软优先，非 keystone 仍可竞争"
            ),
        )
        size = render_figure_size({"width_mm": 460, "height_mm": 320},
                                  prefix="cy")

    if st.button("生成循环图", type="primary", key="cy_go"):
        def _get(n):
            return st.session_state.files[n]["df"] if n != "（无）" else None
        ko = st.session_state.files[ko_name]["df"]
        t = _get(t_name)
        if t is not None and t.shape[1] == 2:
            cols_lower = [c.lower() for c in t.columns]
            has_header = ("classification" in cols_lower
                          or "taxonomy" in cols_lower
                          or any("mag" in c for c in cols_lower))
            if not has_header:
                # 首行是 MAG id + 分类字符串（无 header）→ 作为数据行恢复
                import pandas as _pd
                original_header = list(t.columns)
                t = _pd.concat(
                    [_pd.DataFrame([original_header], columns=["MAG", "Taxonomy"]),
                     t.rename(columns={t.columns[0]: "MAG",
                                       t.columns[1]: "Taxonomy"})],
                    ignore_index=True,
                )
        params = {
            "completeness_threshold": float(comp_thresh),
            "top_n_contributors": top_n,
            "env_rho_min": rho_min,
            "env_p_max": p_max,
            "show_env_panel": show_env,
            "show_couplings": show_couplings,
            "max_cells_per_element": max_cells,
            "group_filter": None if group_filter == "All" else group_filter,
            "annotate_cross_group": annotate_cg,
            "contributor_ranking": ranking,
            "hide_regulator_only_cells": hide_regulator,
            **size,
        }
        try:
            result = cycle_mod.analyze(
                ko, t, _get(k_name), _get(a_name), _get(e_name), _get(m_name),
                params=params,
            )
        except ValueError as e:
            st.error(str(e))
            st.stop()
        st.session_state["_cy_last"] = result
        # 登记循环图导出信息 — file_paths 用当前选中的文件
        _cy_files = {"ko_annotation": ko_name}
        if t_name != "（无）": _cy_files["taxonomy"] = t_name
        if k_name != "（无）": _cy_files["keystone"] = k_name
        if a_name != "（无）": _cy_files["abundance"] = a_name
        if e_name != "（无）": _cy_files["env_factors"] = e_name
        if m_name != "（无）": _cy_files["metadata"] = m_name
        _register_export("cycle_diagram", result=result,
                         file_paths=_cy_files, params=params)

        # Q3: per-group cycles（给 HTML 交互导出用于组切换）
        # 仅当 metadata 含 Group 列 + 至少 2 组时计算
        _md = _get(m_name)
        _per_group = {}
        if (_md is not None and not _md.empty
                and "Group" in _md.columns):
            try:
                from envmeta.geocycle.inference import infer as _infer_one
                _groups = sorted(set(_md["Group"].astype(str)))
                # 软 warning：>10 组时量化提示体积 / 导出时间（不阻止）
                if len(_groups) > 10:
                    st.warning(
                        f"检测到 **{len(_groups)} 组**。per_group_cycles 预计"
                        f"增加 HTML 体积 ~{len(_groups) * 100} KB / 导出时间 "
                        f"+{len(_groups) * 2}s。环境微生物研究通常 ≤ 10 组，"
                        f"建议按生态梯度 / 时间点分批导出 HTML。"
                    )
                if len(_groups) >= 2:
                    _inf_params = {k: v for k, v in params.items()
                                   if k != "group_filter"}
                    for _g in _groups:
                        _p = {**_inf_params, "group_filter": _g}
                        _cd_g = _infer_one(
                            ko, t, _get(k_name), _get(a_name),
                            _get(e_name), _md, params=_p,
                        )
                        _per_group[_g] = _cd_g
            except Exception as _e:  # noqa: BLE001
                st.caption(f"⚠️ 按组跑 infer 失败（不影响主循环图）：{_e}")
        st.session_state["_cy_per_group_last"] = _per_group or None

        # HTML 导出：全样本 CycleData（env / full_corr / sensitivity 用它，
        # 不受循环图 group_filter 影响）。仅当用户当前参数过滤到单组时才
        # 额外跑一次；否则 last 就是全样本
        _full_sample_cd = None
        _cur_gf = params.get("group_filter")
        if _cur_gf and str(_cur_gf).lower() not in ("none", "all", ""):
            try:
                from envmeta.geocycle.inference import (
                    infer as _infer_full_one,
                )
                _full_p = {**params, "group_filter": None}
                _full_sample_cd = _infer_full_one(
                    ko, t, _get(k_name), _get(a_name),
                    _get(e_name), _md, params=_full_p,
                )
            except Exception as _e:  # noqa: BLE001
                st.caption(
                    f"⚠️ 全样本 baseline infer 失败（env 面板将退化为"
                    f"当前过滤视图）：{_e}"
                )
        st.session_state["_cy_full_sample_last"] = _full_sample_cd

    last = st.session_state.get("_cy_last")
    if last is not None:
        st.pyplot(last.figure, use_container_width=False)
        d1, d2, d3 = st.columns(3)
        with d1:
            st.download_button("⬇️ PNG", data=export_to_bytes(last.figure, "png"),
                               file_name="cycle.png", mime="image/png", key="cy_png")
        with d2:
            st.download_button("⬇️ PDF（矢量）", data=export_to_bytes(last.figure, "pdf"),
                               file_name="cycle.pdf", mime="application/pdf", key="cy_pdf")
        with d3:
            st.download_button("⬇️ 统计表（TSV）",
                               data=last.stats.to_csv(sep="\t", index=False).encode("utf-8"),
                               file_name="cycle_stats.tsv",
                               mime="text/tab-separated-values", key="cy_tsv")
        _vector_downloads(last.figure, "cycle", "cy")

        # —— S4.5 独立交互 HTML 导出（T2 α/β/γ/δ/ε 全套）—————
        try:
            from envmeta.geocycle.html_exporter import build_interactive_html as _build_html
            _cycle_data = getattr(last, "data", None)
            if _cycle_data is None:
                raise RuntimeError(
                    "AnalysisResult 未携带 CycleData —— 请重新点击『生成循环图』"
                )
            _hyp = st.session_state.get("_hyp_last")
            _cmp = st.session_state.get("_cy_compare_last")
            _hyp_groups = st.session_state.get("_hyp_multi_last")  # DataFrame or None
            # Q3: per-group 切换所需的每组 CycleData
            _per_group = st.session_state.get("_cy_per_group_last")
            # env / full_corr / sensitivity 面板的全样本 baseline
            _full_sample = st.session_state.get("_cy_full_sample_last")

            _html_bytes = _build_html(
                _cycle_data,
                hypothesis=_hyp,
                compare_df=_cmp,
                hypothesis_by_group=_hyp_groups,
                per_group_cycles=_per_group,
                full_sample_cycle_data=_full_sample,
            )
            st.download_button(
                "📦 导出交互 HTML（独立 SI · D3.js 嵌入）",
                data=_html_bytes,
                file_name="envmeta_cycle_interactive.html",
                mime="text/html",
                key="cy_html_export",
                help=(
                    f"独立 HTML 文件（~{len(_html_bytes) // 1024} KB）含 D3.js 嵌入。"
                    "双击浏览器打开即可离线交互：循环图（4 象限 + 化学物节点 +"
                    "组切换 + 拖拽 + 缩放）+ 假说评分（跨组优先）+ 跨组对比 +"
                    "环境相关 5 档可信度 + SVG / JSON 导出。"
                ),
            )
            # Q5: 状态透明面板 — 让用户一眼看到 HTML 里包含/缺失哪些数据
            def _chk(v):
                return "✅" if v is not None and v is not False else "❌"
            _n_groups = (len(_per_group) if isinstance(_per_group, dict) else 0) or 0
            _hyp_multi_n = 0
            if _hyp_groups is not None:
                try: _hyp_multi_n = len(_hyp_groups)
                except Exception: _hyp_multi_n = 0
            # env 面板数据源提示（full_sample 生效时展示循环图过滤 vs env 全样本）
            _cycle_n = _cycle_data.meta.get(
                "n_samples_used", _cycle_data.meta.get("n_samples", "?"))
            _cur_gf = _cycle_data.params.get("group_filter")
            if _full_sample is not None:
                _full_n = _full_sample.meta.get(
                    "n_samples_used", _full_sample.meta.get("n_samples", "?"))
                _env_src_note = (
                    f"🌐 env 基于全样本 n={_full_n}（循环图视图 "
                    f"group={_cur_gf} n={_cycle_n}）"
                )
            else:
                _env_src_note = f"🌐 env 基于全样本 n={_cycle_n}（循环图视图未过滤）"
            st.caption(
                f"📋 **HTML 已注入数据**：  "
                f"{_chk(_cycle_data)} 循环图视图  ·  "
                f"{_chk(_per_group)} per-group 循环（{_n_groups} 组）  ·  "
                f"{_chk(_hyp)} 单组假说  ·  "
                f"{_chk(_hyp_groups)} 跨组假说（{_hyp_multi_n} 行）  ·  "
                f"{_chk(_cmp)} 通路×组对比  ·  "
                f"{_env_src_note}"
            )
            _missing = []
            if _hyp is None:
                _missing.append("**单组假说** → 到「🧪 假说评分」页跑")
            if _hyp_groups is None:
                _missing.append("**跨组假说** → 「🧪 假说评分」页 ✓ 勾选「跨组对比」后再跑")
            if _cmp is None:
                _missing.append("**通路×组对比** → 展开下方「跨组对比表（回答...）」expander 并跑")
            if _per_group is None and _n_groups == 0:
                _missing.append("**per-group 循环切换** → metadata 需有 Group 列 + ≥ 2 组")
            if _missing:
                st.caption("ℹ️ HTML 里空白的模块是因为以下数据未生成：  \n- " +
                           "  \n- ".join(_missing))
        except Exception as _e:  # noqa: BLE001
            import traceback as _tb
            st.caption(f"⚠️ 交互 HTML 导出失败：{_e}")
            with st.expander("详细错误追踪（调试用）", expanded=False):
                st.code(_tb.format_exc())

        # —— 跨组对比小工具（S2.5-7d）——————————————————
        with st.expander("跨组对比表（回答『A/CK/B 差在哪里』）", expanded=False):
            st.caption(
                "同一通路在不同组的 `total_contribution` 反映丰度加权活性差异；"
                "`top_mag_genus` 变化反映承载者身份差异；"
                "`top_mag_is_keystone` 标出承载者是否为关键物种。"
            )
            if st.button("生成跨组对比表", key="cy_compare_go"):
                from envmeta.analysis.cycle_compare import compare_groups as _cg

                def _get_df(n):
                    return (st.session_state.files[n]["df"]
                            if n != "（无）" else None)

                ko_df = st.session_state.files[ko_name]["df"]
                t_df = _get_df(t_name)
                if t_df is not None and t_df.shape[1] == 2:
                    cols_lower = [c.lower() for c in t_df.columns]
                    has_header = ("classification" in cols_lower
                                  or "taxonomy" in cols_lower
                                  or any("mag" in c for c in cols_lower))
                    if not has_header:
                        import pandas as _pd
                        original_header = list(t_df.columns)
                        t_df = _pd.concat(
                            [_pd.DataFrame([original_header],
                                           columns=["MAG", "Taxonomy"]),
                             t_df.rename(columns={t_df.columns[0]: "MAG",
                                                  t_df.columns[1]: "Taxonomy"})],
                            ignore_index=True,
                        )
                try:
                    compare_df = _cg(
                        ko_df, t_df, _get_df(k_name), _get_df(a_name),
                        _get_df(e_name), _get_df(m_name),
                        params={"completeness_threshold": float(comp_thresh),
                                "top_n_contributors": top_n,
                                "env_rho_min": rho_min,
                                "env_p_max": p_max},
                    )
                    st.session_state["_cy_compare_last"] = compare_df
                except Exception as e:
                    st.error(f"跨组对比失败：{e}")
            cmp_last = st.session_state.get("_cy_compare_last")
            if cmp_last is not None:
                st.dataframe(cmp_last, use_container_width=True)
                st.download_button(
                    "⬇️ 跨组对比表（TSV）",
                    data=cmp_last.to_csv(sep="\t", index=False).encode("utf-8"),
                    file_name="cycle_compare.tsv",
                    mime="text/tab-separated-values",
                    key="cy_cmp_tsv",
                )

        st.info(
            "⚠️ **输出是描述性的，不是因果性的。** "
            "Top-completeness contributor 只表示该 MAG 的 KO 覆盖该通路最多，"
            "**不等于驱动**该通路。相关 ≠ 因果。因果解读需要领域专家基于证据自行判断。"
        )
        with st.expander("活跃通路 + 超阈值环境耦合（主结果）"):
            st.dataframe(
                last.stats[last.stats["type"].isin(["pathway", "env_correlation"])],
                use_container_width=True,
            )
        with st.expander("完整环境相关矩阵（不过滤，用于判断 confounding）"):
            st.caption("同一通路同时与多个环境因子相关时，提示潜在共变。")
            st.dataframe(
                last.stats[last.stats["type"] == "full_correlation"],
                use_container_width=True,
            )
        with st.expander("阈值敏感度扫描（robust = Top-1 在 30/50/70% 三档一致）"):
            sens = last.stats[last.stats["type"] == "sensitivity"].copy()
            if not sens.empty:
                sens["robust"] = sens["mean_completeness"].apply(
                    lambda x: "✅ robust" if x >= 0.5 else "⚠️ threshold-sensitive")
                st.dataframe(sens[["element", "pathway_id", "display_name",
                                   "top_mag", "robust"]],
                             use_container_width=True)

        # —— 📦 Fork Bundle (论文复现包, S4) ——————————————————
        with st.expander("📦 Fork Bundle — 论文复现包", expanded=False):
            st.caption(
                "把当前 KB + 假说 YAML + 分析参数打包成单个 .zip，"
                "读者/审稿人一键加载即可复现本论文。"
                "论文-EnvMeta 绑定发布协议的落地实现（L4 层）。"
            )
            st.info(
                "💡 图表 / 复现脚本 / 文档的批量下载请到「**导出中心**」（左侧导航）。"
            )
            bcols = st.columns(2)

            # ── 加载 Bundle（左列）
            with bcols[0]:
                st.markdown("**📂 加载 Bundle**")
                bundle_upload = st.file_uploader(
                    "选择 .zip Bundle", type=None, key="bundle_upload",
                    help="上传别人发布的 Bundle → 自动解包 KB + 假说",
                )
                if bundle_upload is not None and st.button(
                        "加载 Bundle", key="bundle_load_go"):
                    from envmeta.tools.bundle import (
                        inspect_bundle as _inspect_bundle,
                        load_bundle as _load_bundle,
                    )
                    import tempfile
                    with tempfile.NamedTemporaryFile(
                            suffix=".zip", delete=False) as tf:
                        tf.write(bundle_upload.read())
                        tpath = Path(tf.name)
                    try:
                        info = _inspect_bundle(tpath)
                        b = _load_bundle(tpath)
                        st.session_state["_bundle_last"] = {
                            "info": info,
                            "hypothesis_texts": b.hypothesis_texts,
                            "config": b.config,
                            "manifest": b.manifest,
                            "readme": b.readme,
                        }
                        st.success(
                            f"已加载 {info['manifest'].get('name', '?')}"
                        )
                    except Exception as e:  # noqa: BLE001
                        st.error(f"Bundle 加载失败：{e}")
                    finally:
                        tpath.unlink(missing_ok=True)

            # ── 导出 Bundle（右列）
            with bcols[1]:
                st.markdown("**⬇️ 导出当前状态为 Bundle**")
                bundle_name = st.text_input(
                    "Bundle 名", value="my_paper", key="bundle_name",
                )
                bundle_author = st.text_input(
                    "作者", value="", key="bundle_author",
                )
                bundle_doi = st.text_input(
                    "论文 DOI（可选）", value="", key="bundle_doi",
                )
                # 可选内容
                _kegg_default_path = Path(
                    "envmeta/geocycle/kegg_snapshot.json"
                )
                bundle_include_kegg = st.checkbox(
                    "包含 KEGG 快照（溯源，推荐）",
                    value=_kegg_default_path.exists(),
                    key="bundle_include_kegg",
                    help="把 envmeta/geocycle/kegg_snapshot.json "
                         "打进 kb/kegg_snapshot.json，"
                         "让读者能追溯 KB 的 KEGG 数据来源时间点。",
                    disabled=not _kegg_default_path.exists(),
                )
                bundle_include_readme = st.checkbox(
                    "自动生成 README.md",
                    value=True,
                    key="bundle_include_readme",
                    help="基于名字/作者/DOI 生成一份简短加载指引，"
                         "读者拿到 bundle 就知道怎么用。",
                )
                # 上传的假说预览
                _hyp_up_preview = st.session_state.get("hyp_upload")
                if _hyp_up_preview is not None:
                    st.caption(
                        f"📎 将包含假说 YAML: **{_hyp_up_preview.name}** "
                        "（在上方「🧪 假说评分」区已上传）"
                    )
                else:
                    st.caption(
                        "⚠️ 尚未在「🧪 假说评分」区上传 YAML — "
                        "当前 bundle 不会包含假说。"
                        "请先上传再点创建。"
                    )
                if st.button("创建 Bundle", key="bundle_create_go"):
                    from envmeta.tools.bundle import create_bundle as _create_bundle
                    import tempfile
                    # 收集当前状态：KB（内置）+ 假说（若已上传）+ 分析参数
                    kb_default = Path(
                        "envmeta/geocycle/knowledge_base/elements.json"
                    )
                    hyp_paths: list[Path] = []
                    # 把 session 里上传的 YAML 临时落地
                    hyp_up = st.session_state.get("hyp_upload")
                    if hyp_up is not None:
                        th = Path(tempfile.gettempdir()) / hyp_up.name
                        th.write_bytes(hyp_up.getvalue())
                        hyp_paths.append(th)
                    # 当前循环图参数
                    current_cfg = {
                        "completeness_threshold": float(comp_thresh),
                        "top_n_contributors": top_n,
                        "env_rho_min": rho_min,
                        "env_p_max": p_max,
                        "contributor_ranking": ranking,
                        "max_cells_per_element": max_cells,
                        "hide_regulator_only_cells": hide_regulator,
                    }
                    # 自动 README（若勾选）
                    readme_text = None
                    if bundle_include_readme:
                        _hyp_list = ", ".join(p.name for p in hyp_paths) or "（无）"
                        readme_text = (
                            f"# {bundle_name}\n\n"
                            + (f"**作者**: {bundle_author}\n\n" if bundle_author else "")
                            + (f"**DOI**: {bundle_doi}\n\n" if bundle_doi else "")
                            + "## 如何使用本 Bundle\n\n"
                            "1. 启动 EnvMeta：`streamlit run app.py`\n"
                            "2. 循环图页 → 展开「📦 Fork Bundle」→ 上传本 zip\n"
                            "3. 准备原始数据（见论文 Data Availability）\n"
                            "4. 生成循环图 + 假说评分\n\n"
                            "## 内容\n\n"
                            "- `kb/elements.json`: 元素循环知识库\n"
                            + ("- `kb/kegg_snapshot.json`: KEGG 来源快照\n"
                               if bundle_include_kegg and _kegg_default_path.exists() else "")
                            + f"- `hypotheses/`: 机制假说 YAML（{_hyp_list}）\n"
                            "- `config/cycle_params.yaml`: 分析参数\n\n"
                            "## 更多信息\n\n"
                            "见 [EnvMeta 主仓库](https://github.com/redlizzxy/EnvMeta)"
                            " 的 `paper/bundles/README.md`。\n"
                        )
                    kegg_path = (
                        _kegg_default_path
                        if bundle_include_kegg and _kegg_default_path.exists()
                        else None
                    )
                    out_tmp = Path(tempfile.gettempdir()) / f"{bundle_name}.zip"
                    try:
                        _create_bundle(
                            out_tmp,
                            kb_path=kb_default,
                            hypothesis_paths=hyp_paths,
                            config=current_cfg,
                            kegg_snapshot_path=kegg_path,
                            readme_text=readme_text,
                            name=bundle_name,
                            author=bundle_author,
                            paper_doi=bundle_doi,
                            description="Created via EnvMeta app UI",
                        )
                        st.session_state["_bundle_out_bytes"] = out_tmp.read_bytes()
                        st.session_state["_bundle_out_name"] = f"{bundle_name}.zip"
                        # 构造文件清单告知用户
                        manifest_items = [
                            "✅ manifest.yaml",
                            "✅ kb/elements.json",
                        ]
                        if kegg_path:
                            manifest_items.append("✅ kb/kegg_snapshot.json")
                        if hyp_paths:
                            for p in hyp_paths:
                                manifest_items.append(f"✅ hypotheses/{p.name}")
                        else:
                            manifest_items.append("⚠️ hypotheses/ (空)")
                        manifest_items.append("✅ config/cycle_params.yaml")
                        if readme_text:
                            manifest_items.append("✅ README.md")
                        st.success(
                            f"Bundle 创建成功（{out_tmp.stat().st_size} bytes）"
                            + "\n\n" + "\n".join(manifest_items)
                        )
                    except Exception as e:  # noqa: BLE001
                        st.error(f"创建失败：{e}")
                if st.session_state.get("_bundle_out_bytes"):
                    st.download_button(
                        "⬇️ 下载 Bundle",
                        data=st.session_state["_bundle_out_bytes"],
                        file_name=st.session_state.get(
                            "_bundle_out_name", "bundle.zip",
                        ),
                        mime="application/zip",
                        key="bundle_download",
                    )

            # 展示已加载 bundle 的 summary
            _bundle_last = st.session_state.get("_bundle_last")
            if _bundle_last is not None:
                st.markdown("---\n**已加载 Bundle**")
                m = _bundle_last["manifest"]
                info = _bundle_last["info"]
                st.json({
                    "name": m.get("name"),
                    "author": m.get("author"),
                    "paper_doi": m.get("paper_doi"),
                    "created": m.get("created"),
                    "envmeta_version_bundle": m.get("envmeta_version"),
                    "envmeta_version_runtime": info["envmeta_version_runtime"],
                    "version_mismatch": info["version_mismatch"],
                    "hypotheses": [n for n, _ in _bundle_last["hypothesis_texts"]],
                    "config": _bundle_last["config"],
                })
                if info["version_mismatch"]:
                    st.warning(
                        "⚠️ Bundle 的 envmeta_version 与当前运行时不一致；"
                        "加载仍可用，但部分字段可能含未知扩展。"
                    )

        # —— 机制假说 YAML 评分器（S3）——————————————————
        with st.expander("🧪 假说评分 (可选) — 上传机制 YAML", expanded=False):
            st.caption(
                "上传一份机制假说 YAML，对照本次循环图推断结果评分。"
                "评分是**描述性证据加权**，不是因果证实。"
                "Schema 说明见 `paper/hypotheses/README.md`。"
            )
            _render_interpretation_expander("hypothesis_score", key="interp_hypothesis_score")
            _tmpl = (Path("paper") / "hypotheses" / "arsenic_steel_slag.yaml")
            if _tmpl.exists():
                st.download_button(
                    "⬇️ 下载示例 YAML (arsenic_steel_slag)",
                    data=_tmpl.read_bytes(),
                    file_name="arsenic_steel_slag.yaml",
                    mime="application/x-yaml",
                    key="hyp_tmpl",
                )
            # 接受任意类型：避免错传文件时 streamlit 的"not allowed" 悬浮条
            # 遮挡 × 删除键；格式错误在点"评分"时由 load_hypothesis 抛出
            hyp_file = st.file_uploader(
                "假说 YAML（.yaml / .yml）",
                type=None, key="hyp_upload",
                help="选错文件？点文件名后的 × 可直接移除。"
                     "格式错误会在点击「评分」后提示。",
            )
            hyp_multi_group = st.checkbox(
                "📊 对每组分别评分（生成跨组对比表）",
                value=False, key="hyp_multi_group",
                help="勾选后对 metadata.Group 里每个组分别评分，"
                     "生成 group × 指标 的对比表。"
                     "用于比较不同处理/条件下数据对同一假说的支持度差异，"
                     "是论文 Results 里跨组对比叙事的直接素材。"
                     "💡 **勾选后跨组评分会同步注入到 HTML 交互导出，"
                     "在假说评分 tab 顶部显示跨组卡片组。**",
            )
            if hyp_file is not None and st.button("评分", key="hyp_score_go"):
                try:
                    from envmeta.geocycle.hypothesis import (
                        load_hypothesis as _load_hyp,
                        score as _score_hyp,
                    )
                    from envmeta.geocycle.inference import infer as _infer

                    def _get_file(n):
                        return (st.session_state.files[n]["df"]
                                if n != "（无）" else None)

                    ko_df2 = st.session_state.files[ko_name]["df"]
                    t_df2 = _get_file(t_name)
                    if t_df2 is not None and t_df2.shape[1] == 2:
                        cols_lower = [c.lower() for c in t_df2.columns]
                        has_header = ("classification" in cols_lower
                                      or "taxonomy" in cols_lower
                                      or any("mag" in c for c in cols_lower))
                        if not has_header:
                            original_header = list(t_df2.columns)
                            t_df2 = pd.concat(
                                [pd.DataFrame([original_header],
                                              columns=["MAG", "Taxonomy"]),
                                 t_df2.rename(columns={t_df2.columns[0]: "MAG",
                                                       t_df2.columns[1]: "Taxonomy"})],
                                ignore_index=True,
                            )
                    hyp_params = {
                        "completeness_threshold": float(comp_thresh),
                        "env_rho_min": rho_min,
                        "env_p_max": p_max,
                        "group_filter": None if group_filter == "All" else group_filter,
                    }
                    data = _infer(
                        ko_df2, t_df2, _get_file(k_name), _get_file(a_name),
                        _get_file(e_name), _get_file(m_name),
                        params=hyp_params,
                    )
                    hyp_obj = _load_hyp(hyp_file.read().decode("utf-8"))
                    # S3.5: 若 YAML 含 group_contrast claim，自动调
                    # compare_groups 注入 compare_df
                    has_gc = any(c.type == "group_contrast"
                                  for c in hyp_obj.claims)
                    compare_df_arg = None
                    if has_gc:
                        from envmeta.analysis.cycle_compare import (
                            compare_groups as _compare_groups,
                        )
                        cmp_params = dict(hyp_params)
                        cmp_params.pop("group_filter", None)  # compare 遍历所有组
                        try:
                            compare_df_arg = _compare_groups(
                                ko_df2, t_df2, _get_file(k_name),
                                _get_file(a_name), _get_file(e_name),
                                _get_file(m_name), params=cmp_params,
                            )
                        except Exception as ge:  # noqa: BLE001
                            st.warning(
                                f"自动跨组对比失败，group_contrast "
                                f"将 skipped：{ge}"
                            )
                    hyp_result = _score_hyp(
                        hyp_obj, data, compare_df=compare_df_arg,
                    )
                    st.session_state["_hyp_last"] = hyp_result
                    # 跨组对比（可选）
                    if hyp_multi_group:
                        from envmeta.analysis.hypothesis_compare import (
                            score_by_groups as _score_by_groups,
                        )
                        cmp_params = dict(hyp_params)
                        cmp_params.pop("group_filter", None)
                        try:
                            multi_df = _score_by_groups(
                                hyp_obj,
                                ko_df2, t_df2, _get_file(k_name),
                                _get_file(a_name), _get_file(e_name),
                                _get_file(m_name),
                                params=cmp_params, null_n=299,
                            )
                            st.session_state["_hyp_multi_last"] = multi_df
                        except Exception as me:  # noqa: BLE001
                            st.warning(f"跨组评分失败：{me}")
                            st.session_state["_hyp_multi_last"] = None
                    else:
                        st.session_state["_hyp_multi_last"] = None
                except Exception as e:
                    st.error(f"评分失败：{e}")
            hyp_last = st.session_state.get("_hyp_last")
            if hyp_last is not None:
                _color = {
                    "strong": "#27AE60",
                    "suggestive": "#F39C12",
                    "weak": "#95A5A6",
                    "insufficient": "#E74C3C",
                }.get(hyp_last.label, "#888")
                st.markdown(
                    f"### {hyp_last.hypothesis_name}  "
                    f"<span style='background:{_color};color:white;"
                    f"padding:4px 10px;border-radius:6px;font-size:0.85em;'>"
                    f"{hyp_last.label}</span>  "
                    f"**{hyp_last.overall_score:.2f}**",
                    unsafe_allow_html=True,
                )
                st.caption(
                    f"{hyp_last.n_satisfied}/{hyp_last.n_total} claim 满足"
                    f"（{hyp_last.n_skipped} 条 skipped 不计）"
                )
                # S3.5-ui: 综合解读一句话
                _interp_txt, _interp_type = _interpret_hyp_score(hyp_last)
                getattr(st, _interp_type)(_interp_txt)

                # S3.5-ui 扩展：跨组对比表（若勾选）
                _multi_df = st.session_state.get("_hyp_multi_last")
                _has_multi = _multi_df is not None and not _multi_df.empty
                if _has_multi:
                    st.markdown("#### 📊 跨组支持度对比")
                    _display_df = _multi_df.copy()
                    _display_df["null_p"] = _display_df["null_p"].apply(
                        lambda x: "N/A" if pd.isna(x) else f"{x:.3f}"
                    )
                    _display_df["weight_robust"] = _display_df["weight_robust"].map(
                        {True: "✅", False: "⚠️"}
                    ).fillna("N/A")
                    st.dataframe(
                        _display_df[[
                            "group", "overall_score", "label",
                            "null_p", "weight_robust",
                            "n_satisfied", "n_total", "n_veto",
                            "interpretation",
                        ]],
                        use_container_width=True,
                    )
                    st.download_button(
                        "⬇️ 跨组评分表 (TSV)",
                        data=_multi_df.to_csv(sep="\t", index=False)
                              .encode("utf-8"),
                        file_name="hypothesis_score_by_group.tsv",
                        mime="text/tab-separated-values",
                        key="hyp_multi_tsv",
                    )

                # S3.5: 显示 veto / null_p / weight_robust 三项可信度指标
                if hyp_last.veto_reasons:
                    base_label = hyp_last.params.get(
                        "base_label_before_veto", "?",
                    )
                    st.error(
                        "🚫 **VETOED** — 由 required=true 的 claim 触发硬否决；"
                        f"若无 veto 原本应为 **{base_label}** "
                        f"({hyp_last.overall_score:.2f})。"
                        + "\n\n触发原因:\n"
                        + "\n".join(f"- {r}" for r in hyp_last.veto_reasons)
                    )
                # 勾选跨组对比时，对比表已含当前组的 null_p / robust，
                # 无需再单独展示（减冗余）
                if not _has_multi:
                    mcols = st.columns(2)
                    with mcols[0]:
                        if hyp_last.null_p is not None:
                            if hyp_last.null_p < 0.05:
                                st.success(
                                    f"**null_p = {hyp_last.null_p:.3f}** "
                                    f"(n={hyp_last.null_p_samples})  \n"
                                    "📉 越小越好；此处 < 0.05 = "
                                    "权重设计与数据支持显著一致（特异）"
                                )
                            elif hyp_last.null_p < 0.20:
                                st.warning(
                                    f"**null_p = {hyp_last.null_p:.3f}** "
                                    f"(n={hyp_last.null_p_samples})  \n"
                                    "📉 越小越好；0.05-0.20 = 边界，"
                                    "建议增加独立证据"
                                )
                            else:
                                st.info(
                                    f"**null_p = {hyp_last.null_p:.3f}** "
                                    f"(n={hyp_last.null_p_samples})  \n"
                                    "📉 ≥ 0.20 = 通过率运气主导，"
                                    "权重设计未被数据特异支持"
                                )
                        else:
                            # 区分退化式 N/A（好）vs 无信号 N/A（不好）
                            if hyp_last.overall_score >= 0.95:
                                st.success(
                                    "**null_p = N/A**（退化式）  \n"
                                    "🎯 全 claim satisfied → 排列无意义；"
                                    "此 N/A 是『好』的退化，看 robust 兜底"
                                )
                            else:
                                st.info(
                                    "**null_p = N/A**（无信号）  \n"
                                    "claim 数 < 3 / weight 同质 / score 同质 "
                                    "→ 排列统计力不足，**不是『最好』**"
                                )
                    with mcols[1]:
                        if hyp_last.weight_robust is True:
                            st.success(
                                "✅ **weight robust**  \n"
                                "±20% 权重扰动下 label 保持不变，"
                                "结论不靠阈值挑权重"
                            )
                        elif hyp_last.weight_robust is False:
                            n_flip = sum(
                                1 for r in hyp_last.weight_sensitivity_rows
                                if r.get("flipped")
                            )
                            st.warning(
                                f"⚠️ **weight sensitive**  \n"
                                f"{n_flip}/{len(hyp_last.weight_sensitivity_rows)} "
                                f"个扰动下 label 翻转，建议审视权重分配"
                            )
                        else:
                            st.info(
                                "weight sensitivity = N/A  \n"
                                "claim 数不足，非稳健也非敏感"
                            )

                st.dataframe(hyp_last.to_dataframe(), use_container_width=True)

                if hyp_last.weight_sensitivity_rows:
                    with st.expander("weight sensitivity 详表（OAT ±20%）"):
                        st.dataframe(
                            pd.DataFrame(hyp_last.weight_sensitivity_rows),
                            use_container_width=True,
                        )

                with st.expander("📖 三指标怎么看（独立指标 + 组合判断）"):
                    st.markdown(
                        """
**1. overall_score** (0 ~ 1.0) — 越高越好

加权平均"claim 被数据支持的程度"。1.0 = 全通过。默认阈值：
- ≥ `strong_threshold` (default 0.75) → label = **strong**
- ≥ `suggestive_threshold` (default 0.40) → **suggestive**
- \> 0 → **weak**；= 0 或全 skipped → **insufficient**

---

**2. null_p** (排列检验 p 值) — **越低越好**

把观测到的 score 集合在 claim 间随机重新分配 999 次，
算 P(null_overall ≥ observed)。

| null_p | 含义 |
|---|---|
| < 0.05 | ⭐ 权重设计与数据特异一致（strong + specific）|
| 0.05-0.20 | 🟢 边界显著 |
| ≥ 0.20 | ⚠️ 通过率运气主导，权重设计未被数据特异支持 |
| **N/A (退化式)** | 🎯 overall ≥ 0.95，全通过导致排列无意义。**此 N/A 是"好的"** |
| **N/A (无信号)** | ⚠️ claim<3 / weight 同质 / score 同质。**不是"好"，是统计力不足** |

---

**3. weight_robust** (OAT ±20%) — `True` 最好

对每条 claim 独立扰动 weight ±20%，检查 label 是否翻转。
- `True` = label 稳健，结论不靠挑权重
- `False` = 对 ±20% 扰动敏感，审视权重依据
- `None` = claim 数 < 2，无法测

---

**4. 🚫 VETOED** (required claim)

任一 `required: true` 的 claim 失败 → label 强制 `insufficient`。
overall_score 仍显示（透明度），但 label 被否决。

---

**组合判断优先级**（顶部"一句话解读"自动应用这张表）：

| 条件 | 判定 |
|---|---|
| vetoed | 🚫 假说未通过（必要前提失败）|
| strong + (null_p<0.05 或 退化式 N/A) + robust | ⭐ 最强支持 |
| strong + robust + null_p≥0.20 | ⚠️ strong 但不特异 |
| strong + not robust | ⚠️ strong 但权重敏感 |
| strong + 0.05≤null_p<0.20 | 🟢 strong 边界 |
| suggestive + specific | 🟡 中等支持（权重一致）|
| weak | 🟠 证据薄弱 |
| insufficient (未 veto) | 🔴 证据不足 |
""".strip()
                    )

                c1, c2 = st.columns(2)
                with c1:
                    st.download_button(
                        "⬇️ 评分报告 (TSV)",
                        data=hyp_last.to_dataframe()
                              .to_csv(sep="\t", index=False).encode("utf-8"),
                        file_name="hypothesis_score.tsv",
                        mime="text/tab-separated-values",
                        key="hyp_tsv",
                    )
                with c2:
                    st.download_button(
                        "⬇️ 评分报告 (JSON)",
                        data=hyp_last.to_json().encode("utf-8"),
                        file_name="hypothesis_score.json",
                        mime="application/json",
                        key="hyp_json",
                    )

elif page == "导出中心":
    st.title("📦 导出中心")
    st.caption("所有生成物集中导出 · 图表 / Bundle / 复现脚本 / 文档 四合一")

    registry = st.session_state.get("_export_registry", {})

    tab_fig, tab_bundle, tab_code, tab_doc = st.tabs([
        f"📊 图表 ({len(registry)})",
        "📦 Fork Bundle",
        f"🐍 复现脚本 ({len(registry)})",
        "📚 文档下载",
    ])

    # ── tab 1: 图表 ─────────────────────────────────────────
    with tab_fig:
        if not registry:
            st.info(
                "🤔 还没有已生成的图表。\n\n"
                "请先到**任意分析页**点「🎨 生成图表」按钮，然后回到这里批量下载。\n\n"
                "💡 提示：示例数据可在「首页」一键加载。"
            )
        else:
            st.markdown(f"**共 {len(registry)} 个已生成图表**")
            # 批量 zip 下载（T1.3）
            import io as _io
            import zipfile as _zipfile

            def _build_figures_zip() -> bytes:
                buf = _io.BytesIO()
                with _zipfile.ZipFile(buf, "w", _zipfile.ZIP_DEFLATED) as zf:
                    for _aid, _entry in registry.items():
                        _base = _entry["output_base"]
                        _fig = _entry["result"].figure
                        # PNG + PDF + SVG + TIFF 4 格式
                        for _fmt, _ext in [("png", "png"), ("pdf", "pdf"),
                                           ("svg", "svg"), ("tiff", "tiff")]:
                            try:
                                zf.writestr(
                                    f"{_base}/{_base}.{_ext}",
                                    export_to_bytes(_fig, _fmt),
                                )
                            except Exception:  # noqa: BLE001
                                pass
                        # stats TSV
                        _stats = getattr(_entry["result"], "stats", None)
                        if _stats is not None:
                            try:
                                if isinstance(_stats, pd.DataFrame):
                                    zf.writestr(
                                        f"{_base}/{_base}_stats.tsv",
                                        _stats.to_csv(sep="\t", index=False).encode("utf-8"),
                                    )
                            except Exception:  # noqa: BLE001
                                pass
                return buf.getvalue()

            st.download_button(
                "📦 一键批量下载（ZIP — 所有图 × PNG/PDF/SVG/TIFF + TSV）",
                data=_build_figures_zip(),
                file_name="envmeta_figures_bundle.zip",
                mime="application/zip",
                key="export_all_zip",
                type="primary",
            )
            st.markdown("---")
            # 逐条展示
            for aid, entry in registry.items():
                name = ANALYSIS_INPUTS.get(aid, {}).get("name", aid)
                base = entry["output_base"]
                with st.expander(f"📊 {name}（{base}）", expanded=False):
                    st.pyplot(entry["result"].figure, use_container_width=True)
                    cs = st.columns(4)
                    formats = [("PNG", "png", "image/png"),
                               ("PDF", "pdf", "application/pdf"),
                               ("SVG", "svg", "image/svg+xml"),
                               ("TIFF", "tiff", "image/tiff")]
                    for ci, (label, ext, mime) in enumerate(formats):
                        with cs[ci]:
                            try:
                                st.download_button(
                                    f"⬇️ {label}",
                                    data=export_to_bytes(entry["result"].figure, ext),
                                    file_name=f"{base}.{ext}",
                                    mime=mime,
                                    key=f"ec_{aid}_{ext}",
                                )
                            except Exception as e:  # noqa: BLE001
                                st.caption(f"{label} 导出失败：{e}")
                    _stats = getattr(entry["result"], "stats", None)
                    if isinstance(_stats, pd.DataFrame):
                        st.download_button(
                            "⬇️ 统计表（TSV）",
                            data=_stats.to_csv(sep="\t", index=False).encode("utf-8"),
                            file_name=f"{base}_stats.tsv",
                            mime="text/tab-separated-values",
                            key=f"ec_{aid}_tsv",
                        )
                    # 循环图特殊：加交互 HTML 导出按钮
                    if aid == "cycle_diagram":
                        _cd = getattr(entry["result"], "data", None)
                        if _cd is not None:
                            try:
                                from envmeta.geocycle.html_exporter import (
                                    build_interactive_html as _bh,
                                )
                                _html = _bh(
                                    _cd,
                                    hypothesis=st.session_state.get("_hyp_last"),
                                    compare_df=st.session_state.get("_cy_compare_last"),
                                    hypothesis_by_group=st.session_state.get("_hyp_multi_last"),
                                    per_group_cycles=st.session_state.get("_cy_per_group_last"),
                                    full_sample_cycle_data=st.session_state.get("_cy_full_sample_last"),
                                )
                                st.download_button(
                                    f"📦 导出交互 HTML（~{len(_html) // 1024} KB）",
                                    data=_html,
                                    file_name=f"{base}_interactive.html",
                                    mime="text/html",
                                    key=f"ec_{aid}_html",
                                    help="独立 HTML · D3.js 嵌入 · 离线可用",
                                )
                            except Exception as _e:  # noqa: BLE001
                                st.caption(f"⚠️ HTML 导出失败：{_e}")

    # ── tab 2: Fork Bundle ──────────────────────────────────
    with tab_bundle:
        st.markdown(
            "**📦 Fork Bundle — 论文复现包**\n\n"
            "把当前 KB + 假说 YAML + 分析参数打包成单个 .zip，"
            "读者/审稿人一键加载即可复现本论文。"
        )
        st.caption(
            "Bundle 的创建 / 加载功能完整版位于「生物地球化学循环图」页底部 "
            "的 📦 Fork Bundle expander（因为需要循环图 + 假说评分上下文）。"
        )
        st.button(
            "🚀 去「生物地球化学循环图」页操作 Bundle",
            key="ec_goto_cycle",
            on_click=_goto_page_callback, args=("生物地球化学循环图",),
        )
        # 若已加载过 bundle，这里也展示一下
        _bundle_last = st.session_state.get("_bundle_last")
        if _bundle_last:
            st.success(
                f"当前已加载 Bundle：**{_bundle_last['info']['manifest'].get('name', '?')}**"
            )
            st.caption(
                f"hypotheses: {len(_bundle_last['hypothesis_texts'])} · "
                f"config keys: {len(_bundle_last.get('config', {}))} · "
                f"has readme: {_bundle_last.get('readme') is not None}"
            )

    # ── tab 3: 复现脚本 ──────────────────────────────────────
    with tab_code:
        if not registry:
            st.info("还没有已生成的分析，无法生成复现脚本。")
        else:
            st.markdown(
                f"**共 {len(registry)} 个分析可生成独立可运行的 Python 脚本**"
            )
            # 批量 zip 下载所有 .py 脚本
            def _build_scripts_zip() -> bytes:
                buf = _io.BytesIO()
                with _zipfile.ZipFile(buf, "w", _zipfile.ZIP_DEFLATED) as zf:
                    for _aid, _entry in registry.items():
                        try:
                            _src = generate_code(
                                _aid, _entry["file_paths"], _entry["params"],
                                output_base=_entry["output_base"],
                            )
                            zf.writestr(f"{_aid}_reproduce.py", _src)
                        except Exception as e:  # noqa: BLE001
                            zf.writestr(
                                f"{_aid}_ERROR.txt",
                                f"无法生成 {_aid} 的复现脚本：{e}",
                            )
                return buf.getvalue()

            st.download_button(
                "📦 一键批量下载（ZIP — 全部 .py 复现脚本）",
                data=_build_scripts_zip(),
                file_name="envmeta_scripts_bundle.zip",
                mime="application/zip",
                key="export_scripts_zip",
                type="primary",
            )
            st.markdown("---")
            for aid, entry in registry.items():
                name = ANALYSIS_INPUTS.get(aid, {}).get("name", aid)
                with st.expander(f"🐍 {name}", expanded=False):
                    try:
                        src = generate_code(
                            aid, entry["file_paths"], entry["params"],
                            output_base=entry["output_base"],
                        )
                        st.download_button(
                            "⬇️ 下载 .py 脚本",
                            data=src.encode("utf-8"),
                            file_name=f"{aid}_reproduce.py",
                            mime="text/x-python",
                            key=f"ec_code_{aid}",
                        )
                        st.code(src[:1500] + ("\n# ... (truncated)" if len(src) > 1500 else ""),
                                language="python")
                    except Exception as e:  # noqa: BLE001
                        st.caption(f"⚠️ 无法生成复现脚本：{e}")

    # ── tab 4: 文档下载 ──────────────────────────────────────
    with tab_doc:
        st.markdown("**📚 项目文档一键下载**")
        docs = [
            ("docs/data_preparation_zh.md",
             "📥 数据准备指南（中文）",
             "11 种上游工具 → EnvMeta 映射 + 格式样板 + FAQ"),
            ("paper/tool_comparison.md",
             "🏆 工具对比表",
             "EnvMeta vs Shiny-phyloseq / Anvi'o / plotmicrobiome / MicrobiomeAnalyst（18 维度）"),
            ("paper/benchmarks/time_comparison.md",
             "⏱️ 操作效率对比",
             "14 图 × 传统脚本行数/步骤/耗时 vs EnvMeta"),
            ("README.md",
             "📖 README",
             "项目说明 + 快速开始"),
            ("CLAUDE.md",
             "📜 开发日志 + 设计决策",
             "Phase 0 → v0.6 完整演化记录（含产品定位与核心决策）"),
        ]
        for relpath, title, desc in docs:
            fpath = Path(__file__).parent / relpath
            with st.container(border=True):
                c1, c2 = st.columns([3, 1])
                with c1:
                    st.markdown(f"**{title}**")
                    st.caption(desc)
                    st.caption(f"路径：`{relpath}`")
                with c2:
                    if fpath.exists():
                        st.download_button(
                            "⬇️ 下载",
                            data=fpath.read_bytes(),
                            file_name=fpath.name,
                            mime="text/markdown",
                            key=f"ec_doc_{relpath.replace('/', '_')}",
                        )
                    else:
                        st.caption("❌ 未找到")

# ── 页脚 ──────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "EnvMeta · [GitHub](https://github.com/redlizzxy/EnvMeta) · 环境微生物宏基因组可视化分析平台"
)
