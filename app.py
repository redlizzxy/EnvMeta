"""
EnvMeta — 环境微生物宏基因组可视化分析平台
Streamlit 主入口：streamlit run app.py
"""

import pandas as pd
import streamlit as st

from envmeta import __version__
from envmeta.analysis import gene_heatmap, pcoa, stackplot
from envmeta.export.figure_export import export_to_bytes
from envmeta.file_manager.detector import FileType, detect, read_table
from envmeta.params.common import render_figure_size, render_font_controls

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

page = st.sidebar.radio(
    "功能模块",
    [
        "首页",
        "文件管理",
        "Reads-based 分析",
        "MAG-based 分析",
        "生物地球化学循环图",
        "导出中心",
    ],
)

TYPE_BADGES = {
    FileType.METADATA: "🗂️ metadata",
    FileType.ABUNDANCE_WIDE: "📊 abundance (wide)",
    FileType.DISTANCE_MATRIX: "📐 distance matrix",
    FileType.ALPHA_DIVERSITY: "📈 alpha diversity",
    FileType.CHECKM_QUALITY: "🧬 CheckM quality",
    FileType.ENV_FACTORS: "🌡️ env factors",
    FileType.KO_ABUNDANCE_WIDE: "🧪 KO abundance",
    FileType.UNKNOWN: "❓ unknown",
}
TYPE_OPTIONS = [ft.value for ft in FileType]


def _files_of(*types: FileType) -> dict[str, dict]:
    return {n: info for n, info in st.session_state.files.items() if info["type"] in types}


# ══════════════════════════════════════════════════════════
# 首页
# ══════════════════════════════════════════════════════════
if page == "首页":
    st.title("EnvMeta")
    st.subheader("环境微生物宏基因组可视化分析平台")
    st.markdown("---")
    st.markdown(
        """
        **EnvMeta** 帮助你从宏基因组下游分析的数据文件快速生成发表级图形。

        #### 快速开始
        1. 在 **文件管理** 中上传你的数据文件
        2. 选择分析模块，系统自动匹配可用的输入文件
        3. 调整图形参数，实时预览
        4. 一键导出高质量图形（PNG / PDF）和可复现数据

        #### 当前可用功能（Phase 1 迭代 2）
        | 模块 | 支持内容 |
        |------|---------|
        | 文件识别 | metadata、abundance（宽表）、distance matrix、alpha diversity、CheckM、env factors、KO abundance |
        | Reads-based 分析 | 物种组成堆叠图、PCoA + PERMANOVA、元素循环基因热图 |
        | 调参 | 画布尺寸、字号、排序方式、元素过滤等 |
        | 导出 | PNG（300 DPI）/ PDF（矢量）/ 统计 TSV |
        """
    )

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
    )

    # ── 堆叠图 ───────────────────────────────────────────
    if analysis_type == "物种组成堆叠图":
        abundance_files = _files_of(FileType.ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not abundance_files or not metadata_files:
            st.warning(
                f"需要 1 个丰度表 + 1 个 metadata。当前：丰度表 {len(abundance_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        c1, c2 = st.columns(2)
        with c1:
            ab_name = st.selectbox("丰度表", list(abundance_files.keys()))
        with c2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()), key="stack_md")

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        style = st.sidebar.radio("横轴", ["sample", "group"], horizontal=True, key="stack_style")
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
            with st.expander("查看百分比表"):
                st.dataframe(last.stats.round(2), use_container_width=True)

    # ── PCoA ────────────────────────────────────────────
    elif analysis_type == "β多样性 PCoA":
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
            with st.expander("PERMANOVA 两两对比"):
                st.dataframe(last.params["_stats_tables"]["pairwise_permanova"],
                             use_container_width=True)

    # ── 元素循环基因热图 ───────────────────────────────
    elif analysis_type == "功能基因热图":
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
            with st.expander("原始分组均值矩阵"):
                st.dataframe(last.params["_raw_means"].round(2),
                             use_container_width=True)

    # ── α 多样性箱线图 ─────────────────────────────────
    elif analysis_type == "α多样性":
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
            with st.expander("查看统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # 未实现的子页面
    else:
        st.info(f"「{analysis_type}」将在 Phase 1 迭代 3/4 或 Phase 2 实现。")

# ══════════════════════════════════════════════════════════
# 其他占位页面
# ══════════════════════════════════════════════════════════
elif page == "MAG-based 分析":
    st.title("MAG-based 基因组分析")
    analysis_type = st.selectbox(
        "选择分析类型",
        ["MAG 质量评估", "MAG 丰度热图", "代谢通路完整度", "MAG 元素循环基因谱", "共现网络图"],
    )
    st.info(f"「{analysis_type}」模块开发中 — Phase 2 将实现。")

elif page == "生物地球化学循环图":
    st.title("生物地球化学循环图生成器")
    st.info("核心创新模块 — Phase 3/4 将实现通路推断与交互编辑。")

elif page == "导出中心":
    st.title("导出中心")
    st.info("批量导出功能将在 Phase 2 实现。单图导出在各分析页面下方。")

# ── 页脚 ──────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "EnvMeta · [GitHub](https://github.com/redlizzxy/EnvMeta) · 环境微生物宏基因组可视化分析平台"
)
