"""
EnvMeta — 环境微生物宏基因组可视化分析平台
Streamlit 主入口：streamlit run app.py
"""

import pandas as pd
import streamlit as st

from envmeta import __version__
from envmeta.analysis import (
    cycle_diagram as cycle_mod,
    gene_heatmap, gene_profile as gene_profile_mod,
    lefse as lefse_mod,
    mag_quality as mag_quality_mod,
    pathway as pathway_mod,
    pcoa, rda as rda_mod, stackplot,
)
from envmeta.export.code_generator import generate as generate_code
from envmeta.export.figure_export import export_to_bytes
from envmeta.file_manager.detector import FileType, detect, read_table
from envmeta.params.common import render_figure_size, render_font_controls


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
            _reproduce_button("stackplot",
                              {"abundance": ab_name, "metadata": md_name},
                              last.params, key="stack_code",
                              output_base=f"stackplot_{last.params['style']}")
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
            _reproduce_button("pcoa",
                              {"distance": dist_name, "metadata": md_name},
                              last.params, key="pcoa_code")
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
            _reproduce_button("gene_heatmap",
                              {"ko_abundance": ko_name, "metadata": md_name},
                              last.params, key="heat_code")
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
            _reproduce_button("alpha_boxplot",
                              {"alpha": alpha_name, "metadata": md_name},
                              last.params, key="alpha_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # ── 基因差异分析 log2FC ─────────────────────────────
    elif analysis_type == "基因差异分析 (log2FC)":
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
            _reproduce_button("log2fc",
                              {"ko_abundance": ko_name, "metadata": md_name},
                              last.params, key="log2fc_code",
                              output_base=f"log2fc_{suffix}")
            with st.expander("查看统计表（按 padj 排序）"):
                st.dataframe(last.stats.sort_values("padj").round(4),
                             use_container_width=True)

    # ── RDA 排序图 ─────────────────────────────────
    elif analysis_type == "RDA/CCA 排序":
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
        with c1:
            ab_name = st.selectbox("丰度表", list(abundance_files.keys()), key="rda_ab")
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
            _reproduce_button("rda",
                              {"abundance": ab_name, "env_factors": env_name, "metadata": md_name},
                              last.params, key="rda_code")
            with st.expander("查看 RDA 统计表"):
                st.dataframe(last.stats.round(4), use_container_width=True)

    # ── LEfSe 差异分析 ──────────────────────────────────
    elif analysis_type == "LEfSe 差异分析":
        abundance_files = _files_of(FileType.ABUNDANCE_WIDE)
        metadata_files = _files_of(FileType.METADATA)
        if not abundance_files or not metadata_files:
            st.warning(
                f"需要 1 个丰度表（Species/Genus/Phylum）+ 1 个 metadata。"
                f"当前：丰度表 {len(abundance_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()
        ab_name = st.selectbox("丰度表", list(abundance_files.keys()), key="lefse_ab")
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
    )

    if analysis_type == "MAG 质量评估":
        checkm_files = _files_of(FileType.CHECKM_QUALITY)
        if not checkm_files:
            st.warning("需要 1 个 CheckM / CheckM2 质量表（含 Completeness/Contamination/Genome_Size 列）。")
            st.stop()
        q_name = st.selectbox("CheckM 质量表", list(checkm_files.keys()), key="mq_checkm")
        t_name = st.selectbox(
            "GTDB 分类表（可选）",
            ["（无）"] + [n for n in st.session_state.files.keys() if n != q_name],
            key="mq_tax",
        )
        k_name = st.selectbox(
            "Keystone 物种列表（可选）",
            ["（无）"] + [n for n in st.session_state.files.keys()
                         if n not in (q_name, t_name)],
            key="mq_ks",
        )

        with st.sidebar:
            st.subheader("MAG 质量参数")
            hc = st.slider("高质量 Completeness ≥", 50, 100, 90, key="mq_hc")
            hcon = st.slider("高质量 Contamination <", 1, 10, 5, key="mq_hcon")
            mc = st.slider("中质量 Completeness ≥", 30, 80, 50, key="mq_mc")
            mcon = st.slider("中质量 Contamination <", 5, 20, 10, key="mq_mcon")
            size = render_figure_size({"width_mm": 260, "height_mm": 140}, prefix="mq")

        if st.button("生成 MAG 质量图", type="primary", key="mq_go"):
            q = checkm_files[q_name]["df"]
            t = st.session_state.files[t_name]["df"] if t_name != "（无）" else None
            # taxonomy 表若没有标题，需要设列名
            if t is not None and "classification" not in [c.lower() for c in t.columns] \
                    and "taxonomy" not in [c.lower() for c in t.columns]:
                # 约定：无表头的 MAG\tTaxonomy 两列格式
                if t.shape[1] == 2:
                    t = t.copy()
                    t.columns = ["MAG", "Taxonomy"]
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            params = {
                "high_completeness": float(hc), "high_contamination": float(hcon),
                "med_completeness": float(mc), "med_contamination": float(mcon),
                **size,
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
            _reproduce_button("mag_quality", files_map, last.params, key="mq_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "代谢通路完整度":
        # KO 注释表的文件类型目前落在 UNKNOWN 或 ABUNDANCE_WIDE，允许用户手选任何表
        ko_candidates = list(st.session_state.files.keys())
        if not ko_candidates:
            st.warning("需要至少 1 个 KO 注释表（长格式：MAG + KEGG_ko 两列）。")
            st.stop()
        ko_name = st.selectbox("KO 注释表（MAG + KEGG_ko 长表）", ko_candidates, key="pw_ko")
        t_name = st.selectbox(
            "GTDB 分类表（可选）",
            ["（无）"] + [n for n in ko_candidates if n != ko_name],
            key="pw_tax",
        )
        k_name = st.selectbox(
            "Keystone 物种列表（可选）",
            ["（无）"] + [n for n in ko_candidates if n not in (ko_name, t_name)],
            key="pw_ks",
        )
        a_name = st.selectbox(
            "MAG 丰度表（可选，用于 Top-N 过滤 / 气泡大小）",
            ["（无）"] + [n for n in ko_candidates
                         if n not in (ko_name, t_name, k_name)],
            key="pw_ab",
        )

        with st.sidebar:
            st.subheader("通路完整度参数")
            style = st.radio("图样式", ["heatmap", "bubble"], key="pw_style")
            sort_by = st.selectbox(
                "排序方式",
                ["phylum_then_total", "total", "abundance"],
                key="pw_sort",
            )
            elem_opts = ["arsenic", "sulfur", "iron", "nitrogen"]
            elem_sel = st.multiselect("元素过滤（空=全部）", elem_opts,
                                      default=[], key="pw_elem")
            max_mags = st.slider("最多显示 MAG 数（0=全部）", 0, 200, 50,
                                 step=10, key="pw_max")
            size = render_figure_size({"width_mm": 280, "height_mm": 200},
                                      prefix="pw")

        if st.button("生成通路完整度图", type="primary", key="pw_go"):
            ko = st.session_state.files[ko_name]["df"]
            t = st.session_state.files[t_name]["df"] if t_name != "（无）" else None
            if t is not None and t.shape[1] == 2 and "classification" not in [
                    c.lower() for c in t.columns] and "taxonomy" not in [
                    c.lower() for c in t.columns]:
                t = t.copy()
                t.columns = ["MAG", "Taxonomy"]
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            a = st.session_state.files[a_name]["df"] if a_name != "（无）" else None
            params = {
                "style": style,
                "sort_by": sort_by,
                "element_filter": elem_sel or None,
                "max_mags": max_mags or None,
                **size,
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
            if t_name != "（无）": files_map["taxonomy"] = t_name
            if k_name != "（无）": files_map["keystone"] = k_name
            if a_name != "（无）": files_map["abundance"] = a_name
            _reproduce_button("pathway", files_map, last.params, key="pw_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    elif analysis_type == "MAG 元素循环基因谱":
        ko_candidates = list(st.session_state.files.keys())
        if not ko_candidates:
            st.warning("需要至少 1 个 KO 注释表（MAG + KEGG_ko 长表）。")
            st.stop()
        ko_name = st.selectbox("KO 注释表", ko_candidates, key="gp_ko")
        t_name = st.selectbox(
            "GTDB 分类表（可选）",
            ["（无）"] + [n for n in ko_candidates if n != ko_name],
            key="gp_tax",
        )
        k_name = st.selectbox(
            "Keystone 物种列表（可选）",
            ["（无）"] + [n for n in ko_candidates if n not in (ko_name, t_name)],
            key="gp_ks",
        )
        a_name = st.selectbox(
            "MAG 丰度表（可选）",
            ["（无）"] + [n for n in ko_candidates
                         if n not in (ko_name, t_name, k_name)],
            key="gp_ab",
        )

        with st.sidebar:
            st.subheader("基因谱参数")
            sort_by = st.selectbox(
                "排序方式",
                ["phylum_then_count", "count", "abundance"],
                key="gp_sort",
            )
            elem_opts = ["arsenic", "sulfur", "iron", "nitrogen"]
            elem_sel = st.multiselect("元素过滤（空=全部）", elem_opts,
                                      default=[], key="gp_elem")
            max_mags = st.slider("最多显示 MAG 数（0=全部）", 0, 200, 60,
                                 step=10, key="gp_max")
            show_names = st.checkbox("列标签显示基因名", True, key="gp_names")
            size = render_figure_size({"width_mm": 340, "height_mm": 220},
                                      prefix="gp")

        if st.button("生成基因谱图", type="primary", key="gp_go"):
            ko = st.session_state.files[ko_name]["df"]
            t = st.session_state.files[t_name]["df"] if t_name != "（无）" else None
            if t is not None and t.shape[1] == 2 and "classification" not in [
                    c.lower() for c in t.columns] and "taxonomy" not in [
                    c.lower() for c in t.columns]:
                t = t.copy()
                t.columns = ["MAG", "Taxonomy"]
            k = st.session_state.files[k_name]["df"] if k_name != "（无）" else None
            a = st.session_state.files[a_name]["df"] if a_name != "（无）" else None
            params = {
                "sort_by": sort_by,
                "element_filter": elem_sel or None,
                "max_mags": max_mags or None,
                "show_gene_names": show_names,
                **size,
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
            if t_name != "（无）": files_map["taxonomy"] = t_name
            if k_name != "（无）": files_map["keystone"] = k_name
            if a_name != "（无）": files_map["abundance"] = a_name
            _reproduce_button("gene_profile", files_map, last.params, key="gp_code")
            with st.expander("查看统计表"):
                st.dataframe(last.stats, use_container_width=True)
    else:
        st.info(f"「{analysis_type}」模块开发中 — Phase 2 继续。")

elif page == "生物地球化学循环图":
    st.title("生物地球化学循环图生成器")
    st.caption("Phase 3 v1 — 从 MAG × KO + 环境因子自动推断活跃元素循环通路")

    file_names = list(st.session_state.files.keys())
    if not file_names:
        st.warning("请先在「文件管理」上传 KO 注释表（MAG + KEGG_ko 长表）。")
        st.stop()

    ko_name = st.selectbox("KO 注释表（必需）", file_names, key="cy_ko")
    t_name = st.selectbox(
        "GTDB 分类表（可选，用于 Phylum/Genus 标签）",
        ["（无）"] + [n for n in file_names if n != ko_name],
        key="cy_tax",
    )
    k_name = st.selectbox(
        "Keystone 物种列表（可选）",
        ["（无）"] + [n for n in file_names if n not in (ko_name, t_name)],
        key="cy_ks",
    )
    a_name = st.selectbox(
        "MAG 丰度表（可选，用于通路贡献加权）",
        ["（无）"] + [n for n in file_names
                     if n not in (ko_name, t_name, k_name)],
        key="cy_ab",
    )
    e_name = st.selectbox(
        "环境因子表（可选，用于 env-pathway 相关性）",
        ["（无）"] + [n for n in file_names
                     if n not in (ko_name, t_name, k_name, a_name)],
        key="cy_env",
    )
    m_name = st.selectbox(
        "Metadata（env 表需要搭配）",
        ["（无）"] + [n for n in file_names
                     if n not in (ko_name, t_name, k_name, a_name, e_name)],
        key="cy_md",
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
        show_env = st.checkbox("显示环境耦合面板", True, key="cy_env_panel")
        size = render_figure_size({"width_mm": 360, "height_mm": 260},
                                  prefix="cy")

    if st.button("生成循环图", type="primary", key="cy_go"):
        def _get(n):
            return st.session_state.files[n]["df"] if n != "（无）" else None
        ko = st.session_state.files[ko_name]["df"]
        t = _get(t_name)
        if t is not None and t.shape[1] == 2 and "classification" not in [
                c.lower() for c in t.columns] and "taxonomy" not in [
                c.lower() for c in t.columns]:
            t = t.copy(); t.columns = ["MAG", "Taxonomy"]
        params = {
            "completeness_threshold": float(comp_thresh),
            "top_n_contributors": top_n,
            "env_rho_min": rho_min,
            "env_p_max": p_max,
            "show_env_panel": show_env,
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

elif page == "导出中心":
    st.title("导出中心")
    st.info("批量导出功能将在 Phase 2 实现。单图导出在各分析页面下方。")

# ── 页脚 ──────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "EnvMeta · [GitHub](https://github.com/redlizzxy/EnvMeta) · 环境微生物宏基因组可视化分析平台"
)
