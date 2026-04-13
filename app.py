"""
EnvMeta — 环境微生物宏基因组可视化分析平台
Streamlit 主入口：streamlit run app.py
"""

import pandas as pd
import streamlit as st

from envmeta import __version__
from envmeta.analysis import stackplot
from envmeta.export.figure_export import export_to_bytes
from envmeta.file_manager.detector import FileType, detect, read_table

# ── 页面配置 ──────────────────────────────────────────────
st.set_page_config(
    page_title="EnvMeta",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# session state 初始化：存 {filename: {"df": DataFrame, "type": FileType, "result": DetectionResult}}
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
    FileType.METADATA: ("🗂️ metadata", "blue"),
    FileType.ABUNDANCE_WIDE: ("📊 abundance (wide)", "green"),
    FileType.UNKNOWN: ("❓ unknown", "gray"),
}
TYPE_OPTIONS = [ft.value for ft in FileType]


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
        1. 在 **文件管理** 中上传你的数据文件（丰度表、KO 注释表、metadata 等）
        2. 选择分析模块，系统自动匹配可用的输入文件
        3. 调整图形参数，实时预览
        4. 在 **导出中心** 一键导出高质量图形和可复现代码

        #### 当前可用功能（Phase 1 迭代 1）
        - 文件识别：metadata、abundance（宽表）
        - 分析：物种组成堆叠图（sample / group 两种样式）
        - 导出：PNG / PDF（300 DPI）
        """
    )

# ══════════════════════════════════════════════════════════
# 文件管理（模块 A）
# ══════════════════════════════════════════════════════════
elif page == "文件管理":
    st.title("文件管理")

    uploaded = st.file_uploader(
        "上传数据文件（可多选）",
        accept_multiple_files=True,
        type=None,  # 不限制扩展名
        help="支持 .txt / .tsv / .csv。系统会自动识别文件类型。",
    )

    # 处理新上传：覆盖同名文件
    for up in uploaded or []:
        if up.name not in st.session_state.files:
            try:
                df, enc, sep = read_table(up)
                up.seek(0)
                result = detect(up, filename=up.name)
                st.session_state.files[up.name] = {
                    "df": df,
                    "type": result.file_type,
                    "result": result,
                }
            except Exception as e:
                st.error(f"读取 {up.name} 失败：{e}")

    if not st.session_state.files:
        st.info("请上传文件。也可以在 PowerShell 里拖 `tests/sample_data/*` 里的文件做测试。")
    else:
        st.subheader(f"已识别文件（{len(st.session_state.files)} 个）")
        for fname, info in list(st.session_state.files.items()):
            with st.expander(f"📄 {fname} — {TYPE_BADGES[info['type']][0]}", expanded=True):
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
                        "手动修正类型",
                        TYPE_OPTIONS,
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
            "α多样性",
            "β多样性 PCoA",
            "RDA/CCA 排序",
            "LEfSe 差异分析",
            "功能基因热图",
            "基因差异分析 (log2FC)",
        ],
    )

    if analysis_type != "物种组成堆叠图":
        st.info(f"「{analysis_type}」将在 Phase 1 迭代 2 或 Phase 2 实现。")
    else:
        # ── 堆叠图页面 ──────────────────────────────────────
        abundance_files = {
            n: info for n, info in st.session_state.files.items()
            if info["type"] == FileType.ABUNDANCE_WIDE
        }
        metadata_files = {
            n: info for n, info in st.session_state.files.items()
            if info["type"] == FileType.METADATA
        }

        if not abundance_files or not metadata_files:
            st.warning(
                "需要先在 **文件管理** 上传并识别 1 个丰度表（abundance_wide）+ 1 个 metadata。"
                f"当前：丰度表 {len(abundance_files)}，metadata {len(metadata_files)}。"
            )
            st.stop()

        col_sel1, col_sel2 = st.columns(2)
        with col_sel1:
            ab_name = st.selectbox("丰度表", list(abundance_files.keys()))
        with col_sel2:
            md_name = st.selectbox("Metadata", list(metadata_files.keys()))

        st.sidebar.markdown("---")
        st.sidebar.subheader("参数")
        style = st.sidebar.radio("横轴", ["sample", "group"], horizontal=True)
        top_n = st.sidebar.slider("Top-N", 5, 20, 10)
        drop_unc = st.sidebar.checkbox("过滤 unclassified 行", value=True)
        width_mm = st.sidebar.slider("图宽 (mm)", 80, 300, 160, step=10)
        height_mm = st.sidebar.slider("图高 (mm)", 60, 200, 100, step=10)

        params = {
            "style": style, "top_n": top_n, "drop_unclassified": drop_unc,
            "width_mm": width_mm, "height_mm": height_mm,
        }

        if st.button("🎨 生成图表", type="primary"):
            try:
                result = stackplot.analyze(
                    abundance_files[ab_name]["df"],
                    metadata_files[md_name]["df"],
                    params,
                )
                st.session_state["last_stackplot"] = result
            except Exception as e:
                st.error(f"生成失败：{e}")

        last = st.session_state.get("last_stackplot")
        if last is not None:
            st.pyplot(last.figure, use_container_width=True)

            col_d1, col_d2, col_d3 = st.columns(3)
            with col_d1:
                st.download_button(
                    "⬇️ PNG（300 DPI）",
                    data=export_to_bytes(last.figure, "png"),
                    file_name=f"stackplot_{last.params['style']}.png",
                    mime="image/png",
                )
            with col_d2:
                st.download_button(
                    "⬇️ PDF（矢量）",
                    data=export_to_bytes(last.figure, "pdf"),
                    file_name=f"stackplot_{last.params['style']}.pdf",
                    mime="application/pdf",
                )
            with col_d3:
                st.download_button(
                    "⬇️ 百分比表（TSV）",
                    data=last.stats.to_csv(sep="\t").encode("utf-8"),
                    file_name=f"stackplot_{last.params['style']}_percentage.tsv",
                    mime="text/tab-separated-values",
                )

            with st.expander("查看百分比表"):
                st.dataframe(last.stats.round(2), use_container_width=True)

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
    st.info("批量导出功能将在 Phase 2 实现。堆叠图的单图导出在「Reads-based 分析」页面下方。")

# ── 页脚 ──────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "EnvMeta · [GitHub](https://github.com/redlizzxy/EnvMeta) · 环境微生物宏基因组可视化分析平台"
)
