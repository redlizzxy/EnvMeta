"""
EnvMeta — 环境微生物宏基因组可视化分析平台
Streamlit 主入口：streamlit run app.py
"""

import streamlit as st

from envmeta import __version__

# ── 页面配置 ──────────────────────────────────────────────
st.set_page_config(
    page_title="EnvMeta",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── 侧边栏导航 ────────────────────────────────────────────
st.sidebar.title("EnvMeta")
st.sidebar.caption(f"v{__version__} · 环境微生物宏基因组可视化分析平台")

page = st.sidebar.radio(
    "功能模块",
    [
        "首页",
        "文件管理",          # 模块 A
        "Reads-based 分析",  # 模块 B1
        "MAG-based 分析",    # 模块 B2
        "生物地球化学循环图",  # 模块 E
        "导出中心",          # 模块 D
    ],
)

# ── 首页 ──────────────────────────────────────────────────
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

        #### 功能模块
        | 模块 | 说明 |
        |------|------|
        | 文件管理 | 智能识别上传文件类型，自动推荐可用分析 |
        | Reads-based 分析 | 堆叠图 · α多样性 · PCoA · RDA · LEfSe · 热图 · log2FC |
        | MAG-based 分析 | MAG质量 · 丰度热图 · 通路完整度 · 基因谱 · 网络图 |
        | 生物地球化学循环图 | 从 KO 注释自动推断活跃通路，生成可编辑的循环机制图 |
        | 导出中心 | PDF/SVG/TIFF/PNG 图形 + 可运行 Python 脚本 |
        """
    )

# ── 文件管理（模块 A）────────────────────────────────────
elif page == "文件管理":
    st.title("文件管理")
    st.info("模块开发中 — Phase 1 将实现文件上传与智能识别功能。")

# ── Reads-based 分析（模块 B1）────────────────────────────
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
    st.info(f"「{analysis_type}」模块开发中 — Phase 1/2 将逐步实现。")

# ── MAG-based 分析（模块 B2）─────────────────────────────
elif page == "MAG-based 分析":
    st.title("MAG-based 基因组分析")
    analysis_type = st.selectbox(
        "选择分析类型",
        [
            "MAG 质量评估",
            "MAG 丰度热图",
            "代谢通路完整度",
            "MAG 元素循环基因谱",
            "共现网络图",
        ],
    )
    st.info(f"「{analysis_type}」模块开发中 — Phase 2 将实现。")

# ── 生物地球化学循环图（模块 E）──────────────────────────
elif page == "生物地球化学循环图":
    st.title("生物地球化学循环图生成器")
    st.info("核心创新模块 — Phase 3/4 将实现通路推断与交互编辑。")

# ── 导出中心（模块 D）────────────────────────────────────
elif page == "导出中心":
    st.title("导出中心")
    st.info("模块开发中 — Phase 1 将实现基础图形导出功能。")

# ── 页脚 ──────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "EnvMeta · [GitHub](https://github.com/) · 环境微生物宏基因组可视化分析平台"
)
