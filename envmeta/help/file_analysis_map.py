"""文件类型 ↔ 分析模块双向索引。

ANALYSIS_INPUTS: 每个分析需要哪些文件（required / optional）
FILE_TO_ANALYSIS: 每种文件可以跑哪些分析（反向索引）

数据源为 S8-ux 探索阶段对 envmeta/analysis/*.py 和 app.py 的梳理结果。
新增分析时只需在 ANALYSIS_INPUTS 加一条；FILE_TO_ANALYSIS 自动派生。
"""
from __future__ import annotations

from envmeta.file_manager.detector import FileType


# 分析 ID → 展示信息 + 输入清单
# page: 侧边栏对应哪一页；analysis_type: 页面内部的 selectbox 值（若无则 None）
ANALYSIS_INPUTS: dict[str, dict] = {
    # ─── Reads-based ────────────────────────────────────────────
    "stackplot": {
        "name": "物种组成堆叠图",
        "page": "Reads-based 分析",
        "analysis_type": "物种组成堆叠图",
        "required": [FileType.ABUNDANCE_WIDE, FileType.METADATA],
        "optional": [],
    },
    "pcoa": {
        "name": "β多样性 PCoA",
        "page": "Reads-based 分析",
        "analysis_type": "β多样性 PCoA",
        "required": [FileType.DISTANCE_MATRIX, FileType.METADATA],
        "optional": [],
    },
    "gene_heatmap": {
        "name": "功能基因热图",
        "page": "Reads-based 分析",
        "analysis_type": "功能基因热图",
        "required": [FileType.KO_ABUNDANCE_WIDE, FileType.METADATA],
        "optional": [],
    },
    "alpha_boxplot": {
        "name": "α多样性箱线图",
        "page": "Reads-based 分析",
        "analysis_type": "α多样性",
        "required": [FileType.ALPHA_DIVERSITY, FileType.METADATA],
        "optional": [],
    },
    "log2fc": {
        "name": "基因差异分析 (log2FC)",
        "page": "Reads-based 分析",
        "analysis_type": "基因差异分析 (log2FC)",
        "required": [FileType.KO_ABUNDANCE_WIDE, FileType.METADATA],
        "optional": [],
    },
    "rda": {
        "name": "RDA/CCA 排序",
        "page": "Reads-based 分析",
        "analysis_type": "RDA/CCA 排序",
        "required": [FileType.ABUNDANCE_WIDE, FileType.ENV_FACTORS, FileType.METADATA],
        "optional": [],
    },
    "lefse": {
        "name": "LEfSe 差异分析",
        "page": "Reads-based 分析",
        "analysis_type": "LEfSe 差异分析",
        "required": [FileType.ABUNDANCE_WIDE, FileType.METADATA],
        "optional": [],
    },
    # ─── MAG-based ─────────────────────────────────────────────
    "mag_quality": {
        "name": "MAG 质量评估",
        "page": "MAG-based 分析",
        "analysis_type": "MAG 质量评估",
        "required": [FileType.CHECKM_QUALITY],
        "optional": [FileType.MAG_TAXONOMY, FileType.KEYSTONE_SPECIES],
    },
    "mag_heatmap": {
        "name": "MAG 丰度热图",
        "page": "MAG-based 分析",
        "analysis_type": "MAG 丰度热图",
        "required": [FileType.ABUNDANCE_WIDE],
        "optional": [FileType.MAG_TAXONOMY, FileType.KEYSTONE_SPECIES, FileType.METADATA],
    },
    "pathway": {
        "name": "代谢通路完整度",
        "page": "MAG-based 分析",
        "analysis_type": "代谢通路完整度",
        "required": [FileType.KO_ANNOTATION_LONG],
        "optional": [FileType.ABUNDANCE_WIDE, FileType.MAG_TAXONOMY,
                     FileType.KEYSTONE_SPECIES, FileType.METADATA],
    },
    "gene_profile": {
        "name": "MAG 元素循环基因谱",
        "page": "MAG-based 分析",
        "analysis_type": "MAG 元素循环基因谱",
        "required": [FileType.KO_ANNOTATION_LONG],
        "optional": [FileType.ABUNDANCE_WIDE, FileType.MAG_TAXONOMY,
                     FileType.KEYSTONE_SPECIES],
    },
    "network": {
        "name": "共现网络图（Gephi 辅助）",
        "page": "MAG-based 分析",
        "analysis_type": "共现网络图",
        "required": [FileType.GEPHI_NODES, FileType.GEPHI_EDGES],
        "optional": [FileType.MAG_TAXONOMY, FileType.KEYSTONE_SPECIES],
    },
    # ─── 循环图（单页多功能） ─────────────────────────────────
    "cycle_diagram": {
        "name": "生物地球化学循环图",
        "page": "生物地球化学循环图",
        "analysis_type": None,
        "required": [FileType.KO_ANNOTATION_LONG],
        "optional": [FileType.CHECKM_QUALITY, FileType.MAG_TAXONOMY,
                     FileType.KEYSTONE_SPECIES, FileType.METADATA,
                     FileType.ENV_FACTORS, FileType.ABUNDANCE_WIDE],
    },
    "hypothesis_score": {
        "name": "机制假说评分（YAML）",
        "page": "生物地球化学循环图",
        "analysis_type": None,
        "required": [FileType.KO_ANNOTATION_LONG],
        "optional": [FileType.CHECKM_QUALITY, FileType.MAG_TAXONOMY,
                     FileType.KEYSTONE_SPECIES, FileType.METADATA,
                     FileType.ENV_FACTORS, FileType.ABUNDANCE_WIDE],
    },
}


def _build_file_to_analysis() -> dict[FileType, list[tuple[str, str, str]]]:
    """反向索引：FileType → [(analysis_id, 分析名, 'required'|'optional'), ...]"""
    out: dict[FileType, list[tuple[str, str, str]]] = {ft: [] for ft in FileType}
    for aid, spec in ANALYSIS_INPUTS.items():
        for ft in spec["required"]:
            out[ft].append((aid, spec["name"], "required"))
        for ft in spec["optional"]:
            out[ft].append((aid, spec["name"], "optional"))
    return out


FILE_TO_ANALYSIS: dict[FileType, list[tuple[str, str, str]]] = _build_file_to_analysis()


def analyses_ready(available_types: set[FileType]) -> list[str]:
    """已上传文件类型集合 → 可直接运行的 analysis_id 列表（全部 required 满足）。"""
    ready: list[str] = []
    for aid, spec in ANALYSIS_INPUTS.items():
        if all(ft in available_types for ft in spec["required"]):
            ready.append(aid)
    return ready
