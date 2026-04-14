"""把 analysis 的 params 渲染为独立可运行的 Python 脚本。

公开接口：
    generate(analysis_id, file_paths, params, output_base="out") -> str

生成的脚本读取指定文件，调用对应 `envmeta.analysis.XXX.analyze`，保存
`{output_base}.pdf`、`{output_base}_stats.tsv` 两个产物。
"""
from __future__ import annotations

import datetime
import pprint
from typing import Callable

from envmeta import __version__

SUPPORTED = {"stackplot", "pcoa", "gene_heatmap", "alpha_boxplot", "log2fc", "rda", "lefse", "mag_quality", "pathway", "gene_profile"}

_HEADER_TEMPLATE = '''"""
由 EnvMeta v{version} 于 {timestamp} 生成。

复现步骤：
1. 激活环境：conda activate envmeta
2. 确保输入文件路径正确（见下方 FILES 字典）
3. 运行：python {script_name}.py

产物：
- {output_base}.pdf
- {output_base}_stats.tsv
"""
from __future__ import annotations

import pandas as pd

from envmeta.analysis import {module}
from envmeta.export.figure_export import export_figure

# ── 输入文件 ──────────────────────────────────────────────
FILES = {files_repr}

# ── 分析参数（EnvMeta GUI 导出） ──────────────────────────
PARAMS = {params_repr}

# ── 执行 ──────────────────────────────────────────────────
'''


def _clean_params(params: dict) -> dict:
    """剔除以 `_` 开头的私有键（比如中间 DataFrame、纯缓存）。"""
    return {k: v for k, v in params.items() if not k.startswith("_")}


def _format_dict(d: dict) -> str:
    """用 pprint 保留顺序 + 换行美化。"""
    return pprint.pformat(d, sort_dicts=False, width=90)


def _script_name(analysis_id: str) -> str:
    return f"{analysis_id}_reproduce"


def _header(analysis_id: str, module: str, files: dict, params: dict, output_base: str) -> str:
    return _HEADER_TEMPLATE.format(
        version=__version__,
        timestamp=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        script_name=_script_name(analysis_id),
        output_base=output_base,
        module=module,
        files_repr=_format_dict(files),
        params_repr=_format_dict(params),
    )


# ============================================================================
# 各 analysis 模板
# ============================================================================

def _tpl_stackplot(files: dict, params: dict, output_base: str) -> str:
    return _header("stackplot", "stackplot", files, params, output_base) + f'''
abundance_df = pd.read_csv(FILES["abundance"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = stackplot.analyze(abundance_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t")
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_pcoa(files: dict, params: dict, output_base: str) -> str:
    return _header("pcoa", "pcoa", files, params, output_base) + f'''
distance_df = pd.read_csv(FILES["distance"], sep="\\t", index_col=0)
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = pcoa.analyze(distance_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_gene_heatmap(files: dict, params: dict, output_base: str) -> str:
    return _header("gene_heatmap", "gene_heatmap", files, params, output_base) + f'''
ko_df = pd.read_csv(FILES["ko_abundance"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = gene_heatmap.analyze(ko_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t")
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_alpha_boxplot(files: dict, params: dict, output_base: str) -> str:
    return _header("alpha_boxplot", "alpha_boxplot", files, params, output_base) + f'''
alpha_df = pd.read_csv(FILES["alpha"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = alpha_boxplot.analyze(alpha_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_log2fc(files: dict, params: dict, output_base: str) -> str:
    return _header("log2fc", "log2fc", files, params, output_base) + f'''
ko_df = pd.read_csv(FILES["ko_abundance"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = log2fc.analyze(ko_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_rda(files: dict, params: dict, output_base: str) -> str:
    return _header("rda", "rda", files, params, output_base) + f'''
abundance_df = pd.read_csv(FILES["abundance"], sep="\\t")
env_df = pd.read_csv(FILES["env_factors"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = rda.analyze(abundance_df, env_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_lefse(files: dict, params: dict, output_base: str) -> str:
    return _header("lefse", "lefse", files, params, output_base) + f'''
abundance_df = pd.read_csv(FILES["abundance"], sep="\\t")
metadata_df = pd.read_csv(FILES["metadata"], sep="\\t")

result = lefse.analyze(abundance_df, metadata_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_mag_quality(files: dict, params: dict, output_base: str) -> str:
    return _header("mag_quality", "mag_quality", files, params, output_base) + f'''
quality_df = pd.read_csv(FILES["quality"], sep="\\t")
taxonomy_df = pd.read_csv(FILES["taxonomy"], sep="\\t", header=None, names=["MAG", "Taxonomy"]) if FILES.get("taxonomy") else None
keystone_df = pd.read_csv(FILES["keystone"], sep="\\t") if FILES.get("keystone") else None

result = mag_quality.analyze(quality_df, taxonomy_df, keystone_df, PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_gene_profile(files: dict, params: dict, output_base: str) -> str:
    return _header("gene_profile", "gene_profile", files, params, output_base) + f'''
ko_annotation_df = pd.read_csv(FILES["ko_annotation"], sep="\\t")
taxonomy_df = pd.read_csv(FILES["taxonomy"], sep="\\t", header=None, names=["MAG", "Taxonomy"]) if FILES.get("taxonomy") else None
keystone_df = pd.read_csv(FILES["keystone"], sep="\\t") if FILES.get("keystone") else None
abundance_df = pd.read_csv(FILES["abundance"], sep="\\t") if FILES.get("abundance") else None

result = gene_profile.analyze(ko_annotation_df, taxonomy_df, keystone_df, abundance_df, params=PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


def _tpl_pathway(files: dict, params: dict, output_base: str) -> str:
    return _header("pathway", "pathway", files, params, output_base) + f'''
ko_annotation_df = pd.read_csv(FILES["ko_annotation"], sep="\\t")
taxonomy_df = pd.read_csv(FILES["taxonomy"], sep="\\t", header=None, names=["MAG", "Taxonomy"]) if FILES.get("taxonomy") else None
keystone_df = pd.read_csv(FILES["keystone"], sep="\\t") if FILES.get("keystone") else None
abundance_df = pd.read_csv(FILES["abundance"], sep="\\t") if FILES.get("abundance") else None

result = pathway.analyze(ko_annotation_df, taxonomy_df, keystone_df, abundance_df, params=PARAMS)
export_figure(result.figure, "{output_base}.pdf", "pdf")
result.stats.to_csv("{output_base}_stats.tsv", sep="\\t", index=False)
print(f"Saved: {output_base}.pdf + {output_base}_stats.tsv")
'''


_TEMPLATES: dict[str, Callable[[dict, dict, str], str]] = {
    "stackplot": _tpl_stackplot,
    "pcoa": _tpl_pcoa,
    "gene_heatmap": _tpl_gene_heatmap,
    "alpha_boxplot": _tpl_alpha_boxplot,
    "log2fc": _tpl_log2fc,
    "rda": _tpl_rda,
    "lefse": _tpl_lefse,
    "mag_quality": _tpl_mag_quality,
    "pathway": _tpl_pathway,
    "gene_profile": _tpl_gene_profile,
}


def generate(analysis_id: str, file_paths: dict[str, str],
             params: dict, output_base: str = "out") -> str:
    """为指定分析生成独立的 Python 脚本。

    参数
    -----
    analysis_id: 分析模块名（stackplot / pcoa / gene_heatmap / alpha_boxplot / log2fc / rda）
    file_paths : {"abundance": "data/Phylum.txt", "metadata": "data/metadata.txt"}
    params     : GUI 导出的参数字典（会剔除 `_` 开头的私有键）
    output_base: 产物文件名前缀
    """
    if analysis_id not in SUPPORTED:
        raise ValueError(
            f"不支持的 analysis_id={analysis_id!r}，支持：{sorted(SUPPORTED)}")
    clean = _clean_params(params)
    return _TEMPLATES[analysis_id](file_paths, clean, output_base)
