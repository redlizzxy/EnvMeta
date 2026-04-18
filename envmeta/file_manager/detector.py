"""文件类型识别引擎。

基于表头规则匹配 + 数据特征检测。通过向 _RULES 追加新规则即可扩展支持的类型。
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, Union

import chardet
import numpy as np
import pandas as pd


class FileType(str, Enum):
    METADATA = "metadata"
    ABUNDANCE_WIDE = "abundance_wide"
    DISTANCE_MATRIX = "distance_matrix"
    ALPHA_DIVERSITY = "alpha_diversity"
    CHECKM_QUALITY = "checkm_quality"
    ENV_FACTORS = "env_factors"
    KO_ABUNDANCE_WIDE = "ko_abundance_wide"
    KO_ANNOTATION_LONG = "ko_annotation_long"
    KEYSTONE_SPECIES = "keystone_species"
    MAG_TAXONOMY = "mag_taxonomy"
    GEPHI_NODES = "gephi_nodes"
    GEPHI_EDGES = "gephi_edges"
    UNKNOWN = "unknown"


@dataclass
class DetectionResult:
    file_type: FileType
    confidence: float                    # 0.0 ~ 1.0
    reasons: list[str] = field(default_factory=list)
    preview_df: pd.DataFrame | None = None
    encoding: str = "utf-8"
    separator: str = "\t"


# ============================================================================
# 底层 IO：编码 + 分隔符嗅探
# ============================================================================

def _sniff_encoding(raw: bytes) -> str:
    """用前 32KB 猜编码，置信度低于 0.5 时回退 utf-8。"""
    sample = raw[:32_000]
    guess = chardet.detect(sample)
    enc = (guess.get("encoding") or "utf-8").lower()
    if guess.get("confidence", 0) < 0.5:
        return "utf-8"
    return {"ascii": "utf-8", "gb2312": "gbk"}.get(enc, enc)


def _sniff_separator(text_sample: str) -> str:
    first_line = text_sample.split("\n", 1)[0]
    return "\t" if first_line.count("\t") >= first_line.count(",") else ","


def read_table(source: Union[str, Path, io.BytesIO], n_preview: int = 5) -> tuple[pd.DataFrame, str, str]:
    """从路径或上传的 BytesIO 读取表格，返回 (df, encoding, separator)。"""
    if hasattr(source, "read"):
        raw = source.read()
        source.seek(0)
    else:
        raw = Path(source).read_bytes()

    enc = _sniff_encoding(raw)
    text = raw.decode(enc, errors="replace")
    sep = _sniff_separator(text)
    df = pd.read_csv(io.StringIO(text), sep=sep, dtype=str, keep_default_na=False)
    return df, enc, sep


# ============================================================================
# 规则
# ============================================================================

Rule = Callable[[pd.DataFrame, str], tuple[bool, float, str]]

_KO_PATTERN = re.compile(r"^K\d{5}$")
_ALPHA_METRICS = {"shannon", "simpson", "chao1", "observed_species", "invsimpson",
                   "ace", "richness", "evenness"}


def _has_col(df: pd.DataFrame, *names: str) -> bool:
    lower = {c.lower() for c in df.columns}
    return all(n.lower() in lower for n in names)


def _any_col(df: pd.DataFrame, *names: str) -> bool:
    lower = {c.lower() for c in df.columns}
    return any(n.lower() in lower for n in names)


def _numeric_columns(df: pd.DataFrame) -> list[str]:
    """返回能整体转为数值的列。"""
    out = []
    for c in df.columns:
        try:
            nan_ratio = pd.to_numeric(df[c], errors="coerce").isna().mean()
            if nan_ratio < 0.05:
                out.append(c)
        except Exception:
            pass
    return out


def _rule_metadata(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    if _has_col(df, "SampleID", "Group") or _has_col(df, "Sample_ID", "Group"):
        return True, 0.95, "列含 SampleID + Group"
    return False, 0.0, ""


def _rule_distance_matrix(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """方阵 + 对称 + 首列与列头除首项外一致。"""
    if df.shape[0] != df.shape[1] - 1:
        return False, 0.0, ""

    sample_cols = df.columns[1:].tolist()
    first_col_values = df.iloc[:, 0].astype(str).tolist()
    if sample_cols != first_col_values:
        return False, 0.0, ""

    try:
        mat = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy()
    except Exception:
        return False, 0.0, ""
    if np.isnan(mat).any():
        return False, 0.0, ""
    if not np.allclose(mat, mat.T, atol=1e-6):
        return False, 0.0, ""
    if np.abs(np.diag(mat)).max() > 1e-6:
        return False, 0.0, ""   # 对角线必须为 0

    return True, 0.95, f"方阵 ({mat.shape[0]}×{mat.shape[0]}) + 对称 + 对角零"


def _rule_alpha_diversity(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    has_sample = _any_col(df, "SampleID", "Sample_ID", "#SampleID")
    metrics_hit = [c for c in df.columns if c.lower() in _ALPHA_METRICS]
    if has_sample and metrics_hit:
        return True, 0.9, f"含 Sample ID 列 + α 指数：{', '.join(metrics_hit)}"
    return False, 0.0, ""


def _rule_checkm_quality(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    if _has_col(df, "Completeness", "Contamination"):
        return True, 0.95, "列含 Completeness + Contamination"
    return False, 0.0, ""


def _rule_env_factors(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    if not _any_col(df, "SampleID", "Sample_ID"):
        return False, 0.0, ""
    if not _any_col(df, "Group"):
        return False, 0.0, ""
    num_cols = _numeric_columns(df)
    # 排除 Group（它也是数值-like 时别误计）
    non_meta = [c for c in num_cols if c.lower() not in ("group", "replicate")]
    if len(non_meta) < 2:
        return False, 0.0, ""
    return True, 0.85, f"SampleID + Group + {len(non_meta)} 数值环境因子列"


def _rule_ko_annotation_long(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """MAG × KO 长表（MAG, Gene_ID, KEGG_ko, Description 等）。"""
    has_mag = _any_col(df, "MAG", "Bin", "Bin_ID", "Genome")
    has_ko = _any_col(df, "KEGG_ko", "KEGG ko", "KO", "KO_id")
    if has_mag and has_ko:
        # 区别于宽表：应该**没有**多个数值样本列
        num_cols = _numeric_columns(df)
        # 丰度宽表会有 ≥3 个样本列；长表即使有丰度一般 ≤1
        if len(num_cols) <= 2:
            return True, 0.95, "长表：MAG + KEGG_ko（MAG × KO 注释）"
    return False, 0.0, ""


def _rule_keystone_species(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """Keystone 物种列表：MAG + 网络度量/选取理由。"""
    if not _any_col(df, "MAG", "Bin", "Bin_ID"):
        return False, 0.0, ""
    hints = ("Degree", "Betweenness", "Closeness", "Selection_Reason",
             "Keystone", "hub")
    hits = [h for h in hints if _any_col(df, h)]
    if hits:
        return True, 0.9, f"MAG + keystone 指标列（{', '.join(hits)}）"
    return False, 0.0, ""


def _rule_mag_taxonomy(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """MAG 分类表（GTDB 格式）。支持有/无 header 两种写法。"""
    # 有 header 版本
    if (_any_col(df, "MAG", "Bin", "Genome")
            and _any_col(df, "classification", "taxonomy", "lineage",
                         "Taxonomy")):
        return True, 0.95, "列含 MAG + classification（有 header）"
    # 无 header 版本：2 列，列名本身就是数据（首行当 header 了）
    if df.shape[1] == 2:
        col0_name = str(df.columns[0])
        col1_name = str(df.columns[1])
        # 列 1 名字或数据里多数含 d__/p__/... GTDB 前缀
        tax_re = re.compile(r"[dpocfgs]__")
        col1_vals = df.iloc[:, 0].astype(str).tolist() + df.iloc[:, 1].astype(str).tolist()
        tax_hit = sum(1 for v in col1_vals if tax_re.search(v))
        header_hit = 1 if tax_re.search(col0_name) or tax_re.search(col1_name) else 0
        if tax_hit / max(len(col1_vals), 1) > 0.5 or header_hit > 0:
            return True, 0.9, "2 列含 GTDB 分类前缀（无 header）"
    return False, 0.0, ""


def _rule_ko_abundance_wide(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """宽格式 KO 丰度表：某列叫 KEGG_ko，或首列（或前 2 列中某列）全是 K\d{5} 模式。"""
    if _any_col(df, "KEGG_ko", "KEGG ko"):
        # 需要至少有数值样本列
        num_cols = _numeric_columns(df)
        if len(num_cols) >= 2:
            return True, 0.95, "列含 KEGG_ko + 多个数值样本列"
    # 次级：首列或第二列值多数是 K##### 形式
    for col_idx in (0, 1):
        if col_idx >= df.shape[1]:
            continue
        col = df.iloc[:, col_idx].astype(str)
        ratio = col.str.match(_KO_PATTERN.pattern).mean()
        if ratio > 0.9:
            num_cols = _numeric_columns(df)
            if len(num_cols) >= 2:
                return True, 0.9, f"列 {df.columns[col_idx]!r} 含 KO 标识 ({ratio*100:.0f}% 匹配)"
    return False, 0.0, ""


def _rule_abundance_wide(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """丰度宽表：首列标识 + 其余列纯数值。

    置信度分层（避免 Genus.txt / Species.txt 和 abundance.tsv 冲突时默认选错）：
    - MAG 级（首列 Genome/MAG/Bin，或行值看起来像 Mx_*/MAG_*/bin_*）→ 0.95（高）
    - 分类级（首列 Taxonomy/OTU/Taxon）→ 0.88（中）
    - 都不是 → 不匹配
    """
    if df.shape[1] < 3 or df.shape[0] < 2:
        return False, 0.0, ""
    first_col = df.columns[0].lower()

    # 分类 Taxonomy 关键字
    is_taxon_header = first_col in ("taxonomy", "#otu id", "#otuid", "taxon")
    # MAG 关键字
    is_mag_header = first_col in ("genome", "mag", "bin", "mag_id", "genome_id")
    if not (is_taxon_header or is_mag_header):
        return False, 0.0, ""

    sample_cols = df.columns[1:]
    try:
        numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        nan_ratio = numeric.isna().to_numpy().mean()
    except Exception:
        return False, 0.0, ""

    if nan_ratio >= 0.05:
        return False, 0.0, ""

    # 用行值判断 MAG vs 分类：
    # - 分类级：值以 GTDB/SILVA 前缀开头（k__/p__/c__/o__/f__/g__/s__/d__），
    #   或含分号分隔的 lineage
    # - MAG 级：值以 Mx_/MAG_/bin./GCF_/GCA_/bin_ 开头
    first_vals = df.iloc[:min(20, len(df)), 0].astype(str).tolist()
    taxon_prefixes = ("k__", "p__", "c__", "o__", "f__", "g__", "s__", "d__")
    taxon_like = sum(
        1 for v in first_vals
        if v.startswith(taxon_prefixes) or ";" in v
    )
    mag_prefixes = ("Mx_", "MAG_", "MAG.", "bin_", "bin.", "Bin", "GCF_", "GCA_")
    mag_like = sum(1 for v in first_vals if v.startswith(mag_prefixes))

    n = len(first_vals) or 1
    if is_mag_header or mag_like >= n * 0.5:
        looks_like_mag = True
    elif taxon_like >= n * 0.5:
        looks_like_mag = False
    else:
        # 无法断定：沿用旧行为（低置信度标记为 ABUNDANCE_WIDE）
        looks_like_mag = False

    conf = 0.95 if looks_like_mag else 0.88
    hint = ("MAG 级" if looks_like_mag else "分类级")
    return True, conf, (
        f"{hint}丰度表：首列 {df.columns[0]!r} + "
        f"{len(sample_cols)} 个数值样本列"
    )


def _rule_gephi_nodes(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """Gephi 节点表：有 Id/MAG + Degree + Betweenness。"""
    cols_lower = {c.lower() for c in df.columns}
    has_id = bool(cols_lower & {"id", "mag", "name", "genome"})
    has_degree = "degree" in cols_lower
    has_between = bool(cols_lower & {"betweenness", "betweenness centrality"})
    if has_id and has_degree and has_between:
        return True, 0.95, "Gephi 节点表（Id + Degree + Betweenness）"
    return False, 0.0, ""


def _rule_gephi_edges(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    """Gephi 边表：有 Source + Target（+ 可选 Weight）。"""
    cols_lower = {c.lower() for c in df.columns}
    has_src = bool(cols_lower & {"source", "from", "mag1"})
    has_tgt = bool(cols_lower & {"target", "to", "mag2"})
    if has_src and has_tgt:
        return True, 0.95, "Gephi 边表（Source + Target）"
    return False, 0.0, ""


# 规则顺序：更严格/更独特的先匹配。metadata 和 env_factors 都要 SampleID+Group，
# env_factors 更严格（需要额外数值列），所以它应先于 metadata 优先级考察。
_RULES: list[tuple[FileType, Rule]] = [
    (FileType.ENV_FACTORS, _rule_env_factors),
    (FileType.CHECKM_QUALITY, _rule_checkm_quality),
    (FileType.ALPHA_DIVERSITY, _rule_alpha_diversity),
    (FileType.DISTANCE_MATRIX, _rule_distance_matrix),
    # Keystone / MAG 分类 / KO 长表要排在通用 KO 宽表和丰度宽表之前
    (FileType.KEYSTONE_SPECIES, _rule_keystone_species),
    # Gephi 在 keystone 之后（keystone 表也有 Degree+Betweenness 列）
    (FileType.GEPHI_NODES, _rule_gephi_nodes),
    (FileType.GEPHI_EDGES, _rule_gephi_edges),
    (FileType.KO_ANNOTATION_LONG, _rule_ko_annotation_long),
    (FileType.MAG_TAXONOMY, _rule_mag_taxonomy),
    (FileType.KO_ABUNDANCE_WIDE, _rule_ko_abundance_wide),
    (FileType.METADATA, _rule_metadata),
    (FileType.ABUNDANCE_WIDE, _rule_abundance_wide),
]


# ============================================================================
# 公共入口
# ============================================================================

def detect(source: Union[str, Path, io.BytesIO], filename: str | None = None) -> DetectionResult:
    if filename is None and isinstance(source, (str, Path)):
        filename = Path(source).name
    filename = filename or ""

    df, enc, sep = read_table(source)

    # 按 _RULES 的顺序走：第一个匹配即采纳（更具体的类型排在前）
    best_type = FileType.UNKNOWN
    best_conf = 0.0
    best_reason = ""
    for ftype, rule in _RULES:
        ok, conf, reason = rule(df, filename)
        if ok:
            best_type, best_conf, best_reason = ftype, conf, reason
            break

    return DetectionResult(
        file_type=best_type,
        confidence=best_conf,
        reasons=[best_reason] if best_reason else [],
        preview_df=df.head(5),
        encoding=enc,
        separator=sep,
    )
