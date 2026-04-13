"""文件类型识别引擎。

基于表头规则匹配 + 数据特征检测。通过向 _RULES 追加新规则即可扩展支持的类型。
"""
from __future__ import annotations

import io
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, Union

import chardet
import pandas as pd


class FileType(str, Enum):
    METADATA = "metadata"
    ABUNDANCE_WIDE = "abundance_wide"
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
    # 常见别名归一
    return {"ascii": "utf-8", "gb2312": "gbk"}.get(enc, enc)


def _sniff_separator(text_sample: str) -> str:
    """只看第一行，比较 tab 和逗号数量。"""
    first_line = text_sample.split("\n", 1)[0]
    return "\t" if first_line.count("\t") >= first_line.count(",") else ","


def read_table(source: Union[str, Path, io.BytesIO], n_preview: int = 5) -> tuple[pd.DataFrame, str, str]:
    """从路径或上传的 BytesIO 读取表格，返回 (df, encoding, separator)。

    允许上游传文件路径字符串、Path，或 Streamlit 的 UploadedFile（鸭子类型 .read）。
    """
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
# 规则：每条规则接收 df + filename，返回 (是否匹配, 置信度, 解释)
# ============================================================================

Rule = Callable[[pd.DataFrame, str], tuple[bool, float, str]]


def _rule_metadata(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    cols = {c.lower() for c in df.columns}
    has_sample = any(c in cols for c in ("sampleid", "sample_id", "sample id", "#sampleid"))
    has_group = "group" in cols
    if has_sample and has_group:
        return True, 0.95, "列含 SampleID + Group"
    return False, 0.0, ""


def _rule_abundance_wide(df: pd.DataFrame, filename: str) -> tuple[bool, float, str]:
    if df.shape[1] < 3 or df.shape[0] < 2:
        return False, 0.0, ""
    first_col = df.columns[0].lower()
    first_col_ok = first_col in ("taxonomy", "#otu id", "#otuid", "taxon", "genome")
    if not first_col_ok:
        return False, 0.0, ""

    # 后续列应该全部能解析为数值（样本丰度列）
    sample_cols = df.columns[1:]
    try:
        numeric = df[sample_cols].apply(pd.to_numeric, errors="coerce")
        nan_ratio = numeric.isna().to_numpy().mean()
    except Exception:
        return False, 0.0, ""

    if nan_ratio < 0.05:
        return True, 0.9, f"首列 {df.columns[0]!r} + {len(sample_cols)} 个数值样本列"
    return False, 0.0, ""


_RULES: dict[FileType, Rule] = {
    FileType.METADATA: _rule_metadata,
    FileType.ABUNDANCE_WIDE: _rule_abundance_wide,
}


# ============================================================================
# 公共入口
# ============================================================================

def detect(source: Union[str, Path, io.BytesIO], filename: str | None = None) -> DetectionResult:
    """识别文件类型。

    参数：
        source   — 文件路径、Path 对象，或带 .read() 的文件对象（如 Streamlit UploadedFile）。
        filename — 如果 source 不是路径（如 BytesIO），用这个名字做兜底识别（当前未使用，保留扩展）。

    返回：DetectionResult（含 preview_df 供 UI 显示）。
    """
    if filename is None and isinstance(source, (str, Path)):
        filename = Path(source).name
    filename = filename or ""

    df, enc, sep = read_table(source)

    best_type = FileType.UNKNOWN
    best_conf = 0.0
    best_reason = ""
    for ftype, rule in _RULES.items():
        ok, conf, reason = rule(df, filename)
        if ok and conf > best_conf:
            best_type, best_conf, best_reason = ftype, conf, reason

    return DetectionResult(
        file_type=best_type,
        confidence=best_conf,
        reasons=[best_reason] if best_reason else [],
        preview_df=df.head(5),
        encoding=enc,
        separator=sep,
    )
