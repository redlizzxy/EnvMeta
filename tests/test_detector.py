"""文件识别引擎测试。"""
from pathlib import Path

import pytest

from envmeta.file_manager.detector import FileType, detect

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.mark.parametrize("filename, expected", [
    ("metadata.txt", FileType.METADATA),
    ("Phylum.txt", FileType.ABUNDANCE_WIDE),
    ("Genus.txt", FileType.ABUNDANCE_WIDE),
    ("Species.txt", FileType.ABUNDANCE_WIDE),
    ("alpha.txt", FileType.ALPHA_DIVERSITY),
    ("beta_bray.txt", FileType.DISTANCE_MATRIX),
    ("quality_report.tsv", FileType.CHECKM_QUALITY),
    ("env_factors.txt", FileType.ENV_FACTORS),
    ("abundance.tsv", FileType.ABUNDANCE_WIDE),
    ("kegg_target_only.tsv", FileType.KO_ANNOTATION_LONG),   # S2.5-6
    ("keystone_species.txt", FileType.KEYSTONE_SPECIES),     # S2.5-6
    ("mag_taxonomy_labels.tsv", FileType.MAG_TAXONOMY),      # S2.5-6
    ("ko_tpm.spf", FileType.KO_ABUNDANCE_WIDE),
])
def test_detect_sample_files(filename, expected):
    r = detect(SAMPLE / filename)
    assert r.file_type == expected, (
        f"{filename}: 期望 {expected.value}，得到 {r.file_type.value}（{r.reasons}）"
    )


def test_detect_confidence_for_known():
    r = detect(SAMPLE / "metadata.txt")
    assert r.confidence >= 0.9
    assert r.preview_df is not None
    assert r.preview_df.shape[0] <= 5


def test_abundance_wide_confidence_split_mag_vs_taxon():
    """回归测试：MAG 级丰度表 conf=0.95，TAXON 级 conf=0.88。

    用于 `_first_mag_abundance()` / `_first_taxon_abundance()` 的语义判定。
    跨平台一致性（Windows / Linux）不依赖 session_state 插入顺序。
    """
    mag = detect(SAMPLE / "abundance.tsv")
    assert mag.file_type == FileType.ABUNDANCE_WIDE
    assert mag.confidence >= 0.92, f"MAG 级 conf 应 ≥0.92，实际 {mag.confidence}"

    for taxon_file in ("Genus.txt", "Phylum.txt", "Species.txt"):
        t = detect(SAMPLE / taxon_file)
        assert t.file_type == FileType.ABUNDANCE_WIDE
        assert t.confidence < 0.92, (
            f"{taxon_file} 应为 TAXON 级 (conf<0.92)，实际 {t.confidence}"
        )
