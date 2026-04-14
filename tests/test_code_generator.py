"""代码生成器测试。"""
import ast

import pytest

from envmeta.export.code_generator import SUPPORTED, generate

_ALL_FILES = {
    "abundance": "data/abundance.tsv",
    "metadata": "data/metadata.tsv",
    "distance": "data/beta_bray.txt",
    "ko_abundance": "data/ko_tpm.spf",
    "alpha": "data/alpha.txt",
    "env_factors": "data/env_factors.txt",
}


@pytest.mark.parametrize("analysis_id", sorted(SUPPORTED))
def test_generate_is_valid_python(analysis_id):
    src = generate(analysis_id, _ALL_FILES, {"top_n": 10, "width_mm": 160})
    ast.parse(src)


@pytest.mark.parametrize("analysis_id", sorted(SUPPORTED))
def test_generate_contains_expected_calls(analysis_id):
    src = generate(analysis_id, _ALL_FILES, {"top_n": 10})
    assert f"{analysis_id}.analyze(" in src
    assert "export_figure(" in src
    assert "to_csv(" in src


def test_generate_strips_private_keys():
    src = generate("stackplot", _ALL_FILES,
                   {"top_n": 10, "_raw_means": "should_not_appear"})
    assert "_raw_means" not in src
    assert "top_n" in src


def test_generate_unsupported_raises():
    with pytest.raises(ValueError, match="不支持"):
        generate("nope", _ALL_FILES, {})


def test_generate_output_base_used():
    src = generate("pcoa", _ALL_FILES, {}, output_base="my_pcoa")
    assert 'my_pcoa.pdf' in src
    assert 'my_pcoa_stats.tsv' in src
