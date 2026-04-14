"""RDA 排序图测试。"""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from envmeta.analysis import rda
from envmeta.export.figure_export import export_to_bytes

SAMPLE = Path(__file__).parent / "sample_data"


@pytest.fixture
def rda_inputs():
    ab = pd.read_csv(SAMPLE / "Species.txt", sep="\t")
    env = pd.read_csv(SAMPLE / "env_factors.txt", sep="\t")
    md = pd.read_csv(SAMPLE / "metadata.txt", sep="\t")
    return ab, env, md


def test_rda_smoke(rda_inputs):
    ab, env, md = rda_inputs
    r = rda.analyze(ab, env, md, {"n_permutations": 99})
    assert r.figure is not None
    assert not r.stats.empty
    assert {"variable", "type"}.issubset(r.stats.columns)


def test_rda_explained_variance_present(rda_inputs):
    ab, env, md = rda_inputs
    r = rda.analyze(ab, env, md, {"n_permutations": 99})
    evar = r.stats[r.stats["type"] == "explained_variance"]
    assert len(evar) == 2
    # 前两轴总和应 > 0（基本合理性）
    assert evar["value"].sum() > 0.0


def test_rda_mantel_factors_match_env(rda_inputs):
    ab, env, md = rda_inputs
    r = rda.analyze(ab, env, md, {"n_permutations": 99})
    mantel_rows = r.stats[r.stats["type"] == "mantel"]
    # env_factors.txt 有 pH, Eh, TOC, Total_As 四个数值因子
    assert set(mantel_rows["variable"]) == {"pH", "Eh", "TOC", "Total_As"}
    assert mantel_rows["mantel_p"].between(0, 1).all()


def test_rda_export(rda_inputs):
    ab, env, md = rda_inputs
    r = rda.analyze(ab, env, md, {"n_permutations": 99})
    png = export_to_bytes(r.figure, "png")
    pdf = export_to_bytes(r.figure, "pdf")
    assert png.startswith(b"\x89PNG")
    assert pdf.startswith(b"%PDF")
