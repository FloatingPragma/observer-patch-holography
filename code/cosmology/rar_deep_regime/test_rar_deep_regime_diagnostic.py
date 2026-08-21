"""Tests for the SPARC deep-regime diagnostic."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import rar_deep_regime_diagnostic as rar  # noqa: E402


def test_tables_parse() -> None:
    t1 = rar.read_table1()
    t2 = rar.read_table2()
    assert len(t1) == 175
    assert len(t2) == 3391
    assert t1[0]["name"] == "CamB"
    assert t2[0]["rad_kpc"] == pytest.approx(0.16)


def test_synthetic_deep_law_is_recovered() -> None:
    rng = np.random.default_rng(1)
    a0 = 1.2e-10
    g_bar = 10 ** rng.uniform(-13, -9, 4000)
    g_obs = np.sqrt(g_bar * a0) * 10 ** rng.normal(0, 0.02, g_bar.size)
    fit = rar.deep_fit(g_obs, g_bar)
    assert fit["a0_fixed_exponent_m_s2"] == pytest.approx(a0, rel=0.02)
    assert fit["free_exponent"] == pytest.approx(0.5, abs=0.02)


def test_committed_receipt_matches_code() -> None:
    path = HERE / "receipts" / "sparc_deep_regime_diagnostic.json"
    receipt = json.loads(path.read_text(encoding="utf-8"))
    fresh = rar.run()
    assert receipt["physical_claim"] is False
    assert receipt["source_derived_output"] is False
    assert receipt["sample"]["data_sha256"] == fresh["sample"]["data_sha256"]
    for key in ("deep_radial_acceleration", "baryonic_tully_fisher"):
        assert receipt[key]["a0_fixed_exponent_m_s2"] == pytest.approx(
            fresh[key]["a0_fixed_exponent_m_s2"], rel=1e-9
        )
    assert receipt["additive_all_gradient_extension"]["additive_form_passes_solar_system"] is False


def test_quadrupole_crosscheck_matches_published_values() -> None:
    import qumond_quadrupole_crosscheck as qc

    out = qc.run(epsrel=1e-8)
    fl = out["field_law"]
    assert fl["simple_function_blanchet_novak_inputs"]["relative_difference"] < 0.05
    assert fl["rar_function_park_inputs"]["relative_difference"] < 1e-3
    assert fl["rar_function_park_inputs"]["pull_vs_cassini_sigma"] > 10
    dens = out["density_formulation"]
    assert dens["quadrupole_s2"] == 0.0
    assert dens["tidal_scale_s2"] < 1e-30
    path = HERE / "receipts" / "qumond_quadrupole_crosscheck.json"
    receipt = json.loads(path.read_text(encoding="utf-8"))
    assert receipt["physical_claim"] is False
    assert receipt["field_law"]["rar_function_park_inputs"]["Q2_s2"] == pytest.approx(
        fl["rar_function_park_inputs"]["Q2_s2"], rel=1e-4
    )
