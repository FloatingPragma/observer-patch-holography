#!/usr/bin/env python3
"""Validate the exact compare-only neutrino two-parameter adapter."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_two_parameter_exact_adapter.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"


def test_two_parameter_exact_adapter_hits_both_representative_central_splittings() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_neutrino_two_parameter_exact_adapter"
    assert payload["scope"] == "compare_only_two_parameter_segment_adapter"
    assert payload["promotable"] is False
    assert payload["exact_solution"]["tau_nu"] == pytest.approx(0.49717065800809745, abs=1.0e-15)
    assert payload["exact_solution"]["lambda_nu"] == pytest.approx(1.723958214811294, abs=1.0e-12)
    assert payload["exact_outputs"]["masses_eV"] == pytest.approx(
        [0.01745663294772044, 0.019484199595350048, 0.053081390655025595], abs=1.0e-15
    )
    assert payload["exact_outputs"]["delta_m_sq_eV2"]["21"] == pytest.approx(
        payload["reference_central_values"]["delta_m21_sq_eV2"], abs=1.0e-18
    )
    assert payload["exact_outputs"]["delta_m_sq_eV2"]["32"] == pytest.approx(
        payload["reference_central_values"]["delta_m32_sq_eV2"], abs=1.0e-18
    )
    assert payload["exact_outputs"]["delta_m_sq_eV2"]["31"] == pytest.approx(
        payload["reference_central_values"]["delta_m31_sq_eV2"], abs=1.0e-18
    )
    assert payload["exact_fit_residuals_eV2"]["21"] == pytest.approx(0.0, abs=1.0e-18)
    assert payload["exact_fit_residuals_eV2"]["32"] == pytest.approx(0.0, abs=1.0e-18)
    assert payload["exact_outputs"]["ratio_21_over_32"] == pytest.approx(
        payload["reference_central_values"]["ratio_21_over_32"], abs=1.0e-16
    )
