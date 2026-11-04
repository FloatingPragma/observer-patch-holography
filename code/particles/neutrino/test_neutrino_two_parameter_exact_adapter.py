#!/usr/bin/env python3
"""Validate the exact compare-only neutrino two-parameter adapter."""

from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_neutrino_two_parameter_exact_adapter.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_two_parameter_exact_adapter.json"


def _load_module():
    spec = importlib.util.spec_from_file_location("derive_neutrino_two_parameter_exact_adapter", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _solved_surface():
    """Solve the selector coordinate the way the producer does, without the file."""
    module = _load_module()
    inputs = module.load_segment_inputs(
        certificate=json.loads(module.DEFAULT_CERTIFICATE.read_text(encoding="utf-8")),
        cocycle=json.loads(module.DEFAULT_COCYCLE.read_text(encoding="utf-8")),
        phase_source=json.loads(module.DEFAULT_PHASE_SOURCE.read_text(encoding="utf-8")),
    )
    target_ratio = (
        module.PDG_2025_NO_CENTRAL["delta_m21_sq_eV2"] / module.PDG_2025_NO_CENTRAL["delta_m32_sq_eV2"]
    )
    solved = module._solve_tau_exact_ratio(**inputs, target_ratio=target_ratio)

    def ratio_at(tau_nu: float) -> float:
        return module._surface_at_tau(tau_nu=tau_nu, **inputs)["ratio_21_over_32"]

    return module, solved, ratio_at


def test_two_parameter_exact_adapter_hits_both_representative_central_splittings() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    module = _load_module()
    assert payload["artifact"] == "oph_neutrino_two_parameter_exact_adapter"
    assert payload["scope"] == "compare_only_two_parameter_segment_adapter"
    assert payload["promotable"] is False
    assert payload["proof_chain_role"] == "diagnostic_target_fit_only"
    assert payload["must_not_feed_back"] is True
    # tau_nu is the root of a ratio built from eigenvalue differences; the
    # producer declares the resolution it can reach instead of claiming digits
    # below the binary64 noise floor of the objective.
    assert payload["exact_solution"]["tau_nu"] == pytest.approx(0.49717065800809745, abs=module.TAU_RESOLUTION)
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


def test_steps_thousands_of_times_below_the_declared_tolerance_are_invisible() -> None:
    """The old 1e-15 assertion measured nothing: the objective cannot see it."""
    module, solved, ratio_at = _solved_surface()
    tau = solved["tau_nu"]
    reference = ratio_at(tau)
    for step in (1.0e-15, 1.0e-14, 1.0e-13):
        assert abs(ratio_at(tau + step) - reference) <= module.RATIO_NOISE
        assert abs(ratio_at(tau - step) - reference) <= module.RATIO_NOISE
    # The declared tolerance is not slack: at that step the objective has moved
    # by no more than a hundred noise floors.
    assert abs(ratio_at(tau + module.TAU_RESOLUTION) - reference) <= 1.0e-13
    assert module.TAU_RESOLUTION >= 1.0e-12


def test_declared_tau_tolerance_follows_the_measured_slope() -> None:
    """TAU_RESOLUTION is the noise floor divided by the local objective slope."""
    module, solved, ratio_at = _solved_surface()
    tau = solved["tau_nu"]
    probe = 1.0e-11
    slope = abs(ratio_at(tau + probe) - ratio_at(tau - probe)) / (2.0 * probe)
    assert 1.0e-4 < slope < 1.0e-3
    assert module.RATIO_NOISE / slope <= module.TAU_RESOLUTION
