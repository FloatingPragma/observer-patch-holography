#!/usr/bin/env python3
"""Validate the D11 fixed-ray point-statement artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
FORWARD_SEED_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py"
FORWARD_CERT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed_promotion_certificate.py"
EXACT_ADAPTER_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_reference_exact_adapter.py"
NO_GO_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_fixed_ray_no_go_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_fixed_ray_no_go_theorem.json"


def test_d11_fixed_ray_no_go_theorem_closes_cleanly() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SEED_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_CERT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_ADAPTER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(NO_GO_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_fixed_ray_no_go_theorem"
    assert payload["proof_status"] == "central_pair_off_fixed_ray_on_declared_linear_surface_within_one_sigma"
    assert payload["current_fixed_ray_branch"]["w_HT"] == 0.0
    assert payload["exact_compare_witness"]["w_HT_exact"] == pytest.approx(-0.00248687922025298, abs=1.0e-18)
    assert payload["fixed_ray_point_test"]["central_pair_on_fixed_ray"] is False
    assert payload["fixed_ray_point_test"]["central_pair_within_one_sigma_of_fixed_ray"] is True
    pull = payload["uncertainty_pull"]
    assert pull["sigma_w_HT"] == pytest.approx(0.00366, abs=2.0e-5)
    assert pull["abs_pull_sigma"] == pytest.approx(0.68, abs=0.01)
    assert pull["mt_on_fixed_ray_at_mH_central_gev"] == pytest.approx(172.52, abs=0.01)
    codomain = payload["codomain_dependence"]
    assert codomain["primary_codomain"]["mt_pole_summary_id"] == "Q007TP4"
    assert codomain["auxiliary_codomain"]["mt_summary_id"] == "Q007TP"
    assert abs(codomain["auxiliary_codomain"]["pull_sigma"]) < 1.0
    extension = payload["smallest_supported_extension"]
    assert extension["one_extra_scalar_beyond_fixed_ray"] is True
    assert extension["extension_role"] == "chosen_extension"
    assert extension["forced_by_data_on_declared_surface"] is False
