#!/usr/bin/env python3
"""Validate the D11 fixed-ray no-go theorem artifact."""

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
    assert payload["proof_status"] == "closed_no_go_on_current_one_scalar_fixed_ray"
    assert payload["current_fixed_ray_branch"]["w_HT"] == 0.0
    assert payload["exact_compare_witness"]["w_HT_exact"] == pytest.approx(-0.0003857630977715052, abs=1.0e-18)
    assert payload["fixed_ray_obstruction"]["current_exact_pair_reachable"] is False
    assert payload["smallest_honest_extension"]["one_extra_scalar_beyond_fixed_ray"] is True

