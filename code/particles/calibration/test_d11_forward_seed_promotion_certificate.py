#!/usr/bin/env python3
"""Guard the non-promoting D11 fixed-ray algebra certificate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
FORWARD_SEED_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py"
CERTIFICATE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed_promotion_certificate.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"
MISSING_CERT = ROOT / "particles" / "runs" / "calibration" / "_missing_d11_forward_seed_promotion_certificate.json"


def test_d11_forward_seed_certificate_closes_algebra_not_source_path() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SEED_SCRIPT), "--promotion-certificate", str(MISSING_CERT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(CERTIFICATE_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_forward_seed_promotion_certificate"
    assert payload["status"] == "fixed_ray_algebra_closed"
    assert payload["proof_status"] == "fixed_ray_algebra_closed_on_declared_surface"
    assert payload["fixed_ray_algebra_closed"] is True
    assert payload["fixed_ray_branch_closed"] is False
    assert payload["source_surface_promotable"] is False
    assert payload["predictive_promotion_allowed"] is False
    assert payload["public_surface_candidate_allowed"] is False
    assert payload["forward_path_closed"] is False
    assert payload["seed_equality_certificate"]["residual_abs"] == 0.0
    assert payload["fixed_ray_wedge_vanishing_certificate"]["wedge_value"] == 0.0
    assert payload["smallest_predictive_missing_object"] == "source_emitted_higgs_yukawa_fj_packet"
