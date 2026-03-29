#!/usr/bin/env python3
"""Guard the D11 forward-seed promotion certificate."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FORWARD_SEED_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py"
CERTIFICATE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed_promotion_certificate.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed_promotion_certificate.json"
MISSING_CERT = ROOT / "particles" / "runs" / "calibration" / "_missing_d11_forward_seed_promotion_certificate.json"


def test_d11_forward_seed_promotion_certificate_closes_live_forward_path() -> None:
    subprocess.run([sys.executable, str(FORWARD_SEED_SCRIPT), "--promotion-certificate", str(MISSING_CERT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(CERTIFICATE_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_forward_seed_promotion_certificate"
    assert payload["status"] == "closed"
    assert payload["predictive_promotion_allowed"] is True
    assert payload["seed_equality_certificate"]["residual_abs"] == 0.0
    assert payload["fixed_ray_wedge_vanishing_certificate"]["wedge_value"] == 0.0
    assert payload["smallest_predictive_missing_object"] is None
