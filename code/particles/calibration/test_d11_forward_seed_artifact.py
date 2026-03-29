#!/usr/bin/env python3
"""Guard the compact D11 forward-seed artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_forward_seed.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_forward_seed.json"
MISSING_CERT = ROOT / "particles" / "runs" / "calibration" / "_missing_d11_forward_seed_promotion_certificate.json"


def test_d11_forward_seed_exports_single_scalar_candidate() -> None:
    subprocess.run([sys.executable, str(SCRIPT), "--promotion-certificate", str(MISSING_CERT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_forward_seed"
    assert payload["forward_seed_object"] == "sigma_D11_HT"
    assert payload["seed_status"] == "predictive_forward_law_candidate"
    assert payload["predictive_promotion_allowed"] is False
    assert payload["public_surface_candidate_allowed"] is True
    assert payload["source_seed_law"]["sigma_D11_HT_formula"] == "alpha_u * cos(2*theta_W0) / sqrt(pi)"
