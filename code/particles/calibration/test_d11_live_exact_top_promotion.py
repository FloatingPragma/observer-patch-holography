#!/usr/bin/env python3
"""Validate the exact D11 top-side promotion theorem artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECLARED_SURFACE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_declared_calibration_surface.py"
TOP_SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_live_exact_top_promotion.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_live_exact_top_promotion.json"


def test_d11_live_exact_top_promotion_hits_exact_top_without_inverse_adapter() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TOP_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_live_exact_top_promotion"
    assert payload["theorem_id"] == "D11LiveForwardExactTopPromotion"
    assert payload["proof_status"] == "closed_target_anchored_live_exact_top_promotion"
    assert payload["mass_readout"]["exact_residual_gev"] == pytest.approx(0.0, abs=1.0e-12)
    assert "promotion_of_the_old_fixed_ray_as_exact_pair" in payload["strictly_not_claimed"]

