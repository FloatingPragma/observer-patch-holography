#!/usr/bin/env python3
"""Validate the compare-only D11 top exactifier artifact."""

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


def test_d11_top_exactifier_preserves_fit_without_promotion() -> None:
    subprocess.run([sys.executable, str(DECLARED_SURFACE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(TOP_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_live_exact_top_promotion"
    assert payload["theorem_id"] == "D11TargetAnchoredTopExactifier"
    assert payload["proof_status"] == "compare_only_target_anchored_exactifier"
    assert payload["status"] == "compare_only"
    assert payload["target_ancestry"]["target_values_consumed"] is True
    assert payload["native_source_emission"] is False
    assert payload["source_surface_promotable"] is False
    assert payload["predictive_promotion_allowed"] is False
    assert payload["public_surface_candidate_allowed"] is False
    assert payload["comparison_surface_allowed"] is True
    assert payload["mass_readout"]["exact_residual_gev"] == pytest.approx(0.0, abs=1.0e-12)
    assert "promotion_of_the_old_fixed_ray_as_exact_pair" in payload["strictly_not_claimed"]
    assert "predictive_top_promotion" in payload["strictly_not_claimed"]
