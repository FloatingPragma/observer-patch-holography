#!/usr/bin/env python3
"""Validate the compare-only exact D11 Higgs/top reference adapter."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_reference_exact_adapter.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_reference_exact_adapter.json"


def test_d11_reference_exact_adapter_hits_canonical_targets() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d11_reference_exact_adapter"
    assert payload["scope"] == "compare_only_inverse_slice"
    assert payload["promotable"] is False
    assert payload["predicted_outputs"]["mH_gev"] == pytest.approx(payload["exact_reference_targets"]["mH_gev"], abs=1.0e-12)
    assert payload["predicted_outputs"]["mt_pole_gev"] == pytest.approx(
        payload["exact_reference_targets"]["mt_pole_gev"], abs=1.0e-12
    )
    assert payload["exact_fit_residuals_gev"]["mH_gev"] == pytest.approx(0.0, abs=1.0e-12)
    assert payload["exact_fit_residuals_gev"]["mt_pole_gev"] == pytest.approx(0.0, abs=1.0e-12)
