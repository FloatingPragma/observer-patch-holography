#!/usr/bin/env python3
"""Validate the closed current-family charged quadratic readout theorem."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_quadratic_readout_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_quadratic_readout_theorem.json"


def test_lepton_current_family_quadratic_readout_theorem_closes_exact_readout() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_lepton_current_family_quadratic_readout_theorem"
    assert payload["proof_status"] == "closed_current_family_exact_readout"
    assert payload["quadratic_basis"]["invertible"] is True
    assert payload["quadratic_basis"]["minor_determinant"] != pytest.approx(0.0, abs=1.0e-12)
    assert payload["smallest_constructive_missing_object"] is None
    assert payload["predicted_singular_values_abs"] == pytest.approx(payload["reference_targets"], rel=1.0e-12, abs=1.0e-15)
