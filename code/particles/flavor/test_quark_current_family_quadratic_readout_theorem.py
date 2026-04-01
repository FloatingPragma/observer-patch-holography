#!/usr/bin/env python3
"""Validate the closed current-family quark quadratic readout theorem."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_quadratic_readout_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_quadratic_readout_theorem.json"


def test_quark_current_family_quadratic_readout_theorem_closes_exact_readout() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_quadratic_readout_theorem"
    assert payload["proof_status"] == "closed_current_family_exact_readout"
    assert payload["quadratic_basis"]["invertible"] is True
    assert payload["quadratic_basis"]["minor_determinant"] != pytest.approx(0.0, abs=1.0e-12)
    assert payload["smallest_constructive_missing_object"] is None
    assert payload["predicted_singular_values_u"] == pytest.approx(payload["reference_targets_u"], rel=1.0e-12, abs=1.0e-15)
    assert payload["predicted_singular_values_d"] == pytest.approx(payload["reference_targets_d"], rel=1.0e-12, abs=1.0e-15)
