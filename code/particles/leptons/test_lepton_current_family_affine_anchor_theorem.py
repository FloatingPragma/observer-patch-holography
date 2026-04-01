#!/usr/bin/env python3
"""Validate the exact same-family charged affine-anchor theorem artifact."""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
QUADRATIC_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_quadratic_readout_theorem.py"
EXACT_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_exact_readout.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_affine_anchor_theorem.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_affine_anchor_theorem.json"


def test_lepton_current_family_affine_anchor_theorem() -> None:
    subprocess.run([sys.executable, str(QUADRATIC_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_lepton_current_family_affine_anchor_theorem"
    assert payload["proof_status"] == "closed_current_family_affine_anchor"
    assert payload["theorem_scope"] == "current_family_only"
    assert payload["current_family_geometric_mean"]["value"] == pytest.approx(
        math.exp(payload["current_family_affine_anchor"]["value"]),
        rel=1.0e-12,
        abs=1.0e-15,
    )
    assert payload["reconstructed_from_affine_anchor"] == pytest.approx(
        payload["predicted_singular_values_abs"],
        rel=1.0e-12,
        abs=1.0e-15,
    )
