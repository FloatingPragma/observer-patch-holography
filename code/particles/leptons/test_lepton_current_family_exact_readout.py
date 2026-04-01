#!/usr/bin/env python3
"""Validate the exact current-family charged-lepton readout witness."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
THEOREM_SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_quadratic_readout_theorem.py"
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_current_family_exact_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_current_family_exact_readout.json"


def test_lepton_current_family_exact_readout_hits_reference_targets() -> None:
    subprocess.run([sys.executable, str(THEOREM_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_lepton_current_family_exact_readout"
    assert payload["readout_chain_status"] == "closed_within_current_family_scope"
    assert payload["supporting_readout_theorem"] == "oph_lepton_current_family_quadratic_readout_theorem"
    assert payload["smallest_constructive_missing_object"] is None
    assert payload["predicted_singular_values_abs"] == pytest.approx(payload["reference_targets"], rel=1.0e-12, abs=1.0e-15)
