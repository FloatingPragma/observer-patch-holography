#!/usr/bin/env python3
"""Validate the exact current-family quark readout witness."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys

import pytest


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
THEOREM_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_quadratic_readout_theorem.py"
READOUT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"


def test_quark_current_family_exact_readout_hits_reference_targets() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(THEOREM_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_exact_readout"
    assert payload["proof_status"] == "current_family_exact_witness"
    assert payload["readout_chain_status"] == "closed_within_current_family_scope"
    assert payload["supporting_readout_theorem"] == "oph_quark_current_family_quadratic_readout_theorem"
    assert payload["smallest_constructive_missing_object"] is None
    assert payload["predicted_singular_values_u"] == pytest.approx(payload["reference_targets_u"], rel=1.0e-12, abs=1.0e-15)
    assert payload["predicted_singular_values_d"] == pytest.approx(payload["reference_targets_d"], rel=1.0e-12, abs=1.0e-15)
