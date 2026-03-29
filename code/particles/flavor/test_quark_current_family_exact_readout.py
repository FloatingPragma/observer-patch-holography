#!/usr/bin/env python3
"""Validate the exact current-family quark readout witness."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPREAD_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_spread_map.py"
MEAN_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_sector_mean_split.py"
READOUT_SCRIPT = ROOT / "particles" / "flavor" / "derive_quark_current_family_exact_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "flavor" / "quark_current_family_exact_readout.json"


def test_quark_current_family_exact_readout_hits_reference_targets() -> None:
    subprocess.run([sys.executable, str(SPREAD_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(MEAN_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_quark_current_family_exact_readout"
    assert payload["proof_status"] == "current_family_exact_witness"
    assert payload["predicted_singular_values_u"] == payload["reference_targets_u"]
    assert payload["predicted_singular_values_d"] == payload["reference_targets_d"]
