#!/usr/bin/env python3
"""Guard that the local charged softmax consumer stays diagnostic-only."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_ordered_shape_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_ordered_shape_readout.json"


def test_ordered_shape_readout_is_diagnostic_only() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["proof_status"] == "diagnostic_only_softmax_consumer"
    assert payload["predictive_promotion_allowed"] is False
