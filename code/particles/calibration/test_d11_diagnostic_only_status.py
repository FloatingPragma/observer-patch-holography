#!/usr/bin/env python3
"""Guard that the local D11 sidecar remains diagnostic-only."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "calibration" / "derive_d11_critical_surface_readout.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d11_critical_surface_readout.json"


def test_d11_sidecar_is_diagnostic_only() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["predictive_status"] == "diagnostic_only_until_forward_seed_artifact_exists"
    assert payload["predictive_promotion_allowed"] is False
    assert payload["readout_kernel"]["exact_center_promotion"]["status"] == "diagnostic_only"
