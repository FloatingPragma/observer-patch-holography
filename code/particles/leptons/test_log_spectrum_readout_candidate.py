#!/usr/bin/env python3
"""Guard the new charged ordered-log spectrum artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_lepton_log_spectrum_readout.py"
FORWARD_SCRIPT = ROOT / "particles" / "leptons" / "build_forward_charged_leptons.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "lepton_log_spectrum_readout.json"
FORWARD_OUTPUT = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"


def test_log_spectrum_readout_exports_non_softmax_candidate() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FORWARD_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    forward = json.loads(FORWARD_OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_lepton_log_spectrum_readout"
    assert payload["proof_status"] == "current_family_ordered_package_readback_candidate"
    assert payload["predictive_promotion_allowed"] is False
    assert payload["shape_log_shift_e"] == 0.0
    assert payload["shape_closed"] is True
    assert payload["smallest_constructive_missing_object"] == "eta_source_support_extension_log_per_side"
    assert len(payload["E_e_log_centered"]) == 3
    assert abs(sum(float(value) for value in payload["E_e_log_centered"])) <= 1.0e-12
    assert payload["sigma_e_total_log_per_side"] > 0.0
    assert forward["closure_state"] == "support_extension_missing"
    assert forward["shape_shift_missing"] is False
