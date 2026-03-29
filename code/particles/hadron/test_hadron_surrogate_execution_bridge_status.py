#!/usr/bin/env python3
"""Validate the surrogate hadron execution bridge status artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "hadron" / "derive_hadron_surrogate_execution_bridge_status.py"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "hadron_surrogate_execution_bridge_status.json"


def test_hadron_surrogate_bridge_stays_diagnostic_only() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_hadron_surrogate_execution_bridge_status"
    assert payload["status"] == "surrogate_hmc_execution_bridge_complete"
    assert payload["public_promotion_allowed"] is False
    assert payload["surrogate_execution"] is True
    assert payload["runtime_receipt_frozen"] == {"N_therm": 2048, "N_sep": 512}
    assert payload["surrogate_error_summary"]["all_six_channel_limits_exist"] is True
