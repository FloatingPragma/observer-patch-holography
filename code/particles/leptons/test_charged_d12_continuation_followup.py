#!/usr/bin/env python3
"""Validate the charged D12 continuation follow-up artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_d12_continuation_followup.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_d12_continuation_followup.json"


def test_charged_d12_continuation_followup_is_marked_diagnostic_only() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_charged_d12_continuation_followup"
    assert payload["status"] == "diagnostic_only_d12_continuation"
    assert payload["public_promotion_allowed"] is False
    assert payload["d12_continuation_pair"]["eta_source_support_extension_log_per_side"] < 0.0
    assert payload["d12_continuation_pair"]["sigma_source_support_extension_total_log_per_side"] > 0.0
    assert payload["immediate_followup_required"]["live_scalar_gap"].startswith("eta_source_support_extension")

