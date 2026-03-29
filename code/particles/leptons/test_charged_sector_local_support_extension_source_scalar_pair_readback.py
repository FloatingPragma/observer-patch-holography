#!/usr/bin/env python3
"""Validate the charged support-extension source-scalar pair readback artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ETA_READBACK_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_eta_source_readback.py"
ENDPOINT_RATIO_BREAKER_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_endpoint_ratio_breaker.py"
COMPLETION_LAW_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_completion_law.py"
PAIR_READBACK_SCRIPT = ROOT / "particles" / "leptons" / "derive_charged_sector_local_support_extension_source_scalar_pair_readback.py"
OUTPUT = ROOT / "particles" / "runs" / "leptons" / "charged_sector_local_support_extension_source_scalar_pair_readback.json"


def test_charged_support_extension_source_scalar_pair_readback_exposes_eta_then_sigma_order() -> None:
    subprocess.run([sys.executable, str(COMPLETION_LAW_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ETA_READBACK_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(ENDPOINT_RATIO_BREAKER_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(PAIR_READBACK_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_charged_sector_local_support_extension_source_scalar_pair_readback"
    assert payload["proof_status"] == "source_scalar_pair_readback_formula_closed_values_open"
    assert payload["scalar_order"] == [
        "eta_source_support_extension_log_per_side",
        "sigma_source_support_extension_total_log_per_side",
    ]
    assert payload["next_single_residual_object"] == "eta_source_support_extension_log_per_side"
    assert payload["next_single_residual_object_after_eta"] == "sigma_source_support_extension_total_log_per_side"
    assert payload["eta_source_support_extension_log_per_side"] is None
    assert payload["sigma_source_support_extension_total_log_per_side"] is None

