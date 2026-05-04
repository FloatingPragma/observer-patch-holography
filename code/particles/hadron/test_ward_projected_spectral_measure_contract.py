#!/usr/bin/env python3
"""Smoke tests for the Ward-projected spectral-measure contract."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "hadron" / "derive_ward_projected_spectral_measure_contract.py"
OUTPUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_contract.json"
SCHEMA = ROOT / "particles" / "hadron" / "ward_projected_spectral_measure.schema.json"


def test_contract_emits_constructive_spectral_measure_target() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_ward_projected_spectral_measure_contract"
    assert payload["constructive_next_artifact"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    assert payload["promotion_allowed"] is False
    assert "stable_channel_only_backend_export" in payload["forbidden_promotions"]
    assert schema["properties"]["artifact"]["const"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    required = set(schema["required"])
    assert "finite_volume_levels" in required
    assert "ward_projected_residues" in required
    assert "systematics" in required
