#!/usr/bin/env python3
"""Validate the intrinsic neutrino exact eta-map artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "neutrino" / "derive_intrinsic_neutrino_exact_eta_map.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_eta_map.json"


def test_intrinsic_eta_map_is_exact_once_eta_is_given() -> None:
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_intrinsic_neutrino_exact_eta_map"
    assert payload["theorem_surface_status"] == "intrinsic_builder_complete_exact"
    assert payload["selector_common_scale_invariant"] is True
    assert payload["pmns_status"] == "blocked_pending_shared_charged_lepton_left_basis"
    assert len(payload["masses_gev_sorted"]) == 3
    assert payload["cubic_root_audit_max_abs_gev2"] < 1.0e-30

