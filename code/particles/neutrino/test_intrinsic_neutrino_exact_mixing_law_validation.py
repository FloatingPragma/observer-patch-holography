#!/usr/bin/env python3
"""Validate the intrinsic neutrino exact mixing-law audit artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
ETA_SCRIPT = ROOT / "particles" / "neutrino" / "derive_intrinsic_neutrino_exact_eta_map.py"
VALIDATION_SCRIPT = ROOT / "particles" / "neutrino" / "derive_intrinsic_neutrino_exact_mixing_law_validation.py"
OUTPUT = ROOT / "particles" / "runs" / "neutrino" / "intrinsic_neutrino_exact_mixing_law_validation.json"


def test_intrinsic_mixing_validation_tracks_corrected_atmospheric_statement() -> None:
    subprocess.run([sys.executable, str(ETA_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(VALIDATION_SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_intrinsic_neutrino_exact_mixing_law_validation"
    assert payload["ordered_atmospheric_shift_is_linear"] is True
    assert payload["centroid_gap_is_first_order_invariant"] is True
    assert payload["collective_vector_overlap_from_delta"] > 0.999

