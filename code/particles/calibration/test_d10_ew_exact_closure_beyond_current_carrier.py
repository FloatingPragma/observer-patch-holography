#!/usr/bin/env python3
"""Validate the D10 exact-closure artifact beyond the selected carrier."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
READOUT_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_readout.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
EXACT_CLOSURE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_closure_beyond_current_carrier.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_closure_beyond_current_carrier.json"


def test_d10_exact_closure_emits_split_readout_on_selected_carrier() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(READOUT_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_CLOSURE_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_exact_closure_beyond_current_carrier"
    assert payload["status"] == "closed"
    assert payload["exactness_surface_kind"] == "split_mass_neutral_readout"
    assert payload["derived_readout_compensator"]["kind"] == "derived_not_free"
    assert payload["exact_outputs"]["MW_pole"] > 0.0
    assert payload["exact_outputs"]["MZ_pole"] > 0.0
    assert payload["exact_outputs"]["alpha_em_eff_inv"] > 0.0
    assert payload["exact_outputs"]["sin2w_eff"] > 0.0
    assert payload["stronger_residual_object"] == "EWSinglePostTransportTreeIdentity_D10"
