#!/usr/bin/env python3
"""Validate the exact D10 mass-pair chart artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
FIBERWISE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_tau2_current_carrier_obstruction.py"
EXACT_WZ_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py"
CHART_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_mass_pair_chart_current_carrier.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_mass_pair_chart_current_carrier.json"


def test_d10_exact_mass_pair_chart_current_carrier_reduces_to_selector_gap() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FIBERWISE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(EXACT_WZ_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(CHART_SCRIPT)], check=True, cwd=ROOT)

    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_exact_mass_pair_chart_current_carrier"
    assert payload["status"] == "closed_smaller_primitive"
    assert payload["next_single_residual_object"] == "EWExactMassPairSelector_D10"
    assert payload["current_selector_pullback"]["zero_set"]["tau2_tree_exact"] == 0.0
    assert payload["current_selector_pullback"]["zero_set"]["delta_n_tree_exact"] == 0.0
    assert payload["local_bijectivity_certificate"]["third_coordinate_needed"] is False
