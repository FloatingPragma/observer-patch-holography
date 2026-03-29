#!/usr/bin/env python3
"""Validate the reduced D10 electroweak source-pair artifact."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
FAMILY_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_observable_family.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_source_transport_pair.json"
def test_d10_source_pair_records_two_scalar_family_without_reference_slice() -> None:
    subprocess.run([sys.executable, str(FAMILY_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["artifact"] == "oph_d10_ew_source_transport_pair"
    assert payload["status"] == "selected_two_scalar_family"
    assert payload["family_coordinates"]["tau_Y_formula"] == "sigma_EW - eta_EW"
    assert payload["family_coordinates"]["tau_2_formula"] == "sigma_EW + eta_EW"
    assert payload["special_slices"]["current_one_seed_slice"]["sigma_EW"] == 0.0
    assert payload["predictive_population_closed"] is True
    assert payload["predictive_population_verdict"] == "closed_current_carrier_nonexact_running_quintet"
    assert payload["population_selector_status"] == "closed"
    assert payload["population_evaluator_object"] == "EWGaugeSourceTransportPairPopulationEvaluator_D10"
    assert payload["population_minimality_certificate"]["determinant_nonzero"] is True
    assert payload["population_minimality_certificate"]["third_scalar_needed"] is False
    assert payload["one_scalar_reduction_certificate"]["verdict"] == "no_one_scalar_residual_on_reopened_family"
    assert payload["population_nonuniqueness_certificate"]["verdict"] == "population_selector_now_closed_current_carrier"
    assert payload["smallest_constructive_missing_object"] == "EWSinglePostTransportTreeIdentity_D10"
    assert "exact_wz_candidate" not in payload
    assert payload["first_nonzero_oph_seed_trial"]["coherent_output_quintet"]["MW_pole"] > 0.0
    assert payload["first_nonzero_oph_seed_trial"]["coherent_output_quintet"]["MZ_pole"] > 0.0
