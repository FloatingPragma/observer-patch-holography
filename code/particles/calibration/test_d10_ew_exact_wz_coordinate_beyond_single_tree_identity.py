#!/usr/bin/env python3
"""Guard the D10 exact-W/Z coordinate shell beneath the unsplit tree identity."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE_PAIR_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_source_transport_pair.py"
POPULATION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_population_evaluator.py"
FIBERWISE_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.py"
OBSTRUCTION_SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_tau2_current_carrier_obstruction.py"
SCRIPT = ROOT / "particles" / "calibration" / "derive_d10_ew_exact_wz_coordinate_beyond_single_tree_identity.py"
OUTPUT = ROOT / "particles" / "runs" / "calibration" / "d10_ew_exact_wz_coordinate_beyond_single_tree_identity.json"


def test_d10_exact_wz_coordinate_shell_is_emitted() -> None:
    subprocess.run([sys.executable, str(SOURCE_PAIR_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(POPULATION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(FIBERWISE_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(OBSTRUCTION_SCRIPT)], check=True, cwd=ROOT)
    subprocess.run([sys.executable, str(SCRIPT)], check=True, cwd=ROOT)
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))

    assert payload["artifact"] == "oph_d10_ew_exact_wz_coordinate_beyond_single_tree_identity"
    assert payload["status"] == "open_current_carrier_insufficient"
    assert payload["depends_on_object"] == "EWFiberwisePopulationTreeLaw_D10"
    assert payload["current_carrier_obstruction_artifact"] == "oph_d10_ew_tau2_current_carrier_obstruction"
    assert payload["coordinate_symbol"] == "tau2_tree_exact"
    assert payload["tau2_tree_exact"] is None
    assert payload["next_residual_object_if_open"] == "delta_n_tree_exact"
    assert payload["direct_tau2_emission_blocked"] is True
    assert payload["tauY_from_single_tree_identity_formula"] == "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)"
    beta = payload["carrier_basis_scalar"]["beta_EW"]
    alpha_y = payload["carrier_basis_scalar"]["alphaY_mz"]
    alpha2 = payload["carrier_basis_scalar"]["alpha2_mz"]
    assert abs(beta - ((alpha2 - alpha_y) / (alpha2 + alpha_y))) < 1.0e-15
