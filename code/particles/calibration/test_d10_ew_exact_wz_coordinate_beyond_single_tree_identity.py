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


def test_d10_exact_wz_coordinate_shell_is_emitted(tmp_path: pathlib.Path) -> None:
    canonical_inputs = (
        SOURCE_PAIR_SCRIPT.parent.parent / "runs/calibration/d10_ew_source_transport_pair.json",
        SOURCE_PAIR_SCRIPT.parent.parent / "runs/calibration/d10_ew_population_evaluator.json",
        SOURCE_PAIR_SCRIPT.parent.parent
        / "runs/calibration/d10_ew_fiberwise_population_tree_law_beneath_single_tree_identity.json",
        SOURCE_PAIR_SCRIPT.parent.parent
        / "runs/calibration/d10_ew_tau2_current_carrier_obstruction.json",
        OUTPUT,
    )
    canonical_bytes = {path: path.read_bytes() for path in canonical_inputs}

    source_pair = tmp_path / "source_pair.json"
    population = tmp_path / "population.json"
    fiberwise = tmp_path / "fiberwise.json"
    obstruction = tmp_path / "obstruction.json"
    output = tmp_path / "exact_wz.json"
    subprocess.run(
        [sys.executable, str(SOURCE_PAIR_SCRIPT), "--output", str(source_pair)],
        check=True,
        cwd=ROOT,
    )
    subprocess.run(
        [
            sys.executable,
            str(POPULATION_SCRIPT),
            "--source-pair",
            str(source_pair),
            "--output",
            str(population),
        ],
        check=True,
        cwd=ROOT,
    )
    subprocess.run(
        [
            sys.executable,
            str(FIBERWISE_SCRIPT),
            "--source-pair",
            str(source_pair),
            "--population",
            str(population),
            "--output",
            str(fiberwise),
        ],
        check=True,
        cwd=ROOT,
    )
    subprocess.run(
        [
            sys.executable,
            str(OBSTRUCTION_SCRIPT),
            "--source-pair",
            str(source_pair),
            "--population",
            str(population),
            "--fiberwise-tree-law",
            str(fiberwise),
            "--output",
            str(obstruction),
        ],
        check=True,
        cwd=ROOT,
    )
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--source-pair",
            str(source_pair),
            "--population",
            str(population),
            "--fiberwise-tree-law",
            str(fiberwise),
            "--tau2-obstruction",
            str(obstruction),
            "--output",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert all(path.read_bytes() == before for path, before in canonical_bytes.items())

    assert payload["artifact"] == "oph_d10_ew_exact_wz_coordinate_beyond_single_tree_identity"
    assert payload["status"] == "open_current_carrier_insufficient"
    assert payload["depends_on_object"] == "EWFiberwisePopulationTreeLaw_D10"
    assert payload["current_carrier_obstruction_artifact"] == "oph_d10_ew_tau2_current_carrier_obstruction"
    assert payload["coordinate_symbol"] == "tau2_tree_exact"
    assert payload["tau2_tree_exact"] is None
    assert payload["next_residual_object_if_open"] == "delta_n_tree_exact"
    assert payload["direct_tau2_emission_blocked"] is True
    assert "central W/Z pair" in payload["direct_tau2_emission_blocked_scope"]
    assert "additional neutral coordinates" in payload["direct_tau2_emission_blocked_scope"]
    assert payload["proof_status"] == "exact_interval_excludes_single_tau2_central_WZ_pair_on_current_carrier"
    assert "exact interval certificate" in payload["minimality_certificate"]["why_one_scalar_suffices_after_single_tree_identity"]
    assert payload["tauY_from_single_tree_identity_formula"] == "-(tau2_tree_exact + 2*eta_source) / (1 + 4*tau2_tree_exact^2)"
    beta = payload["carrier_basis_scalar"]["beta_EW"]
    alpha_y = payload["carrier_basis_scalar"]["alphaY_mz"]
    alpha2 = payload["carrier_basis_scalar"]["alpha2_mz"]
    assert abs(beta - ((alpha2 - alpha_y) / (alpha2 + alpha_y))) < 1.0e-15
