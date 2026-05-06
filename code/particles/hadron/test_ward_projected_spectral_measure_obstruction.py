#!/usr/bin/env python3
"""Tests for the Ward-projected hadronic spectral-measure obstruction."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys
import tempfile

from derive_ward_projected_spectral_measure_obstruction import build_obstruction


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "hadron" / "derive_ward_projected_spectral_measure_obstruction.py"


def test_current_local_artifacts_do_not_construct_rho_q() -> None:
    payload = build_obstruction()

    assert payload["artifact"] == "oph_ward_projected_hadronic_spectral_measure_obstruction"
    assert payload["status"] == "rho_Q_not_constructible_from_current_local_artifacts"
    assert payload["promotion_allowed"] is False
    assert payload["external_inputs_used"] is False

    state = payload["local_artifact_state"]
    assert state["rho_levels"]["level_points_count"] == 0
    assert state["rho_phase_shift"]["level_points_count"] == 0
    assert state["rho_phase_shift"]["resonance_populated"] is False
    assert state["full_unquenched_rho_scattering"]["correlation_matrices_populated"] is False
    assert state["full_unquenched_rho_scattering"]["principal_correlators_populated"] is False
    assert state["production_closure"]["production_dump_present"] is False


def test_obstruction_names_the_minimal_missing_source_object() -> None:
    payload = build_obstruction()
    missing = payload["minimal_source_object_still_missing"]

    assert missing["id"] == "oph_qcd_ward_projected_hadronic_spectral_measure"
    required = set(missing["strict_minimum_fields"])
    assert "nonempty finite-volume vector-channel level support" in required
    assert "Ward-projected residues/weights with current normalization" in required
    assert "pushforward rule from finite-volume data to rho_Q(s;P)" in required
    assert "quadrature and OPE/tail bound for the Thomson kernel" in required


def test_nonidentifiability_witness_changes_thomson_moment() -> None:
    payload = build_obstruction()
    proof = payload["theorem_grade_obstruction"]["proof_by_nonidentifiability"]
    completions = proof["two_positive_completions_same_projection"]

    assert completions[0]["moment"] == "1/6"
    assert completions[1]["moment"] == "1/12"
    assert proof["moments_differ"] is True


def test_cli_writes_obstruction_payload_without_external_input() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_path = pathlib.Path(tmp_dir) / "obstruction.json"
        subprocess.run(
            [sys.executable, str(SCRIPT), "--output", str(out_path)],
            cwd=ROOT,
            check=True,
        )
        payload = json.loads(out_path.read_text(encoding="utf-8"))

    assert payload["external_inputs_used"] is False
    assert payload["forbidden_shortcuts"]["codata_or_nist_endpoint"].startswith("not inspected")
