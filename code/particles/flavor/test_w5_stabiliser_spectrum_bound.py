#!/usr/bin/env python3
"""Tests for the W5 residual-rotation stabiliser-spectrum receipt."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_w5_stabiliser_spectrum_bound.py"
ARTIFACT = ROOT / "particles" / "runs" / "flavor" / "w5_stabiliser_spectrum_bound.json"


def _run_to(path: pathlib.Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(path)],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    return json.loads(path.read_text(encoding="utf-8"))


def test_w5_stabiliser_loader_is_byte_exact(tmp_path: pathlib.Path) -> None:
    rebuilt = tmp_path / ARTIFACT.name
    _run_to(rebuilt)
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()
    assert b"\r\n" not in rebuilt.read_bytes()


def test_c3_and_c5_invariance_force_axial_double_degeneracy(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "w5.json")
    checks = payload["canonical_rotation_checks"]
    for group in ("C3", "C5"):
        record = checks[group]
        assert record["linear_constraint_rank_in_W5"] == 4
        assert record["fixed_locus_dimension"] == 1
        assert record["nonzero_projective_dimension"] == 0
        assert record["general_invariant_form"] == "alpha*(n*n^T-I/3)"
        assert record["spectrum"] == ["2*alpha/3", "-alpha/3", "-alpha/3"]
        assert record["eigenvalue_multiplicities"] == [1, 2]
        assert record["simple_spectrum_possible"] is False
        assert (
            record["fixed_block_determinants"]["transverse_vector_det_R_minus_I"]["numeric"]
            > 0
        )
        assert (
            record["fixed_block_determinants"]["planar_spin2_det_R2_minus_I"]["numeric"]
            > 0
        )


def test_c2_fixed_locus_has_three_linear_and_two_projective_dimensions(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "w5.json")
    record = payload["canonical_rotation_checks"]["C2"]
    assert record["linear_constraint_rank_in_W5"] == 2
    assert record["fixed_locus_dimension"] == 3
    assert record["nonzero_projective_dimension"] == 2
    assert record["free_coordinates"] == ["a", "b", "d"]
    assert record["simple_spectrum_possible"] is True
    assert record["simple_spectrum_witness"]["spectrum"] == [1, 2, -3]
    assert (
        record["fixed_block_determinants"]["planar_spin2_det_R2_minus_I"]["exact"]
        == "0"
    )


def test_klein_four_fixed_locus_has_one_projective_parameter(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "w5.json")
    record = payload["canonical_rotation_checks"]["V4"]
    assert record["linear_constraint_rank_in_W5"] == 3
    assert record["fixed_locus_dimension"] == 2
    assert record["nonzero_projective_dimension"] == 1
    assert record["general_invariant_form"] == "diag(a,d,-a-d)"
    assert record["simple_spectrum_possible"] is True
    lattice = payload["subgroup_lattice_argument"]
    assert lattice["exhaustive"] is True
    assert set(lattice["conjugacy_classes_of_subgroups"]) == {
        "1", "C2", "C3", "V4", "C5", "S3", "D5", "A4", "A5"
    }
    dims = lattice["fixed_locus_dimension_by_subgroup"]
    assert dims["V4"] == 2
    assert all(dims[g] <= 1 for g in ("C3", "C5", "S3", "D5", "A4", "A5"))


def test_conclusion_requires_potential_selection_without_universal_no_go(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "w5.json")
    result = payload["w5_orb_reclassification"]
    scope = payload["scope"]
    assert result["current_gap_type"] == "screen_derived_potential_selection_required"
    assert result["symmetry_alone_sufficient"] is False
    assert result["screen_derived_potential_required"] is True
    assert scope["universal_impossibility_claimed"] is False
    assert "simple_spectrum_points_with_C2_or_trivial_stabiliser" in scope["not_excluded"]
    assert "selection_by_a_specific_screen_derived_A5_invariant_potential" in scope["not_excluded"]
    assert scope["literature_claims_used"] is False
