#!/usr/bin/env python3
"""Tests for the exact icosahedral residual-axis angle-spectrum receipt."""

from __future__ import annotations

import json
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "particles" / "flavor" / "derive_icosahedral_axis_angle_spectrum.py"
ARTIFACT = ROOT / "particles" / "runs" / "flavor" / "icosahedral_axis_angle_spectrum.json"


def _run_to(path: pathlib.Path) -> dict:
    subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(path)],
        check=True,
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    return json.loads(path.read_text(encoding="utf-8"))


def test_icosahedral_axis_angle_loader_is_byte_exact(tmp_path: pathlib.Path) -> None:
    rebuilt = tmp_path / ARTIFACT.name
    _run_to(rebuilt)
    assert rebuilt.read_bytes() == ARTIFACT.read_bytes()
    assert b"\r\n" not in rebuilt.read_bytes()


def test_axis_construction_collapses_to_6_10_15_unoriented_axes(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "axis-spectrum.json")
    assert payload["construction"]["oriented_object_counts"] == {
        "vertices": 12,
        "triangular_faces": 20,
        "edges": 30,
    }
    assert payload["construction"]["unoriented_axis_counts"] == {
        "fivefold_vertex_axes": 6,
        "threefold_face_axes": 10,
        "twofold_edge_axes": 15,
        "total": 31,
    }


def test_full_distinct_acute_angle_spectrum_and_multiplicities_are_frozen(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "axis-spectrum.json")
    spectra = payload["pairwise_acute_angle_spectra"]
    expected_angles = {
        "fivefold_x_fivefold": [63.434948822922],
        "fivefold_x_threefold": [37.37736814065, 79.187683036428],
        "fivefold_x_twofold": [31.717474411461, 58.282525588539, 90.0],
        "threefold_x_threefold": [41.810314895779, 70.528779365509],
        "threefold_x_twofold": [
            20.905157447889,
            54.735610317245,
            69.094842552111,
            90.0,
        ],
        "twofold_x_twofold": [36.0, 60.0, 72.0, 90.0],
    }
    expected_pair_counts = {
        "fivefold_x_fivefold": 15,
        "fivefold_x_threefold": 60,
        "fivefold_x_twofold": 90,
        "threefold_x_threefold": 45,
        "threefold_x_twofold": 150,
        "twofold_x_twofold": 105,
    }
    for row_name, angles in expected_angles.items():
        row = spectra[row_name]
        assert [entry["acute_angle_degrees"] for entry in row["angles"]] == angles
        assert row["unoriented_pair_count"] == expected_pair_counts[row_name]
        assert sum(entry["multiplicity"] for entry in row["angles"]) == expected_pair_counts[row_name]
        assert all("cosine_squared_exact" in entry for entry in row["angles"])
        assert all(entry["acute_angle_exact"].startswith("acos(sqrt(") for entry in row["angles"])


def test_golden_angle_self_test_and_direct_axis_no_go_are_narrow(
    tmp_path: pathlib.Path,
) -> None:
    payload = _run_to(tmp_path / "axis-spectrum.json")
    method = payload["method_self_test"]
    minimum = payload["minimum_nonzero_axis_angle"]
    comparison = payload["compare_only_cabibbo_readback"]
    conclusion = payload["conclusion"]

    assert method["identity"] == "angle = arctan(1/phi)"
    assert method["expected_family_pair"] == "fivefold_x_twofold"
    assert method["hard_assertion_passed"] is True
    assert method["expected_angle_degrees"] == 31.717474411461
    assert minimum["angle_degrees"] == 20.905157447889
    assert minimum["family_pair"] == "threefold_x_twofold"
    assert minimum["sine"] == 0.356822089773
    assert comparison["input"]["value"] == 0.225
    assert comparison["input"]["display_value"] == "0.2250"
    assert comparison["input"]["standard_uncertainty"] == 0.0004
    assert comparison["input"]["published_notation"] == "0.2250 +/- 0.0004"
    assert comparison["input"]["determination"] == "Kmu2_decay_constant_ratio"
    assert comparison["input"]["fixture_path"] == (
        "code/particles/data/pdg_2024_vus_kmu2_fixture.json"
    )
    assert len(comparison["input"]["fixture_sha256"]) == 64
    assert comparison["input"]["source"]["publisher"] == "Particle Data Group"
    assert "revised April 2024" in comparison["input"]["source"]["edition"]
    assert comparison["input"]["role"] == "compare_only_not_used_to_construct_or_select_axes"
    assert comparison["one_standard_uncertainty_interval_abs_Vus"] == [
        0.2246,
        0.2254,
    ]
    assert comparison["plus_five_standard_uncertainties_abs_Vus"] == 0.227
    assert comparison[
        "direct_match_available_at_plus_five_standard_uncertainties"
    ] is False
    assert comparison["residual_mismatch_cabibbo_available"] is False
    assert conclusion["scope_is_narrow"] is True
    assert conclusion["not_excluded"] == [
        "all_A5_flavor_models",
        "spinorial_or_other_representations",
        "higher_order_symmetry_breaking",
        "arbitrary_overlap_geometry",
    ]
    assert conclusion["literature_claims_used_to_construct_axis_spectrum"] is False
    assert conclusion["external_comparison_coordinate_used"] is True
