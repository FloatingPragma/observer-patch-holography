from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import fz12_full_symbol_remainder_certificate as certificate
import verify_fz12_full_symbol_remainder_independent as independent


VERIFIER = certificate.HERE / "verify_fz12_full_symbol_remainder_independent.py"


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(independent.canonical_json_bytes(value))


def resign(value: dict, field: str) -> None:
    value[field] = independent.self_digest(value, field)


def independent_run(
    *,
    base_geometry: Path = independent.DEFAULT_BASE_GEOMETRY,
    source: Path = independent.DEFAULT_SOURCE,
    custody: Path = independent.DEFAULT_CUSTODY,
    remainder_lean: Path = independent.DEFAULT_REMAINDER_LEAN,
    receipt: Path = independent.DEFAULT_RECEIPT,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--base-geometry",
            str(base_geometry),
            "--source",
            str(source),
            "--custody",
            str(custody),
            "--remainder-lean",
            str(remainder_lean),
            "--receipt",
            str(receipt),
        ],
        cwd=certificate.HERE,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def mutated_receipt(tmp_path: Path, mutate) -> Path:
    receipt = json.loads(independent.DEFAULT_RECEIPT.read_text())
    mutate(receipt)
    resign(receipt, "receipt_sha256")
    path = tmp_path / "fz12_full_symbol_remainder_receipt.json"
    write_canonical(path, receipt)
    return path


def test_exact_geometry_moments_symbol_and_remainder() -> None:
    receipt = certificate.build_receipt()
    geometry = receipt["geometry_contract"]
    assert geometry == independent.expected_geometry_contract()
    assert geometry["source_direction_count"] == 30
    assert geometry["unoriented_axis_count"] == 15
    assert geometry["raw_carrier_difference_norm_squared"] == "4"
    assert geometry["unit_direction_norm_squared"] == "1"
    assert geometry["equal_weight"] == "1/5"
    assert geometry["q_domain"] == {
        "minimum": "0",
        "maximum": "1",
        "inclusive": True,
    }

    assert receipt["exact_moments"] == independent.expected_exact_moments()
    assert receipt["exact_symbol"] == independent.expected_exact_symbol()
    assert receipt["taylor_remainders"] == independent.expected_taylor_remainders()
    assert receipt["exact_moments"]["M8"] == "10/3 - (8/15) I6(n)"
    assert (
        receipt["taylor_remainders"]["R6"]["global_absolute_coefficient"]
        == "7/388800"
    )
    assert (
        receipt["taylor_remainders"]["R8"]["global_coefficient"]
        == "7/34992000"
    )


def test_independent_exact_geometry_reconstruction() -> None:
    result = independent.independent_geometry_checks()
    assert result == {
        "direction_count": 30,
        "projective_axis_count": 15,
        "raw_difference_norm_squared": "4",
        "M2": "10",
        "M4": "6",
        "M6": "30/7-(2/7)I6(n)",
        "M8": "10/3-(8/15)I6(n)",
        "controls": {
            "double_directed_labels_60": True,
            "projective_axes_15": True,
            "raw_norm2_differences": True,
        },
    }


def test_sharp_range_has_direct_exact_parent_evidence() -> None:
    receipt = certificate.build_receipt()
    assert [pin["path"] for pin in receipt["parent_pins"]] == [
        "code/a5_fingerprint/runtime/a5_multipole_fixed_point_receipt.json",
        "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json",
        "code/a5_fingerprint/runtime/fz12_custody_projection.json",
        "Lean/Screen/SeamCurrentEdge30Remainder.lean",
    ]
    base = json.loads(independent.DEFAULT_BASE_GEOMETRY.read_text())
    independent.verify_base_range_evidence(base)
    assert receipt["exact_moments"]["I6_unit_sphere_range"] == {
        "minimum": "-5/9",
        "maximum": "1",
    }


def test_lean_parent_is_fixed_and_admission_free() -> None:
    raw = independent.load_fixed_lean_parent(independent.DEFAULT_REMAINDER_LEAN)
    source = raw.decode("utf-8")
    assert "theorem seamMoment8_eq (k : Vec3)" in source
    assert "sorry" not in source
    assert "admit" not in source
    lean_pin = certificate.build_receipt()["parent_pins"][-1]
    assert lean_pin == {
        "path": "Lean/Screen/SeamCurrentEdge30Remainder.lean",
        "role": "kernel-checked exact eighth seam moment",
        "bytes": 2968,
        "sha256": (
            "sha256:e5d90cf9d84a4b4394355f9704114937268b2bfd32c478c4f93e6c53cd493928"
        ),
        "theorem": "OPH.SeamCurrentEdge30Remainder.seamMoment8_eq",
        "sorry_free": True,
    }


def test_positivity_is_spatial_and_physical_promotion_stays_open() -> None:
    receipt = certificate.build_receipt()
    assert (
        receipt["positivity_and_monotonicity"]
        == independent.expected_positivity()
    )
    assert receipt["positivity_and_monotonicity"]["lower_bound"] == (
        "(19/20) q^2"
    )
    assert receipt["positivity_and_monotonicity"]["group_velocity_claimed"] is False
    assert receipt["physical_boundary"] == independent.expected_physical_boundary()
    assert all(value is False for value in receipt["physical_boundary"].values())


def test_exposure_is_empty_and_no_score_or_verdict_exists() -> None:
    exposure = certificate.build_receipt()["exposure_boundary"]
    assert exposure == independent.expected_exposure_boundary()
    assert exposure["comparison_inputs"] == []
    assert all(value is False for key, value in exposure.items() if key != "comparison_inputs")


def test_required_negative_controls_are_fail_closed() -> None:
    controls = certificate.build_receipt()["negative_controls"]
    assert controls == independent.expected_negative_controls()
    assert controls["all_detectors_fired"] is True
    assert [row["control"] for row in controls["controls"]] == [
        "double_directed_labels_60",
        "projective_axes_15",
        "raw_norm2_differences",
        "fz11_rank_six_sign",
        "mutated_M8",
        "q_above_one",
        "unequal_weights",
        "extra_hop_shell",
    ]
    assert all(row["detector_fired"] is True for row in controls["controls"])


def test_committed_receipt_is_canonical_byte_exact() -> None:
    assert certificate.verify_committed_receipt() == certificate.build_receipt()
    assert independent.DEFAULT_RECEIPT.read_bytes() == independent.canonical_json_bytes(
        certificate.build_receipt()
    )


def test_independent_verifier_is_separate_and_passes() -> None:
    source = VERIFIER.read_text(encoding="utf-8")
    assert "import fz12_full_symbol_remainder_certificate" not in source
    assert "from fz12_full_symbol_remainder_certificate" not in source
    result = independent_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FZ12_FULL_SYMBOL_REMAINDER_VERIFY_PASS"


@pytest.mark.parametrize(
    ("parent_name", "self_field", "argument"),
    [
        ("a5_multipole_fixed_point_receipt.json", "receipt_sha256", "base_geometry"),
        ("seam_current_edge_prediction_receipt.json", "receipt_sha256", "source"),
        ("fz12_custody_projection.json", "projection_sha256", "custody"),
    ],
)
def test_resigned_json_parent_mutation_fails_fixed_raw_pin(
    tmp_path: Path, parent_name: str, self_field: str, argument: str
) -> None:
    original = independent.RUNTIME / parent_name
    parent = json.loads(original.read_text())
    parent["audit_mutation"] = True
    resign(parent, self_field)
    path = tmp_path / parent_name
    write_canonical(path, parent)
    kwargs = {argument: path}
    result = independent_run(**kwargs)
    assert result.returncode == 1
    assert "parent raw hash drift" in result.stderr


def test_lean_parent_text_mutation_fails_fixed_raw_pin(tmp_path: Path) -> None:
    source = independent.DEFAULT_REMAINDER_LEAN.read_text(encoding="utf-8")
    source = source.replace("(8 / 15 : ℝ)", "(7 / 15 : ℝ)", 1)
    path = tmp_path / "SeamCurrentEdge30Remainder.lean"
    # Preserve LF bytes on Windows so this mutation reaches the hash gate rather
    # than failing earlier because pathlib translated newlines to CRLF.
    path.write_bytes(source.encode("utf-8"))
    result = independent_run(remainder_lean=path)
    assert result.returncode == 1
    assert "remainder Lean raw hash drift" in result.stderr


def test_range_evidence_semantic_mutation_fails() -> None:
    parent = json.loads(independent.DEFAULT_BASE_GEOMETRY.read_text())
    mutated = copy.deepcopy(parent)
    next(
        row
        for row in mutated["critical_points"]["orbits"]
        if row["orbit"] == "face_high"
    )["value"] = "-1/2+0*sqrt5"
    with pytest.raises(independent.VerificationError, match="minimum evidence drift"):
        independent.verify_base_range_evidence(mutated)


@pytest.mark.parametrize(
    ("mutation", "expected_error"),
    [
        ("double_directed_labels_60", "geometry contract drift"),
        ("projective_axes_15", "geometry contract drift"),
        ("raw_norm2_differences", "geometry contract drift"),
        ("fz11_rank_six_sign", "exact moments drift"),
        ("mutated_M8", "exact moments drift"),
        ("q_above_one", "geometry contract drift"),
        ("unequal_weights", "geometry contract drift"),
        ("extra_hop_shell", "geometry contract drift"),
        ("sharp_R6", "Taylor remainder drift"),
        ("sharp_R8", "Taylor remainder drift"),
        ("physical_frequency", "physical boundary drift"),
        ("comparison_input", "exposure boundary drift"),
        ("extra_top_level", "top-level keys drift"),
    ],
)
def test_resigned_semantic_mutations_fail_closed(
    tmp_path: Path, mutation: str, expected_error: str
) -> None:
    def mutate(receipt: dict) -> None:
        if mutation == "double_directed_labels_60":
            receipt["geometry_contract"]["source_direction_count"] = 60
        elif mutation == "projective_axes_15":
            receipt["geometry_contract"]["source_direction_count"] = 15
        elif mutation == "raw_norm2_differences":
            receipt["geometry_contract"]["unit_direction_definition"] = (
                "w_e = carrierSeamDifference(e)"
            )
        elif mutation == "fz11_rank_six_sign":
            receipt["exact_moments"]["M6"] = "30/7 + (2/7) I6(n)"
        elif mutation == "mutated_M8":
            receipt["exact_moments"]["M8"] = "10/3 + (8/15) I6(n)"
        elif mutation == "q_above_one":
            receipt["geometry_contract"]["q_domain"]["maximum"] = "2"
        elif mutation == "unequal_weights":
            receipt["geometry_contract"]["equal_weight"] = "1/6"
        elif mutation == "extra_hop_shell":
            receipt["geometry_contract"]["extra_hop_shell"] = True
        elif mutation == "sharp_R6":
            receipt["taylor_remainders"]["R6"][
                "global_absolute_coefficient"
            ] = "1/100000"
        elif mutation == "sharp_R8":
            receipt["taylor_remainders"]["R8"]["global_coefficient"] = (
                "1/100000000"
            )
        elif mutation == "physical_frequency":
            receipt["physical_boundary"][
                "spatial_symbol_to_frequency_squared_proved"
            ] = True
        elif mutation == "comparison_input":
            receipt["exposure_boundary"]["comparison_inputs"] = ["public data"]
        else:
            receipt["unknown_promotion"] = True

    path = mutated_receipt(tmp_path, mutate)
    result = independent_run(receipt=path)
    assert result.returncode == 1
    assert expected_error in result.stderr


@pytest.mark.parametrize(
    "control",
    [
        "double_directed_labels_60",
        "projective_axes_15",
        "raw_norm2_differences",
        "fz11_rank_six_sign",
        "mutated_M8",
        "q_above_one",
        "unequal_weights",
        "extra_hop_shell",
    ],
)
def test_declared_control_cannot_be_disabled(tmp_path: Path, control: str) -> None:
    def mutate(receipt: dict) -> None:
        row = next(
            item
            for item in receipt["negative_controls"]["controls"]
            if item["control"] == control
        )
        row["detector_fired"] = False
        receipt["negative_controls"]["all_detectors_fired"] = False

    path = mutated_receipt(tmp_path, mutate)
    result = independent_run(receipt=path)
    assert result.returncode == 1
    assert "negative controls drift" in result.stderr
