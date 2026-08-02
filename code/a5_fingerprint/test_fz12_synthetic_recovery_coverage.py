from __future__ import annotations

import json
import os
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

import fz12_synthetic_recovery_coverage as producer
import verify_fz12_synthetic_recovery_coverage_independent as independent


VERIFIER = producer.HERE / "verify_fz12_synthetic_recovery_coverage_independent.py"


@pytest.fixture(scope="module")
def receipt() -> dict:
    return producer.build_receipt()


def write_canonical(path: Path, value: dict) -> None:
    path.write_bytes(independent.canonical_json_bytes(value))


def resign(value: dict, field: str) -> None:
    value[field] = independent.self_digest(value, field)


def verifier_run(
    *,
    source: Path = independent.DEFAULT_SOURCE,
    custody: Path = independent.DEFAULT_CUSTODY,
    remainder: Path = independent.DEFAULT_REMAINDER,
    receipt: Path = independent.DEFAULT_RECEIPT,
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [
            sys.executable,
            str(VERIFIER),
            "--source",
            str(source),
            "--custody",
            str(custody),
            "--remainder",
            str(remainder),
            "--receipt",
            str(receipt),
        ],
        cwd=producer.HERE,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def mutated_receipt(tmp_path: Path, mutation) -> Path:
    value = json.loads(independent.DEFAULT_RECEIPT.read_text())
    mutation(value)
    resign(value, "receipt_sha256")
    path = tmp_path / independent.DEFAULT_RECEIPT.name
    write_canonical(path, value)
    return path


def test_exact_exhaustive_recovery_and_coverage(receipt: dict) -> None:
    rows = receipt["coordinate_results"]
    assert [(row["coefficient"], row["coverage"]) for row in rows] == [
        ("C4_over_a2", "3904/4096"),
        ("B0_over_a4", "3892/4096"),
        ("B6_over_a4", "3892/4096"),
    ]
    assert [row["detection_power"] for row in rows] == [
        "4096/4096",
        "234/4096",
        "195/4096",
    ]
    assert receipt["joint_dimension_eight_result"]["coverage"] == "3892/4096"
    assert (
        receipt["joint_dimension_eight_result"]["detection_power"] == "209/4096"
    )
    assert receipt["recovery_conclusion"] == {
        "leading_coefficient_resolved": True,
        "leading_detection_power": "4096/4096",
        "linked_dimension_eight_terms_resolved": False,
        "linked_dimension_eight_joint_power": "209/4096",
        "interpretation": (
            "the exact synthetic regime rejects a zero leading coefficient "
            "in every error vector; its nominal 95 percent interval covers "
            "the injected value in 3904/4096 vectors, while the linked "
            "higher-order pair has power close to its calibrated "
            "false-positive rate"
        ),
    }


def test_design_is_full_rank_and_inside_certified_domain(receipt: dict) -> None:
    design = receipt["synthetic_design"]
    assert design["row_count"] == 12
    assert design["q_values"] == ["1/4", "1/2", "3/4", "1"]
    assert design["I6_values"] == ["1", "-5/9", "-5/16"]
    assert design["q_domain"] == "0 < q <= 1"
    assert design["full_column_rank"] is True
    study = independent.rebuild_study()
    assert design == independent.expected_design(study)


def test_exact_linked_signal_hierarchy(receipt: dict) -> None:
    hierarchy = receipt["exact_signal_hierarchy"]
    assert hierarchy == {
        "domain": "0 < q <= 1 and -5/9 <= I6 <= 1",
        "leading_absolute_term": "q^4/20",
        "full_linked_dimension_eight_absolute_ratio_to_leading": {
            "minimum": "q^2/45",
            "maximum": "2 q^2/81",
        },
        "anisotropic_dimension_eight_absolute_ratio_to_leading_upper": "q^2/630",
        "derivation": (
            "1/900 <= 1/840 - I6/12600 <= 1/810, divided by 1/20"
        ),
    }
    minimum = Fraction(1, 840) - Fraction(1, 12600)
    maximum = Fraction(1, 840) + Fraction(5, 9 * 12600)
    assert minimum / Fraction(1, 20) == Fraction(1, 45)
    assert maximum / Fraction(1, 20) == Fraction(2, 81)
    assert Fraction(1, 12600) / Fraction(1, 20) == Fraction(1, 630)


def test_dimensionless_sensitivity_grid_keeps_experimental_power_open(
    receipt: dict,
) -> None:
    grid = receipt["synthetic_sensitivity_grid"]
    assert grid["sigma_values"] == ["1/50", "1/100", "1/200", "1/500", "1/1000"]
    assert [row["leading_detection_power"] for row in grid["rows"]] == [
        "1127/4096",
        "3053/4096",
        "4096/4096",
        "4096/4096",
        "4096/4096",
    ]
    assert [
        row["joint_dimension_eight_detection_power"] for row in grid["rows"]
    ] == [
        "202/4096",
        "201/4096",
        "209/4096",
        "228/4096",
        "373/4096",
    ]
    assert "dataset-specific" in grid["scope"]
    assert "remain open" in grid["scope"]


def test_remainder_bias_cannot_change_synthetic_leading_result(receipt: dict) -> None:
    robustness = receipt["full_symbol_remainder_robustness"]
    assert robustness["source_bound"] == (
        "abs(lambda_hat - P6) <= (7/388800) q^8"
    )
    assert robustness["response_bound"] == (
        "for 0 < q <= 1, abs((lambda_hat - P6)/q^2) <= (7/388800) q^6"
    )
    fractions = [
        Fraction(row["fraction_of_interval_half_width"])
        for row in robustness["worst_case_ols_bias_bounds"]
    ]
    assert max(fractions) < Fraction(1, 300)
    assert robustness["minimum_exhaustive_leading_detection_margin"] == (
        "373/196900"
    )
    assert robustness["leading_remainder_bias_to_detection_margin"] == (
        "2370977/185628672"
    )
    assert robustness["worst_case_bias_less_than_detection_margin"] is True
    assert Fraction(robustness["leading_remainder_bias_to_detection_margin"]) < 1
    assert robustness["leading_recovery_changed"] is False
    assert robustness["experimental_remainder_model_supplied"] is False


def test_no_physical_or_experimental_promotion(receipt: dict) -> None:
    boundary = receipt["exposure_and_physical_boundary"]
    assert boundary["synthetic_inputs_only"] is True
    assert boundary["comparison_inputs"] == []
    assert all(
        value is False
        for key, value in boundary.items()
        if key not in {"synthetic_inputs_only", "comparison_inputs"}
    )
    assert receipt["study_scope"]["physical_interpretation"] is False
    assert receipt["study_scope"]["experimental_sensitivity_claim"] is False
    assert boundary["continuous_direction_coverage_calibrated"] is False
    assert boundary["orientation_profiled_over_SO3_mod_A5"] is False
    assert boundary["nuisance_model_coverage_calibrated"] is False
    assert boundary["gaussian_or_detector_response_calibrated"] is False
    assert boundary["five_sigma_tail_calibrated"] is False


def test_negative_controls_are_fail_closed(receipt: dict) -> None:
    controls = receipt["negative_controls"]
    assert controls["scope"] == (
        "contract and mutation guards, not producer-executed "
        "counterfactual coverage studies"
    )
    assert controls["all_guards_passed"] is True
    assert [row["control"] for row in controls["controls"]] == [
        "leading_sign_flip",
        "unlinked_isotropic_term",
        "rank_six_sign_flip",
        "rank_deficient_direction_design",
        "q_above_certified_domain",
        "partial_noise_enumeration",
        "public_comparison_input",
    ]
    assert all(row["guard_passed"] is True for row in controls["controls"])


def test_committed_receipt_is_canonical_and_byte_exact(receipt: dict) -> None:
    assert producer.verify_committed_receipt() == receipt
    assert independent.DEFAULT_RECEIPT.read_bytes() == (
        independent.canonical_json_bytes(receipt)
    )


def test_independent_verifier_is_separate_and_passes() -> None:
    source = VERIFIER.read_text(encoding="utf-8")
    assert "import fz12_synthetic_recovery_coverage" not in source
    assert "from fz12_synthetic_recovery_coverage" not in source
    result = verifier_run()
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "FZ12_SYNTHETIC_COVERAGE_VERIFY_PASS"


@pytest.mark.parametrize(
    ("parent", "self_field", "argument"),
    [
        ("seam_current_edge_prediction_receipt.json", "receipt_sha256", "source"),
        ("fz12_custody_projection.json", "projection_sha256", "custody"),
        ("fz12_full_symbol_remainder_receipt.json", "receipt_sha256", "remainder"),
    ],
)
def test_resigned_parent_mutation_fails_fixed_raw_pin(
    tmp_path: Path, parent: str, self_field: str, argument: str
) -> None:
    original = independent.RUNTIME / parent
    value = json.loads(original.read_text())
    value["audit_mutation"] = True
    resign(value, self_field)
    path = tmp_path / parent
    write_canonical(path, value)
    result = verifier_run(**{argument: path})
    assert result.returncode == 1
    assert "parent raw hash drift" in result.stderr


@pytest.mark.parametrize(
    ("mutation", "expected_error"),
    [
        ("schema", "synthetic schema drift"),
        ("coefficient", "synthetic injected ray drift"),
        ("hierarchy", "exact signal hierarchy drift"),
        ("q_design", "design drift"),
        ("noise", "noise calibration drift"),
        ("sensitivity_grid", "synthetic sensitivity grid drift"),
        ("coverage", "coordinate coverage or power drift"),
        ("joint_power", "joint dimension-eight result drift"),
        ("conclusion", "synthetic recovery conclusion drift"),
        ("remainder_bias", "full-symbol remainder robustness drift"),
        ("public_measurement", "exposure or physical boundary drift"),
        ("physical_sector", "exposure or physical boundary drift"),
        ("control", "negative control guard disabled"),
        ("extra_key", "synthetic top-level keys drift"),
    ],
)
def test_resigned_semantic_mutations_fail_closed(
    tmp_path: Path, mutation: str, expected_error: str
) -> None:
    def mutate(value: dict) -> None:
        if mutation == "schema":
            value["schema"] = "oph.fz12.synthetic_recovery_coverage.v2"
        elif mutation == "coefficient":
            value["injected_frozen_ray"]["C4_over_a2"] = "1/20"
        elif mutation == "hierarchy":
            value["exact_signal_hierarchy"][
                "anisotropic_dimension_eight_absolute_ratio_to_leading_upper"
            ] = "q^2/600"
        elif mutation == "q_design":
            value["synthetic_design"]["q_values"][-1] = "5/4"
        elif mutation == "noise":
            value["noise_and_calibration"]["sigma"] = "1/100"
        elif mutation == "sensitivity_grid":
            value["synthetic_sensitivity_grid"]["rows"][0][
                "joint_dimension_eight_detection_power"
            ] = "4096/4096"
        elif mutation == "coverage":
            value["coordinate_results"][0]["covered_replicates"] += 1
        elif mutation == "joint_power":
            value["joint_dimension_eight_result"]["signal_detections"] += 1
        elif mutation == "conclusion":
            value["recovery_conclusion"]["linked_dimension_eight_terms_resolved"] = True
        elif mutation == "remainder_bias":
            value["full_symbol_remainder_robustness"][
                "leading_recovery_changed"
            ] = True
        elif mutation == "public_measurement":
            value["exposure_and_physical_boundary"]["public_measurement_read"] = True
        elif mutation == "physical_sector":
            value["exposure_and_physical_boundary"]["physical_sector_identified"] = True
        elif mutation == "control":
            value["negative_controls"]["controls"][0]["guard_passed"] = False
        else:
            value["unknown_promotion"] = True

    path = mutated_receipt(tmp_path, mutate)
    result = verifier_run(receipt=path)
    assert result.returncode == 1
    assert expected_error in result.stderr


def test_noncanonical_receipt_is_rejected(tmp_path: Path) -> None:
    value = json.loads(independent.DEFAULT_RECEIPT.read_text())
    path = tmp_path / independent.DEFAULT_RECEIPT.name
    path.write_text(json.dumps(value, indent=2), encoding="ascii")
    result = verifier_run(receipt=path)
    assert result.returncode == 1
    assert "noncanonical synthetic receipt" in result.stderr


def test_independent_study_exact_counts() -> None:
    study = independent.rebuild_study()
    assert [row["covered_replicates"] for row in study["coordinate_rows"]] == [
        3904,
        3892,
        3892,
    ]
    assert [row["signal_detections"] for row in study["coordinate_rows"]] == [
        4096,
        234,
        195,
    ]
    assert study["joint_covered"] == 3892
    assert study["joint_detections"] == 209


def test_declared_controls_cannot_be_deleted(tmp_path: Path) -> None:
    def mutate(value: dict) -> None:
        value["negative_controls"]["controls"].pop()

    path = mutated_receipt(tmp_path, mutate)
    result = verifier_run(receipt=path)
    assert result.returncode == 1
    assert "negative control list drift" in result.stderr
