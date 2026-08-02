#!/usr/bin/env python3
"""Independently verify the target-free FZ-12 synthetic coverage packet.

This verifier does not import the producer.  It rebuilds the rational design,
the ordinary-least-squares map, all 4096 independent Rademacher error vectors,
the coordinate and joint 95 percent calibrations, the injected-signal power,
and the certified-remainder bias bounds from separate code.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
DEFAULT_SOURCE = RUNTIME / "seam_current_edge_prediction_receipt.json"
DEFAULT_CUSTODY = RUNTIME / "fz12_custody_projection.json"
DEFAULT_REMAINDER = RUNTIME / "fz12_full_symbol_remainder_receipt.json"
DEFAULT_RECEIPT = RUNTIME / "fz12_synthetic_recovery_coverage_receipt.json"

SCHEMA = "oph.fz12.synthetic_recovery_coverage.v1"
STATUS = (
    "EXACT_EXHAUSTIVE_SYNTHETIC_LEADING_RECOVERY__"
    "LINKED_DIMENSION_EIGHT_UNRESOLVED__PHYSICAL_COMPARISON_UNARMED"
)

PARENT_SPECS = (
    {
        "relative_path": (
            "code/a5_fingerprint/runtime/seam_current_edge_prediction_receipt.json"
        ),
        "bytes": 9296,
        "sha256": (
            "0b8f0f7573f556ef0f47158fe07eca002c5b35790231d1ef2b75518057d12915"
        ),
        "schema": "oph.seam_current_edge_prediction_candidate.v1",
        "status": (
            "EXACT_SOURCE_NATIVE_EDGE_RAY__PROSPECTIVE_PHYSICAL_BRANCH_UNARMED__"
            "PHYSICAL_PRODUCER_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "frozen FZ-12 normalized spatial coefficient ray",
    },
    {
        "relative_path": "code/a5_fingerprint/runtime/fz12_custody_projection.json",
        "bytes": 3624,
        "sha256": (
            "dfc6574a5bde9b576df6aeef2e03aba8602ec2723a0a4606e7942c405a57c643"
        ),
        "schema": "oph.fz12.custody_projection.v1",
        "status": "FROZEN_FZ12_SOURCE_AND_CUSTODY_PINNED__NO_COMPARISON_DATA",
        "self_field": "projection_sha256",
        "role": "data-free FZ-12 freeze and custody projection",
    },
    {
        "relative_path": (
            "code/a5_fingerprint/runtime/fz12_full_symbol_remainder_receipt.json"
        ),
        "bytes": 5367,
        "sha256": (
            "ce4052f633c891515ebc319e9d7ff2bc0044bb8bb9e76a890e8e22c8e882dcfa"
        ),
        "schema": "oph.fz12.full_symbol_remainder.v1",
        "status": (
            "EXACT_TARGET_FREE_EDGE_SYMBOL_Q_LE_ONE_REMAINDER__"
            "PHYSICAL_FREQUENCY_AND_COMPARISON_OPEN"
        ),
        "self_field": "receipt_sha256",
        "role": "target-free q at most one spatial-symbol remainder contract",
    },
)

Q_VALUES = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1))
I6_VALUES = (Fraction(1), Fraction(-5, 9), Fraction(-5, 16))
NOISE = Fraction(1, 200)
NOISE_GRID = (
    Fraction(1, 50),
    Fraction(1, 100),
    Fraction(1, 200),
    Fraction(1, 500),
    Fraction(1, 1000),
)
COVERAGE = Fraction(19, 20)
BETAS = (Fraction(-1, 20), Fraction(1, 840), Fraction(-1, 12600))
NAMES = ("C4_over_a2", "B0_over_a4", "B6_over_a4")

Matrix = list[list[Fraction]]


class VerificationError(ValueError):
    """A synthetic coverage contract check failed."""


def check(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def self_digest(value: dict[str, Any], field: str) -> str:
    body = {key: item for key, item in value.items() if key != field}
    return "sha256:" + sha256(canonical_json_bytes(body))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def matrix_text(value: Matrix) -> list[list[str]]:
    return [[fraction_text(item) for item in row] for row in value]


def load_parent(
    path: Path, spec: dict[str, Any]
) -> tuple[bytes, dict[str, Any]]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError(f"invalid parent JSON: {path.name}") from error
    check(isinstance(value, dict), f"parent root is not an object: {path.name}")
    check(raw == canonical_json_bytes(value), f"noncanonical parent: {path.name}")
    check(sha256(raw) == spec["sha256"], f"parent raw hash drift: {path.name}")
    check(len(raw) == spec["bytes"], f"parent byte-count drift: {path.name}")
    check(value.get("schema") == spec["schema"], f"parent schema drift: {path.name}")
    check(value.get("status") == spec["status"], f"parent status drift: {path.name}")
    check(
        value.get(spec["self_field"]) == self_digest(value, spec["self_field"]),
        f"parent self-digest drift: {path.name}",
    )
    return raw, value


def transposed(value: Matrix) -> Matrix:
    return [list(column) for column in zip(*value, strict=True)]


def product(left: Matrix, right: Matrix) -> Matrix:
    check(left and right and len(left[0]) == len(right), "matrix dimension mismatch")
    right_t = transposed(right)
    return [
        [sum(x * y for x, y in zip(row, column, strict=True)) for column in right_t]
        for row in left
    ]


def inverse(value: Matrix) -> Matrix:
    size = len(value)
    check(size > 0 and all(len(row) == size for row in value), "nonsquare matrix")
    work = [
        list(row) + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(value)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        check(pivot is not None, "singular matrix")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        work[column] = [entry / divisor for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column], strict=True)
            ]
    return [row[size:] for row in work]


def design() -> Matrix:
    return [
        [q * q, q**4, i6 * q**4]
        for q in Q_VALUES
        for i6 in I6_VALUES
    ]


def rebuild_study() -> dict[str, Any]:
    x = design()
    xt = transposed(x)
    gram = product(xt, x)
    inverse_gram = inverse(gram)
    estimator = product(inverse_gram, xt)
    check(
        product(estimator, x)
        == [[Fraction(int(i == j)) for j in range(3)] for i in range(3)],
        "independent OLS left inverse failed",
    )

    def errors_at(noise: Fraction) -> list[tuple[Fraction, ...]]:
        result = []
        for signs in itertools.product((-1, 1), repeat=12):
            result.append(
                tuple(
                    sum(
                        weight * noise * sign
                        for weight, sign in zip(row, signs, strict=True)
                    )
                    for row in estimator
                )
            )
        check(len(result) == 4096, "independent error universe size drift")
        check(
            all(sum(error[index] for error in result) == 0 for index in range(3)),
            "independent error universe is not centered",
        )
        return result

    errors = errors_at(NOISE)

    def threshold(values: list[Fraction]) -> Fraction:
        ordered = sorted(values)
        rank = (COVERAGE.numerator * len(ordered) + COVERAGE.denominator - 1)
        rank //= COVERAGE.denominator
        return ordered[rank - 1]

    coordinate_rows = []
    thresholds = []
    for index, (name, beta) in enumerate(zip(NAMES, BETAS, strict=True)):
        absolute_errors = [abs(error[index]) for error in errors]
        cutoff = threshold(absolute_errors)
        thresholds.append(cutoff)
        covered = sum(error <= cutoff for error in absolute_errors)
        detections = sum(abs(beta + error[index]) > cutoff for error in errors)
        variance = NOISE**2 * inverse_gram[index][index]
        coordinate_rows.append(
            {
                "coefficient": name,
                "injected_value": fraction_text(beta),
                "interval_half_width": fraction_text(cutoff),
                "covered_replicates": covered,
                "coverage": f"{covered}/4096",
                "signal_detections": detections,
                "detection_power": f"{detections}/4096",
                "estimator_variance": fraction_text(variance),
                "squared_signal_to_noise": fraction_text(beta**2 / variance),
            }
        )

    covariance = [
        [NOISE**2 * inverse_gram[i][j] for j in (1, 2)] for i in (1, 2)
    ]
    precision = inverse(covariance)

    def statistic(vector: tuple[Fraction, Fraction]) -> Fraction:
        return sum(
            vector[i] * precision[i][j] * vector[j]
            for i in range(2)
            for j in range(2)
        )

    null_statistics = [statistic((error[1], error[2])) for error in errors]
    joint_cutoff = threshold(null_statistics)
    joint_covered = sum(value <= joint_cutoff for value in null_statistics)
    joint_detections = sum(
        statistic((BETAS[1] + error[1], BETAS[2] + error[2])) > joint_cutoff
        for error in errors
    )

    q_by_row = [q for q in Q_VALUES for _ in I6_VALUES]
    remainder_bias = [
        sum(
            abs(estimator[index][row])
            * Fraction(7, 388800)
            * q_by_row[row] ** 6
            for row in range(12)
        )
        for index in range(3)
    ]
    leading_detection_margin = min(
        abs(BETAS[0] + error[0]) - thresholds[0] for error in errors
    )
    check(leading_detection_margin > 0, "leading recovery has no strict margin")
    check(
        remainder_bias[0] < leading_detection_margin,
        "certified remainder can erase the leading synthetic detection",
    )

    grid_rows = []
    for noise in NOISE_GRID:
        grid_errors = errors_at(noise)
        grid_thresholds = [
            threshold([abs(error[index]) for error in grid_errors])
            for index in range(3)
        ]
        grid_detections = [
            sum(
                abs(BETAS[index] + error[index]) > grid_thresholds[index]
                for error in grid_errors
            )
            for index in range(3)
        ]
        grid_covariance = [
            [noise**2 * inverse_gram[i][j] for j in (1, 2)] for i in (1, 2)
        ]
        grid_precision = inverse(grid_covariance)

        def grid_statistic(vector: tuple[Fraction, Fraction]) -> Fraction:
            return sum(
                vector[i] * grid_precision[i][j] * vector[j]
                for i in range(2)
                for j in range(2)
            )

        grid_null = [grid_statistic((error[1], error[2])) for error in grid_errors]
        grid_joint_cutoff = threshold(grid_null)
        grid_joint_detections = sum(
            grid_statistic((BETAS[1] + error[1], BETAS[2] + error[2]))
            > grid_joint_cutoff
            for error in grid_errors
        )
        grid_rows.append(
            {
                "sigma": fraction_text(noise),
                "maximum_row_leading_signal_over_sigma": fraction_text(
                    Fraction(1, 20) / noise
                ),
                "maximum_row_full_dimension_eight_signal_over_sigma": fraction_text(
                    Fraction(1, 810) / noise
                ),
                "maximum_row_anisotropic_signal_over_sigma": fraction_text(
                    Fraction(1, 12600) / noise
                ),
                "leading_detection_power": f"{grid_detections[0]}/4096",
                "B0_detection_power": f"{grid_detections[1]}/4096",
                "B6_detection_power": f"{grid_detections[2]}/4096",
                "joint_dimension_eight_detection_power": (
                    f"{grid_joint_detections}/4096"
                ),
            }
        )
    return {
        "design": x,
        "gram": gram,
        "inverse_gram": inverse_gram,
        "coordinate_rows": coordinate_rows,
        "thresholds": thresholds,
        "covariance": covariance,
        "joint_cutoff": joint_cutoff,
        "joint_covered": joint_covered,
        "joint_detections": joint_detections,
        "remainder_bias": remainder_bias,
        "leading_detection_margin": leading_detection_margin,
        "sensitivity_grid": grid_rows,
    }


def expected_parent_pins(
    loaded: list[tuple[bytes, dict[str, Any]]],
) -> list[dict[str, Any]]:
    return [
        {
            "path": spec["relative_path"],
            "role": spec["role"],
            "bytes": len(raw),
            "sha256": "sha256:" + sha256(raw),
            "schema": spec["schema"],
            "status": spec["status"],
            "self_digest": value[spec["self_field"]],
        }
        for spec, (raw, value) in zip(PARENT_SPECS, loaded, strict=True)
    ]


def expected_design(study: dict[str, Any]) -> dict[str, Any]:
    rows = []
    for index, (q, i6, features) in enumerate(
        zip(
            (q for q in Q_VALUES for _ in I6_VALUES),
            (i6 for _ in Q_VALUES for i6 in I6_VALUES),
            study["design"],
            strict=True,
        )
    ):
        rows.append(
            {
                "row": index,
                "q": fraction_text(q),
                "I6": fraction_text(i6),
                "features": [fraction_text(item) for item in features],
            }
        )
    return {
        "row_count": 12,
        "q_values": ["1/4", "1/2", "3/4", "1"],
        "I6_values": ["1", "-5/9", "-5/16"],
        "q_domain": "0 < q <= 1",
        "rows": rows,
        "feature_order": list(NAMES),
        "gram_matrix": matrix_text(study["gram"]),
        "inverse_gram_matrix": matrix_text(study["inverse_gram"]),
        "full_column_rank": True,
    }


def expected_remainder(study: dict[str, Any]) -> dict[str, Any]:
    rows = []
    for name, bound, cutoff in zip(
        NAMES,
        study["remainder_bias"],
        study["thresholds"],
        strict=True,
    ):
        rows.append(
            {
                "coefficient": name,
                "worst_case_bias_bound": fraction_text(bound),
                "fraction_of_interval_half_width": fraction_text(bound / cutoff),
            }
        )
    return {
        "source_bound": "abs(lambda_hat - P6) <= (7/388800) q^8",
        "response_bound": (
            "for 0 < q <= 1, abs((lambda_hat - P6)/q^2) "
            "<= (7/388800) q^6"
        ),
        "worst_case_ols_bias_bounds": rows,
        "minimum_exhaustive_leading_detection_margin": fraction_text(
            study["leading_detection_margin"]
        ),
        "leading_remainder_bias_to_detection_margin": fraction_text(
            study["remainder_bias"][0] / study["leading_detection_margin"]
        ),
        "worst_case_bias_less_than_detection_margin": True,
        "leading_recovery_changed": False,
        "experimental_remainder_model_supplied": False,
    }


def verify_receipt(
    path: Path,
    loaded: list[tuple[bytes, dict[str, Any]]],
) -> dict[str, Any]:
    raw = path.read_bytes()
    try:
        receipt = json.loads(raw)
    except json.JSONDecodeError as error:
        raise VerificationError("invalid synthetic receipt JSON") from error
    check(isinstance(receipt, dict), "synthetic receipt root is not an object")
    check(raw == canonical_json_bytes(receipt), "noncanonical synthetic receipt")
    check(receipt.get("schema") == SCHEMA, "synthetic schema drift")
    check(receipt.get("status") == STATUS, "synthetic status drift")
    check(receipt.get("issue") == 667, "synthetic issue drift")
    check(
        receipt.get("receipt_sha256") == self_digest(receipt, "receipt_sha256"),
        "synthetic self-digest drift",
    )
    expected_keys = {
        "schema",
        "status",
        "issue",
        "parent_pins",
        "study_scope",
        "injected_frozen_ray",
        "exact_signal_hierarchy",
        "synthetic_design",
        "noise_and_calibration",
        "synthetic_sensitivity_grid",
        "coordinate_results",
        "joint_dimension_eight_result",
        "recovery_conclusion",
        "full_symbol_remainder_robustness",
        "exposure_and_physical_boundary",
        "negative_controls",
        "receipt_sha256",
    }
    check(set(receipt) == expected_keys, "synthetic top-level keys drift")
    check(
        receipt.get("parent_pins") == expected_parent_pins(loaded),
        "synthetic parent pins drift",
    )

    source, custody, remainder = (value for _, value in loaded)
    check(
        source.get("conditional_physical_candidate", {}).get("coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "frozen source ray drift",
    )
    check(
        custody.get("projection_scope", {}).get("includes_measurement_values") is False,
        "custody contains measurement values",
    )
    check(
        remainder.get("taylor_remainders", {}).get("R6", {}).get(
            "global_absolute_coefficient"
        )
        == "7/388800",
        "remainder coefficient drift",
    )

    # Cheap semantic preflight keeps adversarial mutation tests bounded.  A
    # valid receipt is subsequently recomputed from all 4096 sign vectors.
    design_preflight = receipt.get("synthetic_design", {})
    check(
        design_preflight.get("q_values") == ["1/4", "1/2", "3/4", "1"]
        and design_preflight.get("I6_values") == ["1", "-5/9", "-5/16"]
        and design_preflight.get("row_count") == 12
        and design_preflight.get("full_column_rank") is True,
        "design drift",
    )
    grid_preflight = receipt.get("synthetic_sensitivity_grid", {})
    grid_rows = grid_preflight.get("rows", [])
    check(
        grid_preflight.get("sigma_values")
        == ["1/50", "1/100", "1/200", "1/500", "1/1000"]
        and len(grid_rows) == 5
        and [row.get("joint_dimension_eight_detection_power") for row in grid_rows]
        == ["202/4096", "201/4096", "209/4096", "228/4096", "373/4096"],
        "synthetic sensitivity grid drift",
    )
    coordinate_preflight = receipt.get("coordinate_results", [])
    check(
        len(coordinate_preflight) == 3
        and [row.get("covered_replicates") for row in coordinate_preflight]
        == [3904, 3892, 3892]
        and [row.get("signal_detections") for row in coordinate_preflight]
        == [4096, 234, 195],
        "coordinate coverage or power drift",
    )
    joint_preflight = receipt.get("joint_dimension_eight_result", {})
    check(
        joint_preflight.get("covered_replicates") == 3892
        and joint_preflight.get("signal_detections") == 209
        and joint_preflight.get("resolved_in_declared_synthetic_regime") is False,
        "joint dimension-eight result drift",
    )
    remainder_preflight = receipt.get("full_symbol_remainder_robustness", {})
    check(
        remainder_preflight.get("leading_recovery_changed") is False
        and remainder_preflight.get("worst_case_bias_less_than_detection_margin")
        is True
        and remainder_preflight.get("experimental_remainder_model_supplied") is False
        and len(remainder_preflight.get("worst_case_ols_bias_bounds", [])) == 3,
        "full-symbol remainder robustness drift",
    )
    check(
        receipt.get("injected_frozen_ray", {}).get("C4_over_a2") == "-1/20"
        and receipt.get("injected_frozen_ray", {}).get("B0_over_a4") == "1/840"
        and receipt.get("injected_frozen_ray", {}).get("B6_over_a4")
        == "-1/12600",
        "synthetic injected ray drift",
    )
    check(
        receipt.get("exact_signal_hierarchy", {}).get(
            "anisotropic_dimension_eight_absolute_ratio_to_leading_upper"
        )
        == "q^2/630",
        "exact signal hierarchy drift",
    )
    check(
        receipt.get("noise_and_calibration", {}).get("sigma") == "1/200",
        "noise calibration drift",
    )
    check(
        receipt.get("recovery_conclusion", {}).get(
            "linked_dimension_eight_terms_resolved"
        )
        is False,
        "synthetic recovery conclusion drift",
    )
    check(
        receipt.get("exposure_and_physical_boundary", {}).get(
            "public_measurement_read"
        )
        is False
        and receipt.get("exposure_and_physical_boundary", {}).get(
            "physical_sector_identified"
        )
        is False,
        "exposure or physical boundary drift",
    )
    check(
        all(
            row.get("guard_passed") is True
            for row in receipt.get("negative_controls", {}).get("controls", [])
        ),
        "negative control guard disabled",
    )

    study = rebuild_study()
    check(
        receipt.get("study_scope")
        == {
            "type": "exact exhaustive target-free synthetic injection and recovery",
            "response": "y = (lambda_hat - q^2) / q^2",
            "injection_model": (
                "y = C4_over_a2 q^2 + (B0_over_a4 + "
                "B6_over_a4 I6) q^4"
            ),
            "fit_model": "unconstrained three-coefficient ordinary least squares",
            "dimension_eight_label": (
                "formal EFT derivative order only; no physical operator is identified"
            ),
            "physical_interpretation": False,
            "experimental_sensitivity_claim": False,
        },
        "synthetic study scope drift",
    )
    check(
        receipt.get("injected_frozen_ray")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
            "B0_over_C4_squared": "10/21",
            "B6_over_C4_squared": "-2/63",
            "B6_over_B0": "-1/15",
        },
        "synthetic injected ray drift",
    )
    check(
        receipt.get("exact_signal_hierarchy")
        == {
            "domain": "0 < q <= 1 and -5/9 <= I6 <= 1",
            "leading_absolute_term": "q^4/20",
            "full_linked_dimension_eight_absolute_ratio_to_leading": {
                "minimum": "q^2/45",
                "maximum": "2 q^2/81",
            },
            "anisotropic_dimension_eight_absolute_ratio_to_leading_upper": (
                "q^2/630"
            ),
            "derivation": (
                "1/900 <= 1/840 - I6/12600 <= 1/810, divided by 1/20"
            ),
        },
        "exact signal hierarchy drift",
    )
    check(receipt.get("synthetic_design") == expected_design(study), "design drift")
    check(
        receipt.get("noise_and_calibration")
        == {
            "law": "independent Rademacher errors, each exactly plus or minus sigma",
            "sigma": "1/200",
            "enumeration": "all 2^12 sign vectors",
            "replicate_count": 4096,
            "nominal_coordinate_coverage": "19/20",
            "interval_rule": (
                "beta_hat plus or minus the smallest exhaustive absolute-error "
                "quantile attaining at least 19/20 coverage"
            ),
            "joint_dimension_eight_rule": (
                "exact covariance quadratic form with its exhaustive 19/20 null quantile"
            ),
            "noise_scope": (
                "dimensionless synthetic stress law, not an experimental response model"
            ),
            "sampling_error": "0",
        },
        "noise calibration drift",
    )
    check(
        receipt.get("synthetic_sensitivity_grid")
        == {
            "q_values": ["1/4", "1/2", "3/4", "1"],
            "sigma_values": ["1/50", "1/100", "1/200", "1/500", "1/1000"],
            "rows": study["sensitivity_grid"],
            "scope": (
                "dimensionless synthetic grid only; dataset-specific response, "
                "covariance, and power remain open"
            ),
            "grid_conclusion": (
                "the leading term reaches full exhaustive power by sigma 1/200, "
                "while joint higher-order power stays at or below 373/4096 on "
                "the declared grid"
            ),
        },
        "synthetic sensitivity grid drift",
    )
    check(
        receipt.get("coordinate_results") == study["coordinate_rows"],
        "coordinate coverage or power drift",
    )
    check(
        receipt.get("joint_dimension_eight_result")
        == {
            "coefficients": ["B0_over_a4", "B6_over_a4"],
            "covariance": matrix_text(study["covariance"]),
            "null_statistic_threshold": fraction_text(study["joint_cutoff"]),
            "covered_replicates": study["joint_covered"],
            "coverage": f"{study['joint_covered']}/4096",
            "signal_detections": study["joint_detections"],
            "detection_power": f"{study['joint_detections']}/4096",
            "resolved_in_declared_synthetic_regime": False,
        },
        "joint dimension-eight result drift",
    )
    check(
        receipt.get("recovery_conclusion")
        == {
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
        },
        "synthetic recovery conclusion drift",
    )
    check(
        receipt.get("full_symbol_remainder_robustness") == expected_remainder(study),
        "full-symbol remainder robustness drift",
    )
    boundary = receipt.get("exposure_and_physical_boundary")
    check(
        boundary
        == {
            "synthetic_inputs_only": True,
            "comparison_inputs": [],
            "public_measurement_read": False,
            "target_measurement_read": False,
            "comparison_permitted": False,
            "comparison_budget_consumed": False,
            "score_emitted": False,
            "evidence_claimed": False,
            "verdict_emitted": False,
            "physical_position_identified": False,
            "physical_sector_identified": False,
            "photon_attachment_proved": False,
            "clock_frame_and_boost_proved": False,
            "wave_packet_and_readout_proved": False,
            "present_experimental_sensitivity_established": False,
            "continuous_direction_coverage_calibrated": False,
            "orientation_profiled_over_SO3_mod_A5": False,
            "nuisance_model_coverage_calibrated": False,
            "gaussian_or_detector_response_calibrated": False,
            "five_sigma_tail_calibrated": False,
        },
        "exposure or physical boundary drift",
    )
    controls = receipt.get("negative_controls")
    check(isinstance(controls, dict), "negative controls missing")
    check(
        controls.get("scope")
        == (
            "contract and mutation guards, not producer-executed "
            "counterfactual coverage studies"
        ),
        "negative control scope drift",
    )
    check(controls.get("all_guards_passed") is True, "negative controls disabled")
    check(
        [row.get("control") for row in controls.get("controls", [])]
        == [
            "leading_sign_flip",
            "unlinked_isotropic_term",
            "rank_six_sign_flip",
            "rank_deficient_direction_design",
            "q_above_certified_domain",
            "partial_noise_enumeration",
            "public_comparison_input",
        ],
        "negative control list drift",
    )
    check(
        all(row.get("guard_passed") is True for row in controls["controls"]),
        "negative control guard disabled",
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--custody", type=Path, default=DEFAULT_CUSTODY)
    parser.add_argument("--remainder", type=Path, default=DEFAULT_REMAINDER)
    parser.add_argument("--receipt", type=Path, default=DEFAULT_RECEIPT)
    args = parser.parse_args(argv)
    try:
        loaded = [
            load_parent(path, spec)
            for path, spec in zip(
                (args.source, args.custody, args.remainder),
                PARENT_SPECS,
                strict=True,
            )
        ]
        verify_receipt(args.receipt, loaded)
    except (OSError, VerificationError) as error:
        print(f"FZ12_SYNTHETIC_COVERAGE_VERIFY_FAIL: {error}", file=sys.stderr)
        return 1
    print("FZ12_SYNTHETIC_COVERAGE_VERIFY_PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
