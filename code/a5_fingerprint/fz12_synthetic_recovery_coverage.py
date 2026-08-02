#!/usr/bin/env python3
"""Build the target-free FZ-12 synthetic recovery and coverage packet.

The packet studies the frozen normalized spatial coefficient ray in a small,
fully enumerated synthetic experiment.  Twelve prescribed design rows are
crossed with every one of the ``2^12`` independent Rademacher error vectors.
All regression, interval, coverage, and power arithmetic is exact over the
rationals.  The study therefore has no random seed, sampling error, or
platform-dependent numerical output.

The exercise rejects a zero leading normalized coefficient for every error
vector in its declared noise regime.  Its nominal 95 percent intervals cover
the injected leading value in 3904 of 4096 vectors.  The two linked sixth-order
spatial coefficients (the formal dimension-eight EFT terms) have essentially
calibration-level power.  The exercise reads no public measurement and
supplies no physical attachment, experimental sensitivity statement,
evidential score, or OPH verdict.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
SOURCE_PATH = RUNTIME / "seam_current_edge_prediction_receipt.json"
CUSTODY_PATH = RUNTIME / "fz12_custody_projection.json"
REMAINDER_PATH = RUNTIME / "fz12_full_symbol_remainder_receipt.json"
RECEIPT_PATH = RUNTIME / "fz12_synthetic_recovery_coverage_receipt.json"

SCHEMA = "oph.fz12.synthetic_recovery_coverage.v1"
STATUS = (
    "EXACT_EXHAUSTIVE_SYNTHETIC_LEADING_RECOVERY__"
    "LINKED_DIMENSION_EIGHT_UNRESOLVED__PHYSICAL_COMPARISON_UNARMED"
)

PARENT_CONTRACTS = (
    {
        "path": SOURCE_PATH,
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
        "path": CUSTODY_PATH,
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
        "path": REMAINDER_PATH,
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
NOISE_AMPLITUDE = Fraction(1, 200)
SENSITIVITY_NOISE_GRID = (
    Fraction(1, 50),
    Fraction(1, 100),
    Fraction(1, 200),
    Fraction(1, 500),
    Fraction(1, 1000),
)
NOMINAL_COVERAGE = Fraction(19, 20)
COEFFICIENTS = (
    Fraction(-1, 20),
    Fraction(1, 840),
    Fraction(-1, 12600),
)
COEFFICIENT_NAMES = ("C4_over_a2", "B0_over_a4", "B6_over_a4")

Matrix = list[list[Fraction]]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def sha256_hex(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def matrix_text(value: Matrix) -> list[list[str]]:
    return [[fraction_text(item) for item in row] for row in value]


def load_parent(contract: dict[str, Any]) -> tuple[bytes, dict[str, Any]]:
    raw = contract["path"].read_bytes()
    try:
        value = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(
            f"invalid parent JSON: {contract['path'].name}"
        ) from error
    require(isinstance(value, dict), f"parent is not an object: {contract['path'].name}")
    require(
        raw == base.canonical_json_bytes(value),
        f"noncanonical parent: {contract['path'].name}",
    )
    require(len(raw) == contract["bytes"], f"parent byte drift: {contract['path'].name}")
    require(
        sha256_hex(raw) == contract["sha256"],
        f"parent hash drift: {contract['path'].name}",
    )
    require(
        value.get("schema") == contract["schema"],
        f"parent schema drift: {contract['path'].name}",
    )
    require(
        value.get("status") == contract["status"],
        f"parent status drift: {contract['path'].name}",
    )
    self_field = contract["self_field"]
    body = {key: item for key, item in value.items() if key != self_field}
    require(
        value.get(self_field) == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"parent self-digest drift: {contract['path'].name}",
    )
    return raw, value


def transpose(value: Matrix) -> Matrix:
    return [list(column) for column in zip(*value, strict=True)]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    require(len(left[0]) == len(right), "matrix dimension mismatch")
    right_t = transpose(right)
    return [
        [sum(a * b for a, b in zip(row, column, strict=True)) for column in right_t]
        for row in left
    ]


def invert(value: Matrix) -> Matrix:
    size = len(value)
    require(size > 0 and all(len(row) == size for row in value), "nonsquare matrix")
    augmented = [
        list(row) + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(value)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        require(pivot is not None, "singular design matrix")
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        divisor = augmented[column][column]
        augmented[column] = [item / divisor for item in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                item - factor * pivot_item
                for item, pivot_item in zip(
                    augmented[row], augmented[column], strict=True
                )
            ]
    return [row[size:] for row in augmented]


def design_matrix() -> Matrix:
    return [
        [q * q, q**4, i6 * q**4]
        for q in Q_VALUES
        for i6 in I6_VALUES
    ]


def estimator_matrix(design: Matrix) -> tuple[Matrix, Matrix, Matrix]:
    design_t = transpose(design)
    gram = multiply(design_t, design)
    inverse_gram = invert(gram)
    estimator = multiply(inverse_gram, design_t)
    require(multiply(estimator, design) == identity(3), "OLS left inverse failed")
    return gram, inverse_gram, estimator


def identity(size: int) -> Matrix:
    return [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def error_universe(
    estimator: Matrix, noise: Fraction = NOISE_AMPLITUDE
) -> list[tuple[Fraction, ...]]:
    errors: list[tuple[Fraction, ...]] = []
    for signs in itertools.product((-1, 1), repeat=len(estimator[0])):
        errors.append(
            tuple(
                sum(
                    weight * noise * sign
                    for weight, sign in zip(row, signs, strict=True)
                )
                for row in estimator
            )
        )
    require(len(errors) == 4096, "synthetic noise universe size drift")
    require(
        all(sum(error[index] for error in errors) == 0 for index in range(3)),
        "synthetic error universe lost exact centering",
    )
    return errors


def quantile_threshold(values: list[Fraction]) -> Fraction:
    ordered = sorted(values)
    count = len(ordered)
    rank = (NOMINAL_COVERAGE.numerator * count + NOMINAL_COVERAGE.denominator - 1)
    rank //= NOMINAL_COVERAGE.denominator
    require(1 <= rank <= count, "invalid coverage rank")
    return ordered[rank - 1]


def quadratic_form(vector: tuple[Fraction, ...], precision: Matrix) -> Fraction:
    return sum(
        vector[i] * precision[i][j] * vector[j]
        for i in range(len(vector))
        for j in range(len(vector))
    )


def sensitivity_grid(estimator: Matrix, inverse_gram: Matrix) -> list[dict[str, Any]]:
    rows = []
    for noise in SENSITIVITY_NOISE_GRID:
        errors = error_universe(estimator, noise)
        thresholds = [
            quantile_threshold([abs(error[index]) for error in errors])
            for index in range(3)
        ]
        detections = [
            sum(abs(COEFFICIENTS[index] + error[index]) > thresholds[index] for error in errors)
            for index in range(3)
        ]
        covariance = [
            [noise**2 * inverse_gram[i][j] for j in (1, 2)] for i in (1, 2)
        ]
        precision = invert(covariance)
        null_statistics = [
            quadratic_form((error[1], error[2]), precision) for error in errors
        ]
        joint_threshold = quantile_threshold(null_statistics)
        joint_detections = sum(
            quadratic_form(
                (COEFFICIENTS[1] + error[1], COEFFICIENTS[2] + error[2]),
                precision,
            )
            > joint_threshold
            for error in errors
        )
        rows.append(
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
                "leading_detection_power": f"{detections[0]}/4096",
                "B0_detection_power": f"{detections[1]}/4096",
                "B6_detection_power": f"{detections[2]}/4096",
                "joint_dimension_eight_detection_power": (
                    f"{joint_detections}/4096"
                ),
            }
        )
    return rows


def exact_study() -> dict[str, Any]:
    design = design_matrix()
    gram, inverse_gram, estimator = estimator_matrix(design)
    errors = error_universe(estimator)
    universe_size = len(errors)

    coordinate_rows = []
    thresholds: list[Fraction] = []
    for index, (name, coefficient) in enumerate(
        zip(COEFFICIENT_NAMES, COEFFICIENTS, strict=True)
    ):
        absolute_errors = [abs(error[index]) for error in errors]
        threshold = quantile_threshold(absolute_errors)
        thresholds.append(threshold)
        covered = sum(value <= threshold for value in absolute_errors)
        detections = sum(
            abs(coefficient + error[index]) > threshold for error in errors
        )
        variance = NOISE_AMPLITUDE**2 * inverse_gram[index][index]
        coordinate_rows.append(
            {
                "coefficient": name,
                "injected_value": fraction_text(coefficient),
                "interval_half_width": fraction_text(threshold),
                "covered_replicates": covered,
                "coverage": f"{covered}/{universe_size}",
                "signal_detections": detections,
                "detection_power": f"{detections}/{universe_size}",
                "estimator_variance": fraction_text(variance),
                "squared_signal_to_noise": fraction_text(coefficient**2 / variance),
            }
        )

    dimension_eight_covariance = [
        [NOISE_AMPLITUDE**2 * inverse_gram[i][j] for j in (1, 2)]
        for i in (1, 2)
    ]
    dimension_eight_precision = invert(dimension_eight_covariance)
    null_statistics = [
        quadratic_form((error[1], error[2]), dimension_eight_precision)
        for error in errors
    ]
    joint_threshold = quantile_threshold(null_statistics)
    joint_covered = sum(value <= joint_threshold for value in null_statistics)
    dimension_eight_signal = (COEFFICIENTS[1], COEFFICIENTS[2])
    joint_detections = sum(
        quadratic_form(
            (
                dimension_eight_signal[0] + error[1],
                dimension_eight_signal[1] + error[2],
            ),
            dimension_eight_precision,
        )
        > joint_threshold
        for error in errors
    )

    remainder_coefficient = Fraction(7, 388800)
    q_by_row = [q for q in Q_VALUES for _ in I6_VALUES]
    remainder_bias_bounds = [
        sum(
            abs(estimator[index][row]) * remainder_coefficient * q_by_row[row] ** 6
            for row in range(len(design))
        )
        for index in range(3)
    ]
    leading_detection_margin = min(
        abs(COEFFICIENTS[0] + error[0]) - thresholds[0] for error in errors
    )
    require(leading_detection_margin > 0, "leading recovery has no strict margin")
    require(
        remainder_bias_bounds[0] < leading_detection_margin,
        "certified remainder can erase the leading synthetic detection",
    )

    return {
        "design_matrix": design,
        "gram_matrix": gram,
        "inverse_gram_matrix": inverse_gram,
        "estimator_matrix": estimator,
        "errors": errors,
        "coordinate_rows": coordinate_rows,
        "thresholds": thresholds,
        "dimension_eight_covariance": dimension_eight_covariance,
        "dimension_eight_precision": dimension_eight_precision,
        "joint_threshold": joint_threshold,
        "joint_covered": joint_covered,
        "joint_detections": joint_detections,
        "remainder_bias_bounds": remainder_bias_bounds,
        "leading_detection_margin": leading_detection_margin,
        "sensitivity_grid": sensitivity_grid(estimator, inverse_gram),
    }


def build_receipt() -> dict[str, Any]:
    loaded = [load_parent(contract) for contract in PARENT_CONTRACTS]
    source, custody, remainder = (item[1] for item in loaded)
    require(
        source.get("conditional_physical_candidate", {}).get("coefficients")
        == {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
        },
        "frozen FZ-12 coefficient ray drift",
    )
    require(
        custody.get("projection_scope", {}).get("includes_measurement_values") is False,
        "custody projection contains measurement values",
    )
    require(
        custody.get("source_receipt", {}).get("sha256")
        == PARENT_CONTRACTS[0]["sha256"],
        "custody projection lost the frozen source pin",
    )
    require(
        remainder.get("geometry_contract", {}).get("q_domain")
        == {"minimum": "0", "maximum": "1", "inclusive": True},
        "remainder q domain drift",
    )
    require(
        remainder.get("taylor_remainders", {}).get("R6", {}).get(
            "global_absolute_coefficient"
        )
        == "7/388800",
        "remainder coefficient drift",
    )
    require(
        remainder.get("exposure_boundary", {}).get("comparison_inputs") == [],
        "remainder parent contains comparison input",
    )

    study = exact_study()
    parent_pins = []
    for contract, (raw, value) in zip(PARENT_CONTRACTS, loaded, strict=True):
        parent_pins.append(
            {
                "path": contract["relative_path"],
                "role": contract["role"],
                "bytes": len(raw),
                "sha256": f"sha256:{sha256_hex(raw)}",
                "schema": contract["schema"],
                "status": contract["status"],
                "self_digest": value[contract["self_field"]],
            }
        )

    design_rows = []
    for index, (q, i6, features) in enumerate(
        zip(
            (q for q in Q_VALUES for _ in I6_VALUES),
            (i6 for _ in Q_VALUES for i6 in I6_VALUES),
            study["design_matrix"],
            strict=True,
        )
    ):
        design_rows.append(
            {
                "row": index,
                "q": fraction_text(q),
                "I6": fraction_text(i6),
                "features": [fraction_text(item) for item in features],
            }
        )

    remainder_rows = []
    for name, bound, threshold in zip(
        COEFFICIENT_NAMES,
        study["remainder_bias_bounds"],
        study["thresholds"],
        strict=True,
    ):
        remainder_rows.append(
            {
                "coefficient": name,
                "worst_case_bias_bound": fraction_text(bound),
                "fraction_of_interval_half_width": fraction_text(bound / threshold),
            }
        )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 667,
        "parent_pins": parent_pins,
        "study_scope": {
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
        "injected_frozen_ray": {
            "C4_over_a2": "-1/20",
            "B0_over_a4": "1/840",
            "B6_over_a4": "-1/12600",
            "B0_over_C4_squared": "10/21",
            "B6_over_C4_squared": "-2/63",
            "B6_over_B0": "-1/15",
        },
        "exact_signal_hierarchy": {
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
        "synthetic_design": {
            "row_count": 12,
            "q_values": [fraction_text(value) for value in Q_VALUES],
            "I6_values": [fraction_text(value) for value in I6_VALUES],
            "q_domain": "0 < q <= 1",
            "rows": design_rows,
            "feature_order": list(COEFFICIENT_NAMES),
            "gram_matrix": matrix_text(study["gram_matrix"]),
            "inverse_gram_matrix": matrix_text(study["inverse_gram_matrix"]),
            "full_column_rank": True,
        },
        "noise_and_calibration": {
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
        "synthetic_sensitivity_grid": {
            "q_values": [fraction_text(value) for value in Q_VALUES],
            "sigma_values": [
                fraction_text(value) for value in SENSITIVITY_NOISE_GRID
            ],
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
        "coordinate_results": study["coordinate_rows"],
        "joint_dimension_eight_result": {
            "coefficients": ["B0_over_a4", "B6_over_a4"],
            "covariance": matrix_text(study["dimension_eight_covariance"]),
            "null_statistic_threshold": fraction_text(study["joint_threshold"]),
            "covered_replicates": study["joint_covered"],
            "coverage": f"{study['joint_covered']}/4096",
            "signal_detections": study["joint_detections"],
            "detection_power": f"{study['joint_detections']}/4096",
            "resolved_in_declared_synthetic_regime": False,
        },
        "recovery_conclusion": {
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
        "full_symbol_remainder_robustness": {
            "source_bound": "abs(lambda_hat - P6) <= (7/388800) q^8",
            "response_bound": (
                "for 0 < q <= 1, abs((lambda_hat - P6)/q^2) "
                "<= (7/388800) q^6"
            ),
            "worst_case_ols_bias_bounds": remainder_rows,
            "minimum_exhaustive_leading_detection_margin": fraction_text(
                study["leading_detection_margin"]
            ),
            "leading_remainder_bias_to_detection_margin": fraction_text(
                study["remainder_bias_bounds"][0]
                / study["leading_detection_margin"]
            ),
            "worst_case_bias_less_than_detection_margin": True,
            "leading_recovery_changed": False,
            "experimental_remainder_model_supplied": False,
        },
        "exposure_and_physical_boundary": {
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
        "negative_controls": {
            "scope": (
                "contract and mutation guards, not producer-executed "
                "counterfactual coverage studies"
            ),
            "all_guards_passed": True,
            "controls": [
                {
                    "control": "leading_sign_flip",
                    "guard_passed": True,
                    "expected_failure": "frozen negative C4 coefficient",
                },
                {
                    "control": "unlinked_isotropic_term",
                    "guard_passed": True,
                    "expected_failure": "frozen B0/C4^2 ratio",
                },
                {
                    "control": "rank_six_sign_flip",
                    "guard_passed": True,
                    "expected_failure": "frozen negative B6 coefficient",
                },
                {
                    "control": "rank_deficient_direction_design",
                    "guard_passed": True,
                    "expected_failure": "three-coefficient identifiability",
                },
                {
                    "control": "q_above_certified_domain",
                    "guard_passed": True,
                    "expected_failure": "q at most one remainder contract",
                },
                {
                    "control": "partial_noise_enumeration",
                    "guard_passed": True,
                    "expected_failure": "exact 4096-replicate coverage",
                },
                {
                    "control": "public_comparison_input",
                    "guard_passed": True,
                    "expected_failure": "synthetic-only exposure boundary",
                },
            ],
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def verify_committed_receipt() -> dict[str, Any]:
    raw = RECEIPT_PATH.read_bytes()
    try:
        committed = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError("invalid committed synthetic receipt") from error
    require(isinstance(committed, dict), "synthetic receipt root is not an object")
    require(raw == base.canonical_json_bytes(committed), "noncanonical synthetic receipt")
    require(committed.get("schema") == SCHEMA, "synthetic receipt schema drift")
    require(committed.get("status") == STATUS, "synthetic receipt status drift")
    body = {key: item for key, item in committed.items() if key != "receipt_sha256"}
    require(
        committed.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        "synthetic receipt self-digest drift",
    )
    rebuilt = build_receipt()
    require(raw == base.canonical_json_bytes(rebuilt), "synthetic receipt result drift")
    return committed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    if args.write and args.verify:
        parser.error("choose --write or --verify")
    receipt = verify_committed_receipt() if args.verify else build_receipt()
    if args.write:
        RECEIPT_PATH.parent.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["recovery_conclusion"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
