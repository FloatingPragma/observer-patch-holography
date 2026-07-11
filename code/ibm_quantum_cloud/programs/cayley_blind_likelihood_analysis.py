#!/usr/bin/env python3
"""Frozen blinded likelihood analysis for the record-gated Cayley benchmark.

The analysis compares fixed process hypotheses using conditional likelihood
ratios.  The calibration channels are frozen from independent data, smoothed
with a positive preregistered Dirichlet pseudocount, and accompanied by frozen
sensitivity channels.  Their finite-sample uncertainty is not analytically
integrated, so this module deliberately does *not* call the resulting ratios
Bayes factors.

The fixed hypotheses are:

* the record-gated repair process;
* the matched lazy-heat process;
* delayed-record, shuffled-record, and inverted-record controllers;
* a globally shared, multiplicity-marginalized label/layout model; and
* a calibration-only noise model.

Every hypothesis supplies a *latent* joint law for each submitted circuit's
``(heated state, decision record, final state)``.  A calibration channel frozen
before held-out execution is then applied without fitting:

    p_observed(y | H, calibration) = sum_x C(y | x) p_latent(x | H).

Leakage is not a conditioning event.  It is a deterministic field of every
joint observed outcome and remains in the multinomial likelihood.  The input
contract also binds declared, submitted, retrieved, and counted shots so a
missing or postselected packet fails before any scientific score is emitted.

This module is deliberately independent of Qiskit and never reads credentials
or submits work.  It evaluates a programmed measurement-and-feedback channel;
it does not compare OPH with unrestricted ordinary quantum mechanics.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


LOCK_SCHEMA_VERSION = "oph.cayley-blind-analysis-lock.v1"
DATA_SCHEMA_VERSION = "oph.cayley-blind-heldout-data.v1"
REPORT_SCHEMA_VERSION = "oph.cayley-blind-likelihood-report.v1"

STATE_WIDTH = 3
STATE_CARDINALITY = 2**STATE_WIDTH
DECISION_CARDINALITY = 2
JOINT_CARDINALITY = STATE_CARDINALITY * DECISION_CARDINALITY * STATE_CARDINALITY

PRIMARY_ENDPOINT = "primary_s3"
REPAIR_MODEL = "record_gated_repair"
LABEL_MODEL = "label_only"
REQUIRED_NULL_MODELS = (
    "lazy_heat",
    "delayed_record",
    "shuffled_record",
    "inverted_record",
    LABEL_MODEL,
    "calibrated_noise",
)
REQUIRED_MODELS = (REPAIR_MODEL, *REQUIRED_NULL_MODELS)
ROW_POINT_MODELS = tuple(model for model in REQUIRED_MODELS if model != LABEL_MODEL)

PER_BACKEND_LR_THRESHOLD = 10.0
POOLED_LR_THRESHOLD = 100.0
SIMULTANEOUS_LEVEL = 0.99
SECONDARY_FAMILY_ALPHA = 0.01
CALIBRATION_CONTROL_MIN_P_VALUE = 0.01
DEFAULT_MAX_LEAKAGE_FRACTION = 0.10
CALIBRATION_UNCERTAINTY_MODE = "conditional_plugin_positive_dirichlet_sensitivity"

_OUTCOME_RE = re.compile(r"^h([0-7])\|d([01])\|f([0-7])\|l([01])$")
_QISKIT_JOINED_RE = re.compile(r"^[01]{7}$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_FORBIDDEN_BLIND_KEYS = {
    "api_key",
    "arm_mapping",
    "label_map",
    "reveal_key",
    "reveal_salt",
    "semantic_mapping",
    "token",
    "unblinded_labels",
}


class AnalysisValidationError(ValueError):
    """Raised when the frozen analysis contract cannot be satisfied."""


def canonical_json_bytes(value: Any) -> bytes:
    """Canonical JSON encoding used by every lock and report hash."""

    try:
        text = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise AnalysisValidationError(f"value is not canonical-JSON serializable: {exc}") from exc
    return text.encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def analysis_code_sha256() -> str:
    """Hash the exact analysis source that is executing."""

    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _json_copy(value: Any) -> Any:
    return json.loads(canonical_json_bytes(value).decode("utf-8"))


def _require_sha256(value: Any, field: str) -> str:
    if not isinstance(value, str) or _SHA256_RE.fullmatch(value) is None:
        raise AnalysisValidationError(f"{field} must be a lowercase SHA-256 digest")
    return value


def _reject_reveal_material(value: Any, path: str = "root") -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if str(key).lower() in _FORBIDDEN_BLIND_KEYS:
                raise AnalysisValidationError(
                    f"blind input contains forbidden reveal/credential field {path}.{key}"
                )
            _reject_reveal_material(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_reveal_material(child, f"{path}[{index}]")


def joint_index(heated: int, decision: int, final: int) -> int:
    if heated not in range(STATE_CARDINALITY):
        raise AnalysisValidationError("heated state is outside the three-qubit code space")
    if decision not in range(DECISION_CARDINALITY):
        raise AnalysisValidationError("decision record is not a bit")
    if final not in range(STATE_CARDINALITY):
        raise AnalysisValidationError("final state is outside the three-qubit code space")
    return (heated * DECISION_CARDINALITY + decision) * STATE_CARDINALITY + final


def joint_tuple(index: int) -> tuple[int, int, int]:
    if index not in range(JOINT_CARDINALITY):
        raise AnalysisValidationError("joint outcome index is outside the frozen space")
    heated_decision, final = divmod(index, STATE_CARDINALITY)
    heated, decision = divmod(heated_decision, DECISION_CARDINALITY)
    return heated, decision, final


def leakage_bit(heated: int, final: int, valid_codes: Sequence[int]) -> int:
    valid = set(int(code) for code in valid_codes)
    return int(heated not in valid or final not in valid)


def outcome_key(
    heated: int,
    decision: int,
    final: int,
    valid_codes: Sequence[int],
) -> str:
    leak = leakage_bit(heated, final, valid_codes)
    return f"h{heated}|d{decision}|f{final}|l{leak}"


def parse_outcome_key(key: str, valid_codes: Sequence[int]) -> tuple[int, int, int, int]:
    match = _OUTCOME_RE.fullmatch(str(key))
    if match is None:
        raise AnalysisValidationError(f"invalid joint outcome key {key!r}")
    heated, decision, final, leak = (int(value) for value in match.groups())
    expected_leak = leakage_bit(heated, final, valid_codes)
    if leak != expected_leak:
        raise AnalysisValidationError(
            f"outcome {key!r} has leakage={leak}, expected {expected_leak}"
        )
    return heated, decision, final, leak


def all_outcome_keys(valid_codes: Sequence[int]) -> tuple[str, ...]:
    return tuple(
        outcome_key(heated, decision, final, valid_codes)
        for heated in range(STATE_CARDINALITY)
        for decision in range(DECISION_CARDINALITY)
        for final in range(STATE_CARDINALITY)
    )


def qiskit_joined_key_to_outcome_key(
    joined_key: str,
    valid_codes: Sequence[int],
) -> str:
    """Decode Qiskit ``join_data()`` bits in strict final|decision|heated order.

    The circuit declares three classical registers in the order ``heated[3]``,
    ``decision[1]``, ``final[3]``.  Qiskit's joined count strings display the
    highest classical bits first, hence the seven compact bits are
    ``final[3] | decision[1] | heated[3]``.  Spaced ``get_counts`` keys, hex
    strings, shortened keys, and any other representation are rejected rather
    than guessed.
    """

    key = str(joined_key)
    if _QISKIT_JOINED_RE.fullmatch(key) is None:
        raise AnalysisValidationError(
            "joined Qiskit outcome must be exactly seven binary bits in "
            "final|decision|heated order"
        )
    final = int(key[0:3], 2)
    decision = int(key[3], 2)
    heated = int(key[4:7], 2)
    return outcome_key(heated, decision, final, valid_codes)


def qiskit_joined_counts_to_analysis_counts(
    joined_counts: Mapping[str, Any],
    valid_codes: Sequence[int],
    *,
    expected_shots: int | None = None,
) -> dict[str, int]:
    """Convert a complete Qiskit joined-count mapping without dropping shots."""

    if not isinstance(joined_counts, Mapping) or not joined_counts:
        raise AnalysisValidationError("joined Qiskit counts must be a nonempty mapping")
    converted: dict[str, int] = {}
    total = 0
    for joined_key, raw_count in joined_counts.items():
        if isinstance(raw_count, bool) or not isinstance(raw_count, int) or raw_count < 0:
            raise AnalysisValidationError("joined Qiskit counts must be nonnegative integers")
        key = qiskit_joined_key_to_outcome_key(str(joined_key), valid_codes)
        if key in converted:
            raise AnalysisValidationError("joined Qiskit count keys are not one-to-one")
        converted[key] = raw_count
        total += raw_count
    if expected_shots is not None:
        if isinstance(expected_shots, bool) or not isinstance(expected_shots, int):
            raise AnalysisValidationError("expected_shots must be an integer")
        if total != expected_shots:
            raise AnalysisValidationError(
                f"joined Qiskit counts contain {total} of {expected_shots} shots"
            )
    return converted


def probability_mapping_to_vector(
    probabilities: Mapping[str, Any],
    valid_codes: Sequence[int],
) -> np.ndarray:
    if not isinstance(probabilities, Mapping):
        raise AnalysisValidationError("candidate probabilities must be an outcome mapping")
    vector = np.zeros(JOINT_CARDINALITY, dtype=float)
    for key, raw_probability in probabilities.items():
        heated, decision, final, _ = parse_outcome_key(str(key), valid_codes)
        if isinstance(raw_probability, bool):
            raise AnalysisValidationError("probabilities cannot be booleans")
        try:
            probability = float(raw_probability)
        except (TypeError, ValueError) as exc:
            raise AnalysisValidationError(f"invalid probability for {key!r}") from exc
        if not math.isfinite(probability) or probability < 0.0:
            raise AnalysisValidationError("probabilities must be finite and nonnegative")
        vector[joint_index(heated, decision, final)] = probability
    total = float(np.sum(vector))
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise AnalysisValidationError(f"candidate probabilities sum to {total}, not one")
    return vector


def vector_to_probability_mapping(
    vector: Sequence[float],
    valid_codes: Sequence[int],
    *,
    include_zeros: bool = False,
) -> dict[str, float]:
    array = np.asarray(vector, dtype=float)
    if array.shape != (JOINT_CARDINALITY,):
        raise AnalysisValidationError("joint probability vector has the wrong dimension")
    result: dict[str, float] = {}
    for index, probability in enumerate(array):
        if include_zeros or probability > 0.0:
            heated, decision, final = joint_tuple(index)
            result[outcome_key(heated, decision, final, valid_codes)] = float(probability)
    return result


def _validate_stochastic_matrix(
    values: Any,
    shape: tuple[int, int],
    field: str,
) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    if matrix.shape != shape:
        raise AnalysisValidationError(f"{field} must have shape {shape}")
    if np.any(~np.isfinite(matrix)) or np.any(matrix < 0.0):
        raise AnalysisValidationError(f"{field} must be finite and nonnegative")
    if not np.allclose(matrix.sum(axis=1), 1.0, atol=1e-12, rtol=0.0):
        raise AnalysisValidationError(f"{field} rows must sum to one")
    return matrix


def _validate_positive_probability_vector(values: Any, size: int, field: str) -> np.ndarray:
    vector = np.asarray(values, dtype=float)
    if vector.shape != (size,):
        raise AnalysisValidationError(f"{field} must have shape ({size},)")
    if np.any(~np.isfinite(vector)) or np.any(vector <= 0.0):
        raise AnalysisValidationError(f"{field} must be finite and strictly positive")
    if not math.isclose(float(np.sum(vector)), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise AnalysisValidationError(f"{field} must sum to one")
    return vector


def _validate_calibration_channel(
    channel: Mapping[str, Any],
    field: str,
    *,
    require_positive: bool = True,
) -> str:
    if not isinstance(channel, Mapping):
        raise AnalysisValidationError(f"{field} must be a mapping")
    mode = channel.get("channel_mode")
    if mode == "full_joint_assignment":
        matrix = _validate_stochastic_matrix(
            channel.get("joint_assignment"),
            (JOINT_CARDINALITY, JOINT_CARDINALITY),
            f"{field}.joint_assignment",
        )
        if require_positive and np.any(matrix <= 0.0):
            raise AnalysisValidationError(
                f"{field}.joint_assignment must be strictly positive after pseudocount smoothing"
            )
        return mode
    if mode == "contamination_mixture":
        raw_epsilon = channel.get("contamination_probability")
        if isinstance(raw_epsilon, bool):
            raise AnalysisValidationError(f"{field}.contamination_probability is invalid")
        try:
            epsilon = float(raw_epsilon)
        except (TypeError, ValueError) as exc:
            raise AnalysisValidationError(
                f"{field}.contamination_probability is invalid"
            ) from exc
        if not math.isfinite(epsilon) or not 0.0 < epsilon < 1.0:
            raise AnalysisValidationError(
                f"{field}.contamination_probability must lie strictly between zero and one"
            )
        _validate_positive_probability_vector(
            channel.get("contamination_distribution"),
            JOINT_CARDINALITY,
            f"{field}.contamination_distribution",
        )
        return mode
    if mode == "factorized_assignment_diagnostic_only":
        matrices = (
            _validate_stochastic_matrix(
                channel.get("heated_assignment"),
                (STATE_CARDINALITY, STATE_CARDINALITY),
                f"{field}.heated_assignment",
            ),
            _validate_stochastic_matrix(
                channel.get("decision_assignment"),
                (DECISION_CARDINALITY, DECISION_CARDINALITY),
                f"{field}.decision_assignment",
            ),
            _validate_stochastic_matrix(
                channel.get("final_assignment"),
                (STATE_CARDINALITY, STATE_CARDINALITY),
                f"{field}.final_assignment",
            ),
        )
        if require_positive and any(np.any(matrix <= 0.0) for matrix in matrices):
            raise AnalysisValidationError(
                f"{field} assignment matrices must be strictly positive after pseudocount smoothing"
            )
        return mode
    raise AnalysisValidationError(f"{field}.channel_mode is unsupported")


def validate_calibration(calibration_id: str, calibration: Mapping[str, Any]) -> str:
    if not isinstance(calibration, Mapping):
        raise AnalysisValidationError(f"calibration {calibration_id!r} is not a mapping")
    if calibration.get("frozen_before_heldout") is not True:
        raise AnalysisValidationError(f"calibration {calibration_id!r} was not frozen")
    if calibration.get("uses_heldout_outputs") is not False:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} does not exclude held-out outputs"
        )
    _require_sha256(
        calibration.get("calibration_receipt_sha256"),
        f"calibrations.{calibration_id}.calibration_receipt_sha256",
    )
    _require_sha256(
        calibration.get("calibration_derivation_sha256"),
        f"calibrations.{calibration_id}.calibration_derivation_sha256",
    )
    if calibration.get("uncertainty_mode") != CALIBRATION_UNCERTAINTY_MODE:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} lacks the frozen conditional uncertainty contract"
        )
    raw_pseudocount = calibration.get("dirichlet_pseudocount")
    if isinstance(raw_pseudocount, bool):
        raise AnalysisValidationError(f"calibration {calibration_id!r} pseudocount is invalid")
    try:
        pseudocount = float(raw_pseudocount)
    except (TypeError, ValueError) as exc:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} pseudocount is invalid"
        ) from exc
    if not math.isfinite(pseudocount) or pseudocount <= 0.0:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} requires a positive Dirichlet pseudocount"
        )

    max_age = calibration.get("max_age_seconds")
    if isinstance(max_age, bool) or not isinstance(max_age, (int, float)):
        raise AnalysisValidationError(f"calibration {calibration_id!r} max age is invalid")
    if not math.isfinite(float(max_age)) or float(max_age) <= 0.0:
        raise AnalysisValidationError(f"calibration {calibration_id!r} max age is invalid")

    control = calibration.get("control_validation")
    if not isinstance(control, Mapping):
        raise AnalysisValidationError(f"calibration {calibration_id!r} lacks control validation")
    try:
        gof_p_value = float(control.get("gof_p_value"))
    except (TypeError, ValueError) as exc:
        raise AnalysisValidationError("calibration control p-value is invalid") from exc
    min_count = control.get("minimum_count_per_prepared_state")
    required_count = control.get("minimum_required_count")
    if (
        not math.isfinite(gof_p_value)
        or not 0.0 <= gof_p_value <= 1.0
        or isinstance(min_count, bool)
        or not isinstance(min_count, int)
        or isinstance(required_count, bool)
        or not isinstance(required_count, int)
        or min_count < 0
        or required_count <= 0
    ):
        raise AnalysisValidationError("calibration control validation fields are invalid")
    computed_control_pass = bool(
        gof_p_value >= CALIBRATION_CONTROL_MIN_P_VALUE and min_count >= required_count
    )
    if control.get("passed") is not computed_control_pass:
        raise AnalysisValidationError("calibration control pass flag is inconsistent")

    mode = _validate_calibration_channel(
        calibration,
        f"calibrations.{calibration_id}",
    )
    primary_eligible = bool(
        mode == "contamination_mixture"
        or (
            mode == "full_joint_assignment"
            and calibration.get("branch_matched_dynamic") is True
        )
    )
    if calibration.get("primary_eligible") is not primary_eligible:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} primary-eligibility flag is inconsistent"
        )

    sensitivity_channels = calibration.get("sensitivity_channels")
    if not isinstance(sensitivity_channels, list) or not sensitivity_channels:
        raise AnalysisValidationError(
            f"calibration {calibration_id!r} needs at least one frozen sensitivity channel"
        )
    sensitivity_ids: set[str] = set()
    for index, channel in enumerate(sensitivity_channels):
        if not isinstance(channel, Mapping):
            raise AnalysisValidationError("calibration sensitivity channels must be mappings")
        sensitivity_id = channel.get("sensitivity_id")
        if (
            not isinstance(sensitivity_id, str)
            or not sensitivity_id
            or sensitivity_id in sensitivity_ids
        ):
            raise AnalysisValidationError("calibration sensitivity IDs must be unique")
        sensitivity_ids.add(sensitivity_id)
        _require_sha256(
            channel.get("derivation_sha256"),
            f"calibrations.{calibration_id}.sensitivity_channels[{index}].derivation_sha256",
        )
        _validate_calibration_channel(
            channel,
            f"calibrations.{calibration_id}.sensitivity_channels[{index}]",
        )
    return mode


def convolve_calibration(
    latent_probabilities: Sequence[float],
    calibration: Mapping[str, Any],
) -> np.ndarray:
    """Apply the frozen true-row/observed-column calibration channel."""

    latent = np.asarray(latent_probabilities, dtype=float)
    if latent.shape != (JOINT_CARDINALITY,):
        raise AnalysisValidationError("latent probability vector has the wrong dimension")
    mode = calibration.get("channel_mode")
    if mode == "full_joint_assignment":
        assignment = _validate_stochastic_matrix(
            calibration["joint_assignment"],
            (JOINT_CARDINALITY, JOINT_CARDINALITY),
            "joint_assignment",
        )
        observed = latent @ assignment
    elif mode == "contamination_mixture":
        epsilon = float(calibration["contamination_probability"])
        contamination = _validate_positive_probability_vector(
            calibration["contamination_distribution"],
            JOINT_CARDINALITY,
            "contamination_distribution",
        )
        observed = (1.0 - epsilon) * latent + epsilon * contamination
    elif mode == "factorized_assignment_diagnostic_only":
        heated_assignment = _validate_stochastic_matrix(
            calibration["heated_assignment"],
            (STATE_CARDINALITY, STATE_CARDINALITY),
            "heated_assignment",
        )
        decision_assignment = _validate_stochastic_matrix(
            calibration["decision_assignment"],
            (DECISION_CARDINALITY, DECISION_CARDINALITY),
            "decision_assignment",
        )
        final_assignment = _validate_stochastic_matrix(
            calibration["final_assignment"],
            (STATE_CARDINALITY, STATE_CARDINALITY),
            "final_assignment",
        )
        cube = latent.reshape(
            STATE_CARDINALITY,
            DECISION_CARDINALITY,
            STATE_CARDINALITY,
        )
        observed_cube = np.einsum(
            "hdf,ha,db,fc->abc",
            cube,
            heated_assignment,
            decision_assignment,
            final_assignment,
            optimize=True,
        )
        observed = observed_cube.reshape(JOINT_CARDINALITY)
    else:
        raise AnalysisValidationError("calibration channel mode is unsupported")
    observed = np.maximum(observed, 0.0)
    total = float(np.sum(observed))
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise AnalysisValidationError(f"calibration convolution has mass {total}")
    return observed / total


def symmetric_assignment(size: int, error_probability: float) -> list[list[float]]:
    """Build a frozen symmetric calibration matrix for synthetic preflight only."""

    if size < 2 or not 0.0 <= error_probability < 1.0:
        raise AnalysisValidationError("invalid symmetric assignment request")
    off_diagonal = error_probability / (size - 1)
    return [
        [
            1.0 - error_probability if row == column else off_diagonal
            for column in range(size)
        ]
        for row in range(size)
    ]


def factorized_calibration_packet(
    *,
    state_error: float,
    decision_error: float,
    receipt_sha256: str,
    derivation_sha256: str | None = None,
    dirichlet_pseudocount: float = 0.5,
) -> dict[str, Any]:
    """Create a diagnostic-only factorized object for simulation/tests.

    It is intentionally ineligible for a primary pass because an independent
    readout convolution cannot represent errors that alter a feedback branch.
    """

    _require_sha256(receipt_sha256, "receipt_sha256")
    if derivation_sha256 is None:
        derivation_sha256 = receipt_sha256
    _require_sha256(derivation_sha256, "derivation_sha256")
    lower_state = max(1e-6, state_error * 0.75)
    lower_decision = max(1e-6, decision_error * 0.75)
    upper_state = min(0.99, max(1e-6, state_error * 1.25))
    upper_decision = min(0.99, max(1e-6, decision_error * 1.25))
    return {
        "frozen_before_heldout": True,
        "uses_heldout_outputs": False,
        "calibration_receipt_sha256": receipt_sha256,
        "calibration_derivation_sha256": derivation_sha256,
        "uncertainty_mode": CALIBRATION_UNCERTAINTY_MODE,
        "dirichlet_pseudocount": dirichlet_pseudocount,
        "max_age_seconds": 86400.0,
        "control_validation": {
            "gof_p_value": 1.0,
            "minimum_count_per_prepared_state": 512,
            "minimum_required_count": 128,
            "passed": True,
        },
        "channel_mode": "factorized_assignment_diagnostic_only",
        "primary_eligible": False,
        "heated_assignment": symmetric_assignment(STATE_CARDINALITY, state_error),
        "decision_assignment": symmetric_assignment(
            DECISION_CARDINALITY,
            decision_error,
        ),
        "final_assignment": symmetric_assignment(STATE_CARDINALITY, state_error),
        "sensitivity_channels": [
            {
                "sensitivity_id": "lower_error",
                "derivation_sha256": derivation_sha256,
                "channel_mode": "factorized_assignment_diagnostic_only",
                "heated_assignment": symmetric_assignment(STATE_CARDINALITY, lower_state),
                "decision_assignment": symmetric_assignment(
                    DECISION_CARDINALITY, lower_decision
                ),
                "final_assignment": symmetric_assignment(STATE_CARDINALITY, lower_state),
            },
            {
                "sensitivity_id": "higher_error",
                "derivation_sha256": derivation_sha256,
                "channel_mode": "factorized_assignment_diagnostic_only",
                "heated_assignment": symmetric_assignment(STATE_CARDINALITY, upper_state),
                "decision_assignment": symmetric_assignment(
                    DECISION_CARDINALITY, upper_decision
                ),
                "final_assignment": symmetric_assignment(STATE_CARDINALITY, upper_state),
            },
        ],
    }


def contamination_calibration_packet(
    *,
    contamination_probability: float,
    sensitivity_probabilities: Sequence[float],
    receipt_sha256: str,
    derivation_sha256: str,
    dirichlet_pseudocount: float = 0.5,
    control_gof_p_value: float = 1.0,
    minimum_count_per_prepared_state: int = 512,
    minimum_required_count: int = 128,
    max_age_seconds: float = 86400.0,
) -> dict[str, Any]:
    """Build a positive, explicit contamination prior for a tractable run.

    This is not claimed to be a branch-matched hardware channel.  It is a
    preregistered conditional noise model whose width is exposed through the
    mandatory sensitivity alternatives.
    """

    _require_sha256(receipt_sha256, "receipt_sha256")
    _require_sha256(derivation_sha256, "derivation_sha256")
    uniform = [1.0 / JOINT_CARDINALITY] * JOINT_CARDINALITY
    channels = []
    for index, probability in enumerate(sensitivity_probabilities):
        channels.append(
            {
                "sensitivity_id": f"contamination_{index}",
                "derivation_sha256": derivation_sha256,
                "channel_mode": "contamination_mixture",
                "contamination_probability": float(probability),
                "contamination_distribution": uniform,
            }
        )
    packet = {
        "frozen_before_heldout": True,
        "uses_heldout_outputs": False,
        "calibration_receipt_sha256": receipt_sha256,
        "calibration_derivation_sha256": derivation_sha256,
        "uncertainty_mode": CALIBRATION_UNCERTAINTY_MODE,
        "dirichlet_pseudocount": dirichlet_pseudocount,
        "max_age_seconds": max_age_seconds,
        "control_validation": {
            "gof_p_value": control_gof_p_value,
            "minimum_count_per_prepared_state": minimum_count_per_prepared_state,
            "minimum_required_count": minimum_required_count,
            "passed": bool(
                control_gof_p_value >= CALIBRATION_CONTROL_MIN_P_VALUE
                and minimum_count_per_prepared_state >= minimum_required_count
            ),
        },
        "channel_mode": "contamination_mixture",
        "primary_eligible": True,
        "contamination_probability": float(contamination_probability),
        "contamination_distribution": uniform,
        "sensitivity_channels": channels,
    }
    validate_calibration("synthetic_contamination", packet)
    return packet


def default_model_metadata() -> dict[str, Any]:
    return {
        REPAIR_MODEL: {
            "kind": "record_conditioned_repair",
            "fixed_without_heldout_fitting": True,
        },
        "lazy_heat": {
            "kind": "depth_matched_open_loop_lazy_heat",
            "fixed_without_heldout_fitting": True,
        },
        "delayed_record": {
            "kind": "record_intervention_null",
            "fixed_without_heldout_fitting": True,
        },
        "shuffled_record": {
            "kind": "record_intervention_null",
            "fixed_without_heldout_fitting": True,
        },
        "inverted_record": {
            "kind": "record_intervention_null",
            "fixed_without_heldout_fitting": True,
        },
        LABEL_MODEL: {
            "kind": "label_layout_mixture_null",
            "fixed_without_heldout_fitting": True,
            "global_shared_mapping_marginal": True,
        },
        "calibrated_noise": {
            "kind": "controller_independent_calibrated_noise_null",
            "fixed_without_heldout_fitting": True,
        },
    }


def _validate_valid_codes(raw_codes: Any, row_id: str) -> tuple[int, ...]:
    if not isinstance(raw_codes, list) or not raw_codes:
        raise AnalysisValidationError(f"row {row_id!r} needs nonempty valid_codes")
    if any(isinstance(code, bool) or not isinstance(code, int) for code in raw_codes):
        raise AnalysisValidationError(f"row {row_id!r} valid_codes must be integers")
    codes = tuple(raw_codes)
    if len(set(codes)) != len(codes) or any(
        code not in range(STATE_CARDINALITY) for code in codes
    ):
        raise AnalysisValidationError(f"row {row_id!r} has invalid or duplicate codes")
    return codes


def build_candidate_provenance(
    candidate_probabilities: Mapping[str, Mapping[str, Any]],
    derivation_sha256_by_model: Mapping[str, str],
) -> dict[str, dict[str, str]]:
    """Bind every frozen point-hypothesis table to data and derivation hashes."""

    if set(candidate_probabilities) != set(derivation_sha256_by_model):
        raise AnalysisValidationError("candidate provenance model set is incomplete")
    result: dict[str, dict[str, str]] = {}
    for model, probabilities in candidate_probabilities.items():
        result[model] = {
            "probability_table_sha256": sha256_json(probabilities),
            "derivation_sha256": _require_sha256(
                derivation_sha256_by_model[model],
                f"candidate derivation for {model}",
            ),
        }
    return result


def _validate_candidate_provenance(
    candidates: Mapping[str, Any],
    provenance: Any,
    *,
    field: str,
) -> None:
    if not isinstance(provenance, Mapping) or set(provenance) != set(candidates):
        raise AnalysisValidationError(f"{field} does not match the candidate table set")
    for model, probabilities in candidates.items():
        record = provenance[model]
        if not isinstance(record, Mapping):
            raise AnalysisValidationError(f"{field}.{model} must be a mapping")
        _require_sha256(record.get("derivation_sha256"), f"{field}.{model}.derivation_sha256")
        claimed_table_hash = _require_sha256(
            record.get("probability_table_sha256"),
            f"{field}.{model}.probability_table_sha256",
        )
        if claimed_table_hash != sha256_json(probabilities):
            raise AnalysisValidationError(f"{field}.{model} probability-table hash mismatch")


def _validate_expected_row(
    row: Mapping[str, Any],
    calibrations: Mapping[str, Any],
    models: Sequence[str],
) -> None:
    row_id = row.get("row_id")
    if not isinstance(row_id, str) or not row_id:
        raise AnalysisValidationError("every expected row needs a nonempty row_id")
    opaque_id = row.get("opaque_id")
    if not isinstance(opaque_id, str) or not opaque_id or opaque_id != row_id:
        raise AnalysisValidationError(
            f"row {row_id!r} must use its unique opaque circuit ID as row_id"
        )
    for field in (
        "endpoint",
        "backend_role",
        "backend_name",
        "layout_id",
        "calibration_id",
    ):
        if not isinstance(row.get(field), str) or not row[field]:
            raise AnalysisValidationError(f"row {row_id!r} needs {field}")
    _require_sha256(
        row.get("logical_qpy_sha256"),
        f"expected_rows.{row_id}.logical_qpy_sha256",
    )
    physical_layout = row.get("physical_layout")
    if (
        not isinstance(physical_layout, list)
        or len(physical_layout) != 4
        or any(isinstance(qubit, bool) or not isinstance(qubit, int) or qubit < 0 for qubit in physical_layout)
        or len(set(physical_layout)) != 4
    ):
        raise AnalysisValidationError(
            f"row {row_id!r} physical_layout must contain four distinct qubit indices"
        )
    if row["calibration_id"] not in calibrations:
        raise AnalysisValidationError(f"row {row_id!r} names an unknown calibration")
    shots = row.get("shots")
    if isinstance(shots, bool) or not isinstance(shots, int) or shots <= 0:
        raise AnalysisValidationError(f"row {row_id!r} shots must be a positive integer")
    raw_leakage_limit = row.get("max_leakage_fraction")
    if isinstance(raw_leakage_limit, bool):
        raise AnalysisValidationError(f"row {row_id!r} leakage limit is invalid")
    try:
        leakage_limit = float(raw_leakage_limit)
    except (TypeError, ValueError) as exc:
        raise AnalysisValidationError(f"row {row_id!r} leakage limit is invalid") from exc
    if not math.isfinite(leakage_limit) or not 0.0 <= leakage_limit < 1.0:
        raise AnalysisValidationError(f"row {row_id!r} leakage limit is invalid")
    valid_codes = _validate_valid_codes(row.get("valid_codes"), row_id)
    candidates = row.get("candidate_probabilities")
    if not isinstance(candidates, Mapping):
        raise AnalysisValidationError(f"row {row_id!r} lacks candidate probabilities")
    if set(candidates) != set(ROW_POINT_MODELS):
        raise AnalysisValidationError(
            f"row {row_id!r} candidate set does not match the frozen model set"
        )
    for model in ROW_POINT_MODELS:
        probability_mapping_to_vector(candidates[model], valid_codes)
    _validate_candidate_provenance(
        candidates,
        row.get("candidate_provenance"),
        field=f"expected_rows.{row_id}.candidate_provenance",
    )


def _validate_secondary_specs(
    specs: Any,
    row_ids: set[str],
    models: set[str],
) -> None:
    if not isinstance(specs, list) or not specs:
        raise AnalysisValidationError("at least one frozen secondary test is required")
    test_ids: set[str] = set()
    for spec in specs:
        if not isinstance(spec, Mapping):
            raise AnalysisValidationError("secondary test specifications must be mappings")
        test_id = spec.get("test_id")
        if not isinstance(test_id, str) or not test_id or test_id in test_ids:
            raise AnalysisValidationError("secondary test IDs must be unique and nonempty")
        test_ids.add(test_id)
        selected_rows = spec.get("row_ids")
        if not isinstance(selected_rows, list) or not selected_rows:
            raise AnalysisValidationError(f"secondary test {test_id!r} has no rows")
        if len(set(selected_rows)) != len(selected_rows) or not set(selected_rows) <= row_ids:
            raise AnalysisValidationError(f"secondary test {test_id!r} has invalid row IDs")
        if spec.get("model", REPAIR_MODEL) not in models:
            raise AnalysisValidationError(f"secondary test {test_id!r} has unknown model")


def build_label_layout_component(
    *,
    component_id: str,
    prior_weight: float,
    row_probabilities: Mapping[str, Mapping[str, Any]],
    component_derivation_sha256: str,
    row_derivation_sha256: Mapping[str, str],
) -> dict[str, Any]:
    if set(row_probabilities) != set(row_derivation_sha256):
        raise AnalysisValidationError("label component row provenance is incomplete")
    return {
        "component_id": component_id,
        "prior_weight": prior_weight,
        "component_derivation_sha256": _require_sha256(
            component_derivation_sha256,
            f"label component {component_id} derivation",
        ),
        "row_probabilities": _json_copy(row_probabilities),
        "row_provenance": {
            row_id: {
                "probability_table_sha256": sha256_json(probabilities),
                "derivation_sha256": _require_sha256(
                    row_derivation_sha256[row_id],
                    f"label component {component_id} row {row_id} derivation",
                ),
            }
            for row_id, probabilities in row_probabilities.items()
        },
    }


def _validate_label_layout_model(
    label_model: Any,
    rows_by_id: Mapping[str, Mapping[str, Any]],
) -> None:
    if not isinstance(label_model, Mapping):
        raise AnalysisValidationError("analysis lock lacks the global label/layout model")
    if label_model.get("mapping_scope") != "global_shared_across_primary_rows":
        raise AnalysisValidationError("label/layout mapping scope is not globally shared")
    _require_sha256(
        label_model.get("component_set_derivation_sha256"),
        "label_layout_model.component_set_derivation_sha256",
    )
    primary_ids = {
        row_id for row_id, row in rows_by_id.items() if row["endpoint"] == PRIMARY_ENDPOINT
    }
    components = label_model.get("components")
    if not isinstance(components, list) or len(components) < 2:
        raise AnalysisValidationError("label/layout model needs at least two frozen components")
    component_ids: set[str] = set()
    weight_total = 0.0
    for index, component in enumerate(components):
        if not isinstance(component, Mapping):
            raise AnalysisValidationError("label/layout components must be mappings")
        component_id = component.get("component_id")
        if (
            not isinstance(component_id, str)
            or not component_id
            or component_id in component_ids
        ):
            raise AnalysisValidationError("label/layout component IDs must be unique")
        component_ids.add(component_id)
        _require_sha256(
            component.get("component_derivation_sha256"),
            f"label_layout_model.components[{index}].component_derivation_sha256",
        )
        raw_weight = component.get("prior_weight")
        if isinstance(raw_weight, bool):
            raise AnalysisValidationError("label/layout component weight is invalid")
        try:
            weight = float(raw_weight)
        except (TypeError, ValueError) as exc:
            raise AnalysisValidationError("label/layout component weight is invalid") from exc
        if not math.isfinite(weight) or weight <= 0.0:
            raise AnalysisValidationError("label/layout component weights must be positive")
        weight_total += weight

        row_probabilities = component.get("row_probabilities")
        row_provenance = component.get("row_provenance")
        if not isinstance(row_probabilities, Mapping) or set(row_probabilities) != primary_ids:
            raise AnalysisValidationError(
                f"label/layout component {component_id!r} must cover every primary circuit"
            )
        if not isinstance(row_provenance, Mapping) or set(row_provenance) != primary_ids:
            raise AnalysisValidationError(
                f"label/layout component {component_id!r} provenance is incomplete"
            )
        for row_id, probabilities in row_probabilities.items():
            probability_mapping_to_vector(probabilities, rows_by_id[row_id]["valid_codes"])
            provenance = row_provenance[row_id]
            if not isinstance(provenance, Mapping):
                raise AnalysisValidationError("label/layout row provenance must be a mapping")
            _require_sha256(
                provenance.get("derivation_sha256"),
                f"label_layout_model.{component_id}.{row_id}.derivation_sha256",
            )
            claimed_hash = _require_sha256(
                provenance.get("probability_table_sha256"),
                f"label_layout_model.{component_id}.{row_id}.probability_table_sha256",
            )
            if claimed_hash != sha256_json(probabilities):
                raise AnalysisValidationError(
                    f"label/layout component {component_id!r} row {row_id!r} hash mismatch"
                )
    if not math.isclose(weight_total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise AnalysisValidationError("label/layout component prior weights must sum to one")


def _validate_lock_body(lock: Mapping[str, Any], *, verify_code_hash: bool) -> None:
    _reject_reveal_material(lock)
    if lock.get("schema_version") != LOCK_SCHEMA_VERSION:
        raise AnalysisValidationError("analysis lock schema version mismatch")
    if lock.get("blinded") is not True or lock.get("revealed") is not False:
        raise AnalysisValidationError("analysis lock must remain blinded and unrevealed")
    if lock.get("lock_state") != "frozen_before_heldout":
        raise AnalysisValidationError("analysis lock was not frozen before held-out execution")
    if lock.get("primary_endpoint") != PRIMARY_ENDPOINT:
        raise AnalysisValidationError("the one primary endpoint must be pooled S3")
    _require_sha256(
        lock.get("catalog_precommitment_sha256"),
        "catalog_precommitment_sha256",
    )
    locked_code_hash = _require_sha256(
        lock.get("analysis_code_sha256"),
        "analysis_code_sha256",
    )
    if verify_code_hash and locked_code_hash != analysis_code_sha256():
        raise AnalysisValidationError("analysis source changed after the lock was frozen")

    thresholds = lock.get("thresholds")
    expected_thresholds = {
        "per_backend_likelihood_ratio": PER_BACKEND_LR_THRESHOLD,
        "pooled_likelihood_ratio": POOLED_LR_THRESHOLD,
        "simultaneous_prediction_level": SIMULTANEOUS_LEVEL,
        "secondary_holm_family_alpha": SECONDARY_FAMILY_ALPHA,
        "calibration_control_min_p_value": CALIBRATION_CONTROL_MIN_P_VALUE,
    }
    if thresholds != expected_thresholds:
        raise AnalysisValidationError("analysis thresholds differ from the frozen protocol")

    models = lock.get("models")
    if not isinstance(models, list) or tuple(models) != REQUIRED_MODELS:
        raise AnalysisValidationError("the frozen model order or required null set changed")
    metadata = lock.get("model_metadata")
    if not isinstance(metadata, Mapping) or set(metadata) != set(models):
        raise AnalysisValidationError("model metadata does not match the candidate set")
    for model in models:
        if metadata[model].get("fixed_without_heldout_fitting") is not True:
            raise AnalysisValidationError("all hypotheses must be fixed without held-out fitting")
    if metadata[LABEL_MODEL].get("global_shared_mapping_marginal") is not True:
        raise AnalysisValidationError("label/layout multiplicity is not a global marginal")

    execution_contract = lock.get("execution_contract")
    if execution_contract != {
        "row_granularity": "one_locked_row_per_opaque_circuit",
        "qiskit_joined_bit_order": "final[3]|decision[1]|heated[3]",
        "calibration_uncertainty": CALIBRATION_UNCERTAINTY_MODE,
        "evidence_measure": "conditional_likelihood_ratio_not_bayes_factor",
    }:
        raise AnalysisValidationError("analysis execution contract changed")

    calibrations = lock.get("calibrations")
    if not isinstance(calibrations, Mapping) or not calibrations:
        raise AnalysisValidationError("analysis lock has no calibration channels")
    for calibration_id, calibration in calibrations.items():
        validate_calibration(str(calibration_id), calibration)

    rows = lock.get("expected_rows")
    if not isinstance(rows, list) or not rows:
        raise AnalysisValidationError("analysis lock has no expected rows")
    row_ids: set[str] = set()
    primary_backends: set[str] = set()
    role_bindings: dict[str, tuple[str, str, tuple[int, ...]]] = {}
    for row in rows:
        _validate_expected_row(row, calibrations, models)
        row_id = row["row_id"]
        if row_id in row_ids:
            raise AnalysisValidationError(f"duplicate expected row {row_id!r}")
        row_ids.add(row_id)
        binding = (
            row["backend_name"],
            row["layout_id"],
            tuple(row["physical_layout"]),
        )
        previous_binding = role_bindings.setdefault(row["backend_role"], binding)
        if previous_binding != binding:
            raise AnalysisValidationError("one backend role maps to multiple backend/layout bindings")
        if row["endpoint"] == PRIMARY_ENDPOINT:
            primary_backends.add(row["backend_name"])
    if len(primary_backends) < 2:
        raise AnalysisValidationError("primary S3 endpoint needs two independent backends")
    if len({row["shots"] for row in rows}) != 1:
        raise AnalysisValidationError("all locked benchmark circuits must use the same shot count")
    rows_by_id = {row["row_id"]: row for row in rows}
    _validate_label_layout_model(lock.get("label_layout_model"), rows_by_id)
    _validate_secondary_specs(lock.get("secondary_tests"), row_ids, set(models))


def validate_analysis_lock(lock: Mapping[str, Any], *, verify_code_hash: bool = True) -> None:
    if not isinstance(lock, Mapping):
        raise AnalysisValidationError("analysis lock must be a mapping")
    claimed_hash = _require_sha256(lock.get("analysis_lock_sha256"), "analysis_lock_sha256")
    unhashed = dict(lock)
    unhashed.pop("analysis_lock_sha256", None)
    if sha256_json(unhashed) != claimed_hash:
        raise AnalysisValidationError("analysis lock hash mismatch")
    _validate_lock_body(unhashed, verify_code_hash=verify_code_hash)


def build_analysis_lock(
    *,
    expected_rows: Sequence[Mapping[str, Any]],
    calibrations: Mapping[str, Mapping[str, Any]],
    secondary_tests: Sequence[Mapping[str, Any]],
    catalog_precommitment_sha256: str,
    label_layout_model: Mapping[str, Any],
    created_utc: str | None = None,
    model_metadata: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Create and self-verify the immutable pre-held-out analysis lock."""

    if created_utc is None:
        created_utc = datetime.now(timezone.utc).isoformat()
    body = {
        "schema_version": LOCK_SCHEMA_VERSION,
        "created_utc": created_utc,
        "lock_state": "frozen_before_heldout",
        "blinded": True,
        "revealed": False,
        "catalog_precommitment_sha256": catalog_precommitment_sha256,
        "analysis_code_sha256": analysis_code_sha256(),
        "primary_endpoint": PRIMARY_ENDPOINT,
        "models": list(REQUIRED_MODELS),
        "model_metadata": _json_copy(
            default_model_metadata() if model_metadata is None else model_metadata
        ),
        "thresholds": {
            "per_backend_likelihood_ratio": PER_BACKEND_LR_THRESHOLD,
            "pooled_likelihood_ratio": POOLED_LR_THRESHOLD,
            "simultaneous_prediction_level": SIMULTANEOUS_LEVEL,
            "secondary_holm_family_alpha": SECONDARY_FAMILY_ALPHA,
            "calibration_control_min_p_value": CALIBRATION_CONTROL_MIN_P_VALUE,
        },
        "execution_contract": {
            "row_granularity": "one_locked_row_per_opaque_circuit",
            "qiskit_joined_bit_order": "final[3]|decision[1]|heated[3]",
            "calibration_uncertainty": CALIBRATION_UNCERTAINTY_MODE,
            "evidence_measure": "conditional_likelihood_ratio_not_bayes_factor",
        },
        "calibrations": _json_copy(calibrations),
        "expected_rows": _json_copy(list(expected_rows)),
        "label_layout_model": _json_copy(label_layout_model),
        "secondary_tests": _json_copy(list(secondary_tests)),
        "likelihood_contract": {
            "family": "stratified_complete_joint_conditional_likelihood_ratios",
            "calibration_convolution": "p_obs(y)=sum_x C(y|x)*p_latent(x)",
            "multinomial_constants_included": True,
            "calibration_refit_on_heldout": False,
            "positive_dirichlet_pseudocount_required": True,
            "sensitivity_bounds_required": True,
            "leakage_conditioning": False,
        },
        "claim_boundary": (
            "Identifies a programmed record-gated channel against frozen controller nulls; "
            "not OPH against unrestricted quantum mechanics."
        ),
    }
    _validate_lock_body(body, verify_code_hash=True)
    locked = dict(body)
    locked["analysis_lock_sha256"] = sha256_json(body)
    validate_analysis_lock(locked)
    return locked


def _counts_to_vector(
    counts: Mapping[str, Any],
    valid_codes: Sequence[int],
) -> np.ndarray:
    if not isinstance(counts, Mapping) or not counts:
        raise AnalysisValidationError("joint counts must be a nonempty mapping")
    vector = np.zeros(JOINT_CARDINALITY, dtype=np.int64)
    for key, raw_count in counts.items():
        heated, decision, final, _ = parse_outcome_key(str(key), valid_codes)
        if isinstance(raw_count, bool) or not isinstance(raw_count, int) or raw_count < 0:
            raise AnalysisValidationError("joint counts must be nonnegative integers")
        vector[joint_index(heated, decision, final)] = raw_count
    return vector


def _validate_data_row(data_row: Mapping[str, Any], expected_row: Mapping[str, Any]) -> np.ndarray:
    row_id = expected_row["row_id"]
    if data_row.get("row_id") != row_id:
        raise AnalysisValidationError(f"data row ID does not match expected row {row_id!r}")
    for field in (
        "opaque_id",
        "logical_qpy_sha256",
        "backend_role",
        "backend_name",
        "layout_id",
        "physical_layout",
    ):
        if data_row.get(field) != expected_row[field]:
            raise AnalysisValidationError(
                f"row {row_id!r} execution provenance does not match locked {field}"
            )
    for field in ("job_id", "group_id"):
        if not isinstance(data_row.get(field), str) or not data_row[field]:
            raise AnalysisValidationError(f"row {row_id!r} lacks {field}")
    for field in (
        "submission_event_sha256",
        "harvest_event_sha256",
        "raw_joined_counts_sha256",
    ):
        _require_sha256(data_row.get(field), f"data row {row_id}.{field}")
    raw_age = data_row.get("calibration_age_seconds")
    if isinstance(raw_age, bool):
        raise AnalysisValidationError(f"row {row_id!r} calibration age is invalid")
    try:
        age = float(raw_age)
    except (TypeError, ValueError) as exc:
        raise AnalysisValidationError(f"row {row_id!r} calibration age is invalid") from exc
    if not math.isfinite(age) or age < 0.0:
        raise AnalysisValidationError(f"row {row_id!r} calibration age is invalid")
    expected_shots = int(expected_row["shots"])
    for field in ("declared_shots", "submitted_shots", "retrieved_shots"):
        value = data_row.get(field)
        if isinstance(value, bool) or not isinstance(value, int) or value != expected_shots:
            raise AnalysisValidationError(
                f"row {row_id!r} {field} does not equal locked shots={expected_shots}"
            )
    if data_row.get("excluded_shots") != 0:
        raise AnalysisValidationError(f"row {row_id!r} excludes shots")
    if data_row.get("postselected") is not False:
        raise AnalysisValidationError(f"row {row_id!r} is postselected")
    if data_row.get("all_outcomes_included") is not True:
        raise AnalysisValidationError(f"row {row_id!r} does not include all outcomes")
    counts = _counts_to_vector(data_row.get("counts"), expected_row["valid_codes"])
    counted_shots = int(np.sum(counts))
    if counted_shots != expected_shots:
        raise AnalysisValidationError(
            f"row {row_id!r} counted {counted_shots} of {expected_shots} shots"
        )
    return counts


def seal_data_packet(
    *,
    analysis_lock_sha256: str,
    rows: Sequence[Mapping[str, Any]],
    manifest_sha256: str,
    submission_journal_sha256: str,
    harvest_journal_sha256: str,
    source_kind: str,
    created_utc: str | None = None,
) -> dict[str, Any]:
    if created_utc is None:
        created_utc = datetime.now(timezone.utc).isoformat()
    packet = {
        "schema_version": DATA_SCHEMA_VERSION,
        "created_utc": created_utc,
        "analysis_lock_sha256": analysis_lock_sha256,
        "manifest_sha256": _require_sha256(manifest_sha256, "manifest_sha256"),
        "submission_journal_sha256": _require_sha256(
            submission_journal_sha256, "submission_journal_sha256"
        ),
        "harvest_journal_sha256": _require_sha256(
            harvest_journal_sha256, "harvest_journal_sha256"
        ),
        "source_kind": source_kind,
        "blinded": True,
        "revealed": False,
        "jobs_complete": True,
        "all_submitted_jobs_included": True,
        "rows": _json_copy(list(rows)),
    }
    _reject_reveal_material(packet)
    packet["data_packet_sha256"] = sha256_json(packet)
    return packet


def validate_data_packet(
    lock: Mapping[str, Any],
    data: Mapping[str, Any],
) -> dict[str, np.ndarray]:
    if not isinstance(data, Mapping):
        raise AnalysisValidationError("held-out data packet must be a mapping")
    _reject_reveal_material(data)
    if data.get("schema_version") != DATA_SCHEMA_VERSION:
        raise AnalysisValidationError("held-out data schema version mismatch")
    if data.get("blinded") is not True or data.get("revealed") is not False:
        raise AnalysisValidationError("held-out packet is not blinded")
    if data.get("jobs_complete") is not True:
        raise AnalysisValidationError("held-out jobs are not marked complete")
    if data.get("all_submitted_jobs_included") is not True:
        raise AnalysisValidationError("held-out packet omits submitted jobs")
    if data.get("analysis_lock_sha256") != lock["analysis_lock_sha256"]:
        raise AnalysisValidationError("held-out packet refers to a different analysis lock")
    for field in ("manifest_sha256", "submission_journal_sha256", "harvest_journal_sha256"):
        _require_sha256(data.get(field), field)
    if data.get("source_kind") not in ("ibm_qpu_hardware", "synthetic_preflight"):
        raise AnalysisValidationError("held-out packet source kind is invalid")
    claimed_hash = _require_sha256(data.get("data_packet_sha256"), "data_packet_sha256")
    unhashed = dict(data)
    unhashed.pop("data_packet_sha256", None)
    if sha256_json(unhashed) != claimed_hash:
        raise AnalysisValidationError("held-out data packet hash mismatch")

    expected = {row["row_id"]: row for row in lock["expected_rows"]}
    rows = data.get("rows")
    if not isinstance(rows, list):
        raise AnalysisValidationError("held-out rows must be a list")
    observed: dict[str, Mapping[str, Any]] = {}
    for row in rows:
        if not isinstance(row, Mapping) or not isinstance(row.get("row_id"), str):
            raise AnalysisValidationError("malformed held-out row")
        if row["row_id"] in observed:
            raise AnalysisValidationError(f"duplicate held-out row {row['row_id']!r}")
        observed[row["row_id"]] = row
    if set(observed) != set(expected):
        missing = sorted(set(expected) - set(observed))
        extra = sorted(set(observed) - set(expected))
        raise AnalysisValidationError(f"held-out row mismatch: missing={missing}, extra={extra}")
    return {
        row_id: _validate_data_row(observed[row_id], expected_row)
        for row_id, expected_row in expected.items()
    }


def multinomial_log_likelihood(counts: Sequence[int], probabilities: Sequence[float]) -> float:
    count_array = np.asarray(counts, dtype=np.int64)
    probability_array = np.asarray(probabilities, dtype=float)
    if count_array.shape != (JOINT_CARDINALITY,) or probability_array.shape != (
        JOINT_CARDINALITY,
    ):
        raise AnalysisValidationError("multinomial vectors have the wrong dimension")
    if np.any(count_array < 0):
        raise AnalysisValidationError("multinomial counts must be nonnegative")
    if np.any(probability_array < 0.0) or not math.isclose(
        float(np.sum(probability_array)),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-10,
    ):
        raise AnalysisValidationError("multinomial probabilities are invalid")
    positive = count_array > 0
    if np.any((probability_array <= 0.0) & positive):
        return -math.inf
    shots = int(np.sum(count_array))
    coefficient = math.lgamma(shots + 1.0) - sum(
        math.lgamma(int(count) + 1.0) for count in count_array
    )
    score = coefficient + float(
        np.sum(count_array[positive] * np.log(probability_array[positive]))
    )
    return float(score)


def _candidate_observed_probabilities(
    expected_row: Mapping[str, Any],
    calibration: Mapping[str, Any],
    model: str,
) -> np.ndarray:
    latent = probability_mapping_to_vector(
        expected_row["candidate_probabilities"][model],
        expected_row["valid_codes"],
    )
    return convolve_calibration(latent, calibration)


def _safe_bayes_factor(log_bayes_factor: float) -> float | str:
    if math.isnan(log_bayes_factor):
        raise AnalysisValidationError("Bayes-factor comparison produced NaN")
    if log_bayes_factor > 709.0:
        return "Infinity"
    if log_bayes_factor < -745.0:
        return 0.0
    return float(math.exp(log_bayes_factor))


def simultaneous_hoeffding_envelope(
    row_records: Sequence[tuple[str, np.ndarray, np.ndarray]],
    *,
    level: float = SIMULTANEOUS_LEVEL,
) -> dict[str, Any]:
    """Bonferroni-Hoeffding simultaneous prediction envelope.

    Each multinomial cell is marginally binomial.  Hoeffding plus a union
    bound over every locked row and all 128 cells gives family coverage at
    least ``level`` without asymptotic or independence assumptions between
    cells.
    """

    if not row_records or not 0.0 < level < 1.0:
        raise AnalysisValidationError("invalid simultaneous envelope request")
    family_size = len(row_records) * JOINT_CARDINALITY
    alpha = 1.0 - level
    log_term = math.log((2.0 * family_size) / alpha)
    violations: list[dict[str, Any]] = []
    max_absolute_error = 0.0
    for row_id, counts, probabilities in row_records:
        shots = int(np.sum(counts))
        epsilon = math.sqrt(log_term / (2.0 * shots))
        frequencies = counts / shots
        max_absolute_error = max(
            max_absolute_error,
            float(np.max(np.abs(frequencies - probabilities))),
        )
        for index, (count, probability) in enumerate(zip(counts, probabilities)):
            lower_probability = max(0.0, float(probability) - epsilon)
            upper_probability = min(1.0, float(probability) + epsilon)
            lower_count = max(0, math.ceil(shots * lower_probability - 1e-12))
            upper_count = min(shots, math.floor(shots * upper_probability + 1e-12))
            if int(count) < lower_count or int(count) > upper_count:
                heated, decision, final = joint_tuple(index)
                violations.append(
                    {
                        "row_id": row_id,
                        "heated": heated,
                        "decision": decision,
                        "final": final,
                        "observed_count": int(count),
                        "predicted_probability": float(probability),
                        "allowed_count": [lower_count, upper_count],
                    }
                )
    return {
        "method": "Bonferroni-Hoeffding simultaneous multinomial-cell envelope",
        "level": level,
        "family_size": family_size,
        "pass": not violations,
        "violation_count": len(violations),
        "violations": violations,
        "max_absolute_fraction_error": max_absolute_error,
    }


def conservative_hoeffding_gof_p_value(
    row_records: Sequence[tuple[str, np.ndarray, np.ndarray]],
) -> tuple[float, float, int]:
    """Return a union-bound p-value for the largest calibrated cell residual."""

    if not row_records:
        raise AnalysisValidationError("secondary test has no rows")
    statistic = 0.0
    cells = len(row_records) * JOINT_CARDINALITY
    for _, counts, probabilities in row_records:
        shots = int(np.sum(counts))
        frequencies = counts / shots
        row_statistic = float(
            np.max(np.sqrt(2.0 * shots) * np.abs(frequencies - probabilities))
        )
        statistic = max(statistic, row_statistic)
    p_value = min(1.0, 2.0 * cells * math.exp(-(statistic**2)))
    return p_value, statistic, cells


def holm_adjust(
    tests: Sequence[tuple[str, float]],
    *,
    alpha: float = SECONDARY_FAMILY_ALPHA,
) -> list[dict[str, Any]]:
    if not tests or not 0.0 < alpha < 1.0:
        raise AnalysisValidationError("invalid Holm family")
    seen: set[str] = set()
    normalized: list[tuple[str, float]] = []
    for test_id, raw_p_value in tests:
        if test_id in seen:
            raise AnalysisValidationError("Holm test IDs must be unique")
        seen.add(test_id)
        p_value = float(raw_p_value)
        if not math.isfinite(p_value) or not 0.0 <= p_value <= 1.0:
            raise AnalysisValidationError("Holm p-values must lie in [0,1]")
        normalized.append((test_id, p_value))

    ordered = sorted(normalized, key=lambda item: (item[1], item[0]))
    adjusted: dict[str, float] = {}
    running = 0.0
    total = len(ordered)
    for rank, (test_id, p_value) in enumerate(ordered):
        running = max(running, min(1.0, (total - rank) * p_value))
        adjusted[test_id] = running
    return [
        {
            "test_id": test_id,
            "raw_p_value": p_value,
            "holm_adjusted_p_value": adjusted[test_id],
            "holm_reject_at_family_alpha": adjusted[test_id] <= alpha,
            "family_alpha": alpha,
        }
        for test_id, p_value in normalized
    ]


def _row_shot_audit(
    expected_row: Mapping[str, Any],
    counts: np.ndarray,
) -> dict[str, Any]:
    valid_codes = expected_row["valid_codes"]
    joint_counts = {
        outcome_key(*joint_tuple(index), valid_codes): int(count)
        for index, count in enumerate(counts)
        if int(count) > 0
    }
    leakage_shots = sum(
        int(count)
        for index, count in enumerate(counts)
        if leakage_bit(joint_tuple(index)[0], joint_tuple(index)[2], valid_codes)
    )
    decision_counts = [
        int(
            sum(
                counts[joint_index(heated, decision, final)]
                for heated in range(STATE_CARDINALITY)
                for final in range(STATE_CARDINALITY)
            )
        )
        for decision in range(DECISION_CARDINALITY)
    ]
    return {
        "row_id": expected_row["row_id"],
        "endpoint": expected_row["endpoint"],
        "backend_id": expected_row["backend_id"],
        "layout_id": expected_row["layout_id"],
        "declared_submitted_retrieved_counted_shots": int(np.sum(counts)),
        "leakage_shots": leakage_shots,
        "decision_counts": decision_counts,
        "nonzero_joint_outcomes": len(joint_counts),
        "joint_counts": joint_counts,
        "joint_counts_sha256": sha256_json(joint_counts),
    }


def run_blind_analysis(
    lock: Mapping[str, Any],
    data: Mapping[str, Any],
) -> dict[str, Any]:
    """Run the immutable blind analysis or raise before scoring."""

    validate_analysis_lock(lock, verify_code_hash=True)
    counts_by_row = validate_data_packet(lock, data)
    expected_rows = {row["row_id"]: row for row in lock["expected_rows"]}
    models = tuple(lock["models"])

    probabilities_by_row: dict[str, dict[str, np.ndarray]] = {}
    likelihoods_by_row: dict[str, dict[str, float]] = {}
    calibration_receipts: dict[str, dict[str, str]] = {}
    for row_id, expected_row in expected_rows.items():
        calibration_id = expected_row["calibration_id"]
        calibration = lock["calibrations"][calibration_id]
        calibration_receipts[calibration_id] = {
            "mode": validate_calibration(calibration_id, calibration),
            "calibration_receipt_sha256": calibration["calibration_receipt_sha256"],
            "frozen_calibration_object_sha256": sha256_json(calibration),
        }
        probabilities_by_row[row_id] = {}
        likelihoods_by_row[row_id] = {}
        for model in models:
            probabilities = _candidate_observed_probabilities(
                expected_row,
                calibration,
                model,
            )
            probabilities_by_row[row_id][model] = probabilities
            likelihoods_by_row[row_id][model] = multinomial_log_likelihood(
                counts_by_row[row_id],
                probabilities,
            )

    primary_rows = [
        row for row in lock["expected_rows"] if row["endpoint"] == PRIMARY_ENDPOINT
    ]
    primary_by_backend: dict[str, list[Mapping[str, Any]]] = {}
    for row in primary_rows:
        primary_by_backend.setdefault(row["backend_id"], []).append(row)

    pooled_scores = {
        model: float(sum(likelihoods_by_row[row["row_id"]][model] for row in primary_rows))
        for model in models
    }
    per_backend_scores: dict[str, dict[str, float]] = {}
    for backend_id, rows in primary_by_backend.items():
        per_backend_scores[backend_id] = {
            model: float(sum(likelihoods_by_row[row["row_id"]][model] for row in rows))
            for model in models
        }

    pooled_bayes_factors: dict[str, dict[str, Any]] = {}
    for null in models[1:]:
        log_bf = pooled_scores[REPAIR_MODEL] - pooled_scores[null]
        pooled_bayes_factors[null] = {
            "log_bayes_factor_repair_over_null": log_bf,
            "bayes_factor_repair_over_null": _safe_bayes_factor(log_bf),
            "passes_pooled_threshold": log_bf > math.log(POOLED_BF_THRESHOLD),
        }

    per_backend_bayes_factors: dict[str, dict[str, dict[str, Any]]] = {}
    for backend_id, scores in per_backend_scores.items():
        per_backend_bayes_factors[backend_id] = {}
        for null in models[1:]:
            log_bf = scores[REPAIR_MODEL] - scores[null]
            per_backend_bayes_factors[backend_id][null] = {
                "log_bayes_factor_repair_over_null": log_bf,
                "bayes_factor_repair_over_null": _safe_bayes_factor(log_bf),
                "passes_per_backend_threshold": log_bf
                > math.log(PER_BACKEND_BF_THRESHOLD),
            }

    global_envelope_records = [
        (
            row["row_id"],
            counts_by_row[row["row_id"]],
            probabilities_by_row[row["row_id"]][REPAIR_MODEL],
        )
        for row in primary_rows
    ]
    global_envelope = simultaneous_hoeffding_envelope(global_envelope_records)
    backend_envelopes = {
        backend_id: simultaneous_hoeffding_envelope(
            [
                (
                    row["row_id"],
                    counts_by_row[row["row_id"]],
                    probabilities_by_row[row["row_id"]][REPAIR_MODEL],
                )
                for row in rows
            ]
        )
        for backend_id, rows in primary_by_backend.items()
    }

    secondary_raw: list[tuple[str, float]] = []
    secondary_diagnostics: dict[str, dict[str, Any]] = {}
    for spec in lock["secondary_tests"]:
        model = spec.get("model", REPAIR_MODEL)
        records = [
            (
                row_id,
                counts_by_row[row_id],
                probabilities_by_row[row_id][model],
            )
            for row_id in spec["row_ids"]
        ]
        raw_p_value, statistic, cells = conservative_hoeffding_gof_p_value(records)
        secondary_raw.append((spec["test_id"], raw_p_value))
        secondary_diagnostics[spec["test_id"]] = {
            "model": model,
            "row_ids": list(spec["row_ids"]),
            "method": "union-bound Hoeffding max-cell goodness-of-fit",
            "max_scaled_cell_residual": statistic,
            "cell_count": cells,
        }
    holm_results = holm_adjust(secondary_raw)
    for result in holm_results:
        result.update(secondary_diagnostics[result["test_id"]])

    every_backend_bf_passes = all(
        result["passes_per_backend_threshold"]
        for nulls in per_backend_bayes_factors.values()
        for result in nulls.values()
    )
    every_pooled_bf_passes = all(
        result["passes_pooled_threshold"] for result in pooled_bayes_factors.values()
    )
    primary_pass = bool(
        every_backend_bf_passes and every_pooled_bf_passes and global_envelope["pass"]
    )

    decisive_null_failure = any(
        all(
            per_backend_bayes_factors[backend_id][null][
                "log_bayes_factor_repair_over_null"
            ]
            < -math.log(100.0)
            for backend_id in per_backend_bayes_factors
        )
        for null in models[1:]
    )
    both_system_envelope_failure = all(
        not envelope["pass"] for envelope in backend_envelopes.values()
    )
    kernel_failure = bool(decisive_null_failure or both_system_envelope_failure)
    if primary_pass:
        verdict = "passes_frozen_reduced_repair_kernel_gate"
    elif kernel_failure:
        verdict = "fails_frozen_reduced_repair_kernel"
    else:
        verdict = "valid_but_inconclusive"

    report_body = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "blinded": True,
        "revealed": False,
        "analysis_lock_sha256": lock["analysis_lock_sha256"],
        "analysis_code_sha256": analysis_code_sha256(),
        "blind_manifest_commitment": lock["blind_manifest_commitment"],
        "heldout_data_packet_sha256": data["data_packet_sha256"],
        "validity": {
            "fail_closed_validation_passed": True,
            "all_expected_rows_present": True,
            "all_submitted_jobs_included": True,
            "shot_conservation_passed": True,
            "postselection_detected": False,
            "calibration_frozen_before_heldout": True,
        },
        "likelihood": {
            "formula": "Multinomial(n; p_obs), p_obs(y)=sum_x C(y|x)*p_latent(x)",
            "point_hypotheses": True,
            "equal_prior_odds": True,
            "multinomial_constants_included": True,
            "probability_floor": None,
            "calibration_receipts": calibration_receipts,
            "primary_pooled_log_likelihoods": pooled_scores,
            "primary_per_backend_log_likelihoods": per_backend_scores,
        },
        "primary_endpoint": {
            "name": PRIMARY_ENDPOINT,
            "backend_count": len(primary_by_backend),
            "row_count": len(primary_rows),
            "required_nulls": list(models[1:]),
            "per_backend_bayes_factors": per_backend_bayes_factors,
            "pooled_bayes_factors": pooled_bayes_factors,
            "per_backend_threshold": PER_BACKEND_BF_THRESHOLD,
            "pooled_threshold": POOLED_BF_THRESHOLD,
            "global_99_percent_simultaneous_envelope": global_envelope,
            "per_backend_99_percent_simultaneous_envelopes": backend_envelopes,
            "label_layout_multiplicity_marginalized": True,
            "passes_gate": primary_pass,
        },
        "secondary_family": {
            "correction": "Holm",
            "family_alpha": SECONDARY_FAMILY_ALPHA,
            "does_not_select_primary": True,
            "tests": holm_results,
        },
        "shot_audit": [
            _row_shot_audit(row, counts_by_row[row["row_id"]])
            for row in lock["expected_rows"]
        ],
        "candidate_probability_hashes": {
            row_id: {
                model: sha256_json(
                    vector_to_probability_mapping(
                        probabilities,
                        expected_rows[row_id]["valid_codes"],
                        include_zeros=True,
                    )
                )
                for model, probabilities in candidates.items()
            }
            for row_id, candidates in probabilities_by_row.items()
        },
        "decision": {
            "verdict": verdict,
            "primary_pass": primary_pass,
            "kernel_failure": kernel_failure,
            "decisive_null_failure": decisive_null_failure,
            "both_system_envelope_failure": both_system_envelope_failure,
        },
        "claim_boundary": (
            "A pass supports a blinded finite self-reading repair implementation against the "
            "frozen controller nulls only. Standard quantum mechanics predicts the programmed "
            "dynamic circuit."
        ),
    }
    report = dict(report_body)
    report["blind_report_sha256"] = sha256_json(report_body)
    return report


def simulate_data_packet(
    lock: Mapping[str, Any],
    *,
    generating_model: str,
    seed: int,
    created_utc: str = "synthetic-preflight",
) -> dict[str, Any]:
    """Generate a clearly labeled synthetic packet for tests and power checks."""

    validate_analysis_lock(lock, verify_code_hash=True)
    if generating_model not in lock["models"]:
        raise AnalysisValidationError("unknown synthetic generating model")
    rng = np.random.default_rng(seed)
    rows = []
    for expected_row in lock["expected_rows"]:
        calibration = lock["calibrations"][expected_row["calibration_id"]]
        probabilities = _candidate_observed_probabilities(
            expected_row,
            calibration,
            generating_model,
        )
        shots = int(expected_row["shots"])
        counts = rng.multinomial(shots, probabilities)
        count_mapping = {
            outcome_key(*joint_tuple(index), expected_row["valid_codes"]): int(count)
            for index, count in enumerate(counts)
            if int(count) > 0
        }
        rows.append(
            {
                "row_id": expected_row["row_id"],
                "declared_shots": shots,
                "submitted_shots": shots,
                "retrieved_shots": shots,
                "excluded_shots": 0,
                "postselected": False,
                "all_outcomes_included": True,
                "counts": count_mapping,
                "synthetic_generating_model": generating_model,
            }
        )
    return seal_data_packet(
        analysis_lock_sha256=lock["analysis_lock_sha256"],
        rows=rows,
        created_utc=created_utc,
    )


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AnalysisValidationError(f"cannot read JSON {path}: {exc}") from exc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the frozen blinded Cayley repair likelihood analysis."
    )
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = run_blind_analysis(_load_json(args.lock), _load_json(args.data))
    except AnalysisValidationError as exc:
        print(f"FAIL_CLOSED: {exc}", file=sys.stderr)
        return 2
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"blind_report_sha256": report["blind_report_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
