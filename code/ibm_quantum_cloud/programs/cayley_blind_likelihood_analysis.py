#!/usr/bin/env python3
"""Frozen blinded likelihood analysis for the record-gated Cayley benchmark.

The analysis compares fixed process hypotheses using conditional likelihood
ratios.  The calibration channels are frozen from independent data, smoothed
with a positive preregistered Dirichlet pseudocount, and accompanied by frozen
sensitivity channels.  Their finite-sample uncertainty is not analytically
integrated, so this module reports conditional likelihood ratios only.

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
from itertools import product
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
STATE_PREPARATION_MODEL = "state_preparation_only"
MH_REFERENCE_MODEL = "kappa_1"
MH_POINT_MODELS = (
    MH_REFERENCE_MODEL,
    "kappa_0",
    "kappa_2",
    "mh_calibrated_noise",
)
REQUIRED_NULL_MODELS = (
    "lazy_heat",
    "delayed_record",
    "shuffled_record",
    "inverted_record",
    STATE_PREPARATION_MODEL,
    LABEL_MODEL,
    "calibrated_noise",
)
REQUIRED_MODELS = (REPAIR_MODEL, *REQUIRED_NULL_MODELS)
CAYLEY_POINT_MODELS = tuple(model for model in REQUIRED_MODELS if model != LABEL_MODEL)

PER_BACKEND_LR_THRESHOLD = 10.0
POOLED_LR_THRESHOLD = 100.0
SIMULTANEOUS_LEVEL = 0.99
SECONDARY_FAMILY_ALPHA = 0.01
CALIBRATION_CONTROL_MIN_P_VALUE = 0.01
CALIBRATION_MAX_BASIS_ERROR_FRACTION = 0.15
DEFAULT_MAX_LEAKAGE_FRACTION = 0.10
INDIVIDUAL_CIRCUIT_CATASTROPHIC_LEAKAGE_FRACTION = 0.25
CALIBRATION_UNCERTAINTY_MODE = "conditional_plugin_positive_dirichlet_sensitivity"

_OUTCOME_RE = re.compile(r"^h([0-7])\|d([01])\|f([0-7])\|l([01])$")
_QISKIT_JOINED_RE = re.compile(r"^[01]{7}$")
_QISKIT_CALIBRATION_RE = re.compile(r"^[01]{4}$")
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


def _binomial_upper_tail(trials: int, errors: int, probability: float) -> float:
    """Exact P[Binomial(trials, probability) >= errors]."""

    if trials < 0 or errors < 0 or errors > trials or not 0.0 < probability < 1.0:
        raise AnalysisValidationError("invalid binomial calibration request")
    if errors == 0:
        return 1.0
    return min(
        1.0,
        float(
            sum(
                math.comb(trials, count)
                * (probability**count)
                * ((1.0 - probability) ** (trials - count))
                for count in range(errors, trials + 1)
            )
        ),
    )


def basis_calibration_control_result(
    *,
    diagnostic_counts_by_opaque_id: Mapping[str, Mapping[str, Any]],
    control_rule: Mapping[str, Any],
    provider_job_ids: Sequence[str],
    calibration_receipt_sha256: str,
) -> dict[str, Any]:
    """Build the job/count-bound 4-bit basis-calibration validity result.

    Each prepared code receives an exact one-sided binomial test of the null
    that its wrong-code probability is no greater than the frozen maximum.
    Bonferroni correction across every prepared code is valid without assuming
    independence between code-specific tests.
    """

    _require_sha256(calibration_receipt_sha256, "calibration_receipt_sha256")
    if not isinstance(diagnostic_counts_by_opaque_id, Mapping):
        raise AnalysisValidationError("diagnostic calibration counts must be a mapping")
    expected_codes = control_rule.get("expected_basis_code_by_opaque_id")
    if not isinstance(expected_codes, Mapping) or set(expected_codes) != set(
        diagnostic_counts_by_opaque_id
    ):
        raise AnalysisValidationError("diagnostic calibration circuit set is incomplete")
    if list(control_rule.get("diagnostic_opaque_ids", ())) != sorted(expected_codes):
        raise AnalysisValidationError("diagnostic calibration order differs from frozen rule")
    maximum_error = float(control_rule.get("maximum_basis_error_fraction"))
    if maximum_error != CALIBRATION_MAX_BASIS_ERROR_FRACTION:
        raise AnalysisValidationError("basis calibration error threshold changed")
    if control_rule.get("multiple_test_correction") != "Bonferroni":
        raise AnalysisValidationError("basis calibration correction changed")
    if (
        not isinstance(provider_job_ids, Sequence)
        or isinstance(provider_job_ids, (str, bytes))
        or not provider_job_ids
        or len(set(provider_job_ids)) != len(provider_job_ids)
        or any(not isinstance(job_id, str) or not job_id for job_id in provider_job_ids)
    ):
        raise AnalysisValidationError("calibration provider job IDs are invalid")

    per_code: list[dict[str, Any]] = []
    raw_counts: dict[str, dict[str, int]] = {}
    for opaque_id in sorted(expected_codes):
        expected_code = expected_codes[opaque_id]
        if isinstance(expected_code, bool) or not isinstance(expected_code, int) or expected_code not in range(16):
            raise AnalysisValidationError("expected calibration basis code is invalid")
        supplied_counts = diagnostic_counts_by_opaque_id[opaque_id]
        if not isinstance(supplied_counts, Mapping) or not supplied_counts:
            raise AnalysisValidationError("diagnostic count table is empty")
        normalized: dict[str, int] = {}
        for key, raw_count in supplied_counts.items():
            bit_string = str(key)
            if _QISKIT_CALIBRATION_RE.fullmatch(bit_string) is None:
                raise AnalysisValidationError(
                    "basis diagnostic outcomes must be exactly four binary bits"
                )
            if isinstance(raw_count, bool) or not isinstance(raw_count, int) or raw_count < 0:
                raise AnalysisValidationError("basis diagnostic counts must be integers")
            normalized[bit_string] = raw_count
        shots = sum(normalized.values())
        if shots <= 0:
            raise AnalysisValidationError("basis diagnostic has no shots")
        expected_key = f"{expected_code:04b}"
        correct = normalized.get(expected_key, 0)
        errors = shots - correct
        raw_p_value = _binomial_upper_tail(shots, errors, maximum_error)
        per_code.append(
            {
                "opaque_id": opaque_id,
                "expected_basis_code": expected_code,
                "shots": shots,
                "correct_count": correct,
                "error_count": errors,
                "error_fraction": errors / shots,
                "raw_upper_tail_p_value": raw_p_value,
            }
        )
        raw_counts[opaque_id] = normalized
    family_size = len(per_code)
    adjusted_p_value = min(
        1.0,
        family_size * min(record["raw_upper_tail_p_value"] for record in per_code),
    )
    summary = {
        "method": "exact one-sided binomial upper-tail with Bonferroni correction",
        "maximum_basis_error_fraction": maximum_error,
        "family_size": family_size,
        "per_code": per_code,
    }
    return {
        "calibration_receipt_sha256": calibration_receipt_sha256,
        "diagnostic_counts_sha256": sha256_json(raw_counts),
        "diagnostic_summary_sha256": sha256_json(summary),
        "diagnostic_opaque_ids": sorted(expected_codes),
        "provider_job_ids": list(provider_job_ids),
        "all_diagnostic_jobs_complete": True,
        "all_diagnostic_shots_included": True,
        "postselected": False,
        "gof_p_value": adjusted_p_value,
        "minimum_count_per_prepared_state": min(record["shots"] for record in per_code),
        "summary": summary,
    }


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


def state_preparation_only_joint_null(
    repair_probabilities: Mapping[str, Any],
    valid_codes: Sequence[int],
) -> dict[str, float]:
    """Match marginals while destroying heated/decision-to-final dependence.

    The null retains the repair hypothesis' ``p(h,d)`` and ``p(f)`` exactly but
    factorizes them as ``p(h,d,f)=p(h,d)p(f)``.  It therefore cannot gain
    evidence merely by matching the final-state histogram.
    """

    repair = probability_mapping_to_vector(repair_probabilities, valid_codes).reshape(
        STATE_CARDINALITY,
        DECISION_CARDINALITY,
        STATE_CARDINALITY,
    )
    heated_decision = np.sum(repair, axis=2)
    final = np.sum(repair, axis=(0, 1))
    factorized = heated_decision[:, :, np.newaxis] * final[np.newaxis, np.newaxis, :]
    return vector_to_probability_mapping(
        factorized.reshape(JOINT_CARDINALITY),
        valid_codes,
    )


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
        calibration.get("calibration_protocol_sha256"),
        f"calibrations.{calibration_id}.calibration_protocol_sha256",
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

    control_rule = calibration.get("control_rule")
    if not isinstance(control_rule, Mapping):
        raise AnalysisValidationError(f"calibration {calibration_id!r} lacks a control rule")
    if control_rule.get("minimum_gof_p_value") != CALIBRATION_CONTROL_MIN_P_VALUE:
        raise AnalysisValidationError("calibration control p-value threshold changed")
    required_count = control_rule.get("minimum_count_per_prepared_state")
    diagnostic_ids = control_rule.get("diagnostic_opaque_ids")
    expected_codes = control_rule.get("expected_basis_code_by_opaque_id")
    if (
        isinstance(required_count, bool)
        or not isinstance(required_count, int)
        or required_count <= 0
        or not isinstance(diagnostic_ids, list)
        or not diagnostic_ids
        or len(set(diagnostic_ids)) != len(diagnostic_ids)
        or any(not isinstance(value, str) or not value for value in diagnostic_ids)
        or diagnostic_ids != sorted(diagnostic_ids)
        or not isinstance(expected_codes, Mapping)
        or set(expected_codes) != set(diagnostic_ids)
        or any(
            isinstance(code, bool) or not isinstance(code, int) or code not in range(16)
            for code in expected_codes.values()
        )
        or len(set(expected_codes.values())) != len(expected_codes)
        or control_rule.get("maximum_basis_error_fraction")
        != CALIBRATION_MAX_BASIS_ERROR_FRACTION
        or control_rule.get("multiple_test_correction") != "Bonferroni"
    ):
        raise AnalysisValidationError("calibration control rule is invalid")

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
    diagnostic_opaque_ids: Sequence[str] = ("synthetic-calibration",),
    diagnostic_expected_codes: Mapping[str, int] | None = None,
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
    diagnostic_ids = sorted(str(value) for value in diagnostic_opaque_ids)
    if diagnostic_expected_codes is None:
        diagnostic_expected_codes = {
            opaque_id: index for index, opaque_id in enumerate(diagnostic_ids)
        }
    return {
        "frozen_before_heldout": True,
        "uses_heldout_outputs": False,
        "calibration_protocol_sha256": receipt_sha256,
        "calibration_derivation_sha256": derivation_sha256,
        "uncertainty_mode": CALIBRATION_UNCERTAINTY_MODE,
        "dirichlet_pseudocount": dirichlet_pseudocount,
        "max_age_seconds": 86400.0,
        "control_rule": {
            "minimum_gof_p_value": CALIBRATION_CONTROL_MIN_P_VALUE,
            "minimum_count_per_prepared_state": 128,
            "diagnostic_opaque_ids": diagnostic_ids,
            "expected_basis_code_by_opaque_id": dict(diagnostic_expected_codes),
            "maximum_basis_error_fraction": CALIBRATION_MAX_BASIS_ERROR_FRACTION,
            "multiple_test_correction": "Bonferroni",
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
    minimum_required_count: int = 128,
    max_age_seconds: float = 86400.0,
    diagnostic_opaque_ids: Sequence[str] = ("synthetic-calibration",),
    diagnostic_expected_codes: Mapping[str, int] | None = None,
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
    diagnostic_ids = sorted(str(value) for value in diagnostic_opaque_ids)
    if diagnostic_expected_codes is None:
        diagnostic_expected_codes = {
            opaque_id: index for index, opaque_id in enumerate(diagnostic_ids)
        }
    packet = {
        "frozen_before_heldout": True,
        "uses_heldout_outputs": False,
        "calibration_protocol_sha256": receipt_sha256,
        "calibration_derivation_sha256": derivation_sha256,
        "uncertainty_mode": CALIBRATION_UNCERTAINTY_MODE,
        "dirichlet_pseudocount": dirichlet_pseudocount,
        "max_age_seconds": max_age_seconds,
        "control_rule": {
            "minimum_gof_p_value": CALIBRATION_CONTROL_MIN_P_VALUE,
            "minimum_count_per_prepared_state": minimum_required_count,
            "diagnostic_opaque_ids": diagnostic_ids,
            "expected_basis_code_by_opaque_id": dict(diagnostic_expected_codes),
            "maximum_basis_error_fraction": CALIBRATION_MAX_BASIS_ERROR_FRACTION,
            "multiple_test_correction": "Bonferroni",
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
        STATE_PREPARATION_MODEL: {
            "kind": "final_marginal_matched_joint_factorization_null",
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
    logical_circuit_sha256: str,
) -> dict[str, dict[str, str]]:
    """Bind every frozen point-hypothesis table to data and derivation hashes."""

    if set(candidate_probabilities) != set(derivation_sha256_by_model):
        raise AnalysisValidationError("candidate provenance model set is incomplete")
    logical_digest = _require_sha256(
        logical_circuit_sha256,
        "candidate logical_circuit_sha256",
    )
    result: dict[str, dict[str, str]] = {}
    for model, probabilities in candidate_probabilities.items():
        result[model] = {
            "probability_table_sha256": sha256_json(probabilities),
            "derivation_sha256": _require_sha256(
                derivation_sha256_by_model[model],
                f"candidate derivation for {model}",
            ),
            "logical_circuit_sha256": logical_digest,
        }
    return result


def _validate_candidate_provenance(
    candidates: Mapping[str, Any],
    provenance: Any,
    *,
    field: str,
    logical_circuit_sha256: str,
) -> None:
    if not isinstance(provenance, Mapping) or set(provenance) != set(candidates):
        raise AnalysisValidationError(f"{field} does not match the candidate table set")
    for model, probabilities in candidates.items():
        record = provenance[model]
        if not isinstance(record, Mapping):
            raise AnalysisValidationError(f"{field}.{model} must be a mapping")
        if set(record) != {
            "probability_table_sha256",
            "derivation_sha256",
            "logical_circuit_sha256",
        }:
            raise AnalysisValidationError(f"{field}.{model} has an unexpected schema")
        _require_sha256(record.get("derivation_sha256"), f"{field}.{model}.derivation_sha256")
        bound_logical_digest = _require_sha256(
            record.get("logical_circuit_sha256"),
            f"{field}.{model}.logical_circuit_sha256",
        )
        if bound_logical_digest != logical_circuit_sha256:
            raise AnalysisValidationError(
                f"{field}.{model} binds a different logical circuit"
            )
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
    if "logical_qpy_sha256" in row:
        raise AnalysisValidationError(
            f"row {row_id!r} uses removed logical_qpy_sha256 identity"
        )
    family = row.get("family")
    if family not in ("cayley", "mh"):
        raise AnalysisValidationError(f"row {row_id!r} has an unsupported analysis family")
    for field in (
        "endpoint",
        "backend_role",
        "backend_name",
        "layout_id",
        "calibration_id",
    ):
        if not isinstance(row.get(field), str) or not row[field]:
            raise AnalysisValidationError(f"row {row_id!r} needs {field}")
    logical_circuit_sha256 = _require_sha256(
        row.get("logical_circuit_sha256"),
        f"expected_rows.{row_id}.logical_circuit_sha256",
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
    if (
        not math.isfinite(leakage_limit)
        or leakage_limit != INDIVIDUAL_CIRCUIT_CATASTROPHIC_LEAKAGE_FRACTION
    ):
        raise AnalysisValidationError(f"row {row_id!r} leakage limit is invalid")
    valid_codes = _validate_valid_codes(row.get("valid_codes"), row_id)
    candidates = row.get("candidate_probabilities")
    if not isinstance(candidates, Mapping):
        raise AnalysisValidationError(f"row {row_id!r} lacks candidate probabilities")
    expected_candidate_models = (
        CAYLEY_POINT_MODELS if family == "cayley" else MH_POINT_MODELS
    )
    if set(candidates) != set(expected_candidate_models):
        raise AnalysisValidationError(
            f"row {row_id!r} candidate set does not match its frozen family model set"
        )
    for model in expected_candidate_models:
        probability_mapping_to_vector(candidates[model], valid_codes)
    if family == "cayley":
        stratum_id = row.get("state_preparation_stratum_id")
        if not isinstance(stratum_id, str) or not stratum_id:
            raise AnalysisValidationError(
                f"row {row_id!r} lacks an opaque state-preparation stratum"
            )
    if family == "mh":
        identifiability_role = row.get("identifiability_role")
        if identifiability_role not in ("identifiable", "abelian_negative_control"):
            raise AnalysisValidationError(f"row {row_id!r} has an invalid MH role")
        kappa_hashes = {
            sha256_json(candidates[model])
            for model in ("kappa_1", "kappa_0", "kappa_2")
        }
        if identifiability_role == "abelian_negative_control" and len(kappa_hashes) != 1:
            raise AnalysisValidationError(
                f"row {row_id!r} abelian control must be exactly nonidentifiable"
            )
    _validate_candidate_provenance(
        candidates,
        row.get("candidate_provenance"),
        field=f"expected_rows.{row_id}.candidate_provenance",
        logical_circuit_sha256=logical_circuit_sha256,
    )


def _validate_secondary_specs(
    specs: Any,
    rows_by_id: Mapping[str, Mapping[str, Any]],
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
        if len(set(selected_rows)) != len(selected_rows) or not set(selected_rows) <= set(
            rows_by_id
        ):
            raise AnalysisValidationError(f"secondary test {test_id!r} has invalid row IDs")
        selected_model = spec.get("model", REPAIR_MODEL)
        if selected_model == LABEL_MODEL or any(
            selected_model not in rows_by_id[row_id]["candidate_probabilities"]
            for row_id in selected_rows
        ):
            raise AnalysisValidationError(
                f"secondary test {test_id!r} model does not match every selected row family"
            )


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
    analysis_row_ids = {
        row_id for row_id, row in rows_by_id.items() if row["endpoint"] == PRIMARY_ENDPOINT
    }
    components = label_model.get("components")
    if not isinstance(components, list) or len(components) < 2:
        raise AnalysisValidationError("label/layout model needs at least two frozen components")
    component_ids: set[str] = set()
    component_table_hashes: set[str] = set()
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
        if not isinstance(row_probabilities, Mapping) or set(row_probabilities) != analysis_row_ids:
            raise AnalysisValidationError(
                f"label/layout component {component_id!r} must cover every primary circuit"
            )
        if not isinstance(row_provenance, Mapping) or set(row_provenance) != analysis_row_ids:
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
        table_set_hash = sha256_json(row_probabilities)
        if table_set_hash in component_table_hashes:
            raise AnalysisValidationError("label/layout components contain duplicate table sets")
        component_table_hashes.add(table_set_hash)
    if not math.isclose(weight_total, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise AnalysisValidationError("label/layout component prior weights must sum to one")
    if label_model.get("reference_component_id") not in component_ids:
        raise AnalysisValidationError("label/layout reference component is absent")


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
        "pooled_backend_family_max_leakage_fraction": DEFAULT_MAX_LEAKAGE_FRACTION,
        "individual_circuit_catastrophic_leakage_fraction": (
            INDIVIDUAL_CIRCUIT_CATASTROPHIC_LEAKAGE_FRACTION
        ),
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
        "evidence_measure": "conditional_likelihood_ratio_only",
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
            if row["family"] != "cayley":
                raise AnalysisValidationError("the primary S3 endpoint must contain Cayley rows")
            primary_backends.add(row["backend_name"])
    if {row["calibration_id"] for row in rows} != set(calibrations):
        raise AnalysisValidationError("calibration set must exactly match referenced circuit rows")
    if len(primary_backends) < 2:
        raise AnalysisValidationError("primary S3 endpoint needs two independent backends")
    mh_by_endpoint: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        if row["family"] == "mh":
            mh_by_endpoint.setdefault(row["endpoint"], []).append(row)
    for endpoint, endpoint_rows in mh_by_endpoint.items():
        roles = {row["identifiability_role"] for row in endpoint_rows}
        if len(roles) != 1:
            raise AnalysisValidationError(f"MH endpoint {endpoint!r} mixes identifiability roles")
        if next(iter(roles)) == "identifiable" and not any(
            len(
                {
                    sha256_json(row["candidate_probabilities"][model])
                    for model in ("kappa_1", "kappa_0", "kappa_2")
                }
            )
            > 1
            for row in endpoint_rows
        ):
            raise AnalysisValidationError(
                f"MH endpoint {endpoint!r} is declared identifiable but has zero separation"
            )
    state_preparation_strata: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        if row["family"] == "cayley":
            state_preparation_strata.setdefault(
                row["state_preparation_stratum_id"], []
            ).append(row)
    for stratum_id, stratum_rows in state_preparation_strata.items():
        if len({row["backend_role"] for row in stratum_rows}) != 1:
            raise AnalysisValidationError(
                f"state-preparation stratum {stratum_id!r} crosses backend roles"
            )
        valid_code_sets = {tuple(row["valid_codes"]) for row in stratum_rows}
        if len(valid_code_sets) != 1:
            raise AnalysisValidationError(
                f"state-preparation stratum {stratum_id!r} crosses physical encodings"
            )
        valid_codes = stratum_rows[0]["valid_codes"]
        total_shots = sum(int(row["shots"]) for row in stratum_rows)
        aggregate_repair = sum(
            int(row["shots"])
            * probability_mapping_to_vector(
                row["candidate_probabilities"][REPAIR_MODEL], valid_codes
            )
            for row in stratum_rows
        ) / total_shots
        expected_null = probability_mapping_to_vector(
            state_preparation_only_joint_null(
                vector_to_probability_mapping(aggregate_repair, valid_codes),
                valid_codes,
            ),
            valid_codes,
        )
        for row in stratum_rows:
            supplied = probability_mapping_to_vector(
                row["candidate_probabilities"][STATE_PREPARATION_MODEL],
                valid_codes,
            )
            if not np.allclose(supplied, expected_null, atol=1e-12, rtol=0.0):
                raise AnalysisValidationError(
                    f"state-preparation stratum {stratum_id!r} does not preserve/factorize marginals"
                )
    if len({row["shots"] for row in rows}) != 1:
        raise AnalysisValidationError("all locked benchmark circuits must use the same shot count")
    rows_by_id = {row["row_id"]: row for row in rows}
    diagnostic_ids = [
        opaque_id
        for calibration in calibrations.values()
        for opaque_id in calibration["control_rule"]["diagnostic_opaque_ids"]
    ]
    if len(set(diagnostic_ids)) != len(diagnostic_ids):
        raise AnalysisValidationError("diagnostic calibration circuit IDs must be globally unique")
    if set(diagnostic_ids) & row_ids:
        raise AnalysisValidationError("diagnostic and dynamic circuit IDs must be disjoint")
    expected_coverage = {
        "dynamic_analysis_opaque_ids": sorted(row_ids),
        "diagnostic_calibration_opaque_ids": sorted(diagnostic_ids),
        "all_catalog_circuits_classified": True,
    }
    if lock.get("catalog_coverage") != expected_coverage:
        raise AnalysisValidationError("catalog coverage does not match dynamic and diagnostic rows")
    _validate_label_layout_model(lock.get("label_layout_model"), rows_by_id)
    _validate_secondary_specs(lock.get("secondary_tests"), rows_by_id)


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
            "pooled_backend_family_max_leakage_fraction": DEFAULT_MAX_LEAKAGE_FRACTION,
            "individual_circuit_catastrophic_leakage_fraction": (
                INDIVIDUAL_CIRCUIT_CATASTROPHIC_LEAKAGE_FRACTION
            ),
        },
        "execution_contract": {
            "row_granularity": "one_locked_row_per_opaque_circuit",
            "qiskit_joined_bit_order": "final[3]|decision[1]|heated[3]",
            "calibration_uncertainty": CALIBRATION_UNCERTAINTY_MODE,
            "evidence_measure": "conditional_likelihood_ratio_only",
        },
        "calibrations": _json_copy(calibrations),
        "expected_rows": _json_copy(list(expected_rows)),
        "catalog_coverage": {
            "dynamic_analysis_opaque_ids": sorted(
                str(row["opaque_id"]) for row in expected_rows
            ),
            "diagnostic_calibration_opaque_ids": sorted(
                str(opaque_id)
                for calibration in calibrations.values()
                for opaque_id in calibration["control_rule"]["diagnostic_opaque_ids"]
            ),
            "all_catalog_circuits_classified": True,
        },
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
    if "logical_qpy_sha256" in data_row:
        raise AnalysisValidationError(
            f"data row {row_id!r} uses removed logical_qpy_sha256 identity"
        )
    for field in (
        "opaque_id",
        "logical_circuit_sha256",
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
    calibration_results: Mapping[str, Mapping[str, Any]],
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
        "calibration_results": _json_copy(calibration_results),
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
    calibration_results = data.get("calibration_results")
    if not isinstance(calibration_results, Mapping) or set(calibration_results) != set(
        lock["calibrations"]
    ):
        raise AnalysisValidationError("held-out packet calibration-result set is incomplete")
    for calibration_id, calibration in lock["calibrations"].items():
        result = calibration_results[calibration_id]
        if not isinstance(result, Mapping):
            raise AnalysisValidationError("calibration result must be a mapping")
        _require_sha256(
            result.get("calibration_receipt_sha256"),
            f"calibration_results.{calibration_id}.calibration_receipt_sha256",
        )
        _require_sha256(
            result.get("diagnostic_counts_sha256"),
            f"calibration_results.{calibration_id}.diagnostic_counts_sha256",
        )
        claimed_summary_hash = _require_sha256(
            result.get("diagnostic_summary_sha256"),
            f"calibration_results.{calibration_id}.diagnostic_summary_sha256",
        )
        summary = result.get("summary")
        if not isinstance(summary, Mapping) or sha256_json(summary) != claimed_summary_hash:
            raise AnalysisValidationError("calibration diagnostic summary hash mismatch")
        if result.get("diagnostic_opaque_ids") != calibration["control_rule"][
            "diagnostic_opaque_ids"
        ]:
            raise AnalysisValidationError("calibration result covers different diagnostic circuits")
        provider_job_ids = result.get("provider_job_ids")
        if (
            not isinstance(provider_job_ids, list)
            or not provider_job_ids
            or len(set(provider_job_ids)) != len(provider_job_ids)
            or any(not isinstance(job_id, str) or not job_id for job_id in provider_job_ids)
        ):
            raise AnalysisValidationError("calibration result lacks provider job IDs")
        if result.get("all_diagnostic_jobs_complete") is not True:
            raise AnalysisValidationError("calibration diagnostic jobs are incomplete")
        if result.get("all_diagnostic_shots_included") is not True:
            raise AnalysisValidationError("calibration diagnostic shots are incomplete")
        if result.get("postselected") is not False:
            raise AnalysisValidationError("calibration result is postselected")
        try:
            gof_p_value = float(result.get("gof_p_value"))
        except (TypeError, ValueError) as exc:
            raise AnalysisValidationError("calibration result p-value is invalid") from exc
        minimum_count = result.get("minimum_count_per_prepared_state")
        if (
            not math.isfinite(gof_p_value)
            or not 0.0 <= gof_p_value <= 1.0
            or isinstance(minimum_count, bool)
            or not isinstance(minimum_count, int)
            or minimum_count < 0
        ):
            raise AnalysisValidationError("calibration result validation fields are invalid")
        per_code = summary.get("per_code")
        if not isinstance(per_code, list) or len(per_code) != len(
            calibration["control_rule"]["diagnostic_opaque_ids"]
        ):
            raise AnalysisValidationError("calibration diagnostic summary is incomplete")
        recomputed_raw_p_values: list[float] = []
        recomputed_minimum_count: int | None = None
        for record in per_code:
            if not isinstance(record, Mapping):
                raise AnalysisValidationError("calibration per-code summary is malformed")
            shots = record.get("shots")
            errors = record.get("error_count")
            if (
                isinstance(shots, bool)
                or not isinstance(shots, int)
                or shots <= 0
                or isinstance(errors, bool)
                or not isinstance(errors, int)
                or errors < 0
                or errors > shots
            ):
                raise AnalysisValidationError("calibration per-code counts are invalid")
            recomputed_raw_p_values.append(
                _binomial_upper_tail(
                    shots,
                    errors,
                    calibration["control_rule"]["maximum_basis_error_fraction"],
                )
            )
            recomputed_minimum_count = (
                shots
                if recomputed_minimum_count is None
                else min(recomputed_minimum_count, shots)
            )
        recomputed_gof = min(
            1.0,
            len(per_code) * min(recomputed_raw_p_values),
        )
        if (
            not math.isclose(gof_p_value, recomputed_gof, rel_tol=0.0, abs_tol=1e-12)
            or minimum_count != recomputed_minimum_count
        ):
            raise AnalysisValidationError("calibration result does not match frozen exact test")
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


def _safe_likelihood_ratio(log_likelihood_ratio: float) -> float | str:
    if math.isnan(log_likelihood_ratio):
        raise AnalysisValidationError("likelihood-ratio comparison produced NaN")
    if log_likelihood_ratio > 709.0:
        return "Infinity"
    if log_likelihood_ratio < -745.0:
        return 0.0
    return float(math.exp(log_likelihood_ratio))


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        raise AnalysisValidationError("cannot marginalize an empty likelihood family")
    maximum = max(float(value) for value in values)
    if not math.isfinite(maximum):
        raise AnalysisValidationError("likelihood marginal contains a non-finite score")
    return maximum + math.log(sum(math.exp(float(value) - maximum) for value in values))


def _calibration_variants(calibration: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    return [calibration, *calibration["sensitivity_channels"]]


def _row_table_log_likelihood(
    expected_row: Mapping[str, Any],
    counts: np.ndarray,
    latent_table: Mapping[str, Any],
    calibration: Mapping[str, Any],
) -> float:
    latent = probability_mapping_to_vector(latent_table, expected_row["valid_codes"])
    observed = convolve_calibration(latent, calibration)
    return multinomial_log_likelihood(counts, observed)


def _point_log_likelihood_ratio_sensitivity_bounds(
    rows: Sequence[Mapping[str, Any]],
    counts_by_row: Mapping[str, np.ndarray],
    calibrations: Mapping[str, Mapping[str, Any]],
    reference_model: str,
    null_model: str,
) -> tuple[float, float]:
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["calibration_id"], []).append(row)
    lower = 0.0
    upper = 0.0
    for calibration_id, selected_rows in grouped.items():
        differences = []
        for channel in _calibration_variants(calibrations[calibration_id]):
            reference_score = sum(
                _row_table_log_likelihood(
                    row,
                    counts_by_row[row["row_id"]],
                    row["candidate_probabilities"][reference_model],
                    channel,
                )
                for row in selected_rows
            )
            null_score = sum(
                _row_table_log_likelihood(
                    row,
                    counts_by_row[row["row_id"]],
                    row["candidate_probabilities"][null_model],
                    channel,
                )
                for row in selected_rows
            )
            differences.append(float(reference_score - null_score))
        lower += min(differences)
        upper += max(differences)
    return float(lower), float(upper)


def _label_component_log_scores(
    rows: Sequence[Mapping[str, Any]],
    counts_by_row: Mapping[str, np.ndarray],
    calibrations: Mapping[str, Mapping[str, Any]],
    label_model: Mapping[str, Any],
) -> dict[str, float]:
    scores: dict[str, float] = {}
    for component in label_model["components"]:
        scores[component["component_id"]] = float(
            sum(
                _row_table_log_likelihood(
                    row,
                    counts_by_row[row["row_id"]],
                    component["row_probabilities"][row["row_id"]],
                    calibrations[row["calibration_id"]],
                )
                for row in rows
            )
        )
    return scores


def _label_log_marginal(
    component_scores: Mapping[str, float],
    label_model: Mapping[str, Any],
) -> float:
    return _logsumexp(
        [
            math.log(float(component["prior_weight"]))
            + float(component_scores[component["component_id"]])
            for component in label_model["components"]
        ]
    )


def _label_log_likelihood_ratio_sensitivity_bounds(
    rows: Sequence[Mapping[str, Any]],
    counts_by_row: Mapping[str, np.ndarray],
    calibrations: Mapping[str, Mapping[str, Any]],
    label_model: Mapping[str, Any],
    reference_model: str,
) -> tuple[float, float]:
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(row["calibration_id"], []).append(row)
    calibration_ids = sorted(grouped)
    channel_families = [
        _calibration_variants(calibrations[calibration_id])
        for calibration_id in calibration_ids
    ]
    combination_count = math.prod(len(family) for family in channel_families)
    if combination_count > 64:
        raise AnalysisValidationError(
            "label sensitivity grid exceeds the frozen exact Cartesian-product limit"
        )
    paired_log_likelihood_ratios: list[float] = []
    for selected_channels in product(*channel_families):
        reference_score = 0.0
        component_scores = {
            component["component_id"]: 0.0 for component in label_model["components"]
        }
        for calibration_id, channel in zip(calibration_ids, selected_channels):
            for row in grouped[calibration_id]:
                reference_score += _row_table_log_likelihood(
                    row,
                    counts_by_row[row["row_id"]],
                    row["candidate_probabilities"][reference_model],
                    channel,
                )
                for component in label_model["components"]:
                    component_scores[component["component_id"]] += _row_table_log_likelihood(
                        row,
                        counts_by_row[row["row_id"]],
                        component["row_probabilities"][row["row_id"]],
                        channel,
                    )
        label_marginal = _label_log_marginal(component_scores, label_model)
        paired_log_likelihood_ratios.append(reference_score - label_marginal)
    return (
        float(min(paired_log_likelihood_ratios)),
        float(max(paired_log_likelihood_ratios)),
    )


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
        "opaque_id": expected_row["opaque_id"],
        "logical_circuit_sha256": expected_row["logical_circuit_sha256"],
        "endpoint": expected_row["endpoint"],
        "backend_role": expected_row["backend_role"],
        "backend_name": expected_row["backend_name"],
        "layout_id": expected_row["layout_id"],
        "physical_layout": list(expected_row["physical_layout"]),
        "declared_submitted_retrieved_counted_shots": int(np.sum(counts)),
        "leakage_shots": leakage_shots,
        "leakage_fraction": leakage_shots / int(np.sum(counts)),
        "max_leakage_fraction": float(expected_row["max_leakage_fraction"]),
        "leakage_gate_pass": leakage_shots / int(np.sum(counts))
        <= float(expected_row["max_leakage_fraction"]),
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
    point_models = tuple(model for model in models if model != LABEL_MODEL)

    probabilities_by_row: dict[str, dict[str, np.ndarray]] = {}
    likelihoods_by_row: dict[str, dict[str, float]] = {}
    calibration_receipts: dict[str, dict[str, Any]] = {}
    for row_id, expected_row in expected_rows.items():
        calibration_id = expected_row["calibration_id"]
        calibration = lock["calibrations"][calibration_id]
        calibration_receipts[calibration_id] = {
            "mode": validate_calibration(calibration_id, calibration),
            "calibration_protocol_sha256": calibration["calibration_protocol_sha256"],
            "calibration_derivation_sha256": calibration[
                "calibration_derivation_sha256"
            ],
            "frozen_calibration_object_sha256": sha256_json(calibration),
            "uncertainty_mode": calibration["uncertainty_mode"],
            "dirichlet_pseudocount": float(calibration["dirichlet_pseudocount"]),
            "primary_eligible": calibration["primary_eligible"],
            "control_rule": _json_copy(calibration["control_rule"]),
            "calibration_result": _json_copy(data["calibration_results"][calibration_id]),
        }
        probabilities_by_row[row_id] = {}
        likelihoods_by_row[row_id] = {}
        for model in expected_row["candidate_probabilities"]:
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
        primary_by_backend.setdefault(row["backend_name"], []).append(row)

    pooled_scores = {
        model: float(sum(likelihoods_by_row[row["row_id"]][model] for row in primary_rows))
        for model in point_models
    }
    pooled_label_component_scores = _label_component_log_scores(
        primary_rows,
        counts_by_row,
        lock["calibrations"],
        lock["label_layout_model"],
    )
    pooled_scores[LABEL_MODEL] = _label_log_marginal(
        pooled_label_component_scores,
        lock["label_layout_model"],
    )
    per_backend_scores: dict[str, dict[str, float]] = {}
    per_backend_label_component_scores: dict[str, dict[str, float]] = {}
    for backend_name, rows in primary_by_backend.items():
        per_backend_scores[backend_name] = {
            model: float(sum(likelihoods_by_row[row["row_id"]][model] for row in rows))
            for model in point_models
        }
        component_scores = _label_component_log_scores(
            rows,
            counts_by_row,
            lock["calibrations"],
            lock["label_layout_model"],
        )
        per_backend_label_component_scores[backend_name] = component_scores
        per_backend_scores[backend_name][LABEL_MODEL] = _label_log_marginal(
            component_scores,
            lock["label_layout_model"],
        )

    pooled_likelihood_ratios: dict[str, dict[str, Any]] = {}
    for null in models[1:]:
        log_lr = pooled_scores[REPAIR_MODEL] - pooled_scores[null]
        if null == LABEL_MODEL:
            sensitivity_bounds = _label_log_likelihood_ratio_sensitivity_bounds(
                primary_rows,
                counts_by_row,
                lock["calibrations"],
                lock["label_layout_model"],
                REPAIR_MODEL,
            )
        else:
            sensitivity_bounds = _point_log_likelihood_ratio_sensitivity_bounds(
                primary_rows,
                counts_by_row,
                lock["calibrations"],
                REPAIR_MODEL,
                null,
            )
        pooled_likelihood_ratios[null] = {
            "conditional_log_likelihood_ratio_repair_over_null": log_lr,
            "conditional_likelihood_ratio_repair_over_null": _safe_likelihood_ratio(log_lr),
            "sensitivity_log_likelihood_ratio_bounds": list(sensitivity_bounds),
            "passes_pooled_threshold": bool(
                log_lr > math.log(POOLED_LR_THRESHOLD)
                and sensitivity_bounds[0] > math.log(POOLED_LR_THRESHOLD)
            ),
        }

    per_backend_likelihood_ratios: dict[str, dict[str, dict[str, Any]]] = {}
    for backend_name, scores in per_backend_scores.items():
        per_backend_likelihood_ratios[backend_name] = {}
        for null in models[1:]:
            log_lr = scores[REPAIR_MODEL] - scores[null]
            backend_rows = primary_by_backend[backend_name]
            if null == LABEL_MODEL:
                sensitivity_bounds = _label_log_likelihood_ratio_sensitivity_bounds(
                    backend_rows,
                    counts_by_row,
                    lock["calibrations"],
                    lock["label_layout_model"],
                    REPAIR_MODEL,
                )
            else:
                sensitivity_bounds = _point_log_likelihood_ratio_sensitivity_bounds(
                    backend_rows,
                    counts_by_row,
                    lock["calibrations"],
                    REPAIR_MODEL,
                    null,
                )
            per_backend_likelihood_ratios[backend_name][null] = {
                "conditional_log_likelihood_ratio_repair_over_null": log_lr,
                "conditional_likelihood_ratio_repair_over_null": _safe_likelihood_ratio(log_lr),
                "sensitivity_log_likelihood_ratio_bounds": list(sensitivity_bounds),
                "passes_per_backend_threshold": bool(
                    log_lr > math.log(PER_BACKEND_LR_THRESHOLD)
                    and sensitivity_bounds[0] > math.log(PER_BACKEND_LR_THRESHOLD)
                ),
            }

    mh_rows_by_endpoint: dict[str, list[Mapping[str, Any]]] = {}
    for row in lock["expected_rows"]:
        if row["family"] == "mh":
            mh_rows_by_endpoint.setdefault(row["endpoint"], []).append(row)
    mh_secondary_families: dict[str, dict[str, Any]] = {}
    for endpoint, rows in mh_rows_by_endpoint.items():
        scores = {
            model: float(
                sum(likelihoods_by_row[row["row_id"]][model] for row in rows)
            )
            for model in MH_POINT_MODELS
        }
        ratios: dict[str, dict[str, Any]] = {}
        for null in MH_POINT_MODELS[1:]:
            log_lr = scores[MH_REFERENCE_MODEL] - scores[null]
            bounds = _point_log_likelihood_ratio_sensitivity_bounds(
                rows,
                counts_by_row,
                lock["calibrations"],
                MH_REFERENCE_MODEL,
                null,
            )
            ratios[null] = {
                "conditional_log_likelihood_ratio_kappa_1_over_null": log_lr,
                "conditional_likelihood_ratio_kappa_1_over_null": _safe_likelihood_ratio(
                    log_lr
                ),
                "sensitivity_log_likelihood_ratio_bounds": list(bounds),
            }
        role = rows[0]["identifiability_role"]
        mh_secondary_families[endpoint] = {
            "identifiability_role": role,
            "row_count": len(rows),
            "log_likelihoods": scores,
            "conditional_likelihood_ratios": ratios,
            "negative_control_exact_zero_check": bool(
                role != "abelian_negative_control"
                or all(
                    math.isclose(
                        ratios[null][
                            "conditional_log_likelihood_ratio_kappa_1_over_null"
                        ],
                        0.0,
                        rel_tol=0.0,
                        abs_tol=1e-12,
                    )
                    for null in ("kappa_0", "kappa_2")
                )
            ),
            "does_not_select_primary": True,
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
        backend_name: simultaneous_hoeffding_envelope(
            [
                (
                    row["row_id"],
                    counts_by_row[row["row_id"]],
                    probabilities_by_row[row["row_id"]][REPAIR_MODEL],
                )
                for row in rows
            ]
        )
        for backend_name, rows in primary_by_backend.items()
    }

    secondary_raw: list[tuple[str, float]] = []
    secondary_diagnostics: dict[str, dict[str, Any]] = {}
    for spec in lock["secondary_tests"]:
        model = spec.get("model", REPAIR_MODEL)
        if model == LABEL_MODEL:
            raise AnalysisValidationError("secondary tests cannot select the global label marginal")
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

    shot_audit = [
        _row_shot_audit(row, counts_by_row[row["row_id"]])
        for row in lock["expected_rows"]
    ]
    individual_catastrophic_leakage_pass = all(
        row["leakage_gate_pass"] for row in shot_audit
    )
    leakage_groups: dict[tuple[str, str], dict[str, int]] = {}
    for row in shot_audit:
        key = (row["backend_name"], expected_rows[row["row_id"]]["family"])
        group = leakage_groups.setdefault(key, {"shots": 0, "leakage_shots": 0})
        group["shots"] += row["declared_submitted_retrieved_counted_shots"]
        group["leakage_shots"] += row["leakage_shots"]
    pooled_leakage_gates = {
        f"{backend_name}|{family}": {
            "backend_name": backend_name,
            "family": family,
            "shots": values["shots"],
            "leakage_shots": values["leakage_shots"],
            "leakage_fraction": values["leakage_shots"] / values["shots"],
            "maximum_fraction": DEFAULT_MAX_LEAKAGE_FRACTION,
            "pass": values["leakage_shots"] / values["shots"]
            <= DEFAULT_MAX_LEAKAGE_FRACTION,
        }
        for (backend_name, family), values in sorted(leakage_groups.items())
    }
    pooled_leakage_pass = all(details["pass"] for details in pooled_leakage_gates.values())
    leakage_gate_pass = bool(
        individual_catastrophic_leakage_pass and pooled_leakage_pass
    )
    data_rows_by_id = {row["row_id"]: row for row in data["rows"]}
    calibration_gate_details: dict[str, dict[str, Any]] = {}
    for calibration_id, calibration in lock["calibrations"].items():
        selected_rows = [
            row for row in lock["expected_rows"] if row["calibration_id"] == calibration_id
        ]
        maximum_observed_age = max(
            float(data_rows_by_id[row["row_id"]]["calibration_age_seconds"])
            for row in selected_rows
        )
        calibration_result = data["calibration_results"][calibration_id]
        control_pass = bool(
            float(calibration_result["gof_p_value"])
            >= float(calibration["control_rule"]["minimum_gof_p_value"])
            and calibration_result["minimum_count_per_prepared_state"]
            >= calibration["control_rule"]["minimum_count_per_prepared_state"]
        )
        age_pass = maximum_observed_age <= float(calibration["max_age_seconds"])
        eligible = calibration["primary_eligible"] is True
        primary_referenced = any(row["endpoint"] == PRIMARY_ENDPOINT for row in selected_rows)
        calibration_gate_details[calibration_id] = {
            "control_pass": control_pass,
            "control_gof_p_value": float(calibration_result["gof_p_value"]),
            "minimum_count_per_prepared_state": calibration_result[
                "minimum_count_per_prepared_state"
            ],
            "calibration_receipt_sha256": calibration_result[
                "calibration_receipt_sha256"
            ],
            "primary_eligible": eligible,
            "primary_referenced": primary_referenced,
            "maximum_observed_age_seconds": maximum_observed_age,
            "maximum_allowed_age_seconds": float(calibration["max_age_seconds"]),
            "age_pass": age_pass,
            "pass": bool(
                control_pass and age_pass and (eligible or not primary_referenced)
            ),
        }
    calibration_gate_pass = all(
        details["pass"] for details in calibration_gate_details.values()
    )
    validity_gate_pass = bool(calibration_gate_pass and leakage_gate_pass)

    every_backend_lr_passes = all(
        result["passes_per_backend_threshold"]
        for nulls in per_backend_likelihood_ratios.values()
        for result in nulls.values()
    )
    every_pooled_lr_passes = all(
        result["passes_pooled_threshold"] for result in pooled_likelihood_ratios.values()
    )
    label_weighted_scores = sorted(
        (
            (
                component["component_id"],
                math.log(float(component["prior_weight"]))
                + pooled_label_component_scores[component["component_id"]],
            )
            for component in lock["label_layout_model"]["components"]
        ),
        key=lambda item: (-item[1], item[0]),
    )
    label_top_gap = label_weighted_scores[0][1] - label_weighted_scores[1][1]
    label_unique_preference = bool(
        label_top_gap > 0.0
        and label_weighted_scores[0][0]
        == lock["label_layout_model"]["reference_component_id"]
    )
    primary_pass = bool(
        validity_gate_pass
        and every_backend_lr_passes
        and every_pooled_lr_passes
        and label_unique_preference
        and global_envelope["pass"]
    )

    decisive_null_failure = any(
        all(
            per_backend_likelihood_ratios[backend_name][null][
                "sensitivity_log_likelihood_ratio_bounds"
            ][1] < -math.log(100.0)
            for backend_name in per_backend_likelihood_ratios
        )
        for null in models[1:]
    )
    both_system_envelope_failure = all(
        not envelope["pass"] for envelope in backend_envelopes.values()
    )
    label_specificity_failure = not label_unique_preference
    kernel_failure = bool(
        validity_gate_pass
        and (
            decisive_null_failure
            or both_system_envelope_failure
            or label_specificity_failure
        )
    )
    if not validity_gate_pass:
        verdict = "invalid_calibration_or_leakage_gate"
    elif primary_pass:
        verdict = "passes_frozen_reduced_repair_kernel_gate"
    elif kernel_failure:
        verdict = "fails_frozen_reduced_repair_kernel"
    else:
        verdict = "valid_but_inconclusive"

    label_ranked = sorted(
        (
            {
                "component_id": component["component_id"],
                "prior_weight": float(component["prior_weight"]),
                "log_likelihood": pooled_label_component_scores[component["component_id"]],
                "log_weighted_likelihood": math.log(float(component["prior_weight"]))
                + pooled_label_component_scores[component["component_id"]],
            }
            for component in lock["label_layout_model"]["components"]
        ),
        key=lambda item: (-item["log_weighted_likelihood"], item["component_id"]),
    )
    for rank, item in enumerate(label_ranked, start=1):
        item["rank"] = rank
        item["posterior_within_label_model"] = math.exp(
            item["log_weighted_likelihood"] - pooled_scores[LABEL_MODEL]
        )
    best_label_score = max(pooled_label_component_scores.values())

    report_body = {
        "schema_version": REPORT_SCHEMA_VERSION,
        "blinded": True,
        "revealed": False,
        "analysis_lock_sha256": lock["analysis_lock_sha256"],
        "analysis_code_sha256": analysis_code_sha256(),
        "catalog_precommitment_sha256": lock["catalog_precommitment_sha256"],
        "heldout_data_packet_sha256": data["data_packet_sha256"],
        "manifest_sha256": data["manifest_sha256"],
        "submission_journal_sha256": data["submission_journal_sha256"],
        "harvest_journal_sha256": data["harvest_journal_sha256"],
        "validity": {
            "fail_closed_validation_passed": True,
            "all_expected_rows_present": True,
            "all_submitted_jobs_included": True,
            "shot_conservation_passed": True,
            "postselection_detected": False,
            "calibration_frozen_before_heldout": True,
            "calibration_gate_pass": calibration_gate_pass,
            "calibration_gates": calibration_gate_details,
            "leakage_gate_pass": leakage_gate_pass,
            "individual_catastrophic_leakage_pass": (
                individual_catastrophic_leakage_pass
            ),
            "pooled_backend_family_leakage_pass": pooled_leakage_pass,
            "pooled_backend_family_leakage_gates": pooled_leakage_gates,
            "overall_validity_gate_pass": validity_gate_pass,
        },
        "conditional_likelihood": {
            "formula": "Multinomial(n; p_obs), p_obs(y)=sum_x C(y|x)*p_latent(x)",
            "evidence_measure": "conditional_likelihood_ratio_only",
            "calibration_uncertainty_integrated": False,
            "positive_dirichlet_pseudocount_required": True,
            "frozen_sensitivity_bounds_required": True,
            "multinomial_constants_included": True,
            "calibration_receipts": calibration_receipts,
            "primary_pooled_log_likelihoods": pooled_scores,
            "primary_per_backend_log_likelihoods": per_backend_scores,
        },
        "primary_endpoint": {
            "name": PRIMARY_ENDPOINT,
            "backend_count": len(primary_by_backend),
            "row_count": len(primary_rows),
            "required_nulls": list(models[1:]),
            "per_backend_conditional_likelihood_ratios": per_backend_likelihood_ratios,
            "pooled_conditional_likelihood_ratios": pooled_likelihood_ratios,
            "per_backend_threshold": PER_BACKEND_LR_THRESHOLD,
            "pooled_threshold": POOLED_LR_THRESHOLD,
            "global_99_percent_simultaneous_envelope": global_envelope,
            "per_backend_99_percent_simultaneous_envelopes": backend_envelopes,
            "label_layout_multiplicity": {
                "mapping_scope": "global_shared_across_primary_rows",
                "component_count": len(label_ranked),
                "log_marginal_likelihood": pooled_scores[LABEL_MODEL],
                "ranked_components": label_ranked,
                "reference_component_id": lock["label_layout_model"][
                    "reference_component_id"
                ],
                "reference_component_is_unique_top": label_unique_preference,
                "top_to_second_log_weighted_likelihood_gap": label_top_gap,
                "repair_log_likelihood_gap_to_best_label_component": pooled_scores[
                    REPAIR_MODEL
                ]
                - best_label_score,
                "per_backend_component_log_likelihoods": per_backend_label_component_scores,
            },
            "passes_gate": primary_pass,
        },
        "secondary_family": {
            "correction": "Holm",
            "family_alpha": SECONDARY_FAMILY_ALPHA,
            "does_not_select_primary": True,
            "tests": holm_results,
        },
        "mh_secondary_families": mh_secondary_families,
        "shot_audit": shot_audit,
        "candidate_probability_hashes": {
            row_id: {
                model: {
                    "logical_circuit_sha256": expected_rows[row_id][
                        "logical_circuit_sha256"
                    ],
                    "latent_probability_table_sha256": expected_rows[row_id][
                        "candidate_provenance"
                    ][model]["probability_table_sha256"],
                    "derivation_sha256": expected_rows[row_id]["candidate_provenance"][
                        model
                    ]["derivation_sha256"],
                    "conditional_observed_prediction_sha256": sha256_json(
                        {
                            "logical_circuit_sha256": expected_rows[row_id][
                                "logical_circuit_sha256"
                            ],
                            "latent_probability_table_sha256": expected_rows[row_id][
                                "candidate_provenance"
                            ][model]["probability_table_sha256"],
                            "calibration_sha256": sha256_json(
                                lock["calibrations"][expected_rows[row_id]["calibration_id"]]
                            ),
                        }
                    ),
                }
                for model in candidates
            }
            for row_id, candidates in probabilities_by_row.items()
        },
        "label_component_probability_hashes": {
            component["component_id"]: {
                row_id: sha256_json(probabilities)
                for row_id, probabilities in component["row_probabilities"].items()
            }
            for component in lock["label_layout_model"]["components"]
        },
        "decision": {
            "verdict": verdict,
            "primary_pass": primary_pass,
            "kernel_failure": kernel_failure,
            "decisive_null_failure": decisive_null_failure,
            "both_system_envelope_failure": both_system_envelope_failure,
            "label_specificity_failure": label_specificity_failure,
            "invalid_calibration_or_leakage_gate": not validity_gate_pass,
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
    rng = np.random.default_rng(seed)
    rows = []
    for expected_row in lock["expected_rows"]:
        selected_generating_model = generating_model
        if (
            expected_row["family"] == "mh"
            and generating_model == REPAIR_MODEL
        ):
            selected_generating_model = MH_REFERENCE_MODEL
        if selected_generating_model not in expected_row["candidate_probabilities"]:
            raise AnalysisValidationError(
                "synthetic generating model does not match an expected row family"
            )
        calibration = lock["calibrations"][expected_row["calibration_id"]]
        probabilities = _candidate_observed_probabilities(
            expected_row,
            calibration,
            selected_generating_model,
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
                "opaque_id": expected_row["opaque_id"],
                "logical_circuit_sha256": expected_row["logical_circuit_sha256"],
                "backend_role": expected_row["backend_role"],
                "backend_name": expected_row["backend_name"],
                "layout_id": expected_row["layout_id"],
                "physical_layout": list(expected_row["physical_layout"]),
                "job_id": f"synthetic-job-{expected_row['backend_role']}",
                "group_id": f"synthetic-group-{expected_row['backend_role']}",
                "submission_event_sha256": "c" * 64,
                "harvest_event_sha256": "d" * 64,
                "raw_joined_counts_sha256": sha256_json(count_mapping),
                "calibration_age_seconds": 0.0,
                "declared_shots": shots,
                "submitted_shots": shots,
                "retrieved_shots": shots,
                "excluded_shots": 0,
                "postselected": False,
                "all_outcomes_included": True,
                "counts": count_mapping,
            }
        )
    return seal_data_packet(
        analysis_lock_sha256=lock["analysis_lock_sha256"],
        rows=rows,
        manifest_sha256="e" * 64,
        submission_journal_sha256="c" * 64,
        harvest_journal_sha256="d" * 64,
        source_kind="synthetic_preflight",
        calibration_results={
            calibration_id: basis_calibration_control_result(
                diagnostic_counts_by_opaque_id={
                    opaque_id: {f"{basis_code:04b}": 512}
                    for opaque_id, basis_code in calibration["control_rule"][
                        "expected_basis_code_by_opaque_id"
                    ].items()
                },
                control_rule=calibration["control_rule"],
                provider_job_ids=[f"synthetic-calibration-job-{calibration_id}"],
                calibration_receipt_sha256=sha256_json(
                    {"synthetic_calibration": calibration_id}
                ),
            )
            for calibration_id, calibration in lock["calibrations"].items()
        },
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
