#!/usr/bin/env python3
"""Strict fail-closed validator for the Ward-projected spectral-measure export.

This module is the single typed validator/construction boundary for the
production export contract (`ward_projected_spectral_measure.schema.json`).
It positively enforces every semantic requirement of the contract; anything
unknown, missing, malformed, contradictory, or non-finite fails closed.
Schema conformance alone is never acceptance: the semantic layer runs on
top of the JSON Schema and each check appends a typed reason, and a payload
is accepted only when the reason list is empty.

Positively enforced requirements (ids in SEMANTIC_REQUIREMENTS):

- finite_numeric_values: every number anywhere in the payload is finite
  (NaN/Inf rejected); booleans are never accepted where numbers are typed.
- nonempty_level_support: at least one ensemble block with at least one level.
- unique_level_identifiers: level_id values are unique across the payload.
- residue_level_reference_integrity: residue rows and levels are in exact
  one-to-one correspondence (every residue references an existing level,
  no dangling or duplicate references, every level carries a residue).
- positive_energies: E_n > 0 for every level.
- s_equals_energy_squared: s_n = E_n^2 within the declared relative
  tolerance S_EQUALS_ENERGY_SQUARED_REL_TOL.
- nonnegative_weights_and_residues: w_n >= 0 and residues >= 0.
- weight_residue_consistency: the declared relation between level weights
  and Ward-projected residues (residue_n = weight_n within
  WEIGHT_RESIDUE_REL_TOL; both are the coefficient Z_n of
  delta(s - E_n^2) in the zero-momentum transverse channel).
- complete_typed_provenance: a closed provenance manifest with at least one
  typed source input (kind from the fixed allowlist) and the two
  no-leak booleans exactly false.
- no_measured_hvp_or_target_inputs: external_targets_used must be an empty
  list, and a recursive scan over every key and string value rejects any
  occurrence of a forbidden target token or a downstream-forbidden key
  (the downstream FORBIDDEN_SOURCE_KEYS set is imported, not copied).
- covariance_dimension_and_symmetry: square matrix matching the declared
  dimension and row basis, symmetric within COVARIANCE_SYMMETRY_REL_TOL.
- covariance_positive_semidefinite: min eigenvalue >= -COVARIANCE_PSD_REL_TOL
  relative to the matrix scale.
- complete_budget_rows: exactly the seven required budget rows.
- finite_ordered_bound_intervals: every budget carries bound_interval with
  finite decimals 0 <= lo <= hi.
- transport_moment_certificate_complete: kernel, Delta_had_image interval,
  quadrature_error_bound and tail_bound, finite and correctly ordered
  (this is the block the downstream source-transport validator requires).
- strict_branch_typing: flavors exactly "2+1", qed from the closed enum;
  quenched, surrogate, target-calibrated and compare-only payloads are
  excluded from acceptance.
- guard_flags_exactly_false: the three guard booleans are literally False
  (strings, 0, None, or missing fail).
- ward_projection_required: projection.ward_projected is literally True and
  the lane is U(1)_Q.
- closed_object_boundaries: unknown keys at any object boundary are
  rejected (additionalProperties: false everywhere in the schema).
- positivity_status_allowlisted: rho_had_or_measure.positivity_status must
  be exactly "certified_positive"; unknown status strings fail closed.

The gate is downstream-compatible by construction: an accepted payload
embeds as `source_measure` into the source-transport payload contract of
code/P_derivation/thomson_spectral_transport.py and passes its
_validate_source_measure checks; verify_issue_317_spectral_measure_packet.py
executes that implication as an integration witness.
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
P_DERIVATION = ROOT / "P_derivation"
if str(P_DERIVATION) not in sys.path:
    sys.path.insert(0, str(P_DERIVATION))

from thomson_spectral_transport import (  # noqa: E402
    FORBIDDEN_SOURCE_KEYS as DOWNSTREAM_FORBIDDEN_KEYS,
    SOURCE_SYSTEMATICS_BUDGETS,
)

SCHEMA_PATH = HERE / "ward_projected_spectral_measure.schema.json"

FORBIDDEN_TARGETS = (
    "CODATA_ALPHA",
    "MUON_G_MINUS_2",
    "EE_TO_HADRONS",
    "RARE_DECAY_DATA",
    "HADRON_MASS_TARGETS",
    "PDG_QCD_FITS",
)

# Declared numerical tolerances of the contract.
S_EQUALS_ENERGY_SQUARED_REL_TOL = 1e-12
WEIGHT_RESIDUE_REL_TOL = 1e-12
COVARIANCE_SYMMETRY_REL_TOL = 1e-12
COVARIANCE_PSD_REL_TOL = 1e-12

ALLOWED_SOURCE_INPUT_KINDS = (
    "bare_lattice_parameter",
    "source_ensemble",
    "oph_source_artifact",
    "declared_convention",
    "backend_provenance",
)

ALLOWED_POSITIVITY_STATUS = ("certified_positive",)

ALLOWED_FLAVORS = ("2+1",)
ALLOWED_QED = ("off", "explicitly_superseded")

SEMANTIC_REQUIREMENTS = (
    "finite_numeric_values",
    "nonempty_level_support",
    "unique_level_identifiers",
    "residue_level_reference_integrity",
    "positive_energies",
    "s_equals_energy_squared",
    "nonnegative_weights_and_residues",
    "weight_residue_consistency",
    "complete_typed_provenance",
    "no_measured_hvp_or_target_inputs",
    "covariance_dimension_and_symmetry",
    "covariance_positive_semidefinite",
    "complete_budget_rows",
    "finite_ordered_bound_intervals",
    "transport_moment_certificate_complete",
    "strict_branch_typing",
    "guard_flags_exactly_false",
    "ward_projection_required",
    "closed_object_boundaries",
    "positivity_status_allowlisted",
)


@dataclass(frozen=True)
class ValidationResult:
    accepted: bool
    reasons: tuple[str, ...] = field(default_factory=tuple)

    def to_json(self) -> dict[str, Any]:
        return {"accepted": self.accepted, "reasons": list(self.reasons)}


def load_schema() -> dict[str, Any]:
    return json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))


def _is_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool)


def _finite_number(value: Any, where: str, reasons: list[str]) -> float | None:
    if not _is_number(value):
        reasons.append(f"not_a_number:{where}")
        return None
    if not math.isfinite(value):
        reasons.append(f"nonfinite_number:{where}")
        return None
    return float(value)


def _scan_finiteness(payload: Any, where: str, reasons: list[str]) -> None:
    """Reject NaN/Inf anywhere in the payload, at any nesting depth."""
    if isinstance(payload, dict):
        for key, value in payload.items():
            _scan_finiteness(value, f"{where}.{key}" if where else str(key), reasons)
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            _scan_finiteness(value, f"{where}[{index}]", reasons)
    elif _is_number(payload) and not math.isfinite(payload):
        reasons.append(f"nonfinite_number:{where}")


def _scan_forbidden_tokens(payload: Any, where: str, reasons: list[str]) -> None:
    """Recursive forbidden-token scan over every key and string value.

    Rejects (a) any key or string value containing a forbidden target token
    (case-insensitive), and (b) any key equal to a downstream-forbidden
    source key. This closes the nested-provenance leak class, e.g.
    {"measured_hvp_input": {"source": "EE_TO_HADRONS"}}, independently of
    whether the schema boundary already rejects the unknown key.
    """
    if isinstance(payload, dict):
        for key, value in payload.items():
            key_text = str(key)
            key_where = f"{where}.{key_text}" if where else key_text
            upper = key_text.upper()
            for token in FORBIDDEN_TARGETS:
                if token in upper:
                    reasons.append(f"TARGET_LEAK_DETECTED:key:{key_where}")
            if key_text in DOWNSTREAM_FORBIDDEN_KEYS:
                reasons.append(f"downstream_forbidden_key:{key_where}")
            _scan_forbidden_tokens(value, key_where, reasons)
    elif isinstance(payload, list):
        for index, value in enumerate(payload):
            _scan_forbidden_tokens(value, f"{where}[{index}]", reasons)
    elif isinstance(payload, str):
        upper = payload.upper()
        for token in FORBIDDEN_TARGETS:
            if token in upper:
                reasons.append(f"TARGET_LEAK_DETECTED:value:{where}")


def _decimal_bound(value: Any, where: str, reasons: list[str]) -> Decimal | None:
    """Parse a nonnegative finite decimal from a number or string."""
    if isinstance(value, bool) or value is None:
        reasons.append(f"invalid_decimal:{where}")
        return None
    if _is_number(value) and not math.isfinite(value):
        reasons.append(f"nonfinite_number:{where}")
        return None
    try:
        parsed = Decimal(str(value))
    except (InvalidOperation, ValueError):
        reasons.append(f"invalid_decimal:{where}")
        return None
    if not parsed.is_finite():
        reasons.append(f"nonfinite_number:{where}")
        return None
    if parsed < 0:
        reasons.append(f"negative_decimal:{where}")
        return None
    return parsed


def _check_interval(value: Any, where: str, reasons: list[str]) -> None:
    if not isinstance(value, dict):
        reasons.append(f"required_interval_missing:{where}")
        return
    lo = _decimal_bound(value.get("lo"), f"{where}.lo", reasons)
    hi = _decimal_bound(value.get("hi"), f"{where}.hi", reasons)
    if lo is not None and hi is not None and lo > hi:
        reasons.append(f"invalid_interval_order:{where}")


def _check_schema(payload: dict[str, Any], schema: dict[str, Any], reasons: list[str]) -> None:
    import jsonschema

    validator = jsonschema.Draft202012Validator(schema)
    for error in sorted(validator.iter_errors(payload), key=lambda e: e.json_path):
        reasons.append(f"schema_violation:{error.json_path}:{error.message}")


def _check_provenance(payload: dict[str, Any], reasons: list[str]) -> None:
    targets = payload.get("external_targets_used")
    if not isinstance(targets, list):
        reasons.append("external_targets_used_not_a_list")
    elif targets:
        reasons.append(
            "TARGET_LEAK_DETECTED:external_targets_used:" + ",".join(str(t) for t in targets)
        )
    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        reasons.append("provenance_manifest_missing")
        return
    source_inputs = provenance.get("source_inputs")
    if not isinstance(source_inputs, list) or not source_inputs:
        reasons.append("provenance_source_inputs_missing_or_empty")
    else:
        for index, item in enumerate(source_inputs):
            if not isinstance(item, dict):
                reasons.append(f"provenance_source_input_not_typed:[{index}]")
                continue
            if item.get("kind") not in ALLOWED_SOURCE_INPUT_KINDS:
                reasons.append(
                    f"provenance_source_input_kind_not_allowlisted:[{index}]:{item.get('kind')!r}"
                )
            identifier = item.get("identifier")
            if not isinstance(identifier, str) or not identifier:
                reasons.append(f"provenance_source_input_identifier_missing:[{index}]")
    if provenance.get("measured_hvp_input_present") is not False:
        reasons.append("provenance_measured_hvp_input_present_not_false")
    if provenance.get("target_calibration_present") is not False:
        reasons.append("provenance_target_calibration_present_not_false")


def _check_levels_and_residues(payload: dict[str, Any], reasons: list[str]) -> None:
    blocks = payload.get("finite_volume_levels")
    if not isinstance(blocks, list) or not blocks:
        reasons.append("finite_volume_levels_missing_or_empty")
        blocks = []
    level_weights: dict[str, float] = {}
    seen_ids: set[str] = set()
    for b_index, block in enumerate(blocks):
        if not isinstance(block, dict):
            reasons.append(f"level_block_not_typed:[{b_index}]")
            continue
        levels = block.get("levels")
        if not isinstance(levels, list) or not levels:
            reasons.append(f"levels_missing_or_empty:[{b_index}]")
            continue
        for l_index, level in enumerate(levels):
            where = f"finite_volume_levels[{b_index}].levels[{l_index}]"
            if not isinstance(level, dict):
                reasons.append(f"level_not_typed:{where}")
                continue
            level_id = level.get("level_id")
            if not isinstance(level_id, str) or not level_id:
                reasons.append(f"level_id_missing:{where}")
                level_id = None
            elif level_id in seen_ids:
                reasons.append(f"duplicate_level_id:{level_id}")
            else:
                seen_ids.add(level_id)
            energy = _finite_number(level.get("energy"), f"{where}.energy", reasons)
            s_value = _finite_number(level.get("s"), f"{where}.s", reasons)
            weight = _finite_number(level.get("weight"), f"{where}.weight", reasons)
            if energy is not None and energy <= 0.0:
                reasons.append(f"nonpositive_energy:{where}")
            if weight is not None and weight < 0.0:
                reasons.append(f"negative_level_weight:{where}")
            if energy is not None and s_value is not None and energy > 0.0:
                expected = energy * energy
                if abs(s_value - expected) > S_EQUALS_ENERGY_SQUARED_REL_TOL * max(1.0, expected):
                    reasons.append(f"s_not_equal_energy_squared:{where}")
            for error_key in ("energy_jackknife_error", "weight_jackknife_error"):
                if error_key in level:
                    err = _finite_number(level[error_key], f"{where}.{error_key}", reasons)
                    if err is not None and err < 0.0:
                        reasons.append(f"negative_error:{where}.{error_key}")
            if level_id is not None and weight is not None:
                level_weights[level_id] = weight

    residues = payload.get("ward_projected_residues")
    if not isinstance(residues, list) or not residues:
        reasons.append("ward_projected_residues_missing_or_empty")
        residues = []
    residue_ids: set[str] = set()
    for r_index, row in enumerate(residues):
        where = f"ward_projected_residues[{r_index}]"
        if not isinstance(row, dict):
            reasons.append(f"residue_row_not_typed:{where}")
            continue
        level_id = row.get("level_id")
        if not isinstance(level_id, str) or not level_id:
            reasons.append(f"residue_level_id_missing:{where}")
            level_id = None
        elif level_id in residue_ids:
            reasons.append(f"duplicate_residue_reference:{level_id}")
        else:
            residue_ids.add(level_id)
        if level_id is not None and level_id not in seen_ids:
            reasons.append(f"residue_references_unknown_level:{level_id}")
        residue = _finite_number(row.get("residue"), f"{where}.residue", reasons)
        if residue is not None and residue < 0.0:
            reasons.append(f"negative_residue:{where}")
        normalization = row.get("current_normalization")
        if not isinstance(normalization, str) or not normalization:
            reasons.append(f"current_normalization_missing:{where}")
        if (
            residue is not None
            and level_id is not None
            and level_id in level_weights
        ):
            weight = level_weights[level_id]
            if abs(residue - weight) > WEIGHT_RESIDUE_REL_TOL * max(1.0, abs(weight)):
                reasons.append(f"weight_residue_inconsistent:{level_id}")
    for level_id in sorted(seen_ids - residue_ids):
        reasons.append(f"level_without_residue:{level_id}")


def _check_measure_block(payload: dict[str, Any], reasons: list[str]) -> None:
    rho = payload.get("rho_had_or_measure")
    if not isinstance(rho, dict):
        reasons.append("rho_had_or_measure_missing")
        return
    if rho.get("support_variable") != "s":
        reasons.append("support_variable_not_s")
    if not rho.get("pushforward_rule"):
        reasons.append("pushforward_rule_missing")
    if rho.get("positivity_status") not in ALLOWED_POSITIVITY_STATUS:
        reasons.append(
            f"positivity_status_not_allowlisted:{rho.get('positivity_status')!r}"
        )


def _check_covariance(payload: dict[str, Any], reasons: list[str]) -> None:
    covariance = payload.get("covariance")
    if not isinstance(covariance, dict):
        reasons.append("covariance_missing")
        return
    dimension = covariance.get("dimension")
    matrix = covariance.get("matrix")
    row_basis = covariance.get("row_basis")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 1:
        reasons.append("covariance_dimension_invalid")
        return
    if not isinstance(row_basis, list) or len(row_basis) != dimension:
        reasons.append("covariance_row_basis_dimension_mismatch")
    if not isinstance(matrix, list) or len(matrix) != dimension:
        reasons.append("covariance_matrix_dimension_mismatch")
        return
    rows: list[list[float]] = []
    for r_index, row in enumerate(matrix):
        if not isinstance(row, list) or len(row) != dimension:
            reasons.append(f"covariance_matrix_dimension_mismatch:row[{r_index}]")
            return
        parsed_row: list[float] = []
        for c_index, value in enumerate(row):
            entry = _finite_number(value, f"covariance.matrix[{r_index}][{c_index}]", reasons)
            if entry is None:
                return
            parsed_row.append(entry)
        rows.append(parsed_row)
    arr = np.array(rows, dtype=float)
    scale = float(np.max(np.abs(arr))) if arr.size else 0.0
    tol_scale = max(1.0, scale)
    if float(np.max(np.abs(arr - arr.T))) > COVARIANCE_SYMMETRY_REL_TOL * tol_scale:
        reasons.append("covariance_not_symmetric")
        return
    min_eig = float(np.linalg.eigvalsh(0.5 * (arr + arr.T)).min())
    if min_eig < -COVARIANCE_PSD_REL_TOL * tol_scale:
        reasons.append("covariance_not_positive_semidefinite")


def _check_transport_certificate(payload: dict[str, Any], reasons: list[str]) -> None:
    certificate = payload.get("transport_moment_certificate")
    if not isinstance(certificate, dict):
        reasons.append("transport_moment_certificate_missing")
        return
    kernel = certificate.get("kernel")
    if not isinstance(kernel, str) or not kernel:
        reasons.append("transport_moment_certificate_kernel_missing")
    _check_interval(
        certificate.get("Delta_had_image"),
        "transport_moment_certificate.Delta_had_image",
        reasons,
    )
    _decimal_bound(
        certificate.get("quadrature_error_bound"),
        "transport_moment_certificate.quadrature_error_bound",
        reasons,
    )
    _decimal_bound(
        certificate.get("tail_bound"),
        "transport_moment_certificate.tail_bound",
        reasons,
    )


def _check_budgets(payload: dict[str, Any], reasons: list[str]) -> None:
    systematics = payload.get("systematics")
    if not isinstance(systematics, dict):
        reasons.append("systematics_missing")
        return
    for name in SOURCE_SYSTEMATICS_BUDGETS:
        budget = systematics.get(name)
        if not isinstance(budget, dict):
            reasons.append(f"budget_missing:{name}")
            continue
        if "bound_interval" not in budget:
            reasons.append(f"budget_bound_interval_missing:{name}")
            continue
        _check_interval(budget["bound_interval"], f"systematics.{name}.bound_interval", reasons)
    for name in systematics:
        if name not in SOURCE_SYSTEMATICS_BUDGETS:
            reasons.append(f"unknown_budget_row:{name}")


def _check_typing_and_guards(payload: dict[str, Any], reasons: list[str]) -> None:
    branch = payload.get("branch")
    if not isinstance(branch, dict):
        reasons.append("branch_missing")
    else:
        if branch.get("flavors") not in ALLOWED_FLAVORS:
            reasons.append(f"branch_flavors_not_allowlisted:{branch.get('flavors')!r}")
        if branch.get("qed") not in ALLOWED_QED:
            reasons.append(f"branch_qed_not_allowlisted:{branch.get('qed')!r}")
    projection = payload.get("projection")
    if not isinstance(projection, dict):
        reasons.append("projection_missing")
    else:
        if projection.get("lane") != "U(1)_Q":
            reasons.append("projection_lane_not_u1q")
        if projection.get("ward_projected") is not True:
            reasons.append("ward_projection_dropped")
    guards = payload.get("guards")
    if not isinstance(guards, dict):
        reasons.append("guards_missing")
        return
    for name in ("stable_channel_only", "surrogate_hadron_artifact", "compare_only_external_endpoint"):
        if guards.get(name) is not False:
            reasons.append(f"guard_not_literal_false:{name}")


def validate_production_payload(
    payload: Any, schema: dict[str, Any] | None = None
) -> ValidationResult:
    """Strict fail-closed acceptance decision for a production export payload."""
    if not isinstance(payload, dict):
        return ValidationResult(False, ("payload_not_an_object",))
    schema = schema or load_schema()
    reasons: list[str] = []
    _check_schema(payload, schema, reasons)
    _scan_finiteness(payload, "", reasons)
    _scan_forbidden_tokens(payload, "", reasons)
    _check_provenance(payload, reasons)
    _check_levels_and_residues(payload, reasons)
    _check_measure_block(payload, reasons)
    _check_covariance(payload, reasons)
    _check_transport_certificate(payload, reasons)
    _check_budgets(payload, reasons)
    _check_typing_and_guards(payload, reasons)
    deduped = tuple(dict.fromkeys(reasons))
    return ValidationResult(accepted=not deduped, reasons=deduped)
