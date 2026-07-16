#!/usr/bin/env python3
"""Fail-closed comparator for a sealed source artifact and corrective target.

``run_bracket.py`` never reads a target. This separate process reads both
objects after emission, verifies the corrective target's coordinate algebra,
and refuses a verdict unless a future detached successor is externally
activated and the artifact is a certified function or enclosure over its
registered P domain.

The 2026-07-16 v3 contract is a permanently inactive post-target-access
erratum scaffold. The current sampled singleton bracket is intentionally not
certified. Therefore a real invocation currently returns ``NOT_EVALUABLE``.
No point-diagnostic containment shortcut is implemented.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import string
import sys
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
from typing import Any

TARGET_ARTIFACT = "oph_hadronic_closure_corrective_target_contract"
EMISSION_ARTIFACT = "oph_ward_projected_payload_source_bracket"
EMISSION_SCHEMA_VERSION = 3
MACHINE_CONTRACT_SCHEMA_VERSION = 1
PI_DECIMAL = Decimal(
    "3.141592653589793238462643383279502884197169399375105820974944592307816406286"
)
CURRENT_INACTIVE_CORRECTIVE_ID = "hadronic_closure_target_2026-07-16_v3"

TOTAL_COORDINATE = "delta_source_total_alpha_inv"
RESIDUAL_COORDINATE = "delta_source_residual_vs_implemented_alpha_inv"
S_QEW_COORDINATE = "s_qew_effective"
S_HADRONIC_COORDINATE = "s_hadronic"
COORDINATES = (
    TOTAL_COORDINATE,
    RESIDUAL_COORDINATE,
    S_QEW_COORDINATE,
    S_HADRONIC_COORDINATE,
)

REQUIRED_COORDINATE_TYPES = {
    TOTAL_COORDINATE: ("total", "inverse_alpha", "map_input_only"),
    RESIDUAL_COORDINATE: ("residual", "inverse_alpha", "diagnostic_only"),
    S_QEW_COORDINATE: (
        "screening_ratio_qew",
        "dimensionless",
        "diagnostic_only",
    ),
    S_HADRONIC_COORDINATE: (
        "screening_ratio_hadronic",
        "dimensionless",
        "diagnostic_only",
    ),
}

REQUIRED_RECEIPT_FIELDS = {
    "target_free_dependency_dag",
    "forbidden_input_scan",
    "source_commit",
    "source_tree_sha256",
    "generator_argv",
    "environment_lock",
    "python_version",
    "numeric_library_versions",
    "precision_and_cutoffs",
    "executable_dependency_sha256",
    "canonical_json_sha256",
    "external_timestamp",
    "payload_work_started_utc",
}


class ScoringError(ValueError):
    """Fail-closed protocol error with a stable machine code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


def _mapping(value: Any, path: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ScoringError("schema_mismatch", f"{path} must be an object")
    return value


def _field(mapping: dict[str, Any], key: str, path: str) -> Any:
    if key not in mapping:
        raise ScoringError("schema_mismatch", f"missing {path}.{key}")
    return mapping[key]


def _decimal(value: Any, path: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise ScoringError("schema_mismatch", f"{path} must be a finite decimal")
    try:
        result = Decimal(str(value))
    except InvalidOperation as exc:
        raise ScoringError("schema_mismatch", f"{path} is not a decimal") from exc
    if not result.is_finite():
        raise ScoringError("schema_mismatch", f"{path} must be finite")
    return result


def _decimal_string(value: Any, path: str) -> Decimal:
    if not isinstance(value, str):
        raise ScoringError("schema_mismatch", f"{path} must be a decimal string")
    return _decimal(value, path)


def _sha256(value: Any, path: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in string.hexdigits for character in value)
    ):
        raise ScoringError("schema_mismatch", f"{path} must be a SHA-256 hex digest")
    return value.lower()


def _nonzero_sha256(value: Any, path: str) -> str:
    digest = _sha256(value, path)
    if digest == "0" * 64:
        raise ScoringError("schema_mismatch", f"{path} cannot be an all-zero digest")
    return digest


def _nonempty_string(value: Any, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ScoringError("schema_mismatch", f"{path} must be a nonempty string")
    return value.strip()


def _utc_timestamp(value: Any, path: str) -> datetime:
    raw = _nonempty_string(value, path)
    try:
        parsed = datetime.fromisoformat(raw.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ScoringError(
            "schema_mismatch", f"{path} is not an ISO-8601 timestamp"
        ) from exc
    if parsed.tzinfo is None:
        raise ScoringError("schema_mismatch", f"{path} must include a timezone")
    return parsed.astimezone(timezone.utc)


def _verified_receipt(value: Any, path: str) -> dict[str, Any]:
    receipt = _mapping(value, path)
    if receipt.get("status") != "PASS" or receipt.get("verified") is not True:
        raise ScoringError(
            "artifact_provenance_mismatch",
            f"{path} must be a verified PASS receipt",
        )
    _nonzero_sha256(_field(receipt, "sha256", path), f"{path}.sha256")
    return receipt


def _difference_equals(left: Decimal, right: Decimal, expected: Decimal) -> bool:
    with localcontext() as context:
        context.prec = 100
        return left - right == expected


def _s_separation_matches(
    left: Decimal,
    right: Decimal,
    alpha_u: Decimal,
    q_naive: Decimal,
) -> bool:
    with localcontext() as context:
        context.prec = 100
        return abs((left - right) - alpha_u / q_naive) <= Decimal("1e-36")


def _decimal_close(left: Decimal, right: Decimal, tolerance: Decimal) -> bool:
    with localcontext() as context:
        context.prec = 100
        return abs(left - right) <= tolerance


def _at_path(mapping: dict[str, Any], dotted_path: str) -> Any:
    value: Any = mapping
    walked: list[str] = []
    for key in dotted_path.split("."):
        parent = ".".join(walked) or "root"
        value = _field(_mapping(value, parent), key, parent)
        walked.append(key)
    return value


def artifact_content_sha256(artifact: dict[str, Any]) -> str:
    """Recompute the source emitter's deterministic embedded hash."""
    digest_source = {
        key: value
        for key, value in artifact.items()
        if key not in {"content_sha256", "wall_time_seconds"}
    }
    # A future certified artifact stores this same digest in its provenance
    # block.  Exclude that redundant copy to avoid a self-referential hash.
    if isinstance(digest_source.get("receipts"), dict):
        digest_source["receipts"] = dict(digest_source["receipts"])
        digest_source["receipts"].pop("canonical_json_sha256", None)
    try:
        canonical = json.dumps(
            digest_source,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ScoringError(
            "artifact_schema_mismatch", "artifact is not canonical finite JSON"
        ) from exc
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _validate_interval(value: Any, path: str) -> dict[str, Decimal]:
    interval = _mapping(value, path)
    # A rigorous enclosure cannot use binary JSON floats as purported
    # outward-rounded bounds.  The recorded width is redundant, but when
    # present it must equal the Decimal endpoint difference exactly.
    for field in ("lo", "hi", "width"):
        if not isinstance(_field(interval, field, path), str):
            raise ScoringError(
                "coordinate_schema_mismatch",
                f"{path}.{field} must be a decimal string",
            )
    lo = _decimal(interval["lo"], f"{path}.lo")
    hi = _decimal(interval["hi"], f"{path}.hi")
    width = _decimal(interval["width"], f"{path}.width")
    if lo > hi:
        raise ScoringError("coordinate_schema_mismatch", f"{path}.lo exceeds hi")
    if width < 0 or width != hi - lo:
        raise ScoringError(
            "coordinate_schema_mismatch",
            f"{path}.width does not equal hi - lo",
        )
    return {"lo": lo, "hi": hi, "width": width}


def _validate_corrective_target(target: dict[str, Any]) -> None:
    """Validate the algebra and epistemic gates present in corrective v3."""
    if target.get("artifact") != TARGET_ARTIFACT:
        raise ScoringError("target_schema_mismatch", "unexpected target artifact")
    _nonempty_string(target.get("id"), "target.id")
    _nonzero_sha256(
        _field(target, "historical_v2_sha256", "target"),
        "target.historical_v2_sha256",
    )

    definitions = _mapping(
        _field(target, "coordinate_definitions", "target"),
        "target.coordinate_definitions",
    )
    required_definitions = {
        "Delta_source_total",
        "Delta_source_residual_vs_implemented",
        "S_QEW_effective",
        "S_hadronic",
        "type_rule",
    }
    if set(definitions) != required_definitions:
        raise ScoringError(
            "target_coordinate_mismatch",
            "target must define total, residual, S_QEW, and S_hadronic separately",
        )

    shared = _mapping(
        _field(target, "shared_payload_contract", "target"),
        "target.shared_payload_contract",
    )
    if shared.get("singleton_P_payload_eligible") is not False:
        raise ScoringError(
            "target_schema_mismatch", "singleton-P payloads must be ineligible"
        )
    if shared.get("sampled_grid_extrema_interval_certificate") is not False:
        raise ScoringError(
            "target_schema_mismatch",
            "sampled grid extrema cannot be an interval certificate",
        )
    if shared.get("Delta_EW_zero_status") != "declared_zero_branch_unproven":
        raise ScoringError(
            "target_coordinate_mismatch", "corrective target must leave Delta_EW open"
        )

    primary = _mapping(
        _field(target, "primary_scoring_rule", "target"),
        "target.primary_scoring_rule",
    )
    if primary.get("producer_and_comparator_separate") is not True:
        raise ScoringError(
            "target_schema_mismatch", "producer and comparator must be separate"
        )
    if primary.get("map_verdicts_independent") is not True:
        raise ScoringError(
            "target_schema_mismatch", "CL-1 and CL-2 verdicts must be independent"
        )
    if primary.get("point_diagnostics_can_decide_verdict") is not False:
        raise ScoringError(
            "target_schema_mismatch", "point diagnostics cannot decide a verdict"
        )
    expected_policy_status = (
        "not_yet_registered"
        if target.get("id") == CURRENT_INACTIVE_CORRECTIVE_ID
        else "registered_in_machine_scoring_contract"
    )
    if primary.get("decision_policy_status") != expected_policy_status:
        raise ScoringError(
            "target_schema_mismatch",
            "decision-policy status is inconsistent with the target version",
        )

    measurement = _mapping(
        _field(target, "measurement_coordinate", "target"),
        "target.measurement_coordinate",
    )
    a_target = _decimal(
        _field(measurement, "A_target_inverse_alpha", "target.measurement_coordinate"),
        "target.measurement_coordinate.A_target_inverse_alpha",
    )
    u_a = _decimal(
        _field(measurement, "u_A_inverse_alpha", "target.measurement_coordinate"),
        "target.measurement_coordinate.u_A_inverse_alpha",
    )
    if measurement.get("u_A_semantics") != (
        "CODATA quoted standard uncertainty (1 sigma); not by itself a hard falsification bound"
    ):
        raise ScoringError(
            "target_coordinate_mismatch",
            "measurement uncertainty semantics are missing or ambiguous",
        )
    p_target = _decimal(
        _field(
            measurement, "P_target_point_diagnostic", "target.measurement_coordinate"
        ),
        "target.measurement_coordinate.P_target_point_diagnostic",
    )
    a0 = _decimal(
        _field(measurement, "a0_at_P_target_point", "target.measurement_coordinate"),
        "target.measurement_coordinate.a0_at_P_target_point",
    )
    lepton = _decimal(
        _field(
            measurement,
            "Delta_lepton_at_P_target_point",
            "target.measurement_coordinate",
        ),
        "target.measurement_coordinate.Delta_lepton_at_P_target_point",
    )
    implemented = _decimal(
        _field(
            measurement,
            "Delta_source_implemented_at_P_target_point",
            "target.measurement_coordinate",
        ),
        "target.measurement_coordinate.Delta_source_implemented_at_P_target_point",
    )
    alpha_u = _decimal(
        _field(
            measurement,
            "alpha_U_at_P_target_point",
            "target.measurement_coordinate",
        ),
        "target.measurement_coordinate.alpha_U_at_P_target_point",
    )
    q_naive = _decimal(
        _field(
            measurement,
            "Delta_quark_naive_at_P_target_point",
            "target.measurement_coordinate",
        ),
        "target.measurement_coordinate.Delta_quark_naive_at_P_target_point",
    )
    if min(a_target, u_a, p_target, q_naive) <= 0:
        raise ScoringError(
            "target_coordinate_mismatch",
            "A, u_A, P, and Delta_quark_naive must be positive",
        )
    if measurement.get("P_of_A_formula") != "P(A)=phi+sqrt(pi)/A":
        raise ScoringError(
            "target_coordinate_mismatch", "unexpected P-of-A diagnostic formula"
        )
    with localcontext() as context:
        context.prec = 100
        phi = (Decimal(1) + Decimal(5).sqrt()) / Decimal(2)
        expected_p = phi + PI_DECIMAL.sqrt() / a_target
    if not _decimal_close(p_target, expected_p, Decimal("1e-36")):
        raise ScoringError(
            "target_coordinate_mismatch",
            "P target is inconsistent with A_target and the declared adapter",
        )

    maps = _mapping(_field(target, "map_targets", "target"), "target.map_targets")
    if set(maps) != {"CL-1", "CL-2"}:
        raise ScoringError(
            "target_schema_mismatch", "target must contain separate CL-1 and CL-2 maps"
        )

    totals: dict[str, Decimal] = {}
    s_qew: dict[str, Decimal] = {}
    expected_map_contracts = {
        "CL-1": (
            "CL1_FULL_THOMSON",
            "G1(alpha)=1/(a0(P(alpha))+Delta_source_total(P(alpha)))",
        ),
        "CL-2": (
            "CL2_FULL_THOMSON_PLUS_ALPHA_U",
            "G2(alpha)=1/(a0(P(alpha))+Delta_source_total(P(alpha))+alpha_U(P(alpha)))",
        ),
    }
    for map_name in ("CL-1", "CL-2"):
        map_target = _mapping(maps[map_name], f"target.map_targets.{map_name}")
        if map_target.get("closure_rows") != [map_name]:
            raise ScoringError(
                "target_schema_mismatch", f"{map_name} must govern only its own row"
            )
        expected_kind, expected_formula = expected_map_contracts[map_name]
        if (
            map_target.get("map_kind") != expected_kind
            or map_target.get("map_formula") != expected_formula
        ):
            raise ScoringError(
                "target_coordinate_mismatch",
                f"{map_name} does not name its exact completed-map contract",
            )
        point = _mapping(
            _field(
                map_target, "point_diagnostics_only", f"target.map_targets.{map_name}"
            ),
            f"target.map_targets.{map_name}.point_diagnostics_only",
        )
        total = _decimal(
            _field(
                point,
                "Delta_source_total_target",
                f"target.map_targets.{map_name}.point_diagnostics_only",
            ),
            f"target.map_targets.{map_name}.point_diagnostics_only.Delta_source_total_target",
        )
        residual = _decimal(
            _field(
                point,
                "Delta_source_residual_vs_implemented",
                f"target.map_targets.{map_name}.point_diagnostics_only",
            ),
            f"target.map_targets.{map_name}.point_diagnostics_only.Delta_source_residual_vs_implemented",
        )
        if not _difference_equals(total, implemented, residual):
            raise ScoringError(
                "target_coordinate_mismatch",
                f"{map_name} total and residual use inconsistent baselines",
            )
        if point.get("S_hadronic_target") is not None:
            raise ScoringError(
                "target_coordinate_mismatch",
                f"{map_name} cannot name S_QEW as a hadronic target while Delta_EW is open",
            )
        with localcontext() as context:
            context.prec = 100
            expected_total = (
                a_target - a0 - (alpha_u if map_name == "CL-2" else Decimal(0))
            )
        if total != expected_total:
            raise ScoringError(
                "target_coordinate_mismatch",
                f"{map_name} total does not follow from A_target, a0, and its map formula",
            )
        totals[map_name] = total
        s_qew[map_name] = _decimal(
            _field(
                point,
                "S_QEW_effective_target",
                f"target.map_targets.{map_name}.point_diagnostics_only",
            ),
            f"target.map_targets.{map_name}.point_diagnostics_only.S_QEW_effective_target",
        )
        with localcontext() as context:
            context.prec = 100
            expected_s = (total - lepton) / q_naive
        if not _decimal_close(s_qew[map_name], expected_s, Decimal("5e-28")):
            raise ScoringError(
                "target_coordinate_mismatch",
                f"{map_name} S_QEW diagnostic is inconsistent with its total",
            )

    if not _difference_equals(totals["CL-1"], totals["CL-2"], alpha_u):
        raise ScoringError(
            "target_coordinate_mismatch",
            "CL-1 and CL-2 total targets must differ by alpha_U",
        )
    if not _s_separation_matches(s_qew["CL-1"], s_qew["CL-2"], alpha_u, q_naive):
        raise ScoringError(
            "target_coordinate_mismatch",
            "CL-1 and CL-2 S_QEW diagnostics have inconsistent separation",
        )


def _validate_machine_contract(contract: dict[str, Any]) -> dict[str, Any]:
    if contract.get("schema_version") != MACHINE_CONTRACT_SCHEMA_VERSION:
        raise ScoringError(
            "target_schema_mismatch", "unexpected machine_scoring_contract schema"
        )
    if contract.get("payload_artifact") != EMISSION_ARTIFACT:
        raise ScoringError(
            "target_schema_mismatch", "machine contract names the wrong artifact"
        )
    if contract.get("payload_schema_version") != EMISSION_SCHEMA_VERSION:
        raise ScoringError(
            "target_schema_mismatch", "machine contract names the wrong payload schema"
        )
    if contract.get("payload_object") != (
        "target_blind_function_or_certified_interval_enclosure"
    ):
        raise ScoringError(
            "target_schema_mismatch",
            "machine contract permits an ineligible payload object",
        )
    for key in ("source_family_id", "scheme_id", "current_definition_id"):
        if not isinstance(contract.get(key), str) or not contract[key]:
            raise ScoringError(
                "target_schema_mismatch", f"missing machine contract {key}"
            )
    if contract.get("source_method_frozen_before_target_access") is not True:
        raise ScoringError(
            "target_schema_mismatch",
            "eligible source method must be frozen before target access",
        )
    if contract.get("target_access_policy") not in {
        "withheld_future_data_release",
        "audited_clean_room_operator",
    }:
        raise ScoringError(
            "target_schema_mismatch",
            "machine contract lacks a target-blind method-selection policy",
        )
    _verified_receipt(
        _field(
            contract,
            "source_method_freeze_receipt",
            "target.machine_scoring_contract",
        ),
        "target.machine_scoring_contract.source_method_freeze_receipt",
    )
    decision = _mapping(
        _field(contract, "decision_policy", "target.machine_scoring_contract"),
        "target.machine_scoring_contract.decision_policy",
    )
    if decision.get("outcomes") != ["COMPATIBLE", "FAIL", "INCONCLUSIVE"]:
        raise ScoringError(
            "target_schema_mismatch",
            "decision policy must distinguish COMPATIBLE, FAIL, and INCONCLUSIVE",
        )
    confidence_raw = _field(decision, "measurement_confidence_level", "decision_policy")
    multiplier_raw = _field(decision, "measurement_sigma_multiplier", "decision_policy")
    if not isinstance(confidence_raw, str) or not isinstance(multiplier_raw, str):
        raise ScoringError(
            "target_schema_mismatch",
            "confidence level and sigma multiplier must be decimal strings",
        )
    confidence = _decimal(
        confidence_raw,
        "target.machine_scoring_contract.decision_policy.measurement_confidence_level",
    )
    multiplier = _decimal(
        multiplier_raw,
        "target.machine_scoring_contract.decision_policy.measurement_sigma_multiplier",
    )
    if not 0 < confidence < 1 or multiplier <= 0:
        raise ScoringError(
            "target_schema_mismatch",
            "measurement confidence and sigma multiplier are invalid",
        )
    if decision.get("combined_uncertainty_rule") != (
        "interval_minkowski_sum_with_frozen_confidence_multiplier"
    ):
        raise ScoringError(
            "target_schema_mismatch", "combined uncertainty rule is not recognized"
        )
    maximum_width = _mapping(
        _field(decision, "maximum_prediction_width", "decision_policy"),
        "target.machine_scoring_contract.decision_policy.maximum_prediction_width",
    )
    if maximum_width.get("coordinate") != "inverse_alpha":
        raise ScoringError(
            "target_schema_mismatch", "maximum prediction width has wrong coordinate"
        )
    width_raw = _field(maximum_width, "value", "maximum_prediction_width")
    if not isinstance(width_raw, str):
        raise ScoringError(
            "target_schema_mismatch",
            "maximum prediction width must be a decimal string",
        )
    width_value = _decimal(
        width_raw,
        "target.machine_scoring_contract.decision_policy.maximum_prediction_width.value",
    )
    if width_value <= 0:
        raise ScoringError(
            "target_schema_mismatch", "maximum prediction width must be positive"
        )
    expected_predicates = {
        "FAIL": "completed_map_interval_disjoint_from_combined_measurement_interval",
        "COMPATIBLE": "eligible_narrow_completed_map_interval_overlaps_combined_measurement_interval",
        "INCONCLUSIVE": "all_other_cases_including_excess_width_or_open_gate",
    }
    if decision.get("predicates") != expected_predicates:
        raise ScoringError(
            "target_schema_mismatch", "decision predicates are incomplete or ambiguous"
        )

    domain = _mapping(
        _field(contract, "p_domain", "target.machine_scoring_contract"),
        "target.machine_scoring_contract.p_domain",
    )
    lo = _decimal_string(
        _field(domain, "lo", "target.machine_scoring_contract.p_domain"),
        "target.machine_scoring_contract.p_domain.lo",
    )
    hi = _decimal_string(
        _field(domain, "hi", "target.machine_scoring_contract.p_domain"),
        "target.machine_scoring_contract.p_domain.hi",
    )
    if lo >= hi or domain.get("kind") != "registered_P_basin":
        raise ScoringError(
            "target_schema_mismatch",
            "registered P domain must be a non-singleton basin",
        )

    schema = _mapping(
        _field(contract, "coordinate_schema", "target.machine_scoring_contract"),
        "target.machine_scoring_contract.coordinate_schema",
    )
    if set(schema) != set(COORDINATES):
        raise ScoringError(
            "target_coordinate_mismatch",
            "machine contract must distinguish total, residual, S_QEW, and S_hadronic",
        )
    for name, (kind, units, role) in REQUIRED_COORDINATE_TYPES.items():
        coordinate = _mapping(
            schema[name], f"target.machine_scoring_contract.coordinate_schema.{name}"
        )
        if (
            coordinate.get("kind") != kind
            or coordinate.get("units") != units
            or coordinate.get("scoring_role") != role
            or not isinstance(coordinate.get("artifact_path"), str)
        ):
            raise ScoringError(
                "target_coordinate_mismatch", f"machine contract misclassifies {name}"
            )
    required_receipts = contract.get("required_receipt_fields")
    if (
        not isinstance(required_receipts, list)
        or set(required_receipts) != REQUIRED_RECEIPT_FIELDS
    ):
        raise ScoringError(
            "target_schema_mismatch",
            "machine contract must require the full provenance receipt set",
        )
    return contract


def _require_activated_target(target: dict[str, Any]) -> dict[str, Any]:
    if target.get("id") == CURRENT_INACTIVE_CORRECTIVE_ID:
        raise ScoringError(
            "target_not_activated",
            "the 2026-07-16 v3 file is an erratum scaffold and can never be activated; issue a detached, externally anchored successor",
        )
    activation = _mapping(
        _field(target, "activation_requirements", "target"),
        "target.activation_requirements",
    )
    inactive = (
        "not_scorable" in str(target.get("registration_status"))
        or target.get("frozen_utc") is None
        or target.get("promotion_or_falsification_allowed") is not True
    )
    if inactive:
        raise ScoringError(
            "target_not_activated",
            "corrective target is not externally frozen or activated; no payload may score",
        )
    required_activation_fields = (
        "governs_payloads_started_after",
        "external_timestamp_receipt",
        "first_eligible_payload_commit_definition",
        "canonical_artifact_digest",
        "detached_activation_manifest",
    )
    for field in required_activation_fields:
        if not activation.get(field):
            raise ScoringError(
                "target_not_activated", f"activation receipt {field} is missing"
            )
    start_after = _utc_timestamp(
        activation["governs_payloads_started_after"],
        "target.activation_requirements.governs_payloads_started_after",
    )
    frozen = _utc_timestamp(target["frozen_utc"], "target.frozen_utc")
    if frozen >= start_after:
        raise ScoringError(
            "target_not_activated",
            "eligible payload work must start strictly after the target freeze",
        )
    commit_definition = _mapping(
        activation["first_eligible_payload_commit_definition"],
        "target.activation_requirements.first_eligible_payload_commit_definition",
    )
    if commit_definition.get("kind") != "exact_source_method_git_commit":
        raise ScoringError(
            "target_not_activated",
            "eligible source method must be bound to an exact Git commit",
        )
    eligible_commit = _nonempty_string(
        _field(
            commit_definition,
            "object_id",
            "target.activation_requirements.first_eligible_payload_commit_definition",
        ),
        "target.activation_requirements.first_eligible_payload_commit_definition.object_id",
    )
    if (
        len(eligible_commit) not in {40, 64}
        or any(character not in string.hexdigits for character in eligible_commit)
        or set(eligible_commit) == {"0"}
    ):
        raise ScoringError(
            "target_not_activated", "eligible source method object id is invalid"
        )
    _verified_receipt(
        _field(
            commit_definition,
            "verification_receipt",
            "target.activation_requirements.first_eligible_payload_commit_definition",
        ),
        "target.activation_requirements.first_eligible_payload_commit_definition.verification_receipt",
    )
    canonical_digest = _nonzero_sha256(
        activation["canonical_artifact_digest"],
        "target.activation_requirements.canonical_artifact_digest",
    )
    timestamp_receipt = _verified_receipt(
        activation["external_timestamp_receipt"],
        "target.activation_requirements.external_timestamp_receipt",
    )
    receipt_time = _utc_timestamp(
        _field(
            timestamp_receipt,
            "timestamp_utc",
            "target.activation_requirements.external_timestamp_receipt",
        ),
        "target.activation_requirements.external_timestamp_receipt.timestamp_utc",
    )
    if timestamp_receipt.get("target_sha256") != canonical_digest:
        raise ScoringError(
            "target_not_activated",
            "external timestamp receipt does not bind the canonical target digest",
        )
    if not frozen <= receipt_time < start_after:
        raise ScoringError(
            "target_not_activated",
            "external timestamp must bind the frozen target before eligible payload work",
        )
    detached = _verified_receipt(
        activation["detached_activation_manifest"],
        "target.activation_requirements.detached_activation_manifest",
    )
    if detached.get("target_sha256") != canonical_digest:
        raise ScoringError(
            "target_not_activated",
            "detached activation manifest does not bind the canonical target digest",
        )
    contract = _mapping(
        _field(target, "machine_scoring_contract", "target"),
        "target.machine_scoring_contract",
    )
    validated = dict(_validate_machine_contract(contract))
    validated["_activation_cutoff_utc"] = start_after.isoformat()
    validated["_canonical_target_digest"] = canonical_digest
    validated["_eligible_source_commit"] = eligible_commit.lower()
    return validated


def _validate_artifact(
    artifact: dict[str, Any], contract: dict[str, Any]
) -> dict[str, dict[str, Decimal]]:
    if artifact.get("artifact") != EMISSION_ARTIFACT:
        raise ScoringError(
            "artifact_schema_mismatch",
            "legacy or unknown emission artifact; historical V1 is never scoreable",
        )
    if artifact.get("schema_version") != EMISSION_SCHEMA_VERSION:
        raise ScoringError(
            "artifact_schema_mismatch", "unexpected emission schema_version"
        )

    recorded_hash = artifact.get("content_sha256")
    _nonzero_sha256(recorded_hash, "artifact.content_sha256")
    if recorded_hash != artifact_content_sha256(artifact):
        raise ScoringError(
            "artifact_hash_mismatch", "content_sha256 verification failed"
        )
    if artifact.get("promotion_allowed") is not False:
        raise ScoringError(
            "artifact_schema_mismatch", "a source emitter cannot authorize promotion"
        )
    if artifact.get("target_or_measurement_inputs_used_in_computation") is not False:
        raise ScoringError(
            "artifact_schema_mismatch",
            "artifact does not attest target-free production",
        )
    scheme = _mapping(_field(artifact, "scheme", "artifact"), "artifact.scheme")
    if scheme.get("same_subtraction_as_a0") is not True:
        raise ScoringError(
            "artifact_schema_mismatch",
            "artifact does not close the same-subtraction scheme lock",
        )
    if artifact.get("payload_object") != contract["payload_object"]:
        raise ScoringError(
            "artifact_not_certified",
            "payload is not a target-blind function or certified P-domain enclosure",
        )
    if artifact.get("source_family_id") != contract["source_family_id"]:
        raise ScoringError("artifact_schema_mismatch", "source_family_id mismatch")
    if scheme.get("scheme_id") != contract["scheme_id"]:
        raise ScoringError("artifact_schema_mismatch", "scheme_id mismatch")
    if artifact.get("current_definition_id") != contract["current_definition_id"]:
        raise ScoringError("artifact_schema_mismatch", "current_definition_id mismatch")

    target_domain = _mapping(
        contract["p_domain"], "target.machine_scoring_contract.p_domain"
    )
    artifact_domain = _mapping(
        _field(artifact, "p_domain", "artifact"), "artifact.p_domain"
    )
    artifact_lo = _decimal_string(
        _field(artifact_domain, "lo", "artifact.p_domain"), "artifact.p_domain.lo"
    )
    artifact_hi = _decimal_string(
        _field(artifact_domain, "hi", "artifact.p_domain"), "artifact.p_domain.hi"
    )
    target_lo = _decimal_string(
        target_domain["lo"], "target.machine_scoring_contract.p_domain.lo"
    )
    target_hi = _decimal_string(
        target_domain["hi"], "target.machine_scoring_contract.p_domain.hi"
    )
    if (
        artifact_domain.get("kind") != "registered_P_basin"
        or artifact_lo != target_lo
        or artifact_hi != target_hi
    ):
        raise ScoringError(
            "evaluation_point_mismatch",
            "artifact P domain differs from the registered target basin",
        )

    actual_schema = _mapping(
        _field(artifact, "coordinate_schema", "artifact"),
        "artifact.coordinate_schema",
    )
    expected_schema = _mapping(
        contract["coordinate_schema"],
        "target.machine_scoring_contract.coordinate_schema",
    )
    if set(actual_schema) != set(COORDINATES):
        raise ScoringError(
            "coordinate_schema_mismatch",
            "artifact must distinguish total, residual, S_QEW, and S_hadronic",
        )
    for name in COORDINATES:
        actual = _mapping(actual_schema[name], f"artifact.coordinate_schema.{name}")
        expected = _mapping(
            expected_schema[name],
            f"target.machine_scoring_contract.coordinate_schema.{name}",
        )
        for field in ("kind", "units", "artifact_path", "scoring_role"):
            if actual.get(field) != expected.get(field):
                raise ScoringError(
                    "coordinate_schema_mismatch",
                    f"{name}.{field} does not match the target contract",
                )

    # The source emitter is allowed to publish a finite-float diagnostic
    # envelope, but only a certified artifact may claim Decimal enclosure
    # semantics. Check the generic P and coordinate schema first so mismatches
    # retain their precise fail-closed codes, then stop uncertified artifacts
    # before parsing certificate-only interval endpoints or receipts.
    certification = _mapping(
        _field(artifact, "certification", "artifact"), "artifact.certification"
    )
    if certification.get("status") != "certified":
        raise ScoringError(
            "artifact_not_certified", "source object lacks a certified enclosure"
        )
    if certification.get("sampled_grid_extrema_interval_certificate") is not False:
        raise ScoringError(
            "artifact_not_certified", "sampled extrema cannot certify an interval"
        )

    _verified_receipt(
        _field(scheme, "scheme_lock_receipt", "artifact.scheme"),
        "artifact.scheme.scheme_lock_receipt",
    )

    intervals: dict[str, dict[str, Decimal]] = {}
    for name in COORDINATES:
        expected = _mapping(
            expected_schema[name],
            f"target.machine_scoring_contract.coordinate_schema.{name}",
        )
        intervals[name] = _validate_interval(
            _at_path(artifact, expected["artifact_path"]),
            f"artifact.{expected['artifact_path']}",
        )

    if certification.get("delta_EW_gate") != "closed":
        raise ScoringError(
            "open_delta_EW_gate",
            "Delta_EW is open, so completed-map scoring is unavailable",
        )
    if certification.get("delta_EW_theorem_status") != "closed":
        raise ScoringError(
            "open_delta_EW_gate",
            "a literal gate label cannot replace a closed Delta_EW theorem",
        )
    _verified_receipt(
        _field(certification, "delta_EW_source_receipt", "artifact.certification"),
        "artifact.certification.delta_EW_source_receipt",
    )
    _validate_interval(
        _field(certification, "delta_EW_interval", "artifact.certification"),
        "artifact.certification.delta_EW_interval",
    )
    numerical_error = _validate_interval(
        _field(certification, "numerical_error_interval", "artifact.certification"),
        "artifact.certification.numerical_error_interval",
    )
    theory_error = _validate_interval(
        _field(certification, "theory_error_interval", "artifact.certification"),
        "artifact.certification.theory_error_interval",
    )
    for name, interval in (
        ("numerical", numerical_error),
        ("theory", theory_error),
    ):
        if not interval["lo"] <= 0 <= interval["hi"]:
            raise ScoringError(
                "artifact_not_certified",
                f"{name} error interval must contain zero",
            )
    derivative_raw = _field(
        certification,
        "derivative_or_lipschitz_bound_over_P_domain",
        "artifact.certification",
    )
    if not isinstance(derivative_raw, str):
        raise ScoringError(
            "artifact_not_certified",
            "derivative or Lipschitz bound must be a decimal string",
        )
    derivative_bound = _decimal(
        derivative_raw,
        "artifact.certification.derivative_or_lipschitz_bound_over_P_domain",
    )
    if derivative_bound < 0:
        raise ScoringError(
            "artifact_not_certified",
            "derivative or Lipschitz bound must be non-negative",
        )
    contraction = _mapping(
        _field(
            certification,
            "completed_map_contraction_bounds",
            "artifact.certification",
        ),
        "artifact.certification.completed_map_contraction_bounds",
    )
    if set(contraction) != {"CL-1", "CL-2"}:
        raise ScoringError(
            "artifact_not_certified",
            "both completed maps require independent contraction bounds",
        )
    for map_name, value in contraction.items():
        bound = _decimal(
            value, f"artifact.certification.completed_map_contraction_bounds.{map_name}"
        )
        if not 0 <= bound < 1:
            raise ScoringError(
                "artifact_not_certified",
                f"{map_name} contraction bound must be in [0,1)",
            )
    _verified_receipt(
        _field(
            certification,
            "component_enclosure_identity_receipt",
            "artifact.certification",
        ),
        "artifact.certification.component_enclosure_identity_receipt",
    )

    receipts = _mapping(_field(artifact, "receipts", "artifact"), "artifact.receipts")
    if set(receipts) != set(contract["required_receipt_fields"]):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "artifact provenance receipt set is incomplete",
        )
    for receipt_field in (
        "target_free_dependency_dag",
        "forbidden_input_scan",
        "environment_lock",
        "external_timestamp",
    ):
        _verified_receipt(receipts[receipt_field], f"artifact.receipts.{receipt_field}")
    for hash_field in (
        "source_tree_sha256",
        "canonical_json_sha256",
    ):
        _nonzero_sha256(receipts[hash_field], f"artifact.receipts.{hash_field}")
    if receipts["canonical_json_sha256"].lower() != recorded_hash:
        raise ScoringError(
            "artifact_provenance_mismatch",
            "canonical JSON receipt does not equal the embedded artifact digest",
        )
    commit = _nonempty_string(
        receipts["source_commit"], "artifact.receipts.source_commit"
    )
    if (
        len(commit) not in {40, 64}
        or any(character not in string.hexdigits for character in commit)
        or set(commit) == {"0"}
    ):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "source_commit must be a nonzero full Git object id",
        )
    if commit.lower() != contract.get("_eligible_source_commit"):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "artifact source commit is not the target-eligible frozen method commit",
        )
    argv = receipts["generator_argv"]
    if (
        not isinstance(argv, list)
        or not argv
        or any(not isinstance(item, str) or not item for item in argv)
    ):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "generator_argv must be a nonempty string array",
        )
    _nonempty_string(receipts["python_version"], "artifact.receipts.python_version")
    versions = _mapping(
        receipts["numeric_library_versions"],
        "artifact.receipts.numeric_library_versions",
    )
    if not versions or any(
        not isinstance(name, str)
        or not name
        or not isinstance(version, str)
        or not version
        for name, version in versions.items()
    ):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "numeric library versions are incomplete",
        )
    precision = _mapping(
        receipts["precision_and_cutoffs"],
        "artifact.receipts.precision_and_cutoffs",
    )
    if not precision or any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in precision.values()
    ):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "precision and cutoffs must be positive integers",
        )
    dependencies = _mapping(
        receipts["executable_dependency_sha256"],
        "artifact.receipts.executable_dependency_sha256",
    )
    if not dependencies:
        raise ScoringError(
            "artifact_provenance_mismatch", "dependency hash receipt is empty"
        )
    for name, digest in dependencies.items():
        _nonempty_string(name, "artifact.receipts.executable_dependency_sha256 key")
        _nonzero_sha256(
            digest, f"artifact.receipts.executable_dependency_sha256.{name}"
        )
    timestamp_receipt = _mapping(
        receipts["external_timestamp"], "artifact.receipts.external_timestamp"
    )
    if timestamp_receipt.get("source_tree_sha256") != receipts["source_tree_sha256"]:
        raise ScoringError(
            "artifact_provenance_mismatch",
            "external timestamp does not bind the emitted source tree",
        )
    dag_receipt = _mapping(
        receipts["target_free_dependency_dag"],
        "artifact.receipts.target_free_dependency_dag",
    )
    if timestamp_receipt.get("dependency_dag_sha256") != dag_receipt.get("sha256"):
        raise ScoringError(
            "artifact_provenance_mismatch",
            "external timestamp does not bind the target-free dependency DAG",
        )
    timestamp = _utc_timestamp(
        _field(
            timestamp_receipt,
            "timestamp_utc",
            "artifact.receipts.external_timestamp",
        ),
        "artifact.receipts.external_timestamp.timestamp_utc",
    )
    work_started = _utc_timestamp(
        receipts["payload_work_started_utc"],
        "artifact.receipts.payload_work_started_utc",
    )
    activation_cutoff = _utc_timestamp(
        _field(contract, "_activation_cutoff_utc", "target.machine_scoring_contract"),
        "target.machine_scoring_contract._activation_cutoff_utc",
    )
    if work_started < activation_cutoff:
        raise ScoringError(
            "artifact_provenance_mismatch",
            "payload work began before this target's eligibility cutoff",
        )
    if timestamp >= work_started:
        raise ScoringError(
            "artifact_provenance_mismatch",
            "source/dependency timestamp must strictly precede payload work",
        )
    return intervals


def score_artifact(artifact: dict[str, Any], target: dict[str, Any]) -> dict[str, Any]:
    """Validate a sealed source artifact, then fail closed before scalar scoring.

    The function raises :class:`ScoringError` for every non-evaluable state.
    A future completed-map solver must replace the final error; adding scalar
    point containment here would violate the corrective target.
    """
    artifact = _mapping(artifact, "artifact")
    target = _mapping(target, "target")
    _validate_corrective_target(target)
    contract = _require_activated_target(target)
    _validate_artifact(artifact, contract)
    raise ScoringError(
        "completed_map_solver_unavailable",
        "validated payload still requires independent interval solves of CL-1 and CL-2",
    )


def _load_json(path: Path, role: str) -> dict[str, Any]:
    def reject_nonfinite(value: str) -> None:
        raise ValueError(f"non-finite JSON number {value}")

    try:
        value = json.loads(
            path.read_text(encoding="utf-8"), parse_constant=reject_nonfinite
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        raise ScoringError(f"{role}_read_error", f"cannot read {role}: {exc}") from exc
    return _mapping(value, role)


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    try:
        artifact = _load_json(args.artifact, "artifact")
        target = _load_json(args.target, "target")
        score_artifact(artifact, target)
        raise AssertionError("score_artifact must return a verdict or fail closed")
    except ScoringError as exc:
        result = {
            "artifact": "oph_ward_projected_payload_score_attempt",
            "schema_version": 2,
            "status": "NOT_EVALUABLE",
            "error_code": exc.code,
            "error": exc.message,
            "map_results": {"CL-1": "NOT_EVALUABLE", "CL-2": "NOT_EVALUABLE"},
            "closure_allowed": False,
            "falsification_allowed": False,
            "promotion_allowed": False,
        }
        if args.artifact.exists():
            result["source_file_sha256"] = _file_sha256(args.artifact)
        if args.target.exists():
            result["target_file_sha256"] = _file_sha256(args.target)

    text = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        output = args.output.resolve()
        if output in {args.artifact.resolve(), args.target.resolve()}:
            print("refusing to overwrite artifact or target", file=sys.stderr)
            return 2
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
