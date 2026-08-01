#!/usr/bin/env python3
"""Issue #649: finite-cut survival versus global-capacity attachment.

The consumed target-free projection of the finite receipt types ``q`` as the
expectation of one normalized class projector and ``1-q`` as the corresponding
one-cut survival reading.  It does not supply an action of that reading on a
global-capacity object.  Independent replay checks the projection against the
full parent receipt outside the producer cone.

This certificate constructs two exact global completions of the same local
datum.  Both are positive, permutation invariant, and compose on disconnected
cuts and under regrouping:

* ``neutral`` leaves global capacity unchanged;
* ``multiplicative`` acts by ``(1-q)^m`` for ``m`` disconnected cuts.

They disagree for every positive capacity and every ``0 < q < 1``.  Even a
monoid-action requirement therefore does not select the attachment.  A second
exact construction exhibits three positive compositional actions on the same
source datum when ``0 < q < 1/6``: no capacity action, one-class blocking, and
six-class-total blocking.  Their one-cut effects are pairwise distinct.  The
declared finite source therefore does not select among those three actions.

No cosmological quantity, measured target, horizon area, or desired capacity
is read.  The result is the bounded non-identifiability exit declared by
issue #649; it emits no numeric cosmic value.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCE_PROJECTION_PATH = (
    HERE / "manifests" / "global_capacity_attachment_source_projection.json"
)
OUTPUT_PATH = HERE / "runtime" / "global_capacity_attachment_receipt.json"

SCHEMA = "oph.global_capacity_attachment.v1"
STATUS = "BOUNDED_NONIDENTIFIABLE_GLOBAL_COMPOSITION"
ISSUE = 649
FINITE_SCHEMA = "oph.edge_center_clock_certificate.v3"
FINITE_STATUS = "conditional_edge_center_arithmetic_with_open_source_and_clock_gates"
SOURCE_PROJECTION_SCHEMA = "oph.global_capacity_attachment_source_projection.v1"


class CertificateError(ValueError):
    """Fail-closed error with a stable mutation-test code."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise CertificateError(code, message)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged_sha256(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError("SOURCE_READ", f"cannot read {path}: {exc}") from exc
    require(isinstance(value, dict), "SOURCE_TYPE", f"{path} is not an object")
    return value


def _pin(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": len(raw),
        "sha256": tagged_sha256(raw),
    }


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _fraction(value: Any, label: str) -> Fraction:
    require(isinstance(value, str), "PAYLOAD_TYPE", f"{label} must be a string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError("ARITHMETIC", f"{label} is not rational") from exc
    return result


def _validate_source_projection() -> None:
    finite = _load(SOURCE_PROJECTION_PATH)
    require(
        finite.get("schema") == SOURCE_PROJECTION_SCHEMA,
        "SOURCE_STATUS",
        "projection schema drift",
    )
    require(
        finite.get("finite_source_schema") == FINITE_SCHEMA,
        "SOURCE_STATUS",
        "finite-source schema drift",
    )
    require(
        finite.get("finite_source_status") == FINITE_STATUS,
        "SOURCE_STATUS",
        "finite-source status drift",
    )
    statement = finite.get("presence_statement")
    require(
        isinstance(statement, str)
        and "finite one-step survival is the presence value 1 - P/24" in statement,
        "SOURCE_SEMANTICS",
        "finite source fails to type one-step presence survival",
    )
    require(
        finite.get("operational_clock_bound") is False,
        "PROMOTION",
        "finite source unexpectedly promotes a physical clock",
    )
    reserve_input = finite.get("reserve_trace_branch", {})
    require(
        reserve_input.get("binding_theorem_present") is False
        and reserve_input.get("derived_from_collar_counts") is False,
        "SOURCE_SEMANTICS",
        "finite source fails to expose the declared reserve branch boundary",
    )
    p_interval = finite.get("P_certified_enclosure", {})
    p_lo = _fraction(p_interval.get("lo"), "P_certified_enclosure.lo")
    p_hi = _fraction(p_interval.get("hi"), "P_certified_enclosure.hi")
    require(
        0 < p_lo <= p_hi < 4,
        "DOMAIN",
        "finite source does not certify 0 < q=P/24 < 1/6",
    )


def _neutral(_survival: Fraction, _cuts: int, capacity: Fraction) -> Fraction:
    return capacity


def _multiplicative(survival: Fraction, cuts: int, capacity: Fraction) -> Fraction:
    return survival**cuts * capacity


def _factor_action(factor: Fraction, cuts: int, capacity: Fraction) -> Fraction:
    return factor**cuts * capacity


def _composition_and_refinement(
    factor: Fraction, capacity: Fraction
) -> tuple[bool, list[dict[str, Any]]]:
    composition_checks = []
    for left in range(4):
        for right in range(4 - left):
            direct = _factor_action(factor, left + right, capacity)
            nested = _factor_action(
                factor, left, _factor_action(factor, right, capacity)
            )
            composition_checks.append(direct == nested)
    refinement_partitions = ([3], [1, 2], [2, 1], [1, 1, 1])
    refinement_checks = []
    for pieces in refinement_partitions:
        refined = capacity
        for piece in reversed(pieces):
            refined = _factor_action(factor, piece, refined)
        coarse = _factor_action(factor, sum(pieces), capacity)
        refinement_checks.append(
            {
                "pieces": list(pieces),
                "coarse_cut_count": sum(pieces),
                "coarse_value": _fraction_text(coarse),
                "iterated_refined_value": _fraction_text(refined),
                "commutes": refined == coarse,
            }
        )
    return all(composition_checks), refinement_checks


def _completion_rows(q: Fraction, capacity: Fraction) -> list[dict[str, Any]]:
    survival = 1 - q
    rows = []
    for name, action in (("neutral", _neutral), ("multiplicative", _multiplicative)):
        values = [action(survival, cuts, capacity) for cuts in range(4)]
        composition_holds, refinement_checks = _composition_and_refinement(
            Fraction(1) if name == "neutral" else survival, capacity
        )
        rows.append(
            {
                "completion": name,
                "formula": (
                    "G_m(C)=C" if name == "neutral" else "G_m(C)=(1-q)^m C"
                ),
                "values_m_0_to_3": [_fraction_text(value) for value in values],
                "identity_at_zero_cuts": values[0] == capacity,
                "positive": all(value > 0 for value in values),
                "permutation_invariant": True,
                "disconnected_cut_composition": composition_holds,
                "cut_count_source_object": "free commutative monoid on one abstract cut",
                "refinement_model": (
                    "a coarse cut count is replaced by a finite partition with "
                    "the same total count"
                ),
                "refinement_partition_checks": refinement_checks,
                "refinement_regrouping_invariant": all(
                    row["commutes"] for row in refinement_checks
                ),
            }
        )
    return rows


def _blocked_event_completion_rows(
    q: Fraction, capacity: Fraction
) -> list[dict[str, Any]]:
    factors = (
        ("no_capacity_action", Fraction(1)),
        ("one_class_projector", 1 - q),
        ("six_class_total_projector", 1 - 6 * q),
    )
    rows: list[dict[str, Any]] = []
    for name, factor in factors:
        values = [_factor_action(factor, cuts, capacity) for cuts in range(4)]
        composition_holds, refinement_checks = _composition_and_refinement(
            factor, capacity
        )
        rows.append(
            {
                "completion": name,
                "formula": f"G_m(C)=({_fraction_text(factor)})^m C",
                "factor": _fraction_text(factor),
                "values_m_0_to_3": [_fraction_text(value) for value in values],
                "identity_at_zero_cuts": values[0] == capacity,
                "positive": all(value > 0 for value in values),
                "permutation_invariant": True,
                "disconnected_cut_composition": composition_holds,
                "refinement_partition_checks": refinement_checks,
                "refinement_regrouping_invariant": all(
                    check["commutes"] for check in refinement_checks
                ),
            }
        )
    return rows


def _witness_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for q in (Fraction(1, 60), Fraction(1, 24), Fraction(1, 12)):
        capacity = Fraction(120)
        survival = 1 - q
        neutral = _neutral(survival, 1, capacity)
        multiplicative = _multiplicative(survival, 1, capacity)
        one_class = 1 - q
        six_class = 1 - 6 * q
        no_action = Fraction(1)
        blocked_completions = _blocked_event_completion_rows(q, capacity)
        rows.append(
            {
                "q": _fraction_text(q),
                "capacity": _fraction_text(capacity),
                "survival": _fraction_text(survival),
                "completions": _completion_rows(q, capacity),
                "different_at_one_cut": neutral != multiplicative,
                "blocked_event_factors": {
                    "one_class": _fraction_text(one_class),
                    "six_mutually_exclusive_classes": _fraction_text(six_class),
                    "no_capacity_action": _fraction_text(no_action),
                },
                "blocked_event_factors_pairwise_distinct": (
                    len({one_class, six_class, no_action}) == 3
                ),
                "blocked_event_completions": blocked_completions,
                "blocked_event_actions_pairwise_distinct_at_one_cut": (
                    len(
                        {
                            Fraction(item["values_m_0_to_3"][1])
                            for item in blocked_completions
                        }
                    )
                    == 3
                ),
            }
        )
    return rows


def build() -> dict[str, Any]:
    """Build the deterministic generic countermodel receipt."""
    _validate_source_projection()
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": STATUS,
        "source_pin": _pin(SOURCE_PROJECTION_PATH),
        "producer_cone": {
            "consumed_source_fields": [
                "finite_source_schema",
                "finite_source_status",
                "presence_statement",
                "reserve_trace_branch.binding_theorem_present",
                "reserve_trace_branch.derived_from_collar_counts",
                "operational_clock_bound",
                "P_certified_enclosure",
            ],
            "cosmological_target_payload_used_in_countermodel": False,
            "desired_capacity_used_in_countermodel": False,
            "completion_selected_from_target": False,
            "full_parent_receipt_not_in_producer_cone": True,
            "projection_checked_against_parent_by_independent_replay": True,
        },
        "source_boundary": {
            "established": (
                "on the declared reserve-trace branch, one normalized class "
                "expectation q=P/24 with 0<q<1/6 and the conditional local "
                "survival 1-q"
            ),
            "not_established": [
                "an action of local survival on global capacity",
                "a disconnected-cut product law acting on capacity",
                "identification of the physical blocked event",
                "a Poisson or projective-limit carrier",
            ],
        },
        "typed_countermodel_class": {
            "connected_cut_source_object": (
                "the free commutative cut-count monoid (Nat,+,0) carrying the "
                "same normalized local reading q"
            ),
            "global_capacity_object": (
                "the positive real ray; the executable receipt serializes an "
                "exact rational witness subfamily"
            ),
            "refinement_relation": (
                "finite partitions of one coarse cut count whose entries have "
                "the same total"
            ),
            "scope": (
                "an abstract source-indistinguishability class; no physical "
                "universe carrier or cut-to-edge identification is asserted, "
                "and the consumed reserve-trace branch remains conditional"
            ),
        },
        "generic_theorem": {
            "domain": "0 < q < 1 and C > 0",
            "neutral_completion": "G_m(C)=C",
            "multiplicative_completion": "G_m(C)=(1-q)^m C",
            "both_are_positive_monoid_actions": True,
            "both_are_permutation_and_refinement_regrouping_invariant": True,
            "they_disagree_at_one_cut": True,
            "conclusion": (
                "the local finite-cut datum plus composition and regrouping laws "
                "does not select a global-capacity attachment"
            ),
            "lean_module": "Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean",
        },
        "exact_witnesses": _witness_rows(),
        "blocked_event_verdict": {
            "one_class_selected": False,
            "six_class_total_selected": False,
            "no_capacity_action_selected": False,
            "declared_completion_set": [
                "no_capacity_action",
                "one_class_projector",
                "six_class_total_projector",
            ],
            "same_source_datum": "one normalized local class expectation q",
            "status": "SOURCE_DOES_NOT_SELECT_AMONG_DECLARED_BLOCKED_EVENT_ACTIONS",
        },
        "composition_boundary": {
            "disconnected_cut_composition_classified": True,
            "per_edge_multiplication_authorized": False,
            "reason": (
                "a factor per edge requires an independently constructed product "
                "decomposition and an edge-to-capacity action"
            ),
        },
        "comparison_boundary": {
            "cosmological_payload_read": False,
            "numeric_cosmic_capacity_emitted": False,
            "forecast_or_comparison_permitted": False,
            "finite_presence_branch_selected": False,
            "poisson_branch_selected": False,
        },
        "consumer": {
            "issue": 648,
            "effect": (
                "consume as a bounded nonidentifiability verdict; no reserve branch "
                "may be promoted without a stronger source law selecting its global action"
            ),
        },
        "reopen_condition": (
            "a source-derived typed action of the connected-cut survival object on "
            "global capacity, including the physical blocked-event projector and "
            "its disconnected-cut/refinement composition law"
        ),
    }
    payload["receipt_sha256"] = tagged_sha256(canonical_json_bytes(payload))
    validate(payload, verify_source=True)
    return payload


def validate(payload: Mapping[str, Any], *, verify_source: bool = True) -> None:
    require(isinstance(payload, Mapping), "PAYLOAD_TYPE", "payload is not an object")
    require(payload.get("schema") == SCHEMA, "SCHEMA", "schema drift")
    require(payload.get("issue") == ISSUE, "SCHEMA", "issue drift")
    require(payload.get("status") == STATUS, "STATUS", "status drift")
    require(
        payload.get("source_pin") == _pin(SOURCE_PROJECTION_PATH),
        "SOURCE_PIN",
        "pin drift",
    )
    if verify_source:
        _validate_source_projection()

    expected_rows = _witness_rows()
    require(payload.get("exact_witnesses") == expected_rows, "ARITHMETIC", "witness drift")
    for row in expected_rows:
        q = _fraction(row["q"], "q")
        require(0 < q < Fraction(1, 6), "DOMAIN", "witness q outside menu domain")
        require(row["different_at_one_cut"] is True, "NONIDENTIFIABILITY", "completions agree")
        require(
            row["blocked_event_factors_pairwise_distinct"] is True,
            "NONIDENTIFIABILITY",
            "blocked-event menu collapsed",
        )
        require(
            row["blocked_event_actions_pairwise_distinct_at_one_cut"] is True,
            "NONIDENTIFIABILITY",
            "blocked-event actions agree after one cut",
        )
        for completion in row["completions"]:
            require(
                completion["identity_at_zero_cuts"] is True
                and completion["positive"] is True
                and completion["permutation_invariant"] is True
                and completion["disconnected_cut_composition"] is True
                and completion["refinement_regrouping_invariant"] is True,
                "COMPOSITION",
                "completion law failed",
            )
            refinement_checks = completion.get("refinement_partition_checks")
            require(
                isinstance(refinement_checks, list)
                and len(refinement_checks) == 4
                and all(item.get("commutes") is True for item in refinement_checks),
                "COMPOSITION",
                "refinement partition check failed",
            )
        blocked_completions = row.get("blocked_event_completions")
        require(
            isinstance(blocked_completions, list)
            and [item.get("completion") for item in blocked_completions]
            == [
                "no_capacity_action",
                "one_class_projector",
                "six_class_total_projector",
            ],
            "NONIDENTIFIABILITY",
            "blocked-event completion set drift",
        )
        for completion in blocked_completions:
            require(
                completion.get("identity_at_zero_cuts") is True
                and completion.get("positive") is True
                and completion.get("permutation_invariant") is True
                and completion.get("disconnected_cut_composition") is True
                and completion.get("refinement_regrouping_invariant") is True,
                "COMPOSITION",
                "blocked-event completion law failed",
            )

    theorem = payload.get("generic_theorem")
    require(isinstance(theorem, Mapping), "PAYLOAD_TYPE", "theorem missing")
    require(
        theorem.get("both_are_positive_monoid_actions") is True
        and theorem.get("both_are_permutation_and_refinement_regrouping_invariant") is True
        and theorem.get("they_disagree_at_one_cut") is True,
        "NONIDENTIFIABILITY",
        "generic theorem flags drift",
    )
    require(
        theorem.get("domain") == "0 < q < 1 and C > 0"
        and theorem.get("neutral_completion") == "G_m(C)=C"
        and theorem.get("multiplicative_completion") == "G_m(C)=(1-q)^m C"
        and theorem.get("lean_module")
        == "Lean/ObserverPatchHolography/GlobalCapacityAttachment.lean",
        "NONIDENTIFIABILITY",
        "generic theorem statement drift",
    )
    source_boundary = payload.get("source_boundary")
    require(
        isinstance(source_boundary, Mapping)
        and source_boundary.get("established")
        == (
            "on the declared reserve-trace branch, one normalized class "
            "expectation q=P/24 with 0<q<1/6 and the conditional local "
            "survival 1-q"
        )
        and source_boundary.get("not_established")
        == [
            "an action of local survival on global capacity",
            "a disconnected-cut product law acting on capacity",
            "identification of the physical blocked event",
            "a Poisson or projective-limit carrier",
        ],
        "SOURCE_SEMANTICS",
        "source boundary drift",
    )
    producer_cone = payload.get("producer_cone")
    require(
        producer_cone
        == {
            "consumed_source_fields": [
                "finite_source_schema",
                "finite_source_status",
                "presence_statement",
                "reserve_trace_branch.binding_theorem_present",
                "reserve_trace_branch.derived_from_collar_counts",
                "operational_clock_bound",
                "P_certified_enclosure",
            ],
            "cosmological_target_payload_used_in_countermodel": False,
            "desired_capacity_used_in_countermodel": False,
            "completion_selected_from_target": False,
            "full_parent_receipt_not_in_producer_cone": True,
            "projection_checked_against_parent_by_independent_replay": True,
        },
        "PROMOTION",
        "producer-cone boundary drift",
    )
    typed_class = payload.get("typed_countermodel_class")
    require(
        isinstance(typed_class, Mapping)
        and typed_class.get("connected_cut_source_object")
        == (
            "the free commutative cut-count monoid (Nat,+,0) carrying the "
            "same normalized local reading q"
        )
        and typed_class.get("global_capacity_object")
        == (
            "the positive real ray; the executable receipt serializes an "
            "exact rational witness subfamily"
        )
        and typed_class.get("refinement_relation")
        == (
            "finite partitions of one coarse cut count whose entries have "
            "the same total"
        )
        and "remains conditional" in str(typed_class.get("scope")),
        "SOURCE_SEMANTICS",
        "typed countermodel class drift",
    )
    blocked = payload.get("blocked_event_verdict")
    require(isinstance(blocked, Mapping), "PAYLOAD_TYPE", "blocked-event verdict missing")
    for key in (
        "one_class_selected",
        "six_class_total_selected",
        "no_capacity_action_selected",
    ):
        require(blocked.get(key) is False, "SELECTION", f"{key} was promoted")
    require(
        blocked.get("declared_completion_set")
        == [
            "no_capacity_action",
            "one_class_projector",
            "six_class_total_projector",
        ],
        "SELECTION",
        "blocked-event declared completion set drift",
    )
    require(
        blocked.get("same_source_datum")
        == "one normalized local class expectation q"
        and blocked.get("status")
        == "SOURCE_DOES_NOT_SELECT_AMONG_DECLARED_BLOCKED_EVENT_ACTIONS",
        "SELECTION",
        "blocked-event bounded verdict drift",
    )
    boundary = payload.get("comparison_boundary")
    require(isinstance(boundary, Mapping), "PAYLOAD_TYPE", "comparison boundary missing")
    require(
        boundary and all(value is False for value in boundary.values()),
        "PROMOTION",
        "comparison or capacity promotion enabled",
    )
    composition = payload.get("composition_boundary")
    require(
        isinstance(composition, Mapping)
        and composition.get("disconnected_cut_composition_classified") is True
        and composition.get("per_edge_multiplication_authorized") is False,
        "COMPOSITION",
        "composition boundary drift",
    )
    recorded_hash = payload.get("receipt_sha256")
    require(isinstance(recorded_hash, str), "HASH", "receipt hash missing")
    unhashed = dict(payload)
    unhashed.pop("receipt_sha256", None)
    require(
        recorded_hash == tagged_sha256(canonical_json_bytes(unhashed)),
        "HASH",
        "receipt hash drift",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.validate_only:
        payload = _load(args.output)
        validate(payload, verify_source=True)
        print("GLOBAL_CAPACITY_ATTACHMENT_VALID")
        return 0
    payload = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(canonical_json_bytes(payload))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
