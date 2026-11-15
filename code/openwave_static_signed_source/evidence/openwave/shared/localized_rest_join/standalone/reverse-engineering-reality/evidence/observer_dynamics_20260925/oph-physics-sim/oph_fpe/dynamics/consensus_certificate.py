from __future__ import annotations

import math
from typing import Any, Iterable

from oph_fpe.claims import FINITE_CONSENSUS_THEOREM_RECEIPT, RECOVERED_CORE, with_claim_metadata


def finite_consensus_theorem_certificate(
    trace: Iterable[dict[str, Any]],
    evidence: dict[str, Any] | None = None,
    *,
    strict_tol: float = 1.0e-12,
) -> dict[str, Any]:
    """Fail-closed C0b consensus theorem certificate.

    Ordinary simulator traces only show that the current search settled.  C0b
    needs theorem-phase strict descent plus replay/confluence evidence.  Missing
    evidence therefore fails the theorem receipt instead of promoting settling.
    """

    rows = [dict(row) for row in trace if isinstance(row, dict)]
    evidence = dict(evidence or {})
    try:
        tolerance = float(strict_tol)
    except (TypeError, ValueError, OverflowError):
        tolerance = math.nan
    tolerance_valid = bool(math.isfinite(tolerance) and tolerance >= 0.0)
    theorem_rows = [row for row in rows if str(row.get("phase", "")) == "theorem"]
    accepted_theorem_rows = [row for row in theorem_rows if row.get("accepted") is True]
    strict_descent_violations: list[dict[str, Any]] = []
    phi_increase_violations: list[dict[str, Any]] = []
    for index, row in enumerate(accepted_theorem_rows):
        touched_delta = _finite_float(row.get("delta_touched_phi", row.get("delta_phi")))
        global_delta = _finite_float(row.get("delta_global_phi", row.get("delta_phi")))
        if touched_delta is None or not tolerance_valid or touched_delta >= -tolerance:
            strict_descent_violations.append(
                _event_violation(row, index=index, key="delta_touched_phi", value=touched_delta)
            )
        if global_delta is None or not tolerance_valid or global_delta > tolerance:
            phi_increase_violations.append(
                _event_violation(row, index=index, key="delta_global_phi", value=global_delta)
            )

    required_replay_fields = {
        "disjoint_commutation_violation_count": 0,
        "local_diamond_violation_count": 0,
        "repair_completeness_violation_count": 0,
        "unique_terminal_quotient_hash_count": 1,
    }
    # Trace-derived fields are never overridden by caller summaries.  A
    # truncated/imported summary is useful diagnostics, but it is not a C0b
    # computed receipt.  Production array replays use the separate in-memory
    # verifier in ``bw_array`` which computes every count from the port arrays.
    theorem_event_count = len(theorem_rows)
    accepted_theorem_move_count = len(accepted_theorem_rows)
    strict_descent_violation_count = len(strict_descent_violations)
    phi_increase_violation_count = len(phi_increase_violations)
    declared_mismatches: list[str] = []
    for name, computed in (
        ("theorem_phase_event_count", theorem_event_count),
        ("accepted_theorem_move_count", accepted_theorem_move_count),
        ("strict_descent_violation_count", strict_descent_violation_count),
        ("accepted_phi_increase_violation_count", phi_increase_violation_count),
    ):
        if name in evidence and _strict_nonnegative_int(evidence.get(name)) != computed:
            declared_mismatches.append(name)
    missing_evidence: list[str] = []
    invalid_evidence: list[str] = []
    if not tolerance_valid:
        invalid_evidence.append("strict_descent_tolerance")
    if theorem_event_count <= 0:
        missing_evidence.append("theorem_phase_repair_events")
    for name in required_replay_fields:
        if name not in evidence:
            missing_evidence.append(name)

    replay_values: dict[str, int | None] = {
        name: _strict_nonnegative_int(evidence.get(name)) for name in required_replay_fields
    }
    disjoint_violations = replay_values["disjoint_commutation_violation_count"]
    diamond_violations = replay_values["local_diamond_violation_count"]
    completeness_violations = replay_values["repair_completeness_violation_count"]
    terminal_hash_count = replay_values["unique_terminal_quotient_hash_count"]
    schedule_replay_count = _strict_nonnegative_int(evidence.get("schedule_replay_count"))
    requested_schedule_replays = _strict_positive_int(evidence.get("requested_schedule_replays", 16))
    for name, value in (*replay_values.items(), ("schedule_replay_count", schedule_replay_count)):
        if name in evidence and value is None:
            invalid_evidence.append(name)
    if requested_schedule_replays is None:
        invalid_evidence.append("requested_schedule_replays")
    if "schedule_replay_count" not in evidence:
        missing_evidence.append("schedule_replay_count")
    if (
        schedule_replay_count is not None
        and requested_schedule_replays is not None
        and schedule_replay_count < requested_schedule_replays
    ):
        missing_evidence.append("sufficient_schedule_replays")

    imported_summary_checks_pass = bool(
        theorem_event_count > 0
        and not missing_evidence
        and not invalid_evidence
        and not declared_mismatches
        and strict_descent_violation_count == 0
        and phi_increase_violation_count == 0
        and disjoint_violations == required_replay_fields["disjoint_commutation_violation_count"]
        and diamond_violations == required_replay_fields["local_diamond_violation_count"]
        and completeness_violations == required_replay_fields["repair_completeness_violation_count"]
        and terminal_hash_count == required_replay_fields["unique_terminal_quotient_hash_count"]
        and schedule_replay_count is not None
        and requested_schedule_replays is not None
        and schedule_replay_count >= requested_schedule_replays
    )
    # Count-only replay evidence is declarative and therefore cannot promote a
    # strong computed theorem receipt.  The fully computed array replay path
    # emits that receipt directly after replaying schedules and AB/BA checks.
    passed = False
    missing_evidence.append("computed_replay_artifact")
    report = {
        "mode": "finite_consensus_theorem_certificate_v1",
        FINITE_CONSENSUS_THEOREM_RECEIPT: passed,
        "finite_consensus_theorem_receipt": passed,
        "receipt": passed,
        "theorem_phase_event_count": theorem_event_count,
        "accepted_theorem_move_count": accepted_theorem_move_count,
        "strict_descent_violation_count": strict_descent_violation_count,
        "accepted_phi_increase_violation_count": phi_increase_violation_count,
        "disjoint_commutation_violation_count": disjoint_violations,
        "local_diamond_violation_count": diamond_violations,
        "repair_completeness_violation_count": completeness_violations,
        "unique_terminal_quotient_hash_count": terminal_hash_count,
        "schedule_replay_count": schedule_replay_count,
        "requested_schedule_replays": requested_schedule_replays,
        "missing_evidence": sorted(set(missing_evidence)),
        "invalid_evidence": sorted(set(invalid_evidence)),
        "declared_summary_mismatches": sorted(set(declared_mismatches)),
        "imported_summary_checks_pass": imported_summary_checks_pass,
        "computed_replay_artifact_present": False,
        "strict_descent_violations": strict_descent_violations[:16],
        "accepted_phi_increase_violations": phi_increase_violations[:16],
        "claim_boundary": (
            "C0b finite consensus theorem receipt. This is stricter than final_phi == 0: it requires "
            "theorem-phase strict touched-overlap descent, disjoint/local replay checks, repair "
            "completeness, and schedule-confluent terminal quotient hashes. Count-only imported summaries "
            "remain diagnostic and cannot emit the strong computed receipt."
        ),
    }
    return with_claim_metadata(
        report,
        claim_level=RECOVERED_CORE,
        receipt=FINITE_CONSENSUS_THEOREM_RECEIPT,
        physical_claim=False,
        observable_id="theorem_phase_repair_replay_and_terminal_hashes",
        fit_objective="strict_descent_confluence_and_repair_completeness",
    )


def _event_violation(
    row: dict[str, Any], *, index: int, key: str, value: float | None
) -> dict[str, Any]:
    return {
        "index": index,
        "cycle": row.get("cycle"),
        "node": row.get("node"),
        key: value,
        "finite_delta": value is not None,
        "reason": row.get("reason"),
    }


def _finite_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return parsed if math.isfinite(parsed) else None


def _strict_nonnegative_int(value: Any) -> int | None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        return None
    return int(value)


def _strict_positive_int(value: Any) -> int | None:
    parsed = _strict_nonnegative_int(value)
    return parsed if parsed is not None and parsed > 0 else None
