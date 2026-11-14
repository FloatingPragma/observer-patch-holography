from __future__ import annotations

from typing import Any

from oph_fpe.dynamics.repair import RepairEvent


def verify_local_law(
    repair_events: list[RepairEvent],
    record_events: list[Any],
    commit_cycles: int,
) -> list[dict[str, Any]]:
    receipts: list[dict[str, Any]] = []
    for event in repair_events:
        ok = event.delta_phi <= 1e-12 or event.reason == "metropolis_hot" or not event.accepted
        receipts.append(
            {
                "verifier": "local_repair_acceptance",
                "cycle": event.cycle,
                "node": event.node,
                "ok": ok,
                "delta_phi": event.delta_phi,
                "reason": event.reason,
            }
        )
    for event in record_events:
        receipts.append(
            {
                "verifier": "record_stability",
                "cycle": event.cycle,
                "node": event.node,
                "ok": event.stable_count >= commit_cycles,
                "stable_count": event.stable_count,
                "required": commit_cycles,
            }
        )
    return receipts
