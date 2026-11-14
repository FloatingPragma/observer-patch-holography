"""Reference-free policy checks for particle inputs used by simulator lanes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable


GENERATIVE_DESTINATIONS = frozenset(
    {
        "source_law",
        "repair_kernel",
        "initial_state",
        "carrier_action",
        "species_assignment",
        "promotion_gate",
    }
)

NON_GENERATIVE_CLASSES = frozenset(
    {
        "measured_endpoint_comparison",
        "external_data_hadron_closure",
        "empirical_external_data_closure",
        "conditional_compare_only_envelope",
        "rejected_target_informed_candidate_status",
        "retrospective_rejection_evidence",
        "compare_only_scheme_analysis",
    }
)


@dataclass(frozen=True)
class ParticleInputRecord:
    name: str
    epistemic_class: str
    destination: str
    sha256: str | None = None
    used: bool = True


def particle_input_non_circularity_report(
    records: Iterable[ParticleInputRecord | dict[str, Any]],
) -> dict[str, Any]:
    normalized = [_record(row) for row in records]
    blockers: list[str] = []
    for row in normalized:
        if not row.used:
            continue
        if row.epistemic_class in NON_GENERATIVE_CLASSES and row.destination in GENERATIVE_DESTINATIONS:
            blockers.append(f"non_generative_input_entered_{row.destination}:{row.name}")
        if not _valid_hash(row.sha256):
            blockers.append(f"input_hash_missing_or_invalid:{row.name}")
    return {
        "schema": "oph_particle_input_non_circularity_v1",
        "records": [row.__dict__ for row in normalized],
        "reference_free_receipt": not blockers,
        "no_target_leak_receipt": not blockers,
        "blockers": blockers,
        "policy": {
            "generative_destinations": sorted(GENERATIVE_DESTINATIONS),
            "non_generative_classes": sorted(NON_GENERATIVE_CLASSES),
        },
    }


def _record(value: ParticleInputRecord | dict[str, Any]) -> ParticleInputRecord:
    if isinstance(value, ParticleInputRecord):
        return value
    if not isinstance(value, dict):
        raise TypeError("particle input record must be a ParticleInputRecord or mapping")
    return ParticleInputRecord(
        name=str(value.get("name") or value.get("key") or "unnamed"),
        epistemic_class=str(value.get("epistemic_class") or "undeclared"),
        destination=str(value.get("destination") or "undeclared"),
        sha256=value.get("sha256") if isinstance(value.get("sha256"), str) else None,
        used=value.get("used", True) is True,
    )


def _valid_hash(value: str | None) -> bool:
    if not isinstance(value, str):
        return False
    raw = value.removeprefix("sha256:")
    return len(raw) == 64 and all(char in "0123456789abcdefABCDEF" for char in raw)
