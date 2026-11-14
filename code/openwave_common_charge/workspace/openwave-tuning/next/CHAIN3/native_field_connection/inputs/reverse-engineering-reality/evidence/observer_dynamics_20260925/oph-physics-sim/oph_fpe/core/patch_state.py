from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Hashable


@dataclass
class PatchState:
    """Finite state owned by one observer patch."""

    hidden: int
    ports: dict[int, Hashable] = field(default_factory=dict)
    gauges: dict[int, Hashable] = field(default_factory=dict)
    scalar: float = 0.0
    phase: float = 0.0
    modular_time: float = 0.0
    modular_depth: float = 0.0
    repair_load: float = 0.0
    capacity: int = 0
    record: tuple[Any, ...] | None = None
    candidate_record: tuple[Any, ...] | None = None
    stable_count: int = 0
    commit_count: int = 0

    def observer_packet(self) -> tuple[Any, ...]:
        port_packet = tuple((neighbor, self.ports[neighbor]) for neighbor in sorted(self.ports))
        scalar_bin = round(self.scalar, 3)
        modular_bin = round(self.modular_depth, 3)
        return (self.hidden, scalar_bin, modular_bin, port_packet)

    def copy(self) -> "PatchState":
        return PatchState(
            hidden=self.hidden,
            ports=dict(self.ports),
            gauges=dict(self.gauges),
            scalar=self.scalar,
            phase=self.phase,
            modular_time=self.modular_time,
            modular_depth=self.modular_depth,
            repair_load=self.repair_load,
            capacity=self.capacity,
            record=self.record,
            candidate_record=self.candidate_record,
            stable_count=self.stable_count,
            commit_count=self.commit_count,
        )

    def to_jsonable(self, group_label) -> dict[str, Any]:
        return {
            "hidden": self.hidden,
            "scalar": self.scalar,
            "phase": self.phase,
            "modular_time": self.modular_time,
            "modular_depth": self.modular_depth,
            "repair_load": self.repair_load,
            "capacity": self.capacity,
            "ports": {str(k): group_label(v) for k, v in sorted(self.ports.items())},
            "gauges": {str(k): group_label(v) for k, v in sorted(self.gauges.items())},
            "record": _packet_to_jsonable(self.record, group_label),
            "stable_count": self.stable_count,
            "commit_count": self.commit_count,
        }


def _packet_to_jsonable(packet: tuple[Any, ...] | None, group_label) -> Any:
    if packet is None:
        return None
    hidden, scalar_bin, modular_bin, ports = packet
    return {
        "hidden": hidden,
        "scalar_bin": scalar_bin,
        "modular_bin": modular_bin,
        "ports": [(neighbor, group_label(value)) for neighbor, value in ports],
    }
