#!/usr/bin/env python3
"""Exact reversible PUBLIC_CHECKPOINT_PACKET reference implementation.

The packet is a finite schema/control built on the twelve ports of an
icosahedral screen.  It is deliberately target-free: no measured Lambda,
electroweak value, operational scale, or desired capacity enters its producer.
It demonstrates the fast theorem branch

    every checkpoint generator injective  =>  M_0(q) = |X_reach(q)|.

It is not the missing source-derived physical packet or the finite-size slack
law.  Those remain the theorem-producing deliverables tracked by the repo.
"""
from __future__ import annotations

import json
from typing import Any, Mapping

from correctable_public_record_capacity import (
    _channel_rows,
    evaluate_terminal,
    public_global_sections,
    reachable_public_sections,
)


PORTS = ("north", "south", *(f"upper_{i}" for i in range(5)), *(f"lower_{i}" for i in range(5)))


def icosahedral_edges() -> list[tuple[str, str]]:
    """Return the 30 undirected edges of a standard twelve-vertex icosahedron."""
    edges: set[tuple[str, str]] = set()

    def add(left: str, right: str) -> None:
        edges.add(tuple(sorted((left, right))))

    for i in range(5):
        upper = f"upper_{i}"
        lower = f"lower_{i}"
        add("north", upper)
        add("south", lower)
        add(upper, f"upper_{(i + 1) % 5}")
        add(lower, f"lower_{(i + 1) % 5}")
        add(upper, lower)
        add(upper, f"lower_{(i - 1) % 5}")
    assert len(edges) == 30
    return sorted(edges)


def build_reference_packet(capacity_dimension: int = 4) -> dict[str, Any]:
    """Build a connected twelve-port packet with reversible checkpoint generators."""
    if capacity_dimension < 2:
        raise ValueError("reference packet requires at least two record labels")
    atoms = [f"record_{i}" for i in range(capacity_dimension)]
    observers = {port: atoms for port in PORTS}
    interfaces = [
        {
            "interface_id": f"{left}--{right}",
            "left_observer": left,
            "right_observer": right,
            "left_readout": {atom: atom for atom in atoms},
            "right_readout": {atom: atom for atom in atoms},
        }
        for left, right in icosahedral_edges()
    ]
    sections = public_global_sections(observers, interfaces)
    section_ids = ["|".join(f"{port}={section[port]}" for port in sorted(section)) for section in sections]
    label_to_section = {
        section["north"]: sid for section, sid in zip(sections, section_ids, strict=True)
    }

    def rows(permutation: Mapping[int, int]) -> dict[str, dict[str, float]]:
        return {
            label_to_section[f"record_{i}"]: {label_to_section[f"record_{permutation[i]}"]: 1.0}
            for i in range(capacity_dimension)
        }

    identity = {i: i for i in range(capacity_dimension)}
    rotate = {i: (i + 1) % capacity_dimension for i in range(capacity_dimension)}
    reflect = {i: (-i) % capacity_dimension for i in range(capacity_dimension)}
    all_ports = list(PORTS)
    return {
        "schema": "PUBLIC_CHECKPOINT_PACKET/v1-reversible",
        "status": "REFERENCE_CONTROL_NOT_PHYSICAL_RECEIPT",
        "regulator_id": "icosahedral-12-port-reference-r0",
        "carrier_type_id": "public-record-carrier",
        "terminal_id": "icosahedral-reversible-q0",
        "capacity_dimension": capacity_dimension,
        "observers": observers,
        "interfaces": interfaces,
        "reachability_witnesses": {
            sid: ["source_seed", f"endogenous_write_{i}", "public_seam_check"]
            for i, sid in enumerate(section_ids)
        },
        "publicness_policy": [all_ports],
        "continuation_manifest_complete": True,
        "global_checkpoint_kernels": [
            {
                "authorized_observers": all_ports,
                "continuation_id": "identity",
                "rows": rows(identity),
            },
            {
                "authorized_observers": all_ports,
                "continuation_id": "cyclic_relabel",
                "rows": rows(rotate),
            },
            {
                "authorized_observers": all_ports,
                "continuation_id": "orientation_reflection",
                "rows": rows(reflect),
            },
        ],
        "local_marginal_consistency_passed": True,
        "projection_supports": {sid: [i] for i, sid in enumerate(section_ids)},
        "refinement_maps": [
            {
                "from_regulator": "icosahedral-12-port-reference-r0",
                "to_regulator": "icosahedral-12-port-reference-r1",
                "record_embedding": {sid: sid for sid in section_ids},
                "no_new_confusability": True,
            }
        ],
        "self_read_predicate_injected": False,
        "supplied_capacity_metadata_read_by_producer": False,
        "lambda_used": False,
        "ew_bridge_used": False,
        "rho_used": False,
    }


def certify_reversible_packet(packet: Mapping[str, Any]) -> dict[str, Any]:
    """Verify injective deterministic generators and the fast capacity identity."""
    sections = public_global_sections(packet["observers"], packet["interfaces"])
    reachable = reachable_public_sections(sections, packet["reachability_witnesses"])
    generator_receipts: list[dict[str, Any]] = []
    for channel in packet["global_checkpoint_kernels"]:
        normalized = _channel_rows(channel, reachable)
        image: list[str] = []
        for source in reachable:
            positive = [(output, probability) for output, probability in normalized[source].items() if probability > 0]
            if len(positive) != 1 or positive[0][1] != 1.0:
                return {
                    "status": "NONDETERMINISTIC_CHECKPOINT_GENERATOR",
                    "continuation_id": channel.get("continuation_id"),
                    "source": source,
                }
            image.append(positive[0][0])
        if len(set(image)) != len(reachable):
            return {
                "status": "NONINJECTIVE_CHECKPOINT_GENERATOR",
                "continuation_id": channel.get("continuation_id"),
            }
        if set(image) != set(reachable):
            return {
                "status": "CHECKPOINT_NOT_CLOSED_ON_REACHABLE_RECORDS",
                "continuation_id": channel.get("continuation_id"),
            }
        generator_receipts.append(
            {
                "continuation_id": channel["continuation_id"],
                "deterministic": True,
                "injective": True,
                "permutation_of_reachable_records": True,
            }
        )

    evaluation = evaluate_terminal(packet)
    if evaluation["status"] != "PASS":
        return evaluation
    count = len(reachable)
    if evaluation["exact_zero_error_capacity"] != count:
        return {
            "status": "REVERSIBLE_CAPACITY_IDENTITY_FAILED",
            "reachable_count": count,
            "capacity": evaluation["exact_zero_error_capacity"],
        }
    return {
        "status": "PASS",
        "packet_status": packet.get("status"),
        "regulator_id": packet.get("regulator_id"),
        "port_count": len(packet["observers"]),
        "interface_count": len(packet["interfaces"]),
        "public_global_section_count": len(sections),
        "reachable_public_record_count": count,
        "checkpoint_generator_receipts": generator_receipts,
        "exact_zero_error_capacity": evaluation["exact_zero_error_capacity"],
        "capacity_dimension": evaluation["capacity_dimension"],
        "fast_branch_identity": "M_0(q)=|X_reach(q)|",
        "robust_terminal_saturation": evaluation["saturation_passed"],
        "rank_one_complete": evaluation["saturation_rank_one_complete"],
        "remaining_physical_obligation": "derive the source packet family and unique finite-size slack zero",
    }


def main() -> None:
    packet = build_reference_packet()
    print(json.dumps({"packet": packet, "certificate": certify_reversible_packet(packet)}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
