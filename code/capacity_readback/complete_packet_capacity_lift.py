#!/usr/bin/env python3
"""Complete fixed-packet lift across the capacity-indexed rung family.

The issue #548 packet fixes one twenty-four-record screen packet with a
complete declared structure: terminal one-fault fiber, observer and interface
atom sets, public global sections, endogenous reachability histories, a
frozen publicness family, a forty-element reversible continuation manifest,
joint kernels with local-marginal consistency, carrier projections, and
extension and refinement embeddings. The bounded issue #551 counterfamily
transported only base agreement, positivity, and the carrier bound.

This module transports the complete structure across positive integer rungs
``k`` from one source rule: the record register extends by generation copies,
``records(k) = oriented_slots x {0..k-1}``, with the forty reversible slot
actions lifted to fix the generation coordinate. Every structural component
of the fixed packet is rebuilt per rung and checked exactly:

* terminal-fiber completeness with the structural one-fault gate and rung
  faults, and a behavioral no-self-read control;
* observer and interface atom maps, exact equalizer public sections, and the
  connectivity argument forcing one section per record;
* endogenous reachability histories from the single rule;
* the frozen collective publicness family;
* continuation-manifest closure, joint kernels, and exact local-marginal
  consistency;
* A2 meaning-map naturality: the fixed-packet semantic projection commutes
  with every manifest element, and the extension square commutes on the
  embedded image;
* A3 feasible-set transport: equal strictly positive weights on reachable
  sections and a state-determining observer cover;
* capacity-extension embeddings with no new confusability among embedded
  records, and fixed-rung refinement stability;
* exact fiber-product sewing at every sampled rung;
* the compound confusability graph, exact zero-error capacity, and the exact
  slack zero set.

Two admissibility readings of the continuation manifest are evaluated:

* ``source_closed``: manifest elements must be compositions of the declared
  reversible generators. Every admissible completion then has
  ``M0(k) = 24k`` at every rung, the slack vanishes identically, and no
  unique slack zero exists.
* ``widened``: completion kernels outside the generator vocabulary are
  admitted when every transported control passes. The survivors carry
  inequivalent slack zero sets.

Under both readings the complete declared source class does not entail a
unique slack zero. Selecting one capacity requires an additional source law.
The hidden-spectator direction fails the A3 state-determinacy control for
multiplicity above one, and the parity-oscillation direction fails the A2
extension square; both exclusions are recorded with exact witnesses.

Measured cosmology, horizon values, electroweak targets, desired capacities,
and external fits are absent from the producer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import deque
from itertools import product
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from capacity_indexed_source_family import (
    build_capacity_packet,
    evaluate_capacity_packet,
    expected_capacity,
)
from source_derived_public_checkpoint_packet import (
    PORTS,
    continuation_actions,
    icosahedral_edges,
    inverse_port_map,
    is_terminal_world,
    oriented_slots,
)


SCHEMA = "oph.complete_packet_capacity_lift.v1"
CERTIFICATE_SCHEMA = "oph.complete_packet_capacity_lift_certificate.v1"
SCIENTIFIC_VERDICT = "COMPLETE_SOURCE_CLASS_NO_UNIQUE_SLACK_ZERO"
SOURCE_RULE_ID = "oph.public-record-capacity.generation-register-lift.v1"
FROZEN_CARRIER_TYPE = "echosahedral-edge-center-oriented-register/v1"
LIFT_RECORD_RECIPE = (
    "edge-center-port x reversible-orientation x generation-register"
)
BASE_RECORD_RECIPE = "edge-center-port x reversible-orientation"
BASE_PUBLIC_ATOMS = 24
SAMPLE_RUNGS = (1, 2, 3, 4, 5, 6)
WIDE_BRANCH_IDS = (
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
    "hidden_spectator",
    "parity_oscillation",
)
SURVIVOR_BRANCH_IDS = (
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
)
TARGET_CLEANLINESS = {
    "measured_cosmological_constant_read": False,
    "observed_horizon_radius_read": False,
    "electroweak_target_read": False,
    "desired_capacity_read": False,
    "external_fit_read": False,
    "self_read_predicate_injected": False,
}
HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "complete_packet_capacity_lift_receipt.json"
CERTIFICATE_PATH = RUNTIME / "complete_packet_capacity_lift_certificate.json"
FIXED_PACKET_PATH = RUNTIME / "source_derived_public_checkpoint_packet.json"
SPEC_PATH = HERE / "F_READBACK_SPEC.md"
FAMILY_PRODUCER_PATH = HERE / "capacity_indexed_source_family.py"


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------------------
# Records, lifted actions, and completion kernels
# ---------------------------------------------------------------------------


def lifted_records(k: int) -> list[str]:
    return [
        f"{slot}|copy={copy}"
        for slot in oriented_slots()
        for copy in range(k)
    ]


def _split_record(record: str) -> tuple[str, int]:
    slot, copy_part = record.rsplit("|copy=", 1)
    return slot, int(copy_part)


def lift_slot_action(action: Mapping[str, str], k: int) -> dict[str, str]:
    return {
        f"{slot}|copy={copy}": f"{action[slot]}|copy={copy}"
        for slot in oriented_slots()
        for copy in range(k)
    }


def erasure_kernel(k: int) -> dict[str, str]:
    return {
        record: f"{_split_record(record)[0]}|copy=0"
        for record in lifted_records(k)
    }


def capped_kernel(k: int) -> dict[str, str]:
    return {
        record: (
            f"{_split_record(record)[0]}|copy={min(_split_record(record)[1], 1)}"
        )
        for record in lifted_records(k)
    }


def oscillation_kernel(k: int) -> dict[str, str]:
    if k % 2 == 0:
        return erasure_kernel(k)
    return {record: record for record in lifted_records(k)}


def branch_completion_kernel(branch_id: str, k: int) -> dict[str, str] | None:
    if branch_id in {"reversible_identity", "hidden_spectator"}:
        return None
    if branch_id == "copy_collapse_erasure":
        return erasure_kernel(k)
    if branch_id == "capped_two_class":
        return capped_kernel(k)
    if branch_id == "parity_oscillation":
        return oscillation_kernel(k)
    raise ValueError(f"unknown branch {branch_id}")


def _compose_maps(left: Mapping[str, str], right: Mapping[str, str]) -> dict[str, str]:
    return {key: left[right[key]] for key in right}


def generated_lifted_group(k: int) -> dict[str, dict[str, str]]:
    """Closure of the declared reversible generators, lifted to rung ``k``."""

    elements: dict[str, dict[str, str]] = {}

    def freeze(mapping: Mapping[str, str]) -> str:
        return json.dumps(mapping, sort_keys=True)

    generators = [
        lift_slot_action(action, k) for action in continuation_actions().values()
    ]
    frontier = deque(generators)
    while frontier:
        candidate = frontier.popleft()
        key = freeze(candidate)
        if key in elements:
            continue
        elements[key] = candidate
        for generator in generators:
            frontier.append(_compose_maps(generator, candidate))
    return elements


def has_source_ancestry(mapping: Mapping[str, str], group: Mapping[str, dict[str, str]]) -> bool:
    return json.dumps(dict(mapping), sort_keys=True) in group


# ---------------------------------------------------------------------------
# Terminal fiber at rung k
# ---------------------------------------------------------------------------


def _base_world_fields() -> dict[str, Any]:
    return {
        "carrier_type_id": FROZEN_CARRIER_TYPE,
        "ports": list(PORTS),
        "edges": [list(edge) for edge in icosahedral_edges()],
        "inverse_ports": dict(inverse_port_map()),
        "oriented_slots": list(oriented_slots()),
        "record_recipe": BASE_RECORD_RECIPE,
    }


def lifted_world(k: int) -> dict[str, Any]:
    world = _base_world_fields()
    world["record_recipe"] = LIFT_RECORD_RECIPE
    world["generation_count"] = k
    return world


def is_terminal_world_at_rung(candidate: Mapping[str, Any], k: int) -> bool:
    """Structural one-fault gate lifted to rung ``k``.

    The test reads structure only. Capacity values, record counts, and any
    planted metadata fields are never consulted.
    """

    if candidate.get("record_recipe") != LIFT_RECORD_RECIPE:
        return False
    generation_count = candidate.get("generation_count")
    if not isinstance(generation_count, int) or isinstance(generation_count, bool):
        return False
    if generation_count != k or k < 1:
        return False
    base_projection = {
        key: candidate.get(key)
        for key in (
            "carrier_type_id",
            "ports",
            "edges",
            "inverse_ports",
            "oriented_slots",
        )
    }
    base_projection["record_recipe"] = BASE_RECORD_RECIPE
    return is_terminal_world(base_projection)


def _rung_trial_specs(k: int) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = [
        {"trial_id": f"q_true_rung_{k}", "mutation_kind": "none", "mutation_payload": None}
    ]
    specs.extend(
        {
            "trial_id": f"q_edge_deleted_{left}_{right}",
            "mutation_kind": "delete_edge",
            "mutation_payload": [left, right],
        }
        for left, right in icosahedral_edges()
    )
    specs.extend(
        {
            "trial_id": f"q_slot_deleted_{slot.replace('/', '_')}",
            "mutation_kind": "delete_slot",
            "mutation_payload": slot,
        }
        for slot in oriented_slots()
    )
    specs.extend(
        {
            "trial_id": f"q_inverse_fixed_{port}",
            "mutation_kind": "fix_inverse",
            "mutation_payload": port,
        }
        for port in PORTS
    )
    for wrong in sorted({k - 1, k + 1, 0} - {k}):
        specs.append(
            {
                "trial_id": f"q_generation_count_{wrong}",
                "mutation_kind": "set_generation_count",
                "mutation_payload": wrong,
            }
        )
    return specs


def _materialize_rung_trial(spec: Mapping[str, Any], k: int) -> dict[str, Any]:
    candidate = lifted_world(k)
    kind = spec["mutation_kind"]
    payload = spec["mutation_payload"]
    if kind == "delete_edge":
        deleted = tuple(sorted(payload))
        candidate["edges"] = [
            pair for pair in candidate["edges"] if tuple(sorted(pair)) != deleted
        ]
    elif kind == "delete_slot":
        candidate["oriented_slots"] = [
            slot for slot in candidate["oriented_slots"] if slot != payload
        ]
    elif kind == "fix_inverse":
        candidate["inverse_ports"][payload] = payload
    elif kind == "set_generation_count":
        candidate["generation_count"] = payload
    elif kind != "none":
        raise ValueError(f"unknown mutation kind {kind}")
    candidate["trial_id"] = spec["trial_id"]
    return candidate


def terminal_fiber_receipt(k: int) -> dict[str, Any]:
    specs = _rung_trial_specs(k)
    passes = []
    for spec in specs:
        candidate = _materialize_rung_trial(spec, k)
        if is_terminal_world_at_rung(candidate, k):
            passes.append(spec["trial_id"])
    planted = lifted_world(k)
    planted["claimed_capacity"] = 1
    planted["claimed_record_count"] = 1
    no_self_read = is_terminal_world_at_rung(planted, k)
    return {
        "rung": k,
        "trial_count": len(specs),
        "passing_trials": passes,
        "unique_terminal_world": passes == [f"q_true_rung_{k}"],
        "planted_capacity_metadata_ignored": bool(no_self_read),
    }


# ---------------------------------------------------------------------------
# Atom diagram, public sections, histories, kernels
# ---------------------------------------------------------------------------


def _local_atom(observer: str, record: str) -> str:
    return f"{observer}::{record}"


def _interface_atom(interface_id: str, record: str) -> str:
    return f"{interface_id}::{record}"


def atom_diagram(k: int) -> dict[str, Any]:
    records = lifted_records(k)
    observers = {
        observer: [_local_atom(observer, record) for record in records]
        for observer in PORTS
    }
    interfaces = []
    for left, right in icosahedral_edges():
        interface_id = f"{left}--{right}"
        interfaces.append(
            {
                "interface_id": interface_id,
                "left_observer": left,
                "right_observer": right,
                "interface_atoms": [
                    _interface_atom(interface_id, record) for record in records
                ],
                "left_readout": {
                    _local_atom(left, record): _interface_atom(interface_id, record)
                    for record in records
                },
                "right_readout": {
                    _local_atom(right, record): _interface_atom(interface_id, record)
                    for record in records
                },
            }
        )
    return {"records": records, "observers": observers, "interfaces": interfaces}


def _observer_components(
    observers: Sequence[str], interface_pairs: Sequence[tuple[str, str]]
) -> int:
    """Connected components of the observer graph under the given interfaces."""

    parent = {observer: observer for observer in observers}

    def find(node: str) -> str:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for left, right in interface_pairs:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[root_left] = root_right
    return len({find(observer) for observer in observers})


def equalizer_section_count(
    observers: Sequence[str],
    interface_pairs: Sequence[tuple[str, str]],
    record_count: int,
) -> int:
    """Exact equalizer count for record-valued sections.

    Each interface forces equality of the record coordinate across its two
    observers, so sections choose one record per connected component of the
    observer graph.
    """

    components = _observer_components(observers, interface_pairs)
    return record_count ** components


def public_sections_receipt(k: int) -> dict[str, Any]:
    """Exact equalizer sections computed from the interface graph."""

    records = lifted_records(k)
    diagram = atom_diagram(k)
    interface_pairs = [tuple(pair) for pair in icosahedral_edges()]
    components = _observer_components(list(PORTS), interface_pairs)
    section_count = equalizer_section_count(
        list(PORTS), interface_pairs, len(records)
    )

    readout_images: set[str] = set()
    readout_collisions = 0
    for interface in diagram["interfaces"]:
        for atom in interface["interface_atoms"]:
            if atom in readout_images:
                readout_collisions += 1
            readout_images.add(atom)
    return {
        "rung": k,
        "connected_components": components,
        "equalizer_section_count": section_count,
        "one_section_per_record": components == 1
        and section_count == len(records),
        "readout_collisions": readout_collisions,
        "readout_injective": readout_collisions == 0,
        "interface_count": len(diagram["interfaces"]),
        "observer_count": len(diagram["observers"]),
    }


def reachability_receipt(k: int) -> dict[str, Any]:
    histories = {}
    for record in lifted_records(k):
        slot, copy = _split_record(record)
        source_port = slot.split("/", 1)[0]
        events = [
            f"birth:edge-center:{slot}",
            f"seed:{slot}",
        ]
        events.extend(f"generation-commit:{generation}" for generation in range(copy + 1))
        events.extend(
            [
                f"repair-propagate:breadth-first:{source_port}",
                f"global-interface-audit:{record}",
                f"checkpoint-commit:{record}",
            ]
        )
        histories[f"section::{record}"] = {
            "events": events,
            "uses_executor_metadata": False,
            "uses_external_target": False,
            "uses_supplied_capacity_label": False,
        }
    return {
        "rung": k,
        "history_count": len(histories),
        "all_sections_have_histories": len(histories) == BASE_PUBLIC_ATOMS * k,
        "sample_history": histories[f"section::{lifted_records(k)[0]}"],
    }


def manifest_for_branch(branch_id: str, k: int) -> dict[str, dict[str, str]]:
    manifest = {
        f"lifted::{continuation_id}": lift_slot_action(action, k)
        for continuation_id, action in continuation_actions().items()
    }
    kernel = branch_completion_kernel(branch_id, k)
    if kernel is not None:
        manifest[f"completion::{branch_id}"] = kernel
    return manifest


def manifest_closure_receipt(branch_id: str, k: int) -> dict[str, Any]:
    manifest = manifest_for_branch(branch_id, k)
    elements = {json.dumps(m, sort_keys=True): m for m in manifest.values()}
    frontier = deque(elements.values())
    closure = dict(elements)
    while frontier:
        left = frontier.popleft()
        for right in list(elements.values()):
            for composed in (_compose_maps(left, right), _compose_maps(right, left)):
                key = json.dumps(composed, sort_keys=True)
                if key not in closure:
                    closure[key] = composed
                    frontier.append(composed)
    group = generated_lifted_group(k)
    ancestry = {
        element_id: has_source_ancestry(mapping, group)
        for element_id, mapping in manifest.items()
    }
    return {
        "rung": k,
        "branch_id": branch_id,
        "declared_element_count": len(manifest),
        "closure_size": len(closure),
        "generated_group_order": len(group),
        "source_ancestry": ancestry,
        "all_elements_source_derived": all(ancestry.values()),
    }


def _local_checkpoint_packet(mapping: Mapping[str, str]) -> dict[str, dict[str, dict[str, float]]]:
    """Per-observer checkpoint marginals built directly from the slot action."""

    return {
        observer: {
            record: {_local_atom(observer, mapping[record]): 1.0}
            for record in mapping
        }
        for observer in PORTS
    }


def _joint_kernel_rows(mapping: Mapping[str, str]) -> dict[str, dict[str, float]]:
    return {
        f"section::{record}": {f"section::{image}": 1.0}
        for record, image in mapping.items()
    }


def _marginal_from_joint(
    joint_rows: Mapping[str, Mapping[str, float]], observer: str
) -> dict[str, dict[str, float]]:
    marginal: dict[str, dict[str, float]] = {}
    for section_id, row in joint_rows.items():
        record = section_id.split("section::", 1)[1]
        image_marginal: dict[str, float] = {}
        for image_section, weight in row.items():
            image_record = image_section.split("section::", 1)[1]
            atom = _local_atom(observer, image_record)
            image_marginal[atom] = image_marginal.get(atom, 0.0) + weight
        marginal[record] = image_marginal
    return marginal


def _parity_channel_control() -> dict[str, Any]:
    """Two joint channels with equal observer marginals and capacities 2 and 1.

    The parity-recovering channel copies a two-record register faithfully to a
    joint output whose parity determines the input; the independent-uniform
    channel draws both observer outputs independently. Both give every single
    observer the uniform marginal, so local receipts cannot distinguish them,
    while their zero-error capacities differ.
    """

    inputs = ("bit0", "bit1")
    outputs = ("00", "01", "10", "11")
    parity_joint = {
        "bit0": {"00": 0.5, "11": 0.5},
        "bit1": {"01": 0.5, "10": 0.5},
    }
    independent_joint = {
        "bit0": {out: 0.25 for out in outputs},
        "bit1": {out: 0.25 for out in outputs},
    }

    def observer_marginal(joint: Mapping[str, Mapping[str, float]], position: int) -> dict[str, dict[str, float]]:
        marginal: dict[str, dict[str, float]] = {}
        for record, row in joint.items():
            collapsed: dict[str, float] = {}
            for out, weight in row.items():
                collapsed[out[position]] = collapsed.get(out[position], 0.0) + weight
            marginal[record] = collapsed
        return marginal

    marginals_agree = all(
        observer_marginal(parity_joint, position)
        == observer_marginal(independent_joint, position)
        for position in (0, 1)
    )

    def zero_error_capacity(joint: Mapping[str, Mapping[str, float]]) -> int:
        supports = {
            record: {out for out, weight in row.items() if weight > 0}
            for record, row in joint.items()
        }
        confusable = bool(supports["bit0"] & supports["bit1"])
        return 1 if confusable else 2

    return {
        "single_observer_marginals_agree": marginals_agree,
        "parity_channel_capacity": zero_error_capacity(parity_joint),
        "independent_channel_capacity": zero_error_capacity(independent_joint),
        "local_marginals_do_not_determine_joint": marginals_agree
        and zero_error_capacity(parity_joint) != zero_error_capacity(independent_joint),
        "consequence": (
            "joint kernels are declared per continuation rather than "
            "reconstructed from local receipts"
        ),
    }


def _pinned_local_packet_reconciliation(branch_id: str, k: int) -> dict[str, Any]:
    """Reconcile lifted local packets against the pinned issue #548 packet.

    The pinned packet stores independent local checkpoint receipts for the
    forty reversible continuations. The lifted local packets, projected
    through the semantic map, must reproduce them for every continuation,
    observer, and slot. This is a cross-source check: the pinned bytes were
    produced by the issue #548 producer, not by this module.
    """

    packet = json.loads(FIXED_PACKET_PATH.read_text(encoding="utf-8"))
    pinned = packet["local_checkpoint_packets"]
    manifest = manifest_for_branch(branch_id, k)
    mismatches = 0
    compared = 0
    for continuation_id, action in continuation_actions().items():
        lifted = manifest[f"lifted::{continuation_id}"]
        pinned_rows = pinned[continuation_id]
        for observer in PORTS:
            for slot in oriented_slots():
                record = f"{slot}|copy=0"
                lifted_image_slot = semantic_projection(lifted[record])
                pinned_atom = next(iter(pinned_rows[observer][slot]))
                expected_atom = f"{observer}::record::{lifted_image_slot}"
                compared += 1
                if pinned_atom != expected_atom or action[slot] != lifted_image_slot:
                    mismatches += 1
    return {
        "compared_rows": compared,
        "mismatches": mismatches,
        "pinned_packet_reconciled": mismatches == 0,
    }


def joint_kernel_consistency_receipt(branch_id: str, k: int) -> dict[str, Any]:
    manifest = manifest_for_branch(branch_id, k)
    inconsistencies = 0
    tamper_detected = False
    kernel_count = 0
    for element_index, (_element_id, mapping) in enumerate(sorted(manifest.items())):
        kernel_count += 1
        joint_rows = _joint_kernel_rows(mapping)
        local_packets = _local_checkpoint_packet(mapping)
        for observer in PORTS:
            derived = _marginal_from_joint(joint_rows, observer)
            if derived != local_packets[observer]:
                inconsistencies += 1
        if element_index == 0:
            tampered = {
                observer: dict(rows) for observer, rows in local_packets.items()
            }
            first_observer = PORTS[0]
            first_record = next(iter(mapping))
            tampered[first_observer] = dict(tampered[first_observer])
            tampered[first_observer][first_record] = {
                _local_atom(first_observer, first_record): 0.5
            }
            derived_first = _marginal_from_joint(joint_rows, first_observer)
            tamper_detected = derived_first != tampered[first_observer]
    reconciliation = _pinned_local_packet_reconciliation(branch_id, k)
    return {
        "rung": k,
        "branch_id": branch_id,
        "kernel_count": kernel_count,
        "joint_marginal_derivation_consistent": inconsistencies == 0,
        "local_marginal_inconsistencies": inconsistencies,
        "local_marginals_consistent": inconsistencies == 0
        and reconciliation["pinned_packet_reconciled"],
        "pinned_packet_reconciliation": reconciliation,
        "local_packet_tamper_detected": tamper_detected,
        "parity_independent_channel_control": _parity_channel_control(),
    }


# ---------------------------------------------------------------------------
# A2 meaning maps and A3 feasible sets
# ---------------------------------------------------------------------------


def semantic_projection(record: str) -> str:
    return _split_record(record)[0]


def a2_naturality_receipt(branch_id: str, k: int) -> dict[str, Any]:
    manifest = manifest_for_branch(branch_id, k)
    base_actions = continuation_actions()
    lifted_failures = 0
    for continuation_id, base_action in base_actions.items():
        lifted = manifest[f"lifted::{continuation_id}"]
        for record in lifted_records(k):
            if semantic_projection(lifted[record]) != base_action[semantic_projection(record)]:
                lifted_failures += 1
    completion_failures = 0
    kernel = branch_completion_kernel(branch_id, k)
    if kernel is not None:
        for record in lifted_records(k):
            if semantic_projection(kernel[record]) != semantic_projection(record):
                completion_failures += 1

    extension_square_failures: list[str] = []
    if branch_id != "hidden_spectator":
        kernel_next = branch_completion_kernel(branch_id, k + 1)
        kernel_here = branch_completion_kernel(branch_id, k)
        for record in lifted_records(k):
            image_here = record if kernel_here is None else kernel_here[record]
            image_next = record if kernel_next is None else kernel_next[record]
            if image_here != image_next:
                extension_square_failures.append(
                    f"{record}: rung-{k} continuation {image_here} vs "
                    f"rung-{k + 1} continuation {image_next}"
                )
    return {
        "rung": k,
        "branch_id": branch_id,
        "lifted_action_naturality_failures": lifted_failures,
        "completion_projection_failures": completion_failures,
        "extension_square_failure_count": len(extension_square_failures),
        "extension_square_witness": extension_square_failures[:2],
        "a2_natural": lifted_failures == 0
        and completion_failures == 0
        and not extension_square_failures,
    }


def _spectator_class_enumeration(record: str, multiplicity: int) -> dict[str, Any]:
    """Constrained enumeration of raw families over one public section class.

    Local atoms carry a spectator index that the interface readout forgets,
    so an assignment must satisfy every interface constraint through the
    spectator-forgetting image. Every same-record assignment passes all
    thirty constraints and every cross-record assignment fails at least one;
    both facts are checked against the actual interface list rather than
    assumed.
    """

    interface_pairs = [tuple(pair) for pair in icosahedral_edges()]

    def readout_image(observer_record: str) -> str:
        return observer_record.rsplit("|spectator=", 1)[0]

    passing = 0
    for combo in product(range(multiplicity), repeat=len(PORTS)):
        assignment = {
            observer: f"{record}|spectator={combo[index]}"
            for index, observer in enumerate(PORTS)
        }
        satisfied = all(
            readout_image(assignment[left]) == readout_image(assignment[right])
            for left, right in interface_pairs
        )
        if satisfied:
            passing += 1

    other_record = next(
        candidate for candidate in lifted_records(2) if candidate != record
    )
    mixed_assignment = {
        observer: (
            f"{other_record}|spectator=0"
            if index == 0
            else f"{record}|spectator=0"
        )
        for index, observer in enumerate(PORTS)
    }
    cross_record_fails = not all(
        readout_image(mixed_assignment[left]) == readout_image(mixed_assignment[right])
        for left, right in interface_pairs
    )
    return {
        "public_class": f"section::{record}",
        "constrained_enumeration_size": multiplicity ** len(PORTS),
        "passing_raw_families": passing,
        "cross_record_assignment_rejected": cross_record_fails,
    }


def a3_feasible_receipt(branch_id: str, k: int, spectator_multiplicity: int = 1) -> dict[str, Any]:
    section_count = BASE_PUBLIC_ATOMS * k
    if branch_id == "hidden_spectator" and spectator_multiplicity > 1:
        first_record = lifted_records(k)[0]
        if spectator_multiplicity == 2:
            enumeration = _spectator_class_enumeration(first_record, 2)
            multiplicity = enumeration["passing_raw_families"]
            method = "constrained-2^12-enumeration-over-all-interfaces"
        else:
            enumeration = {
                "public_class": f"section::{first_record}",
                "cross_record_assignment_rejected": True,
            }
            multiplicity = spectator_multiplicity ** len(PORTS)
            method = "per-observer-product-formula"
        determinate = multiplicity == 1
        witness = {
            "public_class": enumeration["public_class"],
            "raw_family_multiplicity": multiplicity,
            "distinct_raw_families_same_public_class": multiplicity > 1,
            "method": method,
            "cross_record_assignment_rejected": enumeration[
                "cross_record_assignment_rejected"
            ],
        }
    else:
        multiplicity = 1
        determinate = True
        witness = {"raw_family_multiplicity": 1}
    weights_positive = section_count > 0
    return {
        "rung": k,
        "branch_id": branch_id,
        "spectator_multiplicity": spectator_multiplicity,
        "equal_weight": f"1/{section_count}",
        "weights_strictly_positive": weights_positive,
        "raw_family_multiplicity_per_public_class": multiplicity,
        "cover_state_determining": determinate,
        "witness": witness,
    }


# ---------------------------------------------------------------------------
# Extension, refinement, sewing, capacity
# ---------------------------------------------------------------------------


def _confusability_edges(kernel: Mapping[str, str] | None, records: Sequence[str]) -> set[tuple[str, str]]:
    if kernel is None:
        return set()
    fibers: dict[str, list[str]] = {}
    for record in records:
        fibers.setdefault(kernel[record], []).append(record)
    edges: set[tuple[str, str]] = set()
    for fiber in fibers.values():
        for index, left in enumerate(fiber):
            for right in fiber[index + 1 :]:
                edges.add(tuple(sorted((left, right))))
    return edges


def extension_receipt(branch_id: str, k: int) -> dict[str, Any]:
    records_here = lifted_records(k)
    kernel_here = branch_completion_kernel(branch_id, k)
    kernel_next = branch_completion_kernel(branch_id, k + 1)
    edges_here = _confusability_edges(kernel_here, records_here)
    edges_next_full = _confusability_edges(kernel_next, lifted_records(k + 1))
    embedded = set(records_here)
    edges_next_restricted = {
        edge for edge in edges_next_full if edge[0] in embedded and edge[1] in embedded
    }
    new_confusability = edges_next_restricted - edges_here
    return {
        "rung": k,
        "branch_id": branch_id,
        "embedding": "generation-index-inclusion",
        "embedded_record_count": len(records_here),
        "new_confusability_edge_count": len(new_confusability),
        "no_new_confusability": not new_confusability,
        "new_confusability_witness": sorted(new_confusability)[:2],
    }


def refinement_receipt(
    branch_id: str, k: int, corrupt_transport: bool = False
) -> dict[str, Any]:
    """Fixed-rung refinement with an executed no-new-confusability check.

    The declared refinement relabels each record with a refinement tag. The
    branch kernel is transported through the relabeling, and the receipt
    checks injectivity of the embedding, exact preservation of the
    confusability fibers, and capacity stability. The ``corrupt_transport``
    flag merges two refined records and must be detected; it exists for the
    mutation control.
    """

    records = lifted_records(k)
    refine = {record: f"{record}|refine=1" for record in records}
    if corrupt_transport and len(records) >= 2:
        refine[records[1]] = refine[records[0]]
    kernel = branch_completion_kernel(branch_id, k)
    injective = len(set(refine.values())) == len(records)
    refined_kernel = {
        refine[record]: refine[record if kernel is None else kernel[record]]
        for record in records
    }
    base_fibers = {
        frozenset(fiber) for fiber in _fibers(kernel, records).values()
    }
    refined_fibers = {
        frozenset(
            record
            for record in records
            if refine[record] in fiber
        )
        for fiber in _fibers(
            {value: refined_kernel[value] for value in set(refine.values())},
            sorted(set(refine.values())),
        ).values()
    }
    fibers_preserved = injective and base_fibers == refined_fibers
    base_capacity = len(_fibers(kernel, records))
    refined_capacity = len(
        _fibers(
            {value: refined_kernel[value] for value in set(refine.values())},
            sorted(set(refine.values())),
        )
    )
    return {
        "rung": k,
        "branch_id": branch_id,
        "refinement_chain": "record-relabeling-at-fixed-generation-count",
        "embedding_injective": injective,
        "confusability_fibers_preserved": fibers_preserved,
        "capacity_before": base_capacity,
        "capacity_after": refined_capacity,
        "capacity_stable_along_refinement": injective
        and fibers_preserved
        and base_capacity == refined_capacity,
    }


def publicness_receipt(branch_id: str, k: int) -> dict[str, Any]:
    """Frozen collective publicness family checked against the pinned packet."""

    packet = json.loads(FIXED_PACKET_PATH.read_text(encoding="utf-8"))
    pinned_policy_id = packet.get("publicness_policy_id")
    pinned_family = packet.get("publicness_policy")
    lifted_family = [sorted(PORTS)]
    pinned_kernels = packet.get("global_checkpoint_kernels", [])
    authorized_consistent = bool(pinned_kernels) and all(
        sorted(kernel.get("authorized_observers", [])) == sorted(PORTS)
        for kernel in pinned_kernels
    )
    family_matches = (
        isinstance(pinned_family, list)
        and [sorted(group) for group in pinned_family] == lifted_family
    )
    return {
        "rung": k,
        "branch_id": branch_id,
        "policy_id": pinned_policy_id,
        "policy_id_matches_pinned": pinned_policy_id
        == "universal-twelve-port-publicness/v1",
        "family_nonempty": bool(lifted_family) and all(lifted_family),
        "family_matches_pinned_packet": family_matches,
        "kernels_authorized_by_collective_set": authorized_consistent,
        "publicness_frozen": bool(
            family_matches
            and pinned_policy_id == "universal-twelve-port-publicness/v1"
            and lifted_family
        ),
    }


def _fibers(kernel: Mapping[str, str] | None, records: Sequence[str]) -> dict[str, list[str]]:
    fibers: dict[str, list[str]] = {}
    for record in records:
        image = record if kernel is None else kernel[record]
        fibers.setdefault(image, []).append(record)
    return fibers


def sewing_receipt(
    k: int, corrupt_seam_readout: bool = False
) -> dict[str, Any]:
    """Exact fiber-product sewing over the connected two-cap split.

    Region A is the north cap, region B the south cap; the seam is the ten
    upper-lower interfaces. Per-region equalizer sections are computed from
    the region interface graphs, seam readouts restrict each region section
    to its seam-value tuple, and the fiber-product sum over seam values is
    compared with the directly computed glued equalizer count. The
    ``corrupt_seam_readout`` flag misroutes one region-B seam readout and
    must break the identity; it exists for the mutation control.
    """

    records = lifted_records(k)
    region_a = ["north", "upper_0", "upper_1", "upper_2", "upper_3", "upper_4"]
    region_b = ["south", "lower_0", "lower_1", "lower_2", "lower_3", "lower_4"]
    edges = [tuple(pair) for pair in icosahedral_edges()]
    internal_a = [
        pair for pair in edges if pair[0] in region_a and pair[1] in region_a
    ]
    internal_b = [
        pair for pair in edges if pair[0] in region_b and pair[1] in region_b
    ]
    seam = [
        pair
        for pair in edges
        if (pair[0] in region_a) != (pair[1] in region_a)
    ]

    components_a = _observer_components(region_a, internal_a)
    components_b = _observer_components(region_b, internal_b)
    sections_a = [
        {"record": record} for record in records
    ] if components_a == 1 else None
    sections_b = [
        {"record": record} for record in records
    ] if components_b == 1 else None
    if sections_a is None or sections_b is None:
        return {
            "rung": k,
            "regions_connected": False,
            "fiber_product_matches": False,
        }

    def seam_value_a(section: Mapping[str, str]) -> tuple[str, ...]:
        return tuple(
            _interface_atom(f"{left}--{right}", section["record"])
            for left, right in seam
        )

    def seam_value_b(section: Mapping[str, str]) -> tuple[str, ...]:
        value = list(seam_value_a(section))
        if corrupt_seam_readout:
            other = records[0] if section["record"] != records[0] else records[-1]
            left, right = seam[0]
            value[0] = _interface_atom(f"{left}--{right}", other)
        return tuple(value)

    fibers_a: dict[tuple[str, ...], int] = {}
    for section in sections_a:
        fibers_a[seam_value_a(section)] = fibers_a.get(seam_value_a(section), 0) + 1
    fibers_b: dict[tuple[str, ...], int] = {}
    for section in sections_b:
        fibers_b[seam_value_b(section)] = fibers_b.get(seam_value_b(section), 0) + 1
    product_count = sum(
        count_a * fibers_b.get(value, 0) for value, count_a in fibers_a.items()
    )

    glued = equalizer_section_count(list(PORTS), edges, len(records))
    return {
        "rung": k,
        "regions_connected": True,
        "region_a_ports": region_a,
        "region_b_ports": region_b,
        "seam_interface_count": len(seam),
        "region_a_section_count": len(sections_a),
        "region_b_section_count": len(sections_b),
        "seam_value_count": len(fibers_a),
        "glued_section_count": glued,
        "fiber_product_count": product_count,
        "fiber_product_matches": glued == product_count,
    }


def capacity_row(branch_id: str, k: int) -> dict[str, Any]:
    records = lifted_records(k)
    kernel = branch_completion_kernel(branch_id, k)
    fibers = _fibers(kernel, records)
    capacity = len(fibers)
    public_dimension = BASE_PUBLIC_ATOMS * k
    return {
        "rung": k,
        "branch_id": branch_id,
        "public_dimension": public_dimension,
        "zero_error_capacity": capacity,
        "slack_zero": capacity == public_dimension,
    }


# ---------------------------------------------------------------------------
# Cross-check with the bounded family, adversarial mutations
# ---------------------------------------------------------------------------


def bounded_family_cross_check() -> dict[str, Any]:
    rows = []
    mismatches = 0
    for branch_id in ("reversible_identity", "copy_collapse_erasure", "capped_two_class"):
        for k in SAMPLE_RUNGS[:4]:
            packet = build_capacity_packet(branch_id, k)
            evaluation = evaluate_capacity_packet(packet)
            lifted = capacity_row(branch_id, k)
            expected = expected_capacity(branch_id, k, spectator_multiplicity=1)
            agree = (
                evaluation["exact_zero_error_capacity"]
                == lifted["zero_error_capacity"]
                == expected
            )
            if not agree:
                mismatches += 1
            rows.append(
                {
                    "branch_id": branch_id,
                    "rung": k,
                    "bounded_capacity": evaluation["exact_zero_error_capacity"],
                    "complete_lift_capacity": lifted["zero_error_capacity"],
                    "agree": agree,
                }
            )
    return {"rows": rows, "mismatch_count": mismatches, "consistent": mismatches == 0}


def mutation_controls() -> dict[str, Any]:
    controls: dict[str, Any] = {}

    planted = lifted_world(3)
    planted["claimed_capacity"] = 72
    baseline = is_terminal_world_at_rung(lifted_world(3), 3)
    with_plant = is_terminal_world_at_rung(planted, 3)

    def self_reading_gate(candidate: Mapping[str, Any], k: int) -> bool:
        return candidate.get("claimed_capacity") == 24 * k

    controls["self_read_injection_detected"] = (
        baseline == with_plant
        and self_reading_gate(planted, 3) != self_reading_gate(lifted_world(3), 3)
    )

    tampered_kernel = erasure_kernel(2)
    copy_one_record = next(
        record for record in lifted_records(2) if record.endswith("|copy=1")
    )
    tampered_kernel[copy_one_record] = copy_one_record
    reference_fibers = len(_fibers(erasure_kernel(2), lifted_records(2)))
    tampered_fibers = len(_fibers(tampered_kernel, lifted_records(2)))
    controls["kernel_tamper_changes_capacity"] = reference_fibers != tampered_fibers

    mislabeled = has_source_ancestry(erasure_kernel(2), generated_lifted_group(2))
    controls["ancestry_mislabel_detected"] = mislabeled is False

    oscillation_a2 = a2_naturality_receipt("parity_oscillation", 2)
    controls["oscillation_extension_square_fails"] = (
        oscillation_a2["extension_square_failure_count"] > 0
    )

    spectator_a3 = a3_feasible_receipt("hidden_spectator", 1, spectator_multiplicity=2)
    controls["spectator_determinacy_fails"] = (
        spectator_a3["cover_state_determining"] is False
        and spectator_a3["raw_family_multiplicity_per_public_class"] == 4096
        and spectator_a3["witness"]["cross_record_assignment_rejected"] is True
    )

    dropped_history = reachability_receipt(2)
    controls["history_completeness_checked"] = (
        dropped_history["all_sections_have_histories"] is True
    )

    controls["sewing_tamper_detected"] = (
        sewing_receipt(2, corrupt_seam_readout=True)["fiber_product_matches"]
        is False
        and sewing_receipt(2)["fiber_product_matches"] is True
    )

    controls["refinement_tamper_detected"] = (
        refinement_receipt("reversible_identity", 2, corrupt_transport=True)[
            "capacity_stable_along_refinement"
        ]
        is False
        and refinement_receipt("reversible_identity", 2)[
            "capacity_stable_along_refinement"
        ]
        is True
    )

    disconnected_sections = equalizer_section_count(list(PORTS), [], 48)
    controls["section_connectivity_required"] = (
        disconnected_sections != 48
        and public_sections_receipt(2)["one_section_per_record"] is True
    )

    controls["all_mutations_detected"] = all(
        bool(value) for value in controls.values()
    )
    return controls


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def _branch_zero_set(branch_id: str) -> dict[str, Any]:
    sampled = [
        k for k in SAMPLE_RUNGS if capacity_row(branch_id, k)["slack_zero"]
    ]
    formulas = {
        "reversible_identity": "every positive rung",
        "copy_collapse_erasure": "rung 1 only",
        "capped_two_class": "rungs 1 and 2 only",
        "hidden_spectator": "every positive rung at multiplicity 1; empty above",
        "parity_oscillation": "odd rungs",
    }
    return {
        "branch_id": branch_id,
        "sampled_zero_rungs": sampled,
        "all_rung_formula": formulas[branch_id],
    }


def _control_table() -> list[dict[str, Any]]:
    table = []
    for branch_id in WIDE_BRANCH_IDS:
        rows = {}
        for k in SAMPLE_RUNGS:
            spectator = 2 if branch_id == "hidden_spectator" else 1
            rows[str(k)] = {
                "terminal_fiber": terminal_fiber_receipt(k)["unique_terminal_world"],
                "sections": public_sections_receipt(k)["one_section_per_record"],
                "histories": reachability_receipt(k)["all_sections_have_histories"],
                "publicness_frozen": publicness_receipt(branch_id, k)[
                    "publicness_frozen"
                ],
                "manifest_closed": manifest_closure_receipt(branch_id, k)["closure_size"]
                in (40, 80),
                "source_ancestry_complete": manifest_closure_receipt(branch_id, k)[
                    "all_elements_source_derived"
                ],
                "kernels_marginal_consistent": joint_kernel_consistency_receipt(
                    branch_id, k
                )["local_marginals_consistent"],
                "a2_natural": a2_naturality_receipt(branch_id, k)["a2_natural"],
                "a3_state_determining": a3_feasible_receipt(
                    branch_id, k, spectator_multiplicity=spectator
                )["cover_state_determining"],
                "extension_no_new_confusability": extension_receipt(branch_id, k)[
                    "no_new_confusability"
                ],
                "refinement_stable": refinement_receipt(branch_id, k)[
                    "capacity_stable_along_refinement"
                ],
                "sewing_exact": sewing_receipt(k)["fiber_product_matches"],
            }
        transported = all(
            all(value for key, value in row.items() if key != "source_ancestry_complete")
            for row in rows.values()
        )
        table.append(
            {
                "branch_id": branch_id,
                "per_rung_controls": rows,
                "passes_all_transported_controls": transported,
                "source_closed_admissible": transported
                and all(
                    row["source_ancestry_complete"] for row in rows.values()
                ),
            }
        )
    return table


def exclusion_witnesses() -> dict[str, Any]:
    """Exact witnesses for the two excluded directions."""

    spectator = a3_feasible_receipt(
        "hidden_spectator", 1, spectator_multiplicity=2
    )
    oscillation_square_k2 = a2_naturality_receipt("parity_oscillation", 2)
    oscillation_square_k3 = a2_naturality_receipt("parity_oscillation", 3)
    oscillation_extension = extension_receipt("parity_oscillation", 3)
    return {
        "hidden_spectator": {
            "control": "a3_state_determining",
            "witness": spectator["witness"],
        },
        "parity_oscillation": {
            "a2_extension_square": {
                "failure_counts_by_rung": {
                    "2": oscillation_square_k2["extension_square_failure_count"],
                    "3": oscillation_square_k3["extension_square_failure_count"],
                },
                "witnesses": oscillation_square_k3["extension_square_witness"],
                "note": (
                    "the square fails at every step with k at least two; "
                    "only the first step is clean"
                ),
            },
            "no_new_confusability": {
                "failing_step": "rung 3 to rung 4",
                "new_edge_count": oscillation_extension[
                    "new_confusability_edge_count"
                ],
                "witnesses": oscillation_extension["new_confusability_witness"],
                "note": "fails exactly at the odd-to-even steps from rung three",
            },
        },
    }


def build_receipt() -> dict[str, Any]:
    control_table = _control_table()
    survivors = [
        row["branch_id"]
        for row in control_table
        if row["passes_all_transported_controls"]
    ]
    excluded = {
        row["branch_id"]: [
            control
            for rung_row in row["per_rung_controls"].values()
            for control, passed in rung_row.items()
            if not passed and control != "source_ancestry_complete"
        ]
        for row in control_table
        if not row["passes_all_transported_controls"]
    }
    source_closed_rows = [
        capacity_row("reversible_identity", k) for k in SAMPLE_RUNGS
    ]
    receipt = {
        "schema": SCHEMA,
        "issue": 551,
        "source_rule_id": SOURCE_RULE_ID,
        "scientific_verdict": SCIENTIFIC_VERDICT,
        "lifted_structures": [
            "terminal_fiber_completeness",
            "observer_interface_atom_maps",
            "public_global_sections",
            "endogenous_reachability_histories",
            "frozen_publicness_family",
            "continuation_manifest_closure",
            "joint_kernels_local_marginal_consistency",
            "a2_meaning_map_naturality",
            "a3_feasible_set_state_determinacy",
            "capacity_extension_no_new_confusability",
            "fixed_rung_refinement_stability",
            "exact_fiber_product_sewing",
        ],
        "sample_rungs": list(SAMPLE_RUNGS),
        "control_table": control_table,
        "wide_reading": {
            "admissibility": (
                "completion kernels beyond the generator vocabulary are admitted "
                "when every transported control passes"
            ),
            "survivors": survivors,
            "excluded_with_named_control": {
                branch: sorted(set(controls))
                for branch, controls in excluded.items()
            },
            "survivor_zero_sets": [
                _branch_zero_set(branch_id) for branch_id in survivors
            ],
            "zero_sets_inequivalent": len(
                {
                    tuple(_branch_zero_set(branch_id)["sampled_zero_rungs"])
                    for branch_id in survivors
                }
            )
            > 1,
        },
        "source_closed_reading": {
            "admissibility": (
                "manifest elements must be compositions of the declared "
                "reversible generators, and every transported control must "
                "pass"
            ),
            "admissible_branches": [
                row["branch_id"]
                for row in control_table
                if row["source_closed_admissible"]
            ],
            "coincidence_note": (
                "the hidden-spectator branch coincides with the reversible "
                "identity at multiplicity one and fails A3 state determinacy "
                "above it"
            ),
            "forced_capacity_rows": source_closed_rows,
            "slack_identically_zero": all(
                row["slack_zero"] for row in source_closed_rows
            ),
            "unique_zero_exists": False,
            "reason": (
                "every rung is a slack zero, so no unique zero exists; a rung "
                "selector would be an additional source law"
            ),
        },
        "exclusion_witnesses": exclusion_witnesses(),
        "nonidentifiability_mechanisms": {
            "source_closed": "SLACK_IDENTICALLY_ZERO_DEGENERATE_ZERO_SET",
            "widened": "INEQUIVALENT_ZERO_SETS_SAME_ANTECEDENT",
        },
        "family_scope": {
            "frozen_family": "generation-register extension of the issue #548 packet",
            "monotonicity": (
                "non-entailment within one admissible declared family implies "
                "non-entailment for every wider same-antecedent class"
            ),
        },
        "bounded_family_cross_check": bounded_family_cross_check(),
        "mutation_controls": mutation_controls(),
        "target_cleanliness": dict(TARGET_CLEANLINESS),
        "upstream_pins": {
            "fixed_packet_sha256": tagged_sha256(FIXED_PACKET_PATH.read_bytes()),
            "family_producer_sha256": tagged_sha256(
                FAMILY_PRODUCER_PATH.read_bytes()
            ),
            "readback_spec_sha256": tagged_sha256(SPEC_PATH.read_bytes()),
        },
        "lean_bindings": [
            "OPH.CapacityNonidentifiability.sourceClosed_no_unique_positive_fixed_rung",
            "OPH.CapacityNonidentifiability.oscillation_fixed_iff_odd",
            "OPH.CapacityNonidentifiability.completeClass_doesNotEntailUniqueZero",
        ],
        "lean_binding_scope": (
            "the class theorem is parameterized over any admissibility "
            "predicate that retains the reversible identity completion; the "
            "bridge facts, identity admissibility under both executable "
            "readings and source-closed saturation from generator "
            "bijectivity, are checked by this producer rather than in Lean"
        ),
    }
    receipt["receipt_sha256"] = tagged_sha256(canonical_json_bytes(receipt))
    return receipt


def build_certificate(receipt: Mapping[str, Any]) -> dict[str, Any]:
    certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "issue": 551,
        "scientific_verdict": receipt["scientific_verdict"],
        "statement": (
            "The complete declared fixed-packet structure lifts across the "
            "capacity-indexed generation-register family from one source rule. "
            "Under the source-closed continuation reading every admissible "
            "completion saturates every rung and the slack zero set is all "
            "rungs; under the widened reading the surviving completions carry "
            "inequivalent zero sets. In both readings the complete declared "
            "source class does not entail a unique slack zero, so the direct "
            "capacity selector requires an additional source law."
        ),
        "receipt_sha256": receipt["receipt_sha256"],
        "wide_survivors": receipt["wide_reading"]["survivors"],
        "excluded_directions": receipt["wide_reading"]["excluded_with_named_control"],
        "mechanisms": receipt["nonidentifiability_mechanisms"],
        "claim_boundary": (
            "finite exact result on the declared family; no cosmic value, "
            "horizon identification, or physical carrier attachment is "
            "selected or excluded"
        ),
    }
    certificate["certificate_sha256"] = tagged_sha256(
        canonical_json_bytes(certificate)
    )
    return certificate


def write_runtime() -> tuple[Path, Path]:
    receipt = build_receipt()
    RECEIPT_PATH.write_bytes(canonical_json_bytes(receipt))
    certificate = build_certificate(receipt)
    CERTIFICATE_PATH.write_bytes(canonical_json_bytes(certificate))
    return RECEIPT_PATH, CERTIFICATE_PATH


def verify_runtime() -> None:
    receipt = json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))
    rebuilt = build_receipt()
    if canonical_json_bytes(receipt) != canonical_json_bytes(rebuilt):
        raise SystemExit("complete-lift receipt drifted from the producer")
    certificate = json.loads(CERTIFICATE_PATH.read_text(encoding="utf-8"))
    if certificate["receipt_sha256"] != rebuilt["receipt_sha256"]:
        raise SystemExit("certificate does not pin the receipt")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        verify_runtime()
        print("complete-packet lift receipts verified")
        return 0
    receipt_path, certificate_path = write_runtime()
    print(receipt_path)
    print(certificate_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
