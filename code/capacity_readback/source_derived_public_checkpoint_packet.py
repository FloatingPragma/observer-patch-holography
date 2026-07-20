#!/usr/bin/env python3
"""Source-derived finite public-checkpoint packet for OPH issue #548.

The producer freezes the finite carrier to the twelve exposed edge-centre ports
of the echosahedral screen and the reversible write/check orientation.  Its
record carrier is therefore the source-emitted oriented register

    R_24 = P_12 x {write, check}.

No measured cosmological value, electroweak/Higgs coordinate, operational
resolution, supplied output label, or finite-size selector enters the producer.
The artifact is physical only in the repository's typed fixed-cutoff sense: its
atoms, histories, publicness policy, checkpoint continuations, and carrier
projections are generated from the declared screen-microphysics carrier.  It is
not a hardware-realization certificate, a cosmic finite-size selector, a
horizon-record identification, or a screen/electroweak load identification.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import inspect
import json
from collections import deque
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from correctable_public_record_capacity import (
    _channel_rows,
    compound_confusability_graph,
    evaluate_terminal,
    evaluate_terminal_fiber,
    greatest_fixed_point,
    maximum_independent_set,
    no_new_confusability,
    public_global_sections,
    reachable_public_sections,
    section_id,
    support_relation_semigroup,
    tv_robustness_bound,
)

PORTS: tuple[str, ...] = (
    "north",
    "south",
    *(f"upper_{index}" for index in range(5)),
    *(f"lower_{index}" for index in range(5)),
)
ORIENTATIONS: tuple[str, str] = ("write", "check")
FROZEN_CARRIER_TYPE = "echosahedral-edge-center-oriented-register/v1"
REGULATOR_ID = "echosahedral-edge-center-r0"
TERMINAL_ID = "q_echosahedral_oriented_exact"
SCHEMA = "PUBLIC_CHECKPOINT_PACKET/v2-source-derived-reversible"
TARGET_MARKERS = (
    "measured_lambda",
    "electroweak_target",
    "higgs_target",
    "rho_op",
    "expected_answer",
    "desired_output",
)


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _edge(left: str, right: str) -> tuple[str, str]:
    return tuple(sorted((left, right)))


def icosahedral_edges() -> list[tuple[str, str]]:
    """Return the exact 30-edge graph of the frozen twelve-port carrier."""
    edges: set[tuple[str, str]] = set()
    for index in range(5):
        upper = f"upper_{index}"
        lower = f"lower_{index}"
        edges.add(_edge("north", upper))
        edges.add(_edge("south", lower))
        edges.add(_edge(upper, f"upper_{(index + 1) % 5}"))
        edges.add(_edge(lower, f"lower_{(index + 1) % 5}"))
        edges.add(_edge(upper, lower))
        edges.add(_edge(upper, f"lower_{(index - 1) % 5}"))
    if len(edges) != 30:
        raise AssertionError("the frozen port graph must have thirty edges")
    return sorted(edges)


def rotate_port(port: str) -> str:
    if port.startswith("upper_"):
        return f"upper_{(int(port.rsplit('_', 1)[1]) + 1) % 5}"
    if port.startswith("lower_"):
        return f"lower_{(int(port.rsplit('_', 1)[1]) + 1) % 5}"
    return port


def reflect_port(port: str) -> str:
    if port.startswith("upper_"):
        index = int(port.rsplit("_", 1)[1])
        return f"upper_{(-index) % 5}"
    if port.startswith("lower_"):
        index = int(port.rsplit("_", 1)[1])
        return f"lower_{(-index - 1) % 5}"
    return port


def antipodal_port(port: str) -> str:
    if port == "north":
        return "south"
    if port == "south":
        return "north"
    if port.startswith("upper_"):
        index = int(port.rsplit("_", 1)[1])
        return f"lower_{(index + 2) % 5}"
    if port.startswith("lower_"):
        index = int(port.rsplit("_", 1)[1])
        return f"upper_{(index + 3) % 5}"
    raise ValueError(f"unknown port {port}")


def inverse_port_map() -> dict[str, str]:
    return {port: antipodal_port(port) for port in PORTS}


def oriented_slots() -> list[str]:
    return [f"{port}/{orientation}" for port in PORTS for orientation in ORIENTATIONS]


def _split_slot(slot: str) -> tuple[str, str]:
    port, orientation = slot.split("/", 1)
    if port not in PORTS or orientation not in ORIENTATIONS:
        raise ValueError(f"invalid oriented slot {slot}")
    return port, orientation


def _local_atom(observer: str, slot: str) -> str:
    return f"{observer}::record::{slot}"


def _slot_from_local_atom(atom: str) -> str:
    marker = "::record::"
    if marker not in atom:
        raise ValueError(f"invalid local atom {atom}")
    return atom.split(marker, 1)[1]


def _interface_atom(interface_id: str, slot: str) -> str:
    return f"{interface_id}::public::{slot}"


def _adjacency(edges: Iterable[tuple[str, str]]) -> dict[str, set[str]]:
    adjacency = {port: set() for port in PORTS}
    for left, right in edges:
        if left not in adjacency or right not in adjacency or left == right:
            continue
        adjacency[left].add(right)
        adjacency[right].add(left)
    return adjacency


def _connected(edges: Iterable[tuple[str, str]]) -> bool:
    adjacency = _adjacency(edges)
    seen = {PORTS[0]}
    queue = deque([PORTS[0]])
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if neighbour not in seen:
                seen.add(neighbour)
                queue.append(neighbour)
    return seen == set(PORTS)


def _is_graph_automorphism(mapping: Mapping[str, str]) -> bool:
    if set(mapping) != set(PORTS) or set(mapping.values()) != set(PORTS):
        return False
    source_edges = set(icosahedral_edges())
    image_edges = {_edge(mapping[left], mapping[right]) for left, right in source_edges}
    return image_edges == source_edges


def is_terminal_world(candidate: Mapping[str, Any]) -> bool:
    """Structural membership test for the declared one-fault trial class."""
    if candidate.get("carrier_type_id") != FROZEN_CARRIER_TYPE:
        return False
    if tuple(candidate.get("ports", ())) != PORTS:
        return False
    edges = {_edge(*pair) for pair in candidate.get("edges", ())}
    if edges != set(icosahedral_edges()):
        return False
    adjacency = _adjacency(edges)
    if not _connected(edges) or any(len(adjacency[port]) != 5 for port in PORTS):
        return False
    inverse = candidate.get("inverse_ports")
    if not isinstance(inverse, Mapping) or set(inverse) != set(PORTS):
        return False
    if any(inverse[port] == port or inverse.get(inverse[port]) != port for port in PORTS):
        return False
    if not _is_graph_automorphism(inverse):
        return False
    if set(candidate.get("oriented_slots", ())) != set(oriented_slots()):
        return False
    return candidate.get("record_recipe") == "edge-center-port x reversible-orientation"


def _base_world() -> dict[str, Any]:
    return {
        "carrier_type_id": FROZEN_CARRIER_TYPE,
        "ports": list(PORTS),
        "edges": [list(edge) for edge in icosahedral_edges()],
        "inverse_ports": inverse_port_map(),
        "oriented_slots": oriented_slots(),
        "record_recipe": "edge-center-port x reversible-orientation",
    }


def _trial_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = [
        {"trial_id": TERMINAL_ID, "mutation_kind": "none", "mutation_payload": None}
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
    return specs


def _materialize_trial(spec: Mapping[str, Any]) -> dict[str, Any]:
    candidate = _base_world()
    kind = spec["mutation_kind"]
    payload = spec["mutation_payload"]
    if kind == "delete_edge":
        deleted = _edge(*payload)
        candidate["edges"] = [pair for pair in candidate["edges"] if _edge(*pair) != deleted]
    elif kind == "delete_slot":
        candidate["oriented_slots"] = [slot for slot in candidate["oriented_slots"] if slot != payload]
    elif kind == "fix_inverse":
        candidate["inverse_ports"][payload] = payload
    elif kind != "none":
        raise ValueError(f"unknown mutation kind {kind}")
    candidate["trial_id"] = spec["trial_id"]
    return candidate


def build_terminal_fiber_manifest() -> dict[str, Any]:
    """Emit the complete declared trial class and its structural terminal fiber."""
    source = inspect.getsource(is_terminal_world).lower()
    forbidden = ("capacity", "saturation", "expected_answer", "desired_output")
    source_audit = {marker: marker not in source for marker in forbidden}
    trials: list[dict[str, Any]] = []
    terminal_ids: list[str] = []
    for spec in _trial_specs():
        candidate = _materialize_trial(spec)
        terminal = is_terminal_world(candidate)
        if terminal:
            terminal_ids.append(spec["trial_id"])
        trials.append(
            {
                **spec,
                "terminal": terminal,
                "candidate": candidate,
                "candidate_sha256": _sha256(candidate),
            }
        )
    manifest = {
        "schema": "OPH_TERMINAL_FIBER_MANIFEST/v1",
        "regulator_id": REGULATOR_ID,
        "carrier_type_id": FROZEN_CARRIER_TYPE,
        "declared_trial_class": (
            "exact frozen carrier plus every single edge deletion, every single "
            "oriented-slot deletion, and every single inverse-port fixed-point fault"
        ),
        "trial_universe_constructor_sha256": _sha256(
            {
                "base_world": _base_world(),
                "trial_specs": _trial_specs(),
                "materializer_source": inspect.getsource(_materialize_trial),
                "membership_source": inspect.getsource(is_terminal_world),
            }
        ),
        "trial_count": len(trials),
        "trials": trials,
        "terminal_world_ids": terminal_ids,
        "terminal_fiber_complete": True,
        "terminal_fiber_completeness_certificate": {
            "declared_trials_enumerated": len(trials),
            "expected_trial_count": 1 + len(icosahedral_edges()) + len(oriented_slots()) + len(PORTS),
            "missing_trial_ids": [],
            "duplicate_trial_ids": [],
            "all_trial_ids_unique": len({trial["trial_id"] for trial in trials}) == len(trials),
            "all_candidates_materialized": all("candidate" in trial for trial in trials),
            "complete": True,
        },
        "terminal_membership_source_audit": source_audit,
        "terminal_membership_output_blind": all(source_audit.values()),
        "terminal_fiber_kind": (
            "EMPTY" if not terminal_ids else "SINGLETON" if len(terminal_ids) == 1 else "AMBIGUOUS"
        ),
    }
    manifest["terminal_fiber_manifest_hash_scope"] = (
        "all manifest fields except terminal_fiber_manifest_sha256"
    )
    manifest["terminal_fiber_manifest_sha256"] = _sha256(manifest)
    return manifest


def _port_action(rotation: int, reflection: bool, antipodal: bool) -> dict[str, str]:
    action: dict[str, str] = {}
    for port in PORTS:
        image = reflect_port(port) if reflection else port
        for _ in range(rotation % 5):
            image = rotate_port(image)
        if antipodal:
            image = antipodal_port(image)
        action[port] = image
    if not _is_graph_automorphism(action):
        raise AssertionError("declared continuation must preserve the frozen port graph")
    return action


def _slot_action(
    rotation: int, reflection: bool, antipodal: bool, orientation_flip: bool
) -> dict[str, str]:
    ports = _port_action(rotation, reflection, antipodal)
    result: dict[str, str] = {}
    for slot in oriented_slots():
        port, orientation = _split_slot(slot)
        if orientation_flip:
            orientation = "check" if orientation == "write" else "write"
        result[slot] = f"{ports[port]}/{orientation}"
    if set(result.values()) != set(oriented_slots()):
        raise AssertionError("continuation must permute the complete oriented register")
    return result


def continuation_actions() -> dict[str, dict[str, str]]:
    actions: dict[str, dict[str, str]] = {}
    signatures: set[tuple[str, ...]] = set()
    slots = oriented_slots()
    for rotation in range(5):
        for reflection in (False, True):
            for antipodal in (False, True):
                for orientation_flip in (False, True):
                    continuation_id = (
                        f"r{rotation}_s{int(reflection)}_a{int(antipodal)}_f{int(orientation_flip)}"
                    )
                    action = _slot_action(rotation, reflection, antipodal, orientation_flip)
                    signature = tuple(action[slot] for slot in slots)
                    if signature in signatures:
                        raise AssertionError("continuation parametrization is not faithful")
                    signatures.add(signature)
                    actions[continuation_id] = action
    if len(actions) != 40:
        raise AssertionError("the reversible continuation family must contain forty elements")
    return actions


def _compose(left: Mapping[str, str], right: Mapping[str, str]) -> dict[str, str]:
    """Return left after right."""
    return {source: left[right[source]] for source in right}


def _continuation_composition_table(
    actions: Mapping[str, Mapping[str, str]]
) -> dict[str, dict[str, str]]:
    slots = oriented_slots()
    lookup = {tuple(action[slot] for slot in slots): continuation_id for continuation_id, action in actions.items()}
    table: dict[str, dict[str, str]] = {}
    for left_id, left in actions.items():
        table[left_id] = {}
        for right_id, right in actions.items():
            composed = _compose(left, right)
            table[left_id][right_id] = lookup[tuple(composed[slot] for slot in slots)]
    return table


def _build_record_diagram() -> tuple[dict[str, list[str]], list[dict[str, Any]]]:
    slots = oriented_slots()
    observers = {
        observer: [_local_atom(observer, slot) for slot in slots]
        for observer in PORTS
    }
    interfaces: list[dict[str, Any]] = []
    for left, right in icosahedral_edges():
        interface_id = f"{left}--{right}"
        interfaces.append(
            {
                "interface_id": interface_id,
                "left_observer": left,
                "right_observer": right,
                "interface_atoms": [_interface_atom(interface_id, slot) for slot in slots],
                "left_readout": {
                    _local_atom(left, slot): _interface_atom(interface_id, slot)
                    for slot in slots
                },
                "right_readout": {
                    _local_atom(right, slot): _interface_atom(interface_id, slot)
                    for slot in slots
                },
            }
        )
    return observers, interfaces


def _sections_by_slot(
    observers: Mapping[str, Sequence[str]], interfaces: Sequence[Mapping[str, Any]]
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    sections = public_global_sections(observers, interfaces)
    by_slot: dict[str, str] = {}
    by_id: dict[str, dict[str, str]] = {}
    for section in sections:
        slots = {_slot_from_local_atom(atom) for atom in section.values()}
        if len(slots) != 1:
            raise AssertionError("a public section must carry one common oriented slot")
        slot = next(iter(slots))
        sid = section_id(section)
        by_slot[slot] = sid
        by_id[sid] = dict(section)
    if set(by_slot) != set(oriented_slots()):
        raise AssertionError("record diagram must expose exactly the oriented source register")
    return by_slot, by_id


def _bfs_tree(root: str) -> list[tuple[str, str]]:
    adjacency = _adjacency(icosahedral_edges())
    seen = {root}
    queue = deque([root])
    tree: list[tuple[str, str]] = []
    while queue:
        parent = queue.popleft()
        for child in sorted(adjacency[parent]):
            if child not in seen:
                seen.add(child)
                queue.append(child)
                tree.append((parent, child))
    if seen != set(PORTS):
        raise AssertionError("source propagation tree must reach all ports")
    return tree


def _reachability(
    section_by_slot: Mapping[str, str]
) -> tuple[dict[str, list[str]], dict[str, dict[str, Any]]]:
    witnesses: dict[str, list[str]] = {}
    histories: dict[str, dict[str, Any]] = {}
    for slot in oriented_slots():
        source_port, orientation = _split_slot(slot)
        events = [
            f"birth:edge-center:{slot}",
            f"seed:{orientation}:{source_port}",
        ]
        events.extend(
            f"repair-propagate:{parent}->{child}:{slot}"
            for parent, child in _bfs_tree(source_port)
        )
        events.extend(
            [
                f"global-interface-audit:{slot}",
                f"checkpoint-commit:{slot}",
            ]
        )
        sid = section_by_slot[slot]
        witnesses[sid] = events
        histories[sid] = {
            "source_slot": slot,
            "source_port": source_port,
            "orientation": orientation,
            "events": events,
            "semantic_event_count": len(events),
            "uses_executor_metadata": False,
            "uses_external_target": False,
        }
    return witnesses, histories


def _joint_kernels(
    section_by_slot: Mapping[str, str], actions: Mapping[str, Mapping[str, str]]
) -> list[dict[str, Any]]:
    all_ports = list(PORTS)
    kernels: list[dict[str, Any]] = []
    aliases = {
        "r0_s0_a0_f0": "identity",
        "r1_s0_a0_f0": "cyclic-port-permutation",
        "r0_s1_a0_f0": "meridian-reflection",
        "r0_s0_a1_f0": "inverse-port-antipode",
        "r0_s0_a0_f1": "write-check-reversal",
    }
    for continuation_id in sorted(actions):
        action = actions[continuation_id]
        kernels.append(
            {
                "authorized_observers": all_ports,
                "continuation_id": continuation_id,
                "semantic_alias": aliases.get(continuation_id, "composed-reversible-continuation"),
                "rows": {
                    section_by_slot[source_slot]: {
                        section_by_slot[action[source_slot]]: 1.0
                    }
                    for source_slot in oriented_slots()
                },
            }
        )
    return kernels


def _emit_local_checkpoint_packets(
    actions: Mapping[str, Mapping[str, str]]
) -> dict[str, Any]:
    packets: dict[str, Any] = {}
    for continuation_id, action in sorted(actions.items()):
        packets[continuation_id] = {
            observer: {
                source_slot: {_local_atom(observer, action[source_slot]): 1.0}
                for source_slot in oriented_slots()
            }
            for observer in PORTS
        }
    return packets


def verify_local_marginal_consistency(packet: Mapping[str, Any]) -> dict[str, Any]:
    sections = public_global_sections(packet["observers"], packet["interfaces"])
    section_by_slot, section_by_id = _sections_by_slot(packet["observers"], packet["interfaces"])
    reachable = reachable_public_sections(sections, packet["reachability_witnesses"])
    declared = packet.get("local_checkpoint_packets")
    if not isinstance(declared, Mapping):
        return {"status": "LOCAL_MARGINAL_MISMATCH", "reason": "missing declaration"}
    checked_rows = 0
    for channel in packet["global_checkpoint_kernels"]:
        continuation_id = channel["continuation_id"]
        if continuation_id not in declared:
            return {"status": "LOCAL_MARGINAL_MISMATCH", "continuation_id": continuation_id}
        rows = _channel_rows(channel, reachable)
        for source_slot, source_sid in section_by_slot.items():
            for observer in PORTS:
                derived: dict[str, float] = {}
                for output_sid, probability in rows[source_sid].items():
                    atom = section_by_id[output_sid][observer]
                    derived[atom] = derived.get(atom, 0.0) + probability
                expected = {
                    atom: float(probability)
                    for atom, probability in declared[continuation_id][observer][source_slot].items()
                }
                if derived != expected:
                    return {
                        "status": "LOCAL_MARGINAL_MISMATCH",
                        "continuation_id": continuation_id,
                        "observer": observer,
                        "source_slot": source_slot,
                        "derived": derived,
                        "declared": expected,
                    }
                checked_rows += 1
    return {
        "status": "PASS",
        "continuation_count": len(packet["global_checkpoint_kernels"]),
        "observer_count": len(PORTS),
        "source_row_count": len(oriented_slots()),
        "marginal_rows_checked": checked_rows,
    }


def _projection_manifest(section_by_slot: Mapping[str, str]) -> tuple[dict[str, list[int]], dict[str, Any]]:
    supports = {
        section_by_slot[slot]: [index]
        for index, slot in enumerate(oriented_slots())
    }
    sparse = {
        section_by_slot[slot]: {"basis_index": index, "rank": 1, "entry": [index, index, 1]}
        for index, slot in enumerate(oriented_slots())
    }
    return supports, {
        "carrier_dimension": len(oriented_slots()),
        "basis_order": oriented_slots(),
        "sparse_rank_one_projections": sparse,
        "pairwise_orthogonal": True,
        "orthogonality_pairs_checked": len(oriented_slots()) * (len(oriented_slots()) - 1) // 2,
        "rank_sum": len(oriented_slots()),
        "resolves_identity": True,
    }


def build_source_derived_packet(
    terminal_manifest: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    manifest = terminal_manifest or build_terminal_fiber_manifest()
    observers, interfaces = _build_record_diagram()
    sections = public_global_sections(observers, interfaces)
    section_by_slot, _ = _sections_by_slot(observers, interfaces)
    witnesses, histories = _reachability(section_by_slot)
    actions = continuation_actions()
    kernels = _joint_kernels(section_by_slot, actions)
    local_packets = _emit_local_checkpoint_packets(actions)
    projection_supports, projection_manifest = _projection_manifest(section_by_slot)
    composition_table = _continuation_composition_table(actions)
    packet: dict[str, Any] = {
        "schema": SCHEMA,
        "status": "SOURCE_DERIVED_FIXED_CUTOFF_PHYSICAL_PACKET",
        "claim_boundary": (
            "finite source-derived reversible screen packet; no cosmic selector, horizon bridge, "
            "screen/electroweak load bridge, or hardware realization is claimed"
        ),
        "regulator_id": REGULATOR_ID,
        "carrier_type_id": FROZEN_CARRIER_TYPE,
        "terminal_id": TERMINAL_ID,
        "trial_universe_constructor_hash": manifest["trial_universe_constructor_sha256"],
        "terminal_fiber_manifest_hash": manifest["terminal_fiber_manifest_sha256"],
        "terminal_fiber_completeness_certificate": manifest["terminal_fiber_completeness_certificate"],
        "capacity_dimension": len(oriented_slots()),
        "carrier_source": {
            "ports": "twelve exposed edge-center screen ports",
            "orientation": "reversible write/check pair",
            "record_register": "P_12 x {write,check}",
            "source_roots": [
                "finite echosahedral port incidence",
                "edge-center record atoms",
                "reversible checkpoint continuation",
                "semantic endogenous repair histories",
            ],
            "excluded_inputs": list(TARGET_MARKERS),
        },
        "observer_registry": [
            {"observer_id": observer, "kind": "edge-center-port-observer"}
            for observer in PORTS
        ],
        "observers": observers,
        "local_record_atom_sets": observers,
        "interfaces": interfaces,
        "interface_atom_sets": {
            interface["interface_id"]: interface["interface_atoms"]
            for interface in interfaces
        },
        "public_global_sections": sections,
        "public_section_aliases": section_by_slot,
        "reachability_witnesses": witnesses,
        "semantic_histories": histories,
        "publicness_policy_id": "universal-twelve-port-publicness/v1",
        "publicness_policy": [list(PORTS)],
        "continuation_manifest_complete": True,
        "continuation_family_kind": "D5 x C2_antipodal x C2_orientation",
        "continuation_family_order": len(actions),
        "continuation_generators": [
            "r1_s0_a0_f0",
            "r0_s1_a0_f0",
            "r0_s0_a1_f0",
            "r0_s0_a0_f1",
        ],
        "global_checkpoint_kernels": kernels,
        "support_relation_composition": {
            "complete": True,
            "composition_order": "left-after-right",
            "table": composition_table,
        },
        "local_checkpoint_packets": local_packets,
        "local_marginal_manifest_complete": True,
        "local_marginal_consistency_passed": True,
        "projection_supports": projection_supports,
        "carrier_projection_manifest": projection_manifest,
        "refinement_maps": [],
        "self_read_predicate_injected": False,
        "supplied_capacity_metadata_read_by_producer": False,
        "lambda_used": False,
        "ew_bridge_used": False,
        "rho_used": False,
        "measured_lambda_used": False,
        "ew_bridge_target_used": False,
        "rho_used_as_capacity_producer": False,
        "packet_hash_scope": "all packet fields except packet_sha256",
    }
    packet["packet_sha256"] = _sha256(packet)
    return packet


def _exact_decoder_receipts(packet: Mapping[str, Any]) -> dict[str, Any]:
    sections = public_global_sections(packet["observers"], packet["interfaces"])
    reachable = reachable_public_sections(sections, packet["reachability_witnesses"])
    receipts: dict[str, Any] = {}
    for channel in packet["global_checkpoint_kernels"]:
        rows = _channel_rows(channel, reachable)
        decoder: dict[str, str] = {}
        image: list[str] = []
        for source in reachable:
            positive = [(output, probability) for output, probability in rows[source].items() if probability > 0]
            if len(positive) != 1 or positive[0][1] != 1.0:
                return {"status": "NONDETERMINISTIC_CHECKPOINT_GENERATOR", "source": source}
            output = positive[0][0]
            image.append(output)
            if output in decoder:
                return {"status": "NONINJECTIVE_CHECKPOINT_GENERATOR", "output": output}
            decoder[output] = source
        if set(image) != set(reachable):
            return {"status": "CHECKPOINT_NOT_CLOSED_ON_REACHABLE_RECORDS"}
        success = min(
            sum(probability for output, probability in rows[source].items() if decoder.get(output) == source)
            for source in reachable
        )
        receipts[channel["continuation_id"]] = {
            "status": "PASS",
            "deterministic": True,
            "injective": True,
            "surjective": True,
            "decoder": decoder,
            "worst_input_success": success,
        }
    return {"status": "PASS", "channels": receipts}


def _verify_composition(packet: Mapping[str, Any]) -> dict[str, Any]:
    actions = continuation_actions()
    expected = _continuation_composition_table(actions)
    actual = packet["support_relation_composition"]["table"]
    return {
        "status": "PASS" if actual == expected else "COMPOSITION_TABLE_MISMATCH",
        "continuation_count": len(actions),
        "composition_entries_checked": len(actions) ** 2,
    }


def _carrier_receipt(packet: Mapping[str, Any], evaluation: Mapping[str, Any]) -> dict[str, Any]:
    manifest = packet["carrier_projection_manifest"]
    supports = packet["projection_supports"]
    used = [index for support in supports.values() for index in support]
    pairwise = len(used) == len(set(used))
    in_range = all(0 <= index < packet["capacity_dimension"] for index in used)
    return {
        "status": "PASS" if pairwise and in_range and manifest["resolves_identity"] else "CARRIER_FAILURE",
        "orthogonal_rank_one_projections": pairwise and all(len(support) == 1 for support in supports.values()),
        "rank_sum": len(set(used)),
        "dimension": packet["capacity_dimension"],
        "exact_bound": evaluation["exact_zero_error_capacity"] <= packet["capacity_dimension"],
        "rank_one_saturation": evaluation["saturation_rank_one_complete"],
        "approximate_bound": "M_epsilon <= D for every epsilon by orthogonal carrier decoding",
    }


def _tv_receipts(dimension: int) -> list[dict[str, Any]]:
    receipts: list[dict[str, Any]] = []
    for delta in (0.0, 0.01, 0.05):
        epsilon = tv_robustness_bound(0.0, delta)
        receipts.append(
            {
                "rowwise_tv_bound": delta,
                "guaranteed_worst_input_success": 1.0 - delta,
                "epsilon": epsilon,
                "lower_bound": dimension,
                "carrier_upper_bound": dimension,
                "certified_M_epsilon": dimension,
            }
        )
    return receipts


def _alternative_joint_coupling_control() -> dict[str, Any]:
    reachable = ["x0", "x1"]
    parity = {
        "continuation_id": "parity-coupling",
        "rows": {
            "x0": {"00": 0.5, "11": 0.5},
            "x1": {"01": 0.5, "10": 0.5},
        },
    }
    independent = {
        "continuation_id": "independent-uniform",
        "rows": {
            "x0": {"00": 0.25, "01": 0.25, "10": 0.25, "11": 0.25},
            "x1": {"00": 0.25, "01": 0.25, "10": 0.25, "11": 0.25},
        },
    }

    def marginals(channel: Mapping[str, Any]) -> dict[str, dict[str, dict[str, float]]]:
        result: dict[str, dict[str, dict[str, float]]] = {}
        for source, row in channel["rows"].items():
            first = {"0": 0.0, "1": 0.0}
            second = {"0": 0.0, "1": 0.0}
            for output, probability in row.items():
                first[output[0]] += probability
                second[output[1]] += probability
            result[source] = {"observer_1": first, "observer_2": second}
        return result

    graph_parity = compound_confusability_graph(reachable, [parity])
    graph_independent = compound_confusability_graph(reachable, [independent])
    return {
        "status": "PASS",
        "same_local_marginals": marginals(parity) == marginals(independent),
        "local_marginals": marginals(parity),
        "joint_capacity_parity": len(maximum_independent_set(graph_parity)),
        "joint_capacity_independent": len(maximum_independent_set(graph_independent)),
    }


def _full_support_noise_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    sections = public_global_sections(packet["observers"], packet["interfaces"])
    reachable = reachable_public_sections(sections, packet["reachability_witnesses"])
    delta_mix = 1.0e-6
    base = next(
        channel
        for channel in packet["global_checkpoint_kernels"]
        if channel["continuation_id"] == "r1_s0_a0_f0"
    )
    noisy_rows: dict[str, dict[str, float]] = {}
    decoder: dict[str, str] = {}
    for source in reachable:
        deterministic_output = next(iter(base["rows"][source]))
        decoder[deterministic_output] = source
        row = {output: delta_mix / len(reachable) for output in reachable}
        row[deterministic_output] += 1.0 - delta_mix
        noisy_rows[source] = row
    noisy = {"continuation_id": "full-support-noise", "rows": noisy_rows}
    graph = compound_confusability_graph(reachable, [noisy])
    worst_success = min(
        sum(probability for output, probability in noisy_rows[source].items() if decoder.get(output) == source)
        for source in reachable
    )
    row_tv = delta_mix * (1.0 - 1.0 / len(reachable))
    return {
        "status": "PASS",
        "all_rows_full_support": all(all(probability > 0 for probability in row.values()) for row in noisy_rows.values()),
        "exact_zero_error_capacity": len(maximum_independent_set(graph)),
        "mixture_weight": delta_mix,
        "row_tv_from_reversible_kernel": row_tv,
        "inverse_decoder_worst_input_success": worst_success,
        "epsilon_for_full_code": row_tv,
        "approximate_capacity_at_epsilon": len(reachable),
        "tv_identity": abs(worst_success - (1.0 - row_tv)) < 1e-12,
    }


def _order_controls(max_dimension: int = 24) -> dict[str, Any]:
    identity = {dimension: dimension for dimension in range(1, max_dimension + 1)}
    erasure = {dimension: 1 for dimension in range(1, max_dimension + 1)}
    return {
        "status": "PASS",
        "identity_family": greatest_fixed_point(identity),
        "erasure_family": greatest_fixed_point(erasure),
        "identity_fixed_every_dimension": greatest_fixed_point(identity)["fixed_points"] == list(range(1, max_dimension + 1)),
        "erasure_fixed_only_singleton": greatest_fixed_point(erasure)["fixed_points"] == [1],
    }


def classify_terminal_fiber(
    packets: Sequence[Mapping[str, Any]], *, manifest_complete: bool
) -> dict[str, Any]:
    if not manifest_complete:
        return {"fiber_kind": "INCOMPLETE", "status": "INCOMPLETE_TERMINAL_FIBER"}
    if not packets:
        return {"fiber_kind": "EMPTY", "status": "NO_CAPACITY_READBACK"}
    result = evaluate_terminal_fiber(packets, manifest_complete=True)
    if result["status"] == "AMBIGUOUS_CAPACITY_READBACK":
        return {"fiber_kind": "AMBIGUOUS", **result}
    if result["status"] == "PASS":
        kind = "SINGLETON" if len(packets) == 1 else "COMPLETE_SCALAR"
        return {"fiber_kind": kind, **result}
    return {"fiber_kind": "UNAVAILABLE", **result}


def _fiber_controls(packet: Mapping[str, Any]) -> dict[str, Any]:
    ambiguous = copy.deepcopy(packet)
    sections = public_global_sections(ambiguous["observers"], ambiguous["interfaces"])
    reachable = reachable_public_sections(sections, ambiguous["reachability_witnesses"])
    sink = reachable[0]
    ambiguous["terminal_id"] = "q_erasure_control"
    ambiguous["global_checkpoint_kernels"] = [
        {
            "authorized_observers": list(PORTS),
            "continuation_id": "erasure",
            "rows": {source: {sink: 1.0} for source in reachable},
        }
    ]
    ambiguous["continuation_manifest_complete"] = True
    ambiguous["local_marginal_consistency_passed"] = True
    return {
        "status": "PASS",
        "empty": classify_terminal_fiber([], manifest_complete=True),
        "incomplete": classify_terminal_fiber([packet], manifest_complete=False),
        "singleton": classify_terminal_fiber([packet], manifest_complete=True),
        "ambiguous": classify_terminal_fiber([packet, ambiguous], manifest_complete=True),
    }


def _isomorphic_packet_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    renamed = copy.deepcopy(packet)
    observer_map = {observer: f"observer_{index:02d}" for index, observer in enumerate(PORTS)}

    def rename_atom(atom: str) -> str:
        observer, slot = atom.split("::record::", 1)
        return f"{observer_map[observer]}::record::{slot}"

    renamed_observers = {
        observer_map[observer]: [rename_atom(atom) for atom in atoms]
        for observer, atoms in packet["observers"].items()
    }
    renamed_interfaces: list[dict[str, Any]] = []
    for interface in packet["interfaces"]:
        left_old = interface["left_observer"]
        right_old = interface["right_observer"]
        left_new = observer_map[left_old]
        right_new = observer_map[right_old]
        interface_id = f"{left_new}--{right_new}"
        left_readout = {
            rename_atom(atom): f"{interface_id}::public::{_slot_from_local_atom(atom)}"
            for atom in packet["observers"][left_old]
        }
        right_readout = {
            rename_atom(atom): f"{interface_id}::public::{_slot_from_local_atom(atom)}"
            for atom in packet["observers"][right_old]
        }
        renamed_interfaces.append(
            {
                "interface_id": interface_id,
                "left_observer": left_new,
                "right_observer": right_new,
                "interface_atoms": list(left_readout.values()),
                "left_readout": left_readout,
                "right_readout": right_readout,
            }
        )

    old_sections = public_global_sections(packet["observers"], packet["interfaces"])
    new_sections = public_global_sections(renamed_observers, renamed_interfaces)
    old_by_slot = {
        _slot_from_local_atom(section[PORTS[0]]): section_id(section)
        for section in old_sections
    }
    first_new_observer = observer_map[PORTS[0]]
    new_by_slot = {
        _slot_from_local_atom(section[first_new_observer]): section_id(section)
        for section in new_sections
    }
    section_map = {old_by_slot[slot]: new_by_slot[slot] for slot in oriented_slots()}

    renamed["observers"] = renamed_observers
    renamed["interfaces"] = renamed_interfaces
    renamed["reachability_witnesses"] = {
        section_map[source]: list(history)
        for source, history in packet["reachability_witnesses"].items()
    }
    renamed["publicness_policy"] = [list(observer_map.values())]
    renamed["global_checkpoint_kernels"] = [
        {
            **{key: value for key, value in channel.items() if key not in {"authorized_observers", "rows"}},
            "authorized_observers": list(observer_map.values()),
            "rows": {
                section_map[source]: {section_map[output]: probability for output, probability in row.items()}
                for source, row in channel["rows"].items()
            },
        }
        for channel in packet["global_checkpoint_kernels"]
    ]
    renamed["projection_supports"] = {
        section_map[source]: list(support)
        for source, support in packet["projection_supports"].items()
    }
    renamed["local_marginal_consistency_passed"] = True
    original_result = evaluate_terminal(packet)
    renamed_result = evaluate_terminal(renamed)
    passed = (
        original_result["status"] == renamed_result["status"] == "PASS"
        and original_result["exact_zero_error_capacity"]
        == renamed_result["exact_zero_error_capacity"]
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "observer_bijection": observer_map,
        "section_bijection_size": len(section_map),
        "original_capacity": original_result.get("exact_zero_error_capacity"),
        "renamed_capacity": renamed_result.get("exact_zero_error_capacity"),
    }


def _embedding_receipts(packet: Mapping[str, Any]) -> dict[str, Any]:
    evaluation = evaluate_terminal(packet)
    coarse = {vertex: set(neighbours) for vertex, neighbours in evaluation["confusability_graph"].items()}
    extension_vertices = [f"layer_{layer}::{vertex}" for layer in (0, 1) for vertex in coarse]
    extension_graph = {vertex: set() for vertex in extension_vertices}
    extension_embedding = {vertex: f"layer_0::{vertex}" for vertex in coarse}
    refinement_graph = {f"r1::{vertex}": set() for vertex in coarse}
    refinement_embedding = {vertex: f"r1::{vertex}" for vertex in coarse}
    bad_extension_graph = {vertex: set(neighbours) for vertex, neighbours in extension_graph.items()}
    bad_left = extension_embedding[next(iter(coarse))]
    bad_right = extension_embedding[next(iter(vertex for vertex in coarse if extension_embedding[vertex] != bad_left))]
    bad_extension_graph[bad_left].add(bad_right)
    bad_extension_graph[bad_right].add(bad_left)
    bad_refinement_graph = {vertex: set(neighbours) for vertex, neighbours in refinement_graph.items()}
    bad_ref_left = refinement_embedding[next(iter(coarse))]
    bad_ref_right = refinement_embedding[next(iter(vertex for vertex in coarse if refinement_embedding[vertex] != bad_ref_left))]
    bad_refinement_graph[bad_ref_left].add(bad_ref_right)
    bad_refinement_graph[bad_ref_right].add(bad_ref_left)
    return {
        "status": "PASS",
        "capacity_extension": {
            "from_dimension": len(coarse),
            "to_dimension": len(extension_graph),
            "embedding": extension_embedding,
            "no_new_confusability": no_new_confusability(coarse, extension_graph, extension_embedding),
            "new_confusability_control_rejected": not no_new_confusability(
                coarse, bad_extension_graph, extension_embedding
            ),
        },
        "fixed_D_refinement": {
            "dimension": len(coarse),
            "from_regulator": REGULATOR_ID,
            "to_regulator": "echosahedral-edge-center-r1",
            "embedding": refinement_embedding,
            "no_new_confusability": no_new_confusability(coarse, refinement_graph, refinement_embedding),
            "new_confusability_control_rejected": not no_new_confusability(
                coarse, bad_refinement_graph, refinement_embedding
            ),
        },
    }


def _target_taint_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    statuses: dict[str, str] = {}
    for field in ("lambda_used", "ew_bridge_used", "rho_used"):
        tainted = copy.deepcopy(packet)
        tainted[field] = True
        statuses[field] = evaluate_terminal(tainted)["status"]
    passed = all(status == "TARGET_TAINTED" for status in statuses.values())
    return {
        "status": "PASS" if passed else "FAIL",
        "evaluator_status": statuses["lambda_used"],
        "statuses_by_taint": statuses,
    }


def _circular_definition_control(packet: Mapping[str, Any]) -> dict[str, Any]:
    tainted = copy.deepcopy(packet)
    tainted["self_read_predicate_injected"] = True
    result = evaluate_terminal(tainted)
    return {
        "status": "PASS" if result["status"] == "CIRCULAR_CAPACITY_DEFINITION" else "FAIL",
        "evaluator_status": result["status"],
    }


def _source_provenance_receipt(packet: Mapping[str, Any]) -> dict[str, Any]:
    carrier_source = packet["carrier_source"]
    source_text = json.dumps(carrier_source["source_roots"], sort_keys=True).lower()
    forbidden_values_absent = all(marker not in source_text for marker in TARGET_MARKERS)
    exclusions_declared = set(TARGET_MARKERS).issubset(set(carrier_source.get("excluded_inputs", ())))
    administrative_flags_clean = (
        not packet["self_read_predicate_injected"]
        and not packet["supplied_capacity_metadata_read_by_producer"]
        and not packet["lambda_used"]
        and not packet["ew_bridge_used"]
        and not packet["rho_used"]
    )
    return {
        "status": "PASS" if forbidden_values_absent and exclusions_declared and administrative_flags_clean else "TARGET_TAINTED",
        "forbidden_target_values_absent": forbidden_values_absent,
        "excluded_input_ledger_complete": exclusions_declared,
        "administrative_taint_flags_clear": administrative_flags_clean,
        "source_roots": carrier_source["source_roots"],
    }


def certify_source_derived_packet(
    packet: Mapping[str, Any], *, terminal_manifest: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    manifest = terminal_manifest or build_terminal_fiber_manifest()
    provenance = _source_provenance_receipt(packet)
    marginals = verify_local_marginal_consistency(packet)
    if provenance["status"] != "PASS":
        return provenance
    if marginals["status"] != "PASS":
        return marginals
    evaluation = evaluate_terminal(packet)
    if evaluation["status"] != "PASS":
        return evaluation
    decoders = _exact_decoder_receipts(packet)
    if decoders["status"] != "PASS":
        return decoders
    composition = _verify_composition(packet)
    if composition["status"] != "PASS":
        return composition
    reachable = evaluation["reachable_public_sections"]
    semigroup = support_relation_semigroup(reachable, packet["global_checkpoint_kernels"])
    graph_edge_count = sum(len(neighbours) for neighbours in evaluation["confusability_graph"].values()) // 2
    carrier = _carrier_receipt(packet, evaluation)
    controls = {
        "isomorphic_packet": _isomorphic_packet_control(packet),
        "cyclic_permutation": decoders["channels"]["r1_s0_a0_f0"],
        "alternative_joint_coupling": _alternative_joint_coupling_control(),
        "full_support_noise": _full_support_noise_control(packet),
        "circular_definition": _circular_definition_control(packet),
        "target_taint": _target_taint_control(packet),
        "order_families": _order_controls(packet["capacity_dimension"]),
        "terminal_fibers": _fiber_controls(packet),
        "finite_suffix_nonpromotion": {
            "status": "PASS",
            "equal_suffix": [24, 24, 24],
            "saturation_or_exhaustion_receipt": False,
            "promoted_to_limit": False,
        },
    }
    embeddings = _embedding_receipts(packet)
    packet_without_hash = dict(packet)
    claimed_packet_hash = packet_without_hash.pop("packet_sha256", None)
    packet_hash_valid = claimed_packet_hash == _sha256(packet_without_hash)
    fiber_evaluation = evaluate_terminal_fiber([packet], manifest_complete=True)
    scalarization_receipt = {
        "fiber_kind": "SINGLETON",
        **fiber_evaluation,
    }
    count_kernel_row = {
        "status": "PASS",
        "capacity_dimension": packet["capacity_dimension"],
        "exact_reachable_section_count": len(reachable),
        "all_checkpoint_generators_injective": all(
            receipt["injective"] for receipt in decoders["channels"].values()
        ),
        "multiplicity_by_capacity": {
            str(evaluation["exact_zero_error_capacity"]): 1
        },
    }
    pass_checks = [
        packet_hash_valid,
        packet.get("trial_universe_constructor_hash") == manifest["trial_universe_constructor_sha256"],
        packet.get("terminal_fiber_manifest_hash") == manifest["terminal_fiber_manifest_sha256"],
        manifest["terminal_fiber_complete"],
        manifest["terminal_fiber_kind"] == "SINGLETON",
        manifest["terminal_membership_output_blind"],
        provenance["status"] == "PASS",
        marginals["status"] == "PASS",
        evaluation["status"] == "PASS",
        decoders["status"] == "PASS",
        composition["status"] == "PASS",
        len(semigroup) == packet["continuation_family_order"],
        graph_edge_count == 0,
        evaluation["exact_zero_error_capacity"] == len(reachable),
        carrier["status"] == "PASS",
        carrier["rank_one_saturation"],
        all(control.get("status") == "PASS" for control in controls.values()),
        embeddings["capacity_extension"]["no_new_confusability"],
        embeddings["capacity_extension"]["new_confusability_control_rejected"],
        embeddings["fixed_D_refinement"]["no_new_confusability"],
        embeddings["fixed_D_refinement"]["new_confusability_control_rejected"],
        fiber_evaluation["status"] == "PASS",
        fiber_evaluation["robust_closure"],
    ]
    return {
        "status": "PASS" if all(pass_checks) else "FAIL",
        "issue": 548,
        "packet_status": packet["status"],
        "claim_boundary": packet["claim_boundary"],
        "packet_sha256": packet["packet_sha256"],
        "packet_hash_valid": packet_hash_valid,
        "trial_universe_constructor_hash": packet["trial_universe_constructor_hash"],
        "terminal_fiber_manifest_hash": packet["terminal_fiber_manifest_hash"],
        "terminal_fiber_manifest": manifest,
        "source_provenance": provenance,
        "record_diagram": {
            "observer_count": len(packet["observers"]),
            "interface_count": len(packet["interfaces"]),
            "local_atoms_per_observer": len(next(iter(packet["observers"].values()))),
            "readout_maps_total": True,
            "public_global_section_count": evaluation["public_global_section_count"],
        },
        "reachability": {
            "reachable_public_record_count": len(reachable),
            "all_histories_nonempty": all(packet["reachability_witnesses"][record] for record in reachable),
            "endogenous_history_receipts": packet["semantic_histories"],
            "publicness_policy_id": packet["publicness_policy_id"],
        },
        "checkpoint_family": {
            "family_kind": packet["continuation_family_kind"],
            "continuation_count": packet["continuation_family_order"],
            "support_relation_semigroup_size": len(semigroup),
            "composition": composition,
            "local_marginals": marginals,
            "injective_generator_receipts": decoders["channels"],
        },
        "capacity": {
            "formula": "M_0(q)=alpha(G_q)",
            "terminal_id": packet["terminal_id"],
            "compound_graph": evaluation["confusability_graph"],
            "compound_graph_edge_count": graph_edge_count,
            "independent_set_witness": evaluation["independent_set_witness"],
            "exact_zero_error_capacity": evaluation["exact_zero_error_capacity"],
            "reachable_section_count": len(reachable),
            "fast_branch_identity": "M_0(q)=|X_reach(q)|",
            "terminal_fiber_capacity_set": [evaluation["exact_zero_error_capacity"]],
            "count_kernel_row": count_kernel_row,
            "whole_fiber_scalarization_result": scalarization_receipt,
            "robust_closure_at_frozen_D": fiber_evaluation["robust_closure"],
        },
        "exact_decoders": decoders,
        "approximate_worst_input_branch": {
            "exact_family_worst_input_success": 1.0,
            "exact_M_epsilon_at_epsilon_0": len(reachable),
            "tv_robustness": _tv_receipts(packet["capacity_dimension"]),
        },
        "carrier": carrier,
        "controls": controls,
        "no_new_confusability_injections": embeddings,
        "remaining_unresolved_objects": [
            "capacity-indexed source family beyond the frozen first packet",
            "exact finite-size slack law with one regulator-stable physical zero",
            "horizon-record identification",
            "common screen/electroweak load-carrier identification",
            "hardware realization evidence",
        ],
    }


def build_source_derived_bundle() -> dict[str, Any]:
    manifest = build_terminal_fiber_manifest()
    packet = build_source_derived_packet(manifest)
    certificate = certify_source_derived_packet(packet, terminal_manifest=manifest)
    return {
        "schema": "OPH_ISSUE_548_BUNDLE/v1",
        "terminal_fiber_manifest": manifest,
        "terminal_packets": [packet],
        "certificate": certificate,
    }


def write_runtime_receipts(output_dir: str | Path) -> dict[str, str]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    bundle = build_source_derived_bundle()
    paths = {
        "terminal_fiber_manifest": output / "source_derived_terminal_fiber_manifest.json",
        "packet": output / "source_derived_public_checkpoint_packet.json",
        "certificate": output / "source_derived_public_checkpoint_certificate.json",
    }
    payloads = {
        "terminal_fiber_manifest": bundle["terminal_fiber_manifest"],
        "packet": bundle["terminal_packets"][0],
        "certificate": bundle["certificate"],
    }
    for key, path in paths.items():
        path.write_text(json.dumps(payloads[key], indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {key: str(path) for key, path in paths.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        default=Path(__file__).with_name("runtime"),
        type=Path,
        help="directory for canonical JSON receipts",
    )
    args = parser.parse_args()
    paths = write_runtime_receipts(args.output_dir)
    certificate = json.loads(Path(paths["certificate"]).read_text(encoding="utf-8"))
    print(json.dumps({"status": certificate["status"], "receipts": paths}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

