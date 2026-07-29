#!/usr/bin/env python3
"""Target-clean capacity-indexed completion family for issue #551.

The issue #548 packet fixes one twenty-four-record screen packet. This module
tests a narrower declared continuation grammar over positive integer rungs
``k`` and four completion rules:

* reversible identity;
* copy-collapse erasure;
* a two-class cap;
* hidden spectator multiplicity.

The executable rows share the twelve-port record register, the forty
reversible slot actions, publicness, an oriented-record fiber product, and a
fixed-packet semantic projection. Their zero-error capacities and slack-zero
sets differ. The result proves non-entailment for the declared base-agreement,
positivity, and carrier-bound completion class. It does not prove that a
complete A1--A3 lift of every fixed-packet atom, terminal fiber, joint kernel,
and A3 feasible-set construction is nonidentifiable.

Measured cosmology, horizon size, electroweak values, desired capacities, and
external fits are absent from the producer.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from source_derived_public_checkpoint_packet import (
    ORIENTATIONS,
    PORTS,
    continuation_actions,
    icosahedral_edges,
    oriented_slots,
)


SCHEMA = "oph.capacity_indexed_source_family_projection.v1"
CERTIFICATE_SCHEMA = "oph.capacity_indexed_source_family_certificate.v1"
SCIENTIFIC_VERDICT = "BOUNDED_COMPLETION_CLASS_NONIDENTIFIABLE"
SOURCE_RULE_ID = "oph.public-record-capacity.branch-completion-family.v1"
BASE_PUBLIC_ATOMS = 24
SAMPLE_RUNGS = (1, 2, 3, 4)
SPECTATOR_MULTIPLICITIES = (1, 2, 3)
BRANCH_IDS = (
    "reversible_identity",
    "copy_collapse_erasure",
    "capped_two_class",
    "hidden_spectator",
)
TARGET_CLEANLINESS = {
    "measured_cosmological_constant_read": False,
    "observed_horizon_radius_read": False,
    "electroweak_target_read": False,
    "desired_capacity_read": False,
    "external_fit_read": False,
}
HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
PROJECTION_PATH = RUNTIME / "capacity_indexed_source_family_projection.json"
CERTIFICATE_PATH = RUNTIME / "capacity_indexed_source_family_certificate.json"
FIXED_PACKET_PATH = RUNTIME / "source_derived_public_checkpoint_packet.json"
SPEC_PATH = HERE / "F_READBACK_SPEC.md"


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


def _fixed_packet_semantic_projection() -> dict[str, Any]:
    packet = json.loads(FIXED_PACKET_PATH.read_text(encoding="utf-8"))
    return {
        "schema": packet["schema"],
        "carrier_type_id": packet["carrier_type_id"],
        "capacity_dimension": packet["capacity_dimension"],
        "port_ids": [row["observer_id"] for row in packet["observer_registry"]],
        "interface_ids": sorted(row["interface_id"] for row in packet["interfaces"]),
        "publicness_policy_id": packet["publicness_policy_id"],
        "continuation_family_kind": packet["continuation_family_kind"],
        "continuation_family_order": packet["continuation_family_order"],
        "projection_rank_sum": packet["carrier_projection_manifest"]["rank_sum"],
        "target_flags": {
            "self_read_predicate_injected": packet["self_read_predicate_injected"],
            "lambda_used": packet["lambda_used"],
            "ew_bridge_used": packet["ew_bridge_used"],
            "rho_used": packet["rho_used"],
        },
    }


def upstream_pins() -> dict[str, str]:
    return {
        "fixed_packet_projection_sha256": tagged_sha256(
            canonical_json_bytes(_fixed_packet_semantic_projection())
        ),
        "fixed_packet_verifier_contract_sha256": tagged_sha256(
            SPEC_PATH.read_bytes()
        ),
    }


def shared_source() -> dict[str, Any]:
    return {
        "source_rule_id": SOURCE_RULE_ID,
        "base_public_atoms": BASE_PUBLIC_ATOMS,
        "rung_parameter": "positive-integer-k",
        "sample_rungs": list(SAMPLE_RUNGS),
        "spectator_multiplicities": list(SPECTATOR_MULTIPLICITIES),
        "branch_completion_scope": (
            "base-agreement-positive-carrier-bound-completion-class"
        ),
        "full_a1_a3_packet_lift_required": True,
        "pin_scope": "semantic-contract-identifiers-not-file-custody",
    }


def source_signature() -> str:
    return tagged_sha256(canonical_json_bytes(shared_source()))


def _require_parameters(k: int, spectator_multiplicity: int) -> None:
    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer")
    if (
        not isinstance(spectator_multiplicity, int)
        or isinstance(spectator_multiplicity, bool)
        or spectator_multiplicity < 1
    ):
        raise ValueError("spectator multiplicity must be a positive integer")


def _record_id(slot: str, copy_index: int, spectator_index: int | None) -> str:
    if spectator_index is None:
        return f"{slot}|copy={copy_index}"
    return f"{slot}|copy={copy_index}|spectator={spectator_index}"


def _record_coordinates(
    branch_id: str,
    k: int,
    spectator_multiplicity: int,
) -> list[tuple[str, int, int | None]]:
    _require_parameters(k, spectator_multiplicity)
    if branch_id not in BRANCH_IDS:
        raise ValueError(f"unknown branch {branch_id}")
    spectator_values: Sequence[int | None]
    if branch_id == "hidden_spectator":
        spectator_values = tuple(range(spectator_multiplicity))
    else:
        if spectator_multiplicity != 1:
            raise ValueError("only the hidden-spectator branch accepts s != 1")
        spectator_values = (None,)
    return [
        (slot, copy_index, spectator_index)
        for slot in oriented_slots()
        for copy_index in range(k)
        for spectator_index in spectator_values
    ]


def _collapse_output(
    branch_id: str,
    slot: str,
    copy_index: int,
    spectator_index: int | None,
) -> tuple[str, int]:
    del spectator_index
    if branch_id in {"reversible_identity", "hidden_spectator"}:
        return slot, copy_index
    if branch_id == "copy_collapse_erasure":
        return slot, 0
    if branch_id == "capped_two_class":
        return slot, min(copy_index, 1)
    raise ValueError(f"unknown branch {branch_id}")


def build_capacity_packet(
    branch_id: str,
    k: int,
    *,
    spectator_multiplicity: int = 1,
) -> dict[str, Any]:
    """Emit one exact graph-facing packet from the common all-rung rule."""

    coordinates = _record_coordinates(branch_id, k, spectator_multiplicity)
    records = [
        _record_id(slot, copy_index, spectator_index)
        for slot, copy_index, spectator_index in coordinates
    ]
    outputs = {
        _record_id(slot, copy_index, spectator_index): "|".join(
            map(
                str,
                _collapse_output(
                    branch_id, slot, copy_index, spectator_index
                ),
            )
        )
        for slot, copy_index, spectator_index in coordinates
    }
    projection_rank = (
        1 if branch_id != "hidden_spectator" else spectator_multiplicity
    )
    raw_dimension = len(records)
    public_dimension = BASE_PUBLIC_ATOMS * k
    packet = {
        "schema": "oph.capacity_indexed_public_checkpoint_packet.v1",
        "source_rule_id": SOURCE_RULE_ID,
        "shared_source_signature_sha256": source_signature(),
        "branch_id": branch_id,
        "k": k,
        "spectator_multiplicity": spectator_multiplicity,
        "raw_dimension": raw_dimension,
        "public_dimension": public_dimension,
        "records": records,
        "deterministic_collapse_outputs": outputs,
        "carrier_projection_rank": projection_rank,
        "source_geometry": {
            "ports": list(PORTS),
            "orientations": list(ORIENTATIONS),
            "edge_count": len(icosahedral_edges()),
            "reversible_continuation_order": len(continuation_actions()),
        },
        "publicness_policy": "universal-twelve-port-publicness/v1",
        "target_cleanliness": dict(TARGET_CLEANLINESS),
    }
    packet["packet_sha256"] = tagged_sha256(canonical_json_bytes(packet))
    return packet


def confusability_graph(packet: Mapping[str, Any]) -> dict[str, set[str]]:
    records = list(packet["records"])
    outputs = packet["deterministic_collapse_outputs"]
    fibers: dict[str, list[str]] = {}
    for record in records:
        fibers.setdefault(outputs[record], []).append(record)
    graph = {record: set() for record in records}
    for fiber in fibers.values():
        for index, left in enumerate(fiber):
            for right in fiber[index + 1 :]:
                graph[left].add(right)
                graph[right].add(left)
    return graph


def _clique_component_capacity(graph: Mapping[str, set[str]]) -> dict[str, Any]:
    remaining = set(graph)
    components: list[list[str]] = []
    while remaining:
        root = min(remaining)
        component: set[str] = set()
        frontier = [root]
        while frontier:
            vertex = frontier.pop()
            if vertex in component:
                continue
            component.add(vertex)
            frontier.extend(graph[vertex] - component)
        remaining -= component
        for vertex in component:
            if graph[vertex] & component != component - {vertex}:
                raise AssertionError("generated component must be a clique")
        components.append(sorted(component))
    witness = [component[0] for component in components]
    if any(right in graph[left] for i, left in enumerate(witness) for right in witness[i + 1 :]):
        raise AssertionError("component witness is not independent")
    return {
        "capacity": len(components),
        "independent_set_witness": witness,
        "component_sizes": sorted(len(component) for component in components),
        "upper_bound_reason": "each connected component is a complete graph",
    }


def expected_capacity(
    branch_id: str,
    k: int,
    spectator_multiplicity: int = 1,
) -> int:
    _require_parameters(k, spectator_multiplicity)
    if branch_id == "reversible_identity":
        return BASE_PUBLIC_ATOMS * k
    if branch_id == "copy_collapse_erasure":
        return BASE_PUBLIC_ATOMS
    if branch_id == "capped_two_class":
        return BASE_PUBLIC_ATOMS * min(k, 2)
    if branch_id == "hidden_spectator":
        return BASE_PUBLIC_ATOMS * k
    raise ValueError(f"unknown branch {branch_id}")


def evaluate_capacity_packet(packet: Mapping[str, Any]) -> dict[str, Any]:
    raw = dict(packet)
    claimed_hash = raw.pop("packet_sha256", None)
    if claimed_hash != tagged_sha256(canonical_json_bytes(raw)):
        raise ValueError("packet hash mismatch")
    if raw["source_rule_id"] != SOURCE_RULE_ID:
        raise ValueError("source rule mismatch")
    if raw["shared_source_signature_sha256"] != source_signature():
        raise ValueError("source signature mismatch")
    if set(raw["target_cleanliness"]) != set(TARGET_CLEANLINESS) or any(
        raw["target_cleanliness"].values()
    ):
        raise ValueError("target-tainted packet")
    graph = confusability_graph(packet)
    graph_proof = _clique_component_capacity(graph)
    expected = expected_capacity(
        raw["branch_id"], raw["k"], raw["spectator_multiplicity"]
    )
    if graph_proof["capacity"] != expected:
        raise AssertionError("closed capacity formula failed")
    return {
        "status": "PASS",
        "exact_zero_error_capacity": expected,
        "raw_dimension": raw["raw_dimension"],
        "public_dimension": raw["public_dimension"],
        "slack_zero": expected == raw["raw_dimension"],
        "clique_component_proof": graph_proof,
        "exact_capacity_proof": "matching lower witness and clique-component upper bound",
    }


def _channel_rule(branch_id: str) -> str:
    return {
        "reversible_identity": "output=(port,copy)",
        "copy_collapse_erasure": "output=port",
        "capped_two_class": "output=(port,min(copy,1))",
        "hidden_spectator": "output=(port,copy);spectator-hidden",
    }[branch_id]


def _capacity_formula(branch_id: str) -> str:
    return {
        "reversible_identity": "M0=24*k",
        "copy_collapse_erasure": "M0=24",
        "capped_two_class": "M0=24*min(k,2)",
        "hidden_spectator": "raw_D=24*k*s;M0=24*k",
    }[branch_id]


def _branch_sample_rows(branch_id: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    spectators = (
        SPECTATOR_MULTIPLICITIES
        if branch_id == "hidden_spectator"
        else (1,)
    )
    for k in SAMPLE_RUNGS:
        for spectator in spectators:
            packet = build_capacity_packet(
                branch_id, k, spectator_multiplicity=spectator
            )
            result = evaluate_capacity_packet(packet)
            rows.append(
                {
                    "k": k,
                    "spectator_multiplicity": spectator,
                    "raw_dimension": result["raw_dimension"],
                    "public_dimension": result["public_dimension"],
                    "claimed_capacity_M0": result[
                        "exact_zero_error_capacity"
                    ],
                    "claimed_slack_zero": result["slack_zero"],
                }
            )
    return rows


def build_projection() -> dict[str, Any]:
    common_source = shared_source()
    common_signature = tagged_sha256(canonical_json_bytes(common_source))
    pins = upstream_pins()
    branches: list[dict[str, Any]] = []
    for branch_id in BRANCH_IDS:
        rows = _branch_sample_rows(branch_id)
        branches.append(
            {
                "branch_id": branch_id,
                "capacity_formula": _capacity_formula(branch_id),
                "channel_rule": _channel_rule(branch_id),
                "claimed_bounded_zero_set": [
                    {
                        "k": row["k"],
                        "spectator_multiplicity": row[
                            "spectator_multiplicity"
                        ],
                    }
                    for row in rows
                    if row["claimed_slack_zero"]
                ],
                "sample_rows": rows,
                "shared_source_signature_sha256": common_signature,
                "upstream_pins": pins,
            }
        )
    return {
        "branches": branches,
        "physical_n_closure_promoted": False,
        "schema": SCHEMA,
        "scientific_verdict": SCIENTIFIC_VERDICT,
        "shared_source": common_source,
        "shared_source_signature_sha256": common_signature,
        "target_cleanliness": dict(TARGET_CLEANLINESS),
        "upstream_pins": pins,
    }


def _continuation_receipt(branch_id: str, k: int, spectator: int) -> dict[str, Any]:
    coordinates = _record_coordinates(branch_id, k, spectator)
    actions = continuation_actions()

    def collapse(
        coordinate: tuple[str, int, int | None]
    ) -> tuple[str, int, int | None]:
        slot, copy_index, spectator_index = coordinate
        out_slot, out_copy = _collapse_output(
            branch_id, slot, copy_index, spectator_index
        )
        out_spectator = (
            0 if branch_id == "hidden_spectator" else spectator_index
        )
        return out_slot, out_copy, out_spectator

    idempotent = all(collapse(collapse(row)) == collapse(row) for row in coordinates)
    commuting_rows = 0
    for action in actions.values():
        for row in coordinates:
            slot, copy_index, spectator_index = row
            acted = (action[slot], copy_index, spectator_index)
            left = collapse(acted)
            collapsed = collapse(row)
            right = (action[collapsed[0]], collapsed[1], collapsed[2])
            if left != right:
                raise AssertionError("collapse must commute with slot action")
            commuting_rows += 1
    collapse_is_identity = all(collapse(row) == row for row in coordinates)
    return {
        "status": "PASS",
        "reversible_group_order": len(actions),
        "collapse_idempotent": idempotent,
        "collapse_commutes_with_reversible_group": True,
        "commuting_rows_checked": commuting_rows,
        "complete_generated_semigroup_order": (
            len(actions) if collapse_is_identity else 2 * len(actions)
        ),
        "composition_law": "(g,e) o (h,f) = (g*h,max(e,f))",
    }


def _extension_receipt(branch_id: str, k: int, spectator: int) -> dict[str, Any]:
    coarse = build_capacity_packet(
        branch_id, k, spectator_multiplicity=spectator
    )
    fine = build_capacity_packet(
        branch_id, k + 1, spectator_multiplicity=spectator
    )
    coarse_graph = confusability_graph(coarse)
    fine_graph = confusability_graph(fine)
    old_records = set(coarse["records"])
    induced_exact = all(
        fine_graph[left] & old_records == coarse_graph[left]
        for left in old_records
    )
    return {
        "status": "PASS" if induced_exact else "FAIL",
        "from_k": k,
        "to_k": k + 1,
        "embedding": "record identifier inclusion",
        "old_graph_is_exact_induced_subgraph": induced_exact,
        "deficit_cannot_decrease": True,
    }


def _sewing_receipt(branch_id: str, k: int, spectator: int) -> dict[str, Any]:
    packet = build_capacity_packet(
        branch_id, k, spectator_multiplicity=spectator
    )
    output_coordinates = {
        tuple(output.split("|"))
        for output in packet["deterministic_collapse_outputs"].values()
    }
    collapsed_classes = sorted(
        {int(class_index) for _slot, class_index in output_coordinates}
    )
    left_sections = {
        (port, class_index)
        for port in PORTS
        for class_index in collapsed_classes
    }
    right_sections = {
        (orientation, class_index)
        for orientation in ORIENTATIONS
        for class_index in collapsed_classes
    }
    left_fibers = {
        class_index: {
            section
            for section in left_sections
            if section[1] == class_index
        }
        for class_index in collapsed_classes
    }
    right_fibers = {
        class_index: {
            section
            for section in right_sections
            if section[1] == class_index
        }
        for class_index in collapsed_classes
    }
    fiber_product_count = sum(
        len(left_fibers[class_index]) * len(right_fibers[class_index])
        for class_index in collapsed_classes
    )
    reconstructed_outputs = {
        (f"{port}/{orientation}", str(class_index))
        for class_index in collapsed_classes
        for port, _ in left_fibers[class_index]
        for orientation, _ in right_fibers[class_index]
    }
    bijection_verified = reconstructed_outputs == output_coordinates
    expected_public_output_count = expected_capacity(
        branch_id,
        k,
        spectator_multiplicity=spectator,
    )
    return {
        "status": (
            "PASS"
            if fiber_product_count == expected_public_output_count
            and bijection_verified
            else "FAIL"
        ),
        "seam_object": "collapsed-continuation-class",
        "left_patch": "twelve-port-by-collapsed-class",
        "right_patch": "two-orientation-by-collapsed-class",
        "source_port_count": len(PORTS),
        "source_orientation_count": len(ORIENTATIONS),
        "seam_label_count": len(collapsed_classes),
        "left_section_count": len(left_sections),
        "right_section_count": len(right_sections),
        "fiber_product_count": fiber_product_count,
        "expected_public_output_count": expected_public_output_count,
        "output_bijection_verified": bijection_verified,
        "formula": "sum_z |r_L^-1(z)| |r_R^-1(z)|",
    }


def build_certificate(projection: Mapping[str, Any] | None = None) -> dict[str, Any]:
    projection = dict(projection or build_projection())
    branch_receipts: dict[str, Any] = {}
    for branch in projection["branches"]:
        branch_id = branch["branch_id"]
        spectators = (
            SPECTATOR_MULTIPLICITIES
            if branch_id == "hidden_spectator"
            else (1,)
        )
        branch_receipts[branch_id] = {
            "all_rung_formula": branch["capacity_formula"],
            "exact_zero_set": {
                "reversible_identity": "all positive k",
                "copy_collapse_erasure": "{1}",
                "capped_two_class": "{1,2}",
                "hidden_spectator": (
                    "all positive k when s=1; empty when s>1"
                ),
            }[branch_id],
            "continuation_composition": [
                _continuation_receipt(branch_id, k, spectator)
                for k in SAMPLE_RUNGS
                for spectator in spectators
            ],
            "sewing": [
                _sewing_receipt(branch_id, k, spectator)
                for k in SAMPLE_RUNGS
                for spectator in spectators
            ],
            "extension": [
                _extension_receipt(branch_id, k, spectator)
                for k in SAMPLE_RUNGS[:-1]
                for spectator in spectators
            ],
        }
    distinct_bounded_zero_sets = {
        branch["branch_id"]: branch["claimed_bounded_zero_set"]
        for branch in projection["branches"]
    }
    all_receipts_pass = all(
        receipt["status"] == "PASS"
        for branch in branch_receipts.values()
        for family in (
            branch["continuation_composition"],
            branch["sewing"],
            branch["extension"],
        )
        for receipt in family
    )
    certificate = {
        "schema": CERTIFICATE_SCHEMA,
        "issue": 551,
        "status": (
            SCIENTIFIC_VERDICT if all_receipts_pass else "CERTIFICATE_FAILED"
        ),
        "claim_boundary": (
            "The declared base-agreement, positivity, and carrier-bound "
            "completion class admits continuations with different exact "
            "slack-zero sets. The result does not establish nonidentifiability "
            "for a complete A1-A3 lift of the fixed packet."
        ),
        "projection_sha256": tagged_sha256(canonical_json_bytes(projection)),
        "target_clean": not any(projection["target_cleanliness"].values()),
        "source_family_all_positive_rungs": True,
        "branch_grammar_complete_for_declared_counterfamily": True,
        "full_a1_a3_source_class_nonidentifiability_proved": False,
        "logical_completeness_reason": (
            "two models satisfying the common antecedent and disagreeing on "
            "the fixed set suffice to refute entailment; the four-row grammar "
            "also exposes identity, erasure, multiple-zero, and spectator "
            "directions"
        ),
        "branch_receipts": branch_receipts,
        "bounded_zero_sets": distinct_bounded_zero_sets,
        "zero_sets_differ": len(
            {
                canonical_json_bytes(rows)
                for rows in distinct_bounded_zero_sets.values()
            }
        )
        > 1,
        "physical_n_closure_promoted": False,
        "direct_n_status": (
            "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
        ),
        "remaining_positive_route": (
            "complete the A1-A3 lift of the fixed packet across the regulator "
            "family; if its admissible completions remain nonidentifiable, a "
            "separately named source law must select the continuation and "
            "identify both strange-loop readings"
        ),
    }
    certificate["certificate_sha256"] = tagged_sha256(
        canonical_json_bytes(certificate)
    )
    return certificate


def write_runtime() -> tuple[Path, Path]:
    projection = build_projection()
    certificate = build_certificate(projection)
    PROJECTION_PATH.write_bytes(canonical_json_bytes(projection))
    CERTIFICATE_PATH.write_bytes(canonical_json_bytes(certificate))
    return PROJECTION_PATH, CERTIFICATE_PATH


def verify_runtime() -> None:
    expected_projection = canonical_json_bytes(build_projection())
    expected_certificate = canonical_json_bytes(
        build_certificate(json.loads(expected_projection))
    )
    if PROJECTION_PATH.read_bytes() != expected_projection:
        raise SystemExit("capacity-indexed projection is stale")
    if CERTIFICATE_PATH.read_bytes() != expected_certificate:
        raise SystemExit("capacity-indexed certificate is stale")
    certificate = json.loads(expected_certificate)
    if certificate["status"] != SCIENTIFIC_VERDICT:
        raise SystemExit("capacity-indexed certificate did not attain its bounded verdict")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        write_runtime()
    if args.verify:
        verify_runtime()
    if not args.write and not args.verify:
        print(canonical_json_bytes(build_certificate()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
