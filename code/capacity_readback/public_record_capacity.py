#!/usr/bin/env python3
"""Superseded Pro4 checkpoint-fixed capacity control.

This evaluator reproduces the earlier shared-boundary-basis branch.  Fixed
projections are not the Pro5 capacity definition: a cyclic permutation can
have a trivial fixed algebra while preserving every label.  The canonical
evaluator is ``correctable_public_record_capacity.py``, which uses global
joint kernels and compound confusability graphs.

The bundled examples are mathematical controls with no physical promotion
weight.  CL-7 remains open until a source-derived trial-universe family emits
the complete receipts specified in F_READBACK_SPEC.md.
"""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


HERE = Path(__file__).resolve().parent

FAIL_STATES = {
    "NO_TERMINAL_NORMAL_FORM",
    "NO_TRIAL_REACHABILITY",
    "INCOMPLETE_TERMINAL_FIBER",
    "NO_STABLE_OBSERVER_SECTOR",
    "INCOMPLETE_OBSERVER_REGISTRY",
    "INCOMPLETE_INTERFACE_MANIFEST",
    "NO_PUBLIC_RECORD_DESCENT",
    "NO_PUBLIC_REACHABILITY_DESCENT",
    "NO_PUBLIC_RECORD_CHANNEL",
    "NO_BOUNDARY_REPRESENTATION",
    "NO_REACHABLE_PUBLIC_RECORD",
    "AMBIGUOUS_CAPACITY_READBACK",
    "CIRCULAR_CAPACITY_DEFINITION",
    "TARGET_TAINTED",
    "EW_BRIDGE_USED_AS_CAPACITY_PRODUCER",
    "RHO_USED_AS_CAPACITY_PRODUCER",
    "CAPACITY_EXTENSION_NATURALITY_FAILED",
    "NON_NATURAL_REFINEMENT",
    "NOT_MONOTONE",
    "NOT_DEFLATIONARY",
}


class UnionFind:
    def __init__(self, items: Sequence[tuple[str, str]]) -> None:
        self.parent = {item: item for item in items}

    def find(self, item: tuple[str, str]) -> tuple[str, str]:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: tuple[str, str], right: tuple[str, str]) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            self.parent[right_root] = left_root
        else:
            self.parent[left_root] = right_root


def _failure(status: str, **extra: Any) -> dict[str, Any]:
    if status not in FAIL_STATES:
        raise ValueError(f"unknown fail-closed status: {status}")
    return {"status": status, **extra}


def _validate_boundary_dimension(boundary_dimension: int) -> None:
    if not isinstance(boundary_dimension, int) or boundary_dimension < 1:
        raise ValueError("boundary_dimension must be a positive integer")


def _atom_table(
    boundary_dimension: int, terminal: Mapping[str, Any]
) -> dict[tuple[str, str], dict[str, Any]]:
    table: dict[tuple[str, str], dict[str, Any]] = {}
    observers = terminal.get("observers", {})
    for observer_id, observer in observers.items():
        atoms = observer.get("atoms", [])
        for atom in atoms:
            atom_id = atom["id"]
            key = (observer_id, atom_id)
            if key in table:
                raise ValueError(f"duplicate local record atom {key}")
            support = atom["boundary_support"]
            if (
                not isinstance(support, list)
                or not support
                or any(not isinstance(value, int) for value in support)
                or len(set(support)) != len(support)
                or any(not 0 <= value < boundary_dimension for value in support)
            ):
                raise ValueError(
                    f"boundary support for {key} must be a nonempty basis subset"
                )
            if "boundary_rank" in atom and atom["boundary_rank"] != len(support):
                raise ValueError(f"boundary rank disagrees with support for {key}")
            if not isinstance(atom["reachable"], bool):
                raise ValueError(f"reachable flag for {key} must be boolean")
            table[key] = dict(atom)
    return table


def evaluate_terminal(
    boundary_dimension: int,
    terminal: Mapping[str, Any],
) -> dict[str, Any]:
    """Compute the exact stable public code on one terminal quotient.

    Local commutative record algebras are supplied by their atoms.  Each
    interface row equates the restrictions of one atom from each endpoint.
    The resulting equivalence classes are the atoms of the finite public
    equalizer on this shared-boundary-basis branch.
    """
    _validate_boundary_dimension(boundary_dimension)
    terminal_id = terminal.get("terminal_id", "UNNAMED_TERMINAL")

    if not terminal.get("normal_form_certified", False):
        return _failure("NO_TERMINAL_NORMAL_FORM", terminal_id=terminal_id)
    if not terminal.get("reachable_from_trial", False):
        return _failure("NO_TRIAL_REACHABILITY", terminal_id=terminal_id)
    if terminal.get("self_read_predicate_injected", False):
        return _failure("CIRCULAR_CAPACITY_DEFINITION", terminal_id=terminal_id)
    if terminal.get("ew_bridge_target_used", False):
        return _failure(
            "EW_BRIDGE_USED_AS_CAPACITY_PRODUCER", terminal_id=terminal_id
        )
    if terminal.get("rho_used_as_capacity_producer", False):
        return _failure("RHO_USED_AS_CAPACITY_PRODUCER", terminal_id=terminal_id)
    if terminal.get("target_metadata_read_by_producer", False):
        return _failure("TARGET_TAINTED", terminal_id=terminal_id)
    if not terminal.get("observer_registry_complete", False):
        return _failure("INCOMPLETE_OBSERVER_REGISTRY", terminal_id=terminal_id)
    if not terminal.get("interface_manifest_complete", False):
        return _failure("INCOMPLETE_INTERFACE_MANIFEST", terminal_id=terminal_id)
    if not terminal.get("public_instrument_faithful", False) or not terminal.get(
        "public_reread_zero_error", False
    ):
        return _failure("NO_PUBLIC_RECORD_CHANNEL", terminal_id=terminal_id)

    observers = terminal.get("observers", {})
    if not observers:
        return _failure("NO_STABLE_OBSERVER_SECTOR", terminal_id=terminal_id)
    atom_table = _atom_table(boundary_dimension, terminal)
    if not atom_table:
        return _failure("NO_STABLE_OBSERVER_SECTOR", terminal_id=terminal_id)
    for observer_id in observers:
        observer_support_entries = [
            value
            for (candidate_id, _), atom in atom_table.items()
            if candidate_id == observer_id
            for value in atom["boundary_support"]
        ]
        if (
            len(observer_support_entries) != len(set(observer_support_entries))
            or set(observer_support_entries) != set(range(boundary_dimension))
        ):
            return _failure(
                "NO_BOUNDARY_REPRESENTATION",
                terminal_id=terminal_id,
                reason=(
                    f"local record atoms for {observer_id} do not resolve the "
                    "boundary identity"
                ),
            )

    equalizer = UnionFind(sorted(atom_table))
    for interface in terminal.get("interfaces", []):
        left_observer = interface["left_observer"]
        right_observer = interface["right_observer"]
        for row in interface.get("atom_pairs", []):
            left = (left_observer, row["left_atom"])
            right = (right_observer, row["right_atom"])
            if left not in atom_table or right not in atom_table:
                raise ValueError(f"interface references an unknown atom: {left}, {right}")
            equalizer.union(left, right)

    members_by_root: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for key in sorted(atom_table):
        members_by_root.setdefault(equalizer.find(key), []).append(key)
    ordered_members = sorted(tuple(sorted(members)) for members in members_by_root.values())
    public_id_by_local: dict[tuple[str, str], str] = {}
    for index, members in enumerate(ordered_members):
        public_id = f"p{index}"
        for member in members:
            public_id_by_local[member] = public_id

    public_rows: dict[str, dict[str, Any]] = {}
    for index, members in enumerate(ordered_members):
        public_id = f"p{index}"
        supports = {
            tuple(sorted(atom_table[member]["boundary_support"])) for member in members
        }
        if len(supports) != 1:
            return _failure(
                "NO_BOUNDARY_REPRESENTATION",
                terminal_id=terminal_id,
                reason=f"support mismatch in {public_id}",
            )
        reachable_flags = {atom_table[member]["reachable"] for member in members}
        if len(reachable_flags) != 1:
            return _failure(
                "NO_PUBLIC_REACHABILITY_DESCENT",
                terminal_id=terminal_id,
                reason=f"reachability mismatch in {public_id}",
            )
        image_ids: set[str] = set()
        for observer_id, atom_id in members:
            checkpoint_next = atom_table[(observer_id, atom_id)]["checkpoint_next"]
            target = (observer_id, checkpoint_next)
            if target not in atom_table:
                raise ValueError(f"checkpoint references an unknown atom: {target}")
            image_ids.add(public_id_by_local[target])
        if len(image_ids) != 1:
            return _failure(
                "NO_PUBLIC_RECORD_DESCENT",
                terminal_id=terminal_id,
                reason=f"checkpoint does not descend on {public_id}",
            )
        support = supports.pop()
        public_rows[public_id] = {
            "members": [f"{observer}:{atom}" for observer, atom in members],
            "boundary_support": list(support),
            "boundary_rank": len(support),
            "reachable": reachable_flags.pop(),
            "checkpoint_next": image_ids.pop(),
        }

    all_support_entries = [
        value
        for row in public_rows.values()
        for value in row["boundary_support"]
    ]
    if (
        len(all_support_entries) != len(set(all_support_entries))
        or set(all_support_entries) != set(range(boundary_dimension))
    ):
        return _failure(
            "NO_BOUNDARY_REPRESENTATION",
            terminal_id=terminal_id,
            reason=(
                "public atom supports are not pairwise orthogonal or do not resolve "
                "the supplied boundary identity"
            ),
        )

    checkpoint_components = UnionFind(
        [(public_id, "") for public_id in sorted(public_rows)]
    )
    for public_id, row in public_rows.items():
        checkpoint_components.union(
            (public_id, ""), (row["checkpoint_next"], "")
        )
    component_members: dict[tuple[str, str], list[str]] = {}
    for public_id in sorted(public_rows):
        root = checkpoint_components.find((public_id, ""))
        component_members.setdefault(root, []).append(public_id)

    stable_code: list[dict[str, Any]] = []
    for component in sorted(tuple(sorted(x)) for x in component_members.values()):
        component_set = set(component)
        preimage = {
            public_id
            for public_id, row in public_rows.items()
            if row["checkpoint_next"] in component_set
        }
        assert preimage == component_set
        reachable = any(public_rows[public_id]["reachable"] for public_id in component)
        if reachable:
            stable_code.append(
                {
                    "projection_atoms": list(component),
                    "boundary_rank": sum(
                        public_rows[public_id]["boundary_rank"]
                        for public_id in component
                    ),
                    "checkpoint_fixed": True,
                    "reachable": True,
                }
            )

    stable_code_size = len(stable_code)
    if stable_code_size == 0:
        return _failure("NO_REACHABLE_PUBLIC_RECORD", terminal_id=terminal_id)
    assert stable_code_size <= len(public_rows) <= boundary_dimension

    saturated = stable_code_size == boundary_dimension
    saturation_rank_one = saturated and all(
        row["boundary_rank"] == 1 for row in stable_code
    )
    if saturated:
        assert saturation_rank_one
        assert sum(row["boundary_rank"] for row in stable_code) == boundary_dimension

    slack = "0" if saturated else f"log({boundary_dimension}/{stable_code_size})"
    return {
        "status": "PASS",
        "terminal_id": terminal_id,
        "normal_form_hash": terminal.get("normal_form_hash"),
        "public_equalizer_atom_count": len(public_rows),
        "public_equalizer_atoms": public_rows,
        "checkpoint_descent_passed": True,
        "stable_public_code": stable_code,
        "stable_public_code_size": stable_code_size,
        "stable_public_log_capacity": f"log({stable_code_size})",
        "capacity_slack": slack,
        "record_capacity_le_boundary_capacity_passed": True,
        "saturated": saturated,
        "saturation_rank_one_complete_code": saturation_rank_one,
        "nontrivial_record_capacity": stable_code_size > 1,
    }


def evaluate_terminal_fiber(
    boundary_dimension: int,
    terminals: Sequence[Mapping[str, Any]],
    *,
    fiber_manifest_complete: bool,
) -> dict[str, Any]:
    """Return the set-valued row and scalarize only on complete-fiber agreement."""
    _validate_boundary_dimension(boundary_dimension)
    if not fiber_manifest_complete:
        return _failure("INCOMPLETE_TERMINAL_FIBER")
    if not terminals:
        return _failure("NO_TERMINAL_NORMAL_FORM")
    terminal_ids = [terminal.get("terminal_id") for terminal in terminals]
    if None in terminal_ids or len(set(terminal_ids)) != len(terminal_ids):
        raise ValueError("a complete terminal fiber needs unique terminal identifiers")

    terminal_results = [
        evaluate_terminal(boundary_dimension, terminal) for terminal in terminals
    ]
    failed = [row for row in terminal_results if row["status"] != "PASS"]
    if failed:
        return {
            "status": failed[0]["status"],
            "boundary_dimension": boundary_dimension,
            "terminal_results": terminal_results,
        }

    capacities = [row["stable_public_code_size"] for row in terminal_results]
    kernel = Counter(capacities)
    support = sorted(kernel)
    base = {
        "boundary_dimension": boundary_dimension,
        "terminal_fiber_size": len(terminals),
        "terminal_fiber_complete": True,
        "terminal_results": terminal_results,
        "set_valued_readback": support,
        "count_kernel_row": {str(value): kernel[value] for value in support},
        "terminal_fiber_capacity_diameter": support[-1] - support[0],
        "existential_self_read_closure": boundary_dimension in support,
    }
    if len(support) != 1:
        return {"status": "AMBIGUOUS_CAPACITY_READBACK", **base}

    readback = support[0]
    return {
        "status": "PASS",
        **base,
        "scalar_readback_dimension": readback,
        "scalar_readback_entropy": f"log({readback})",
        "deterministic_self_read_closure": readback == boundary_dimension,
        "nontriviality_status": (
            "NONTRIVIAL_CAPACITY_WITNESSED"
            if readback > 1
            else "NONTRIVIALITY_NOT_ESTABLISHED"
        ),
    }


def greatest_fixed_point(capacity_map: Mapping[int, int]) -> dict[str, Any]:
    """Certify the greatest fixed point of a finite monotone deflationary map."""
    if not capacity_map:
        raise ValueError("capacity_map must be nonempty")
    if any(not isinstance(key, int) for key in capacity_map) or any(
        not isinstance(value, int) for value in capacity_map.values()
    ):
        raise ValueError("capacity_map keys and values must be integers")
    chain = sorted(capacity_map)
    if chain != list(range(chain[-1] + 1)):
        raise ValueError("capacity_map domain must be the finite chain 0..D_max")
    if any(value not in capacity_map for value in capacity_map.values()):
        raise ValueError("capacity_map must be a self-map of its declared chain")
    for dimension in chain:
        if capacity_map[dimension] > dimension:
            return _failure(
                "NOT_DEFLATIONARY", witness=[dimension, capacity_map[dimension]]
            )
    for left, right in zip(chain, chain[1:], strict=False):
        if capacity_map[left] > capacity_map[right]:
            return _failure(
                "NOT_MONOTONE",
                witness=[left, right, capacity_map[left], capacity_map[right]],
            )

    path = [chain[-1]]
    while capacity_map[path[-1]] != path[-1]:
        path.append(capacity_map[path[-1]])
        assert len(path) <= len(chain)
    fixed_points = [dimension for dimension in chain if capacity_map[dimension] == dimension]
    greatest = path[-1]
    assert greatest == max(fixed_points)
    return {
        "status": "PASS",
        "monotone": True,
        "deflationary": True,
        "iteration_from_top": path,
        "fixed_points": fixed_points,
        "greatest_fixed_point": greatest,
        "fixed_point_unique": len(fixed_points) == 1,
        "positive_nontrivial_fixed_point": greatest > 1,
        "selection_rule": "greatest self-sustaining public-record capacity",
    }


def capacity_extension_naturality(
    capacity_table: Mapping[int, int],
    adjacent_embedding_certificates: Mapping[tuple[int, int], bool],
) -> dict[str, Any]:
    """Check the adjacent embedding receipts that imply capacity monotonicity."""
    if not capacity_table:
        raise ValueError("capacity_table must be nonempty")
    dimensions = sorted(capacity_table)
    if dimensions != list(range(1, dimensions[-1] + 1)):
        raise ValueError("capacity_table domain must be 1..D_max")
    if any(
        not isinstance(value, int) or not 1 <= value <= dimension
        for dimension, value in capacity_table.items()
    ):
        raise ValueError("capacity table values must be integers in 1..D")
    required = list(zip(dimensions, dimensions[1:], strict=False))
    missing = [
        pair for pair in required if not adjacent_embedding_certificates.get(pair, False)
    ]
    if missing:
        return _failure(
            "CAPACITY_EXTENSION_NATURALITY_FAILED",
            reason="missing adjacent stable-code embedding certificate",
            missing_pairs=[list(pair) for pair in missing],
        )
    for left, right in required:
        if capacity_table[left] > capacity_table[right]:
            return _failure(
                "CAPACITY_EXTENSION_NATURALITY_FAILED",
                reason="certified extension does not inject the stable code",
                witness=[left, right, capacity_table[left], capacity_table[right]],
            )
    return {
        "status": "PASS",
        "adjacent_embeddings_complete": True,
        "capacity_monotone": True,
        "certified_pairs": [list(pair) for pair in required],
    }


def refinement_stabilization(
    boundary_dimension: int,
    code_sizes: Sequence[int],
    *,
    embeddings_certified: bool,
    refinement_exhausted: bool = False,
) -> dict[str, Any]:
    """Check the exact monotone bounded-integer refinement theorem."""
    _validate_boundary_dimension(boundary_dimension)
    if not code_sizes:
        raise ValueError("at least one refinement code size is required")
    if not embeddings_certified:
        return _failure("NON_NATURAL_REFINEMENT")
    if any(
        not isinstance(value, int) or not 1 <= value <= boundary_dimension
        for value in code_sizes
    ):
        raise ValueError("refinement code sizes must be integers in 1..D")
    if any(
        left > right
        for left, right in zip(code_sizes, code_sizes[1:], strict=False)
    ):
        return _failure("NON_NATURAL_REFINEMENT", witness=list(code_sizes))

    current = code_sizes[-1]
    suffix_length = 1
    for value in reversed(code_sizes[:-1]):
        if value != current:
            break
        suffix_length += 1
    saturated = current == boundary_dimension
    return {
        "status": "PASS",
        "embedding_naturality_passed": True,
        "monotone_nondecreasing": True,
        "bounded_by_boundary_dimension": True,
        "eventual_stabilization_theorem_applies": True,
        "strict_increases_remaining_upper_bound": boundary_dimension - current,
        "observed_constant_suffix_length": suffix_length,
        "permanent_stabilization_certified": saturated or refinement_exhausted,
        "permanent_stabilization_reason": (
            "capacity bound saturated"
            if saturated
            else "refinement family exhausted"
            if refinement_exhausted
            else "eventual only; no finite stage claim"
        ),
    }


def boundary_basis_trial(
    boundary_dimension: int,
    *,
    checkpoint: str,
    terminal_id: str,
    normal_form_hash: str = "same-terminal-normal-form",
) -> dict[str, Any]:
    """Build a two-observer shared-basis control diagram."""
    _validate_boundary_dimension(boundary_dimension)
    if checkpoint not in {"identity", "erase_to_first"}:
        raise ValueError("checkpoint must be identity or erase_to_first")

    def atoms() -> list[dict[str, Any]]:
        return [
            {
                "id": f"b{index}",
                "boundary_support": [index],
                "boundary_rank": 1,
                "reachable": True,
                "checkpoint_next": (
                    f"b{index}" if checkpoint == "identity" else "b0"
                ),
            }
            for index in range(boundary_dimension)
        ]

    return {
        "terminal_id": terminal_id,
        "normal_form_certified": True,
        "normal_form_hash": normal_form_hash,
        "reachable_from_trial": True,
        "observer_registry_complete": True,
        "interface_manifest_complete": True,
        "public_instrument_faithful": True,
        "public_reread_zero_error": True,
        "self_read_predicate_injected": False,
        "target_metadata_read_by_producer": False,
        "ew_bridge_target_used": False,
        "rho_used_as_capacity_producer": False,
        "observers": {
            "alice": {"atoms": atoms()},
            "bob": {"atoms": atoms()},
        },
        "interfaces": [
            {
                "left_observer": "alice",
                "right_observer": "bob",
                "atom_pairs": [
                    {"left_atom": f"b{index}", "right_atom": f"b{index}"}
                    for index in range(boundary_dimension)
                ],
            }
        ],
    }


def example_receipt() -> dict[str, Any]:
    """Emit exact nonphysical controls for the direct capacity contract."""
    dimension = 4
    identity_a = boundary_basis_trial(
        dimension, checkpoint="identity", terminal_id="identity-a"
    )
    identity_b = boundary_basis_trial(
        dimension, checkpoint="identity", terminal_id="identity-b"
    )
    erasure = boundary_basis_trial(
        dimension, checkpoint="erase_to_first", terminal_id="erasure"
    )
    identity_fiber = evaluate_terminal_fiber(
        dimension, [identity_a, identity_b], fiber_manifest_complete=True
    )
    erasure_result = evaluate_terminal(dimension, erasure)
    ambiguous_fiber = evaluate_terminal_fiber(
        dimension, [identity_a, erasure], fiber_manifest_complete=True
    )
    incomplete_fiber = evaluate_terminal_fiber(
        dimension, [identity_a], fiber_manifest_complete=False
    )
    identity_gfp = greatest_fixed_point({value: value for value in range(5)})
    constant_one_gfp = greatest_fixed_point(
        {0: 0, **{value: 1 for value in range(1, 5)}}
    )
    refinement = refinement_stabilization(
        dimension, [1, 2, 4], embeddings_certified=True
    )
    extension = capacity_extension_naturality(
        {value: value for value in range(1, 5)},
        {(value, value + 1): True for value in range(1, 4)},
    )

    assert identity_fiber["scalar_readback_dimension"] == dimension
    assert erasure_result["stable_public_code_size"] == 1
    assert ambiguous_fiber["status"] == "AMBIGUOUS_CAPACITY_READBACK"
    assert incomplete_fiber["status"] == "INCOMPLETE_TERMINAL_FIBER"
    return {
        "artifact": "public_stable_record_capacity_example_v1",
        "status": "SCHEMA_ONLY",
        "physical_content": False,
        "moves_cl7": False,
        "cl7_status": "open",
        "canonical_F_producer": "stable public record code size",
        "supplied_dimension": dimension,
        "supplied_dimension_factorization": [2, 2],
        "boundary_coordinate": "D = dim(H_boundary)",
        "entropy_coordinate": "F_r(log D) = log(M_pub(D))",
        "identity_boundary_basis_control": identity_fiber,
        "same_normal_form_erasure_control": erasure_result,
        "ambiguous_terminal_fiber_control": ambiguous_fiber,
        "incomplete_terminal_fiber_control": incomplete_fiber,
        "identity_map_greatest_fixed_point_control": identity_gfp,
        "constant_one_greatest_fixed_point_control": constant_one_gfp,
        "capacity_extension_naturality_control": extension,
        "refinement_stabilization_control": refinement,
        "rho_op_role": "independent estimator only; never an input to F",
        "self_read_predicate_injected": False,
        "supplied_capacity_metadata_read_by_producer": False,
        "measured_lambda_used": False,
        "ew_bridge_target_used": False,
        "rho_used_as_capacity_producer": False,
        "fail_states": sorted(FAIL_STATES),
        "open_direct_capacity_gates": [
            "source-derived complete trial-universe and terminal-fiber family",
            "complete public observer, overlap, reachability, and checkpoint receipts",
            "capacity-extension embeddings and naturality",
            "refinement record embeddings",
            "nontrivial fixed-point or greatest-capacity selection receipt",
            "large-screen capacity-density control",
        ],
        "postconstruction_physical_gates": [
            "HORIZON-RECORD-SATURATION",
            "COMMON-EW-LOAD-CARRIER",
            "independent operational rho_op experiment",
        ],
    }


if __name__ == "__main__":
    payload = example_receipt()
    out = HERE / "runtime" / "public_record_capacity_example.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
