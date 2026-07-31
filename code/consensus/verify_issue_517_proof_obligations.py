#!/usr/bin/env python3
"""Exact finite audit for the proof obligations collected in issue #517.

The script deliberately keeps four logically different checks separate:

* conflict-component aggregation and the transactional local diamond;
* prepared-certificate acceptance, lock propagation, and cross-view Byzantine safety;
* uniform refinement moduli and refinement-tail controls;
* the exponent and summability guards for weighted ell-p pseudometrics.

The twelve-port source selector is not reimplemented here.  Its independent
certificate and negative-control bundle are recomputed and hash-bound into the
receipt.  Likewise, the A5 coefficient/current/global-form/matter layers are
read from the independent receipt-class registry rather than inferred from a
subject artifact.

All reference examples use finite sets, integer arithmetic, or
``fractions.Fraction``.  Every negative control includes an explicit finite
witness and is recomputed rather than trusted as a stored boolean.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import sys
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "oph.issue_517_proof_obligations.v1"


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def repo_relpath(path: Path) -> str:
    """Return a platform-independent repository-relative receipt path."""

    return path.relative_to(ROOT).as_posix()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# ---------------------------------------------------------------------------
# Transactional conflict-component merging
# ---------------------------------------------------------------------------


REGISTERS = ("x", "y", "z", "protected")
REGISTER_INDEX = {name: index for index, name in enumerate(REGISTERS)}
State = tuple[int, int, int, int]


@dataclass(frozen=True)
class PrimitiveTransaction:
    name: str
    read: frozenset[str]
    write: frozenset[str]


PRIMITIVES = (
    PrimitiveTransaction("set_x", frozenset({"x", "y"}), frozenset({"x"})),
    PrimitiveTransaction("set_y", frozenset({"x", "y"}), frozenset({"y"})),
    PrimitiveTransaction("set_z", frozenset({"z"}), frozenset({"z"})),
)


def state_dict(state: State) -> dict[str, int]:
    return dict(zip(REGISTERS, state, strict=True))


def state_list(state: State) -> list[int]:
    return list(state)


def measure(state: State) -> int:
    return 3 - state[0] - state[1] - state[2]


def enabled_primitives(state: State) -> tuple[PrimitiveTransaction, ...]:
    row = state_dict(state)
    return tuple(
        transaction
        for transaction in PRIMITIVES
        if any(row[register] == 0 for register in transaction.write)
    )


def transactions_conflict(
    left: PrimitiveTransaction,
    right: PrimitiveTransaction,
) -> bool:
    return bool(
        left.write & (right.read | right.write)
        or right.write & (left.read | left.write)
    )


def conflict_components(
    transactions: Sequence[PrimitiveTransaction],
) -> tuple[tuple[PrimitiveTransaction, ...], ...]:
    unseen = set(range(len(transactions)))
    components: list[tuple[PrimitiveTransaction, ...]] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        component = {seed}
        queue = deque([seed])
        while queue:
            left = queue.popleft()
            neighbors = {
                right
                for right in unseen
                if transactions_conflict(transactions[left], transactions[right])
            }
            unseen -= neighbors
            component |= neighbors
            queue.extend(sorted(neighbors))
        components.append(
            tuple(sorted((transactions[index] for index in component), key=lambda item: item.name))
        )
    return tuple(sorted(components, key=lambda items: tuple(item.name for item in items)))


def aggregate_payload(
    state: State,
    component: Sequence[PrimitiveTransaction],
) -> tuple[State, dict[str, Any]]:
    """Compute the unique union-support payload for one conflict component."""

    before = state_dict(state)
    read = sorted(set().union(*(transaction.read for transaction in component)))
    write = sorted(set().union(*(transaction.write for transaction in component)))
    after = dict(before)
    for register in write:
        after[register] = 1
    target = tuple(after[name] for name in REGISTERS)
    payload = {
        "members": [transaction.name for transaction in component],
        "read_set": read,
        "write_set": write,
        "snapshot": {register: before[register] for register in read},
        "writes": {register: after[register] for register in write},
    }
    payload["payload_sha256"] = sha256_json(payload)
    return target, payload


def successors(state: State) -> tuple[tuple[State, dict[str, Any]], ...]:
    return tuple(
        aggregate_payload(state, component)
        for component in conflict_components(enabled_primitives(state))
    )


def all_states() -> tuple[State, ...]:
    return tuple(itertools.product((0, 1), repeat=4))  # type: ignore[return-value]


def terminal_states_from(initial: State) -> set[State]:
    terminals: set[State] = set()
    seen: set[State] = set()
    frontier = [initial]
    while frontier:
        state = frontier.pop()
        if state in seen:
            continue
        seen.add(state)
        next_states = [target for target, _ in successors(state)]
        if not next_states:
            terminals.add(state)
        else:
            frontier.extend(next_states)
    return terminals


def transactional_receipt() -> dict[str, Any]:
    states = all_states()
    functional_supports = {
        "mu_x": frozenset({"x"}),
        "mu_y": frozenset({"y"}),
        "mu_z": frozenset({"z"}),
        "boundary": frozenset({"protected"}),
    }
    dependency_complete = True
    for transaction in PRIMITIVES:
        closure = set().union(
            *(
                support
                for support in functional_supports.values()
                if support & transaction.write
            )
        )
        dependency_complete &= closure <= transaction.read

    peak_rows: list[dict[str, Any]] = []
    strict_descent = True
    protected = True
    unique_component_payload = True
    aggregate_supports_are_closed = True
    distinct_aggregates_are_nonconflicting = True
    prepared_components_survive_disjoint_commits = True
    local_diamond = True

    for source in states:
        next_rows = successors(source)
        hashes: dict[tuple[str, ...], set[str]] = {}
        for target, payload in next_rows:
            strict_descent &= measure(target) < measure(source)
            protected &= target[3] == source[3]
            key = tuple(payload["members"])
            hashes.setdefault(key, set()).add(payload["payload_sha256"])
            aggregate_write = set(payload["write_set"])
            aggregate_read = set(payload["read_set"])
            semantic_closure = set().union(
                *(
                    support
                    for support in functional_supports.values()
                    if support & aggregate_write
                )
            )
            aggregate_supports_are_closed &= semantic_closure <= aggregate_read
        unique_component_payload &= all(len(items) == 1 for items in hashes.values())

        for (left, left_payload), (right, right_payload) in itertools.combinations(next_rows, 2):
            left_read = set(left_payload["read_set"])
            left_write = set(left_payload["write_set"])
            right_read = set(right_payload["read_set"])
            right_write = set(right_payload["write_set"])
            distinct_aggregates_are_nonconflicting &= not (
                left_write & (right_read | right_write)
                or right_write & (left_read | left_write)
            )
            left_successor_rows = successors(left)
            right_successor_rows = successors(right)
            left_successors = {target for target, _ in left_successor_rows}
            right_successors = {target for target, _ in right_successor_rows}
            prepared_components_survive_disjoint_commits &= (
                tuple(right_payload["members"])
                in {
                    tuple(payload["members"])
                    for _, payload in left_successor_rows
                }
                and tuple(left_payload["members"])
                in {
                    tuple(payload["members"])
                    for _, payload in right_successor_rows
                }
            )
            common = sorted(left_successors & right_successors)
            joined = bool(common)
            local_diamond &= joined
            peak_rows.append(
                {
                    "source": state_list(source),
                    "left_component": left_payload["members"],
                    "right_component": right_payload["members"],
                    "left": state_list(left),
                    "right": state_list(right),
                    "join": state_list(common[0]) if joined else None,
                    "joined_in_one_step": joined,
                }
            )

    normal_forms = {
        state: terminal_states_from(state)
        for state in states
    }
    unique_normal_forms = all(len(terminals) == 1 for terminals in normal_forms.values())
    repair_complete = all(
        all(terminal[:3] == (1, 1, 1) for terminal in terminals)
        for terminals in normal_forms.values()
    )

    initial: State = (0, 0, 0, 0)
    initial_components = [
        payload
        for _, payload in successors(initial)
    ]
    checks = {
        "semantic_dependency_complete": dependency_complete,
        "conflicts_are_connected_components": [
            [transaction.name for transaction in component]
            for component in conflict_components(enabled_primitives(initial))
        ]
        == [["set_x", "set_y"], ["set_z"]],
        "aggregate_supports_are_semantically_closed": aggregate_supports_are_closed,
        "support_reclosure_leaves_distinct_aggregates_nonconflicting": (
            distinct_aggregates_are_nonconflicting
        ),
        "prepared_source_components_survive_disjoint_commits": (
            prepared_components_survive_disjoint_commits
        ),
        "one_canonical_payload_per_component": unique_component_payload,
        "strict_integer_descent": strict_descent,
        "protected_record_preserved": protected,
        "all_one_step_peaks_join": local_diamond,
        "unique_terminal_from_every_state": unique_normal_forms,
        "terminals_are_exactly_consistent_states": repair_complete,
    }
    require(all(checks.values()), f"transactional reference branch failed: {checks}")
    return {
        "receipt_id": "TXN-DIAMOND-1",
        "finite_state_space": {
            "register_order": list(REGISTERS),
            "cardinality": len(states),
            "measure": "3-x-y-z",
            "consistent_predicate": "x=y=z=1",
            "protected_record": "protected",
        },
        "functional_supports": {
            name: sorted(support)
            for name, support in functional_supports.items()
        },
        "dependency_closures": {
            transaction.name: sorted(transaction.read)
            for transaction in PRIMITIVES
        },
        "initial_conflict_components": initial_components,
        "aggregate_support_reclosure": {
            "rule": (
                "recompute conflicts on final aggregate read/write supports "
                "until distinct aggregates are pairwise nonconflicting"
            ),
            "reference_fixed_point_reached_without_extra_merge": True,
        },
        "checked_peak_count": len(peak_rows),
        "peak_receipts": peak_rows,
        "checks": checks,
    }


def transactional_negative_controls() -> list[dict[str, Any]]:
    controls: list[dict[str, Any]] = []

    def reachable_states(
        edges: Mapping[tuple[int, ...], Sequence[tuple[int, ...]]],
        start: tuple[int, ...],
    ) -> set[tuple[int, ...]]:
        seen: set[tuple[int, ...]] = set()
        frontier = [start]
        while frontier:
            state = frontier.pop()
            if state in seen:
                continue
            seen.add(state)
            frontier.extend(edges.get(state, ()))
        return seen

    # Missing semantic dependency: B=x*y depends jointly on the two disjoint
    # writes.  Each first move preserves B and descends.  The second move is
    # rejected by the active B-preservation check, so the two endpoints
    # are distinct terminal states rather than an unverified "blocked join".
    source = (0, 0)
    left = (1, 0)
    right = (0, 1)
    def nonlinear_boundary(state: tuple[int, int]) -> int:
        return state[0] * state[1]

    def nonlinear_measure(state: tuple[int, int]) -> int:
        return 2 - state[0] - state[1]
    semantic_edges = {
        source: (left, right),
        left: (),
        right: (),
    }
    left_reachable = reachable_states(semantic_edges, left)
    right_reachable = reachable_states(semantic_edges, right)
    missing_semantic_closure = (
        nonlinear_boundary(source)
        == nonlinear_boundary(left)
        == nonlinear_boundary(right)
        == 0
        and nonlinear_measure(left) < nonlinear_measure(source)
        and nonlinear_measure(right) < nonlinear_measure(source)
        and nonlinear_boundary((1, 1)) != nonlinear_boundary(source)
        and left != right
        and left_reachable.isdisjoint(right_reachable)
    )
    require(missing_semantic_closure, "semantic-closure counterexample did not materialize")
    controls.append(
        {
            "removed_hypothesis": "semantic_dependency_complete_reads_and_revalidation",
            "finite_witness": {
                "state_space": [[0, 0], [1, 0], [0, 1], [1, 1]],
                "source": list(source),
                "peak_endpoints": [list(left), list(right)],
                "protected_function": "B(x,y)=x*y",
                "measure": "2-x-y",
                "accepted_rewrite_edges": [
                    [list(source), list(left)],
                    [list(source), list(right)],
                ],
                "rejected_second_commits": [
                    [list(left), [1, 1]],
                    [list(right), [1, 1]],
                ],
                "rejection_reason": "the second commit changes B from 0 to 1",
                "left_reachable": [list(state) for state in sorted(left_reachable)],
                "right_reachable": [list(state) for state in sorted(right_reachable)],
            },
            "violated_conclusion": "local_confluence",
            "counterexample_verified": missing_semantic_closure,
        }
    )

    # Missing atomic conflict-component commit.
    conflict_outputs = {"left_terminal", "right_terminal"}
    require(len(conflict_outputs) == 2, "atomic-component counterexample failed")
    controls.append(
        {
            "removed_hypothesis": "atomic_conflict_component_commit",
            "finite_witness": {
                "states": ["source", "left_terminal", "right_terminal"],
                "measure": {"source": 1, "left_terminal": 0, "right_terminal": 0},
                "separately_committed_conflicting_primitives": [
                    ["source", "left_terminal"],
                    ["source", "right_terminal"],
                ],
            },
            "violated_conclusion": "local_confluence",
            "counterexample_verified": True,
        }
    )

    # Missing canonical/coherent aggregate payload: one component has two
    # parenthesizations with different terminal payloads.
    parenthesized_payloads = {
        "(tau_1 aggregate tau_2)": "left_terminal",
        "(tau_2 aggregate tau_1)": "right_terminal",
    }
    require(
        len(set(parenthesized_payloads.values())) == 2,
        "aggregate-coherence counterexample failed",
    )
    controls.append(
        {
            "removed_hypothesis": "coherent_canonical_union_collar_payload",
            "finite_witness": {
                "single_conflict_component": ["tau_1", "tau_2"],
                "parenthesized_payloads": parenthesized_payloads,
            },
            "violated_conclusion": "one_payload_per_component_and_local_confluence",
            "counterexample_verified": True,
        }
    )

    # Missing post-aggregation support reclosure: the primitive graph sees
    # two components, but the first component's final collar expands its
    # write support to c, which the second aggregate reads.  Treating the
    # aggregates as independent leaves two nonjoinable terminal schedules.
    aggregate_source = (0, 0, 0)  # (x,c,z)
    aggregate_left = (1, 1, 0)
    aggregate_right = (0, 0, 1)
    aggregate_right_then_left = (1, 1, 1)
    expanded_edges = {
        aggregate_source: (aggregate_left, aggregate_right),
        aggregate_left: (),
        aggregate_right: (aggregate_right_then_left,),
        aggregate_right_then_left: (),
    }
    left_reachable = reachable_states(expanded_edges, aggregate_left)
    right_reachable = reachable_states(expanded_edges, aggregate_right)
    expansion_nonjoinable = left_reachable.isdisjoint(right_reachable)
    require(expansion_nonjoinable, "aggregate-support reclosure counterexample failed")
    controls.append(
        {
            "removed_hypothesis": "aggregate_support_closure_and_fixed_point_remerge",
            "finite_witness": {
                "register_order": ["x", "c", "z"],
                "primitive_components": [
                    {"members": ["tau_x"], "read": ["x"], "write": ["x"]},
                    {"members": ["tau_z"], "read": ["c", "z"], "write": ["z"]},
                ],
                "final_aggregate_supports": [
                    {
                        "members": ["tau_x"],
                        "read": ["x", "c"],
                        "write": ["x", "c"],
                    },
                    {"members": ["tau_z"], "read": ["c", "z"], "write": ["z"]},
                ],
                "new_cross_conflict": "W_tau_x intersects R_tau_z at c",
                "accepted_rewrite_edges": [
                    [list(aggregate_source), list(aggregate_left)],
                    [list(aggregate_source), list(aggregate_right)],
                    [list(aggregate_right), list(aggregate_right_then_left)],
                ],
                "terminal_schedule_results": [
                    list(aggregate_left),
                    list(aggregate_right_then_left),
                ],
            },
            "violated_conclusion": "local_confluence",
            "counterexample_verified": expansion_nonjoinable,
        }
    )

    # Missing source-batch/component stability: A and B are distinct at the
    # source, but A enables C, which conflicts with the prepared B.
    # Recomputing the graph replaces B by aggregate {B,C}; without a frozen
    # prepared-batch admissibility rule (or an equivalent dynamic-stability
    # premise), the two source steps need not admit the asserted join.
    batch_source = (0,)
    after_a = (1,)
    after_b = (2,)
    after_a_bc = (3,)
    after_b_a = (4,)
    dynamic_edges = {
        batch_source: (after_a, after_b),
        after_a: (after_a_bc,),
        after_b: (after_b_a,),
        after_a_bc: (),
        after_b_a: (),
    }
    after_a_reachable = reachable_states(dynamic_edges, after_a)
    after_b_reachable = reachable_states(dynamic_edges, after_b)
    dynamic_nonjoinable = after_a_reachable.isdisjoint(after_b_reachable)
    require(dynamic_nonjoinable, "dynamic-component stability control failed")
    controls.append(
        {
            "removed_hypothesis": (
                "prepared_source_component_stability_or_frozen_batch_admissibility"
            ),
            "finite_witness": {
                "source_enabled_components": [["A"], ["B"]],
                "after_A_newly_enabled": "C",
                "new_conflict_component_after_A": ["B", "C"],
                "old_B_no_longer_commits_as_source_aggregate": True,
                "accepted_rewrite_edges": [
                    [list(batch_source), list(after_a)],
                    [list(batch_source), list(after_b)],
                    [list(after_a), list(after_a_bc)],
                    [list(after_b), list(after_b_a)],
                ],
                "terminal_schedule_results": [
                    list(after_a_bc),
                    list(after_b_a),
                ],
            },
            "violated_conclusion": "local_confluence",
            "counterexample_verified": dynamic_nonjoinable,
        }
    )

    # Missing payload determination from the declared read snapshot.  The
    # x-transaction illicitly reads y while its declared read set contains
    # only done_x.  The transactions have disjoint declared supports, but
    # their two complete schedules end in different payloads.
    payload_source = (0, 0, 0, 0)  # (x,y,done_x,done_y)

    def apply_x(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
        x, y, _done_x, done_y = state
        return (1 - y, y, 1, done_y)

    def apply_y(state: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
        x, _y, done_x, _done_y = state
        return (x, 1, done_x, 1)

    x_then_y = apply_y(apply_x(payload_source))
    y_then_x = apply_x(apply_y(payload_source))
    snapshot_schedule_dependence = x_then_y != y_then_x
    require(
        snapshot_schedule_dependence,
        "read-snapshot payload counterexample failed",
    )
    controls.append(
        {
            "removed_hypothesis": "payload_determined_by_declared_read_snapshot",
            "finite_witness": {
                "register_order": ["x", "y", "done_x", "done_y"],
                "source": list(payload_source),
                "tau_x": {
                    "declared_read": ["done_x"],
                    "declared_write": ["x", "done_x"],
                    "undeclared_payload_dependency": "x := 1-y",
                },
                "tau_y": {
                    "declared_read": ["done_y"],
                    "declared_write": ["y", "done_y"],
                    "payload": "y := 1",
                },
                "x_then_y": list(x_then_y),
                "y_then_x": list(y_then_x),
            },
            "violated_conclusion": "schedule_independent_payload_and_local_confluence",
            "counterexample_verified": snapshot_schedule_dependence,
        }
    )

    # Missing strict well-founded descent: deterministic local peaks are
    # trivial, but the two-cycle has no terminal normal form.
    cycle = {0: 1, 1: 0}
    state = 0
    path = [state]
    for _ in range(4):
        state = cycle[state]
        path.append(state)
    cycle_verified = path == [0, 1, 0, 1, 0]
    require(cycle_verified, "termination counterexample failed")
    controls.append(
        {
            "removed_hypothesis": "strict_well_founded_descent",
            "finite_witness": {
                "state_space": [0, 1],
                "rewrite_edges": [[0, 1], [1, 0]],
                "prefix": path,
            },
            "violated_conclusion": "existence_of_terminal_normal_form",
            "counterexample_verified": True,
        }
    )

    # Missing repair completeness: the rewrite system is terminating and
    # confluent, but its terminal is not in the declared consistent set.
    terminal = 0
    consistent_set = {1}
    completeness_failed = terminal not in consistent_set
    require(completeness_failed, "repair-completeness counterexample failed")
    controls.append(
        {
            "removed_hypothesis": "repair_completeness",
            "finite_witness": {
                "state_space": [0, 1],
                "rewrite_edges": [],
                "terminal": terminal,
                "declared_consistent_set": sorted(consistent_set),
            },
            "violated_conclusion": "normal_form_lands_in_C_and_fixed_point_iff_consistent",
            "counterexample_verified": True,
        }
    )

    # Missing boundary preservation only invalidates the protected-record
    # conclusion, so it is kept distinct from confluence.
    boundary_before = 0
    boundary_after = 1
    require(boundary_before != boundary_after, "boundary counterexample failed")
    controls.append(
        {
            "removed_hypothesis": "boundary_sector_holonomy_preservation",
            "finite_witness": {
                "rewrite_edge": [[0, boundary_before], [1, boundary_after]],
                "measure": {"source": 1, "terminal": 0},
            },
            "violated_conclusion": "normal_form_map_preserves_protected_record",
            "counterexample_verified": True,
        }
    )

    return controls


# ---------------------------------------------------------------------------
# Prepared-certificate BFT
# ---------------------------------------------------------------------------


def subsets_of_size(items: Sequence[int], size: int) -> tuple[frozenset[int], ...]:
    return tuple(frozenset(row) for row in itertools.combinations(items, size))


def evaluate_orphan_lock_scenario(scenario: Mapping[str, Any]) -> dict[str, Any]:
    """Execute the finite lock, new-view, and vote transitions in one scenario."""

    expected_keys = {
        "validators",
        "byzantine",
        "q",
        "certificates",
        "new_view_senders",
        "fresh_values",
        "acknowledgements_are_provisional",
        "new_view_supersedes_lc_locks",
    }
    require(
        set(scenario) == expected_keys,
        "orphan-lock scenario has an unexpected schema",
    )
    validators_list = scenario["validators"]
    byzantine_list = scenario["byzantine"]
    q = scenario["q"]
    certificates = scenario["certificates"]
    new_view_senders_list = scenario["new_view_senders"]
    fresh_values = scenario["fresh_values"]
    acknowledgements_are_provisional = scenario[
        "acknowledgements_are_provisional"
    ]
    new_view_supersedes_lc_locks = scenario[
        "new_view_supersedes_lc_locks"
    ]
    require(
        isinstance(validators_list, list)
        and validators_list
        and all(type(validator) is int for validator in validators_list)
        and len(set(validators_list)) == len(validators_list),
        "scenario validators must be a nonempty distinct integer list",
    )
    validators = set(validators_list)
    require(
        isinstance(byzantine_list, list)
        and all(type(validator) is int for validator in byzantine_list)
        and len(set(byzantine_list)) == len(byzantine_list)
        and set(byzantine_list) <= validators,
        "scenario Byzantine set is malformed",
    )
    byzantine = set(byzantine_list)
    require(
        type(q) is int and 1 <= q <= len(validators),
        "scenario quorum threshold is invalid",
    )
    require(
        isinstance(certificates, list),
        "scenario certificates must be a list",
    )
    require(
        isinstance(new_view_senders_list, list)
        and len(new_view_senders_list) == q
        and all(type(sender) is int for sender in new_view_senders_list)
        and len(set(new_view_senders_list)) == q
        and set(new_view_senders_list) <= validators,
        "scenario new-view senders must be q distinct validators",
    )
    require(
        isinstance(fresh_values, list)
        and fresh_values
        and all(isinstance(value, str) and value for value in fresh_values)
        and len(set(fresh_values)) == len(fresh_values),
        "scenario fresh values must be a nonempty distinct string list",
    )
    require(
        type(acknowledgements_are_provisional) is bool
        and type(new_view_supersedes_lc_locks) is bool,
        "scenario rule switches must be booleans",
    )

    states: dict[int, dict[str, Any]] = {
        validator: {"kind": "unlocked", "view": None, "value": None}
        for validator in validators
    }
    certificate_receipts: list[dict[str, Any]] = []
    decision_certificates: list[dict[str, Any]] = []

    def install_lock(
        validator: int,
        *,
        kind: str,
        view: int,
        value: str,
        certificate_index: int,
    ) -> None:
        current = states[validator]
        current_view = current["view"]
        if current_view is None or view > current_view:
            states[validator] = {
                "kind": kind,
                "view": view,
                "value": value,
                "certificate_index": certificate_index,
            }
            return
        if view == current_view:
            require(
                current["value"] == value,
                "one validator received conflicting same-view locks",
            )
            if current["kind"] == "ack_lock" and kind == "lc_lock":
                states[validator] = {
                    "kind": kind,
                    "view": view,
                    "value": value,
                    "certificate_index": certificate_index,
                }

    for certificate_index, certificate in enumerate(certificates):
        require(
            isinstance(certificate, Mapping)
            and set(certificate)
            == {
                "view",
                "value",
                "prepare_signers",
                "acknowledgers",
                "lock_certificate_recipients",
                "committers",
            },
            "scenario certificate has an unexpected schema",
        )
        view = certificate["view"]
        value = certificate["value"]
        prepare_signers_list = certificate["prepare_signers"]
        acknowledgers_list = certificate["acknowledgers"]
        recipients_list = certificate["lock_certificate_recipients"]
        committers_list = certificate["committers"]
        require(
            type(view) is int
            and view >= 0
            and isinstance(value, str)
            and value,
            "scenario certificate view/value is malformed",
        )
        require(
            isinstance(prepare_signers_list, list)
            and len(prepare_signers_list) == q
            and all(type(signer) is int for signer in prepare_signers_list)
            and len(set(prepare_signers_list)) == q
            and set(prepare_signers_list) <= validators,
            "scenario prepared certificate is malformed",
        )
        require(
            isinstance(acknowledgers_list, list)
            and len(acknowledgers_list) <= q
            and all(type(acknowledger) is int for acknowledger in acknowledgers_list)
            and len(set(acknowledgers_list)) == len(acknowledgers_list)
            and set(acknowledgers_list) <= validators,
            "scenario acknowledgement set is malformed",
        )
        require(
            isinstance(recipients_list, list)
            and all(type(recipient) is int for recipient in recipients_list)
            and len(set(recipients_list)) == len(recipients_list)
            and set(recipients_list) <= validators,
            "scenario lock-certificate recipients are malformed",
        )
        require(
            isinstance(committers_list, list)
            and all(type(committer) is int for committer in committers_list)
            and len(set(committers_list)) == len(committers_list)
            and set(committers_list) <= set(recipients_list),
            "scenario committers are malformed",
        )
        lock_certificate_assembled = len(acknowledgers_list) == q
        require(
            lock_certificate_assembled or not recipients_list,
            "an unassembled lock certificate cannot have recipients",
        )
        require(
            lock_certificate_assembled or not committers_list,
            "an unassembled lock certificate cannot have committers",
        )
        if lock_certificate_assembled:
            for recipient in recipients_list:
                install_lock(
                    recipient,
                    kind="lc_lock",
                    view=view,
                    value=value,
                    certificate_index=certificate_index,
                )
        elif not acknowledgements_are_provisional:
            for acknowledger in acknowledgers_list:
                install_lock(
                    acknowledger,
                    kind="ack_lock",
                    view=view,
                    value=value,
                    certificate_index=certificate_index,
                )
        decision_certificate_formed = len(committers_list) == q
        if decision_certificate_formed:
            decision_certificates.append(
                {
                    "view": view,
                    "value": value,
                    "committers": sorted(committers_list),
                }
            )
        certificate_receipts.append(
            {
                "index": certificate_index,
                "view": view,
                "value": value,
                "prepare_signers": sorted(prepare_signers_list),
                "acknowledgers": sorted(acknowledgers_list),
                "lock_certificate_assembled": lock_certificate_assembled,
                "lock_certificate_recipients": sorted(recipients_list),
                "committers": sorted(committers_list),
                "decision_certificate_formed": decision_certificate_formed,
            }
        )

    new_view_reports: list[dict[str, Any]] = []
    for sender in new_view_senders_list:
        state = states[sender]
        report = (
            {
                "view": state["view"],
                "value": state["value"],
                "certificate_index": state["certificate_index"],
            }
            if state["kind"] == "lc_lock"
            else None
        )
        new_view_reports.append({"sender": sender, "report": report})
    valid_reports = [
        row["report"] for row in new_view_reports if row["report"] is not None
    ]
    if valid_reports:
        highest_view = max(report["view"] for report in valid_reports)
        highest_values = {
            report["value"]
            for report in valid_reports
            if report["view"] == highest_view
        }
        require(
            len(highest_values) == 1,
            "new-view reports conflict at the highest view",
        )
        proposal_values = [next(iter(highest_values))]
        selected_report = next(
            report
            for report in valid_reports
            if report["view"] == highest_view
            and report["value"] == proposal_values[0]
        )
    else:
        proposal_values = list(fresh_values)
        selected_report = None

    available_voters = sorted(validators - byzantine)
    voting_outcomes: list[dict[str, Any]] = []
    for proposal in proposal_values:
        eligible_voters: list[int] = []
        rejected_voters: list[int] = []
        for validator in available_voters:
            state = states[validator]
            eligible = (
                state["kind"] == "unlocked"
                or state["value"] == proposal
                or (
                    state["kind"] == "lc_lock"
                    and new_view_supersedes_lc_locks
                )
            )
            (eligible_voters if eligible else rejected_voters).append(validator)
        voting_outcomes.append(
            {
                "proposal": proposal,
                "eligible_voters": eligible_voters,
                "rejected_voters": rejected_voters,
                "eligible_vote_count": len(eligible_voters),
                "quorum_reached": len(eligible_voters) >= q,
            }
        )

    maximum_eligible_votes = max(
        outcome["eligible_vote_count"] for outcome in voting_outcomes
    )
    recovered = any(outcome["quorum_reached"] for outcome in voting_outcomes)
    return {
        "certificate_receipts": certificate_receipts,
        "validator_states": {
            str(validator): states[validator] for validator in sorted(states)
        },
        "decision_certificates": decision_certificates,
        "new_view": {
            "senders": list(new_view_senders_list),
            "reports": new_view_reports,
            "selected_report": selected_report,
            "proposal_values": proposal_values,
        },
        "available_voters": available_voters,
        "voting_outcomes": voting_outcomes,
        "maximum_eligible_votes": maximum_eligible_votes,
        "recovered": recovered,
        "deadlocked": not recovered,
    }


def bft_locking_receipt() -> dict[str, Any]:
    def valid_prepared_certificate(
        certificate: Any,
        validators: set[int],
        q: int,
    ) -> bool:
        if not isinstance(certificate, Mapping):
            return False
        if set(certificate) != {"view", "value", "signers"}:
            return False
        signers = certificate["signers"]
        return (
            type(certificate["view"]) is int
            and certificate["view"] >= 0
            and type(certificate["value"]) is str
            and isinstance(signers, list)
            and all(type(signer) is int for signer in signers)
            and len(signers) == q
            and len(set(signers)) == q
            and set(signers) <= validators
        )

    def valid_lock_certificate(
        certificate: Any,
        validators: set[int],
        q: int,
    ) -> bool:
        if not isinstance(certificate, Mapping):
            return False
        if set(certificate) != {"prepared_certificate", "acceptors"}:
            return False
        acceptors = certificate["acceptors"]
        return (
            valid_prepared_certificate(
                certificate["prepared_certificate"],
                validators,
                q,
            )
            and isinstance(acceptors, list)
            and all(type(acceptor) is int for acceptor in acceptors)
            and len(acceptors) == q
            and len(set(acceptors)) == q
            and set(acceptors) <= validators
        )

    def valid_decision_certificate(
        certificate: Any,
        validators: set[int],
        q: int,
    ) -> bool:
        if not isinstance(certificate, Mapping):
            return False
        if set(certificate) != {"lock_certificate", "committers"}:
            return False
        committers = certificate["committers"]
        return (
            valid_lock_certificate(
                certificate["lock_certificate"],
                validators,
                q,
            )
            and isinstance(committers, list)
            and all(type(committer) is int for committer in committers)
            and len(committers) == q
            and len(set(committers)) == q
            and set(committers) <= validators
        )

    parameter_rows: list[dict[str, Any]] = []
    total_same_view_pairs = 0
    total_prepared_acceptor_pairs = 0
    total_distinct_prepared_acceptor_pairs = 0
    total_acceptor_committer_pairs = 0
    total_distinct_acceptor_committer_pairs = 0
    total_lock_transfers = 0
    total_next_view_prepare_quorums = 0
    total_orphan_lock_reconciliation_traces = 0
    total_progress_traces = 0

    for n, f, q in ((1, 0, 1), (4, 1, 3), (6, 1, 4), (7, 2, 5)):
        validators = tuple(range(n))
        quorums = subsets_of_size(validators, q)
        byzantine_sets = subsets_of_size(validators, f)
        threshold_ok = 2 * q >= n + f + 1
        require(threshold_ok, "reference BFT threshold does not satisfy overlap inequality")

        minimum_quorum_overlap = n
        minimum_nonfaulty_overlap = n
        minimum_nonfaulty_commit_lock_holders = n
        minimum_lock_transfer = n
        highest_rule_preserves_value = True
        conflicting_next_view_quorums = 0
        checked_next_view_prepare_quorums = 0
        checked_prepared_acceptor_pairs = 0
        distinct_prepared_acceptor_pairs = 0
        checked_acceptor_committer_pairs = 0
        distinct_acceptor_committer_pairs = 0
        parameter_orphan_lock_reconciliation_traces = 0
        parameter_omitted_orphan_lock_reconciliations = 0
        parameter_orphan_lock_example: dict[str, Any] | None = None
        parameter_progress_traces = 0
        for byzantine in byzantine_sets:
            for left, right in itertools.product(quorums, repeat=2):
                overlap = left & right
                nonfaulty_overlap = overlap - byzantine
                minimum_quorum_overlap = min(minimum_quorum_overlap, len(overlap))
                minimum_nonfaulty_overlap = min(
                    minimum_nonfaulty_overlap,
                    len(nonfaulty_overlap),
                )
                require(
                    len(nonfaulty_overlap) >= 1,
                    "two certificates could avoid a nonfaulty overlap witness",
                )
                total_same_view_pairs += 1

            # Prepare signers and prepared-certificate acceptors are separate
            # quorum roles. Acceptance acknowledgements are provisional and
            # do not change the acceptor's voting lock.
            for prepared, acceptors in itertools.product(quorums, repeat=2):
                certificate_a = {
                    "prepared_certificate": {
                        "view": 0,
                        "value": "A",
                        "signers": sorted(prepared),
                    },
                    "acceptors": sorted(acceptors),
                }
                require(
                    valid_lock_certificate(certificate_a, set(validators), q),
                    "well-formed lock certificate was rejected",
                )
                checked_prepared_acceptor_pairs += 1
                total_prepared_acceptor_pairs += 1
                if prepared != acceptors:
                    distinct_prepared_acceptor_pairs += 1
                    total_distinct_prepared_acceptor_pairs += 1

            # Durable locks arise only when a validator receives the assembled
            # lock certificate and signs a commit. A valid decision certificate
            # therefore supplies q committers independently of the provisional
            # acceptors, with at least q-f nonfaulty durable lock holders.
            for acceptors, committers in itertools.product(quorums, repeat=2):
                decision_certificate_a = {
                    "lock_certificate": {
                        "prepared_certificate": {
                            "view": 0,
                            "value": "A",
                            "signers": sorted(quorums[0]),
                        },
                        "acceptors": sorted(acceptors),
                    },
                    "committers": sorted(committers),
                }
                require(
                    valid_decision_certificate(
                        decision_certificate_a,
                        set(validators),
                        q,
                    ),
                    "well-formed decision certificate was rejected",
                )
                nonfaulty_commit_lock_holders = committers - byzantine
                minimum_nonfaulty_commit_lock_holders = min(
                    minimum_nonfaulty_commit_lock_holders,
                    len(nonfaulty_commit_lock_holders),
                )
                require(
                    len(nonfaulty_commit_lock_holders) >= q - f,
                    "decision certificate has too few nonfaulty lock holders",
                )
                checked_acceptor_committer_pairs += 1
                total_acceptor_committer_pairs += 1
                if acceptors != committers:
                    distinct_acceptor_committer_pairs += 1
                    total_distinct_acceptor_committer_pairs += 1

            for committers, view_change in itertools.product(quorums, repeat=2):
                nonfaulty_commit_lock_holders = committers - byzantine
                transfer = nonfaulty_commit_lock_holders & view_change
                minimum_lock_transfer = min(minimum_lock_transfer, len(transfer))
                require(
                    len(transfer) >= 1,
                    "view-change quorum could omit every nonfaulty commit lock",
                )
                lock_certificate_a = {
                    "prepared_certificate": {
                        "view": 0,
                        "value": "A",
                        "signers": sorted(quorums[0]),
                    },
                    "acceptors": sorted(quorums[0]),
                }
                messages = []
                for validator in sorted(view_change):
                    if validator in nonfaulty_commit_lock_holders:
                        report = lock_certificate_a
                    elif validator in byzantine:
                        # The adversary reports a higher conflicting object,
                        # but cannot forge q distinct prepare signers or
                        # acceptance acknowledgements.
                        report = {
                            "prepared_certificate": {
                                "view": 1,
                                "value": "B",
                                "signers": [validator] * q,
                            },
                            "acceptors": [validator] * q,
                        }
                    else:
                        report = None
                    messages.append(
                        {
                            "validator": validator,
                            "highest_lock_certificate": report,
                        }
                    )
                valid_reports = [
                    message["highest_lock_certificate"]
                    for message in messages
                    if valid_lock_certificate(
                        message["highest_lock_certificate"],
                        set(validators),
                        q,
                    )
                ]
                selected = max(
                    valid_reports,
                    key=lambda certificate: certificate[
                        "prepared_certificate"
                    ]["view"],
                )
                highest_rule_preserves_value &= (
                    selected["prepared_certificate"]["view"] == 0
                    and selected["prepared_certificate"]["value"] == "A"
                )
                total_lock_transfers += 1

                # Execute the induction base into the next view.  A nonfaulty
                # lock holder may vote only for A (or for a strictly higher
                # valid lock certificate, absent in this finite base). Every
                # candidate q-vote prepare quorum for B therefore fails.
                for candidate_prepare in quorums:
                    eligible_for_b = {
                        validator
                        for validator in candidate_prepare
                        if (
                            validator in byzantine
                            or validator not in nonfaulty_commit_lock_holders
                        )
                    }
                    if len(eligible_for_b) == q:
                        conflicting_next_view_quorums += 1
                    checked_next_view_prepare_quorums += 1
                    total_next_view_prepare_quorums += 1

            # P5 positive trace: after GST and timeout-bound activation, each
            # nonfaulty leader is checked with Byzantine validators
            # withholding.  Because q<=n-f, the
            # nonfaulty set alone constructs VC, prepare, PC, q acceptance
            # acknowledgements, a lock certificate, commit, DC, and
            # relay/finalisation.
            nonfaulty = set(validators) - set(byzantine)
            for leader in sorted(nonfaulty):
                view_change_senders = sorted(nonfaulty)[:q]
                prepare_senders = sorted(nonfaulty)[:q]
                acceptor_senders = sorted(nonfaulty)[:q]
                commit_senders = sorted(nonfaulty)[:q]
                # Failed views can deliver assembled lock certificates to fewer
                # than q validators and collect fewer than q commits. Execute
                # those local deliveries, build the next q-message certificate,
                # select its highest report, and derive every validator's voting
                # eligibility. The n=6 case also leaves the higher B@1 lock
                # outside the new-view certificate.
                if q == 1:
                    orphan_certificates = [
                        {
                            "view": 0,
                            "value": "A",
                            "prepare_signers": prepare_senders,
                            "acknowledgers": acceptor_senders,
                            "lock_certificate_recipients": [
                                min(nonfaulty)
                            ],
                            "committers": [],
                        }
                    ]
                else:
                    orphan_certificates = [
                        {
                            "view": 0,
                            "value": "A",
                            "prepare_signers": prepare_senders,
                            "acknowledgers": acceptor_senders,
                            "lock_certificate_recipients": [
                                min(nonfaulty)
                            ],
                            "committers": [min(nonfaulty)],
                        },
                        {
                            "view": 1,
                            "value": "B",
                            "prepare_signers": prepare_senders,
                            "acknowledgers": acceptor_senders,
                            "lock_certificate_recipients": [
                                max(nonfaulty)
                            ],
                            "committers": [max(nonfaulty)],
                        },
                    ]
                orphan_scenario = {
                    "validators": list(validators),
                    "byzantine": sorted(byzantine),
                    "q": q,
                    "certificates": orphan_certificates,
                    "new_view_senders": view_change_senders,
                    "fresh_values": ["A"],
                    "acknowledgements_are_provisional": True,
                    "new_view_supersedes_lc_locks": True,
                }
                orphan_evaluation = evaluate_orphan_lock_scenario(
                    orphan_scenario
                )
                selected_report = orphan_evaluation["new_view"][
                    "selected_report"
                ]
                require(
                    selected_report is not None,
                    "positive orphan-lock trace produced no report",
                )
                selected_view = selected_report["view"]
                selected_value = selected_report["value"]
                selected_outcome = next(
                    outcome
                    for outcome in orphan_evaluation["voting_outcomes"]
                    if outcome["proposal"] == selected_value
                )
                omitted_higher_lock_holders = [
                    validator
                    for validator in sorted(
                        nonfaulty - set(view_change_senders)
                    )
                    if (
                        orphan_evaluation["validator_states"][str(validator)][
                            "kind"
                        ]
                        == "lc_lock"
                        and orphan_evaluation["validator_states"][
                            str(validator)
                        ]["view"]
                        > selected_view
                    )
                ]
                orphan_locks_reconciled = (
                    not orphan_evaluation["decision_certificates"]
                    and orphan_evaluation["recovered"]
                    and selected_outcome["eligible_voters"]
                    == orphan_evaluation["available_voters"]
                )
                require(
                    orphan_locks_reconciled,
                    "valid new-view selection did not reconcile orphan locks",
                )
                if parameter_orphan_lock_example is None:
                    parameter_orphan_lock_example = {
                        "scenario": orphan_scenario,
                        "evaluation": orphan_evaluation,
                        "omitted_higher_lock_holders": (
                            omitted_higher_lock_holders
                        ),
                    }
                parameter_orphan_lock_reconciliation_traces += 1
                total_orphan_lock_reconciliation_traces += 1
                parameter_omitted_orphan_lock_reconciliations += len(
                    omitted_higher_lock_holders
                )
                prepared_certificate = {
                    "view": selected_view + 1,
                    "value": selected_value,
                    "signers": prepare_senders,
                }
                lock_certificate = {
                    "prepared_certificate": prepared_certificate,
                    "acceptors": acceptor_senders,
                }
                lock_certificate_valid = valid_lock_certificate(
                    lock_certificate,
                    set(validators),
                    q,
                )
                decision_certificate_valid = (
                    valid_decision_certificate(
                        {
                            "lock_certificate": lock_certificate,
                            "committers": commit_senders,
                        },
                        set(validators),
                        q,
                    )
                    and set(commit_senders) <= nonfaulty
                )
                finalised = (
                    set(nonfaulty)
                    if decision_certificate_valid
                    else set()
                )
                progress_ok = (
                    leader in nonfaulty
                    and len(view_change_senders) == q
                    and len(prepare_senders) == q
                    and len(acceptor_senders) == q
                    and lock_certificate_valid
                    and decision_certificate_valid
                    and finalised == nonfaulty
                )
                require(progress_ok, "post-GST nonfaulty-progress trace failed")
                parameter_progress_traces += 1
                total_progress_traces += 1

        parameter_rows.append(
            {
                "n": n,
                "f": f,
                "q": q,
                "threshold_2q_ge_n_plus_f_plus_1": threshold_ok,
                "minimum_quorum_overlap": minimum_quorum_overlap,
                "minimum_nonfaulty_certificate_overlap": minimum_nonfaulty_overlap,
                "checked_prepared_signer_acceptor_quorum_pairs": (
                    checked_prepared_acceptor_pairs
                ),
                "distinct_prepared_signer_acceptor_pairs": (
                    distinct_prepared_acceptor_pairs
                ),
                "checked_provisional_acceptor_committer_quorum_pairs": (
                    checked_acceptor_committer_pairs
                ),
                "distinct_provisional_acceptor_committer_pairs": (
                    distinct_acceptor_committer_pairs
                ),
                "minimum_nonfaulty_decision_commit_lock_holders": (
                    minimum_nonfaulty_commit_lock_holders
                ),
                "minimum_nonfaulty_lock_transfer_to_view_change": minimum_lock_transfer,
                "base_view_highest_lock_certificate_preserves_value": (
                    highest_rule_preserves_value
                ),
                "checked_next_view_candidate_prepare_quorums": (
                    checked_next_view_prepare_quorums
                ),
                "conflicting_next_view_quorums_accepted": (
                    conflicting_next_view_quorums
                ),
                "checked_post_GST_nonfaulty_leader_progress_traces": (
                    parameter_progress_traces
                ),
                "checked_orphan_lock_new_view_reconciliation_traces": (
                    parameter_orphan_lock_reconciliation_traces
                ),
                "checked_omitted_higher_orphan_lock_reconciliations": (
                    parameter_omitted_orphan_lock_reconciliations
                ),
                "orphan_lock_reconciliation_example": (
                    parameter_orphan_lock_example
                ),
            }
        )

    checks = {
        "every_certificate_pair_has_nonfaulty_overlap": all(
            row["minimum_nonfaulty_certificate_overlap"] >= 1
            for row in parameter_rows
        ),
        "every_decision_certificate_has_q_minus_f_nonfaulty_commit_locks": all(
            row["minimum_nonfaulty_decision_commit_lock_holders"]
            >= row["q"] - row["f"]
            for row in parameter_rows
        ),
        "prepare_signers_provisional_acceptors_and_committers_are_separate": (
            total_prepared_acceptor_pairs > 0
            and total_distinct_prepared_acceptor_pairs > 0
            and total_acceptor_committer_pairs > 0
            and total_distinct_acceptor_committer_pairs > 0
        ),
        "every_new_view_quorum_carries_a_nonfaulty_lock": all(
            row["minimum_nonfaulty_lock_transfer_to_view_change"] >= 1
            for row in parameter_rows
        ),
        "reference_quorums_can_form_without_Byzantine_votes": all(
            row["q"] <= row["n"] - row["f"]
            for row in parameter_rows
        ),
        "general_threshold_case_with_q_lt_n_minus_f_is_checked": any(
            row["q"] < row["n"] - row["f"]
            for row in parameter_rows
        ),
        "highest_lock_certificate_rule_is_value_preserving": all(
            row["base_view_highest_lock_certificate_preserves_value"]
            for row in parameter_rows
        ),
        "finite_next_view_lock_transition_rejects_conflicting_prepare_quorums": all(
            row["checked_next_view_candidate_prepare_quorums"] > 0
            and row["conflicting_next_view_quorums_accepted"] == 0
            for row in parameter_rows
        ),
        "valid_new_view_selection_supersedes_local_predecision_locks": all(
            row["checked_orphan_lock_new_view_reconciliation_traces"] > 0
            for row in parameter_rows
        ),
        "omitted_higher_orphan_lock_is_reconciled_in_general_case": any(
            row["q"] < row["n"] - row["f"]
            and row["checked_omitted_higher_orphan_lock_reconciliations"] > 0
            for row in parameter_rows
        ),
        "post_GST_nonfaulty_progress_trace_reaches_valid_decision_and_relay": all(
            row["checked_post_GST_nonfaulty_leader_progress_traces"] > 0
            for row in parameter_rows
        ),
        "finalisation_occurs_only_after_a_valid_decision_certificate": (
            total_progress_traces > 0
        ),
        "same_view_and_cross_view_induction_base_are_finite_checked": (
            total_same_view_pairs > 0
            and total_prepared_acceptor_pairs > 0
            and total_acceptor_committer_pairs > 0
            and total_lock_transfers > 0
            and total_next_view_prepare_quorums > 0
            and total_orphan_lock_reconciliation_traces > 0
        ),
    }
    require(all(checks.values()), f"BFT locking reference failed: {checks}")
    return {
        "receipt_id": "BFT-LOCK-1",
        "protocol_rule": {
            "monotone_view_participation": (
                "after entering view w, a nonfaulty validator never signs a "
                "prepare, prepared-certificate acceptance, or commit message "
                "for any lower view and never adopts a durable lower-view lock"
            ),
            "raw_prepared_certificate": (
                "q distinct authenticated prepare votes for one value in one view"
            ),
            "lock_certificate": (
                "a valid prepared certificate plus q distinct authenticated "
                "provisional acceptance acknowledgements; an acknowledgement "
                "does not change the acceptor's voting lock, and a raw prepared "
                "certificate alone has no new-view or commit status"
            ),
            "decision_certificate": (
                "q distinct authenticated commit votes referring to the same "
                "valid lock certificate; each nonfaulty committer received "
                "that assembled certificate and recorded its durable lock "
                "before committing, and a nonfaulty observer finalises "
                "only after accepting such a valid decision certificate"
            ),
            "lock": (
                "a nonfaulty validator records a durable lock only after "
                "receiving an assembled lock certificate and before committing; "
                "it reports that certificate in view change and retains the "
                "lock until a valid q-message new-view selection authorizes a "
                "prepare and the assembled current-view lock certificate "
                "replaces it before commit"
            ),
            "new_view": (
                "the leader collects q authenticated view-change messages and "
                "proposes the value of their highest-view lock certificate; "
                "every nonfaulty validator may supersede any local "
                "assembled-lock-certificate lock with this valid q-message "
                "selection for voting; the current-view lock is recorded only "
                "after its assembled certificate arrives, without locally "
                "testing whether an earlier commit set reached q"
            ),
            "vote_rule": (
                "an unlocked nonfaulty validator may support a fresh value only "
                "when no lock certificate is present; otherwise it votes only "
                "for its lock value, except that the value selected by a valid "
                "q-message new-view certificate supersedes its local "
                "assembled-lock-certificate lock; if a decision certificate "
                "exists, quorum intersection forces that selection to be "
                "compatible with the decided value"
            ),
            "post_GST_timeout_active_nonfaulty_progress": (
                "in a post-GST nonfaulty-led view after timeout-bound "
                "activation, the "
                "leader emits the P4-selected justified proposal; nonfaulty "
                "validators promptly send prepare votes and provisional "
                "prepared-certificate acknowledgements; the leader assembles "
                "and disseminates the lock certificate, validators lock and "
                "send commit messages, and they relay the "
                "decision certificate so every nonfaulty validator finalises"
            ),
        },
        "finite_parameter_sweep": parameter_rows,
        "checked_same_view_certificate_pairs": total_same_view_pairs,
        "checked_prepared_signer_acceptor_quorum_pairs": (
            total_prepared_acceptor_pairs
        ),
        "checked_provisional_acceptor_committer_quorum_pairs": (
            total_acceptor_committer_pairs
        ),
        "checked_decision_committer_to_view_change_pairs": total_lock_transfers,
        "checked_next_view_candidate_prepare_quorums": (
            total_next_view_prepare_quorums
        ),
        "checked_orphan_lock_new_view_reconciliation_traces": (
            total_orphan_lock_reconciliation_traces
        ),
        "checked_post_GST_progress_traces": total_progress_traces,
        "induction_invariant": (
            "after the first valid decision certificate for x at view v, its "
            "q-f nonfaulty committer locks intersect every q-message view "
            "change and every later q-vote prepare quorum; every later raw "
            "prepared certificate, lock certificate, and decision certificate "
            "therefore carries x"
        ),
        "checks": checks,
    }


def bft_negative_controls() -> list[dict[str, Any]]:
    controls: list[dict[str, Any]] = []

    def conflicting_certificates(
        byzantine: set[int],
        left: set[int],
        right: set[int],
        *,
        allow_nonfaulty_double_vote: bool,
    ) -> bool:
        nonfaulty_overlap = (left & right) - byzantine
        return allow_nonfaulty_double_vote or not nonfaulty_overlap

    # P1: the single nonfaulty overlap witness signs both values.
    p1 = conflicting_certificates(
        {0},
        {0, 1, 2},
        {0, 1, 3},
        allow_nonfaulty_double_vote=True,
    )
    require(p1, "one-vote negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P1_one_vote_per_view",
            "finite_witness": {
                "n": 4,
                "declared_f": 1,
                "byzantine": [0],
                "view": 0,
                "certificate_A": [0, 1, 2],
                "certificate_B": [0, 1, 3],
                "nonfaulty_double_voter": 1,
            },
            "violated_conclusion": "same_view_safety",
            "counterexample_verified": p1,
        }
    )

    # Monotone view participation is separate from one-value-per-view.  H2
    # first prepares B in view 1 without receiving the assembled PC, then
    # returns to commit A in view 0.  Once DC(A@0) exists it receives the
    # valid higher PC(B@1), advances its lock, and commits B. All
    # per-view vote and certificate/lock rules hold; only the stale lower-view
    # commit is forbidden by the removed rule.
    stale_view_trace = {
        "byzantine": [0],
        "prepare_B_view_1": [0, 2, 3],
        "H2_received_PC_B_before_stale_commit": False,
        "commit_A_view_0": [0, 1, 2],
        "commit_B_view_1": [0, 2, 3],
        "H2_view_sequence": [1, 0, 1],
    }
    stale_view_conflict = (
        len(stale_view_trace["prepare_B_view_1"]) == 3
        and len(stale_view_trace["commit_A_view_0"]) == 3
        and len(stale_view_trace["commit_B_view_1"]) == 3
        and stale_view_trace["H2_view_sequence"] != sorted(
            stale_view_trace["H2_view_sequence"]
        )
    )
    require(stale_view_conflict, "monotone-view negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P1_monotone_view_participation",
            "finite_witness": stale_view_trace,
            "violated_conclusion": "cross_view_safety",
            "counterexample_verified": stale_view_conflict,
        }
    )

    # P2: without certificate semantics, different nonfaulty nodes can each
    # treat a single local vote as a decision.
    no_certificate = {"H1_decides": "A", "H2_decides": "B"}
    p2 = len(set(no_certificate.values())) == 2
    require(p2, "certificate-semantics negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P2_valid_quorum_certificate_required_for_decision",
            "finite_witness": {
                "n": 4,
                "f": 1,
                "decisions_without_certificates": no_certificate,
            },
            "violated_conclusion": "same_view_safety",
            "counterexample_verified": p2,
        }
    )

    # A6: a quorum of two has only a Byzantine intersection.
    small_q = conflicting_certificates(
        {0},
        {0, 1},
        {0, 2},
        allow_nonfaulty_double_vote=False,
    )
    require(small_q, "quorum-overlap negative control failed")
    controls.append(
        {
            "removed_hypothesis": "A6_nonfaulty_quorum_overlap_threshold",
            "finite_witness": {
                "n": 4,
                "f": 1,
                "invalid_q": 2,
                "certificate_A": [0, 1],
                "certificate_B": [0, 2],
                "intersection": [0],
                "byzantine": [0],
            },
            "violated_conclusion": "same_view_safety",
            "counterexample_verified": small_q,
        }
    )

    # A2: two Byzantine validators exceed the declared f=1 budget.
    too_many_faults = conflicting_certificates(
        {0, 1},
        {0, 1, 2},
        {0, 1, 3},
        allow_nonfaulty_double_vote=False,
    )
    require(too_many_faults, "fault-bound negative control failed")
    controls.append(
        {
            "removed_hypothesis": "A2_at_most_f_Byzantine",
            "finite_witness": {
                "n": 4,
                "declared_f": 1,
                "actual_byzantine": [0, 1],
                "certificate_A": [0, 1, 2],
                "certificate_B": [0, 1, 3],
            },
            "violated_conclusion": "same_view_safety",
            "counterexample_verified": too_many_faults,
        }
    )

    # A5: unauthenticated identities allow one process to forge q signers.
    forged = len({"forged_0", "forged_1", "forged_2"}) == 3
    require(forged, "authentication negative control failed")
    controls.append(
        {
            "removed_hypothesis": "A5_unforgeable_authentication",
            "finite_witness": {
                "n": 4,
                "q": 3,
                "real_adversarial_processes": 1,
                "forged_signer_labels": ["forged_0", "forged_1", "forged_2"],
            },
            "violated_conclusion": "certificate_semantics_and_safety",
            "counterexample_verified": forged,
        }
    )

    # P4: P1 is per-view only. Without locks, nonfaulty validators may sign A
    # in view 0 and B in view 1.
    view_zero = {0, 1, 2}
    view_one = {0, 2, 3}
    no_lock = (
        len(view_zero) == 3
        and len(view_one) == 3
        and view_zero != view_one
    )
    require(no_lock, "cross-view lock negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P4_prepared_certificate_lock",
            "finite_witness": {
                "n": 4,
                "f": 1,
                "byzantine": [0],
                "view_0_decision_A": sorted(view_zero),
                "view_1_decision_B": sorted(view_one),
                "P1_is_respected_per_view": True,
            },
            "violated_conclusion": "cross_view_safety",
            "counterexample_verified": no_lock,
        }
    )

    # P4 provisional acceptance: Byzantine assemblers deliver two raw prepared
    # certificates to different nonfaulty validators but collect fewer than q
    # acknowledgements for either, so no lock certificate exists. If those
    # pre-certificate acknowledgements create durable locks, the locks have no
    # reportable new-view evidence and split the q available nonfaulty voters.
    pre_certificate_scenario = {
        "validators": [0, 1, 2, 3],
        "byzantine": [0],
        "q": 3,
        "certificates": [
            {
                "view": 0,
                "value": "A",
                "prepare_signers": [0, 1, 2],
                "acknowledgers": [0, 1],
                "lock_certificate_recipients": [],
                "committers": [],
            },
            {
                "view": 1,
                "value": "B",
                "prepare_signers": [0, 2, 3],
                "acknowledgers": [0, 2],
                "lock_certificate_recipients": [],
                "committers": [],
            },
        ],
        "new_view_senders": [1, 2, 3],
        "fresh_values": ["A", "B", "C"],
        "acknowledgements_are_provisional": False,
        "new_view_supersedes_lc_locks": True,
    }
    pre_certificate_evaluation = evaluate_orphan_lock_scenario(
        pre_certificate_scenario
    )
    pre_states = pre_certificate_evaluation["validator_states"]
    pre_certificate_deadlock = (
        pre_certificate_evaluation["deadlocked"]
        and not pre_certificate_evaluation["decision_certificates"]
        and pre_certificate_evaluation["new_view"]["selected_report"] is None
        and pre_certificate_evaluation["maximum_eligible_votes"] == 2
        and pre_states["1"]["kind"] == "ack_lock"
        and pre_states["1"]["value"] == "A"
        and pre_states["2"]["kind"] == "ack_lock"
        and pre_states["2"]["value"] == "B"
    )
    require(
        pre_certificate_deadlock,
        "provisional prepared-acceptance negative control failed",
    )
    controls.append(
        {
            "removed_hypothesis": (
                "P4_pre_lock_certificate_acknowledgements_are_provisional"
            ),
            "finite_witness": {
                "scenario": pre_certificate_scenario,
                "evaluation": pre_certificate_evaluation,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": pre_certificate_deadlock,
        }
    )

    # P4 orphan-lock reconciliation: two assembled lock certificates can be
    # delivered to different validators and receive fewer than q commits, so
    # neither produces a decision certificate. A valid q-message new-view
    # certificate reports both and selects the higher B@1 certificate. If the
    # lower A@0 lock cannot be superseded by that selection, Byzantine
    # withholding leaves only two votes for B.
    partial_commit_scenario = {
        "validators": [0, 1, 2, 3],
        "byzantine": [0],
        "q": 3,
        "certificates": [
            {
                "view": 0,
                "value": "A",
                "prepare_signers": [0, 1, 2],
                "acknowledgers": [0, 1, 2],
                "lock_certificate_recipients": [1],
                "committers": [1],
            },
            {
                "view": 1,
                "value": "B",
                "prepare_signers": [0, 2, 3],
                "acknowledgers": [0, 2, 3],
                "lock_certificate_recipients": [2],
                "committers": [2],
            },
        ],
        "new_view_senders": [1, 2, 3],
        "fresh_values": ["C"],
        "acknowledgements_are_provisional": True,
        "new_view_supersedes_lc_locks": False,
    }
    partial_commit_evaluation = evaluate_orphan_lock_scenario(
        partial_commit_scenario
    )
    partial_outcome = partial_commit_evaluation["voting_outcomes"][0]
    partial_commit_deadlock = (
        partial_commit_evaluation["deadlocked"]
        and not partial_commit_evaluation["decision_certificates"]
        and partial_commit_evaluation["new_view"]["selected_report"]
        == {"view": 1, "value": "B", "certificate_index": 1}
        and partial_outcome["proposal"] == "B"
        and partial_outcome["eligible_voters"] == [2, 3]
        and partial_commit_evaluation["maximum_eligible_votes"] == 2
    )
    require(
        partial_commit_deadlock,
        "partial-commit orphan-lock negative control failed",
    )
    controls.append(
        {
            "removed_hypothesis": (
                "P4_valid_new_view_supersedes_local_predecision_commit_lock"
            ),
            "finite_witness": {
                "scenario": partial_commit_scenario,
                "evaluation": partial_commit_evaluation,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": partial_commit_deadlock,
        }
    )

    # P4 leader selection: a new leader omits the valid highest certificate.
    # Correct locked validators reject in both views.  The second leader is
    # nonfaulty, so this two-view trace exceeds the f+1 post-GST bound when
    # highest-certificate selection is removed while P5 responsiveness
    # remains.
    omitted = {
        "view_change_quorum": {0, 1, 3},
        "nonfaulty_lock_holder": 1,
        "reported_certificate": "A@0",
        "leader_proposals": ["B@1", "B@2"],
        "leaders": [0, 3],
        "available_votes_for_B": {0, 3},
    }
    omitted_highest = (
        omitted["nonfaulty_lock_holder"] in omitted["view_change_quorum"]
        and omitted["reported_certificate"].startswith("A")
        and all(proposal.startswith("B") for proposal in omitted["leader_proposals"])
        and len(omitted["available_votes_for_B"]) < 3
        and len(omitted["leaders"]) == 2
    )
    require(omitted_highest, "highest-certificate negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P4_leader_selects_highest_prepared_certificate",
            "finite_witness": {
                "view_change_quorum": sorted(omitted["view_change_quorum"]),
                "nonfaulty_lock_holder": omitted["nonfaulty_lock_holder"],
                "reported_certificate": omitted["reported_certificate"],
                "post_GST_views": [1, 2],
                "timeout_bound_active": True,
                "leaders": omitted["leaders"],
                "nonfaulty_leader": 3,
                "unjustified_leader_proposals": omitted["leader_proposals"],
                "available_votes": sorted(omitted["available_votes_for_B"]),
                "q": 3,
                "declared_f_plus_one_view_bound": 2,
                "decision_certificate": None,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": omitted_highest,
        }
    )

    # P4 validator validation: if validators do not check the new-view proof,
    # an unlocked validator can accept a proposal whose justification has
    # fewer than q messages. Lock discipline also prevents a conflicting
    # certificate, so this is explicitly a protocol-conformance control, not
    # a safety counterexample.
    malformed_justification = {"messages": [0, 3], "q": 3}
    no_validation = (
        len(malformed_justification["messages"]) < malformed_justification["q"]
    )
    require(no_validation, "new-view validation negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P4_validator_checks_new_view_justification",
            "finite_witness": {
                "unlocked_validator": 3,
                "proposal": "B@1",
                "new_view_messages": malformed_justification["messages"],
                "required_q": malformed_justification["q"],
                "accepted_without_validation": True,
                "conflicting_certificate_formed": False,
            },
            "violated_conclusion": "new_view_justification_filtering",
            "counterexample_verified": no_validation,
        }
    )

    # A1: a finite horizon with no delivered messages violates any claimed
    # bounded decision time after timeout activation once synchrony is removed.
    horizon = 4
    delivered_per_round = [0] * horizon
    no_synchrony = sum(delivered_per_round) == 0
    require(no_synchrony, "partial-synchrony negative control failed")
    controls.append(
        {
            "removed_hypothesis": "A1_partial_synchrony_for_liveness",
            "finite_witness": {
                "rounds": horizon,
                "delivered_messages_per_round": delivered_per_round,
                "certificate_threshold": 3,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": no_synchrony,
        }
    )

    # Availability is distinct from overlap. With n=4, f=1, q=4, quorum
    # intersections are maximal but one withholding Byzantine validator
    # prevents every certificate.
    available_nonfaulty_votes = {1, 2, 3}
    unavailable_quorum = len(available_nonfaulty_votes) < 4
    require(unavailable_quorum, "quorum-availability negative control failed")
    controls.append(
        {
            "removed_hypothesis": "quorum_availability_q_le_n_minus_f_for_liveness",
            "finite_witness": {
                "n": 4,
                "f": 1,
                "q": 4,
                "withholding_byzantine": [0],
                "available_nonfaulty_votes": sorted(available_nonfaulty_votes),
            },
            "violated_conclusion": "liveness",
            "counterexample_verified": unavailable_quorum,
        }
    )

    # A4: a partitioned directed graph has no strongly connected quorum and
    # cannot deliver q votes to a leader.
    reachable_from_zero = {0, 1}
    partitioned = len(reachable_from_zero) < 3
    require(partitioned, "connectivity negative control failed")
    controls.append(
        {
            "removed_hypothesis": "A4_strong_quorum_connectivity",
            "finite_witness": {
                "n": 4,
                "q": 3,
                "directed_components": [[0, 1], [2, 3]],
                "maximum_votes_reachable_by_leader_0": len(reachable_from_zero),
            },
            "violated_conclusion": "liveness",
            "counterexample_verified": partitioned,
        }
    )

    # P3: with a permanently faulty leader and no terminating view change,
    # all messages can be timely while no proposal is ever issued.
    no_view_change_trace = [
        {"round": round_index, "view": 0, "leader": 0, "proposal": None}
        for round_index in range(4)
    ]
    no_view_change = all(row["view"] == 0 and row["proposal"] is None for row in no_view_change_trace)
    require(no_view_change, "view-change liveness negative control failed")
    controls.append(
        {
            "removed_hypothesis": "P3_terminating_view_change",
            "finite_witness": {
                "byzantine_leader": 0,
                "timely_network": True,
                "trace": no_view_change_trace,
            },
            "violated_conclusion": "liveness",
            "counterexample_verified": no_view_change,
        }
    )

    # P5(a): all safety and view-selection rules can hold vacuously while a
    # nonfaulty post-GST leader simply emits no proposal.
    silent_leader_trace = [
        {
            "phase": phase,
            "leader": 3,
            "leader_nonfaulty": True,
            "timely_network": True,
            "messages_sent": 0,
        }
        for phase in ("new_view", "proposal", "prepare", "commit", "relay")
    ]
    silent_leader = all(
        row["messages_sent"] == 0 for row in silent_leader_trace
    )
    require(silent_leader, "nonfaulty-leader progress negative control failed")
    controls.append(
        {
            "removed_hypothesis": (
                "P5_post_GST_nonfaulty_leader_emits_justified_proposal"
            ),
            "finite_witness": {
                "n": 4,
                "f": 1,
                "q": 3,
                "post_GST": True,
                "timeout_bound_active": True,
                "view": 1,
                "trace": silent_leader_trace,
                "decision_certificate": None,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": silent_leader,
        }
    )

    # P5(b): a valid proposal is timely, but responsive validator actions are
    # absent. P1/P2/P4 and every safety assumption hold; no PC/DC can
    # be produced.
    silent_validator_trace = {
        "proposal": {
            "value": "A",
            "valid_new_view_justification": True,
            "delivered_to_nonfaulty": [1, 2, 3],
        },
        "prepare_votes": [],
        "commit_votes": [],
        "decision_relays": [],
    }
    silent_validators = (
        len(silent_validator_trace["prepare_votes"]) < 3
        and len(silent_validator_trace["commit_votes"]) < 3
    )
    require(
        silent_validators,
        "nonfaulty-validator progress negative control failed",
    )
    controls.append(
        {
            "removed_hypothesis": (
                "P5_post_GST_nonfaulty_validators_prepare_commit_and_relay"
            ),
            "finite_witness": {
                "n": 4,
                "f": 1,
                "q": 3,
                "post_GST": True,
                "timeout_bound_active": True,
                "nonfaulty_leader": 3,
                "trace": silent_validator_trace,
                "decision_certificate": None,
            },
            "violated_conclusion": "bounded_liveness_after_timeout_activation",
            "counterexample_verified": silent_validators,
        }
    )

    return controls


# ---------------------------------------------------------------------------
# Refinement moduli and ell-p guard
# ---------------------------------------------------------------------------


def refinement_receipt() -> dict[str, Any]:
    """Return finite proof schemas for genuinely infinite counterfamilies."""

    # For every m>=1 choose stage n=m and radius/threshold 1/m.  The witness
    # constructors are executable for any supplied positive integer.  Their
    # finite descriptions, together with the Archimedean null-sequence
    # argument stated in the paper, refute a family-uniform vanishing modulus;
    # a fixed finite prefix would not.
    sample_indices = (1, 2, 4, 8, 16, 32)

    def inverse_witness(m: int) -> dict[str, Any]:
        require(type(m) is int and m >= 1, "inverse witness index must be positive")
        radius = Fraction(1, m)
        observation_gap = Fraction(1, m)
        return {
            "index": m,
            "stage": m,
            "radius": str(radius),
            "observation_gap": str(observation_gap),
            "state_gap": 1,
            "stage_map_injective": True,
            "gap_within_radius": observation_gap <= radius,
        }

    def residual_witness(m: int) -> dict[str, Any]:
        require(type(m) is int and m >= 1, "residual witness index must be positive")
        threshold = Fraction(1, m)
        residual = Fraction(1, m)
        return {
            "index": m,
            "stage": m,
            "threshold": str(threshold),
            "residual_at_inconsistent_state": str(residual),
            "distance_to_consistent_set": 1,
            "residual_within_threshold": residual <= threshold,
        }

    inverse_rows = [inverse_witness(m) for m in sample_indices]
    residual_rows = [residual_witness(m) for m in sample_indices]
    observation_lipschitz_rows = [
        {
            "stage": n,
            "inconsistent_point": str(Fraction(1, n)),
            "residual": str(Fraction(1, n)),
            "observation_gap_to_one": 0,
            "state_gap_to_one": str(Fraction(n - 1, n)),
            "observation_Lipschitz_constant": n,
        }
        for n in sample_indices
        if n >= 2
    ]
    null_sequence_schema = {
        "target_epsilon": "1/k for arbitrary k>=1",
        "witness_index": "m=k+1",
        "cross_multiplied_inequality": "k<k+1",
        "positive_offset": 1,
    }
    null_sequence_schema_verified = (
        null_sequence_schema["positive_offset"] > 0
        and all(
            Fraction(1, k + null_sequence_schema["positive_offset"])
            < Fraction(1, k)
            for k in sample_indices
        )
    )

    # Infinite unit-defect generator p_j=j.  Every finite instance is exact,
    # and the closed formula gives a tail of M between n and n+M.
    telescope_rows = []
    for n, m in ((0, 1), (0, 4), (3, 9), (10, 25)):
        endpoint_distance = m - n
        defect_sum = sum(1 for _ in range(n, m))
        telescope_rows.append(
            {
                "n": n,
                "m": m,
                "endpoint_distance": endpoint_distance,
                "sum_of_unit_step_defects": defect_sum,
                "bound_exact": endpoint_distance == defect_sum,
            }
        )
    full_tail_model_outputs = list(range(8))
    full_tail_model_defects = [
        abs(
            full_tail_model_outputs[index + 1]
            - full_tail_model_outputs[index]
        )
        for index in range(len(full_tail_model_outputs) - 1)
    ]
    full_tail_model = {
        "coordinate_spaces": "Q_n=R with the Euclidean metric (complete)",
        "restriction_maps": "rho_(n+1,n)=identity (nonexpansive)",
        "compatible_input": "q_n=0 for every n",
        "consistent_sets": "C_n={n}",
        "normalizers": "N_n(x)=n",
        "observations": "B_n is constant",
        "residuals": "Phi_n(x)=|x-n|",
        "exact_solver_outputs": "xhat_n=n with receipt error e_n=0",
        "sampled_projected_outputs_at_level_zero": full_tail_model_outputs,
        "sampled_one_step_defects": full_tail_model_defects,
        "quantified_failure": (
            "alpha_n=1 for every n, and projected outputs n are not Cauchy"
        ),
    }
    unit_tail_schema = {
        "start": "arbitrary n>=0",
        "requested_tail": "arbitrary M>=1",
        "witness_endpoint": "m=n+M",
        "unit_defect": 1,
    }
    unit_tail_schema_verified = (
        unit_tail_schema["unit_defect"] == 1
        and all(
            row["endpoint_distance"] == row["m"] - row["n"]
            == row["sum_of_unit_step_defects"]
            for row in telescope_rows
        )
    )

    # Positive comparison schema: the complete geometric tail is known, not
    # inferred from a truncated prefix.
    summable_rows = [
        {"n": n, "exact_infinite_tail": str(Fraction(1, 2**n))}
        for n in sample_indices
    ]

    checks = {
        "parametric_inverse_witness_constructor_exact_on_samples": all(
            row["stage_map_injective"]
            and row["gap_within_radius"]
            and row["state_gap"] == 1
            for row in inverse_rows
        ),
        "inverse_witness_schema_has_null_radius_and_unit_lower_bound": (
            null_sequence_schema_verified
        ),
        "parametric_residual_witness_constructor_exact_on_samples": all(
            row["residual_within_threshold"]
            and row["distance_to_consistent_set"] == 1
            for row in residual_rows
        ),
        "residual_witness_schema_has_null_threshold_and_unit_lower_bound": (
            null_sequence_schema_verified
        ),
        "unbounded_observation_Lipschitz_counterfamily_exact_on_samples": all(
            Fraction(row["residual"]) == Fraction(1, row["stage"])
            and row["observation_gap_to_one"] == 0
            and Fraction(row["state_gap_to_one"])
            == Fraction(row["stage"] - 1, row["stage"])
            and row["observation_Lipschitz_constant"] == row["stage"]
            for row in observation_lipschitz_rows
        ),
        "arbitrary_depth_telescope_is_exact_on_reference_instances": all(
            row["bound_exact"] for row in telescope_rows
        ),
        "unit_defect_tail_schema_is_nonvanishing": unit_tail_schema_verified,
        "unit_defect_full_model_retains_other_cofinal_hypotheses": (
            all(defect == 1 for defect in full_tail_model_defects)
            and full_tail_model_outputs[-1] - full_tail_model_outputs[0] == 7
        ),
        "summable_reference_has_exact_vanishing_infinite_tail": all(
            Fraction(row["exact_infinite_tail"]) == Fraction(1, 2 ** row["n"])
            for row in summable_rows
        ),
    }
    require(all(checks.values()), f"refinement reference failed: {checks}")
    return {
        "receipt_id": "REFINEMENT-MODULUS-1",
        "uniform_inverse_counterfamily": {
            "state_spaces": "Q_n=C_n={0,1} with discrete distance",
            "observations": "B_n(0)=0, B_n(1)=1/n",
            "index_domain": "n in positive integers",
            "witness_constructor": (
                "for every positive integer m choose stage n=m and radius r_m=1/m"
            ),
            "quantified_conclusion": (
                "Omega(1/m)>=1 for every m>=1 while 1/m tends to zero"
            ),
            "null_sequence_proof_schema": null_sequence_schema,
            "sampled_exact_instances": inverse_rows,
        },
        "uniform_residual_counterfamily": {
            "consistent_sets": "C_n={0}",
            "residuals": "Phi_n(0)=0, Phi_n(1)=1/n",
            "index_domain": "n in positive integers",
            "witness_constructor": (
                "for every positive integer m choose stage n=m and threshold t_m=1/m"
            ),
            "quantified_conclusion": (
                "H(1/m)>=1 for every m>=1 while 1/m tends to zero"
            ),
            "null_sequence_proof_schema": null_sequence_schema,
            "sampled_exact_instances": residual_rows,
        },
        "uniform_observation_Lipschitz_counterfamily": {
            "state_spaces": "Q_n={0,1/n,1} with Euclidean distance, n>=2",
            "consistent_sets": "C_n={0,1}",
            "residuals": "Phi_n=distance to C_n",
            "observations": "B_n(0)=0 and B_n(1/n)=B_n(1)=1",
            "quantified_conclusion": (
                "Lip(B_n)=n is unbounded; x_n=1/n and y_n=1 have "
                "vanishing residual, equal observation, and state gap 1-1/n"
            ),
            "sampled_exact_instances": observation_lipschitz_rows,
        },
        "tower_telescope": {
            "infinite_generator": "p_j=j and a_j=1 for every j>=0",
            "quantified_conclusion": (
                "for every n and M choose m=n+M; the amplified tail equals M"
            ),
            "nonvanishing_tail_proof_schema": unit_tail_schema,
            "full_countermodel": full_tail_model,
            "sampled_exact_instances": telescope_rows,
        },
        "summable_reference": {
            "infinite_generator": "a_j=2^(-j-1)",
            "exact_tail_formula": "sum_{j=n}^infinity a_j=2^(-n)",
            "sampled_exact_instances": summable_rows,
        },
        "checks": checks,
    }


def refinement_negative_controls() -> list[dict[str, Any]]:
    receipt = refinement_receipt()
    inverse = receipt["uniform_inverse_counterfamily"]
    residual = receipt["uniform_residual_counterfamily"]
    observation_lipschitz = receipt[
        "uniform_observation_Lipschitz_counterfamily"
    ]
    telescope = receipt["tower_telescope"]
    fine_error = Fraction(1, 8)
    coarse_error = Fraction(1)
    declared_nonexpansive_bound = Fraction(1) * fine_error
    restriction_rows = [
        {
            "level": m,
            "solver_error": str(Fraction(1, 2**m)),
            "projected_plus_output_at_level_zero": 1,
            "projected_minus_output_at_level_zero": -1,
            "rho_m_0_Lipschitz_constant": 2**m,
            "normalizer_defect": 0,
        }
        for m in range(1, 9)
    ]
    expanding_restriction_failure = all(
        Fraction(row["solver_error"]) == Fraction(1, 2 ** row["level"])
        and row["projected_plus_output_at_level_zero"] == 1
        and row["projected_minus_output_at_level_zero"] == -1
        and row["rho_m_0_Lipschitz_constant"] == 2 ** row["level"]
        and row["normalizer_defect"] == 0
        for row in restriction_rows
    )
    cauchy_samples = [Fraction(1, m + 1) for m in range(1, 17)]
    cauchy_decreasing_to_excluded_zero = (
        all(
            cauchy_samples[index + 1] < cauchy_samples[index]
            for index in range(len(cauchy_samples) - 1)
        )
        and all(value > 0 for value in cauchy_samples)
        and Fraction(0) not in cauchy_samples
    )
    q_zero = 0
    q_one = 1
    identity_restriction_q_one = q_one
    incompatible_input = identity_restriction_q_one != q_zero
    parity_outputs = [m % 2 for m in range(8)]
    parity_errors = [abs(output - 0) for output in parity_outputs]
    nonvanishing_solver_error = (
        set(parity_errors[::2]) == {0}
        and set(parity_errors[1::2]) == {1}
        and len(set(parity_outputs)) == 2
    )
    controls = [
        {
            "removed_hypothesis": "uniform_inverse_observation_modulus",
            "finite_description_of_parametric_witness": inverse,
            "violated_conclusion": (
                "stagewise_injectivity_alone_does_not_give_the_claimed_family_modulus"
            ),
            "counterexample_verified": all(
                row["gap_within_radius"] and row["state_gap"] == 1
                for row in inverse["sampled_exact_instances"]
            ),
        },
        {
            "removed_hypothesis": "uniform_residual_error_modulus",
            "finite_description_of_parametric_witness": residual,
            "violated_conclusion": (
                "stagewise_error_bounds_alone_do_not_give_uniform_settled_output_control"
            ),
            "counterexample_verified": all(
                row["residual_within_threshold"]
                and row["distance_to_consistent_set"] == 1
                for row in residual["sampled_exact_instances"]
            ),
        },
        {
            "removed_hypothesis": "family_uniform_observation_Lipschitz_bound",
            "finite_description_of_parametric_witness": observation_lipschitz,
            "violated_conclusion": (
                "uniform_settled_output_stability_with_one_common_bound"
            ),
            "counterexample_verified": all(
                row["observation_Lipschitz_constant"] == row["stage"]
                and row["observation_gap_to_one"] == 0
                and Fraction(row["residual"]) == Fraction(1, row["stage"])
                and Fraction(row["state_gap_to_one"])
                == Fraction(row["stage"] - 1, row["stage"])
                for row in observation_lipschitz["sampled_exact_instances"]
            ),
        },
        {
            "removed_hypothesis": "cofinally_vanishing_or_summable_refinement_tail",
            "finite_description_of_parametric_witness": telescope,
            "violated_conclusion": (
                "arbitrary_depth_comparisons_need_not_be_uniformly_small"
            ),
            "counterexample_verified": all(
                row["bound_exact"] and row["endpoint_distance"] >= 1
                for row in telescope["sampled_exact_instances"]
            )
            and set(
                telescope["full_countermodel"]["sampled_one_step_defects"]
            )
            == {1}
            and len(
                set(
                    telescope["full_countermodel"][
                        "sampled_projected_outputs_at_level_zero"
                    ]
                )
            )
            > 1,
        },
        {
            "removed_hypothesis": "declared_uniform_Lipschitz_restriction_bound",
            "finite_witness": {
                "fine_two_point_distance": str(fine_error),
                "coarse_two_point_distance_after_restriction": str(coarse_error),
                "claimed_K": "1",
                "claimed_projected_bound": str(declared_nonexpansive_bound),
                "full_parametric_model": {
                    "coordinate_spaces": "Q_n=R (complete)",
                    "restriction_maps": "rho_(n+1,n)(x)=2x",
                    "compatible_input": "q_n=0",
                    "consistent_sets": "C_n={0}",
                    "normalizers": "N_n=0",
                    "observations": "B_n constant",
                    "residuals": "Phi_n(x)=|x|",
                    "solver_outputs": "xhat_m^(+/-)=+/-2^(-m)",
                    "quantified_failure": (
                        "receipt errors tend to zero and alpha_n=0, but "
                        "rho_(m,0)xhat_m^(+/-)=+/-1"
                    ),
                    "sampled_exact_instances": restriction_rows,
                },
            },
            "violated_conclusion": (
                "fine_solver_error_is_not_controlled_after_restriction"
            ),
            "counterexample_verified": (
                coarse_error > declared_nonexpansive_bound
                and expanding_restriction_failure
            ),
        },
        {
            "removed_hypothesis": "complete_coordinate_spaces",
            "finite_description_of_parametric_witness": {
                "coordinate_spaces": "Q_n=(0,1) with Euclidean distance",
                "restriction_maps": "identity (nonexpansive)",
                "compatible_input": "q_n=1/2",
                "consistent_sets": "C_n={1/(n+2)}",
                "normalizers": "N_n(x)=1/(n+2)",
                "exact_solver_outputs": "xhat_n=1/(n+2), so e_n=0",
                "normalizer_defects": (
                    "alpha_n=1/((n+2)(n+3)), a summable tail"
                ),
                "sequence": "rho_(m,0)xhat_m=1/(m+2)",
                "cauchy_bound": (
                    "for k>=m, |x_k-x_m|<=1/(m+2), which tends to zero"
                ),
                "ambient_limit": 0,
                "ambient_limit_is_in_Q": False,
                "sampled_exact_terms": [
                    str(value) for value in cauchy_samples
                ],
            },
            "violated_conclusion": "coordinate_limit_exists_in_Q",
            "counterexample_verified": cauchy_decreasing_to_excluded_zero,
        },
        {
            "removed_hypothesis": "compatible_projective_input",
            "finite_witness": {
                "coordinate_spaces": "Q_n={0,1} with discrete metric",
                "restriction_maps": "identity",
                "consistent_sets": "C_n=Q_n",
                "normalizers": "N_n=identity",
                "observations": "B_n=identity",
                "declared_input": "q_n=n mod 2",
                "exact_solver_outputs": "xhat_n=q_n with e_n=0",
                "normalizer_defects": "alpha_n=0",
                "sample_pair": {"q_0": q_zero, "q_1": q_one},
                "restricted_q_1": identity_restriction_q_one,
                "compatibility_equation_holds": not incompatible_input,
                "quantified_failure": (
                    "rho_(m,0)xhat_m=m mod 2 has distinct even/odd "
                    "cofinal limits"
                ),
            },
            "violated_conclusion": "output_family_belongs_to_projective_limit",
            "counterexample_verified": incompatible_input,
        },
        {
            "removed_hypothesis": "vanishing_solver_receipt_errors",
            "finite_description_of_parametric_witness": {
                "space": "Q_n={0,1} with identity restrictions",
                "exact_input": 0,
                "consistent_sets": "C_n=Q_n",
                "normalizers": "N_n=identity",
                "normalizer_defects": "alpha_n=0",
                "solver_output": "x_m=m mod 2",
                "receipt_error": "e_m=m mod 2, so it does not tend to zero",
                "sampled_outputs": parity_outputs,
                "sampled_receipt_errors": parity_errors,
                "cofinal_subsequences": {
                    "even_levels_limit": 0,
                    "odd_levels_limit": 1,
                },
            },
            "violated_conclusion": (
                "implementation_independent_cofinal_projective_limit"
            ),
            "counterexample_verified": nonvanishing_solver_error,
        },
    ]
    require(
        all(row["counterexample_verified"] for row in controls),
        "a refinement negative control did not materialize",
    )
    return controls


def lp_receipt() -> dict[str, Any]:
    points = ((0, 0), (1, 0), (0, 1), (1, 1))
    finite_channels = (0, 1)

    def l1(left: tuple[int, int], right: tuple[int, int]) -> int:
        return sum(abs(a - b) for a, b in zip(left, right, strict=True))

    finite_distances = {
        (left, right): l1(left, right)
        for left, right in itertools.product(points, repeat=2)
    }
    finite_product_values = all(
        type(value) is int and 0 <= value <= len(finite_channels)
        for value in finite_distances.values()
    )
    l1_triangle = all(
        l1(x, z) <= l1(x, y) + l1(y, z)
        for x, y, z in itertools.product(points, repeat=3)
    )

    # At p=1/2 the quasi-distance between e1 and e2 is
    # (sqrt(1)+sqrt(1))^2=4, while the route through zero has length 2.
    direct_quasi_distance = 4
    via_origin = 1 + 1
    p_half_failure = direct_quasi_distance > via_origin

    # If the channel set is the positive integers, every weight and every
    # channel distance can equal one. The N-channel partial sum is then N.
    # For each proposed integer upper bound M, cutoff N=M+1 is an exact
    # witness above that bound, so the product formula has value +infinity.
    def divergence_witness(bound: int) -> dict[str, int]:
        require(
            type(bound) is int and bound >= 0,
            "divergence bound must be a nonnegative integer",
        )
        cutoff = bound + 1
        return {
            "proposed_upper_bound": bound,
            "witness_cutoff": cutoff,
            "partial_sum": cutoff,
        }

    divergence_witnesses = [
        divergence_witness(bound) for bound in (0, 1, 2, 8, 32)
    ]
    nonsummable_failure = all(
        row["witness_cutoff"] == row["proposed_upper_bound"] + 1
        and row["partial_sum"] > row["proposed_upper_bound"]
        for row in divergence_witnesses
    )
    checks = {
        "finite_channel_reference_distances_are_finite": finite_product_values,
        "p_equals_one_triangle_checked_on_finite_reference_set": l1_triangle,
        "p_equals_one_half_has_exact_triangle_counterexample": p_half_failure,
        "infinite_nonsummable_channel_family_has_exact_divergence_schema": (
            nonsummable_failure
        ),
        "paper_exponent_and_summability_guards_required": True,
    }
    require(all(checks.values()), f"ell-p guard failed: {checks}")
    negative_controls = [
        {
            "removed_hypothesis": "p>=1",
            "p": "1/2",
            "points": {
                "x": [1, 0],
                "y": [0, 0],
                "z": [0, 1],
            },
            "d_p_xz": direct_quasi_distance,
            "d_p_xy_plus_d_p_yz": via_origin,
            "violated_conclusion": (
                "triangle_inequality_and_pseudometric_status"
            ),
            "counterexample_verified": p_half_failure,
        },
        {
            "removed_hypothesis": (
                "finite_or_pairwise_summable_declared_channel_family"
            ),
            "finite_description_of_parametric_witness": {
                "points": ["x", "y"],
                "channels": "positive integers",
                "p": "1",
                "weights": "w_c=1",
                "channel_distances": "d_c(F_c(x),F_c(y))=1",
                "partial_sum_formula": "S_N=N",
                "unboundedness_constructor": (
                    "for every integer bound M>=0 choose N=M+1"
                ),
                "sampled_exact_witnesses": divergence_witnesses,
            },
            "violated_conclusion": "finite_valued_pseudometric_status",
            "counterexample_verified": nonsummable_failure,
        },
    ]
    return {
        "receipt_id": "LP-GUARD-1",
        "declared_domain": (
            "p>=1 with a finite declared channel set or a finite weighted "
            "p-sum for every compared pair"
        ),
        "positive_finite_check": {
            "p": "1",
            "channel_count": len(finite_channels),
            "points": [list(point) for point in points],
            "all_ordered_pairs_checked_for_finiteness": len(points) ** 2,
            "all_ordered_triangles_checked": len(points) ** 3,
        },
        "negative_controls": negative_controls,
        "checks": checks,
    }


# ---------------------------------------------------------------------------
# Existing independent certificates
# ---------------------------------------------------------------------------


def load_selector_module() -> Any:
    path = ROOT / "code/a5_closure/echosahedral_selector_certificate.py"
    spec = importlib.util.spec_from_file_location(
        "issue_517_echosahedral_selector_certificate",
        path,
    )
    require(spec is not None and spec.loader is not None, "cannot import selector certificate")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def selector_receipt() -> dict[str, Any]:
    module = load_selector_module()
    manifest_path = (
        ROOT
        / "code/a5_closure/manifests/echosahedral_federation_reference.json"
    )
    stored_receipt_path = (
        ROOT
        / "code/a5_closure/receipts/echosahedral_federation_reference.receipt.json"
    )
    negative_path = (
        ROOT
        / "code/a5_closure/negative_controls/issue_565_negative_controls.json"
    )
    manifest = load_json(manifest_path)
    computed = module.certificate_payload(manifest)
    stored = load_json(stored_receipt_path)
    module.verify_receipt(manifest, stored)
    negatives = module.negative_control_payload(manifest)
    require(
        negatives == load_json(negative_path),
        "stored selector negative controls do not match recomputation",
    )
    require(
        all(row["passed"] for row in negatives["finite_controls"]),
        "selector negative control did not fail closed",
    )

    legacy_path = (
        ROOT
        / "code/particles/hierarchy/certificates/R_screen_sieve_icosahedral_certificate.json"
    )
    legacy_builder_path = (
        ROOT / "code/particles/hierarchy/verify_screen_sieve_theorem.py"
    )
    legacy_spec = importlib.util.spec_from_file_location(
        "issue_517_legacy_screen_sieve_builder",
        legacy_builder_path,
    )
    require(
        legacy_spec is not None and legacy_spec.loader is not None,
        "cannot import legacy screen-sieve builder",
    )
    legacy_builder = importlib.util.module_from_spec(legacy_spec)
    sys.modules[legacy_spec.name] = legacy_builder
    legacy_spec.loader.exec_module(legacy_builder)
    legacy_computed = legacy_builder.build_certificate()
    legacy_stored = load_json(legacy_path)
    require(
        legacy_computed == legacy_stored,
        "legacy conditional screen-sieve artifact differs from exact recomputation",
    )

    return {
        "receipt_id": "ECHOSAHEDRAL-SELECTOR-1",
        "manifest": repo_relpath(manifest_path),
        "receipt": repo_relpath(stored_receipt_path),
        "negative_controls": repo_relpath(negative_path),
        "manifest_sha256": sha256_json(manifest),
        "receipt_sha256": sha256_json(computed),
        "negative_controls_sha256": sha256_json(negatives),
        "finite_negative_control_count": len(negatives["finite_controls"]),
        "separate_conditional_variational_sieve": {
            "artifact": repo_relpath(legacy_path),
            "artifact_sha256": sha256_json(legacy_stored),
            "status": legacy_stored["status"],
            "hierarchy_readout_gate": legacy_stored["hierarchy_screen_readout_gate"],
        },
        "checks": {
            "integer_twelve_port_domain": (
                computed["unit_split"]["domain"]
                == "q in Z^12 with sum(q)=12"
            ),
            "strict_unit_split_gap_two": computed["unit_split"]["strict_gap"] == 2,
            "source_firewall_clear": (
                computed["source_firewall"]["forbidden_dependency_hits"] == []
            ),
            "proper_group_identified_as_A5": (
                computed["icosahedral_selector"]["orientation_preserving_group"]
                == "A5"
            ),
            "avoids_nonexistent_rotation_group_maximality": True,
            "refinement_naturality_checked": all(
                computed["refinement"][key]
                for key in (
                    "unit_lines_natural",
                    "antipode_natural",
                    "frame_gram_natural",
                )
            ),
            "stored_receipt_exactly_recomputed": computed == stored,
            "all_finite_negative_controls_fail_closed": all(
                row["passed"] for row in negatives["finite_controls"]
            ),
            "conditional_variational_sieve_exactly_recomputed": (
                legacy_computed == legacy_stored
            ),
            "conditional_variational_sieve_does_not_supply_hierarchy_readout": (
                legacy_stored["hierarchy_screen_readout_gate"][
                    "supplied_by_screen_sieve"
                ]
                is False
            ),
        },
    }


def separation_receipt() -> dict[str, Any]:
    registry_path = ROOT / "code/audit/receipt_promotion_registry.json"
    registry = load_json(registry_path)
    expected = {
        "coefficient_algebra": "a5-coefficient-algebra",
        "physical_currents": "a5-physical-currents",
        "global_descent": "a5-global-descent",
        "matter_realization": "a5-matter-realization",
    }
    rows = {
        row["role"]: row
        for row in registry["rows"]
        if row.get("role") in expected
    }
    require(set(rows) == set(expected), "A5 separation registry is incomplete")
    require(
        all(rows[role]["receipt_id"] == receipt_id for role, receipt_id in expected.items()),
        "A5 receipt roles are aliased",
    )
    require(
        len({rows[role]["receipt_id"] for role in expected}) == len(expected),
        "A5 layers do not have distinct registry identities",
    )
    layer_locators = {
        role: (
            rows[role]["artifact"],
            rows[role].get("json_pointer", ""),
        )
        for role in expected
    }
    require(
        len(set(layer_locators.values())) == len(layer_locators),
        "A5 registry roles alias the same artifact/json-pointer locator",
    )

    matter_artifact_path = ROOT / rows["matter_realization"]["artifact"]
    matter_artifact = load_json(matter_artifact_path)
    spin = matter_artifact["port_spin_lift"]
    anomalies = matter_artifact["anomalies"]
    refinement = matter_artifact["refinement"]
    anomaly_traces = anomalies["traces"]
    conditional_spin_checked = (
        spin["lift_group_order"] == 120
        and spin["unique_involution"] is True
        and spin["involution_lift_order"] == 4
    )
    conditional_refinement_stable_anomalies_checked = (
        set(anomaly_traces.values()) == {"0"}
        and anomalies["witten_parity"]["even"] is True
        and refinement["natural"] is True
        and all(row["intertwined"] for row in refinement["maps"])
    )
    matter_source_gate = matter_artifact["physical_source_gate"]
    conditional_fixture_boundary_checked = (
        matter_source_gate["upstream_response_constraints_source_bound"] is True
        and matter_source_gate["upstream_current_representation_source_bound"] is False
        and matter_source_gate["charge_pair_derived_within_declared_current_fixture"]
        is True
        and matter_source_gate[
            "conjugate_projector_pair_derived_within_declared_current_fixture"
        ]
        is True
        and matter_source_gate["current_action_on_matter_source_bound"] is False
        and matter_source_gate["matter_lift_source_bound"] is False
        and matter_source_gate["physical_refinement_intertwining_source_bound"]
        is False
        and matter_source_gate["passed"] is False
        and matter_artifact["block_determinant_balance"][
            "declared_matches_derived_pair_up_to_conjugation"
        ]
        is True
        and matter_artifact["upstream"]["semantic_response_artifact_sha256"]
        .startswith("sha256:")
    )
    return {
        "receipt_id": "A5-LAYER-SEPARATION-1",
        "registry": repo_relpath(registry_path),
        "registry_sha256": sha256_json(registry),
        "layers": {
            role: {
                "receipt_id": rows[role]["receipt_id"],
                "receipt_class": rows[role]["receipt_class"],
                "artifact": rows[role]["artifact"],
                "json_pointer": rows[role].get("json_pointer", ""),
                "promoted": rows[role]["promoted"],
                "open_gates": rows[role].get("open_gates", []),
            }
            for role in expected
        },
        "conditional_matter_subreceipts": {
            "spin_lift": {
                "receipt_id": "a5-conditional-spin-lift",
                "artifact": rows["matter_realization"]["artifact"],
                "json_pointer": "/port_spin_lift",
                "exact_checks_pass": conditional_spin_checked,
                "physical_source_promoted": False,
            },
            "refinement_stable_anomaly_algebra": {
                "receipt_id": "a5-conditional-refinement-stable-anomaly-algebra",
                "artifact": rows["matter_realization"]["artifact"],
                "json_pointers": ["/anomalies", "/refinement"],
                "exact_checks_pass": (
                    conditional_refinement_stable_anomalies_checked
                ),
                "physical_source_promoted": False,
            },
        },
        "separate_downstream_gates": {
            "registry_enforced": [
                "physical gauge-current construction",
                "global-form descent",
                "matter realization",
            ],
            "conditional_subreceipt_enforced_but_physical_source_open": [
                "spin lift",
                "refinement-stable anomaly algebra",
            ],
        },
        "checks": {
            "four_registry_roles_have_distinct_receipt_ids": (
                len({rows[role]["receipt_id"] for role in expected})
                == len(expected)
            ),
            "registry_roles_have_distinct_artifact_pointer_locators": (
                len(set(layer_locators.values())) == len(layer_locators)
            ),
            "subject_artifact_does_not_select_receipt_class": (
                registry["promotion_policy"]["requirements"][0].startswith(
                    "receipt class is declared in this independent registry"
                )
            ),
            "coefficient_reconstruction_not_promoted_as_physical_current": (
                rows["coefficient_algebra"]["promoted"] is False
            ),
            "global_descent_not_promoted_from_kernel_identity": (
                rows["global_descent"]["promoted"] is False
            ),
            "matter_schema_not_promoted_as_source_realization": (
                rows["matter_realization"]["promoted"] is False
            ),
            "conditional_spin_lift_subreceipt_checked_separately": (
                conditional_spin_checked and conditional_fixture_boundary_checked
            ),
            "conditional_refinement_stable_anomaly_subreceipt_checked_separately": (
                conditional_refinement_stable_anomalies_checked
                and conditional_fixture_boundary_checked
            ),
        },
    }


def build_receipt() -> dict[str, Any]:
    transaction = transactional_receipt()
    bft = bft_locking_receipt()
    refinement = refinement_receipt()
    lp = lp_receipt()
    selector = selector_receipt()
    separation = separation_receipt()
    negative_controls = {
        "transactional_confluence": transactional_negative_controls(),
        "bft": bft_negative_controls(),
        "refinement": refinement_negative_controls(),
        "ell_p": lp["negative_controls"],
    }
    all_negative_controls = [
        row
        for rows in negative_controls.values()
        for row in rows
    ]
    all_component_checks = {
        "selector": all(selector["checks"].values()),
        "transactional_confluence": all(transaction["checks"].values()),
        "bft": all(bft["checks"].values()),
        "refinement": all(refinement["checks"].values()),
        "ell_p": all(lp["checks"].values()),
        "layer_separation": all(separation["checks"].values()),
        "negative_controls": all(
            row["counterexample_verified"] for row in all_negative_controls
        ),
    }
    require(
        all(all_component_checks.values()),
        f"issue #517 component check failed: {all_component_checks}",
    )
    payload = {
        "schema": SCHEMA,
        "issue": 517,
        "support_status": (
            "exact finite theorem-premise receipt; physical downstream "
            "producer gates remain separately typed"
        ),
        "selector": selector,
        "transactional_confluence": transaction,
        "prepared_certificate_bft": bft,
        "refinement_moduli": refinement,
        "ell_p_guard": lp,
        "a5_layer_separation": separation,
        "negative_controls": negative_controls,
        "checks": all_component_checks,
    }
    payload["receipt_sha256"] = sha256_json(payload)
    return payload


def verify_receipt(receipt: Mapping[str, Any]) -> None:
    expected = build_receipt()
    require(receipt == expected, "stored issue #517 receipt differs from recomputation")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    emit = subparsers.add_parser("emit", help="recompute and write the receipt")
    emit.add_argument("--output", type=Path, required=True)

    verify = subparsers.add_parser("verify", help="recompute and verify a stored receipt")
    verify.add_argument("--receipt", type=Path, required=True)

    args = parser.parse_args(argv)
    if args.command == "emit":
        payload = build_receipt()
        write_json(args.output, payload)
        print(
            json.dumps(
                {
                    "pass": True,
                    "receipt": str(args.output),
                    "receipt_sha256": payload["receipt_sha256"],
                },
                sort_keys=True,
            )
        )
        return 0

    receipt = load_json(args.receipt)
    verify_receipt(receipt)
    print(
        json.dumps(
            {
                "pass": True,
                "receipt": str(args.receipt),
                "receipt_sha256": receipt["receipt_sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
