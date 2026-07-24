#!/usr/bin/env python3
"""Export the verified rooted-tree packet-net repair domain."""

from __future__ import annotations

import argparse
import json
from collections import deque
from datetime import datetime, timezone
from fractions import Fraction
from math import sqrt
from itertools import product
from pathlib import Path
from typing import Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "consensus" / "runs" / "verified_tree_packet_net_domain.json"

ROOT_VERTEX = "r"
PARENTS = {"a": "r", "b": "a", "c": "a"}
CHILDREN = {"r": ["a"], "a": ["b", "c"], "b": [], "c": []}
VERTICES = [ROOT_VERTEX, "a", "b", "c"]
EDGES = [("r", "a"), ("a", "b"), ("a", "c")]
ALPHABET = [0, 1, 2]
GAUGE_LABELS = [0, 1]
WEIGHTS = {("r", "a"): 5, ("a", "b"): 1, ("a", "c"): 1}

DEFAULT_REFERENCE_SIGMA_B = {
    0: Fraction(1, 3),
    1: Fraction(1, 3),
    2: Fraction(1, 3),
}
DEFAULT_PETZ_CONDITIONAL = {
    0: {0: Fraction(1, 2), 1: Fraction(1, 4), 2: Fraction(1, 4)},
    1: {0: Fraction(1, 4), 1: Fraction(1, 2), 2: Fraction(1, 4)},
    2: {0: Fraction(1, 4), 1: Fraction(1, 4), 2: Fraction(1, 2)},
}


PacketState = tuple[tuple[int, ...], tuple[int, ...]]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _edge_key(edge: tuple[str, str]) -> str:
    return f"{edge[0]}-{edge[1]}"


def all_packet_states() -> Iterable[PacketState]:
    for packets in product(ALPHABET, repeat=len(VERTICES)):
        for hidden in product(GAUGE_LABELS, repeat=len(VERTICES)):
            yield packets, hidden


def _idx(vertex: str) -> int:
    return VERTICES.index(vertex)


def inconsistent_edges(state: PacketState) -> list[tuple[str, str]]:
    packets, _hidden = state
    return [(u, v) for u, v in EDGES if packets[_idx(u)] != packets[_idx(v)]]


def potential(state: PacketState) -> int:
    return sum(WEIGHTS[edge] for edge in inconsistent_edges(state))


def enabled_sites(state: PacketState) -> list[str]:
    packets, _hidden = state
    return [
        vertex
        for vertex in VERTICES
        if vertex != ROOT_VERTEX and packets[_idx(vertex)] != packets[_idx(PARENTS[vertex])]
    ]


def repair(state: PacketState, site: str) -> PacketState:
    if site == ROOT_VERTEX:
        return state
    packets, hidden = state
    parent = PARENTS[site]
    updated = list(packets)
    updated[_idx(site)] = packets[_idx(parent)]
    return tuple(updated), hidden


def terminal_states_from(state: PacketState) -> set[PacketState]:
    seen = {state}
    queue = deque([state])
    terminals: set[PacketState] = set()
    while queue:
        current = queue.popleft()
        sites = enabled_sites(current)
        if not sites:
            terminals.add(current)
            continue
        for site in sites:
            nxt = repair(current, site)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return terminals


def tree_checks() -> dict[str, object]:
    total_states = 0
    max_terminals = 0
    max_steps_seen = 0
    for state in all_packet_states():
        total_states += 1
        sites = enabled_sites(state)
        is_consistent = not inconsistent_edges(state)
        if is_consistent != (not sites):
            raise AssertionError(f"repair completeness failed for {state}")
        for site in sites:
            nxt = repair(state, site)
            if potential(nxt) >= potential(state):
                raise AssertionError(f"Lyapunov descent failed for {state} at {site}")
        terminals = terminal_states_from(state)
        max_terminals = max(max_terminals, len(terminals))
        if len(terminals) != 1:
            raise AssertionError(f"nonunique terminal states for {state}: {terminals}")
        terminal = next(iter(terminals))
        if inconsistent_edges(terminal):
            raise AssertionError(f"terminal state is inconsistent: {terminal}")
        max_steps_seen = max(max_steps_seen, potential(state))
    return {
        "total_states_checked": total_states,
        "repair_completeness": True,
        "strict_lyapunov_descent": True,
        "unique_terminal_normal_form": True,
        "max_terminal_count_seen": max_terminals,
        "potential_step_bound": max_steps_seen,
    }


def _as_fraction(value: object) -> Fraction:
    """Parse fixture numbers without introducing binary floating-point error."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    return Fraction(str(value))


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _fraction_matrix_text(matrix: list[list[Fraction]]) -> list[list[str]]:
    return [[_fraction_text(value) for value in row] for row in matrix]


def compute_petz_recovery_witness(
    *,
    reference_sigma_b: Mapping[int, object] | None = None,
    conditional: Mapping[int, Mapping[int, object]] | None = None,
) -> dict[str, object]:
    """Recompute the explicit classical Petz recovery-channel certificate.

    The input algebra is the diagonal three-point algebra on ``B``.  The
    output algebra is the diagonal nine-point algebra on ``B x D``.  In the
    lexicographic bases below the recovery matrix is

        R[(b_out,d), b_in] = 1[b_out=b_in] sigma(d|b_in).

    For this commutative channel its Choi matrix is diagonal, with the entries
    of this 9x3 stochastic matrix on that diagonal.  Complete positivity is
    therefore checked by entrywise nonnegativity, trace preservation by exact
    rational column sums, and trace-norm contractivity by the induced l1 norm.
    The marginal ``Tr_D o R`` is independently evaluated on every basis input.

    No boolean supplied by the subject is trusted.  The optional arguments
    exist so mutation tests can perturb one defining antecedent at a time.
    """
    raw_reference = (
        DEFAULT_REFERENCE_SIGMA_B
        if reference_sigma_b is None
        else reference_sigma_b
    )
    raw_conditional = (
        DEFAULT_PETZ_CONDITIONAL
        if conditional is None
        else conditional
    )

    reference = {
        b: _as_fraction(raw_reference.get(b, 0))
        for b in ALPHABET
    }
    kernel = {
        b: {
            d: _as_fraction(raw_conditional.get(b, {}).get(d, 0))
            for d in ALPHABET
        }
        for b in ALPHABET
    }

    output_basis = [(b, d) for b in ALPHABET for d in ALPHABET]
    recovery_matrix: list[list[Fraction]] = [
        [
            kernel[b_in][d] if b_out == b_in else Fraction(0)
            for b_in in ALPHABET
        ]
        for b_out, d in output_basis
    ]

    # For the commutative measure-and-prepare extension, these matrix entries
    # are exactly the diagonal of the Choi matrix (including structural zeros).
    choi_diagonal = [
        recovery_matrix[row][column]
        for column in range(len(ALPHABET))
        for row in range(len(output_basis))
    ]
    choi_nonnegative = all(value >= 0 for value in choi_diagonal)

    column_sums = [
        sum((row[column] for row in recovery_matrix), Fraction(0))
        for column in range(len(ALPHABET))
    ]
    column_residuals = [value - 1 for value in column_sums]
    trace_preserving = all(value == 0 for value in column_residuals)

    marginal_matrix: list[list[Fraction]] = []
    for b_out in ALPHABET:
        marginal_matrix.append(
            [
                sum(
                    (
                        recovery_matrix[output_basis.index((b_out, d))][b_in]
                        for d in ALPHABET
                    ),
                    Fraction(0),
                )
                for b_in in ALPHABET
            ]
        )
    identity_matrix = [
        [Fraction(int(row == column)) for column in ALPHABET]
        for row in ALPHABET
    ]
    basis_marginal_recovery = marginal_matrix == identity_matrix

    absolute_column_sums = [
        sum((abs(row[column]) for row in recovery_matrix), Fraction(0))
        for column in range(len(ALPHABET))
    ]
    induced_l1_norm = max(absolute_column_sums)
    trace_norm_contractive = choi_nonnegative and trace_preserving and induced_l1_norm <= 1
    cptp = choi_nonnegative and trace_preserving

    reference_sum = sum(reference.values(), Fraction(0))
    reference_normalized = reference_sum == 1
    support_gap = min(reference.values())
    full_support = support_gap > 0

    failures: list[str] = []
    if not choi_nonnegative:
        failures.append("classical_choi_diagonal_has_negative_entry")
    if not trace_preserving:
        failures.append("recovery_column_normalization_failed")
    if not basis_marginal_recovery:
        failures.append("basis_marginal_recovery_failed")
    if induced_l1_norm > 1:
        failures.append("induced_l1_norm_exceeds_one")
    if not reference_normalized:
        failures.append("reference_state_not_normalized")
    if not full_support:
        failures.append("reference_state_not_full_support")

    witness_pass = (
        cptp
        and trace_norm_contractive
        and basis_marginal_recovery
        and reference_normalized
        and full_support
    )
    return {
        "channel": "classical_full_support_petz_for_trace_D",
        "formula": "R(mu_B)(b,d)=mu_B(b) sigma(d|b)",
        "witness_kind": "independently_computed_classical_recovery_channel",
        "audit_issue": 518,
        "arithmetic": "fractions.Fraction exact rational arithmetic",
        "input_basis": [f"b={b}" for b in ALPHABET],
        "output_basis": [f"(b={b},d={d})" for b, d in output_basis],
        "recovery_matrix_shape": [len(output_basis), len(ALPHABET)],
        "recovery_matrix": _fraction_matrix_text(recovery_matrix),
        "classical_choi_diagonal": [_fraction_text(value) for value in choi_diagonal],
        "classical_choi_diagonal_nonnegative": choi_nonnegative,
        "minimum_choi_diagonal": _fraction_text(min(choi_diagonal)),
        "column_sums_exact": {
            str(b): _fraction_text(column_sums[index])
            for index, b in enumerate(ALPHABET)
        },
        "column_residuals_exact": {
            str(b): _fraction_text(column_residuals[index])
            for index, b in enumerate(ALPHABET)
        },
        "trace_preserving": trace_preserving,
        "basis_marginal_recovery_matrix": _fraction_matrix_text(marginal_matrix),
        "basis_marginal_recovery_identity": basis_marginal_recovery,
        "absolute_column_sums_exact": {
            str(b): _fraction_text(absolute_column_sums[index])
            for index, b in enumerate(ALPHABET)
        },
        "induced_l1_norm": _fraction_text(induced_l1_norm),
        "cptp": cptp,
        "trace_norm_contractive": trace_norm_contractive,
        "reference_state_normalized": reference_normalized,
        "support_gap_gamma_sigma": float(support_gap),
        "support_gap_exact": _fraction_text(support_gap),
        "full_support": full_support,
        "inverse_sqrt_bound": 1.0 / sqrt(float(support_gap)) if full_support else None,
        "input_domain": "all diagonal B-packet distributions because gamma_sigma > 0",
        "support_obstruction": "sigma_B(b)=0 with input mass at b makes the inverse sigma_B^{-1/2} undefined unless the domain is restricted",
        "reference_sigma_B": {
            str(b): _fraction_text(value)
            for b, value in reference.items()
        },
        "sample_conditional": {
            str(b): {
                str(d): _fraction_text(kernel[b][d])
                for d in ALPHABET
            }
            for b in ALPHABET
        },
        "sample_column_sums": {
            str(b): float(column_sums[index])
            for index, b in enumerate(ALPHABET)
        },
        "failures": failures,
        "pass": witness_pass,
        "claim_boundary": {
            "closed_here": (
                "CPTP, basis marginal recovery, full-support domain, and induced "
                "l1 contractivity for this explicit classical diagonal channel"
            ),
            "not_closed_here": [
                "a noncommutative Petz recovery theorem for arbitrary packet algebras",
                "a physical quantum channel or hardware realization",
                "the remaining receipt-class, Thomson, hierarchy, dark-sector, or Planck-in/out work in issue #518",
            ],
        },
    }


def petz_checks() -> dict[str, object]:
    witness = compute_petz_recovery_witness()
    if not witness["pass"]:
        raise AssertionError(f"Petz recovery witness failed: {witness['failures']}")
    return witness


def build_payload() -> dict[str, object]:
    return {
        "artifact": "oph_verified_tree_packet_net_domain",
        "object_id": "ConsensusTreePacketRepairDomain_Issue238",
        "generated_utc": _timestamp(),
        "issue": 238,
        "role": (
            "Nontrivial exported packet-net fixture with exhaustively checked repair "
            "mechanics, an exact classical diagonal Petz recovery witness, and "
            "quotient compatibility."
        ),
        "domain": {
            "graph": {
                "root": ROOT_VERTEX,
                "vertices": VERTICES,
                "oriented_edges": [_edge_key(edge) for edge in EDGES],
                "parents": PARENTS,
            },
            "packet_alphabet": ALPHABET,
            "hidden_gauge_labels_per_vertex": GAUGE_LABELS,
            "state_space": "S_i = Z_3 x Z_2; interface projections read only the Z_3 packet",
            "repair": "for non-root i, if x_i != x_parent(i), set x_i := x_parent(i) and leave hidden labels fixed",
            "weights": {_edge_key(edge): weight for edge, weight in WEIGHTS.items()},
        },
        "theorem_checks": tree_checks(),
        "petz_domain": petz_checks(),
        "quotient_compatibility": {
            "gauge_action": "independent permutations of the hidden Z_2 labels at each vertex",
            "reason": "interface projections, potentials, enabled sites, and repairs depend only on packet labels x_i",
            "descends_to_quotient": True,
            "physical_law_use": "gauge-invariant observables factor through the quotient normal form on this verified branch",
        },
        "general_obstructions_isolated": [
            "nonzero cycle holonomy on non-tree affine packet nets makes strict global consistency impossible",
            "packet-closure failure prevents a microscopic instrument from pushing forward to autonomous packet dynamics",
            "Petz recovery without full support requires an explicit support-domain restriction",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Export the verified tree packet-net domain artifact.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
