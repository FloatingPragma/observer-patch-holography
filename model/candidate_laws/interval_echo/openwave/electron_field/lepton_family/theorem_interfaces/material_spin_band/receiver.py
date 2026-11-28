"""Calibrate a receiver coin through local registered terms and check spatial support."""

from fractions import Fraction
from itertools import pairwise

import numpy as np
import sympy as sp

from model.executable_model import rotations


def receiver_control(board, receiver_edge, receive, sign, strength):
    a, b = receiver_edge
    assert b == receive and board.degree[a] == 5 and board.degree[b] == 6
    edges = set(board.seams)
    tip = next(
        p
        for p in range(12 * (b // 12), 12 * (b // 12 + 1))
        if tuple(sorted((a, p))) in edges and tuple(sorted((b, p))) in edges
    )
    positions = (a, b, tip)
    words = [(p, c) for p in positions for c in range(2)]
    match = int(sign == -1)
    sharp = sp.eye(6)
    code = sharp[:, 2:4]
    bright = (
        sharp[:, words.index((b, match))] + sharp[:, words.index((a, 1 - match))]
    ) / sp.sqrt(2)
    written = bright * bright.H
    assert written * written == written

    def mark(phase):
        return written + phase * (sharp - written)

    def blank(left, right):
        addresses = {left: right, right: left}
        return sp.Matrix.hstack(
            *(sharp[:, words.index((addresses.get(p, p), k))] for p, k in words)
        )

    mark_swap = mark(-1)
    contacts = {"mark": mark_swap, "left": blank(a, tip), "right": blank(b, tip)}
    tape = ("mark", "left", "right", "mark", "left", "mark")
    u_x = sharp
    paths = [[b], [b]]
    for term in tape:
        u_x = contacts[term] * u_x
        for c in range(2):
            where = next(i for i, v in enumerate(u_x * code[:, c]) if v != 0)
            paths[c].append(words[where][0])
    assert all(
        tuple(sorted((p, q))) in edges
        for path in paths
        for p, q in pairwise(path)
        if p != q
    )
    # A full marked pulse contributes a minus to each nonmoving column, and
    # each moving column flips its coin. This parity identity is closed under
    # all products of full registered pulses.
    phase_by_coin = sp.diag(*(1 if c == 0 else -1 for _, c in words))
    assert phase_by_coin * mark_swap * phase_by_coin == -abs_entrywise(mark_swap)
    assert all(
        phase_by_coin * contacts[op] * phase_by_coin == contacts[op]
        for op in ("left", "right")
    )

    # The actual private degree-five reference term at a acts on n_a^2.
    # Divide out a common baseline phase to keep the exact matrices small.
    difference = strength * ((1 + sign) ** 2 - 1)
    time = (
        Fraction(3, 2 * abs(difference))
        if difference < 0
        else Fraction(1, 2 * difference)
    )
    period = Fraction(2, strength)
    inv_time = period - time
    assert 0 < time < period and 0 < inv_time < period
    assert (-difference * time - Fraction(-1, 2)) % 2 == 0
    phase_at_a = sp.diag(*(-sp.I if p == a else 1 for p, _ in words))
    prep = phase_at_a * mark(-sp.I)
    u_z = prep.H * mark(-1) * prep
    u_s = prep.H * mark(-sp.I) * prep
    assert u_x.H * u_x == sharp
    assert (u_z.H * u_z).applyfunc(sp.simplify) == sharp
    assert u_x * code == code * sp.Matrix([[0, 1], [1, 0]])
    z = sp.diag(1 if match == 0 else -1, 1 if match == 1 else -1)
    s = sp.diag(1 if match == 0 else -sp.I, 1 if match == 1 else -sp.I)
    assert (u_z * code).applyfunc(sp.simplify) == code * z
    assert (u_s * code).applyfunc(sp.simplify) == code * s
    assert ((u_x * u_z + u_z * u_x) * code).applyfunc(sp.simplify) == sp.zeros(6, 2)
    assert ((u_x * u_z).H * u_z * u_x * code).applyfunc(sp.simplify) == -code
    y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    assert ((u_s * u_x * u_s.H * code) - (-sign) * code * y).applyfunc(
        sp.simplify
    ) == sp.zeros(6, 2)

    fix_gate = [perm for perm in rotations() if perm[b % 12] == b % 12]
    fix_both = [perm for perm in fix_gate if perm[a % 12] == a % 12]
    assert len(fix_gate) == 5 and len(fix_both) == 1
    return (
        {
            "local_positions": positions,
            "x_full_pulse_edges_in_time_order": [
                receiver_edge if h == "mark" else (a, tip) if h == "left" else (b, tip)
                for h in tape
            ],
            "x_material_address_prefixes_for_two_coin_inputs": paths,
            "z_protected_pulse_tape_with_durations_in_pi_units": [
                [receiver_edge, str(Fraction(1, 2 * strength))],
                [f"port {a} reference", str(time)],
                [receiver_edge, str(Fraction(1, strength))],
                [f"port {a} reference", str(inv_time)],
                [receiver_edge, str(Fraction(3, 2 * strength))],
            ],
            "z_diagonal_in_coin_read_basis": [int(z[i, i]) for i in range(2)],
            "phase_tape_to_calibrate_y": "replace middle marked exposure with half and conjugate stable X",
            "equal_order_bright_probability_given_gateway": 1,
            "x_z_against_z_x_dark_probability_given_gateway": 1,
            "perpendicular_physical_turn_orders_preserving_gate_and_marker": 0,
            "proper_cell_rotations_preserving_exact_gate_and_marker": len(fix_both),
            "both_full_monomial_controls_have_commuting_code_restrictions": True,
        },
        np.asarray(u_x.evalf(), complex),
        np.asarray(u_z.evalf(), complex),
    )


def abs_entrywise(matrix):
    return matrix.applyfunc(abs)
