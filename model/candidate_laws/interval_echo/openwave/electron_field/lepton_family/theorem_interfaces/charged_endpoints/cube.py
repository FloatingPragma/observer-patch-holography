"""Exact charged spin transfers on the archived member-zero cube."""

from dataclasses import dataclass, replace
from itertools import pairwise, product

import numpy as np


CELLS = tuple(product((0, 1), repeat=3))
LINKS = tuple(
    (x, tuple(v + int(j == axis) for j, v in enumerate(x)))
    for x in CELLS
    for axis in range(3)
    if x[axis] == 0
)
CENTER = (0, 0, 0)
A, B, C = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
PARTNERS = ((1, 0, 1), (0, 1, 1))
STATES = tuple(product(range(2), repeat=2))
EXCHANGED = (A, CENTER, C), (B, CENTER, A), (C, CENTER, B)
UNCROSSED = (A, CENTER, C, CENTER, A), (B, CENTER, B)
FIRST_DIRECTIONS = {
    tuple(sorted((u, v))): (u, v) for u, v in ((A, CENTER), (B, CENTER), (CENTER, C))
}
MASKS = tuple(m for m in range(32) if m.bit_count() in (2, 4))
Q6 = MASKS.index(24)
BETA_SEA = sum(1 << (60 * x + r) for x in range(8) for r in range(30, 60))


def mode(x, signed, spin):
    return 60 * CELLS.index(x) + 15 * (spin + (2 if signed < 0 else 0)) + Q6


def c_action(mask, j, create):
    occupied = mask & (1 << j) != 0
    if occupied == create:
        raise ArithmeticError("The controlled mode is Pauli blocked")
    phase = (-1) ** ((mask & ((1 << j) - 1)).bit_count() % 2)
    return mask ^ (1 << j), phase


@dataclass(frozen=True)
class Word:
    occupied: int
    characters: tuple


def addressed_character(word, u, v, signed):
    e = tuple(sorted((u, v)))
    j = LINKS.index(e)
    q = signed if u == e[1] else -signed
    links = list(word.characters)
    links[j] += q
    assert links[j] in (-1, 0, 1)
    return replace(word, characters=tuple(links))


def charge(word):
    return {
        x: 6
        * sum(
            ((word.occupied >> mode(x, q, s)) & 1) - int(q == -1)
            for q, s in product((1, -1), range(2))
        )
        for x in CELLS
    }


def gauss(word):
    result = charge(word)
    for (x, y), e in zip(LINKS, word.characters, strict=True):
        result[x] -= 6 * e
        result[y] += 6 * e
    return result


def prepare(signed, spins):
    word = Word(BETA_SEA, (0,) * len(LINKS))
    amplitude = 1
    for tip, holder, s in zip((A, B), PARTNERS, spins, strict=True):
        if signed > 0:
            operations = (mode(holder, -1, 0), False), (mode(tip, 1, s), True)
        else:
            operations = (mode(tip, -1, s), False), (mode(holder, 1, 0), True)
        for j, create in operations:
            mask, parity = c_action(word.occupied, j, create)
            word = replace(word, occupied=mask)
            amplitude *= parity
        word = addressed_character(word, holder, tip, signed)
        assert not any(gauss(word).values())
    return word, amplitude


def flight(word, routes, signed, spins, unsigned=False):
    """Evaluate the pure basis columns of exp[-iπ P H_link P] at six slots."""
    ownership = {A: 0, B: 1}
    history = []
    amplitude = 1 + 0j
    for leg in routes:
        for u, v in pairwise(leg):
            assert u in ownership and v not in ownership
            bond = tuple(sorted((u, v)))
            first = FIRST_DIRECTIONS[bond]
            link = LINKS.index(bond)
            written_orientation = signed if first[0] == bond[1] else -signed
            assert word.characters[link] == (
                0 if (u, v) == first else written_orientation
            )
            label = ownership.pop(u)
            ownership[v] = label
            s = spins[label]
            before = charge(word)
            if signed == 1:
                operations = (mode(u, 1, s), False), (mode(v, 1, s), True)
            else:
                operations = (mode(u, -1, s), True), (mode(v, -1, s), False)
            local_sign = 1
            for j, create in operations:
                mask, parity = c_action(word.occupied, j, create)
                word = replace(word, occupied=mask)
                local_sign *= parity
            word = addressed_character(word, u, v, signed)
            after = charge(word)
            d = [after[cell] - before[cell] for cell in CELLS]
            assert d == [
                6 * signed * (int(cell == v) - int(cell == u)) for cell in CELLS
            ]
            assert not any(gauss(word).values()) and len(ownership) == 2
            amplitude *= 1j * (1 if unsigned else local_sign)
            history.append(
                {
                    "directed_link": [u, v],
                    "native_charge_transfer_by_cell": d,
                    "occupied_mobile_addresses": sorted(ownership),
                    "link_character_word": word.characters,
                }
            )
    assert len(history) == 6
    return word, amplitude, tuple(ownership[p] for p in (A, B)), history


def phase_bank(matrix, reference):
    symmetric = np.array([1, 0, 0, 0])
    antisymmetric = np.array([0, 1, -1, 0])
    result = {}
    for name, v, norm in (
        ("equal_sharp_spins", symmetric, 1),
        ("antisymmetric_pair", antisymmetric, 2),
        ("unresolved_independent_spin_references", np.eye(4), 4),
    ):
        ports = [(reference + phase * matrix) @ v / 2 for phase in (1, -1)]
        result[name] = [float(np.vdot(port, port).real / norm) for port in ports]
        assert sum(result[name]) == 1
    return result


def exchange(signed):
    by_input = dict(zip(STATES, [prepare(signed, s) for s in STATES], strict=True))
    u = np.zeros((4, 4), dtype=complex)
    v = np.zeros_like(u)
    unsigned = np.zeros_like(u)
    history = None
    for i, s in enumerate(STATES):
        initial, factor = by_input[s]
        swapped = STATES.index(s[::-1])
        final = by_input[s[::-1]]
        assert factor == final[1]
        for routes, target, expected, origin in (
            (EXCHANGED, u, swapped, (1, 0)),
            (UNCROSSED, v, i, (0, 1)),
        ):
            word, amplitude, actual_origins, trace = flight(initial, routes, signed, s)
            assert word == by_input[STATES[expected]][0] and actual_origins == origin
            target[expected, i] = amplitude * factor / by_input[STATES[expected]][1]
            if i == 0 and routes == EXCHANGED:
                history = trace
        alternative, phase, _, _ = flight(initial, EXCHANGED, signed, s, unsigned=True)
        assert alternative == final[0]
        unsigned[swapped, i] = phase
    assert len({w.characters for w, _ in by_input.values()}) == 1
    swap = np.eye(4)[:, [0, 2, 1, 3]]
    # Each reference executes six active positive pulses, contributing i⁶.
    np.testing.assert_array_equal(v, -np.eye(4))
    np.testing.assert_array_equal(u, swap)
    np.testing.assert_array_equal(unsigned, -swap)
    crossed = phase_bank(u, v)
    scalar = phase_bank(unsigned, v)
    assert crossed["equal_sharp_spins"] == [0, 1]
    assert scalar["equal_sharp_spins"] == [1, 0]
    assert crossed["unresolved_independent_spin_references"] == [0.25, 0.75]
    # Rows are the two received spin modes, columns the two untouched refs.
    original_bell_bank = np.eye(4) / 2
    crossed_bell_bank = swap @ original_bell_bank
    received = u @ original_bell_bank
    bell_original = abs(np.vdot(original_bell_bank, received)) ** 2
    bell_crossed = abs(np.vdot(crossed_bell_bank, received)) ** 2
    assert (bell_original, bell_crossed) == (1 / 4, 1)
    return {
        "native_signed_charge": 6 * signed,
        "exchange_relative_to_uncrossed_real_matrix": (v.conj().T @ u)
        .real.astype(int)
        .tolist(),
        "unsigned_relative_to_uncrossed_matrix": (v.conj().T @ unsigned)
        .real.astype(int)
        .tolist(),
        "native_bright_dark": crossed,
        "unsigned_control_bright_dark": scalar,
        "two_spin_reference_joint_Bell_matches_after_exchange": {
            "original_endpoint_assignment": [int(4 * bell_original), 4],
            "crossed_endpoint_assignment": int(bell_crossed),
        },
        "matched_reference_bright_dark_from_self_comparison": [1, 0],
        "outcome_other_than_declared_terminal_records_given_stock": 0,
        "single_hop_receiver_increment_at_half_transfer": 3 * signed,
        "single_hop_incoming_current_at_half_transfer": 3 * signed,
        "exchange_single_codeword_route_log": history,
        "original_fixed_compensator_characters": next(iter(by_input.values()))[
            0
        ].characters,
    }
