"""Actual charged hop, two local spin memories and their coincidence probabilities."""

from functools import cache
from itertools import product

import sympy as sp

from ..charged_endpoints.cube import c_action
from ..charged_endpoints.measure import GAMMA, clifford, read_archive


@cache
def wilson_columns():
    gammas = clifford(read_archive(GAMMA))
    beta = gammas[3]
    e = sp.eye(4)
    frame = sp.Matrix.hstack(
        e[:, 0] + e[:, 2], e[:, 1] + e[:, 3], e[:, 2] - e[:, 0], e[:, 3] - e[:, 1]
    ) / sp.sqrt(2)
    assert frame.H * beta * frame == sp.diag(1, 1, -1, -1)
    conjugation = sp.simplify(frame.T * gammas[1] * frame)
    assert conjugation == sp.Matrix(
        [[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]]
    )
    for gamma in gammas[:3]:
        hop = (-beta + beta * gamma) / 2
        for oriented in (hop, hop.H):
            local = sp.simplify(frame.H * oriented * frame)
            assert local[:2, :2] == -sp.eye(2) / 2
            assert local[2:, 2:] == sp.eye(2) / 2
            assert all(
                sp.simplify(local[:, s].dot(local[:, s], hermitian=True))
                == sp.Rational(1, 2)
                for s in range(4)
            )
    return conjugation


def move(word, source, target):
    removed, s = c_action(word, source, False)
    added, t = c_action(removed, target, True)
    return added, s * t


def populations(word):
    # Each two-mode group is the actual q6 beta doublet at holder, origin, receiver.
    return tuple(((word >> (2 * i)) & 3).bit_count() for i in range(3))


def check_conjugate_source_map():
    # In the plus chart the three measured beta groups are negative, positive,
    # positive. The paired minus chart uses their opposite beta groups.
    m = wilson_columns()
    conjugate = sp.zeros(64)
    for word in range(64):
        target, phase = 63, 1
        for slot in reversed(range(6)):
            if word & (1 << slot):
                site, s = divmod(slot, 2)
                beta_offset = 2 if site == 0 else 0
                native = m[beta_offset + s, (2 - beta_offset) + (1 - s)]
                target, fermion_sign = c_action(target, 2 * site + 1 - s, False)
                phase *= native * fermion_sign
        conjugate[target, word] = phase
    assert conjugate.H * conjugate == sp.eye(64)
    plus_in, minus_in = 0b000011, 0b111100
    plus_out, minus_out = 0b010001, 0b011101
    positive, negative = sp.zeros(64), sp.zeros(64)
    positive[plus_out, plus_in], negative[minus_out, minus_in] = 1, 1
    transported = conjugate * positive * conjugate.H
    assert transported in (negative, -negative)
    assert conjugate * (positive.H * positive) * conjugate.H == negative.H * negative


def spin_records(sign, *, mark=False):
    wilson_columns()
    if sign not in (-1, 1):
        raise ValueError("Expected a signed transfer")
    initial = 0b000011 if sign == 1 else 0b111100
    records = sp.zeros(4, 1)
    reset_words = set()
    for s in (0, 1):
        pair_op = (s, 2 + s) if sign == 1 else (2 + s, s)
        first, pair_phase = move(initial, *pair_op)
        created = tuple(
            a - b for a, b in zip(populations(first), populations(initial), strict=True)
        )
        assert created == (-sign, sign, 0)
        origin_code = (first >> 2) & 3
        assert origin_code in (1, 2)
        mark_phase = -1 if mark and origin_code == 2 else 1
        hop_op = (2 + s, 4 + s) if sign == 1 else (4 + s, 2 + s)
        received, hop_phase = move(first, *hop_op)
        returned, reverse_phase = move(received, *reversed(hop_op))
        assert returned == first and reverse_phase * hop_phase == 1
        end_change = tuple(
            a - b
            for a, b in zip(populations(received), populations(first), strict=True)
        )
        assert end_change == (0, -sign, sign)
        # Each local neutral memory swaps its blank spin label for the resident
        # spin in occupation one. Other local occupation sectors carry a refusal.
        assert (received & 3) in (1, 2) and (received >> 4) in (1, 2)
        holder = (received & 3).bit_length() - 1
        receiver = (received >> 4).bit_length() - 1
        reset_words.add((received & 0b001100) | 0b010001)
        # exp(-i*pi*H_hop) in its two-word block has phase +i for a
        # positive-beta excess and -i for a negative-beta vacancy.
        records[2 * holder + receiver] += (
            sp.I * sign * hop_phase * pair_phase * mark_phase / sp.sqrt(2)
        )
    assert reset_words == {0b010001 if sign == 1 else 0b011101}
    output = records * records.H
    singlet = sp.Matrix([0, 1, -1, 0]) / sp.sqrt(2)
    if mark:
        singlet = sp.kronecker_product(sp.eye(2), sp.diag(1, -1)) * singlet
    assert output == singlet * singlet.H
    assert populations(initial)[2] in (
        0,
        2,
    )  # shuttering the hop gives no receiver spin
    return output


def coincidences(density):
    z = sp.diag(1, -1)
    x = sp.Matrix([[0, 1], [1, 0]])
    table = []
    for left, right in product((z, x), ((z + x) / sp.sqrt(2), (z - x) / sp.sqrt(2))):
        probabilities = {}
        for a, b in product((1, -1), repeat=2):
            effect = (
                sp.kronecker_product(sp.eye(2) + a * left, sp.eye(2) + b * right) / 4
            )
            probabilities[a, b] = sp.simplify((effect * density).trace())
        assert sp.simplify(sum(probabilities.values())) == 1
        table.append(sp.simplify(sum(a * b * p for (a, b), p in probabilities.items())))
    return table, sp.simplify(abs(table[0] + table[1] + table[2] - table[3]))


def receiver_phase_ports():
    # A normalized hypercharge pulse on the receiver multiplies its active
    # occupation word by i**N. All unaddressed charge phases cancel on K.
    pulse = sp.diag(*(sp.I ** populations(word)[2] for word in range(64)))
    sine_plus = (sp.eye(2) + sp.Matrix([[0, -sp.I], [sp.I, 0]])) / 2
    readings = []
    for sign, source, output in ((1, 0b000011, 0b010001), (-1, 0b111100, 0b011101)):
        accepted = sp.zeros(64)
        accepted[output, source] = 1
        shifted = pulse * accepted * pulse.H
        assert shifted == sp.I**sign * accepted
        probe = sp.Matrix([1, shifted[output, source]]) / sp.sqrt(2)
        probability = sp.simplify((probe.H * sine_plus * probe)[0])
        assert probability == (1 + sign) // 2
        readings.append(int(probability))
    assert readings == [1, 0]
    return readings


def check():
    check_conjugate_source_map()
    positive, negative = (spin_records(sign) for sign in (1, -1))
    assert positive == negative
    correlations, combined = coincidences(positive)
    dephased = sp.diag(*positive.diagonal())
    controls, classical = coincidences(dephased)
    assert combined == 2 * sp.sqrt(2) and classical == sp.sqrt(2)
    tagged = spin_records(1, mark=True)
    assert tagged == spin_records(-1, mark=True)
    x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    xx, zz = sp.kronecker_product(x, x), sp.kronecker_product(z, z)
    assert [(xx * rho).trace() for rho in (positive, tagged)] == [-1, 1]
    assert [(zz * rho).trace() for rho in (positive, tagged)] == [-1, -1]
    return {
        "held_correlations": [str(x) for x in correlations],
        "combined": str(combined),
        "dephased_correlations": [str(x) for x in controls],
        "dephased_combined": str(classical),
        "origin_mark_off_on_XX": [-1, 1],
        "origin_mark_off_on_ZZ": [-1, -1],
        "receiver_current_sine_plus_for_positive_and_negative": receiver_phase_ports(),
    }
