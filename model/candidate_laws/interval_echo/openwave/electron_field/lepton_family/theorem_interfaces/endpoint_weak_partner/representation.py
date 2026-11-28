"""Exact block and weak generators of the Wilson exterior representation."""

from itertools import combinations

import sympy as sp


MASKS = tuple(m for m in range(32) if m.bit_count() in (2, 4))
CHARGED, NEUTRAL, CONJUGATE = 23, 15, 24
DOUBLET = NEUTRAL, CHARGED
A = tuple(7 ^ (1 << i) for i in range(3))
Q = tuple(tuple((1 << i) | (1 << a) for a in (3, 4)) for i in range(3))


def weights(normalization):
    return {m: sum(normalization[i] for i in range(5) if (1 << i) & m) for m in MASKS}


HYPER = weights((-2, -2, -2, 3, 3))
T3_SIX = weights((0, 0, 0, 3, -3))
PHOTON = {m: HYPER[m] + T3_SIX[m] for m in MASKS}


def unit(row, column):
    return sp.eye(5)[:, row] * sp.eye(5)[column, :]


def exterior(op):
    """Apply one fundamental generator to the literal ascending exterior masks."""
    matrix = sp.zeros(len(MASKS))
    indices = {m: i for i, m in enumerate(MASKS)}
    for m in MASKS:
        for i in range(5):
            if not (1 << i) & m:
                continue
            after = m ^ (1 << i)
            phase = (-1) ** ((m & ((1 << i) - 1)).bit_count())
            for j in range(5):
                if (1 << j) & after or op[j, i] == 0:
                    continue
                signed = phase * (-1) ** ((after & ((1 << j) - 1)).bit_count())
                matrix[indices[after | 1 << j], indices[m]] += signed * op[j, i]
    return matrix


def color():
    pairs = [(unit(i, j), unit(j, i)) for i, j in combinations(range(3), 2)]
    return [
        *(e + f for e, f in pairs),
        *(-sp.I * (e - f) for e, f in pairs),
        unit(0, 0) - unit(1, 1),
        unit(1, 1) - unit(2, 2),
    ]


def check():
    r = unit(3, 4)
    weak = [(r + r.T) / 2, -sp.I * (r - r.T) / 2, (unit(3, 3) - unit(4, 4)) / 2]
    center = sp.diag(-2, -2, -2, 3, 3)
    all_currents = [*color(), *weak, center]
    actual = [exterior(g) for g in all_currents]
    blocks = {}
    for i, m in enumerate(MASKS):
        blocks.setdefault(HYPER[m], []).append(i)
    assert {q: len(ix) for q, ix in blocks.items()} == {-4: 3, -3: 2, 1: 6, 2: 3, 6: 1}
    for q, block in blocks.items():
        p = sp.diag(*(int(i in block) for i in range(len(MASKS))))
        assert all(g * p == p * g for g in actual), q
    col = {
        m: sp.eye(15)[:, MASKS.index(m)] for m in (*DOUBLET, CONJUGATE, *A, *sum(Q, ()))
    }
    assert exterior(r) * col[CHARGED] == col[NEUTRAL]
    assert exterior(r) * col[CONJUGATE] == sp.zeros(15, 1)
    assert all(
        g * col[m] == sp.zeros(15, 1) for g in actual[:8] for m in (*DOUBLET, CONJUGATE)
    )
    weak_casimir = sum((g**2 for g in actual[8:11]), sp.zeros(15))
    for m in DOUBLET:
        assert weak_casimir * col[m] == 3 * col[m] / 4
    assert weak_casimir * col[CONJUGATE] == sp.zeros(15, 1)

    # Form the color-contracted doublet in the actual tensor representation.
    # Matrix multiplication f*Z + Z*g^T is its infinitesimal action.
    h = [
        sum(
            ((-1) ** i * col[A[i]] * col[Q[i][index]].T for i in range(3)), sp.zeros(15)
        )
        for index in range(2)
    ]
    for g in actual:
        for a in range(2):
            expected = sum(
                (
                    g[MASKS.index(DOUBLET[b]), MASKS.index(DOUBLET[a])] * h[b]
                    for b in range(2)
                ),
                sp.zeros(15),
            )
            assert g * h[a] + h[a] * g.T == expected
    # The antisymmetric D⊗D line, times the singlet q=6, is gauge invariant.
    for g in actual:
        weak_block = g.extract(blocks[-3], blocks[-3])
        q6 = g[MASKS.index(CONJUGATE), MASKS.index(CONJUGATE)]
        eps = sp.Matrix([[0, 1], [-1, 0]])
        assert weak_block * eps + eps * weak_block.T + q6 * eps == sp.zeros(2)

    pair_charges = {HYPER[i] - HYPER[j] for i in MASKS for j in MASKS}
    assert 3 not in pair_charges and -3 not in pair_charges
    return {
        "charged_doublet_neutral_partner_and_conjugate_singlet": [
            [m, HYPER[m], T3_SIX[m], PHOTON[m]] for m in (CHARGED, NEUTRAL, CONJUGATE)
        ],
        "native_charge_blocks": [[q, len(ix)] for q, ix in sorted(blocks.items())],
        "all_twelve_generators_preserve_each_block": True,
        "weak_raising_coefficient_in_normalized_internal_basis": 1,
        "weak_raising_annihilates_the_singlet": True,
        "native_color_singlet_composite_weak_doublet": {
            "antitriplet_masks": A,
            "weak_triplet_color_masks_in_two_weak_slots": Q,
            "weak_representation": "dimension two, hypercharge times six is -3",
            "one_creator_one_annihilator_cannot_supply_either_higgs_charge": True,
            "contains_two_native_fermion_defects": True,
        },
        "global_central_phase_exponent_of_bare_negative_charge_conversion": 3,
        "global_central_phase_exponent_of_every_gauge_only_link_operator": 0,
        "converted_finite_code_charge_includes_auxiliary_composite": True,
    }
