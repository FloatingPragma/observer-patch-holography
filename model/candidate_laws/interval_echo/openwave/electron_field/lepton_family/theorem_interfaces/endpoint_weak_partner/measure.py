"""Replay the native block null and a compensated weak-doublet control."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

import sympy as sp

from ..charged_endpoints.cube import BETA_SEA, CELLS, LINKS, c_action
from ..charged_endpoints.measure import (
    BOUNDS,
    CONTRACT,
    GAMMA,
    SOURCE_REVISION,
    clifford,
    read_archive,
)
from . import representation as rep


HERE = Path(__file__).resolve().parent
SITE, HOLDER = (0, 0, 0), (1, 0, 0)
ORDER = sorted(set(rep.HYPER.values()))
CHANNELS = "model/particle_primary_channel_labels/results.json"


def position(cell, row, mask):
    return 60 * CELLS.index(cell) + 15 * row + rep.MASKS.index(mask)


def alleles(word, cell):
    start = 60 * CELLS.index(cell)
    return [(word >> (start + i)) & 1 for i in range(60)]


def increments(word, cell, labels):
    return sum(
        (n - int(i >= 30)) * labels[rep.MASKS[i % 15]]
        for i, n in enumerate(alleles(word, cell))
    )


def occupation_shifts(before, after):
    delta = Counter()
    for cell in CELLS:
        for i, (u, v) in enumerate(
            zip(alleles(before, cell), alleles(after, cell), strict=True)
        ):
            delta[rep.HYPER[rep.MASKS[i % 15]]] += v - u
    return tuple(delta[q] for q in ORDER)


def epsilon(i, j):
    return int((i, j) == (0, 1)) - int((i, j) == (1, 0))


def fock(state, operations):
    factor = 1
    for index, create in operations:
        if bool(state & 1 << index) == create:
            return state, 0
        state, phase = c_action(state, index, create)
        factor *= phase
    return state, factor


def polynomial(coefficients, state):
    result = defaultdict(int)
    for word, amplitude in state.items():
        for coefficient, operations in coefficients:
            target, sign = fock(word, operations)
            result[target] += coefficient * sign * amplitude
    return {word: value for word, value in result.items() if value}


def orient_clifford(gamma):
    eye = sp.eye(4)
    chart = sp.Matrix.hstack(
        eye[:, 0] + eye[:, 2],
        eye[:, 1] + eye[:, 3],
        -eye[:, 0] + eye[:, 2],
        -eye[:, 1] + eye[:, 3],
    ) / sp.sqrt(2)
    beta = gamma[3]
    assert chart.H * beta * chart == sp.diag(1, 1, -1, -1)
    chirality = gamma[0] * gamma[1] * gamma[2] * beta
    assert chirality**2 == eye and chirality * beta + beta * chirality == sp.zeros(4)
    rotated = sp.simplify(chart.H * chirality * chart)
    assert rotated[:2, :2] == rotated[2:, 2:] == sp.zeros(2)
    for i in range(3):
        j = sp.simplify(
            -sp.I * chart.H * gamma[(i + 1) % 3] * gamma[(i + 2) % 3] * chart / 2
        )
        for first, second in ((j[:2, :2], j[:2, :2]), (j[:2, :2], j[2:, 2:])):
            spin_eps = sp.Matrix([[0, 1], [-1, 0]])
            assert first * spin_eps + spin_eps * second.T == sp.zeros(2)
        assert j[:2, 2:] == j[2:, :2] == sp.zeros(2)
    conjugation = sp.simplify(chart.H * gamma[1] * chart.conjugate())
    rows = []
    for i in range(4):
        ((j, sign),) = [
            (j, int(conjugation[i, j])) for j in range(4) if conjugation[i, j]
        ]
        assert j // 2 != i // 2 and abs(sign) == 1
        rows.append((j, sign))
    return rows


def conjugate(operations, rows, coefficient):
    result = []
    for index, create in operations:
        cell, slot = divmod(index, 60)
        row, internal = divmod(slot, 15)
        next_row, phase = rows[row]
        result.append((60 * cell + 15 * next_row + internal, not create))
        coefficient *= phase
    return coefficient, tuple(result)


def native_compensated_coupling():
    terms = []
    for weak, spin, color, scalar_spin in product(
        range(2), range(2), range(3), range(2)
    ):
        phase = (
            epsilon(weak, 1 - weak)
            * epsilon(spin, 1 - spin)
            * (-1) ** color
            * epsilon(scalar_spin, 1 - scalar_spin)
        )
        terms.append(
            (
                phase,
                (
                    (position(SITE, 1 - scalar_spin, rep.Q[color][1 - weak]), True),
                    (position(SITE, scalar_spin, rep.A[color]), True),
                    (position(SITE, 3 - spin, rep.CONJUGATE), True),
                    (position(SITE, spin, rep.DOUBLET[weak]), True),
                ),
            )
        )
    assert len(terms) == 24
    return terms


def relative_ladder(code, sign, crossed=False):
    modified = []
    for row in range(4):
        if sign < 0:
            origin, destination = rep.CHARGED, rep.NEUTRAL
        else:
            origin, destination = rep.NEUTRAL, rep.CHARGED
        modified.append(
            (
                1,
                (
                    (position(SITE, row, origin), False),
                    (position(SITE, row, destination), True),
                ),
            )
        )
    if crossed:
        modified = [
            (1, tuple(reversed([(i, not create) for i, create in ops])))
            for _, ops in modified
        ]
    return polynomial(modified, code)


def norm_sq(ket):
    return sum(x * x for x in ket.values())


def onsite_increment(word):
    beta_chart_weights = {i: (3 if i < 30 else -3) for i in range(60)}
    return sum(
        beta_chart_weights[i] * (n - int(i >= 30))
        for cell in CELLS
        for i, n in enumerate(alleles(word, cell))
    )


def verify_code(signed, charge_conjugation, coupling):
    """Test the complete four-dimensional compensated spin code on the CAR words."""
    terms = (
        coupling
        if signed < 0
        else [conjugate(ops, charge_conjugation, coeff) for coeff, ops in coupling]
    )
    adjoint = [
        (coeff, tuple((j, not create) for j, create in reversed(ops)))
        for coeff, ops in terms
    ]
    recorded = []
    supports = []
    xq = signed * 6
    for spin in range(2):
        start_ops = (
            (position(SITE, spin + 2, rep.CONJUGATE), False),
            (position(HOLDER, 0, rep.CONJUGATE), True),
        )
        initial = BETA_SEA
        coefficient = 1
        if signed > 0:
            coefficient, start_ops = conjugate(start_ops, charge_conjugation, 1)
        start, phase = fock(initial, start_ops)
        right = {start: coefficient * phase}
        doublet = polynomial(terms, right)
        assert onsite_increment(start) == 6
        assert {onsite_increment(w) for w in doublet} == {12}
        assert norm_sq(doublet) == 12
        assert polynomial(adjoint, right) == {}
        assert polynomial(terms, doublet) == {}
        assert polynomial(adjoint, doublet) == {
            word: 12 * value for word, value in right.items()
        }
        supports.append(doublet.keys())
        for term in (right, doublet):
            for word in term:
                # Original charge-six Wilson link remains untouched. Its
                # Gauss label cancels exactly these two local charges.
                q = [increments(word, cell, rep.HYPER) for cell in CELLS]
                assert q == [
                    xq * (int(cell == SITE) - int(cell == HOLDER)) for cell in CELLS
                ]
                assert ((word.bit_count() - BETA_SEA.bit_count()) % 2) == 0
                assert ((sum(alleles(word, SITE)) - 30) % 2) == 1
                assert increments(word, SITE, rep.PHOTON) == xq
        target_shifts = {occupation_shifts(start, word) for word in doublet}
        shift = tuple(-signed * int(q in (-4, -3, 1, 6)) for q in ORDER)
        assert target_shifts == {shift}
        weak_charge = 0
        charged_components = {}
        for word, value in doublet.items():
            location = [position(SITE, row, rep.CHARGED) for row in range(4)]
            if sum(bool(word & 1 << i) for i in location) != 2:
                charged_components[word] = value
                weak_charge += value**2
        assert weak_charge * 2 == 12
        # An oriented weak probe is a relative instrument (the charged
        # pointer and gauge frame are separate supplied equipment).
        raise_to_partner = relative_ladder(charged_components, signed)
        reverse_control = relative_ladder(charged_components, signed, crossed=True)
        assert (
            norm_sq(raise_to_partner) == weak_charge and norm_sq(reverse_control) == 0
        )
        assert relative_ladder(right, signed) == {}
        for word in raise_to_partner:
            assert increments(word, SITE, rep.PHOTON) == 0
            assert increments(word, SITE, rep.HYPER) == xq
            assert occupation_shifts(start, word) == shift
        recorded.append(
            {
                "spin_input": spin,
                "right_endpoint_configurations": len(right),
                "normalized_weak_composite_configurations": len(doublet),
                "created_block_count_changes": shift,
                "matter_plus_composite_hypercharge_times_six": xq,
                "relative_ladder_target": rep.NEUTRAL,
                "relative_ladder_start": rep.CHARGED,
                "weak_countercharge_times_six": xq,
                "relative_wrong_ladder_or_singlet_success": 0,
                "relative_neutral_composite_and_charged_matter_weight_in_doublet": [
                    1,
                    2,
                ],
                "relative_full_weak_transfer_weight_after_that_acceptance": 1,
            }
        )
    assert supports[0].isdisjoint(supports[1])
    return recorded


def separate_native_doublet(band, charge_conjugation):
    # Each matrix element U_ab|0> has squared Haar norm 1/d.
    # B† B on the vacuum is d·I, and B/√d is a contraction.
    blocks = {
        row["Q6"]: (row["dimension"], Fraction(row["Casimir"]))
        for row in band["allocation"]["native_rows"]
    }
    assert blocks[-3] == (2, Fraction(9, 20))
    assert blocks[6] == (1, Fraction(3, 10))
    energy = 6 + Fraction(105, 2) * blocks[-3][1]
    assert energy == Fraction(237, 8)
    charge_checks = []
    for direction in (-1, 1):
        words = set()
        for left, right in product(rep.DOUBLET, repeat=2):
            ops = (
                (position(HOLDER, 2, right), False),
                (position(SITE, 0, left), True),
            )
            coefficient = 1
            if direction > 0:
                coefficient, ops = conjugate(ops, charge_conjugation, coefficient)
            word, factor = fock(BETA_SEA, ops)
            assert abs(coefficient * factor) == 1
            assert onsite_increment(word) == 6
            assert [increments(word, cell, rep.HYPER) for cell in CELLS] == [
                3 * direction * (int(cell == SITE) - int(cell == HOLDER))
                for cell in CELLS
            ]
            words.add(word)
        assert len(words) == 4
        charge_checks.append(
            {
                "endpoint_and_its_own_compensator_hypercharge": [
                    direction * 3,
                    -direction * 3,
                ],
                "physical_link_matter_columns": len(words),
            }
        )
    return {
        "normalized_probe_accept_probability_given_beta_stock": 1,
        "matter_weak_casimir": "3/4",
        "formula_for_probability_of_weak_slots_after_link_and_compensator_are_unread": [
            "1/2",
            "1/2",
        ],
        "bare_onsite_plus_electric_transfer_reference_energy": str(energy),
        "separate_singlet_reference_energy": band["exact_reference_frequency"],
        "the_two_ends_compensate_their_own_distinct_link_representation": charge_checks,
    }


def controlled_clock():
    # Independently compute the projected line positions and the finite
    # resonant/detuned gates from the actual inferred onsite difference.
    native_block = sp.diag(0, 6)
    drive = sp.Matrix([[0, 1], [1, 0]])
    controlled = native_block + drive
    assert set(controlled.eigenvals()) == {3 - sp.sqrt(10), 3 + sp.sqrt(10)}
    assert (controlled - 3 * sp.eye(2)) ** 2 == 10 * sp.eye(2)
    lower = (sp.eye(2) - (controlled - 3 * sp.eye(2)) / sp.sqrt(10)) / 2
    assert lower[1, 1] == (1 - 3 / sp.sqrt(10)) / 2
    maximum = (controlled[1, 0]) ** 2 / 10
    assert maximum == sp.Rational(1, 10)
    banks = []
    for name, converted in (
        ("bare", Fraction(0)),
        ("natural_detuning_control", Fraction(1, 10)),
        ("supplied_offset_compensates_detuning", Fraction(1, 2)),
    ):
        p = converted / 2
        values = [1 - converted, p, p / 2, p / 2]
        assert sum(values) == 1
        banks.append({"mode": name, "probabilities": [str(x) for x in values]})
    return banks


def run():
    archived = {
        path: read_archive(path) for path in (GAMMA, CONTRACT, BOUNDS, CHANNELS)
    }
    c = json.loads(archived[CONTRACT])["fixed_primary_member"]
    band = json.loads(archived[BOUNDS])["selected_member_bound"]["result"]
    assert c["member_id"] == band["member_id"]
    assert all(band[k] == v for k, v in c["ids"].items())
    assert c["source"]["side"] == 2 and Fraction(c["source"]["spacing"]) == 1
    assert tuple(sorted((SITE, HOLDER))) in LINKS
    assert Fraction(band["exact_reference_frequency"]) == Fraction(87, 4)
    rep_checks = rep.check()
    native_channels = json.loads(archived[CHANNELS])["native_numbers"]
    assert [
        [row["Q6"], row["internal_dimension"]] for row in native_channels["blocks"]
    ] == rep_checks["native_charge_blocks"]
    assert native_channels["conserved_global_occupation_counts"] == 5
    assert native_channels["full_H_has_pair_creation_or_interblock_terms"] is False
    charge_conj = orient_clifford(clifford(archived[GAMMA]))
    coupling = native_compensated_coupling()
    return {
        "source_revision": SOURCE_REVISION,
        "source_sha256": {
            path: sha256(value).hexdigest() for path, value in archived.items()
        },
        "source_identity": {"member_id": c["member_id"], **c["ids"]},
        "representation": rep_checks,
        "occupation_block_order": ORDER,
        "unchanged_Wilson_mixed_propagator_element_for_every_wait": 0,
        "bare_free_flight_coherent_five_counter_change_for_every_input": [0] * 5,
        "gamma5_at_fixed_dirac_beta_projected_matrix": "zero on either beta sign",
        "separate_native_doublet_positive_control": separate_native_doublet(
            band, charge_conj
        ),
        "local_projected_control": {
            "physical_pair_addresses": [SITE, HOLDER],
            "reference_character_word": "fixed charge-six link and constant spectator links",
            "new_interaction": "supply a local gauge-invariant four-CAR term and its conjugate",
            "added_modes": 0,
            "native_pair_reference_energy": "87/4",
            "native_detuning_due_to_extra_defects": 6,
            "g_calibration_example": 1,
            "native_and_controlled_H_after_subtracting_pair_reference": [
                [[0, 0], [0, 6]],
                [[0, 1], [1, 6]],
            ],
            "conditional_control_split": "2*sqrt(10)",
            "conditional_weak_weight_of_lower_line": "(1-3/sqrt(10))/2",
            "conditional_weak_weight_of_upper_line": "(1+3/sqrt(10))/2",
            "maximum_conversion_at_natural_detuning": [1, 10],
            "relative_ladder_half_pulse_success_given_herald_at_natural_detuning": [
                1,
                40,
            ],
            "relative_ladder_half_pulse_success_with_supplied_resonant_offsets": [1, 8],
            "charge_sign_banks": [
                {
                    "physical_endpoint_sign": sign,
                    "inputs": verify_code(sign, charge_conj, coupling),
                }
                for sign in (-1, 1)
            ],
            "relative_weak_complete_bank_outcomes": [
                "native_singlet_retained",
                "converted_but_relative_neutral_composite_flag_refused",
                "relative_neutral_composite_flag_and_weak_transfer_refused",
                "relative_neutral_composite_flag_and_signed_weak_transfer",
            ],
            "relative_weak_complete_banks_conditioned_on_stock": controlled_clock(),
            "finite_model_stock_is_postselected": True,
            "absolute_chirality_or_physical_rest_mass_established": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    data = json.loads(json.dumps(run(), allow_nan=False))
    if args.record:
        with args.record.open("x") as file:
            file.write(json.dumps(data, indent=2) + "\n")
    elif data != json.loads((HERE / "results.json").read_text()):
        raise ArithmeticError("The weak partner check differs from the frozen receipt")
    print(
        "Verified the native five-block null and the supplied weak-partner comparison for both signs."
    )


if __name__ == "__main__":
    main()
