"""Independent adversarial oracles for the charge and spin candidate modules.

Charge uses a Laurent-polynomial convolution. Spin uses vector tensor spinor
generators, and exchange uses permutation inversions of ordered wedges.
These differ from the producer's occupation census, gauge-kernel reduction
and iterative Jordan-Wigner implementation.
"""
from collections import Counter
from functools import lru_cache
from itertools import product
import argparse
import json
from pathlib import Path
import sys

import sympy as s

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oph_openwave import charge, spin


def charge_polynomial():
    polynomial = Counter({(0, 0): 1})
    for electric_thirds, centre in [(-1, -2)]*3 + [(3, 3), (0, 3)]:
        updated = polynomial.copy()
        for (q, k), count in polynomial.items():
            updated[q+electric_thirds, (k+centre) % 6] += count
        polynomial = updated
    return polynomial


def charge_oracle(rows):
    actual = Counter((3*r['Q'], r['k']) for r in rows)
    expected = charge_polynomial()
    if actual != expected:
        raise ArithmeticError('charge/centre generating polynomial disagrees')
    return expected


def tensor_spin():
    # Cartesian vector representation, independent of fermion occupation lifts.
    vector = [s.Matrix(3, 3, lambda b, c: -s.I*s.LeviCivita(a, b, c)) for a in range(3)]
    pauli = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    total = [s.kronecker_product(j, s.eye(2)) + s.kronecker_product(s.eye(3), p/2)
             for j, p in zip(vector, pauli)]
    casimir = sum((j*j for j in total), s.zeros(6))
    projector = (s.Rational(15, 4)*s.eye(6)-casimir)/3
    return total, casimir, projector, pauli


def tensor_rotation(quaternion, pauli):
    w, x, y, z = quaternion
    cross = s.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    vector = s.eye(3)+2*w*cross+2*cross*cross
    spinor = w*s.eye(2)-s.I*sum((c*p for c, p in zip((x, y, z), pauli)), s.zeros(2))
    return s.kronecker_product(vector, spinor).applyfunc(s.expand)


def check_turn(matrix):
    if (matrix**4+s.eye(6)).applyfunc(s.simplify) != s.zeros(6):
        raise ArithmeticError('central rotation sign erased')
    if (matrix**8-s.eye(6)).applyfunc(s.simplify) != s.zeros(6):
        raise ArithmeticError('double turn fails')


def wedge_action(creator_modes, occupied, fermionic=True):
    if set(creator_modes) & set(occupied):
        return None, 0
    # Sort the concatenated wedge directly; no iterative CAR/bit-count routine.
    inversions = sum(a > b for a in creator_modes for b in occupied)
    sign = (-1)**inversions if fermionic else 1
    return tuple(sorted((*creator_modes, *occupied))), sign


def exchange_oracle(kernel, masks):
    polynomials = []
    for site in range(2):
        for col in range(6):
            terms = []
            for row, mask in enumerate(masks):
                value = kernel[row, col]
                if value:
                    if not value.is_Integer:
                        raise ArithmeticError('audit requires exact integer kernel coordinates')
                    modes = tuple(i+6*site for i in range(6) if mask & (1 << i))
                    terms.append((modes, int(value)))
            polynomials.append(terms)

    @lru_cache(maxsize=None)
    def act(index, occupied, fermionic=True):
        values = Counter()
        for modes, coefficient in polynomials[index]:
            target, sign = wedge_action(modes, occupied, fermionic)
            if sign:
                values[target] += coefficient*sign
        return tuple((target, coefficient) for target, coefficient in values.items() if coefficient)

    def compose(a, b, occupied, fermionic=True):
        values = Counter()
        for middle, first in act(b, occupied, fermionic):
            for target, second in act(a, middle, fermionic):
                values[target] += first*second
        return {target: coefficient for target, coefficient in values.items() if coefficient}

    nonzero, cases, witness = 0, 0, None
    for a, b, mask in product(range(6), range(6, 12), range(4096)):
        occupied = tuple(i for i in range(12) if mask & (1 << i))
        ab, ba = compose(a, b, occupied), compose(b, a, occupied)
        if ab != {target: -v for target, v in ba.items()}:
            raise ArithmeticError('independent wedge exchange failed')
        cases += 1
        nonzero += bool(ab)
        if ab and witness is None:
            witness = (a, b, occupied)
    if witness is None:
        raise ArithmeticError('independent exchange oracle was vacuous')
    a, b, occupied = witness
    mutated_ab, mutated_ba = compose(a, b, occupied, False), compose(b, a, occupied, False)
    if mutated_ab == {target: -v for target, v in mutated_ba.items()}:
        raise ArithmeticError('commuting mutation escaped independent wedge oracle')
    return cases, nonzero, {'left_column': a, 'right_column': b-6, 'occupied': list(occupied)}


def run():
    charge_result = charge.run()
    rows = charge.enumerate_code()
    census = charge_oracle(rows)
    try:
        charge_oracle(charge.enumerate_code(electric_scale=charge.F(1, 2)))
    except ArithmeticError:
        scale_mutation = True
    else:
        raise ArithmeticError('charge normalization mutation did not fire')
    halves = [j for j in range(6) if 2*j % 6 == 0]
    if halves != [0, 3] or charge_result['electron_charge_and_class']['k'] != 0:
        raise ArithmeticError('electron centre-class counterexample failed')
    if any(q == charge.F(3, 2) or q == charge.F(-3, 2) for q, _ in census):
        raise ArithmeticError('half electric charge entered declared code')

    masks, kernel, gram, left, producer_spin, producer_casimir, p_half, p_three_halves = spin.local_kernel()
    total, casimir, tensor_half, pauli = tensor_spin()
    if producer_casimir.charpoly().as_expr() != casimir.charpoly().as_expr():
        raise ArithmeticError('spin 1 tensor spin-half Casimir disagrees with actual kernel')
    if casimir.eigenvals() != {s.Rational(3, 4): 2, s.Rational(15, 4): 4}:
        raise ArithmeticError('tensor decomposition failed')
    max_character_difference = 0
    for q in spin.quaternion_elements():
        reference = tensor_rotation(q, pauli)
        actual = spin.exterior_rotation(spin.spinor(q))
        for tensor_p, actual_p in ((tensor_half, p_half), (s.eye(6)-tensor_half, p_three_halves)):
            difference = s.simplify(s.trace(tensor_p*reference)-s.trace(actual_p*actual))
            if difference != 0:
                raise ArithmeticError('independent tensor character disagrees')
    q = (s.sqrt(2)/2, 0, 0, s.sqrt(2)/2)
    quarter = tensor_rotation(q, pauli)
    check_turn(quarter)
    try:
        check_turn((1+s.I)/s.sqrt(2)*quarter)
    except ArithmeticError:
        turn_mutation = True
    else:
        raise ArithmeticError('phase-erasure mutation did not fire')
    cases, nonzero, exchange_witness = exchange_oracle(kernel, masks)
    spin_result = spin.run()
    if nonzero != spin_result['nonzero_exchange_comparisons']:
        raise ArithmeticError('producer and independent exchange support disagree')
    for result in (charge_result, spin_result):
        if result['full_criterion_met'] is not False:
            raise ArithmeticError('finite calculation promoted to full physical criterion')

    return {
        'schema': 'oph.openwave.charge_spin_adversarial_review.v1',
        'computational_pass': True,
        'reviewer_implementation': 'Separate agent; Laurent-polynomial charge census, vector tensor spinor representation and permutation-inversion exterior products.',
        'claims': [
            {'claim': 'All 32 charge/centre weights and absence of half-unit electric weights', 'verdict': 'CONFIRMED', 'evidence': 'Independent generating-polynomial census matches every joint charge/centre multiplicity.'},
            {'claim': 'Odd-centre modular no-halving proves electron indivisibility', 'verdict': 'REFUTED', 'origin': 'Initial WORK_ITEMS item 8 shortcut, not the corrected candidate claim', 'evidence': 'Electron-like singlet has k=0; 2j=0 admits j=0,3. Charge fractions and dynamical isolation require separate arguments.'},
            {'claim': 'Full OpenWave charge/quark criteria', 'verdict': 'QUALIFIED', 'evidence': 'Electric generator and normalization declared; no topological integral or confinement/energetic-fractionation experiment.'},
            {'claim': 'Actual six-state kernel decomposes as spin 1/2 plus spin 3/2', 'verdict': 'CONFIRMED', 'evidence': 'Independent spin 1 tensor spin 1/2 Casimir and all 48 projected characters match.'},
            {'claim': 'Full turn minus identity and double turn identity', 'verdict': 'CONFIRMED', 'evidence': 'Independent vector tensor spinor rotations give the central sign; phase-erasure mutation fails.'},
            {'claim': 'All 36 disjoint composite pairs anticommute on all 4096 states', 'verdict': 'CONFIRMED', 'evidence': f'{cases} independent wedge comparisons, {nonzero} nonzero products; commuting mutation fails.'},
            {'claim': 'Derived physical electron spin-statistics', 'verdict': 'QUALIFIED', 'evidence': 'CAR/spinor action are inputs. Both spin sectors occur; no electron assignment, Hamiltonian persistence or continuum spin-statistics theorem.'},
        ],
        'diagnostics': {'joint_charge_centre_census': [{'electric_charge_thirds': str(q), 'centre': k, 'multiplicity': n} for (q, k), n in sorted(census.items())],
                        'electron_equal_class_halves': halves, 'Casimir_multiplicities': {str(k): int(v) for k, v in casimir.eigenvals().items()},
                        'binary_cube_character_comparisons': 96, 'maximum_exact_character_difference': max_character_difference,
                        'exchange_basis_comparisons': cases, 'nonzero_exchange_comparisons': nonzero},
        'mutations': {'halved_electric_normalization_rejected': scale_mutation, 'central_phase_erasure_rejected': turn_mutation,
                      'commuting_creators_rejected': True},
        'mutation_witness': exchange_witness,
        'official_OpenWave_icons_assigned': False,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    rendered = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end='')
