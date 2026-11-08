"""Exact charged bounded functional calculus on the actual selected cube."""
from fractions import Fraction as F
from math import lcm
import numpy as np

from model.protected_defect_catalog import charge
from model.centered_particle_spectrum.certificate import decimal_bounds
from model.spin_operational_position.bounds import generic_instrument_plan


def operators():
    op = charge.operators()
    indices = np.arange(256)
    x = [F(sum(op['signs'][j]*(1-2*((word>>j)&1)) for j in range(8)), 4)
         for word in indices]
    diagonal = [value/(1+value*value) for value in x]
    if any(abs(value) > F(1, 2) for value in diagonal):
        raise ArithmeticError('bounded rational field function exceeds one half')
    denominator = lcm(*(value.denominator for value in diagonal))
    numerator = np.array([int(value*denominator) for value in diagonal], dtype=np.int64)
    w = op['walsh']
    transformed = w.T@(numerator[:, None]*w)
    charges = np.array(op['charges'])
    plus = w@(transformed*(charges[:, None]-charges[None, :] == 1))@w.T
    divisor = denominator*256**2
    minus = plus.T.copy()
    for sign, matrix in ((1, plus), (-1, minus)):
        if not np.array_equal(op['twice_charge']@matrix-matrix@op['twice_charge'], 2*sign*matrix):
            raise ArithmeticError('bounded creator does not carry the actual source charge')
    permutation = op['charge_reversing_permutation']
    if not np.array_equal(plus[np.ix_(permutation, permutation)], -minus):
        raise ArithmeticError('actual cube symmetry does not pair the bounded creators')
    return op, tuple(diagonal), plus, minus, divisor


def build():
    op, diagonal, plus, minus, divisor = operators()
    g, error, proof = charge.ground_input()
    norm = int(g@g)
    rows = []
    probability_sum = F()
    for sign, matrix in ((1, plus), (-1, minus)):
        vector = matrix.astype(object)@g
        probability = F(int(vector@vector), divisor*divisor*norm)
        probability_sum += probability
        # ||B||<=1/2, so the effect has norm at most1/4.
        lower, upper = max(F(), probability-error/2), min(F(1, 4), probability+error/2)
        if lower <= 0:
            raise ArithmeticError('same finite source charged preparation has zero certified weight')
        rows.append({'source_charge': sign, 'original_probability_interval': decimal_bounds(lower, upper, 12),
            'effect_operator_norm_upper': '1/4',
            'all_256_spin_basis_rows_retained': True})
    return {'same_selected_finite_source_ground_certificate': proof,
        'local_test_field': 'x=sum_x epsilon_x Z_x/4; f(x)=x/(1+x^2)',
        'exact_charge_projection': 'sum_q P_(q+1) f(x) P_q, with its adjoint',
        'all_256_diagonal_function_values': [str(x) for x in diagonal],
        'creator_matrix_integer_divisor': divisor,
        'charge_operator_commutator_exact': True,
        'cube_symmetry_pairs_creator_spectra': True,
        'creation_records': rows,
        'complement_original_probability_interval': decimal_bounds(max(F(), 1-probability_sum-error), min(F(1), 1-probability_sum+error), 12),
        'complement_effect_operator_lower': '1/2',
        'finite_native_instrument_plan': generic_instrument_plan(8, 3, F(1, 10**10)),
        'literal_full_8_cube_creator_tape_executed_here': False,
        'finite_cube_substituted_for_ordered_continuum_state': False}
