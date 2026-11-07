"""A reduced native mixing channel cannot become a unitary by static coding.

This concerns unread carrier output only. Retained records, coherent
environments, active correction and noiseless memory are outside the claim.
"""
from collections import Counter
from fractions import Fraction as F

import numpy as np
import sympy as s

from ..interacting_extension._source import ROOT, PIN
from ..interacting_extension.source_selection import native_record_dynamics as native


SOURCE_FILES = (
    'model/particle_identification/passive_clock.py',
    'model/particle_identification/test_passive_clock.py',
    'model/particle_identification/PASSIVE_CLOCK.md',
    'model/interacting_extension/source_selection/native_record_dynamics.py',
    'model/interacting_extension/source_selection/PROOF.md',
    'model/interacting_extension/_source.py',
    'model/interacting_extension/source_manifest.json',
    'model/observer/native.py',
    'model/observer/quantum_mixing.py',
    'model/observer/source_discrimination/circuit.py',
    'model/common_source/slow_geometry.py',
)


def _integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f'{name} must be an integer at least {minimum}')
    return value


def unitary_error_lower(mixing_error, logical_dimension):
    """Half-diamond lower bound given a replacement-channel error upper bound."""
    d = _integer(logical_dimension, 'logical dimension', 2)
    if type(mixing_error) not in (int, F) or mixing_error < 0:
        raise ValueError('mixing error must be an exact nonnegative rational')
    return max(F(), 1-F(1, d*d)-F(mixing_error))


def native_mixing_upper(attempts, copies=1):
    """Outward rational upper bound for independent binary native carriers.

    sqrt(5)<9/4 gives r=(55+sqrt(5))/60 < 229/240. Tensor products
    of copies cost a factor copies by channel telescoping, not by a
    claim about correlated or interacting native carriers.
    """
    n = _integer(attempts, 'attempt count', 0)
    k = _integer(copies, 'independent copy count', 1)
    return min(F(1), 33*k*F(229, 240)**n)


def passive_error_lower(attempts, logical_dimension, copies=1):
    """Static CPTP encode/native evolution/decode versus any target unitary."""
    return unitary_error_lower(native_mixing_upper(attempts, copies), logical_dimension)


def source_contract_audit(repo=ROOT):
    """Recheck full quantum blocks, not just the classical population walk."""
    x = s.Symbol('x')
    exceptional = x*x-110*x+3020
    other_factors = (
        (x-54, 5), (x-52, 5), (x-50, 11), (x-48, 5),
        (x*x-100*x+2488, 1), (x*x-100*x+2490, 4),
        (x*x-100*x+2498, 4), (x**3-150*x*x+7480*x-124024, 5),
    )
    expected = s.Poly((x-60)*exceptional**3 *
                      s.prod(factor**multiplicity for factor, multiplicity in other_factors), x)
    # Exact real-root isolation certifies that no unexamined factor overtakes
    # the displayed exceptional root. Both exceptional roots are positive.
    for factor, _ in other_factors:
        polynomial = s.Poly(factor, x)
        intervals = polynomial.intervals(eps=s.Rational(1, 100))
        if (sum(multiplicity for _, multiplicity in intervals) != polynomial.degree()
                or any(not (0 < low <= high < 57) for (low, high), _ in intervals)):
            raise ArithmeticError('population factor is not confined below the extremal root')
    if not s.Rational(2)**2 < 5 < s.Rational(9, 4)**2:
        raise ArithmeticError('outward radical bound failed')

    rows = {}
    population = None
    for choice in ('C', 'E'):
        instrument = native.native_instrument(repo, choice)
        blocks = native.native_blocks(instrument)
        if instrument.dimension != 66 or sum(len(pairs) for _, pairs, _, _ in blocks) != 4356:
            raise ArithmeticError('native carrier or complete quantum operator domain changed')
        minimum, off_diagonal_maximum = 60, 0
        population_here = None
        census = Counter()
        for difference, pairs, waiting, moving in blocks:
            matrix = waiting+moving
            if matrix.dtype.kind not in 'iu' or (matrix < 0).any() or not np.array_equal(matrix, matrix.T):
                raise ArithmeticError('exact nonnegative self-adjoint native block required')
            census[len(pairs)] += 1
            minimum = min(minimum, int((2*np.diag(matrix)-matrix.sum(axis=1)).min()))
            if any(difference):
                if any(a == b for a, b in pairs):
                    raise ArithmeticError('population direction entered a coherence block')
                off_diagonal_maximum = max(off_diagonal_maximum, int(matrix.sum(axis=1).max()))
            else:
                if any(a != b for a, b in pairs) or population_here is not None:
                    raise ArithmeticError('incorrect population block')
                population_here = matrix
        if minimum <= 0 or off_diagonal_maximum >= 57 or dict(census) != {66: 1, 10: 132, 1: 2970}:
            raise ArithmeticError('full native spectral bounds changed')
        if population_here is None or not np.all(population_here.sum(axis=1) == 60):
            raise ArithmeticError('native population transfer is not stochastic')
        if population is None:
            population = population_here
            if s.Matrix(population.tolist()).charpoly(x).as_poly() != expected:
                raise ArithmeticError('exact native population characteristic polynomial changed')
        elif not np.array_equal(population, population_here):
            raise ArithmeticError('C/E population laws differ')
        rows[choice] = {'quantum_block_size_census': {str(size): count for size, count in sorted(census.items())},
                        'off_diagonal_spectral_radius_numerator_upper_over_60': off_diagonal_maximum,
                        'full_quantum_gershgorin_lower_numerator_over_60': minimum,
                        'complete_operator_directions': 4356}
    return {'source_commit': PIN, 'carrier_dimension': 66, 'candidates': rows,
            'nonconstant_quantum_spectral_radius': '(55+sqrt(5))/60',
            'rational_spectral_radius_upper': '229/240',
            'native_half_diamond_mixing_upper': 'min(1,33*((55+sqrt(5))/60)**n)',
            'population_characteristic_polynomial': str(expected.as_expr()),
            'source_choices': 'binary invariant sector; uniform independent seams and fair coins; un-tilted C or E',
            'retained_native_records_available_to_decoder': False}


def certificate(repo=ROOT):
    lower = passive_error_lower(200, 2)
    if lower <= F(74, 100):
        raise ArithmeticError('the passive-clock example lost its strict error bound')
    return {'source_audit': source_contract_audit(repo),
            'logical_dimension': 2, 'attempts': 200, 'independent_carriers': 1,
            'outward_rational_mixing_upper': str(native_mixing_upper(200)),
            'passive_unitary_error_lower': str(lower),
            'strict_rational_error_lower': '74/100',
            'comparison': 'arbitrary static CPTP encoding and decoding around unread native evolution',
            'record_assisted_active_or_noiseless_memory_routes_excluded': False,
            'source_to_interacting_QFT_bridge_derived': False}
