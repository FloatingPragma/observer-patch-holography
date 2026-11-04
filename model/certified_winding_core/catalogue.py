"""Prove retained interval targets with a fresh integer-certified full basis."""
from fractions import Fraction as F
import hashlib
import math

import numpy as np

from . import core, spectrum

BITS = 26
COUPLINGS = (F(1, 65536), F(1, 2**24), F(1, 2**40))


def certified_case(target):
    radius, degree, mass = target['radius'], target['degree'], int(target['mass'])
    n = 2*(2*radius+1)**3
    values = target['all_eigenvalues']
    if (target['mass'] != str(mass) or target['spacing'] != '1' or target['velocity'] != '1'
            or target['dimension'] != n or len(values) != n
            or any(type(x) is not float or not math.isfinite(x) for x in values)
            or values != sorted(values)):
        raise ValueError('complete ordered finite target spectrum required')
    result = spectrum.calculate(radius, degree, mass, BITS)
    scale = 1 << BITS
    g, r, s0, largest = F(2*n, scale), F(2*(mass+2)*n, scale), F(n, scale), F(mass+2)
    if (F(result['integer_Gram_defect_squared']) > g*g
            or F(result['integer_intertwining_residual_squared']) > r*r
            or F(result['source_Frobenius_error_squared_upper']) > s0*s0):
        raise ArithmeticError('the exact residual exceeds a declared portable bound')
    error = 2*(s0+(r+2*largest*(1+g)*g)/(1-g))
    centers = [F(round(F(x)*scale), scale) for x in values]
    proposed = [[x-error, x+error] for x in centers]
    fresh = [[F(x) for x in pair] for pair in result['ordered_energy_intervals']]
    if any(abs((a+b)/2) > largest for a, b in fresh):
        raise ArithmeticError('the proposed diagonal exceeds the declared magnitude bound')
    for fixed, proved in zip(proposed, fresh, strict=True):
        if not fixed[0] <= proved[0] <= proved[1] <= fixed[1]:
            raise ArithmeticError('the independently certified spectrum escapes a retained target')
    # The exact source construction is independent of the eigensolver.
    hr, hi, _ = spectrum.source(radius, degree, mass, BITS)
    digest = hashlib.sha256()
    for a in (hr, hi):
        for array in (a.indptr, a.indices, a.data):
            digest.update(np.asarray(array, dtype='<i8').tobytes())
    threshold = F(mass, 4)
    inside = lambda a, b: -threshold < a <= b < threshold
    if any(not inside(a, b) and not (b < -threshold or a > threshold) for a, b in proposed):
        raise ArithmeticError('the retained window is not resolved by the declared intervals')
    zero_count = sum(a <= 0 <= b for a, b in proposed)
    exact_two = degree % 2 != 0 and zero_count == 2
    gap = min(min(abs(a), abs(b)) for a, b in proposed if not a <= 0 <= b)
    output = {'radius': radius, 'degree': degree, 'mass': mass, 'dimension': n,
              'spacing': '1', 'velocity': '1', 'dyadic_bits': BITS,
              'source_approximation_sha256': digest.hexdigest(),
              'source_error_cap': str(s0), 'Gram_defect_cap': str(g),
              'integer_residual_cap': str(r), 'retained_target_error': str(error),
              'all_integer_intermediates_proved_int64_safe': True,
              'all_source_eigenvalues_enclosed_with_multiplicity': True,
              'ordered_energy_intervals': [[str(x) for x in pair] for pair in proposed],
              'window_absolute_energy': str(threshold),
              'window_rank_exact': sum(inside(a, b) for a, b in proposed),
              'intervals_containing_zero': zero_count,
              'exact_kernel_dimension': 2 if exact_two else (0 if zero_count == 0 else None),
              'inversion_index': 2 if degree % 2 else None,
              'nonzero_gap_lower': str(gap) if exact_two or zero_count == 0 else None,
              'infinite_volume_kernel_or_gap_certified': False}
    if output['window_rank_exact'] != target['window_rank']:
        raise ArithmeticError('certification disagrees with the retained window count')
    output['full_Fock_contact_cases'] = [core.bounds(output, u) for u in COUPLINGS]
    return output
