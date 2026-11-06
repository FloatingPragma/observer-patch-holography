"""All-quadratic exclusion in the native two-particle occupation dictionary."""
from __future__ import annotations

from fractions import Fraction as Q

from ..native_effective_fields import construction as field
from . import certificate as c

PAIRS = ((0, 1), (2, 3), (0, 2), (1, 3))
SIGNS = (1, 1, -1, -1)


def pair_indices(pairs=PAIRS):
    c.require(type(pairs) is tuple and len(pairs) == 4, 'four occupation pairs required')
    result = []
    for pair in pairs:
        c.require(type(pair) is tuple and len(pair) == 2 and all(type(v) is int for v in pair)
                  and 0 <= pair[0] < pair[1] < 12, 'ordered distinct native ports required')
        state = tuple(int(v in pair) for v in range(12))
        result.append(c.source.INDEX[state])
    return tuple(result)


def contrast(diagonal, pairs=PAIRS):
    c.require(type(diagonal) is tuple and len(diagonal) == c.D, 'complete 66-state diagonal required')
    return sum(sign*c.rational(diagonal[i]) for sign, i in zip(SIGNS, pair_indices(pairs)))


def certificate(order=16):
    c.integer(order, 1, 24, 'logarithm order')
    counts, denominator = c.source.source_counts()
    g = field.generator()
    c.require(all(g[i][j] == denominator*int(i == j)-counts[i][j]
                  for i in range(c.D) for j in range(c.D)),
              'reused local generator differs from native transfer counts')
    local_contrast = contrast(tuple(Q(g[i][i], denominator) for i in range(c.D)))
    c.require(local_contrast == Q(-1, 15), 'native local interaction witness changed')
    polynomial = field.log_polynomial(order)
    r = field.log_remainder(order)
    logarithm_contrast = contrast(tuple(polynomial[i][i] for i in range(c.D)))
    # Each diagonal remainder lies in [0,r]. There are exactly two plus and
    # two minus entries, hence a 2r (not 4r) enclosure for their contrast.
    lower, upper = logarithm_contrast-2*r, logarithm_contrast+2*r
    c.require(upper < 0 or lower > 0, 'logarithm interaction contrast is not separated from zero')
    distance = min(abs(lower), abs(upper))/4
    return {'occupation_pairs': PAIRS, 'basis_indices': pair_indices(),
            'all_native_sector_states_retained': c.D,
            'local_generator_I_minus_P_contrast': local_contrast,
            'local_generator_distance_from_all_quadratics_lower': abs(local_contrast)/4,
            'logarithm_order': order, 'logarithm_contrast_interval': (lower, upper),
            'logarithm_operator_distance_from_all_quadratics_lower': distance,
            'strict_logarithm_distance_simple_lower': Q((distance*10**8).numerator//(distance*10**8).denominator,10**8),
            'arbitrary_one_body_hoppings_and_potentials_allowed_in_excluded_class': True,
            'native_occupation_dictionary_fixed': True,
            'arbitrary_reencoding_excluded': False,
            'relativistic_or_chiral_continuum_claim': False}
