"""Exact all-sector entropy data for the native canonical reference.

Logarithms are retained as rational linear combinations of log(rational).
The source supplies population and occupation; neither is a proper volume.
"""
from fractions import Fraction as F
from math import comb

from model.native_capacity_restriction import restriction


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('an exact integer in the stated domain is required')
    return value


def blocks(population, occupied, retained):
    """Every local number block, including impossible zero-weight outcomes."""
    law = restriction.number_law(population, occupied, retained)
    return tuple({'number': k, 'dimension': comb(retained, k),
                  'probability': probability,
                  'eigenvalue': probability / comb(retained, k)}
                 for k, probability in enumerate(law))


def supported(rows):
    return tuple(row for row in rows if row['probability'])


def log_expression(terms):
    """Canonical rational coefficients; no rounded logarithms or prime factoring."""
    value = {}
    for argument, coefficient in terms:
        if type(argument) is not F or argument <= 0 or type(coefficient) is not F:
            raise ValueError('positive rational log arguments and rational coefficients required')
        if argument != 1:
            value[argument] = value.get(argument, F()) + coefficient
    return tuple((argument, coefficient) for argument, coefficient in sorted(value.items()) if coefficient)


def combine(*weighted):
    return log_expression((argument, F(scale) * coefficient)
                          for scale, terms in weighted for argument, coefficient in terms)


def eigenvalue_log_ratio(population, occupied, retained, lower_number):
    """w_(k+1)/w_k without computing a macroscopic binomial coefficient."""
    rows = blocks(population, occupied, retained)
    integer(lower_number)
    if lower_number >= retained or not all(rows[k]['probability'] for k in (lower_number, lower_number + 1)):
        raise ValueError('two adjacent supported occupation sectors required')
    k = lower_number
    return F(occupied - k, population - retained - occupied + k + 1)


def entropy_expression(rows):
    return log_expression((row['eigenvalue'], -row['probability']) for row in supported(rows))


def first_law_expression(rows, sector_tangent):
    """Pair -log(reference) with a trace-zero Hermitian tangent by block traces."""
    if (type(sector_tangent) is not tuple or len(sector_tangent) != len(rows)
            or any(type(x) is not F for x in sector_tangent) or sum(sector_tangent)):
        raise ValueError('a complete exact trace-zero sector tangent is required')
    if any(x and not row['probability'] for row, x in zip(rows, sector_tangent, strict=True)):
        raise ValueError('the tangent must remain on the faithful reference support')
    return log_expression((row['eigenvalue'], -x)
                          for row, x in zip(rows, sector_tangent, strict=True) if x)


def modular_flow_allowance(minimum_eigenvalue, half_trace_distance, modular_time):
    """Operator-norm flow error per unit norm observable on a common support.

    Inputs must certify the trace distance of two normalized densities and
    the minimum eigenvalue of the sector-scalar reference. This function does
    not certify that premise from source marginals or factorize references.
    """
    m, error, time = minimum_eigenvalue, half_trace_distance, modular_time
    if any(type(x) is not F for x in (m, error, time)) or not 0 < m <= 1 or not 0 <= error < m:
        raise ValueError('exact rational 0 <= error < minimum eigenvalue <= 1 required')
    return min(F(2), 2 * abs(time) * error / (m - error))


def first_law_remainder(minimum_eigenvalue, hs_squared, parameter, tangent_norm_upper):
    """D(sigma+tX||sigma) <= t^2 ||X||_HS^2/(2(m-|t| ||X||))."""
    values = minimum_eigenvalue, hs_squared, parameter, tangent_norm_upper
    if any(type(x) is not F for x in values):
        raise ValueError('exact rational response bounds required')
    m, hs, time, norm = values
    floor = m - abs(time) * norm
    if not 0 < m <= 1 or hs < 0 or norm < 0 or floor <= 0:
        raise ValueError('a positive common spectral floor is required')
    return time * time * hs / (2 * floor)


def finite_entropy_data(population, occupied, retained):
    """Explicit small-state factorization and symbolic entropy identity.

    The cap protects the output from materializing cosmological binomials.
    All-capacity modular data use only ``blocks`` and falling-factorial laws.
    """
    rows = blocks(population, occupied, retained)
    if population > 4096:
        raise ValueError('explicit exterior dimensions are limited to 4096 ports')
    total = comb(population, occupied)
    result = []
    for row in rows:
        k = row['number']
        exterior = comb(population - retained, occupied - k) if 0 <= occupied - k <= population - retained else 0
        if row['probability'] != F(row['dimension'] * exterior, total):
            raise ArithmeticError('native restriction disagrees with tensor multiplicity')
        result.append({**row, 'exterior_dimension': exterior})
    # w_k = e_k / D: use this identity to combine the logarithms exactly.
    bulk = log_expression([(F(total), F(1))] +
                          [(F(row['exterior_dimension']), -row['probability']) for row in supported(result)])
    edge = log_expression((F(row['exterior_dimension']), row['probability']) for row in supported(result))
    full = log_expression([(F(total), F(1))])
    if combine((1, bulk), (1, edge)) != full:
        raise ArithmeticError('same-state bulk/edge entropy split failed')
    return tuple(result), bulk, edge, full


def sector_redistribution(population, occupied, retained, number=0):
    """A positive state tangent changing multiplicity entropy on a fixed cut."""
    rows = blocks(population, occupied, retained)
    ratio = eigenvalue_log_ratio(population, occupied, retained, number)
    tangent = tuple(F(-1 if k == number else 1 if k == number + 1 else 0) for k in range(retained + 1))
    radius = min(rows[number]['probability'], rows[number + 1]['probability']) / 2
    if radius <= 0:
        raise ArithmeticError('redistribution has no two-sided positive interval')
    return {'tangent': tangent, 'positive_parameter_radius': radius,
            'edge_derivative': log_expression([(ratio, F(1))]),
            'bulk_derivative': log_expression([(ratio, F(-1))]),
            'edge_ratio': ratio, 'static_cut_area_derivative': F(0)}
