"""Exact bounds and sparse response moments for the literal fair source.

All generators act on the configuration transfer Hilbert space. An
exponential endpoint comparison does not implement a coherent source clock.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
from math import comb, factorial

from ..native_correlation_sectors import sectors as sector


def integer(value, lower=0):
    if type(value) is not int or value < lower:
        raise ValueError('an exact integer in the stated domain is required')
    return value


def rational(value, lower=F(0)):
    if type(value) not in (int, F) or value < lower:
        raise ValueError('an exact nonnegative rational is required')
    return F(value)


def graph_data(family):
    family = sector.fair_graph(family)
    degrees = [0] * len(family.weights)
    for a, b in family.seams:
        degrees[a] += 1
        degrees[b] += 1
    return len(degrees), len(family.seams), max(degrees)


def generator(family, values, degree):
    """Apply sum_e(I-S_e) exactly, independently of transfer counts_step."""
    ports, _, _ = graph_data(family)
    values, _ = sector.coefficients(values, ports, degree)
    neighbours = [[] for _ in family.weights]
    for a, b in family.seams:
        neighbours[a].append(b)
        neighbours[b].append(a)
    output = defaultdict(int)
    for subset, value in values.items():
        occupied = set(subset)
        for a in subset:
            for b in neighbours[a]:
                if b not in occupied:
                    target = tuple(sorted((occupied - {a}) | {b}))
                    output[subset] += value
                    output[target] -= value
    return {subset: value for subset, value in output.items() if value}


def degree_budget(seams, maximum_degree, correlation_degree):
    """Bound the rescaled logarithm on degrees up to correlation_degree."""
    integer(seams, 1); integer(maximum_degree, 1); integer(correlation_degree)
    product = correlation_degree * maximum_degree
    if product >= seams:
        raise ValueError('the bounded-degree logarithm requires r*Delta<M')
    return {'spectral_ceiling': 2 * product,
            'generator_error_upper': F(product * product, seams - product)}


def native_time(seams, attempts):
    integer(seams, 1); integer(attempts)
    return F(attempts, 2 * seams)


def response_budget(seams, attempts, maximum_degree=None, correlation_degree=None):
    """Endpoint contraction error at exact native time, never a history norm.

    The full-sector term uses the owning federation floor P>=10I/21.
    Optional degree data supplies a second, independent bound.
    """
    time = native_time(seams, attempts)
    full = F(0) if attempts == 0 else min(F(1), F(3, 5 * attempts))
    value = {'native_time': time, 'full_sector_response_error_upper': full}
    if maximum_degree is not None or correlation_degree is not None:
        budget = degree_budget(seams, maximum_degree, correlation_degree)
        bounded = min(F(1), time * budget['generator_error_upper'])
        value.update(budget)
        value['bounded_degree_response_error_upper'] = bounded
        value['best_response_error_upper'] = min(full, bounded)
    else:
        value['best_response_error_upper'] = full
    return value


def rounded_time_budget(seams, time, maximum_degree=None, correlation_degree=None):
    integer(seams, 1); time = rational(time)
    scaled = 2 * seams * time
    attempts = scaled.numerator // scaled.denominator
    remainder = scaled - attempts
    full = (F(0) if time == 0 else F(1) if attempts == 0
            else min(F(1), (F(3, 5) + remainder / 2) / attempts))
    value = {'time': time, 'attempts': attempts, 'fractional_attempt': remainder,
             'full_sector_response_error_upper': full}
    if maximum_degree is not None or correlation_degree is not None:
        budget = degree_budget(seams, maximum_degree, correlation_degree)
        bounded = min(F(1), native_time(seams, attempts) * budget['generator_error_upper']
                      + remainder * F(budget['spectral_ceiling'], 2 * seams))
        value.update(budget)
        value['bounded_degree_response_error_upper'] = bounded
        value['best_response_error_upper'] = min(full, bounded)
    else:
        value['best_response_error_upper'] = full
    return value


def full_generator_obstruction(ports, seams, occupied):
    """A trace-derived lower bound, conditional on a positive source transfer."""
    integer(ports, 2); integer(seams, 1); integer(occupied)
    if occupied > ports:
        raise ValueError('occupation exceeds the port count')
    delta = F(occupied * (ports - occupied), ports * (ports - 1))
    return {'mean_generator_eigenvalue': 2 * seams * delta,
            'rescaled_log_minus_generator_norm_lower': seams * delta * delta,
            'mean_transfer_deficit': delta}


def moments(family, values, degree, order):
    """Normalized exact spectral moments of the local exchange generator."""
    ports, _, _ = graph_data(family)
    values, _ = sector.coefficients(values, ports, degree)
    integer(order)
    if order > 12 or not values:
        raise ValueError('a nonzero probe and at most twelve sparse powers are required')
    norm = sector.dot(values, values)
    current = values
    answer = []
    for k in range(order + 1):
        answer.append(F(sector.dot(values, current), norm))
        if k < order:
            current = generator(family, current, degree)
    return tuple(answer)


def response_intervals(family, values, degree, attempts, order=8):
    """Exact enclosures for both the native lag and local exchange response.

    This evaluates sparse spectral moments. It is not a completed expanded
    interpreter run. The finite-program budget belongs to the source owner.
    """
    ports, seams, maximum_degree = graph_data(family)
    integer(attempts); integer(order, 2)
    if order > 12 or order % 2:
        raise ValueError('an even sparse polynomial order from two to twelve is required')
    values, actual_degree = sector.coefficients(values, ports, degree)
    if not values or actual_degree != degree or sector.lower(values, ports, degree):
        raise ValueError('a nonzero primitive-degree probe is required')
    budget = degree_budget(seams, maximum_degree, degree)
    mu = moments(family, values, degree, order)
    time = native_time(seams, attempts)
    exponential = tuple(sum((-time)**j * mu[j] / factorial(j) for j in range(m + 1))
                        for m in (order - 1, order))
    native = tuple(sum((-1)**j * F(comb(attempts, j), (2 * seams)**j) * mu[j]
                       for j in range(min(m, attempts) + 1)) for m in (order - 1, order))
    # Both exact responses have positive spectral measures on [0,1].
    exponential = (max(F(0), exponential[0]), min(F(1), exponential[1]))
    native = (max(F(0), native[0]), min(F(1), native[1]))
    if exponential[0] > exponential[1] or native[0] > native[1]:
        raise ArithmeticError('positive response enclosure is inconsistent')
    return {'attempts': attempts, 'native_time': time, 'polynomial_order': order,
            'moments': mu, 'local_exchange_response_interval': exponential,
            'native_response_interval': native,
            'response_budget': response_budget(seams, attempts, maximum_degree, degree),
            'exponential_interval_width_upper':
                (time * budget['spectral_ceiling'])**order / factorial(order)}


def preparation_composition(variance, preparation_score_error, response_error):
    """Compose an unnormalized score error with normalized endpoint error."""
    variance = rational(variance); preparation_score_error = rational(preparation_score_error)
    response_error = rational(response_error)
    if not 0 < variance <= 1:
        raise ValueError('a positive source probe variance at most one is required')
    return preparation_score_error / variance + response_error
