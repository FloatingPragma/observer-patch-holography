"""Exact positive spectral filters and their native-record estimators.

The effects act on the reconstructed transfer Hilbert space. Evaluating
their probabilities from several record experiments does not execute a
coherent energy measurement on the native working carrier.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import comb

from ..native_exchange_clock import exchange
from ..native_vacuum_response import response


def integer(value, lower=0):
    return exchange.integer(value, lower)


def rational(value):
    if type(value) not in (int, F):
        raise ValueError('an exact rational is required')
    return F(value)


def plan(seams, maximum_degree, correlation_degree, order):
    """Bernstein effects for Q=P**b on the supplied reducing degree block."""
    integer(seams, 1); integer(maximum_degree, 1)
    integer(correlation_degree, 1); integer(order, 1)
    product = maximum_degree * correlation_degree
    if seams < 2 * product or order > 64:
        raise ValueError('M>=2*r*Delta and order from one through64 required')
    block = seams // (2 * product)
    coefficients = tuple(tuple(
        ((-1)**(j-k) * comb(order, k) * comb(order-k, j-k) if j >= k else 0)
        for j in range(order+1)) for k in range(order+1))
    amplification = tuple(sum(abs(x) for x in row[1:]) for row in coefficients)
    if tuple(sum(row[j] for row in coefficients) for j in range(order+1)) != (1,) + (0,)*order:
        raise ArithmeticError('the complete polynomial family lost its identity sum')
    return {'seams': seams, 'maximum_degree': maximum_degree,
            'correlation_degree': correlation_degree, 'order': order,
            'block_attempts': block, 'transfer_floor': F(1, 2),
            'energy_scale': F(2 * seams, block),
            'lag_attempts': tuple(block*j for j in range(order+1)),
            'power_coefficients': coefficients,
            'bin_amplifications': amplification,
            'total_l1_amplification': sum(amplification),
            'local_exchange_energy_error_upper': exchange.degree_budget(
                seams, maximum_degree, correlation_degree)['generator_error_upper']}


def histogram(coefficients, moments):
    """Unclipped estimator, with the algebraic zeroth moment fixed to one.

    Approximate nonzero-lag moments can yield negative entries. They are
    retained, together with their certified estimation error.
    """
    if (type(coefficients) is not tuple or type(moments) is not tuple
            or not coefficients or len(moments) != len(coefficients)):
        raise ValueError('one literal moment tuple and square coefficient table required')
    width = len(moments)
    values = tuple(rational(x) for x in moments)
    if values[0] != 1:
        raise ValueError('the exact normalized zeroth moment must equal one')
    if any(type(row) is not tuple or len(row) != width
           or any(type(x) is not int for x in row) for row in coefficients):
        raise ValueError('a square exact integer coefficient table required')
    if tuple(sum(row[j] for row in coefficients) for j in range(width)) != (1,)+(0,)*(width-1):
        raise ValueError('the polynomial effects must sum exactly to the identity')
    return tuple(sum((F(c)*v for c,v in zip(row, values, strict=True)), F())
                 for row in coefficients)


def negative_log_interval(value, terms=32):
    """Exact rational enclosure of -log(value), with no floating log."""
    value = rational(value); integer(terms, 1)
    if not 0 < value <= 1 or terms > 4096:
        raise ValueError('a probability in(0,1] and at most4096 terms required')
    deficit = 1-value
    lower = sum((deficit**j/j for j in range(1, terms+1)), F())
    return lower, lower + deficit**(terms+1)/((terms+1)*value)


def energy_labels(seams, block, order, terms=32):
    """All n+1 labels, including the explicit zero-count overflow bin."""
    integer(seams, 1); integer(block, 1); integer(order, 1)
    if order > 64:
        raise ValueError('the bounded detector order is at most64')
    scale = F(2*seams, block)
    return (None,) + tuple(tuple(scale*x for x in negative_log_interval(F(k,order), terms))
                            for k in range(1,order+1))


def resolution_budget(seams, maximum_degree, correlation_degree, order, margin):
    """Uniform probability of a wrong finite transfer-energy label.

    This is the exact Bernstein POVM's spectral resolution, separate from
    estimating its probabilities with a signed record statistic.
    """
    detector = plan(seams, maximum_degree, correlation_degree, order)
    margin = rational(margin)
    if not 0 < margin < F(1,2):
        raise ValueError('a margin strictly between zero and one half required')
    return {'transfer_probability_margin': margin,
            'bad_label_probability_upper': min(F(1), 1/(4*order*margin*margin)),
            'counted_energy_error_on_good_labels': detector['energy_scale']*margin/(F(1,2)-margin),
            'capacity_independent_energy_error_upper':
                8*maximum_degree*correlation_degree*margin/(F(1,2)-margin),
            'zero_count_overflow_probability_upper': F(1,2**order),
            'exchange_energy_error_on_good_labels':
                detector['energy_scale']*margin/(F(1,2)-margin)
                + detector['local_exchange_energy_error_upper']}


def estimation_budget(detector, variance, normalized_systematic_error,
                      total_variation_sampling_error, failure_probability):
    """Complete-vector error for independent native lag repetitions.

    The bound is half the l1 distance to the ideal probability vector;
    the finite-sample vector itself is permitted to be signed.
    """
    variance, eta, error, probability = map(rational, (
        variance, normalized_systematic_error, total_variation_sampling_error,
        failure_probability))
    if not 0 < variance <= 1 or eta < 0 or error <= 0 or not 0 < probability < 1:
        raise ValueError('positive variance/error, nonnegative bias and probability in(0,1) required')
    order = detector['order']; amplification = detector['total_l1_amplification']
    integer(order,1); integer(amplification,1)
    required = order * amplification**2 / (4*probability*variance**2*error**2)
    repeats = (required.numerator+required.denominator-1)//required.denominator
    return {'estimated_lag_settings': order, 'zeroth_moment_is_exact': True,
            'independent_repetitions_per_setting': repeats,
            'total_independent_repetitions': order*repeats,
            'half_l1_systematic_error_upper': F(amplification,2)*eta,
            'half_l1_sampling_error_upper': error,
            'failure_probability_upper': probability,
            'stop_allowance_for_uniform_stop_probability': f'{F(amplification,2)/variance} * p_stop',
            'signed_estimates_are_not_clipped_or_renormalized': True}


def exact_case(family, probe, degree, order):
    ports, seams, maximum_degree = exchange.graph_data(family)
    detector = plan(seams, maximum_degree, degree, order)
    moments = response.exact_native_lags(family, probe, degree, detector['lag_attempts'])
    weights = histogram(detector['power_coefficients'], moments)
    if min(weights) < 0 or sum(weights) != 1:
        raise ArithmeticError('native exact lag moments violate positive spectral completeness')
    return {'ports': ports, 'detector': detector, 'normalized_native_lag_scores': moments,
            'spectral_bin_weights': weights,
            'counted_energy_labels': energy_labels(seams, detector['block_attempts'], order)}
