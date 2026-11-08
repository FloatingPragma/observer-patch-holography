"""Exact continuation coefficients, record budgets and response witnesses.

The unitary group belongs to the reconstructed transfer Hilbert space. Its
matrix elements are statistical reconstructions, not native unitary commands.
"""
from __future__ import annotations

from fractions import Fraction as F
from math import comb, factorial

from ..native_exchange_clock import exchange as exchange
from ..native_correlation_preparation import preparation as preparation


def integer(value, minimum=0):
    return exchange.integer(value, minimum)


def rational(value):
    if type(value) not in (int, F):
        raise ValueError('an exact rational is required')
    return F(value)


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def multiply(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def scale(a, b):
    return a[0]*b, a[1]*b


def absolute_upper(a):
    return abs(a[0])+abs(a[1])


def continuation(seams, maximum_degree, correlation_degree, time, order):
    """A polynomial in P**block for exp(-it[-2M log P]).

    The spectral assertion is conditional on the owning primitive-degree
    bound. The source packet supplies that bound for its actual fair graphs.
    """
    integer(seams, 1); integer(maximum_degree, 1)
    integer(correlation_degree, 1); integer(order)
    time = rational(time)
    product = maximum_degree*correlation_degree
    if seams < 2*product or order > 96:
        raise ValueError('M>=2*r*Delta and polynomial order at most96 required')
    block = seams//(2*product)
    phase = 2*seams*time/block
    exponent = max(1, (abs(phase).numerator+abs(phase).denominator-1)//abs(phase).denominator)
    if exponent > 128:
        raise ValueError('the exact coefficient evaluator bounds ceil(abs(phase)) by128')
    coefficients = [(F(1), F(0))]
    for k in range(order):
        coefficients.append(scale(multiply(coefficients[-1], (F(k), -phase)), F(1, k+1)))
    powers = []
    for j in range(order+1):
        entry = (F(0), F(0))
        for k in range(j, order+1):
            entry = add(entry, scale(coefficients[k], (-1)**j*comb(k, j)))
        powers.append(entry)
    tail = (F(0) if not time else F(2**exponent)-sum(
        (F(comb(exponent+k-1, k), 2**k) for k in range(order+1)), F()))
    amplification = sum((absolute_upper(z) for z in powers), F())
    if tail < 0 or sum((z[0] for z in powers), F()) != 1 or sum((z[1] for z in powers), F()):
        raise ArithmeticError('the exact continuation polynomial lost its constant mode')
    return {'seams': seams, 'maximum_degree': maximum_degree,
            'correlation_degree': correlation_degree, 'time': time,
            'block_attempts': block, 'block_spectral_floor': F(1, 2),
            'phase_parameter': phase, 'coefficient_majorant_integer': exponent,
            'order': order, 'lag_attempts': tuple(block*j for j in range(order+1)),
            'power_coefficients': tuple(powers), 'operator_error_upper': tail,
            'coefficient_l1_upper': amplification,
            'local_exchange_unitary_error_upper': abs(time)*exchange.degree_budget(
                seams, maximum_degree, correlation_degree)['generator_error_upper']}


def reconstruct(coefficients, normalized_lag_responses):
    """Exact complex linear statistic, with no positivity or clipping step."""
    if type(coefficients) is not tuple or type(normalized_lag_responses) is not tuple:
        raise ValueError('literal coefficient and response tuples required')
    if len(coefficients) != len(normalized_lag_responses) or not coefficients:
        raise ValueError('one response for every continuation setting required')
    total = (F(), F())
    for z, value in zip(coefficients, normalized_lag_responses, strict=True):
        if type(z) is not tuple or len(z) != 2:
            raise ValueError('exact complex coefficient pairs required')
        total = add(total, scale((rational(z[0]), rational(z[1])), rational(value)))
    return total


def sampling_budget(amplification, variance, error, failure_probability, settings):
    """Independent repetitions per real lag, Chebyshev plus a union bound."""
    integer(settings, 1)
    amplification, variance, error, failure_probability = map(
        rational, (amplification, variance, error, failure_probability))
    if amplification < 1 or not 0 < variance <= 1 or error <= 0 or not 0 < failure_probability < 1:
        raise ValueError('positive normalization/error and a probability in(0,1) required')
    needed = settings*amplification**2/(failure_probability*variance**2*error**2)
    repetitions = (needed.numerator+needed.denominator-1)//needed.denominator
    return {'settings': settings, 'repetitions_per_setting': repetitions,
            'total_independent_repetitions': settings*repetitions,
            'complex_statistic_error': error,
            'failure_probability_upper': failure_probability,
            'independent_repetitions_required': True}


def source_budget(level, occupied, time=F(1, 64), order=12, precision_bits=40):
    """Authenticate every finite native experiment or retain nonadmission."""
    integer(level); integer(occupied)
    if level > 2:
        raise ValueError('the concrete source attachment evaluates levels zero through two')
    family, first, second, _ = exchange.sector.probe_data(level)
    ports, seams, degree = exchange.graph_data(family)
    preparation.variance(ports, occupied)
    value = continuation(seams, degree, 2, time, order)
    graph = preparation.mixing.graph_certificate(ports, family.seams)
    prepared = preparation.preparation_budget(ports, occupied,
        F(graph['quantum_gap_lower_per_attempt']), precision_bits)
    maximum_attempts = prepared['warmup_attempts']+value['lag_attempts'][-1]+4
    admitted = maximum_attempts <= exchange.sector.source.interpreter.MAX_ATTEMPTS
    programs = (tuple(preparation.prepared_program(level, occupied,
        prepared['warmup_attempts'], lag) for lag in value['lag_attempts']) if admitted else ())
    systematic = value['coefficient_l1_upper']*F(prepared['normalized_score_error_upper'])
    return {**value, 'level': level, 'ports': ports, 'occupation': occupied,
            'probe_seams': (first, second), 'stationary_variance': preparation.variance(ports, occupied),
            'preparation': prepared, 'maximum_attempts_per_experiment': maximum_attempts,
            'native_single_program_attempt_ceiling': exchange.sector.source.interpreter.MAX_ATTEMPTS,
            'all_setting_programs_admitted': admitted,
            'program_digests': tuple(exchange.sector.source.interpreter.digest(p) for p in programs),
            'complete_record_settings': len(value['lag_attempts']),
            'amplified_preparation_error_upper': systematic,
            'transfer_vacuum_systematic_error_upper': systematic+value['operator_error_upper'],
            'local_exchange_systematic_error_upper': systematic+value['operator_error_upper']+
                value['local_exchange_unitary_error_upper'],
            'full_expanded_histories_executed': False,
            'sampling': sampling_budget(value['coefficient_l1_upper'], preparation.variance(ports, occupied),
                                        F(1, 1024), F(1, 100), len(value['lag_attempts']))}


def exact_native_lags(family, probe, degree, lags):
    """Bounded exact transfer moments; not sampling or expanded execution."""
    ports, seams, _ = exchange.graph_data(family)
    probe, _ = exchange.sector.coefficients(probe, ports, degree)
    if (not probe or exchange.sector.lower(probe, ports, degree)
            or type(lags) is not tuple or any(type(k) is not int or not 0 <= k <= 96 for k in lags)):
        raise ValueError('a nonzero primitive probe and bounded exact lag tuple required')
    norm = exchange.sector.dot(probe, probe); current = probe
    collected = {0: F(1)}
    for k in range(1, max(lags, default=0)+1):
        current = exchange.sector.counts_step(family, current, degree)
        collected[k] = F(exchange.sector.dot(probe, current), norm*(2*seams)**k)
    return tuple(collected[k] for k in lags)


def unitary_taylor(moments, time):
    """An exact complex Taylor value; the source ceiling bounds its remainder."""
    time = rational(time)
    if type(moments) is not tuple or len(moments) < 2:
        raise ValueError('exact normalized spectral moments required')
    values = tuple(rational(x) for x in moments)
    if values[0] != 1 or any(x < 0 for x in values):
        raise ValueError('normalized nonnegative energy moments required')
    answer = (F(), F())
    for k, moment in enumerate(values):
        unit = ((F(1), F(0)), (F(0), F(-1)), (F(-1), F(0)), (F(0), F(1)))[k % 4]
        answer = add(answer, scale(unit, time**k*moment/factorial(k)))
    return answer


def response_witness(family, probe, degree=2, occupied=2, time=F(1, 64), kick=F(1, 4096)):
    """A necessary intervention test for an identity transfer dictionary.

    F is the literal product of two occupation contrasts, so ||F||<=1.
    This routine requires exactly that primitive contrast coefficient pattern.
    """
    ports, seams, maximum_degree = exchange.graph_data(family)
    integer(degree, 1); integer(occupied)
    probe, _ = exchange.sector.coefficients(probe, ports, degree)
    vertices = sorted({x for pair in probe for x in pair})
    if degree != 2 or len(vertices) != 4 or len(probe) != 4 or sorted(probe.values()) != [-1, -1, 1, 1]:
        raise ValueError('a literal four-port contrast product required')
    from itertools import combinations
    if not any(probe == exchange.sector.contrast(a, b, ports)
               for a in combinations(vertices, 2) for b in combinations(vertices, 2) if not set(a)&set(b)):
        raise ValueError('the probe must have the original contrast signs')
    time = rational(time); kick = rational(kick)
    if time <= 0 or kick <= 0:
        raise ValueError('positive exact response time and phase parameter required')
    mu = exchange.moments(family, probe, degree, 3)
    variance = preparation.variance(ports, occupied)
    susceptibility = 2*time*mu[1]-time**3*mu[3]/3
    transfer_allowance = 2*time*exchange.degree_budget(seams, maximum_degree, degree)['generator_error_upper']
    transfer_lower = susceptibility-transfer_allowance
    finite_signal = kick*variance*transfer_lower-2*kick*kick
    if finite_signal <= 0:
        raise ArithmeticError('the exact finite intervention bound has no positive margin')
    return {'time': time, 'phase_kick': kick, 'variance': variance,
            'normalized_generator_moments_0_through_3': mu,
            'local_exchange_normalized_susceptibility_lower': susceptibility,
            'transfer_normalized_susceptibility_lower': transfer_lower,
            'transfer_finite_kick_score_difference_lower': finite_signal,
            'native_bare_record_difference': F(0),
            'native_stationary_working_state_difference_for_any_unitary': F(0),
            'dictionary_requires_phase_kick_and_future_readout_intertwining': True,
            'all_encodings_or_coherent_future_interventions_excluded': False}
