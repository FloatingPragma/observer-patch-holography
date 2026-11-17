"""Recover a finite detector Hilbert space from incident-clock correlations.

The recovery function consumes rational correlations only. Native source
enumeration is a separate certificate of those inputs and their continuation.
Exact rank is a law-level result; sampled records have a separate error budget.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path
import random
from typing import NamedTuple

import sympy as sp

from ..observer import native
from ..observer_time_fields import core
from ..observer_time_fields import spatial_dynamics as waves
from ..native_action import incident_refinement as clocks
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*core.SOURCE_PATHS, *clocks.screen.SOURCE_PATHS, *waves.SOURCE_PATHS,
    'model/native_action/incident_refinement.py',
    'model/relative_observer_geometry/temporal_records.py'))))


def _rational(value):
    return F(int(value.p), int(value.q))


def _matrix_rows(matrix):
    return tuple(tuple(_rational(value) for value in row) for row in matrix.tolist())


def _positive(matrix):
    # Sylvester's criterion is exact here; no numerical rank tolerance occurs.
    return all(matrix[:i, :i].det(method='domain-ge') > 0
               for i in range(1, matrix.rows+1))


class Representation(NamedTuple):
    moments: tuple
    dimension: int
    gram: tuple
    shifted_gram: tuple
    transfer: tuple
    polynomial: tuple  # ascending powers, leading coefficient one
    gram_lower: F


def recover(moments):
    """Recover the cyclic positive transfer from an exact flat moment sequence.

    The declared model class is c_k=<f,T**k f>, with 0<T<I on this cyclic
    space. A zero Schur remainder at lag 2r proves the degree-r recurrence
    under that class. Finite noisy moments do not establish this equality.
    """
    if (type(moments) is not tuple or not 3 <= len(moments) <= 31
            or len(moments) % 2 != 1
            or any(type(value) is not F for value in moments)):
        raise ValueError('an odd bounded tuple of exact rational correlations is required')
    values = [sp.Rational(value.numerator, value.denominator) for value in moments]
    if values[0] <= 0:
        raise ValueError('a nonzero positive field variance is required')
    for rank in range(1, (len(values)+1)//2):
        gram = sp.Matrix(rank, rank, lambda i, j: values[i+j])
        if not _positive(gram):
            raise ValueError('the leading temporal Gram is not positive definite')
        next_column = sp.Matrix(values[rank:2*rank])
        coefficients = gram.inv()*next_column
        remainder = values[2*rank]-(next_column.T*coefficients)[0]
        if remainder < 0:
            raise ValueError('the temporal moments have a negative Schur remainder')
        if remainder:
            continue
        shifted = sp.Matrix(rank, rank, lambda i, j: values[i+j+1])
        if not _positive(shifted) or not _positive(gram-shifted):
            raise ValueError('the recovered transfer is not a strict positive contraction')
        transfer = sp.zeros(rank)
        for i in range(rank-1):
            transfer[i+1, i] = 1
        transfer[:, rank-1] = coefficients
        if gram*transfer != shifted:
            raise ArithmeticError('the record metric does not make transfer self-adjoint')
        for index in range(rank, len(values)):
            if values[index] != sum(coefficients[j]*values[index-rank+j] for j in range(rank)):
                raise ValueError('later supplied correlations violate the flat recurrence')
        return Representation(tuple(F(v.numerator, v.denominator) for v in moments), rank,
            _matrix_rows(gram), _matrix_rows(shifted), _matrix_rows(transfer),
            tuple(_rational(-v) for v in coefficients)+(F(1),),
            _rational(1/sp.trace(gram.inv())))
    raise ValueError('the supplied correlations do not certify a finite flat extension')


def predict(representation, lag):
    """Predict a correlation from a previously checked record representation."""
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    if type(lag) is not int or not 0 <= lag <= 128:
        raise ValueError('an exact correlation lag from zero through 128 is required')
    values = list(representation.moments)
    r = representation.dimension
    while len(values) <= lag:
        values.append(-sum(representation.polynomial[j]*values[-r+j] for j in range(r)))
    return values[lag]


@lru_cache(maxsize=1)
def _mean_kernel():
    source = core.source()
    identity, remote, incident = sp.eye(12), sp.zeros(12), sp.zeros(12)
    for a, b in source.seams:
        swap = sp.eye(12)
        swap[a, a] = swap[b, b] = 0
        swap[a, b] = swap[b, a] = 1
        if {a, b}.intersection(source.region):
            incident += (identity+swap)/source.denominator
        else:
            remote += (identity+swap)/source.denominator
    return incident*(identity-remote).inv()


def source_correlations(max_lag=20):
    """All four two-field joint laws, via the independent full 66-state kernel.

    Only initial local bits and incident signed currents are observer inputs.
    Full-state propagation certifies their exact laws; it is not observer data.
    """
    if type(max_lag) is not int or not 14 <= max_lag <= 32:
        raise ValueError('an exact source lag from fourteen through 32 is required')
    source, kernel = core.source(), core.completed()
    indicators = tuple(tuple(state[port] for state in source.states) for port in (0, 1))
    columns = indicators
    denominator = len(source.states)
    rows = []
    mean_kernel, mean_power = _mean_kernel(), sp.eye(12)
    for lag in range(max_lag+1):
        joint = tuple(tuple(F(sum(a*b for a, b in zip(indicators[i], columns[j], strict=True)),
                                  denominator) for j in (0, 1)) for i in (0, 1))
        for i in (0, 1):
            for j in (0, 1):
                independent = sp.Rational(1, 36)+sp.Rational(5, 33)*(mean_power[i, j]-sp.Rational(1, 12))
                if joint[i][j] != _rational(independent):
                    raise ArithmeticError('full native records and the occupation reduction disagree')
        rows.append(joint)
        if lag != max_lag:
            columns = tuple(tuple(sum(a*b for a, b in zip(row, column, strict=True))
                                  for row in kernel.counts) for column in columns)
            denominator *= kernel.denominator
            mean_power = mean_kernel*mean_power
    return tuple(rows)


def _outward(value, upper=False):
    scale = 10**12
    units = -((-value.numerator*scale)//value.denominator) if upper else value.numerator*scale//value.denominator
    return F(units, scale)


def _negative_log_bounds(value, terms=24):
    # -log(x)=2 sum z^(2j+1)/(2j+1), z=(1-x)/(1+x), 0<x<1.
    z = (1-value)/(1+value)
    lower = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    return lower, lower+2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))


def spectrum(representation):
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    x = sp.Symbol('x')
    polynomial = sp.Poly(sum(sp.Rational(c.numerator, c.denominator)*x**i
                             for i, c in enumerate(representation.polynomial)), x)
    intervals = sp.polys.polytools.intervals(polynomial, eps=sp.Rational(1, 10**12))
    result = []
    for (left, right), multiplicity in intervals:
        left, right = _rational(left), _rational(right)
        if not 0 < left <= right < 1 or multiplicity != 1:
            raise ArithmeticError('the detector polynomial lost its distinct positive roots')
        energy = _negative_log_bounds(right)[0], _negative_log_bounds(left)[1]
        result.append({'transfer_interval': [str(_outward(left)), str(_outward(right, True))],
            'energy_interval_per_incident_tick': [str(_outward(energy[0])), str(_outward(energy[1], True))]})
    return result


def sampling_bound(samples, tolerance, settings=4):
    """Two-sided Hoeffding union bound, with rational exponential majorant."""
    if (type(samples) is not int or not 1 <= samples <= 10**9
            or type(settings) is not int or not 1 <= settings <= 64
            or type(tolerance) is not F or not 0 < tolerance < 1):
        raise ValueError('bounded sample count, settings and rational tolerance required')
    exponent = 2*samples*tolerance*tolerance
    denominator = sum((exponent**j/math.factorial(j) for j in range(65)), F(0))
    return min(F(1), 2*settings/denominator)


def perturbation_budget(representation, error):
    """Uniform correlation error gives a certified generalized-eigenvalue bound.

    If every used correlation changes by <=error, both Hankel errors have
    operator norm <=r*error. With g<=lambda_min(H0), the ordered generalized
    eigenvalues change by <=2*r*error/(g-r*error), when r*error<g.
    This does not certify exact rank or the infinite continuation from samples.
    """
    if type(representation) is not Representation or representation != recover(representation.moments):
        raise ValueError('an authentic recovered temporal representation is required')
    if type(error) is not F or error < 0:
        raise ValueError('a nonnegative exact correlation error is required')
    eta = representation.dimension*error
    if eta >= representation.gram_lower:
        return None
    return 2*eta/(representation.gram_lower-eta)


def sampled_records(seed=24017, samples=2048):
    """Independent reset histories; later fields use incident currents only.

    The seed reproduces executions, not an IID-entropy certificate. The bound
    is conditional on the stated independent uniform preparation/scheduler and
    fair native coin services. There is no postselection or capped timeout.
    """
    if type(seed) is not int or not 0 <= seed < 2**32 or type(samples) is not int or not 1 <= samples <= 8192:
        raise ValueError('a bounded exact seed and sample count required')
    source = core.source()
    rng = random.Random(seed)
    rows = []
    for lag in (1, 2, 3, 4):
        successes, attempts, digest = 0, 0, hashlib.sha256()
        for _ in range(samples):
            state = source.states[rng.randrange(66)]
            initial, counter, ticks = state[0], state[0], 0
            digest.update(bytes((initial, state[1])))
            while ticks < lag:
                seam = source.seams[rng.randrange(30)]
                state, (_, current) = native.native(state, seam, rng.randrange(2))
                attempts += 1
                if set(seam).intersection(source.region):
                    ticks += 1
                    counter += current*(int(seam[1] == 0)-int(seam[0] == 0))
                    if counter not in (0, 1):
                        raise ArithmeticError('authentic incident currents lost the local field')
                    digest.update(bytes((*seam, current+1)))
            # This hidden-state equality is an audit, never an observer input.
            if counter != state[0]:
                raise ArithmeticError('record counter disagrees with the native occupation')
            successes += initial*counter
        rows.append({'lag_in_incident_ticks': lag, 'samples': samples,
            'joint_ones_count': successes, 'native_attempts': attempts,
            'accessible_record_sha256': digest.hexdigest()})
    return rows


def screen_compression():
    """Certified finite temporal subspace of both actual screen processes.

    This uses their Bernoulli density and whole-root clock, not the isolated
    total-two source's spectrum. It is a compression, not a closed model.
    """
    ports, edges, incident, observed, exits, _ = clocks.geometry()
    counts = clocks.exit_counts(ports, edges, observed, exits, 48)
    intervals = [(F(2, 9), F(2, 9))]
    for row in clocks.density_response(ports, edges, incident, 3):
        low, high = map(F, row['return_probability_interval'])
        spatial, _ = clocks.exit_bound(counts, row['scheduled_ticks'], len(incident))
        intervals.append((F(2, 9)*low-spatial, F(2, 9)*high+spatial))
    # The certificate uses the exported intervals themselves, including their
    # additional outward rounding, so consumers can reconstruct its budget.
    intervals = [(_outward(lo), _outward(hi, True)) for lo, hi in intervals]
    centers = [(lo+hi)/2 for lo, hi in intervals]
    error = max((hi-lo)/2 for lo, hi in intervals)
    gram = sp.Matrix(2, 2, lambda i, j: sp.Rational(centers[i+j]))
    shifted = sp.Matrix(2, 2, lambda i, j: sp.Rational(centers[i+j+1]))
    lowers = []
    for matrix, allowance in ((gram, 2*error), (shifted, 2*error), (gram-shifted, 4*error)):
        if not _positive(matrix):
            raise ArithmeticError('the source temporal pencil has no positive midpoint certificate')
        lower = _rational(1/sp.trace(matrix.inv()))-allowance
        if lower <= 0:
            raise ArithmeticError('the full-screen temporal positivity is unresolved by its error budget')
        lowers.append(lower)
    # Min-max perturbation of the generalized Rayleigh quotient. Since both
    # pencils are positive contractions, numerator and denominator errors add.
    spectral_error = 4*error/lowers[0]
    x = sp.Symbol('x')
    polynomial = sp.Poly((shifted-x*gram).det(), x).monic()
    roots = sp.polys.polytools.intervals(polynomial, eps=sp.Rational(1, 10**12))
    spectral = []
    for (left, right), multiplicity in roots:
        if multiplicity != 1:
            raise ArithmeticError('the two temporal Ritz values merged')
        spectral.append([str(_outward(_rational(left)-spectral_error)),
                         str(_outward(_rational(right)+spectral_error, True))])
    return {'source_levels': [3, 5], 'preparation': 'Bernoulli product density 1/3',
        'clock_region_ports': 12, 'incident_seams': len(incident),
        'centered_root_field_covariance_intervals': [[str(_outward(lo)), str(_outward(hi, True))]
                                                    for lo, hi in intervals],
        'temporal_vectors': ['f', 'T f'], 'temporal_gram_dimension_lower': 2,
        'gram_positive_lower': str(lowers[0]), 'shifted_gram_positive_lower': str(lowers[1]),
        'gram_minus_shifted_positive_lower': str(lowers[2]),
        'uniform_covariance_midpoint_error': str(error),
        'ordered_compressed_transfer_eigenvalue_intervals': spectral,
        'generalized_eigenvalue_error_upper': str(spectral_error),
        'intervals_hold_separately_for_either_actual_full_screen': True,
        'compressed_transfer_is_closed_or_equals_isolated_carrier_spectrum': False,
        'all_level_continuum_or_physical_time_claimed': False}


def incident_clock_bound(level, incident_seams, horizon=F(1), tolerance=F(1, 8)):
    """Uniform count/inverse bounds on the simultaneous native density scale.

    The event set is fixed before the IID seam schedule is drawn. Dependence
    between different clocks, endpoint outcomes and field records is retained.
    This finite probability bound uses the exact source exposure, not seconds.
    """
    if (type(level) is not int or not 0 <= level <= 100
            or type(incident_seams) is not int
            or not 1 <= incident_seams <= 630*4**level
            or type(horizon) is not F or horizon <= 0
            or type(tolerance) is not F or tolerance <= 0):
        raise ValueError('a bounded source level, nonempty event set and positive rational bounds required')
    seams, scale = 630*4**level, 1260*16**level
    probability = F(incident_seams, seams)
    ticks = scale*probability
    covered = math.ceil(scale*(horizon+tolerance))
    failure = min(F(1), F(covered)*(1-probability)/(scale*ticks*tolerance**2))
    return {'level': level, 'native_attempt_scale': scale,
        'source_ports': 240*4**level, 'source_seams': seams,
        'incident_seams': incident_seams, 'incident_probability': str(probability),
        'mean_incident_ticks_per_continuum_unit': str(ticks),
        'horizon': str(horizon), 'martingale_tolerance': str(tolerance),
        'covered_native_attempts': covered,
        'uniform_count_error_upper': str(tolerance+F(1, scale)),
        'uniform_inverse_error_upper': str(tolerance+1/ticks),
        'failure_probability_upper': str(failure)}


def predictable_port_clock_bound(level, horizon=F(1), tolerance=F(1, 8)):
    """Uniform bounds for passive access chosen before each fresh seam draw.

    The port can depend on earlier field records. Its integer clock credits
    an incident event by 30/degree, including zero current. Access and degree
    disclosure are supplied; autonomous motion and lab time are not derived.
    """
    fixed = incident_clock_bound(level, 5, horizon, tolerance)
    scale, seams = fixed['native_attempt_scale'], fixed['source_seams']
    b = F(scale, seams)
    return {'level': level, 'native_attempt_scale': scale, 'source_seams': seams,
        'source_ports': fixed['source_ports'], 'allowed_port_degrees': [5, 6],
        'integer_credits_per_incident_event': {'5': 6, '6': 5},
        'conditional_credit_mean_per_attempt': str(F(30, seams)),
        'conditional_credit_variance_upper': str(900*(F(1, 5*seams)-F(1, seams**2))),
        'mean_credits_per_continuum_unit': str(30*b),
        'horizon': str(horizon), 'martingale_tolerance': str(tolerance),
        'covered_native_attempts': fixed['covered_native_attempts'],
        'uniform_count_error_upper': fixed['uniform_count_error_upper'],
        'uniform_inverse_error_upper': str(tolerance+F(1, 5)/b),
        'failure_probability_upper': fixed['failure_probability_upper']}


def predictable_native_clock_control():
    """Complete native first histories and conditional second-step clock laws.

    A nonzero incident current moves readout access from a degree-five port
    to an adjacent degree-six port for the next attempt. Every first history
    is retained. The second clock law uses all actual seams, with both coins;
    its reward is state independent, so histories may be grouped by next port.
    """
    graph = clocks.screen.source.federation(0)
    seams = len(graph.seams)
    first = next(p for p, row in enumerate(graph.neighbours) if len(row) == 5)
    second = next(p for p in graph.neighbours[first] if len(graph.neighbours[p]) == 6)
    initial = tuple(int(p % 3 == 0) for p in range(graph.ports))
    laws = {}
    conditional = []
    for port in (first, second):
        degree = len(graph.neighbours[port])
        law = {}
        for seam in graph.seams:
            reward = (30//degree)*int(port in seam)
            law[reward] = law.get(reward, 0)+2
        laws[port] = law
        mean = sum(F(reward*count, 2*seams) for reward, count in law.items())
        variance = sum(F(count, 2*seams)*(reward-mean)**2 for reward, count in law.items())
        if mean != F(30, seams) or variance != 900*(F(1, seams*degree)-F(1, seams**2)):
            raise ArithmeticError('predictable native port compensator changed')
        conditional.append({'port': port, 'degree': degree,
            'bare_tick_probability': str(F(degree, seams)),
            'credit_mean': str(mean), 'credit_variance': str(variance)})
    first_histories = {}
    accepted_credits = 0
    for seam in graph.seams:
        tick = int(first in seam)
        for coin in (0, 1):
            after, (_, current) = native.native(initial, seam, coin)
            decoded = initial[first]+current*(int(seam[1] == first)-int(seam[0] == first))
            if decoded != after[first]:
                raise ArithmeticError('adaptive access failed the incident-current readback')
            selected = second if tick and current else first
            key = 6*tick, selected
            first_histories[key] = first_histories.get(key, 0)+1
            accepted_credits += 6*int(bool(tick and current))
    totals = {}
    cross = 0
    for (reward, selected), count in first_histories.items():
        for next_reward, next_count in laws[selected].items():
            total = reward+next_reward
            totals[total] = totals.get(total, 0)+count*next_count
            cross += reward*next_reward*count*next_count
    branches = (2*seams)**2
    selected_counts = {str(port): sum(count for (_, p), count in first_histories.items() if p == port)
                       for port in (first, second)}
    if (sum(totals.values()) != branches or min(selected_counts.values()) == 0
            or F(cross, branches) != F(30, seams)**2):
        raise ArithmeticError('complete adaptive clock histories lost the martingale property')
    # Choosing an endpoint after the seam is revealed always produces a tick.
    postselected = sum(F(30, len(graph.neighbours[seam[0]])) for seam in graph.seams)/seams
    return {'source_level': 0, 'source_ports': graph.ports, 'source_seams': seams,
        'initial_occupied_ports': sum(initial), 'initial_port': first, 'next_port_if_nonzero_current': second,
        'selection_uses_only_previous_incident_current': True,
        'first_native_histories': 2*seams, 'second_native_histories': branches,
        'next_port_first_history_counts': selected_counts,
        'conditional_second_step_laws': conditional,
        'two_step_credit_histogram': {str(k): v for k, v in sorted(totals.items())},
        'two_step_credit_mean': str(sum(F(k*v, branches) for k, v in totals.items())),
        'successive_credit_covariance': str(F(cross, branches)-F(30, seams)**2),
        'after_seam_endpoint_selection_credit_mean': str(postselected),
        'count_only_nonzero_current_credit_mean': str(F(accepted_credits, 2*seams)),
        'proper_time_or_autonomous_worldline_identified': False}


def incident_wave_clock_bound(level, degree, band=F(1), horizon=F(1),
                              tolerance=F(1, 8), step_error=F(0)):
    """Attach the admitted coherent walk to zero-inclusive native port ticks.

    A protected wave register receives one fixed W per incident clock event.
    Clock repairs pause during these atomic apparatus services. The bound
    does not describe passive repairs interleaved on the wave register.
    """
    fixed = incident_clock_bound(level, degree, horizon, tolerance)
    if (degree not in (5, 6) or type(band) is not F or band <= 0
            or type(step_error) is not F or step_error < 0):
        raise ValueError('a primitive source degree, positive rational band and nonnegative step error required')
    wave = waves.wave_bounds(level, band, horizon)
    eta = F(fixed['uniform_count_error_upper'])
    beta_squared = F(2*degree**2, 315)
    square = beta_squared*band
    grid = 10**12
    units = math.isqrt(square.numerator*grid**2//square.denominator)
    units += int(F(units**2, grid**2) < square)
    prefactor = F(units, grid)
    ideal = prefactor*(eta+(horizon+eta)*wave['relative_wave_frequency_error'])
    steps = math.floor(F(fixed['mean_incident_ticks_per_continuum_unit'])*(horizon+eta))
    accumulated = steps*step_error
    # (1+epsilon)^N-1 <= N epsilon/(1-N epsilon), when N epsilon<1.
    # This avoids an enormous rational power on growing source levels.
    raw = accumulated/(1-accumulated) if accumulated < 1 else None
    return {'level': level, 'port_degree': degree, 'band_upper': str(band),
        'native_attempt_scale': fixed['native_attempt_scale'],
        'squared_wave_normalization_Q': wave['squared_coherent_clock_normalization'],
        'squared_incident_to_wave_parameter_factor_beta': str(beta_squared),
        'beta_times_sqrt_band_upper': str(prefactor),
        'relative_wave_frequency_error': str(wave['relative_wave_frequency_error']),
        'count_error_eta': str(eta), 'horizon': str(horizon),
        'clock_failure_probability_upper': fixed['failure_probability_upper'],
        'uniform_ideal_compressed_amplitude_error_upper': str(min(F(2), ideal)),
        'maximum_triggered_steps_on_clock_good_event': steps,
        'per_step_operator_error': str(step_error),
        'unitary_step_accumulated_error_upper': str(min(F(2), accumulated)),
        'raw_step_accumulated_error_upper': None if raw is None else str(raw),
        'raw_step_exact_error_formula': '(1+epsilon_l)^N-1',
        'growing_compilation_requirement': 'sqrt(Q_l)*epsilon_l -> 0',
        'protected_wave_register_and_paused_atomic_services_required': True}


def active_wave_clock_control():
    """Exact two-attempt native clock histories driving the existing arc walk.

    The timing carrier and protected coherent wave register are separate.
    The latter starts in the admitted normalized occupation amplitude. All
    clock histories and both Ramsey outcomes are included; no source pulse
    compiler for this whole walk is asserted.
    """
    coherent = waves.coherent
    field = tuple(6*state[0]-1 for state in coherent.STATES)
    norm = sum(value*value for value in field)
    probabilities = {}
    for phase in (-1, 1):
        for ticks in range(3):
            outputs = coherent.ramsey_outputs(field, ticks, phase)
            pair = tuple(coherent.norm_squared(output)/norm for output in outputs)
            if sum(pair) != 1:
                raise ArithmeticError('the triggered coherent walk lost a Ramsey outcome')
            probabilities[phase, ticks] = pair
    source = core.source()
    histories = {}
    for initial in source.states:
        for first_seam in source.seams:
            for first_coin in (0, 1):
                first, (_, first_current) = native.native(initial, first_seam, first_coin)
                first_tick = int(0 in first_seam)
                for second_seam in source.seams:
                    for second_coin in (0, 1):
                        _, (_, second_current) = native.native(first, second_seam, second_coin)
                        second_tick = int(0 in second_seam)
                        ticks = first_tick+second_tick
                        movements = int(bool(first_tick and first_current))+int(bool(second_tick and second_current))
                        key = ticks, movements
                        histories[key] = histories.get(key, 0)+1
    branches = len(source.states)*(2*len(source.seams))**2
    if sum(histories.values()) != branches or norm != 330:
        raise ArithmeticError('the complete native timing-source history changed')
    rows = []
    plus = opposite = accepted = reset = joint = mean_movements = F(0)
    for (ticks, movements), count in sorted(histories.items()):
        weight = F(count, branches)
        pair = probabilities[-1, ticks]
        rows.append({'incident_ticks': ticks, 'nonzero_incident_currents': movements,
            'native_histories': count, 'joint_Ramsey_probabilities': [str(weight*p) for p in pair]})
        plus += weight*pair[0]
        opposite += weight*probabilities[1, ticks][0]
        accepted += weight*probabilities[-1, movements][0]
        reset += weight*probabilities[-1, int(ticks > 0)][0]
        joint += weight*movements*pair[0]
        mean_movements += weight*movements
    if plus == opposite or plus == accepted or plus == reset:
        raise ArithmeticError('an active-wave clock or continuation control disappeared')
    return {'elementary_native_states': len(source.states), 'native_attempt_cap': 2,
        'complete_native_histories': branches, 'admitted_wave_arc_dimension': coherent.DIMENSION,
        'wave_preparation': '(6*n_0-1)/sqrt(330) in the separate coherent arc register',
        'Ramsey_by_incident_ticks': [[str(p) for p in probabilities[-1, ticks]] for ticks in range(3)],
        'joint_clock_current_and_Ramsey_law': rows,
        'complete_Ramsey_plus': str(plus), 'opposite_reflection_Ramsey_plus': str(opposite),
        'nonzero_current_only_trigger_Ramsey_plus': str(accepted),
        'reprepare_before_each_step_Ramsey_plus': str(reset),
        'nonzero_current_count_Ramsey_plus_covariance': str(joint-mean_movements*plus),
        'growing_wave_source_pulse_compiler_executed': False,
        'passive_repairs_act_on_protected_wave_register': False}


def source_clock_incidence(level):
    """Check clock exposure against actual refined native seams and ancestry."""
    if type(level) is not int or level not in (0, 1, 2):
        raise ValueError('an executed source level zero, one or two is required')
    graph = clocks.screen.source.federation(level)
    ancestry = clocks.screen.source.cell_parents(0, level)
    patch = {p for p in range(graph.ports) if ancestry[p//12] == 0}
    patch_events = sum(a in patch or b in patch for a, b in graph.seams)
    boundary = sum((a in patch) != (b in patch) for a, b in graph.seams)
    carrier_events = sum(a//12 == 0 or b//12 == 0 for a, b in graph.seams)
    degree_counts = {degree: sum(len(row) == degree for row in graph.neighbours)
                     for degree in (5, 6)}
    expected = (63*4**level+3*2**level)//2
    if (len(patch) != 12*4**level or boundary != 3*2**level
            or patch_events != expected or carrier_events != 33
            or degree_counts != {5: 180*4**level, 6: 60*4**level}):
        raise ArithmeticError('the actual native incident-clock census changed')
    return {'level': level, 'port_degree_census': {str(k): v for k, v in degree_counts.items()},
        'one_carrier_incident_seams': carrier_events,
        'base_face_ports': len(patch), 'base_face_incident_seams': patch_events,
        'base_face_boundary_seams': boundary}


def clock_fluctuation_covariance(level, first, second, shared):
    """Exact covariance per continuum unit after scaling by sqrt(a_l/M_l).

    Inputs count two actual event sets and their intersection. Fixed bounded
    counts give the Brownian overlap limit; finite source counts retain the
    negative common-scheduler correction, including disjoint event sets.
    """
    incident_clock_bound(level, first)
    incident_clock_bound(level, second)
    seams = 630*4**level
    if (type(shared) is not int
            or not max(0, first+second-seams) <= shared <= min(first, second)):
        raise ValueError('the intersection must define two possible source event sets')
    return ((F(first)-F(first*first, seams), F(shared)-F(first*second, seams)),
            (F(shared)-F(first*second, seams), F(second)-F(second*second, seams)))


def source_clock_fluctuations(level):
    """Actual adjacent, disjoint and cross-carrier port-clock comparisons."""
    graph = clocks.screen.source.federation(level)
    five = [p for p, neighbors in enumerate(graph.neighbours) if len(neighbors) == 5]
    first = five[0]
    adjacent = next(p for p in five if p in graph.neighbours[first])
    disjoint = next(p for p in five if p != first and p not in graph.neighbours[first])
    result = []
    for label, pair in (('adjacent_degree_five', (first, adjacent)),
                        ('nonadjacent_degree_five', (first, disjoint)),
                        ('cross_carrier_degree_six', graph.inter_seams[0])):
        a, b = pair
        da, db = len(graph.neighbours[a]), len(graph.neighbours[b])
        shared = sum(a in seam and b in seam for seam in graph.seams)
        covariance = clock_fluctuation_covariance(level, da, db, shared)
        result.append({'comparison': label, 'ports': [a, b],
            'incidence_counts': [da, db, shared],
            'scaled_count_covariance_per_continuum_unit': [[str(x) for x in row] for row in covariance],
            'limiting_correlation_squared_for_fixed_incidence_counts': str(F(shared**2, da*db))})
    return result


def native_clock_field_dependence():
    """Complete native one-step law: a clock mark and two-time field correlate.

    The small elementary source is a finite dependence control, not a sampled
    execution of the large refined screen or its continuum limit.
    """
    source = core.source()
    clock, field, joint, branches, endpoint, clock_endpoint = 0, 0, 0, 0, 0, 0
    for before in source.states:
        for seam in source.seams:
            tick = int(0 in seam)
            for coin in (0, 1):
                after, (_, current) = native.native(before, seam, coin)
                decoded = before[0]+current*(int(seam[1] == 0)-int(seam[0] == 0))
                if decoded != after[0]:
                    raise ArithmeticError('the clocked current record lost the field')
                value = before[0]*decoded
                clock += tick; field += value; joint += tick*value; branches += 1
                endpoint += decoded; clock_endpoint += tick*decoded
    covariance = F(joint, branches)-F(clock*field, branches**2)
    if branches != 3960 or covariance == 0:
        raise ArithmeticError('the native clock/field dependence control changed')
    return {'elementary_source_states': 66, 'native_branches': branches,
        'clock_mean': str(F(clock, branches)),
        'two_time_occupation_product_mean': str(F(field, branches)),
        'joint_clock_field_mean': str(F(joint, branches)),
        'clock_field_covariance': str(covariance),
        'clock_single_endpoint_covariance': str(F(clock_endpoint, branches)-F(clock*endpoint, branches**2)),
        'clock_field_independence_asserted': False}


def run():
    joint = source_correlations(20)
    centered = tuple(row[0][0]-F(1, 36) for row in joint)
    scalar = recover(centered[:15])
    # Both fields are already readable from the same incident record service.
    plus = tuple(sum(row[i][j]-F(1, 36) for i in (0, 1) for j in (0, 1)) for row in joint)
    minus = tuple(sum((1 if i == j else -1)*(row[i][j]-F(1, 36))
                       for i in (0, 1) for j in (0, 1)) for row in joint)
    symmetric, antisymmetric = recover(plus[:9]), recover(minus[:7])
    if (scalar.dimension, symmetric.dimension, antisymmetric.dimension) != (7, 4, 3):
        raise ArithmeticError('the exact detector-visible mode census changed')
    if any(predict(scalar, k) != centered[k] for k in range(15, 21)):
        raise ArithmeticError('the detector recurrence failed held-out native lags')
    if any(row[0][0] != row[1][1] or row[0][1] != row[1][0] for row in joint):
        raise ArithmeticError('the root exchange symmetry changed')
    # Exact spatial observability rank, independent of the correlation pencil.
    kernel = _mean_kernel()
    q = sp.eye(12)-sp.ones(12)/12
    rows = []
    for port in (0, 1):
        vector = q[:, port].T
        for _ in range(12):
            rows.append(vector)
            vector = vector*kernel
    observability = sp.Matrix.vstack(*rows)
    visible_rank = observability.rank()
    if visible_rank != 7:
        raise ArithmeticError('the calibrated occupation-space visibility changed')
    dark = sp.Matrix.hstack(*sp.Matrix.vstack(observability, sp.ones(1, 12)).nullspace())
    dark_transfer, parameters = dark.gauss_jordan_solve(kernel*dark)
    if dark.cols != 4 or parameters.rows or observability*dark != sp.zeros(24, 4):
        raise ArithmeticError('the centered dark occupation space is not invariant')
    sampled = sampled_records()
    tolerance = F(1, 16)
    for row in sampled:
        predicted = joint[row['lag_in_incident_ticks']][0][0]
        discrepancy = abs(F(row['joint_ones_count'], row['samples'])-predicted)
        row.update(exact_joint_probability=str(predicted), absolute_error=str(discrepancy),
                   within_declared_error=discrepancy <= tolerance)

    def representation(value):
        return {'dimension': value.dimension, 'recurrence_polynomial_ascending': list(map(str, value.polynomial)),
            'gram_minimum_eigenvalue_lower': str(value.gram_lower),
            'transfer_spectrum': spectrum(value)}

    return {'schema': 'oph.detector_temporal_records.v1',
        'source': {'binary_states': 66, 'ports': 12, 'clock_region': [0, 1], 'incident_seams': 9,
            'preparation': 'uniform binary total-two source; both local initial bits disclosed',
            'tick': 'each scheduled seam meeting port zero or one, including zero current',
            'physical_seconds_or_modular_time_identified': False},
        'exact_two_field_joint_ones_by_lag': [[[str(p) for p in row] for row in block] for block in joint],
        'centered_single_field': representation(scalar),
        'symmetric_two_field': representation(symmetric),
        'antisymmetric_two_field': representation(antisymmetric),
        'held_out_lags': list(range(15, 21)), 'held_out_recurrence_matches_full_native_source': True,
        'spatial_visibility': {'centered_occupation_dimension': 11, 'visible_dimension': int(visible_rank),
            'dark_dimension': 11-int(visible_rank), 'rank_three_carrier_band_identified_with_temporal_modes': False,
            'dark_transfer_polynomial_ascending': [str(_rational(c)) for c in
                                                  reversed(dark_transfer.charpoly().all_coeffs())],
            'dark_to_linear_two_field_correlations_not_all_record_protocols': True,
            'physical_spatial_coordinates_recovered': False},
        'finite_record_budget': {'histories': sampled, 'simultaneous_joint_probability_tolerance': str(tolerance),
            'conditional_failure_probability_upper': str(_outward(sampling_bound(2048, tolerance), True)),
            'independent_fair_source_services_assumed_not_certified_by_seed': True,
            'raw_seven_mode_spectral_error_certified_at_this_tolerance': perturbation_budget(scalar, tolerance),
            'certified_spectral_error_at_correlation_error_1e_minus_26': str(
                perturbation_budget(scalar, F(1, 10**26))),
            'finite_samples_prove_exact_rank_or_continuation': False},
        'actual_screen_temporal_compression': screen_compression(),
        'simultaneous_refinement_incident_clocks': {
            'native_attempt_scale': 'a_l=1260*16^l',
            'normalized_clock': 'C_l(t)=N_l(floor(a_l*t))/(a_l*p_l)',
            'inverse_clock': 'U_l(t)=T_l(floor(a_l*p_l*t))/a_l; T_l(0)=0',
            'counting_includes_zero_current': True,
            'actual_source_incidence': [source_clock_incidence(level) for level in (0, 1, 2)],
            'actual_source_clock_fluctuations': [source_clock_fluctuations(level) for level in (0, 1, 2)],
            'finite_bounds': [incident_clock_bound(level, degree)
                              for level in (0, 2, 6) for degree in (5, 6, 33)],
            'complete_native_dependence_control': native_clock_field_dependence(),
            'clock_and_field_independence_required': False,
            'centered_clock_fluctuations_asymptotically_independent_of_density_snapshots': True,
            'independence_scope': 'finite vectors, stationary canonical preparation, seam-only clocks, fixed bounded incidence pattern',
            'joint_inverse_clock_transfer_uses_proved_continuous_density_path_limits': True,
            'density_path_topologies': ['C([0,T],H^-3(S2)) for linear interpolation',
                'C([0,T],H^-2(S2)) after the finite Green square-root readout',
                'D with the same continuous limits for stopped step versions'],
            'local_clock_grants_complete_spatial_field_access': False,
            'modular_parameter_proper_time_or_lab_seconds_identified': False},
        'predictable_port_access': {
            'normalized_clock': 'C_l(t)=R_l(floor(a_l*t))/(30*a_l/M_l), R_l(n)=sum_(j<=n) (30/d_j)*I_j',
            'inverse_clock': 'U_l(t)=T_l(floor(30*a_l/M_l*t))/a_l; T_l(k)=inf{n:R_l(n)>=k}; T_l(0)=0',
            'port_selected_before_fresh_seam_draw': True,
            'previous_field_records_may_select_next_port': True,
            'passive_access_path_and_degree_disclosure_supplied': True,
            'finite_bounds': [predictable_port_clock_bound(level) for level in (0, 2, 6)],
            'complete_native_adaptation_control': predictable_native_clock_control(),
            'same_density_and_green_field_path_limits_after_inverse_clock': True,
            'clock_and_field_independence_required': False,
            'fixed_incidence_gaussian_jitter_claim_extended_to_adaptive_ports': False,
            'autonomous_motion_modular_flow_or_proper_time_identified': False},
        'active_incident_wave_clock': {
            'coherent_wave_steps': 'one admitted W_l per primitive port incident event, including zero current',
            'continuum_compressed_wave': 'cos((2*d/sqrt(630))*tau*sqrt(A)) on the proved bounded band',
            'finite_bounds': [incident_wave_clock_bound(level, degree, step_error=F(1, 10**12))
                             for level in (0, 2, 6) for degree in (5, 6)],
            'complete_native_history_control': active_wave_clock_control(),
            'clock_and_wave_history_independence_required': False,
            'protected_wave_register_paused_services_and_common_readout_supplied': True,
            'continuation_selected_by_passive_repair_or_atomic_clock_identified': False},
        'scope': {'exact_laws_reconstruct_finite_cyclic_positive_transfer': True,
            'positive_transfer_class_is_existing_native_theorem': True,
            'coherent_real_time_implementation_or_lab_duration_derived': False,
            'same_spectrum_on_refined_screens_asserted': False},
        'source_bindings': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Detector temporal Hilbert space and native correlation predictions replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
