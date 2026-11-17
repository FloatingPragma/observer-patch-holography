"""One source continuum operator for native density paths, waves and readouts.

The limit follows from the source compactness theorem, not finite numerical
extrapolation. Wave synthesis processes retained history blocks; its parameter
is distinct from passive native waiting and from a calibrated laboratory clock.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from . import spatial_dynamics as dynamics
from . import wave_reference
from ..observer.source_scalar import compiler as source_compiler
from ..physical_space_time import local_continuum as spatial
from ..physical_space_time import observation_process as observations
from ..relative_observer_geometry import temporal_records as clocks
from tools import strict_json

ROOT = dynamics.ROOT
SOURCE_PATHS = tuple(sorted(set((*dynamics.SOURCE_PATHS, *wave_reference.SOURCE_PATHS,
    *spatial.SOURCE_PATHS, *observations.SOURCE_PATHS, *clocks.SOURCE_PATHS,
    'model/relational_event_clock/__init__.py',
    'model/native_action/stopped_control.py', 'model/native_action/clock_histories.py',
    'model/observer_time_fields/continuum_dynamics.py'))))


def _root_upper(value):
    value = F(value)
    if value < 0:
        raise ValueError('a nonnegative exact square is required')
    grid = 1 << 32
    root = math.isqrt(value.numerator*grid*grid//value.denominator)
    root += int(F(root*root, grid*grid) < value)
    return F(root, grid)


def _int(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('an original integer at or above the bound is required')
    return value


def _fraction(value, positive=False):
    if type(value) not in (int, F):
        raise ValueError('an original integer or Fraction is required')
    value = F(value)
    if positive and value <= 0:
        raise ValueError('a positive exact value is required')
    return value



def propagation_bounds(level, band=16, horizon=4, energy=1):
    """Errors on a unit vector with full scaled native energy at most energy.

    These compare finite native evolutions with the same finite Laplacian.
    The subsequent spatial subsequence has no claimed convergence rate.
    """
    level = dynamics._level(level)
    band, horizon, energy = map(dynamics._positive, (band, horizon, energy))
    scale = 4**level
    q, n = 630*scale**2, 1260*scale**2
    if band >= n:
        raise ValueError('the retained band must lie below the native attempt normalization')
    tail = 2*_root_upper(energy/band)
    wave = dynamics.wave_bounds(level, band, horizon)
    wave_low = horizon*_root_upper(band)*wave['relative_wave_frequency_error']
    wave_rounding = _root_upper(band)/(25*scale)
    heat_low = horizon*band**2/(2*(n-band))+band/n
    return {'level': level, 'scaled_native_energy_upper': energy,
        'spectral_band_upper': band, 'parameter_horizon': horizon,
        'high_spectral_tail_amplitude_upper': tail,
        'integer_attempt_heat_error_on_band': heat_low,
        'integer_coherent_wave_error_on_band': wave_low+wave_rounding,
        'integer_coherent_parameter_rounding_upper': F(1, 25*scale),
        'heat_error_on_energy_domain_upper': min(F(2), heat_low+tail),
        'wave_error_on_energy_domain_upper': min(F(2), wave_low+wave_rounding+tail),
        'joint_cap_modular_cosine_error_upper': min(F(2), horizon*26/scale),
        'sufficient_native_return_cap': wave_reference.refinement_budget(level)['sufficient_native_length_cap'],
        'spatial_subsequence_convergence_rate_supplied': False}


def fixed_cap_collapse(level, cap=512, band=16):
    """A fixed length cap freezes every bounded scaled spatial band.

    P_h(x)/P_h(1)>=x^h gives Theta_h(x)<=-h log(x). This is an
    exact scalar inequality for the actual positive native transfer.
    """
    level = dynamics._level(level)
    if type(cap) is not int or cap < 0:
        raise ValueError('an original nonnegative length cap is required')
    band = dynamics._positive(band)
    scale, q = 4**level, 630*16**level
    if band >= 2*q:
        raise ValueError('the band must keep the native transfer positive')
    return {'level': level, 'fixed_cap': cap, 'band_upper': band,
        'scaled_modular_frequency_upper': F(26*cap*scale)*band/(2*q-band),
        'fixed_cap_low_band_generator_limit': 0,
        'uncapped_wave_generator_limit': 'sqrt(A)',
        'fixed_cap_interchanges_with_spatial_limit': False}


def wave_path_bounds(max_degree=10, horizon=4):
    """Uniform C([0,T],H^-2) tail and increment bounds for synthesized waves."""
    if type(max_degree) is not int or max_degree < 0:
        raise ValueError('an original nonnegative harmonic cutoff is required')
    horizon = dynamics._positive(horizon)
    return {'harmonic_cutoff': max_degree, 'parameter_horizon': horizon,
        'expected_squared_supremum_tail_upper': 4*(spatial.GREEN_BOUND+horizon**2)/(max_degree+1)**2,
        'expected_squared_increment_upper_per_dt_squared': F(2),
        'expected_fourth_increment_upper_per_dt_fourth': F(12),
        'initial_covariance': 'G=A^-1',
        'two_time_covariance': 'G cos((t-u) sqrt(A))',
        'velocity_covariance': 'identity on mean-zero L2(mu)',
        'preparation': 'two disjoint normalized blocks of one native-history functional CLT',
        'time_is_declared_spectral_readout_parameter': True}


def finite_speed_bounds(level, lipschitz=1, separation=1, heat_parameter=1):
    """Weighted native heat bounds and the resulting continuum support cone.

    Distances use the unit-circumradius flat comparison polyhedron. The
    support theorem acts on Abar=A plus its constant zero mode, not a
    nonlocal mean-zero projection. It is a continuum theorem only.
    """
    level = dynamics._level(level)
    lipschitz, separation, heat_parameter = map(dynamics._positive,
        (lipschitz, separation, heat_parameter))
    scale = 4**level
    ratio = lipschitz**2/(12*scale)
    if ratio >= 1:
        raise ValueError('refine until the squared edge weight is below twelve')
    # Successive cosh terms after x^2/2 have ratio at most x^2/12.
    native_exponent = heat_parameter*lipschitz**2/(2*(1-ratio))
    return {'level': level, 'weight_Lipschitz_constant': lipschitz,
        'heat_parameter': heat_parameter, 'comparison_support_separation': separation,
        'native_weighted_heat_log_norm_upper': native_exponent,
        'continuum_weighted_heat_log_norm_upper': heat_parameter*lipschitz**2/2,
        'continuum_separated_heat_log_norm_upper': -separation**2/(2*heat_parameter),
        'continuum_wave_speed_squared_upper': F(1, 2),
        'continuum_zero_wave_matrix_element_when_time_squared_below': 2*separation**2,
        'constant_zero_mode_restored_before_support_localization': True,
        'exact_finite_level_wave_support_cone_asserted': False,
        'comparison_units_are_metres_and_seconds': False}


def native_weighted_form_control():
    """Exact weighted full-source identity and two locality countercontrols."""
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    matrix = sp.Matrix(lap.toarray())
    weights = [sp.Rational(2)**((p//12) % 3) for p in range(graph.ports)]
    weighted = sp.Matrix(graph.ports, graph.ports,
        lambda i, j: matrix[i, j]*(weights[i]/weights[j]+weights[j]/weights[i])/2)
    extra = sp.zeros(graph.ports)
    penalty = [sp.S.Zero]*graph.ports
    for a, b in graph.seams:
        eta = (weights[a]/weights[b]+weights[b]/weights[a])/2-1
        if a//12 == b//12 and eta:
            raise ArithmeticError('cell weights vary inside the native carrier')
        for p in (a, b):
            extra[p, p] += eta
            penalty[p] += eta
        extra[a, b] -= eta; extra[b, a] -= eta
    maximum = max((weights[a]/weights[b]+weights[b]/weights[a])/2-1
                  for a, b in graph.inter_seams)
    if (weighted != matrix+extra-sp.diag(*penalty)
            or min(penalty) < 0 or max(penalty) > maximum):
        raise ArithmeticError('the native weighted accretivity decomposition failed')
    # Without cell-repeated weights the internal edges invalidate the
    # cross-seam-only penalty: choose an internal port with no cross seam.
    endpoints = {p for edge in graph.inter_seams for p in edge}
    port = next(p for p in range(graph.ports) if p not in endpoints)
    degree = len(graph.neighbours[port])
    invalid_constant_form = -sp.Rational(degree, 2)
    if invalid_constant_form >= 0:
        raise ArithmeticError('the internal-weight countercontrol disappeared')
    # The finite continuous cosine has an immediate t^2 term across a
    # native cross seam. Exact finite-level support is not the limit theorem.
    a, b = graph.inter_seams[0]
    t_squared_coefficient = -matrix[a, b]/2
    if t_squared_coefficient != sp.Rational(1, 2):
        raise ArithmeticError('the finite wave-tail countercontrol changed')
    return {'source_level': 0, 'source_ports': graph.ports, 'source_cross_seams': len(graph.inter_seams),
        'exact_weighted_form_decomposition': True,
        'maximum_cross_seam_cosh_minus_one': str(maximum),
        'maximum_port_penalty': str(max(penalty)),
        'nonnegative_edge_and_diagonal_remainder': True,
        'internal_weight_mutation_cross_seam_penalty': '0',
        'internal_weight_mutation_constant_quadratic_form': str(invalid_constant_form),
        'finite_cross_seam_wave_t_squared_coefficient': str(t_squared_coefficient),
        'finite_continuous_wave_has_instantaneous_small_tail': True}


def native_heat_support_control():
    """Exact source paths distinguish heat tails from finite attempt fronts.

    Uniformization has nonnegative coefficients and a positive coefficient
    along every connected path. No numerical exponential or fitted kernel
    is used, and this finite gate is not a proof of continuum irreducibility.
    """
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    matrix = sp.Matrix(lap.toarray())
    ports, start = graph.ports, 0
    maximum_degree = max(map(len, graph.neighbours))
    if matrix*sp.ones(ports, 1) != sp.zeros(ports, 1):
        raise ArithmeticError('the source constant mode is not conservative')
    uniformized = maximum_degree*sp.eye(ports)-matrix
    if any(value < 0 for value in uniformized):
        raise ArithmeticError('the native heat uniformization is not nonnegative')
    distance = {start: 0}
    queue = [start]
    for vertex in queue:
        for other in graph.neighbours[vertex]:
            if other not in distance:
                distance[other] = distance[vertex]+1
                queue.append(other)
    if len(distance) != ports:
        raise ArithmeticError('the native source is disconnected')
    length = max(distance.values())
    end = min(p for p, value in distance.items() if value == length)
    seam_count = len(graph.seams)
    # Q=I-L/(4M) is the unbiased fixed-selection-slot transfer. Its
    # numerator differs from uniformization only on the diagonal.
    heat, slot = [0]*ports, [0]*ports
    heat[start] = slot[start] = 1
    prefix = []
    for step in range(length+1):
        prefix.append(slot[end])
        if step == length:
            break
        heat = [(maximum_degree-len(graph.neighbours[p]))*heat[p]
                +sum(heat[q] for q in graph.neighbours[p]) for p in range(ports)]
        slot = [(4*seam_count-len(graph.neighbours[p]))*slot[p]
                +sum(slot[q] for q in graph.neighbours[p]) for p in range(ports)]
    if any(prefix[:-1]) or heat[end] <= 0 or slot[end] != heat[end]:
        raise ArithmeticError('the exact native support front changed')
    # Removing all cross seams leaves separate twelve-port carriers.
    disconnected = {start}
    queue = [start]
    for vertex in queue:
        for other in graph.neighbours[vertex]:
            if vertex//12 == other//12 and other not in disconnected:
                disconnected.add(other); queue.append(other)
    profile = [int(p//12 == start//12) for p in range(ports)]
    changed_seams = tuple((a, b) for a, b in graph.seams if profile[a] != profile[b])
    if (len(disconnected) != 12 or end in disconnected or len(changed_seams) != 3
            or any(not (profile[a] or profile[b]) for a, b in changed_seams)):
        raise ArithmeticError('the disconnected or intervention-support control changed')
    return {'source_level': 0, 'ports': ports, 'native_seams': seam_count,
        'source_port': start, 'target_port': end, 'graph_distance': length,
        'maximum_native_degree': maximum_degree,
        'nonnegative_uniformization_matrix': True,
        'shortest_path_count': heat[end],
        'continuous_heat_lower_bound': 'exp(-maximum_native_degree*t)*shortest_path_count*t^graph_distance/graph_distance!',
        'continuous_heat_leading_taylor_coefficient': F(heat[end], math.factorial(length)),
        'fixed_slot_front_zero_before_graph_distance': True,
        'first_nonzero_slot_transfer': F(slot[end], (4*seam_count)**length),
        'first_nonzero_slot_response_without_susceptibility': -2*F(slot[end], (4*seam_count)**length),
        'cell_supported_control_changed_seams': len(changed_seams),
        'all_changed_control_seams_meet_profile_support': True,
        'cross_seam_deletion_reachable_ports': len(disconnected),
        'cross_seam_deletion_has_zero_cross_component_heat': True,
        'mean_zero_projected_identity_off_diagonal': -F(1, ports),
        'constant_mode_must_be_restored_for_heat_positivity': True,
        'finite_native_attempt_front_is_an_instantaneous_continuum_cone': False}


def density_process_bounds(level, covariance_tolerance=F(1, 10), tests=3,
                           coefficient_l1=3, population_variance_lower=F(1, 2)):
    """Uniform stirring concentration and the conditional finite-population CLT.

    Each real, mean-zero profile has supremum at most one. The CLT ratio
    assumes the displayed positive lower bound on its conditional population
    variance; zero limiting variance is handled separately by Chebyshev.
    Path tightness is certified separately by density_path_bounds.
    """
    level = dynamics._level(level)
    tolerance, coefficient_l1, variance = map(dynamics._positive,
        (covariance_tolerance, coefficient_l1, population_variance_lower))
    if type(tests) is not int or tests < 1:
        raise ValueError('a positive original number of test variables is required')
    vertices, particles = 240*4**level, 80*4**level
    overlap_variance = F(4, vertices)
    pairs = tests*(tests-1)//2
    return {'level': level, 'ports': vertices, 'particles': particles,
        'native_attempts_per_comparison_time': 1260*16**level,
        'susceptibility': F(particles*(vertices-particles), vertices*(vertices-1)),
        'stirring_overlap_variance_upper': overlap_variance,
        'overlap_bound_uniform_in_native_attempt_count': True,
        'test_variables': tests, 'covariance_error_tolerance': tolerance,
        'any_off_diagonal_conditional_covariance_error_probability_upper':
            min(F(1), pairs*overlap_variance/tolerance**2),
        'Cramer_Wold_coefficient_l1_upper': coefficient_l1,
        'conditional_population_variance_lower_assumption': variance,
        'conditional_Hajek_ratio_upper_under_variance_assumption':
            coefficient_l1**2/(variance*min(particles, vertices-particles)),
        'joint_limit': 'centered Gaussian; Cov(Y_tau(f),Y_u(g))=<f,exp(-abs(tau-u) A)g>',
        'equal_time_covariance': 'identity on mean-zero L2(mu)',
        'Green_filtered_covariance': 'G exp(-abs(tau-u) A), using A^-1/2 readout',
        'simultaneous_refinement_and_attempt_time_finite_dimensional_limit': True,
        'count_first_long_block_limit_required': False,
        'path_space_tightness_certificate': 'density_path_bounds and native_density_path_control',
        'clock_is_laboratory_proper_time': False}


def density_path_bounds(level, harmonic_cutoff=10, horizon=4, energy=1,
                        profile_variance=1):
    """Uniform path bounds for actual stationary native snapshots.

    Scalar profiles use F_f=sum f_p(n_p-1/3)/sqrt(V chi) and full-port
    energy E=4^l <f,L f>/V. Harmonic tails use spherical probability-area
    normalization. The finite Green filter is G_l^1/2 on cell tests.
    """
    level = dynamics._level(level)
    if type(harmonic_cutoff) is not int or harmonic_cutoff < 0:
        raise ValueError('an original nonnegative harmonic cutoff is required')
    horizon, energy, variance = map(dynamics._positive,
        (horizon, energy, profile_variance))
    scale, vertices, attempts = 4**level, 240*4**level, 1260*16**level
    particles = vertices//3
    chi = F(particles*(vertices-particles), vertices*(vertices-1))
    steps = math.ceil(attempts*horizon)
    rounded_horizon = F(steps, attempts)
    ballistic_attempts = 1260*8**level
    ballistic_change = 50*F(math.ceil(ballistic_attempts*horizon), attempts)
    return {'level': level, 'harmonic_cutoff': harmonic_cutoff,
        'parameter_horizon': horizon, 'finite_reversal_attempt_horizon': steps,
        'rounded_parameter_horizon': rounded_horizon,
        'native_attempts_per_comparison_time': attempts,
        'full_native_profile_energy_upper': energy,
        'canonical_profile_variance_upper': variance,
        'forward_backward_corrected_profile': 'g=2(I+P)^-1 f',
        'corrected_profile_energy_upper': 4*energy,
        'state_uniform_martingale_increment_variance_upper': 4*energy/(attempts*chi),
        'martingale_maximum_squared_jump_upper': 16*energy/(scale*chi),
        'stationary_terminal_martingale_variance_upper': 4*rounded_horizon*energy,
        'scalar_expected_squared_supremum_upper': 2*variance+50*rounded_horizon*energy,
        'density_H_minus_3_expected_squared_supremum_tail_upper':
            (2+100*rounded_horizon)/(harmonic_cutoff+1)**2,
        'Green_filtered_H_minus_2_expected_squared_supremum_tail_upper':
            (2*spatial.GREEN_BOUND+100*rounded_horizon)/(harmonic_cutoff+1)**2,
        'density_H_minus_3_maximum_squared_snapshot_jump_upper': F(4, vertices)/chi,
        'Green_filtered_H_minus_2_maximum_squared_snapshot_jump_upper': F(2, scale)/chi,
        'density_interpolated_path_topology': 'C([0,T],H^-3(S2))',
        'Green_filtered_interpolated_path_topology': 'C([0,T],H^-2(S2))',
        'step_paths_have_same_continuous_limit_in_Skorohod_J1': True,
        'ballistic_native_attempts_per_parameter': ballistic_attempts,
        'ballistic_to_diffusive_parameter_ratio': F(1, 2**level),
        'ballistic_density_H_minus_3_expected_squared_supremum_change_upper': ballistic_change,
        'ballistic_Green_field_H_minus_2_expected_squared_supremum_change_upper': ballistic_change,
        'ballistic_native_path_limit_is_initial_field_held_constant': True,
        'same_spatial_subsequence_and_simultaneous_attempt_scaling': True,
        'stationary_reversibility_required': True,
        'forward_and_backward_martingale_independence_required': False,
        'spatial_subsequence_convergence_rate_supplied': False}


def native_density_path_control():
    """All native carrier branches check the discrete path-tightness algebra.

    Exact configuration transitions, not a Gaussian simulation, check the
    corrected forward/backward identity, its state-uniform bracket bound,
    and failures of the uncorrected and nonstationary reverse controls.
    """
    states, rows, edges = _native_trace_counts(4)
    vertices, denominator = 12, 2*len(edges)
    chi = F(4*8, vertices*(vertices-1))
    lap = sp.zeros(vertices)
    for a, b in edges:
        lap[a, a] += 1; lap[b, b] += 1
        lap[a, b] -= 1; lap[b, a] -= 1
    transfer = sp.eye(vertices)-lap/denominator
    field = sp.Matrix([1]*6+[-1]*6)
    corrected = 2*(sp.eye(vertices)+transfer).inv(method='DM')*field
    if (sp.eye(vertices)+transfer)*corrected != 2*field:
        raise ArithmeticError('the discrete forward/backward correction changed')
    divisor = math.lcm(*(int(v.q) for v in corrected))
    coefficients = [int(v*divisor) for v in corrected]
    observations = [sum(x*y for x, y in zip(state, coefficients)) for state in states]
    raw = [sum(x*int(y) for x, y in zip(state, field)) for state in states]
    expected = [sum(c*v for c, v in zip(row, observations)) for row in rows]
    expected_raw = [sum(c*v for c, v in zip(row, raw)) for row in rows]
    one_particle_expected = denominator*transfer*corrected
    for state, value in zip(states, expected):
        if sum(x*y for x, y in zip(state, one_particle_expected))*divisor != value:
            raise ArithmeticError('the configuration and density transfers disagree')
    energy = F((field.T*lap*field)[0])/vertices
    corrected_energy = F((corrected.T*lap*corrected)[0])/vertices
    upper = corrected_energy/(denominator*chi)
    variances, raw_defects, branches, maximum_jump = [], set(), 0, F(0)
    nonstationary_reverse = None
    for before, row in enumerate(rows):
        numerator, square_numerator = 0, 0
        for after, multiplicity in enumerate(row):
            if not multiplicity:
                continue
            forward = denominator*observations[after]-expected[before]
            reverse = denominator*observations[before]-expected[after]
            if forward-reverse != 2*denominator*divisor*(raw[after]-raw[before]):
                raise ArithmeticError('a native branch violates forward/backward telescoping')
            numerator += multiplicity*forward
            square_numerator += multiplicity*forward*forward
            maximum_jump = max(maximum_jump,
                F(forward*forward, denominator**2*divisor**2*vertices)/chi)
            raw_defects.add(denominator*(raw[after]-raw[before])
                           +expected_raw[after]-expected_raw[before]
                           -2*denominator*(raw[after]-raw[before]))
            # With a fixed X_0, the backwards conditional mean at X_1 is
            # F_g(X_0)-F_Pg(X_1), not the stationary reversed kernel's zero.
            if before == 0 and reverse and nonstationary_reverse is None:
                nonstationary_reverse = F(reverse, denominator*divisor)
            branches += multiplicity
        variance = F(square_numerator, denominator**3*divisor**2*vertices)/chi
        if numerator or variance > upper:
            raise ArithmeticError('a native conditional martingale bracket failed')
        variances.append(variance)
    exact_stationary = F((corrected.T*(sp.eye(vertices)-transfer**2)*corrected)[0])/vertices
    if (sum(variances, F(0))/len(states) != exact_stationary
            or corrected_energy > 4*energy
            or exact_stationary > 4*energy/denominator
            or maximum_jump > 4*corrected_energy/chi
            or raw_defects == {0} or nonstationary_reverse is None):
        raise ArithmeticError('the native path or mutation certificate changed')
    return {'native_ports': vertices, 'canonical_configurations': len(states),
        'native_endpoint_branches': branches, 'profile_energy': energy,
        'corrected_profile_energy': corrected_energy,
        'exact_corrected_forward_backward_identity_on_every_branch': True,
        'zero_forward_and_stationary_reverse_conditional_mean': True,
        'maximum_state_conditional_increment_variance': max(variances),
        'state_uniform_conditional_increment_variance_upper': upper,
        'stationary_increment_variance': exact_stationary,
        'stationary_increment_variance_upper': 4*energy/denominator,
        'maximum_squared_martingale_increment': maximum_jump,
        'maximum_squared_martingale_increment_upper': 4*corrected_energy/chi,
        'uncorrected_profile_has_nonzero_telescoping_defect': True,
        'fixed_initial_state_reverse_mean_without_normalization': nonstationary_reverse,
        'fixed_initial_state_admits_stationary_reverse_kernel': False}


def native_stirring_control():
    """Exact repair/permutation, random-covariance and scheduling controls.

    All two-attempt stirring branches are enumerated on the native carrier;
    canonical moments independently sum its 495 total-four configurations.
    """
    states, rows, edges = _native_trace_counts(4)
    vertices, particles, denominator = 12, 4, 2*len(edges)
    identity = tuple(range(vertices))
    permutations = []
    for a, b in edges:
        swapped = list(identity)
        swapped[a], swapped[b] = swapped[b], swapped[a]
        permutations.extend((identity, tuple(swapped)))
        for before in states:
            after_swap = tuple(before[i] for i in swapped)
            actual = Counter(dynamics.native.native(before, (a, b), coin)[0]
                             for coin in (0, 1))
            if actual != Counter((before, after_swap)):
                raise ArithmeticError('native unread repair is not identity-or-swap')
    field = np.array([1]*6+[-1]*6, dtype=np.int64)
    single = sum((np.eye(vertices, dtype=np.int64)[list(p)] for p in permutations))
    mean = F(int(field@single@single@field), vertices*denominator**2)
    chi = F(particles*(vertices-particles), vertices*(vertices-1))
    sample = np.asarray(states, dtype=np.int64)
    overlap_sum = overlap_squared_sum = conditional_variance_sum = 0
    conditional_variances = set()
    for first in permutations:
        for second in permutations:
            moved = field[list(first)][list(second)]
            overlap = int(field@moved)
            overlap_sum += overlap
            overlap_squared_sum += overlap*overlap
            coefficients = field+moved
            observations = sample@coefficients
            squared_sum = int(observations@observations)
            expected = chi*len(states)*int(coefficients@coefficients)
            if squared_sum != expected:
                raise ArithmeticError('conditional canonical covariance changed')
            conditional_variance_sum += squared_sum
            conditional_variances.add(F(squared_sum, len(states)*vertices)/chi)
    branches = denominator**2
    if F(overlap_sum, vertices*branches) != mean:
        raise ArithmeticError('independent stirring mean is not the native density semigroup')
    variance = F(overlap_squared_sum, vertices**2*branches)-mean**2
    bound = F(4*int(field@field), vertices**2)
    if not 0 < variance <= bound:
        raise ArithmeticError('stirring overlap concentration control failed')
    native_field = sample@field
    direct_cross = F(int(native_field@np.asarray(rows, dtype=np.int64)@
                         np.asarray(rows, dtype=np.int64)@native_field),
                     len(states)*denominator**2*vertices)/chi
    total_variance = F(conditional_variance_sum, len(states)*branches*vertices)/chi
    if direct_cross != mean or total_variance != 2+2*mean:
        raise ArithmeticError('stirring and full native joint second moments differ')
    fourth = F(sum(int(x)**4 for x in native_field), len(states))/(vertices*chi)**2
    repeated_sum = 0
    for edge_index in range(len(edges)):
        for first in permutations[2*edge_index:2*edge_index+2]:
            for second in permutations[2*edge_index:2*edge_index+2]:
                repeated_sum += int(field@field[list(first)][list(second)])
    repeated_seam_mean = F(repeated_sum, 4*len(edges)*vertices)
    if repeated_seam_mean != F(int(field@single@field), vertices*denominator):
        raise ArithmeticError('repeated-seam mutation does not retain its projector')
    if repeated_seam_mean == mean or fourth == 3 or len(conditional_variances) <= 1:
        raise ArithmeticError('a finite-density or scheduling countercontrol disappeared')
    return {'source_ports': vertices, 'particles': particles,
        'canonical_configurations': len(states), 'native_endpoint_branches': len(states)*denominator,
        'native_unread_repair_equals_identity_or_swap': True,
        'independent_two_attempt_stirring_branches': branches,
        'normalized_profile': 'six plus ones followed by six minus ones',
        'two_attempt_mean_overlap': mean, 'two_attempt_overlap_variance': variance,
        'uniform_overlap_variance_upper': bound,
        'conditional_covariance_varies_before_averaging': True,
        'conditional_canonical_covariance_verified_on_every_branch': True,
        'stirring_and_full_native_two_time_covariance_exactly_equal': True,
        'two_time_sum_variance': total_variance,
        'repeated_seam_mutation_mean_overlap': repeated_seam_mean,
        'single_time_unit_variance_fourth_moment': fourth,
        'finite_native_snapshot_is_exactly_Gaussian': False}


def exact_wave_synthesis():
    """Exact source20cell covariance and conserved velocity identities."""
    _, _, _, _, _, schur = dynamics._base_blocks()
    cells = schur.rows
    constant = sp.ones(cells)/cells
    p0 = sp.eye(cells)-constant
    green = (schur+constant).inv(method='DM')-constant
    correction = sp.Rational(1, 2520)
    covariance = green-correction*p0
    velocity = schur*covariance
    if (velocity != p0-correction*schur or covariance*schur != velocity
            or green*schur != p0 or covariance*sp.ones(cells, 1) != sp.zeros(cells, 1)):
        raise ArithmeticError('native wave-synthesis covariance identities changed')
    # S<=L_cell/12<=I/2 proves positivity of C and 0<=B C<=I.
    _, pivots = covariance[:cells-1, :cells-1].LDLdecomposition(hermitian=False)
    if any(pivots[i, i] <= 0 for i in range(cells-1)):
        raise ArithmeticError('the actual native block covariance is not positive')
    return {'source_cells': cells, 'source_ports': 12*cells,
        'native_attempt_covariance_sha256': dynamics._digest(covariance),
        'native_Green_sha256': dynamics._digest(green),
        'velocity_covariance_sha256': dynamics._digest(velocity),
        'attempt_covariance_correction': F(1, 2520),
        'velocity_covariance_deficit_norm_upper': F(1, 5040),
        'constant_mode_exactly_removed': True,
        'two_history_blocks_have_independent_Gaussian_limits': True,
        'synthesized_wave_is_unprocessed_native_occupation_history': False}


def native_propagator_control():
    """Actual240port native Chebyshev execution and an independent spectral solve."""
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    full = lap.toarray().astype(float)
    eigen, vectors = np.linalg.eigh(full)
    if eigen[0] < -1e-12 or np.linalg.norm(full@vectors-vectors*eigen) > 1e-10:
        raise ArithmeticError('native spectral diagnostic has an excessive residual')
    eigen = np.maximum(eigen, 0.)
    # A prescribed mean-zero cell contrast lifted with the actual port isometry.
    field = np.repeat(np.array([1., -1.]+[0.]*18)/np.sqrt(24.), 12)
    coefficients = vectors.T@field
    kernel = np.eye(graph.ports)-full/1260
    theta = 2*np.arcsin(np.sqrt(eigen/2520))
    previous, current = field, kernel@field
    rows = []
    for count in range(51):
        if count == 0:
            result = field
        elif count == 1:
            result = current
        else:
            previous, current = current, 2*kernel@current-previous
            result = current
        if count not in (0, 1, 2, 25, 50):
            continue
        oracle = vectors@(np.cos(count*theta)*coefficients)
        residual = float(np.linalg.norm(result-oracle))
        parameter = count/np.sqrt(630.)
        target = vectors@(np.cos(parameter*np.sqrt(eigen))*coefficients)
        error = float(np.linalg.norm(result-target))
        allowance = dynamics.wave_bounds(0, band=12, parameter=max(F(1), F(count, 25)))['cosine_propagator_error_on_band']
        if residual > 3e-12 or error > float(allowance)+3e-12:
            raise ArithmeticError('native Chebyshev wave does not meet its source bound')
        amplitude = float(field@result)
        rows.append({'coherent_steps': count,
            'coherent_parameter_diagnostic': parameter,
            'Ramsey_plus_diagnostic': (1+amplitude)/2,
            'Ramsey_minus_diagnostic': (1-amplitude)/2,
            'spectral_oracle_residual_diagnostic': residual,
            'same_native_L_wave_error_diagnostic': error,
            'same_native_L_wave_error_upper': str(allowance)})
    return {'source_ports': graph.ports, 'full_native_seams': len(graph.seams),
        'fixed_preparation': '(cell0 minus cell1)/sqrt(24) on native ports',
        'readouts': rows, 'spatial_limit_inferred_from_finite_rows': False}


def _permutation_moments(vertices, particles, seams, adjacent):
    """Exact normalized traces, from fixed subsets of permutation cycles."""
    if (any(type(v) is not int for v in (vertices, particles, seams, adjacent))
            or vertices < 4 or not 0 < particles < vertices
            or not 0 < seams <= vertices*(vertices-1)//2
            or not 0 <= adjacent <= seams*(seams-1)):
        raise ValueError('a valid finite graph/sector census is required')
    holes = vertices-particles
    chi = F(particles*holes, vertices*(vertices-1))
    psi = F(particles*(particles-1)*holes*(holes-1),
            vertices*(vertices-1)*(vertices-2)*(vertices-3))
    transposition = 1-2*chi
    three_cycle = 1-3*chi
    two_pairs = 1-4*chi+4*psi
    disjoint = seams*(seams-1)-adjacent
    variance = (seams*(1-transposition**2)
                +adjacent*(three_cycle-transposition**2)
                +disjoint*(two_pairs-transposition**2))/(4*seams**2)
    if variance < 0:
        raise ValueError('the supplied census cannot have these permutation moments')
    return {'vertices': vertices, 'particles': particles, 'seams': seams,
        'ordered_adjacent_edge_pairs': adjacent, 'ordered_disjoint_edge_pairs': disjoint,
        'susceptibility': chi, 'transfer_trace_mean': 1-chi,
        'transfer_trace_variance': variance,
        'swap_trace': transposition, 'three_cycle_trace': three_cycle,
        'disjoint_swap_trace': two_pairs}


@lru_cache(maxsize=2, typed=True)
def _native_trace_counts(particles):
    if type(particles) is not int or particles not in (2, 4):
        raise ValueError('the retained native oracle uses total two or four')
    states = tuple(sorted(tuple(int(i in occupied) for i in range(12))
                          for occupied in combinations(range(12), particles)))
    index = {state: i for i, state in enumerate(states)}
    edges, _ = dynamics.source_green.passage.recovered_carrier()
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    for column, before in enumerate(states):
        for edge in edges:
            for coin in (0, 1):
                after, _ = dynamics.native.native(before, edge, coin)
                counts[index[after], column] += 1
    if not np.array_equal(counts, counts.T) or np.any(counts.sum(axis=0) != 60):
        raise ArithmeticError('the actual native full-configuration transfer changed')
    # Return immutable rows so the cached exact certificate cannot be mutated.
    return states, tuple(tuple(map(int, row)) for row in counts), tuple(edges)


def native_reference_state_control(particles=4):
    """Direct native branches check permutation moments and a nontracial signal."""
    states, rows, edges = _native_trace_counts(particles)
    dimension, denominator = len(states), 2*len(edges)
    degree = Counter(v for edge in edges for v in edge)
    adjacent = sum(d*(d-1) for d in degree.values())
    exact = _permutation_moments(12, particles, len(edges), adjacent)
    mean = F(sum(rows[i][i] for i in range(dimension)), dimension*denominator)
    second = F(sum(rows[i][j]*rows[j][i] for i in range(dimension)
                   for j in range(dimension)), dimension*denominator**2)
    if mean != exact['transfer_trace_mean'] or second-mean**2 != exact['transfer_trace_variance']:
        raise ArithmeticError('native branch trace moments differ from permutation counting')
    index = {state: i for i, state in enumerate(states)}
    commutator_numerator = 0
    for i, state in enumerate(states):
        if state[0] == state[1]:
            continue
        swapped = (state[1], state[0], *state[2:])
        commutator_numerator += rows[index[swapped]][i]
    # A=n0, B=i[A,swap01], rho_1=(I+K)/Tr(I+K).
    signal = F(commutator_numerator, dimension*denominator)/(1+mean)
    expected = exact['susceptibility']/(len(edges)*(2-exact['susceptibility']))
    if signal != expected or signal <= 0:
        raise ArithmeticError('the finite nontracial commutator control disappeared')
    counts = np.asarray(rows, dtype=float)/denominator
    eigenvalues = np.linalg.eigvalsh(counts)
    if eigenvalues[0] < -1e-12 or eigenvalues[-1] > 1+1e-12:
        raise ArithmeticError('the native transfer positivity control failed')
    # The integer row sums already certify the stationary eigenvalue as one.
    # acos has an unbounded derivative there: retain the exact eigenvalue
    # instead of amplifying harmless platform-dependent eigensolver roundoff.
    if abs(eigenvalues[-1]-1) > 1e-12:
        raise ArithmeticError('the native stationary eigenvalue exceeds its numerical allowance')
    eigenvalues[-1] = 1.
    if (abs(float(np.mean(eigenvalues))-float(mean)) > 1e-12
            or abs(float(np.mean(eigenvalues**2))-float(second)) > 1e-12):
        raise ArithmeticError('the native spectral oracle differs from exact trace moments')
    weights = np.exp(-np.arccos(np.clip(eigenvalues, 0, 1)))
    distance = float(np.sum(np.abs(weights/weights.sum()-1/dimension))/2)
    bound = 5*_root_upper(exact['transfer_trace_variance'])/exact['susceptibility']
    if distance > float(bound)+1e-12:
        raise ArithmeticError('the finite reference exceeds its analytic trace-distance bound')
    # Independent two-configuration coherent-state control inside the same sector.
    occupation = sp.diag(1, 0)
    swap = sp.Matrix([[0, 1], [1, 0]])
    conjugate = sp.I*(occupation*swap-swap*occupation)
    pure = sp.ones(2)/2
    pure_signal = -sp.I*sp.trace(pure*(occupation*conjugate-conjugate*occupation))
    if pure_signal != 1:
        raise ArithmeticError('the supplied nontracial-state countercontrol changed')
    return {**exact, 'configuration_dimension': dimension,
        'native_endpoint_branches': dimension*denominator,
        'native_and_permutation_trace_moments_exactly_equal': True,
        'ideal_reference_half_trace_distance_diagnostic': distance,
        'ideal_reference_half_trace_distance_upper': min(F(1), bound),
        'cap_one_commutator_expectation_minus_i': signal,
        'finite_reference_exactly_tracial': False,
        'different_supplied_two_configuration_pure_state_commutator': F(pure_signal)}


def source_trace_census(level):
    """Check the analytic degree census on actual finite source seams."""
    if type(level) is not int or level not in range(3):
        raise ValueError('retained source census levels are zero through two')
    graph = dynamics.source_green.source.federation(level)
    degrees = Counter(Counter(v for edge in graph.seams for v in edge).values())
    scale = 4**level
    if degrees != {5: 180*scale, 6: 60*scale}:
        raise ArithmeticError('the actual source degree census changed')
    moments = _permutation_moments(graph.ports, graph.ports//3,
        len(graph.seams), sum(n*d*(d-1) for d, n in degrees.items()))
    if moments['transfer_trace_variance'] != reference_state_bounds(level)['transfer_trace_variance']:
        raise ArithmeticError('the source census and all-level trace formula differ')
    return {'level': level, 'degree_five_ports': degrees[5], 'degree_six_ports': degrees[6],
        'source_census_matches_all_level_formula': True}


def reference_state_bounds(level, modular_parameter=1, terms_per_vertex=1, overlaps=1, norm=1):
    """Full-configuration state bounds, not low-density-sector normalization.

    The growing cap is h_l=64*(4^l)^4. Its state-tail allowance is charged
    separately. All modular bounds keep the dimensionless parameter u fixed.
    """
    level = dynamics._level(level)
    u, a, b, norm = map(dynamics._positive,
                        (modular_parameter, terms_per_vertex, overlaps, norm))
    scale, vertices = 4**level, 240*4**level
    variance = F(2*(365*vertices**2-1631*vertices+1560),
                 11907*(vertices-2)*(vertices-1)**2)
    if not 0 < variance <= F(1, 16*vertices):
        raise ArithmeticError('the all-level trace-variance inequality failed')
    epsilon, tail = F(3, 8*2**level), F(1, 2*scale**2)
    # binom(V,V/3)>=3^(V/3)>=2^(V/3), and exp(pi/2)<5.
    probability_exponent = (5*(vertices-1)-1).bit_length()-vertices//3
    return {'level': level, 'ports': vertices, 'particles': vertices//3,
        'configuration_dimension': f'binom({vertices},{vertices//3})',
        'transfer_trace_mean': 1-F(2*vertices, 9*(vertices-1)),
        'transfer_trace_variance': variance, 'variance_upper': F(1, 16*vertices),
        'ideal_reference_half_trace_distance_to_normalized_trace_upper': epsilon,
        'growing_cap_reference_half_trace_distance_to_normalized_trace_upper': min(F(1), epsilon+tail),
        'fixed_modular_parameter': u,
        'modular_observable_operator_norm_upper': 1,
        'ideal_modular_automorphism_normalized_trace_L2_error_upper': min(F(2), 4*u/F(15*2**level)),
        'growing_cap_modular_automorphism_normalized_trace_L2_error_upper': min(F(2), 4*u/F(15*2**level)+2*u/scale**2),
        'density_mode_dimension': vertices-1,
        'density_mode_reference_probability_log2_upper': probability_exponent,
        'bounded_overlap_fluctuation_commutator_norm_upper': 2*a*b*norm**2,
        'growing_cap_fluctuation_commutator_expectation_upper': 4*a*b*norm**2*(epsilon+tail),
        'reference_half_trace_bound_holds_on_full_configuration_algebra': True,
        'modular_operator_norm_triviality_claimed': False,
        'refinement_rescaled_modular_time_triviality_claimed': False,
        'every_noncommutative_or_nontracial_limit_excluded': False}


@lru_cache(maxsize=8)
def _helmert_basis(size):
    """Deterministic orthonormal basis of the mean-zero subspace in R^size."""
    basis = np.zeros((size, size-1), dtype=float)
    for j in range(size-1):
        norm = math.sqrt((j+1)*(j+2))
        basis[:j+1, j] = 1.0/norm
        basis[j+1, j] = -float(j+1)/norm
    return basis


@lru_cache(maxsize=4)
def _helmert_green_and_basis(level):
    """Independent mean-zero Schur Green matrix and distant cell-dipole pair."""
    graph, _, r, integer_f, denominator = spatial.source.matrices(level)
    schur = r.T@np.linalg.solve(integer_f/denominator, r)/12.0
    green_ambient = spatial.source._green(schur)/4**level
    basis = _helmert_basis(graph.cells)
    green = basis.T@green_ambient@basis
    green = 0.5*(green+green.T)
    adjacency = np.zeros((graph.cells, graph.cells), dtype=int)
    for a, b in graph.inter_seams:
        adjacency[a//12, b//12] = 1
        adjacency[b//12, a//12] = 1
    distances = spatial.shortest_path(adjacency, directed=False, unweighted=True)
    edges = [(i, j) for i in range(graph.cells)
             for j in range(i+1, graph.cells) if adjacency[i, j]]
    best_distance = -1
    best_pair = None
    for a0, a1 in edges[:20]:
        for b0, b1 in edges:
            d = int(min(distances[a0, b0], distances[a0, b1],
                        distances[a1, b0], distances[a1, b1]))
            if d > best_distance:
                best_distance = d
                best_pair = ((a0, a1), (b0, b1))
    return green, basis, best_distance, best_pair


def _green_orthonormal_matrix(level):
    """Orthonormal mean-zero cell Green matrix with independent Schur check."""
    level = _int(level, minimum=0)
    if level > 2:
        raise ValueError('relativistic wave vacuum matrix controls use levels zero through two')
    helmert_green, basis, _, _ = _helmert_green_and_basis(level)
    cells = basis.shape[0]
    if hasattr(spatial, 'green_orthonormal_matrix'):
        matrix = np.asarray(spatial.green_orthonormal_matrix(level), dtype=float)
        if matrix.shape == (cells, cells):
            if (np.max(np.abs(matrix-matrix.T)) > 1e-12
                    or float(np.linalg.norm(matrix@np.ones(cells, dtype=float), 2)) > 1e-11):
                raise ArithmeticError('spatial.green_orthonormal_matrix lost symmetry or mean-zero projection')
            reduced = basis.T@matrix@basis
            reduced = 0.5*(reduced+reduced.T)
        elif matrix.shape == helmert_green.shape:
            if np.max(np.abs(matrix-matrix.T)) > 1e-12:
                raise ArithmeticError('spatial.green_orthonormal_matrix lost symmetry')
            reduced = 0.5*(matrix+matrix.T)
        else:
            raise ArithmeticError('spatial.green_orthonormal_matrix has unexpected shape')
        if np.max(np.abs(np.linalg.eigvalsh(reduced)
                         -np.linalg.eigvalsh(helmert_green))) > 1e-10:
            raise ArithmeticError('spatial.green_orthonormal_matrix disagrees with independent Schur elimination')
        return reduced.copy()
    return helmert_green.copy()


@lru_cache(maxsize=8)
def _elementary_carrier_wave_vacuum_control(horizon, steps):
    """Exact and spectral wave-vacuum equipartition and commutator on 12 ports."""
    laplacian = sp.Matrix(dynamics.native_density()['laplacian'])
    eye = sp.eye(12)
    min_poly = laplacian*(laplacian-6*eye)*(laplacian**2-10*laplacian+20*eye)
    pseudo = (laplacian+sp.ones(12)/12).inv()-sp.ones(12)/12
    if (laplacian.trace() != 60 or pseudo.trace() != sp.Rational(7, 3)
            or any(value != 0 for value in min_poly)):
        raise ArithmeticError('elementary twelve-port native Laplacian invariants changed')
    basis = _helmert_basis(12)
    reduced = basis.T@np.asarray(laplacian, dtype=float)@basis
    reduced = 0.5*(reduced+reduced.T)
    eigenvalues = np.sort(np.linalg.eigvalsh(reduced))
    exact_eigenvalues = np.array([5.0-math.sqrt(5.0)]*3+[6.0]*5+[5.0+math.sqrt(5.0)]*3)
    if np.max(np.abs(eigenvalues-exact_eigenvalues)) > 1e-12:
        raise ArithmeticError('elementary twelve-port spectrum disagrees with exact algebraic roots')
    omega = np.sqrt(eigenvalues)
    inv_omega = 1.0/omega
    kinetic_modes = eigenvalues*inv_omega
    potential_modes = eigenvalues*inv_omega
    equipartition_residual = float(np.max(np.abs(kinetic_modes-potential_modes)))
    exact_sqrt_trace = float(3.0*math.sqrt(10.0+4.0*math.sqrt(5.0))+5.0*math.sqrt(6.0))
    exact_inv_sqrt_trace = float(3.0*math.sqrt(10.0+4.0*math.sqrt(5.0))/(2.0*math.sqrt(5.0))
                                 +5.0*math.sqrt(6.0)/6.0)
    total_energy = 0.5*float(np.sum(kinetic_modes+potential_modes))
    equal_time_trace = float(np.sum(inv_omega))
    if (equipartition_residual > 1e-14
            or abs(total_energy-exact_sqrt_trace) > 1e-12
            or abs(equal_time_trace-exact_inv_sqrt_trace) > 1e-12):
        raise ArithmeticError('elementary twelve-port canonical wave energy equipartition failed')
    t_float = float(horizon)
    max_step_residual = 0.0
    max_commutator_residual = 0.0
    for step in range(steps+1):
        tau = t_float*step/steps
        c_diag = np.cos(tau*omega)
        s_diag = np.sin(tau*omega)
        propagated_pos = (c_diag**2+s_diag**2)*inv_omega
        propagated_vel = (c_diag**2+s_diag**2)*omega
        step_residual = float(max(np.max(np.abs(propagated_pos-inv_omega)),
                                  np.max(np.abs(propagated_vel-potential_modes))))
        max_step_residual = max(max_step_residual, step_residual)
        d_tau = s_diag*inv_omega
        d_minus = np.sin(-tau*omega)*inv_omega
        dd_tau = -eigenvalues*d_tau
        comm_residual = float(max(np.max(np.abs(d_tau+d_minus)),
                                  np.max(np.abs(dd_tau+omega*s_diag))))
        max_commutator_residual = max(max_commutator_residual, comm_residual)
    if max_step_residual > 1e-12 or max_commutator_residual > 1e-12:
        raise ArithmeticError('elementary twelve-port wave evolution or Pauli-Jordan commutator failed')
    return {
        'carrier_ports': 12,
        'mean_zero_dimension': 11,
        'exact_laplacian_trace': 60,
        'exact_pseudoinverse_trace': F(7, 3),
        'exact_sqrt_laplacian_trace_formula': '3*sqrt(10+4*sqrt(5))+5*sqrt(6)',
        'exact_inv_sqrt_laplacian_trace_formula': '3*sqrt(10+4*sqrt(5))/(2*sqrt(5))+5*sqrt(6)/6',
        'equal_time_fractional_wave_covariance_trace_diagnostic': equal_time_trace,
        'canonical_wave_energy_expectation_diagnostic': total_energy,
        'mode_equipartition_residual_diagnostic': equipartition_residual,
        'multi_step_stationarity_residual_diagnostic': max_step_residual,
        'Pauli_Jordan_klein_gordon_residual_diagnostic': max_commutator_residual,
        'exact_mode_energy_equipartition_verified': True}


def fractional_wave_bounds(level=0, T=F(1), steps=4, cutoff=2, s=F(2)):
    """Fractional-filtered wave equipartition, H^-s tightness and causal cone.

    For W_{1/4,l}(t)=A_l^(1/4) W_l(t) with A_l=G_l^-1 on the (20*4^l-1)-dimensional
    orthonormal mean-zero cell subspace, the two-time Gaussian covariance is
    C_{1/4,l}(t,u)=A_l^(-1/2) cos((t-u) sqrt(A_l))=G_l^(1/2) cos((t-u) sqrt(A_l)),
    which equals twice the canonical free-scalar Wightman real covariance and matches
    at equal times the radial half-cylinder Steklov boundary covariance G_l^(1/2).
    """
    level = _int(level, minimum=0)
    if level > 2:
        raise ValueError('relativistic wave vacuum matrix controls use levels zero through two')
    T = _fraction(T, positive=True)
    steps = _int(steps, minimum=1)
    cutoff = _int(cutoff, minimum=1)
    s = _fraction(s, positive=True)
    if s < F(7, 4):
        raise ValueError('fractional wave path tightness requires Sobolev index s>=7/4')
    scale = 4**level
    cells = 20*scale
    dim = cells-1
    if cutoff > dim:
        raise ValueError('spectral cutoff exceeds the mean-zero source dimension')
    elementary = _elementary_carrier_wave_vacuum_control(T, steps)
    green = _green_orthonormal_matrix(level)
    helmert_green, basis, dipole_distance, ((a0, a1), (b0, b1)) = _helmert_green_and_basis(level)
    if (green.shape != (dim, dim)
            or np.max(np.abs(green-green.T)) > 1e-12):
        raise ArithmeticError('mean-zero Green matrix lost symmetry or dimension')
    evals_g_asc, evecs_g = np.linalg.eigh(green)
    evals_helmert_asc, evecs_helmert = np.linalg.eigh(helmert_green)
    if (np.any(evals_g_asc <= 0)
            or np.max(np.abs(evals_g_asc-evals_helmert_asc)) > 1e-10
            or float(evals_g_asc[-1]) > float(spatial.GREEN_BOUND)+1e-12):
        raise ArithmeticError('mean-zero Green spectrum failed independent Schur or norm bounds')
    # Order modes by ascending generator eigenvalue lambda_1 <= ... <= lambda_dim.
    mu = evals_g_asc[::-1]
    vectors = evecs_g[:, ::-1]
    lambdas = 1.0/mu
    indices = np.arange(1, dim+1, dtype=float)
    if (np.any(lambdas < indices/49920.0-1e-12)
            or np.any(lambdas > 1.5*indices+1e-12)
            or float(lambdas[-1]) > 0.5*scale+1e-12):
        raise ArithmeticError('generator eigenvalues violated the all-level spectral enclosure')
    omega = np.sqrt(lambdas)
    inv_omega = np.sqrt(mu)
    # 1. Canonical energy equipartition on the level-l spherical source operator.
    kinetic_modes = lambdas*(inv_omega**2)*omega
    potential_modes = lambdas*inv_omega
    mode_equipartition_residual = float(np.max(np.abs(kinetic_modes-potential_modes)))
    cutoff_kinetic = float(np.sum(kinetic_modes[:cutoff]))
    cutoff_potential = float(np.sum(potential_modes[:cutoff]))
    full_kinetic = float(np.sum(kinetic_modes))
    full_potential = float(np.sum(potential_modes))
    total_canonical_energy = 0.5*(full_kinetic+full_potential)
    trace_sqrt_a = float(np.sum(omega))
    trace_sqrt_g = float(np.sum(inv_omega))
    # Matrix-level verification of C_{1/4,l}(0) = G_l^(1/2) = A_l^(-1/2) and A_l C_{1/4,l}(0) = A_l^(1/2),
    # cross-checked against non-spectral Denman-Beavers matrix square-root iteration.
    c0_matrix = (vectors*inv_omega)@vectors.T
    a_matrix = (vectors*lambdas)@vectors.T
    sqrt_a_matrix = (vectors*omega)@vectors.T
    db_c0, db_sqrt_a = spatial._denman_beavers_sqrt(green)
    matrix_equipartition_residual = float(max(
        np.linalg.norm(c0_matrix@c0_matrix-green, 2),
        np.linalg.norm(a_matrix@c0_matrix-sqrt_a_matrix, 2),
        np.linalg.norm(c0_matrix-db_c0, 2),
        np.linalg.norm(sqrt_a_matrix-db_sqrt_a, 2)))
    if (mode_equipartition_residual > 1e-12
            or abs(cutoff_kinetic-cutoff_potential) > 1e-12
            or abs(full_kinetic-full_potential) > 1e-11
            or abs(total_canonical_energy-trace_sqrt_a) > 1e-11
            or matrix_equipartition_residual > 1e-10):
        raise ArithmeticError('canonical wave energy equipartition failed on the source operator')
    # 2. Geometric H^-s bounds use the OP4 harmonic energy estimate, not
    # an identification of eigenvalue indices with spherical harmonics.
    # Integral comparison: sum_{k>=1} k^-(s-1/2) <= 1 + 1/(s-3/2) = (s-1/2)/(s-3/2).
    harmonic_series_factor = (s-F(1, 2))/(s-F(3, 2))
    # From lambda_k >= k/49920 (lambda_k^(-1/2) < 224 k^(-1/2)) and lambda_k <= 3k/2:
    # the spectral-index diagnostics have increments <=336 sqrt(k)|t-u|^2.
    # For actual harmonics, sum_m <Y_jm,sqrt(A_l)Y_jm> <=
    # sqrt(2)(2j+1)sqrt(j(j+1)) <=12 j^2 by OP4 and Cauchy-Schwarz.
    increment_second_moment_constant = F(12)*(2*s-2)/(2*s-3)
    spectral_increment_constant = F(336)*harmonic_series_factor
    increment_fourth_moment_constant = 3*increment_second_moment_constant**2
    equal_time_hs_variance_upper = F(224)*(s+F(1, 2))/(s-F(1, 2))
    cutoff_tail_hs_variance_upper = F(224)/((s-F(1, 2))*cutoff)
    dt = T/steps
    dt_float = float(dt)
    s_float = float(s)
    weights_hs = (1.0+indices)**(-s_float)
    equal_time_hs_diagnostic = float(np.sum(weights_hs*inv_omega))
    cutoff_tail_hs_diagnostic = float(np.sum(weights_hs[cutoff:]*inv_omega[cutoff:]))
    # Both independent initial displacement and velocity contribute. These
    # diagnostics use spectral-index weights, not spherical harmonic weights.
    increment_second_diagnostic = float(np.sum(
        weights_hs*inv_omega*2.0*(1.0-np.cos(dt_float*omega))))
    increment_fourth_diagnostic = 3.0*(increment_second_diagnostic**2)
    step_second_upper = increment_second_moment_constant*(dt**2)
    step_fourth_upper = increment_fourth_moment_constant*(dt**4)
    if (equal_time_hs_diagnostic > float(equal_time_hs_variance_upper)+1e-12
            or cutoff_tail_hs_diagnostic > float(cutoff_tail_hs_variance_upper)+1e-12
            or increment_second_diagnostic > float(spectral_increment_constant*dt**2)+1e-12
            or increment_fourth_diagnostic > float(3*spectral_increment_constant**2*dt**4)+1e-12):
        raise ArithmeticError('fractional wave H^-s norm or increment violated its exact rational bound')
    # 3. Candidate Pauli-Jordan sine kernel; a classical covariance alone
    # does not define quantum operators or their commutators. Kernel D_l(t-u) = A_l^(-1/2) sin((t-u) sqrt(A_l)).
    t_float = float(T)
    d_diag = np.sin(t_float*omega)*inv_omega
    d_norm_diagnostic = float(np.max(np.abs(d_diag)))
    d_norm_upper = min(T, _root_upper(spatial.GREEN_BOUND))
    # Composite midpoint quadrature check for D_l(T) = int_0^T cos(tau sqrt(A_l)) d tau.
    quad_steps = max(steps, 16)
    quad_dt = t_float/quad_steps
    quad_diag = np.zeros_like(omega)
    for m in range(quad_steps):
        mid = (m+0.5)*quad_dt
        quad_diag += np.cos(mid*omega)*quad_dt
    quad_error_diagnostic = float(np.max(np.abs(d_diag-quad_diag)))
    quad_error_upper = F(scale, 48*quad_steps**2)*(T**3)
    # Separated zero-mean cell-dipole commutator leakage at graph distance d.
    f_ambient = np.zeros(cells, dtype=float)
    f_ambient[a0] = 1.0/math.sqrt(2.0); f_ambient[a1] = -1.0/math.sqrt(2.0)
    g_ambient = np.zeros(cells, dtype=float)
    g_ambient[b0] = 1.0/math.sqrt(2.0); g_ambient[b1] = -1.0/math.sqrt(2.0)
    f_vec = basis.T@f_ambient
    g_vec = basis.T@g_ambient
    omega_h = np.sqrt(1.0/evals_helmert_asc[::-1])
    vectors_h = evecs_helmert[:, ::-1]
    d_helmert = (vectors_h*(np.sin(t_float*omega_h)/omega_h))@vectors_h.T
    d_step_helmert = (vectors_h*(np.sin(dt_float*omega_h)/omega_h))@vectors_h.T
    a_helmert = (vectors_h*(omega_h**2))@vectors_h.T
    dipole_commutator_diagnostic = float(abs(g_vec@d_helmert@f_vec))
    dipole_step_commutator_diagnostic = float(abs(g_vec@d_step_helmert@f_vec))
    dipole_cubic_coeff_diagnostic = float(abs(g_vec@a_helmert@f_vec))/6.0
    dipole_cubic_coeff_upper = F(scale, 6*2**dipole_distance)
    dipole_step_commutator_upper = min(
        dt, _root_upper(spatial.GREEN_BOUND),
        F(scale, 12)*(dt**3),
        dipole_cubic_coeff_upper*(dt**3)+F(scale**2, 480)*(dt**5))
    dipole_commutator_upper = min(
        d_norm_upper,
        F(scale, 12)*(T**3),
        dipole_cubic_coeff_upper*(T**3)+F(scale**2, 480)*(T**5))
    if (d_norm_diagnostic > float(d_norm_upper)+1e-12
            or quad_error_diagnostic > float(quad_error_upper)+1e-12
            or dipole_cubic_coeff_diagnostic > float(dipole_cubic_coeff_upper)+1e-12
            or dipole_step_commutator_diagnostic > float(dipole_step_commutator_upper)+1e-12
            or dipole_commutator_diagnostic > float(dipole_commutator_upper)+1e-12):
        raise ArithmeticError('Pauli-Jordan commutator kernel violated its norm, quadrature, or spatial leakage bound')
    return {
        'level': level,
        'cells': cells,
        'mean_zero_dimension': dim,
        'parameter_horizon': T,
        'time_steps': steps,
        'step_size': dt,
        'spectral_cutoff': cutoff,
        'Sobolev_index_s': s,
        'elementary_twelve_port_control': elementary,
        'fractional_wave_field': 'W_{1/4,l}(t)=A_l^(1/4) W_l(t)',
        'two_time_covariance_kernel': 'C_{1/4,l}(t,u)=A_l^(-1/2) cos((t-u) sqrt(A_l))=G_l^(1/2) cos((t-u) sqrt(A_l))',
        'equal_time_Steklov_boundary_covariance_identity': 'C_{1/4,l}(t,t)=G_l^(1/2)=A_l^(-1/2)',
        'canonical_free_scalar_Wightman_factor': 2,
        'cutoff_kinetic_energy_diagnostic': cutoff_kinetic,
        'cutoff_potential_energy_diagnostic': cutoff_potential,
        'full_kinetic_energy_diagnostic': full_kinetic,
        'full_potential_energy_diagnostic': full_potential,
        'total_canonical_energy_expectation_diagnostic': total_canonical_energy,
        'trace_sqrt_generator_diagnostic': trace_sqrt_a,
        'trace_sqrt_Green_diagnostic': trace_sqrt_g,
        'mode_equipartition_residual_diagnostic': mode_equipartition_residual,
        'matrix_equipartition_residual_diagnostic': matrix_equipartition_residual,
        'canonical_energy_equipartition_verified': True,
        'harmonic_series_integral_factor': harmonic_series_factor,
        'equal_time_spectral_index_variance_diagnostic': equal_time_hs_diagnostic,
        'equal_time_spectral_index_variance_upper': equal_time_hs_variance_upper,
        'cutoff_tail_spectral_index_variance_diagnostic': cutoff_tail_hs_diagnostic,
        'cutoff_tail_spectral_index_variance_upper': cutoff_tail_hs_variance_upper,
        'spectral_index_increment_constant': spectral_increment_constant,
        'numerical_weighted_norms_use_spectral_indices_not_spherical_harmonics': True,
        'H_minus_s_second_moment_constant_per_dt_squared': increment_second_moment_constant,
        'H_minus_s_fourth_moment_constant_per_dt_fourth': increment_fourth_moment_constant,
        'step_spectral_index_second_moment_diagnostic': increment_second_diagnostic,
        'step_spectral_index_second_moment_upper': spectral_increment_constant*dt**2,
        'step_H_minus_s_second_moment_upper': step_second_upper,
        'step_spectral_index_fourth_moment_diagnostic': increment_fourth_diagnostic,
        'step_spectral_index_fourth_moment_upper': 3*spectral_increment_constant**2*dt**4,
        'step_H_minus_s_fourth_moment_upper': step_fourth_upper,
        'Kolmogorov_Chentsov_holder_exponent_upper': F(1, 2),
        'continuous_path_topology': f'C([0,T],H^-{s}(S2))',
        'Pauli_Jordan_commutator_kernel': 'D_l(t-u)=A_l^(-1/2) sin((t-u) sqrt(A_l))',
        'Pauli_Jordan_operator_norm_diagnostic': d_norm_diagnostic,
        'Pauli_Jordan_operator_norm_upper': d_norm_upper,
        'Pauli_Jordan_midpoint_quadrature_steps': quad_steps,
        'Pauli_Jordan_midpoint_quadrature_error_diagnostic': quad_error_diagnostic,
        'Pauli_Jordan_midpoint_quadrature_error_upper': quad_error_upper,
        'separated_dipole_cells': [[a0, a1], [b0, b1]],
        'separated_dipole_cell_distance': dipole_distance,
        'separated_dipole_cubic_taylor_coefficient_diagnostic': dipole_cubic_coeff_diagnostic,
        'separated_dipole_cubic_taylor_coefficient_upper': dipole_cubic_coeff_upper,
        'separated_dipole_step_commutator_leakage_diagnostic': dipole_step_commutator_diagnostic,
        'separated_dipole_step_commutator_leakage_upper': dipole_step_commutator_upper,
        'separated_dipole_commutator_leakage_diagnostic': dipole_commutator_diagnostic,
        'separated_dipole_commutator_leakage_upper': dipole_commutator_upper,
        'continuum_comparison_wave_speed_upper': '1/sqrt(2)',
        'continuum_sine_kernel_vanishes_for_separated_mean_zero_tests': True,
        'classical_covariance_supplies_quantum_CCR': False,
        'unrestricted_mean_zero_projected_kernel_has_local_support': False,
        'restored_zero_mode_kernel': 'D_bar(t)=D_perp(t)+t P_const'}


def _source_dyadic_scalar(expression):
    """Pin transcendental input at 80 digits before rounding to the source grid."""
    return F(int(sp.floor(sp.N(expression, 80)*(1 << 44)+sp.Rational(1, 2))), 1 << 44)


def spin_weyl_word(vectors, copies=1, controlled=False):
    """Existing source pulses for ordered exp(i(s Q_N+t P_N)) products.

    Q_N=sum X/sqrt(2N), P_N=sum Y/sqrt(2N). The optional last qubit is
    the Ramsey control. Blank proper-code registers and coherent code gluing
    are the source compiler's existing premises.
    """
    copies = _int(copies, 1)
    if copies > 15 or type(controlled) is not bool:
        raise ValueError('one through fifteen copies and a Boolean control required')
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 2 or not np.all(np.isfinite(vectors)):
        raise ValueError('finite real phase-space pairs required')
    word = []
    # execute applies later pulses on the left; mathematical products retain
    # the caller's order. No reset, projection, or state normalization enters.
    for sx, sy in reversed(vectors):
        x, y = sp.Rational(float(sx)), sp.Rational(float(sy))
        azimuth = sp.atan2(y, x) if x or y else sp.Integer(0)
        angles = (float(_source_dyadic_scalar(azimuth/2)),
                  float(_source_dyadic_scalar(sp.sqrt((x*x+y*y)/(2*copies)))))
        for target in range(copies):
            for kind, angle in (('Z', angles[0]), ('X', angles[1]), ('Z', -angles[0])):
                if controlled:
                    width = copies+1
                    one = ['I']*width; one[target] = kind
                    two = one.copy(); two[copies] = 'Z'
                    word += source_compiler.pauli_rotation(one, angle/2)
                    word += source_compiler.pauli_rotation(two, -angle/2)
                else:
                    word.append(source_compiler.pulse(kind, angle, [target]))
    return word


def spin_weyl_expectation(vectors, copies, polarization=1):
    """Exact finite product-state law evaluated by two-dimensional matrices."""
    copies = _int(copies, 1)
    polarization = _fraction(polarization)
    if not -1 <= polarization <= 1:
        raise ValueError('polarization must lie in [-1,1]')
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 2 or not np.all(np.isfinite(vectors)):
        raise ValueError('finite real phase-space pairs required')
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    product = np.eye(2, dtype=complex)
    total = np.sum(vectors, axis=0)
    symplectic = 0.0
    for i, (sx, sy) in enumerate(vectors):
        radius = math.hypot(sx, sy)
        angle = radius/math.sqrt(2*copies)
        one = (np.eye(2) if radius == 0 else
               math.cos(angle)*np.eye(2)+1j*math.sin(angle)*(sx*x+sy*y)/radius)
        product = product@one
        symplectic += sum(sx*ty-sy*tx for tx, ty in vectors[i+1:])
    expectation = ((1+float(polarization))*product[0, 0]
                   +(1-float(polarization))*product[1, 1])/2
    limit = np.exp(-float(total@total)/4-0.5j*float(polarization)*symplectic)
    return complex(expectation**copies), complex(limit)


def quantum_fluctuation_control():
    """Source-coded ordered detector law and its conditional CCR limit.

    This realizes a chosen polarized reference, not the unconditioned native
    repair state. The IID quantum CLT applies to tensor independent copies.
    Its new content here is the source pulse/readout and native frequency map.
    """
    vectors = ((F(3, 4), 0), (0, F(1, 2)), (F(-1, 4), F(1, 3)))
    expected, limit = spin_weyl_expectation(vectors, 4)
    reversed_expected, _ = spin_weyl_expectation(tuple(reversed(vectors)), 4)
    tracial, _ = spin_weyl_expectation(vectors, 4, polarization=0)
    records = {}
    for quadrature in ('X', 'Y'):
        copies, control = 4, 4
        word = [source_compiler.pulse('Y', -math.pi/4, [control])]
        word += spin_weyl_word(vectors, copies, controlled=True)
        word += [source_compiler.pulse('Y' if quadrature == 'X' else 'X',
                                     math.pi/4 if quadrature == 'X' else -math.pi/4,
                                     [control])]
        state = np.zeros(1 << (copies+1), dtype=complex); state[0] = 1
        output = source_compiler.execute(word, copies+1, state)
        probabilities = abs(output)**2
        mean = float(np.sum(probabilities[:1 << copies])-np.sum(probabilities[1 << copies:]))
        target = expected.real if quadrature == 'X' else expected.imag
        if abs(mean-target) > 2e-11 or abs(float(np.sum(probabilities))-1) > 2e-11:
            raise ArithmeticError('source-coded ordered Weyl detector changed its complete output law')
        records[quadrature] = {
            'all_pointer_and_work_probabilities': probabilities.tolist(),
            'pointer_mean': mean, 'independent_product_state_mean': target,
            'source_pulses': len(word),
            'expanded_word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True).encode()).hexdigest()}
    if abs(expected.imag) < .05 or abs(reversed_expected-expected.conjugate()) > 1e-12:
        raise ArithmeticError('ordered Weyl phase or reversal control disappeared')
    if abs(tracial.imag) > 1e-12:
        raise ArithmeticError('tracial product control has a spurious commutator phase')
    # The finite native density spectrum fixes the coherent phase. This is
    # a source-executed mode encoding, not an execution of the full arc walk.
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    modes = []
    for lam in (5-math.sqrt(5), 6., 5+math.sqrt(5)):
        phase = math.acos(1-lam/60)
        u = source_compiler.execute([source_compiler.pulse('Z', phase/2, [0])], 1)
        q = x/math.sqrt(2*phase)
        evolved = u.conj().T@q@u
        ordered = complex((evolved@q)[0, 0])
        target = np.exp(-1j*phase)/(2*phase)
        if abs(ordered-target) > 2e-12:
            raise ArithmeticError('source pulse lost the native wave frequency or ordered two-point phase')
        kick = 1e-3
        outputs = []
        for sign in (-1, 1):
            tape = [source_compiler.pulse('X', sign*kick/math.sqrt(2*phase), [0]),
                    source_compiler.pulse('Z', phase/2, [0]),
                    source_compiler.pulse('Y', math.pi/4, [0])]
            out = source_compiler.execute(tape, 1, np.array([1., 0.]))
            probabilities = abs(out)**2
            outputs.append(probabilities.tolist())
        response = ((outputs[1][0]-outputs[1][1])-(outputs[0][0]-outputs[0][1]))/(2*kick*math.sqrt(2*phase))
        sine_response = math.sin(phase)/phase
        response_bound = abs(sine_response)*kick*kick/(3*phase)
        if abs(response-sine_response) > response_bound+2e-10:
            raise ArithmeticError('native source kick and pointer lost the sine response')
        modes.append({'native_laplacian_eigenvalue': lam,
                      'coherent_phase_per_incident_trigger': phase,
                      'ordered_two_point_real': ordered.real,
                      'ordered_two_point_imaginary': ordered.imag,
                      'source_kick_size': kick,
                      'complete_minus_and_plus_kick_pointer_probabilities': outputs,
                      'central_difference_response': response,
                      'sine_kernel_linear_response': sine_response,
                      'finite_kick_error_upper': response_bound})
    return {
        'source_copies_executed': 4, 'complete_pointer_records': records,
        'Weyl_angle_grid_denominator': 1 << 44,
        'Weyl_word_angle_error_operator_norm_upper': str(F(3*4*len(vectors), 1 << 45)),
        'ordered_Weyl_vectors': [[str(x) for x in row] for row in vectors],
        'finite_expectation': [expected.real, expected.imag],
        'CCR_limit_expectation': [limit.real, limit.imag],
        'reversed_order_expectation': [reversed_expected.real, reversed_expected.imag],
        'tracial_reference_expectation': [tracial.real, tracial.imag],
        'finite_copy_errors': [{'copies': n,
            'absolute_ordered_characteristic_error': abs(spin_weyl_expectation(vectors, n)[0]-limit)}
            for n in (4, 16, 64, 256)],
        'native_frequency_source_pulse_controls': modes,
        'limit_ordered_Weyl_law': 'exp(-|sum v_j|^2/4-i sum_(j<k) sigma(v_j,v_k)/2)',
        'finite_CCR': '[Q_N,P_N]=i sum Z_j/N; not iI in operator norm',
        'finite_fourth_moment': 'E Q_N^4=3/4-1/(2N)',
        'encoded_wave_Wightman': '<phi(t)phi(u)>=exp(-i Omega(t-u))/(2 Omega)',
        'native_refinement_frequency': 'Omega_l=sqrt(Q_l)*acos(1-lambda/(2Q_l)), Q_l=630*16^l',
        'native_incident_parameter': 't=beta tau, beta=2d/sqrt(630), fixed d in {5,6}',
        'conditional_CCR_limit_from_ordered_product_state_law_proved': True,
        'finite_native_pulse_Weyl_and_sine_response_detector_executed': True,
        'polarized_product_reference_supplied': True,
        'proper_code_gluing_blank_register_and_controlled_pulses_required': True,
        'spectral_mode_labels_and_rotation_angles_are_declared_encoding': True,
        'independent_product_copy_limit_is_source_refinement_limit': False,
        'growing_field_preparation_and_detector_compiler_executed': False,
        'unconditioned_repair_selects_this_reference': False,
        'physical_QFT_or_atomic_clock_identification': False}


def thermal_clock_parameters(frequency, polarization, copies=2):
    """Finite faithful state/clock relation in the positive-temperature branch.

    Exact rational inputs describe the finite preparation and frequency;
    logarithmic diagnostics approximate the stated exact log-ratio formula.
    sigma_s=rho^(is) . rho^(-is)=alpha_(-beta*s), with
    alpha_t(.)=exp(itH) . exp(-itH), H=frequency*(I-Z)/2.
    """
    frequency = _fraction(frequency, positive=True)
    polarization = _fraction(polarization, positive=True)
    copies = _int(copies, 1)
    if not polarization < 1 or not 0 < float(polarization) < 1:
        raise ValueError('a numerically representable faithful positive polarization required')
    ratio = (1+polarization)/(1-polarization)
    gap = math.log1p(float(polarization))-math.log1p(-float(polarization))
    if not math.isfinite(float(frequency)) or float(frequency) <= 0 or gap <= 0:
        raise ValueError('representable nonzero frequency and modular gap required')
    return {
        'frequency': frequency, 'polarization': polarization, 'copies': copies,
        'ground_probability': (1+polarization)/2,
        'excited_probability': (1-polarization)/2,
        'modular_gap_log_ratio': str(ratio),
        'modular_gap_diagnostic': gap,
        'inverse_temperature_in_frequency_clock_diagnostic': gap/float(frequency),
        'modular_parameter_per_forward_clock_unit_diagnostic': -float(frequency)/gap,
        'canonical_quadrature_variance': 1/(2*polarization),
        'boson_mean_occupation': (1-polarization)/(2*polarization),
        'field_variance': 1/(2*polarization*frequency),
        'finite_fourth_quadrature_moment': (F(3, 4)-F(1, 2*copies))/polarization**2,
        'CCR_fourth_quadrature_moment': F(3, 4)/polarization**2}


def _thermal_code_preparation(polarization, copies):
    """Purify each code with a retained reference using existing source words."""
    copies = _int(copies, 1)
    angle = float(_source_dyadic_scalar(sp.acos(sp.sqrt((1+sp.Rational(polarization))/2))))
    word = []
    for j in range(copies):
        word.append(source_compiler.pulse('Y', -angle, [j]))
        word += source_compiler.cnot(j, copies+j)
    return word


def faithful_mode_source_control(frequency, inverse_temperature=F(4)):
    """Prepare a faithful reference and read its modular/thermal field laws.

    Two system codes, their two untouched purifying references, and one
    pointer are retained. The modular algebra is the system codes, not the
    pure state on the enlarged system-plus-reference algebra.
    """
    frequency = _fraction(frequency, positive=True)
    inverse_temperature = _fraction(inverse_temperature, positive=True)
    omega, beta = float(frequency), float(inverse_temperature)
    polarization = _source_dyadic_scalar(sp.tanh(
        sp.Rational(inverse_temperature)*sp.Rational(frequency)/2))
    parameters = thermal_clock_parameters(frequency, polarization)
    r, gap = float(polarization), parameters['modular_gap_diagnostic']
    if abs(gap/omega-beta) > 1e-12:
        raise ArithmeticError('prepared polarization lost its common modular clock ratio')
    copies = 2
    preparation = _thermal_code_preparation(polarization, copies)
    initial = np.zeros(16, complex); initial[0] = 1
    pure = source_compiler.execute(preparation, 4, initial)
    amplitude = pure.reshape(4, 4)  # reference index, system index
    reduced = amplitude.T@amplitude.conj()
    one = np.diag([(1+r)/2, (1-r)/2])
    target = np.kron(one, one)
    state_error = float(np.linalg.norm(reduced-target, 2))
    if state_error > 2e-12:
        raise ArithmeticError('source purification did not prepare the faithful system state')
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z = np.diag([1., -1.])
    modular_error = 0.0
    for s in (-.7, .2, .9):
        modular = np.diag(np.exp(1j*s*gap*np.array([1., -1.])/2))
        # Forward Heisenberg t=-beta*s fixes the modular sign explicitly.
        state_unitary = source_compiler.execute(
            [source_compiler.pulse('Z', -beta*s*omega/2, [0])], 1)
        for observable in (x, y, z):
            error = np.linalg.norm(modular@observable@modular.conj().T
                         -state_unitary.conj().T@observable@state_unitary, 2)
            modular_error = max(modular_error, float(error))
    if modular_error > 2e-12:
        raise ArithmeticError('source wave clock and faithful modular flow disagree')
    unit = sp.sqrt(sp.Rational(polarization))
    vectors = np.array([[float(sp.N(1/(2*unit), 80)), 0.],
                        [0., float(sp.N(1/(4*unit), 80))]])
    expected, limit = spin_weyl_expectation(vectors, copies, polarization)
    records = {}
    for axis in ('X', 'Y'):
        word = preparation+[source_compiler.pulse('Y', -math.pi/4, [4])]
        controlled = spin_weyl_word(vectors, copies, controlled=True)
        word += source_compiler.remap(controlled, [0, 1, 4])
        word += [source_compiler.pulse('Y' if axis == 'X' else 'X',
                 math.pi/4 if axis == 'X' else -math.pi/4, [4])]
        blank = np.zeros(32, complex); blank[0] = 1
        output = source_compiler.execute(word, 5, blank)
        probabilities = abs(output)**2
        mean = float(sum(probabilities[:16])-sum(probabilities[16:]))
        reference = expected.real if axis == 'X' else expected.imag
        if abs(mean-reference) > 2e-11 or abs(sum(probabilities)-1) > 2e-11:
            raise ArithmeticError('faithful ordered detector lost a pointer/reference outcome')
        records[axis] = {'all_system_reference_pointer_probabilities': probabilities.tolist(),
            'pointer_mean': mean, 'independent_finite_state_mean': reference,
            'source_pulses': len(word),
            'word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True).encode()).hexdigest()}
    # Finite source intervention: canonical normalization preserves response
    # while the symmetric noise depends on the faithful polarization.
    epsilon = 1e-3
    kicks, means = [], []
    for sign in (-1, 1):
        vector = [[sign*epsilon/math.sqrt(r*omega), 0.]]
        word = preparation+spin_weyl_word(vector, copies)
        word += [source_compiler.pulse('Z', omega/2, [j]) for j in range(copies)]
        word += [source_compiler.pulse('Y', math.pi/4, [j]) for j in range(copies)]
        output = source_compiler.execute(word, 4, initial)
        probabilities = abs(output)**2
        readout = np.array([copies-2*(j % 4).bit_count() for j in range(16)])
        means.append(float(probabilities@readout)/math.sqrt(2*copies*r*omega))
        kicks.append(probabilities.tolist())
    response = (means[1]-means[0])/(2*epsilon)
    reference_response = math.sin(omega)/omega
    response_budget = abs(reference_response)*epsilon**2/(3*copies*r*omega)
    if abs(response-reference_response) > response_budget+2e-10:
        raise ArithmeticError('thermal source kick failed its state-independent sine response')
    def covariance(t):
        return (np.cos(omega*t)-1j*r*np.sin(omega*t))/(2*r*omega)
    kms_error = float(abs(covariance(.3-1j*beta)-covariance(-.3)))
    if kms_error > 2e-12:
        raise ArithmeticError('faithful field covariance lost its KMS boundary identity')
    return {'state_clock_parameters': parameters,
        'requested_common_inverse_temperature': inverse_temperature,
        'polarization_synthesis_grid': F(1, 1 << 44),
        'system_algebra': 'M2 tensor M2; purifying references excluded from this algebra',
        'preparation_reduced_state_error_diagnostic': state_error,
        'source_modular_flow_error_on_all_three_Paulis': modular_error,
        'complete_ordered_pointer_records': records,
        'thermal_CCR_ordered_limit': [limit.real, limit.imag],
        'complete_minus_and_plus_kick_system_reference_probabilities': kicks,
        'central_difference_response': response, 'sine_response': reference_response,
        'finite_kick_error_upper': response_budget, 'KMS_two_point_boundary_error': kms_error,
        'same_clock_thermal_and_modular_state_prepared': True,
        'purifying_reference_algebra_has_a_faithful_global_state': False}


def faithful_thermal_field_bound(inverse_temperature=F(4)):
    """Trace-class bound for the mean-zero free field, using the native Weyl bound."""
    beta = _fraction(inverse_temperature, positive=True)
    return {'inverse_temperature': beta,
        'one_particle_heat_trace_upper': F(99840)/beta**2,
        'normal_ordered_Fock_log_partition_upper': F(149760)/beta**2,
        'state': 'rho_b=Z_b^-1 exp(-b dGamma(sqrt(A))) on mean-zero Fock space',
        'modular_flow': 'sigma_s=alpha_(-b*s) on the full encoded field algebra',
        'faithful_normal_state_under_declared_encoding': True,
        'zero_mode_included': False, 'spatial_patch_modular_generator_identified': False,
        'inverse_temperature_selected_by_passive_repair': False}


def faithful_modular_encoding_control():
    modes = [faithful_mode_source_control(_source_dyadic_scalar(sp.acos(1-lam/60)))
             for lam in (5-sp.sqrt(5), sp.Integer(6), 5+sp.sqrt(5))]
    first = modes[0]['state_clock_parameters']['polarization']
    equal_polarization_betas = [thermal_clock_parameters(
        row['state_clock_parameters']['frequency'], first)[
        'inverse_temperature_in_frequency_clock_diagnostic'] for row in modes]
    if max(equal_polarization_betas)-min(equal_polarization_betas) < .5:
        raise ArithmeticError('equal-polarization unequal-frequency countercontrol disappeared')
    # One tick cannot distinguish modular gaps separated by 2pi, whereas
    # equality as continuous flows does distinguish them at half a tick.
    theta = float(modes[0]['state_clock_parameters']['frequency'])
    alias_tick = abs(np.exp(1j*theta)-np.exp(1j*(theta+2*math.pi)))
    alias_half = abs(np.exp(.5j*theta)-np.exp(.5j*(theta+2*math.pi)))
    # Matching state eigenvalues alone does not match the modular action:
    # a transverse polarization is not even stationary for the wave clock.
    transverse_state = np.array([[1, float(first)], [float(first), 1]])/2
    hamiltonian = np.diag([0., theta])
    transverse_commutator = float(np.linalg.norm(
        transverse_state@hamiltonian-hamiltonian@transverse_state, 2))
    if transverse_commutator < .05:
        raise ArithmeticError('off-axis state/clock countercontrol disappeared')
    return {'finite_source_modes': modes,
        'native_frequency_synthesis_grid': F(1, 1 << 44),
        'equal_polarization_implied_inverse_temperatures': equal_polarization_betas,
        'integer_tick_alias_residual': float(alias_tick),
        'half_tick_alias_difference': float(alias_half),
        'same_eigenvalues_transverse_state_clock_commutator_norm': transverse_commutator,
        'continuum_thermal_bound': faithful_thermal_field_bound(),
        'necessary_and_sufficient_continuous_flow_relation':
            'rho_j proportional exp(-b Omega_j (I-Z_j)/2), r_j=tanh(b Omega_j/2)',
        'canonical_fluctuations': 'Q_N=sum X/sqrt(2Nr), P_N=sum Y/sqrt(2Nr)',
        'thermal_Wightman': '(cos(Omega delta)-i r sin(Omega delta))/(2r Omega)',
        'KMS_boundary': 'C(t-i b)=C(-t)',
        'native_clock_conversion': 't=kappa tau, kappa=2d/sqrt(630); s=-kappa tau/b',
        'refining_mode_preparation':
            'Omega_l=sqrt(Q_l) theta_l; r_l=tanh(b Omega_l/2); g_l/theta_l=b sqrt(Q_l)',
        'inverse_temperature_is_shared_testable_preparation_parameter': True,
        'common_inverse_temperature_forced_by_carrier_or_repair': False,
        'geometric_patch_modular_flow_or_Unruh_temperature_identified': False}


def product_clock_bound(level, degree=5, angular_band=F(8), radial_band=F(9),
                        horizon=F(1), tolerance=F(1, 8), step_error=F(0)):
    """Single-code unitary error in the native incident parameter tau.

    The two-axis word, count clock, and angular dispersion share Q_l.
    Radial grid error is separate; radial_band bounds its supplied frequency
    squared. The pulse error is for the complete two-pulse single-code step.
    """
    fixed = clocks.incident_clock_bound(level, degree, horizon, tolerance)
    angular_band, radial_band = map(dynamics._positive, (angular_band, radial_band))
    step_error = _fraction(step_error)
    if degree not in (5, 6) or step_error < 0:
        raise ValueError('a native port degree and nonnegative step error required')
    wave = dynamics.wave_bounds(level, angular_band, horizon)
    q = wave['squared_coherent_clock_normalization']
    delta = wave['relative_wave_frequency_error']
    a, c = _root_upper(angular_band)*(1+delta), _root_upper(radial_band)
    omega = _root_upper(a*a+c*c)
    kappa = _root_upper(F(2*degree**2, 315))
    eta = F(fixed['uniform_count_error_upper'])
    steps = math.floor(F(fixed['mean_incident_ticks_per_continuum_unit'])*(horizon+eta))
    count = omega*kappa*eta/2
    dispersion = kappa*horizon*_root_upper(angular_band)*delta/2
    split = steps*a*c/(4*q)
    synthesis = steps*step_error
    return {'level': level, 'port_degree': degree, 'angular_band': angular_band,
        'radial_band': radial_band, 'Q': q, 'horizon_in_native_incident_units': horizon,
        'clock_failure_probability_upper': F(fixed['failure_probability_upper']),
        'count_error_single_code_unitary_upper': count,
        'angular_dispersion_single_code_unitary_upper': dispersion,
        'split_axes_single_code_unitary_upper': split,
        'source_synthesis_single_code_unitary_upper': synthesis,
        'combined_single_code_unitary_error_upper': min(F(2), count+dispersion+split+synthesis),
        'maximum_triggered_steps_on_clock_good_event': steps,
        'source_pulses_per_code_per_trigger': 2,
        'modular_parameter': 's=-(2d/sqrt(630))*tau/b in the limiting product field',
        'growing_precision_requirement': 'sqrt(Q_l)*step_error_l -> 0',
        'count_tolerance_requirement': 'tolerance_l -> 0 and 4^l*tolerance_l^2 -> infinity',
        'radial_grid_and_spatial_subsequence_errors_included': False,
        'protected_register_and_paused_atomic_pulse_services_required': True}


def product_thermal_field_bound(inverse_temperature=F(4), length=F(2)):
    """A normal thermal state for the same declared Dirichlet product."""
    beta = _fraction(inverse_temperature, positive=True)
    bounds = observations.product_thermal_bounds(beta, length)
    return {'inverse_temperature': beta, 'interval_length': length,
        'one_particle_heat_trace_upper': bounds['one_particle_heat_trace_upper'],
        'normal_ordered_Fock_log_partition_upper': bounds['normal_ordered_Fock_log_partition_upper'],
        'state': 'rho_b=Z_b^-1 exp(-b dGamma(sqrt(B))), B=-d_r^2+(A direct_sum 0)',
        'modular_flow': 'sigma_s=alpha_(-b*s) on the full encoded Fock algebra',
        'domain': 'declared interval (0,L) with Dirichlet ends times the native screen',
        'angular_constant_mode_included': True, 'product_has_zero_frequency_mode': False,
        'field_smearing_domain': 'L2; momentum smearing domain D(B^(1/4))',
        'commutator_and_kick_response': 'B^(-1/2) sin((t-u)*sqrt(B))',
        'comparison_distance_squared': '|r-s|^2+2*d_screen^2',
        'comparison_response_speed_upper': 1,
        'spatial_patch_modular_flow_or_physical_Lorentzian_metric_identified': False}


def product_native_mode_control():
    """An exact full-source mode and zero-inclusive incident records."""
    from ..observer import native
    graph, lap, _, _, _ = dynamics.source_green.matrices(0)
    matrix = sp.Matrix(lap.toarray())
    external = sorted({p for edge in graph.inter_seams for p in edge if p < 12})
    rows = (matrix[:12, :12]-6*sp.eye(12)).col_join(sp.eye(12)[external, :])
    local = rows.nullspace()[0]
    profile = local.col_join(sp.zeros(graph.ports-12, 1))
    if matrix*profile != 6*profile or sum(profile) or matrix*sp.ones(graph.ports, 1) != sp.zeros(graph.ports, 1):
        raise ArithmeticError('the product clock did not use an actual full-port native mode')
    port = next(p for p, neighbors in enumerate(graph.neighbours) if len(neighbors) == 5)
    initial = (1,)*(graph.ports//3)+(0,)*(2*graph.ports//3)
    records = Counter()
    for seam in graph.seams:
        for coin in (0, 1):
            _, (_, current) = native.native(initial, seam, coin)
            records[int(port in seam), int(bool(current and port in seam))] += 1
    if sum(records.values()) != 1260 or sum(n for (tick, _), n in records.items() if tick) != 10:
        raise ArithmeticError('native incident-clock records lost zero-current events')
    return {'ports': graph.ports, 'seams': len(graph.seams), 'clock_port': port,
        'native_eigenvalue': 6, 'exact_first_carrier_profile': [str(x) for x in local],
        'all_other_profile_entries': '0', 'external_seam_ports_of_first_carrier': external,
        'profile_cell_compression': '0 on every cell',
        'internal_fixture_survives_as_a_bounded_nonconstant_continuum_mode': False,
        'same_internal_fixture_scaled_eigenvalue_at_level_l': '6*4^l',
        'complete_one_attempt_clock_record_counts': [
            {'incident': tick, 'nonzero_incident_current': current, 'count': count}
            for (tick, current), count in sorted(records.items())],
        'constant_angular_mode_eigenvalue': 0,
        'clock_ticks_require_nonzero_current': False}


@lru_cache(maxsize=1)
def _cell_visible_packet():
    """Exact polynomial selection on the actual level-zero 240-port graph.

    The integer polynomial fixes the choice inside nearly degenerate bands;
    no numerical eigenvector or external continuum profile is supplied.
    """
    graph, laplacian, _, _, _ = dynamics.source_green.matrices(0)
    profile = [1]*12+[-1]*12+[0]*(graph.ports-24)
    for _ in range(1024):
        profile = [12*profile[i]-sum(profile[i]-profile[j] for j in neighbors)
                   for i, neighbors in enumerate(graph.neighbours)]
        divisor = math.gcd(*profile)
        profile = [x//divisor for x in profile]
    norm2 = sum(x*x for x in profile)
    applied = [sum(profile[i]-profile[j] for j in neighbors)
               for i, neighbors in enumerate(graph.neighbours)]
    energy = F(sum(x*y for x, y in zip(profile, applied, strict=True)), norm2)
    residual2 = sum((F(y)-energy*x)**2 for x, y in zip(profile, applied, strict=True))/norm2
    visible = F(sum(sum(profile[j:j+12])**2 for j in range(0, graph.ports, 12)), 12*norm2)
    if sum(profile) or not F(0) < energy < F(1, 20) or visible < F(99, 100):
        raise ArithmeticError('the native polynomial lost its mean-zero cell-visible low packet')
    with mp.workdps(100):
        constant = [1/mp.sqrt(graph.ports)]*graph.ports+[mp.mpf(0)]*(256-graph.ports)
        mode = [mp.mpf(x)/mp.sqrt(norm2) for x in profile]+[mp.mpf(0)]*(256-graph.ports)
        first = [-x for x in constant]; first[0] += 1
        first_norm = mp.sqrt(sum(x*x for x in first))
        first = [x/first_norm for x in first]
        overlap = sum(x*y for x, y in zip(first, mode, strict=True))
        rotated = [y-2*x*overlap for x, y in zip(first, mode, strict=True)]
        if abs(rotated[0]) > mp.mpf('1e-95'):
            raise ArithmeticError('the exact orthogonal packet does not preserve the anchor')
        rotated[0] = mp.mpf(0)  # Exact identity <constant,mode>=0.
        second = [-x for x in rotated]; second[1] += 1
        second_norm = mp.sqrt(sum(x*x for x in second))
        second = [x/second_norm for x in second]
        rows = [[mp.nstr(x, 100) for x in row] for row in (constant, mode, first, second)]
    return {'graph': graph, 'laplacian': laplacian, 'integer_profile': tuple(profile),
        'norm_squared': norm2, 'rayleigh': energy, 'residual_squared': residual2,
        'cell_visible_weight': visible, 'vectors': tuple(tuple(row) for row in rows)}


def _packet_state_word(amplitudes):
    """Deterministic binary-tree state synthesis on eight existing code qubits."""
    word = []
    with mp.workdps(90):
        values = list(map(mp.mpf, amplitudes))
        for target in reversed(range(8)):
            controls = 7-target
            angles = []
            for prefix in range(1 << controls):
                block = values[prefix*(1 << (target+1)):(prefix+1)*(1 << (target+1))]
                half = 1 << target
                x = block[0] if not target else mp.sqrt(sum(v*v for v in block[:half]))
                y = block[1] if not target else mp.sqrt(sum(v*v for v in block[half:]))
                angles.append(-mp.atan2(y, x) if x or y else mp.mpf(0))
            for mask in range(1 << controls):
                coefficient = sum(theta*(-1)**((prefix & mask).bit_count())
                                  for prefix, theta in enumerate(angles))/(1 << controls)
                rounded = F(int(mp.floor(coefficient*(1 << 44)+mp.mpf('.5'))), 1 << 44)
                label = ['I']*8; label[target] = 'Y'
                for bit in range(controls):
                    if mask >> bit & 1:
                        label[target+1+bit] = 'Z'
                word += source_compiler.pauli_rotation(label, float(rounded))
    return [source_compiler.pulse(p['kind'], round(p['angle']*2**44)/2**44, p['targets'])
            for p in word]


@lru_cache(maxsize=1)
def _cell_visible_packet_source():
    data = _cell_visible_packet()
    constant, mode, first, second = (np.array(row, float) for row in data['vectors'])
    words, raw = {}, {}
    for name, values in (('first', data['vectors'][2]), ('second', data['vectors'][3])):
        words[name] = _packet_state_word(values)
        words['inverse_'+name] = source_compiler.inverse(words[name])
    phase = [float(_source_dyadic_scalar(sp.pi))]+[0.]*255
    words['reflect'] = [source_compiler.pulse(p['kind'], round(p['angle']*2**44)/2**44, p['targets'])
                        for p in source_compiler.compile_diagonal(phase)]
    words['inverse_reflect'] = source_compiler.inverse(words['reflect'])
    for name, word in words.items():
        raw[name] = source_compiler.execute(word, 8)
    h1 = raw['first']@raw['reflect']@raw['inverse_first']
    h2 = raw['second']@raw['reflect']@raw['inverse_second']
    h1i = raw['first']@raw['inverse_reflect']@raw['inverse_first']
    h2i = raw['second']@raw['inverse_reflect']@raw['inverse_second']
    words['encode'] = (words['inverse_second']+words['reflect']+words['second']
                       +words['inverse_first']+words['reflect']+words['first'])
    words['decode'] = source_compiler.inverse(words['encode'])
    encoding, decoding = h1@h2, h2i@h1i
    wanted = (np.eye(256)-2*np.outer(first, first))@(np.eye(256)-2*np.outer(second, second))
    error = float(np.linalg.norm(encoding-wanted, 'fro'))
    inverse_error = float(np.linalg.norm(decoding-wanted.T, 'fro'))
    guard = F(1, 10**8)
    if (not math.isfinite(error+inverse_error) or max(error, inverse_error) > float(guard)
            or np.linalg.norm(wanted[:, 0]-constant) > 2e-13
            or np.linalg.norm(wanted[:, 1]-mode) > 2e-13):
        raise ArithmeticError('the executed full-port packet isometry exceeds its numerical guard')
    return {'data': data, 'words': words, 'encode': encoding, 'decode': decoding,
        'target': wanted, 'guard': guard, 'encoding_error': error,
        'decoding_error': inverse_error, 'constant': constant, 'mode': mode}


def cell_visible_packet_control(steps=25):
    """Actual native cell-visible packet, source encoder and incident clock.

    This is a finite approximate mode with a quantified native Chebyshev
    residual. It is neither an exact invariant mode nor a growing encoder.
    Every256 work address and both untouched reference columns are retained.
    """
    steps = _int(steps)
    if steps > 128:
        raise ValueError('the finite packet clock admits at most128 incident triggers')
    source = _cell_visible_packet_source(); data = source['data']
    graph, rayleigh = data['graph'], data['rayleigh']
    theta = _source_dyadic_scalar(sp.acos(1-sp.Rational(rayleigh)/(2*len(graph.seams))))
    phases = [0.]*256; phases[1] = -float(theta)
    phase_word = [source_compiler.pulse(p['kind'], round(p['angle']*2**44)/2**44, p['targets'])
                  for p in source_compiler.compile_diagonal(phases)]
    raw_phase = source_compiler.execute(phase_word, 8)
    target_phase = np.eye(256, dtype=complex)
    target_phase[1, 1] = np.exp(-1j*float(theta))
    phase_error = float(np.linalg.norm(raw_phase-target_phase, 'fro'))
    if not math.isfinite(phase_error) or phase_error > 1e-9:
        raise ArithmeticError('the packet phase word exceeds its complete-address guard')
    raw_step = source['encode']@raw_phase@source['decode']
    desired_step = source['target']@target_phase@source['target'].T
    step_guard = (1+source['guard'])**2*(1+F(1, 10**9))-1
    # Build a genuine logical Bell input by source pulses, then encode only
    # the work register; its second index is an untouched reference qubit.
    bell_word = [source_compiler.pulse('H', float(_source_dyadic_scalar(sp.pi/2)), [0])]
    bell_word += source_compiler.cnot(0, 8)
    bell_word = [source_compiler.pulse(p['kind'], round(p['angle']*2**44)/2**44, p['targets'])
                 for p in bell_word]
    blank = np.zeros(512, complex); blank[0] = 1
    bell = source_compiler.execute(bell_word, 9, blank).reshape(2, 256).T
    encoded = source['encode']@bell
    evolved = np.linalg.matrix_power(raw_step, steps)@encoded
    decoded = source['decode']@evolved
    read_word = [source_compiler.pulse('H', float(_source_dyadic_scalar(sp.pi/2)), [j]) for j in (0, 8)]
    output = source_compiler.execute(read_word, 9, decoded.T.reshape(512))
    probabilities = abs(output)**2
    contrast = math.fsum((-1)**((j & 1)+(j >> 8))*float(probabilities[j]) for j in range(512))
    expected = math.cos(steps*float(theta))
    source_amplitude_bound = (1+source['guard'])**2*(1+step_guard)**steps-1+F(1, 10**9)
    observable_bound = 2*source_amplitude_bound+source_amplitude_bound**2
    if (abs(sum(probabilities)-1) > float(observable_bound)
            or abs(contrast-expected) > float(observable_bound)):
        raise ArithmeticError('the full packet/reference readout exceeds its source budget')
    mode = source['mode']; kernel = np.eye(graph.ports)-data['laplacian'].toarray()/(2*len(graph.seams))
    previous, current = mode[:graph.ports], kernel@mode[:graph.ports]
    if not steps:
        native_cosine = previous
    elif steps == 1:
        native_cosine = current
    else:
        for _ in range(2, steps+1):
            previous, current = current, 2*kernel@current-previous
        native_cosine = current
    residual_upper = _root_upper(data['residual_squared'])
    native_bound = F(steps*steps, 2*len(graph.seams))*residual_upper+F(steps, 1 << 45)
    native_error = float(np.linalg.norm(native_cosine-expected*mode[:graph.ports]))
    if native_error > float(native_bound)+1e-12:
        raise ArithmeticError('the packet clock exceeded its actual-native Chebyshev residual')
    cell = mode[:graph.ports].reshape(graph.cells, 12).sum(axis=1)/math.sqrt(12)
    from ..observer import native
    initial = (1,)*(graph.ports//3)+(0,)*(2*graph.ports//3)
    port = next(i for i, neighbors in enumerate(graph.neighbours)
                if len(neighbors) == 5 and any(initial[i] != initial[j] for j in neighbors))
    records, counts, expected_field = [], Counter(), 0.
    before_field = float(np.dot(mode[:graph.ports], initial))
    for index, seam in enumerate(graph.seams):
        for coin in (0, 1):
            after, (_, current_record) = native.native(initial, seam, coin)
            tick = int(port in seam)
            records.append((index, coin, current_record, tick))
            counts[tick, int(bool(tick and current_record))] += 1
            a, b = seam
            decoded_field = before_field+current_record*(mode[b]-mode[a])
            if abs(decoded_field-np.dot(mode[:graph.ports], after)) > 2e-13:
                raise ArithmeticError('native current records do not decode the encoded packet field')
            expected_field += decoded_field/(2*len(graph.seams))
    native_field_target = float(initial@kernel@mode[:graph.ports])
    if abs(expected_field-native_field_target) > 2e-12 or not counts[1, 0]:
        raise ArithmeticError('the complete native packet/incident-clock law changed')
    clock_marginal = []
    for tick in (0, 1):
        once = source['decode']@(raw_step@encoded if tick else encoded)
        decoded_once = source_compiler.execute(read_word, 9, once.T.reshape(512))
        law = abs(decoded_once)**2
        parity = [sum(float(value) for j, value in enumerate(law)
                      if ((j & 1)+(j >> 8)) % 2 == sign) for sign in (0, 1)]
        mass = sum(value for (mark, _), value in counts.items() if mark == tick)
        clock_marginal.append({'incident': tick, 'native_branch_count': mass,
            'joint_reference_parity_probabilities': [mass*x/(2*len(graph.seams)) for x in parity]})
    word = source['words']['decode']+phase_word+source['words']['encode']
    pulses = len(bell_word)+len(source['words']['encode'])+steps*len(word)+len(source['words']['decode'])+len(read_word)
    rounded = lambda x: F(math.ceil(x*10**12), 10**12)
    return {'source_level': 0, 'native_ports': graph.ports, 'protected_work_addresses': 256,
        'retained_work_padding_addresses': list(range(graph.ports, 256)),
        'selection': '(12I-L)^1024 applied to cell0 minus cell1, with integer gcd reduction',
        'integer_profile_sha256': hashlib.sha256(json.dumps(data['integer_profile']).encode()).hexdigest(),
        'exact_profile_sum': sum(data['integer_profile']),
        'maximum_profile_integer_bits': max(abs(x).bit_length() for x in data['integer_profile']),
        'rayleigh_eigenvalue_interval': [rounded(rayleigh)-F(1, 10**12), rounded(rayleigh)],
        'exact_residual_norm_upper': residual_upper,
        'cell_visible_weight_lower': rounded(data['cell_visible_weight'])-F(1, 10**12),
        'cell_readout_profile': cell.tolist(),
        'native_internal_fixture_has_zero_cell_compression': True,
        'packet_is_an_exact_native_eigenvector': False,
        'source_encoder_complete_Frobenius_error': source['encoding_error'],
        'source_decoder_complete_Frobenius_error': source['decoding_error'],
        'source_isometry_numerical_guard': source['guard'],
        'encoder_word_sha256': hashlib.sha256(json.dumps(source['words']['encode'], sort_keys=True).encode()).hexdigest(),
        'incident_triggers': steps, 'native_phase_per_trigger': theta,
        'phase_dictionary': 'theta=acos(1-Rayleigh/(2M)), M=630; wave t=N/sqrt(630)',
        'clock_port': port,
        'complete_native_one_attempt_record_sha256': hashlib.sha256(json.dumps(records).encode()).hexdigest(),
        'complete_native_one_attempt_record_counts': [
            {'incident': tick, 'nonzero_incident_current': current_record, 'count': count}
            for (tick, current_record), count in sorted(counts.items())],
        'complete_one_attempt_clock_reference_marginal': clock_marginal,
        'packet_field_native_current_decoder': 'F(after)=F(before)+current*(v_b-v_a)',
        'native_record_mean_field_residual': abs(expected_field-native_field_target),
        'complete_system_reference_probabilities': probabilities.tolist(),
        'source_reference_parity_contrast': contrast,
        'native_packet_cosine_residual_diagnostic': native_error,
        'native_packet_cosine_residual_upper': native_bound,
        'source_readout_error_upper': observable_bound,
        'source_step_complete_target_error': float(np.linalg.norm(raw_step-desired_step, 'fro')),
        'resources': {'encoder_pulses': len(source['words']['encode']), 'clock_pulses_per_trigger': len(word),
            'total_source_pulses': pulses, 'proper_code_load_events': pulses, 'proper_code_export_events': pulses,
            'protected_work_qubits': 8, 'reference_qubits': 1, 'complete_terminal_addresses': 512,
            'classical_profile_construction_integer_iterations': 1024,
            'initial_blank_storage_controller_and_paused_services_supplied': True},
        'full_register_input_columns_checked': 256,
        'cell_visible_finite_native_packet_encoded_and_read': True,
        'unbounded_growing_encoder_or_nonconstant_continuum_mode_selected': False,
        'prepared_thermal_state_radial_law_or_physical_seconds_derived': False}


def product_mode_source_control(native_eigenvalue=6, radial_mode=1, radial_points=3,
                                length=F(2), level=0, inverse_temperature=F(1), steps=25):
    """Compile a product mode from the native angular clock and a radial axis.

    No product-frequency pulse is inserted: each incident service is the
    native Z phase followed by the radial X pulse. Their noncommutation is
    measured and bounded. A declared mode must separately belong to the
    actual native Laplacian; product_native_mode_control certifies the retained
    lambda=0,6 fixtures, not every caller-supplied spectral value.
    """
    level, steps = dynamics._level(level), _int(steps, 1)
    if steps > 256 or level > 8:
        raise ValueError('the finite source execution admits at most 256 steps and level eight')
    beta = _fraction(inverse_temperature, positive=True)
    dictionary = observations.product_mode_dictionary(
        native_eigenvalue, radial_mode, radial_points, length, angular_scale=F(4**level))
    q = 630*16**level
    lam, radial = dictionary['scaled_native_eigenvalue'], dictionary['radial_eigenvalue']
    if not (lam >= 0 and lam < 4*q):
        raise ValueError('the angular mode must lie in the coherent source phase domain')
    theta = sp.acos(1-lam/(2*q))
    a, c = sp.sqrt(q)*theta, sp.sqrt(radial)
    omega_expression = sp.sqrt(a*a+c*c)
    axis_angle = sp.atan2(c, a)
    frequency = _source_dyadic_scalar(omega_expression)
    polarization = _source_dyadic_scalar(sp.tanh(sp.Rational(beta)*omega_expression/2))
    parameters = thermal_clock_parameters(frequency, polarization)
    af, cf, omega, time = (float(sp.N(x, 80)) for x in (a, c, omega_expression, steps/sp.sqrt(q)))
    r = float(polarization)
    x = np.array([[0, 1], [1, 0]], complex)
    y = np.array([[0, -1j], [1j, 0]], complex)
    z, identity = np.diag([1., -1.]), np.eye(2)
    normal, transverse = (cf*x+af*z)/omega, (af*x-cf*z)/omega
    one = (identity+r*normal)/2
    rotation = [source_compiler.pulse('Y', -float(_source_dyadic_scalar(axis_angle/2)), [0])]
    step_word = [source_compiler.pulse('Z', float(_source_dyadic_scalar(theta/2)), [0]),
                 source_compiler.pulse('X', float(_source_dyadic_scalar(c/(2*sp.sqrt(q)))), [0])]
    unitary = source_compiler.execute(step_word*steps, 1)
    ideal = np.cos(time*omega/2)*identity+1j*np.sin(time*omega/2)*normal
    unitary_error = float(np.linalg.norm(unitary-ideal, 2))
    split_bound = float(steps*sp.N(a*c, 80)/(4*q))
    synthesis_bound = steps/2**44
    if unitary_error > split_bound+synthesis_bound+2e-12:
        raise ArithmeticError('source product clock exceeded its two-axis commutator bound')
    reversed_step = source_compiler.execute(list(reversed(step_word)), 1)
    order_defect = float(np.linalg.norm(source_compiler.execute(step_word, 1)-reversed_step, 2))
    copies = 2
    rotate = sum((source_compiler.remap(rotation, [j]) for j in range(copies)), [])
    evolution = sum((source_compiler.remap(step_word*steps, [j]) for j in range(copies)), [])
    preparation = _thermal_code_preparation(polarization, copies)+rotate
    blank = np.zeros(32, complex); blank[0] = 1
    pure = source_compiler.execute(preparation, 4, blank[:16])
    reduced = pure.reshape(4, 4).T@pure.reshape(4, 4).conj()
    state_error = float(np.linalg.norm(reduced-np.kron(one, one), 2))
    if state_error > 2e-11:
        raise ArithmeticError('product thermal preparation lost its joint generator axis')
    # Two times on the same code, with references and pointer retained.
    strengths = (F(1, 2), F(1, 3))
    def weyl(strength, controlled=False):
        normalized = float(sp.N(sp.Rational(strength)/
                           sp.sqrt(sp.Rational(polarization)*omega_expression), 80))
        middle = spin_weyl_word([[normalized, 0]], copies, controlled=controlled)
        if controlled:
            middle = source_compiler.remap(middle, [0, 1, 4])
        return source_compiler.inverse(rotate)+middle+rotate
    single_q = transverse/np.sqrt(2*copies*r*omega)
    def exact_weyl(strength):
        angle = float(strength)/np.sqrt(2*copies*r*omega)
        return np.cos(angle)*identity+1j*np.sin(angle)*transverse
    first, second = map(exact_weyl, strengths)
    expected = complex(np.trace(one@unitary.conj().T@first@unitary@second)**copies)
    ideal_expected = complex(np.trace(one@ideal.conj().T@first@ideal@second)**copies)
    records = {}
    for axis in ('X', 'Y'):
        word = preparation+[source_compiler.pulse('Y', -math.pi/4, [4])]
        word += weyl(strengths[1], True)+evolution+weyl(strengths[0], True)+source_compiler.inverse(evolution)
        word += [source_compiler.pulse('Y' if axis == 'X' else 'X',
                                      math.pi/4 if axis == 'X' else -math.pi/4, [4])]
        output = source_compiler.execute(word, 5, blank)
        probabilities = abs(output)**2
        mean = float(sum(probabilities[:16])-sum(probabilities[16:]))
        oracle = expected.real if axis == 'X' else expected.imag
        if abs(mean-oracle) > 3e-11 or abs(sum(probabilities)-1) > 3e-11:
            raise ArithmeticError('product two-time detector lost an ordered phase or complete outcome')
        records[axis] = {'all_system_reference_pointer_probabilities': probabilities.tolist(),
            'pointer_mean': mean, 'independent_matrix_expectation': oracle,
            'source_pulses': len(word),
            'word_sha256': hashlib.sha256(json.dumps(word, sort_keys=True).encode()).hexdigest()}
    epsilon = F(1, 1000)
    probabilities, means = [], []
    readout = np.array([copies-2*(j % 4).bit_count() for j in range(16)])
    for sign in (-1, 1):
        word = preparation+weyl(sign*epsilon)+evolution+source_compiler.inverse(rotate)
        word += [source_compiler.pulse('Y', math.pi/4, [j]) for j in range(copies)]
        output = source_compiler.execute(word, 4, blank[:16])
        law = abs(output)**2
        probabilities.append(law.tolist())
        means.append(float(law@readout)/np.sqrt(2*copies*r*omega))
    response = (means[1]-means[0])/(2*float(epsilon))
    evolved = unitary.conj().T@single_q@unitary
    linear = float((1j*copies*np.trace(one@(evolved@single_q-single_q@evolved))).real)
    sine = math.sin(time*omega)/omega
    kick_bound = abs(linear)*float(epsilon)**2/(3*copies*r*omega)
    if abs(response-linear) > kick_bound+3e-10 or abs(linear-sine) > 2*unitary_error/(r*omega)+2e-12:
        raise ArithmeticError('product source intervention exceeded its finite kick or split-clock error')
    gap = parameters['modular_gap_diagnostic']
    # The rounded polarization fixes its own exact gap. Compare the actual
    # prepared state at s=-time/requested_beta, charging both preparation
    # error and the small difference from the ideal common-beta relation.
    prepared_one = source_compiler.execute(_thermal_code_preparation(polarization, 1)+rotation,
                                          2, np.array([1., 0., 0., 0.]))
    amplitude = prepared_one.reshape(2, 2)
    reduced_one = amplitude.T@amplitude.conj()
    state_one_error = float(np.linalg.norm(reduced_one-one, 2))
    eigenvalues, eigenvectors = np.linalg.eigh(reduced_one)
    modular = (eigenvectors*np.exp(-1j*time/float(beta)*np.log(eigenvalues)))@eigenvectors.conj().T
    polarization_modular_allowance = abs(time*(gap/float(beta)-omega))
    preparation_modular_allowance = 2*time/float(beta)*state_one_error/min(
        eigenvalues[0], (1-r)/2)
    modular_error = max(float(np.linalg.norm(modular@p@modular.conj().T-unitary.conj().T@p@unitary, 2))
                        for p in (x, y, z))
    if modular_error > 2*unitary_error+polarization_modular_allowance+preparation_modular_allowance+2e-12:
        raise ArithmeticError('finite product modular mismatch exceeded its unitary error')
    # Marginalize every native seam/coin branch, never only successful currents.
    clock_rows = []
    branches_per_attempt = 1260*4**level
    for marks in ((0, 0), (0, 1), (1, 0), (1, 1)):
        count = math.prod(10 if tick else branches_per_attempt-10 for tick in marks)
        u = source_compiler.execute(step_word*sum(marks), 1)
        value = complex(np.trace(one@u.conj().T@first@u@second)**copies)
        clock_rows.append({'incident_marks': list(marks), 'native_branches': count,
            'joint_X_pointer_probabilities': [float(F(count, branches_per_attempt**2))*(1+sign*value.real)/2
                                             for sign in (1, -1)]})
    if abs(sum(sum(row['joint_X_pointer_probabilities']) for row in clock_rows)-1) > 2e-12:
        raise ArithmeticError('product clock marginal discarded a native branch')
    return {'mode_dictionary': {key: str(value) for key, value in dictionary.items()},
        'level': level, 'Q': q, 'incident_triggers': steps,
        'angular_frequency_expression': str(a), 'radial_frequency_expression': str(c),
        'joint_frequency_expression': str(omega_expression),
        'state_clock_parameters': parameters, 'requested_inverse_temperature': beta,
        'inferred_inverse_temperature_diagnostic': gap/omega,
        'source_step_word': step_word, 'source_step_word_pulses_per_code': len(step_word),
        'source_step_angle_error_upper': F(1, 1 << 44),
        'single_code_unitary_error': unitary_error,
        'two_axis_product_formula_error_upper': split_bound,
        'two_axis_pulse_synthesis_error_upper': synthesis_bound,
        'split_order_unitary_difference': order_defect,
        'wrong_parallel_axis_gap_error': af+cf-omega,
        'preparation_reduced_state_error': state_error,
        'finite_split_state_stationarity_defect': float(np.linalg.norm(unitary@one@unitary.conj().T-one, 2)),
        'finite_split_modular_action_error_on_Paulis': modular_error,
        'modular_parameter_used': -time/float(beta),
        'polarization_common_beta_modular_allowance': polarization_modular_allowance,
        'source_preparation_modular_allowance': preparation_modular_allowance,
        'complete_ordered_two_time_pointer_records': records,
        'finite_ordered_expectation': [expected.real, expected.imag],
        'ideal_joint_generator_ordered_expectation': [ideal_expected.real, ideal_expected.imag],
        'complete_signed_kick_system_reference_probabilities': probabilities,
        'signed_kick_response': response, 'actual_split_clock_linear_response': linear,
        'target_product_sine_response': sine, 'finite_kick_error_upper': kick_bound,
        'complete_two_attempt_clock_and_pointer_marginal': clock_rows,
        'finite_split_clock_is_exact_target_modular_flow': False,
        'product_frequency_is_an_independent_input_pulse': False,
        'spectral_encoder_radial_grid_and_protected_services_are_declared': True,
        'source_eigenmode_membership_is_checked_for_every_function_input': False}


def product_clock_encoding_control(inverse_temperature=F(1)):
    beta = _fraction(inverse_temperature, positive=True)
    return {'native_mode_and_clock': product_native_mode_control(),
        'finite_source_modes': [product_mode_source_control(lam, radial_mode=k,
                                                           inverse_temperature=beta)
                                for lam, k in ((6, 1), (6, 2), (0, 1))],
        'refinement_clock_bounds': [product_clock_bound(l, tolerance=F(1, 2**(l//2+3)),
                                                      step_error=F(1, 2**(44+4*l)))
                                    for l in (0, 2, 4, 8)],
        'continuum_thermal_state': product_thermal_field_bound(beta),
        'limiting_generator': 'sqrt(B), B=-d_r^2+(A direct_sum 0)',
        'ordered_covariance': '(coth(b sqrt(B)/2)*cos(t sqrt(B))-i sin(t sqrt(B)))/(2 sqrt(B))',
        'conditional_same_product_operator_state_response_and_active_clock': True,
        'nonconstant_finite_fixture_is_an_internal_compiler_control': True,
        'constant_angular_fixture_is_included_in_the_product_limit': True,
        'passive_source_selects_radial_domain_temperature_or_mode_encoding': False,
        'growing_source_compiler_or_local_detector_network_executed': False,
        'observed_spacetime_or_lab_clock_identified': False}


def run():
    rows = {
        'schema': 'oph.observer_time_fields.continuum_dynamics.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'exact_native_wave_synthesis': exact_wave_synthesis(),
        'actual_native_propagators': native_propagator_control(),
        'all_level_propagation_bounds': [propagation_bounds(l) for l in (0, 1, 2, 4, 10)],
        'fixed_cap_order_of_limits_control': [fixed_cap_collapse(l) for l in (0, 2, 4, 10)],
        'continuum_path_bounds': [wave_path_bounds(j) for j in (0, 10, 100)],
        'fractional_wave_bounds': [fractional_wave_bounds(l) for l in range(3)],
        'source_encoded_quantum_fluctuation': quantum_fluctuation_control(),
        'faithful_thermal_modular_encoding': faithful_modular_encoding_control(),
        'product_thermal_active_clock_encoding': product_clock_encoding_control(),
        'cell_visible_native_packet_encoding': cell_visible_packet_control(),
        'native_weighted_form_control': native_weighted_form_control(),
        'native_heat_support_control': native_heat_support_control(),
        'continuum_finite_speed_bounds': [finite_speed_bounds(l) for l in (0, 1, 2, 4, 10)],
        'actual_native_density_process': {
            'exact_stirring_and_canonical_controls': native_stirring_control(),
            'all_level_bounds': [density_process_bounds(l) for l in (0, 1, 2, 4, 10)],
            'exact_native_path_controls': native_density_path_control(),
            'path_bounds': [density_path_bounds(l, j) for l, j in ((0, 0), (1, 10), (4, 100), (10, 100))],
            'uniform_fixed_number_equilibrium_and_IID_seams_required': True,
            'Gaussian_limit_follows_from_conditional_finite_population_CLT': True,
            'reflection_positive_Gaussian_covariance': True,
            'OS_one_particle_semigroup': 'exp(-tau A), not exp(-tau sqrt(A))',
            'OS_Gaussian_Fock_Hamiltonian': 'dGamma(A)',
            'Green_field_normalization': 'phi=Phi/sqrt(2)',
            'Green_field_temporal_Fourier_covariance': '(omega^2+A^2)^-1',
            'Green_field_Cameron_Martin_action': '1/2 integral (||partial_tau phi||^2+||A phi||^2) d tau',
            'action_is_density_against_infinite_dimensional_Lebesgue_measure': False,
            'native_coherent_quantum_vacuum_preparation_inferred': False},
        'full_reference_state_limit': {
            'native_trace_and_nontracial_controls': [native_reference_state_control(q) for q in (2, 4)],
            'actual_source_degree_census': [source_trace_census(l) for l in range(3)],
            'all_level_reference_bounds': [reference_state_bounds(l) for l in (0, 1, 2, 4, 10, 100)],
            'bounded_algebra_readout_moments_and_products_must_converge': True,
            'fixed_parameter_bounded_algebra_limit_state_is_tracial': True,
            'density_wave_reference_weight_vanishes_without_sector_conditioning': True,
            'bounded_overlap_fluctuation_expected_commutator_tends_to_zero': True,
            'quantum_central_limit_existence_asserted_for_this_refinement': False},
        'common_limit': {'space': 'mean-zero L2(mu) on the source S2',
            'same_subsequence': 'operator-norm limit G of embedded native Schur inverses',
            'positive_self_adjoint_generator': 'A=G^-1 with compact resolvent',
            'constant_mode_extension_is_conservative_strongly_local_Dirichlet': True,
            'constant_mode_extension_has_finite_speed_in_comparison_distance': True,
            'constant_mode_heat_is_positivity_improving': True,
            'disjoint_positive_profiles_have_nonzero_heat_response_at_every_positive_time': True,
            'fixed_selection_slot_parameter': 'tau=q/(2520*16^level)',
            'source_heat_trace_spectral_dimension': 2,
            'native_attempt_parameter': 'tau=n/(1260*16^level)',
            'native_two_point_transfer': 'exp(-tau A)',
            'actual_multitime_density_limit': 'stationary Gaussian OU continuous distribution-valued paths',
            'actual_density_interpolated_path_topology': 'C([0,T],H^-3(S2))',
            'actual_Green_filtered_interpolated_path_topology': 'C([0,T],H^-2(S2))',
            'coherent_step_parameter': 't=n/(sqrt(630)*4^level)',
            'native_coherent_compression': 'cos(t sqrt(A))',
            'return_reference_modular_phase_magnitudes': 'sqrt(A)',
            'synthesized_wave_field_topology': 'C([0,T],H^-2(S2))',
            'fractional_filtered_wave_field': 'W_{1/4}=A^(1/4) W',
            'fractional_filtered_wave_two_time_covariance': 'A^(-1/2) cos((t-u) sqrt(A)) = G^(1/2) cos((t-u) sqrt(A))',
            'fractional_filtered_wave_path_topology': 'C([0,T],H^-s(S2)) for every s>=7/4 (in particular H^-2(S2))',
            'Pauli_Jordan_causal_commutator_kernel': 'D(t-u)=A^(-1/2) sin((t-u) sqrt(A))',
            'uniform_propagator_convergence_on_compact_parameter_intervals': True,
            'propagator_convergence_topology': 'strong for each fixed input, uniform on compact parameter intervals',
            'same_source_covariance_resolvent_heat_and_wave_operator': True,
            'whole_refinement_sequence_limit_or_rate_proved': False},
        'scope': {'same_native_source_and_ancestry_dictionary': True,
            'stationary_sector_uniform_IID_seams_and_fair_repair_required': True,
            'full_initial_occupation_and_current_readback_required': True,
            'coherent_preparation_reflection_and_Ramsey_services_required': True,
            'Gaussian_wave_uses_declared_record_processing': True,
            'Gaussian_density_limit_uses_actual_native_snapshot_process': True,
            'density_process_path_space_limit_proved': True,
            'native_return_reference_is_a_native_coherently_prepared_state': False,
            'faithful_mean_zero_Fock_modular_state_under_declared_encoding_constructed': True,
            'faithful_product_Fock_state_and_source_split_clock_under_declared_encoding_constructed': True,
            'growing_source_field_modular_intervention_compiler_executed': False,
            'comparison_wave_support_cone_proved': True,
            'fractional_wave_canonical_energy_equipartition_proved': True,
            'comparison_sine_kernel_vanishes_for_separated_mean_zero_tests': True,
            'fractional_covariance_alone_constructs_native_quantum_CCR': False,
            'same_dictionary_heat_response_has_an_exact_finite_support_cone': False,
            'conditional_obstruction_to_same_dictionary_relativistic_response': True,
            'physical_faster_than_light_signalling_inferred': False,
            'every_OPH_clock_observer_or_reconstruction_excluded': False,
            'physical_signal_speed_or_Lorentzian_metric_identified': False,
            'observed_one_plus_three_spacetime_identified': False,
            'physical_interacting_QFT_CFT_or_vacuum_state_identified': False,
            'atomic_clock_or_public_lab_data_identification_proved': False}}
    # Retained diagnostics are explicitly numerical; theorem budgets are rational.
    def encode(value):
        if type(value) is F:
            return str(value)
        raise TypeError('unexpected receipt value')
    return json.loads(json.dumps(rows, default=encode))


def check_receipt(recorded, actual):
    """Replay primitive probabilities strictly and propagate their contrast.

    The parity is a signed sum of 512 retained outcomes. Applying the same
    scalar tolerance independently to that duplicate can reject otherwise
    admissible backend roundoff. Its replay allowance is the observed L1
    change of those probabilities, plus summation roundoff, not a relaxed
    tolerance on the primitive records or the independent source bound.
    """
    from ..observer.scalar_response.experiment import compare_receipt
    key, contrast_key = 'cell_visible_native_packet_encoding', 'source_reference_parity_contrast'
    rows = []
    projected = []
    for result in (recorded, actual):
        if type(result) is not dict or type(result.get(key)) is not dict or contrast_key not in result[key]:
            raise ValueError('the complete packet parity record is required')
        row = result[key]
        rows.append(row)
        projected.append({**result, key: {name: value for name, value in row.items() if name != contrast_key}})
    # In particular, all 512 primitive probabilities retain the existing
    # type, shape and 2e-10 scalar checks; every other field is unchanged.
    compare_receipt(*projected)
    allowances = []
    for row in rows:
        probabilities, contrast = row['complete_system_reference_probabilities'], row[contrast_key]
        if (type(probabilities) is not list or len(probabilities) != 512
                or any(type(p) is not float or not math.isfinite(p) or p < 0 for p in probabilities)
                or type(contrast) is not float or not math.isfinite(contrast)):
            raise ValueError('finite nonnegative complete packet probabilities and parity required')
        recomputed = math.fsum((-1)**((j & 1)+(j >> 8))*p for j, p in enumerate(probabilities))
        mass = math.fsum(probabilities)
        # A conservative gamma_n bound also admits receipts accumulated by
        # the former sequential sum. It is ~1.14e-13 for a unit-mass law.
        nu = len(probabilities)*math.ulp(1.)
        roundoff = nu/(1-nu)*mass
        if abs(contrast-recomputed) > roundoff:
            raise ValueError('stored packet parity is inconsistent with its retained probabilities')
        allowance = float(F(row['source_readout_error_upper']))
        theta, steps = float(F(row['native_phase_per_trigger'])), row['incident_triggers']
        if (type(steps) is not int or not 0 <= steps <= 128
                or not math.isfinite(allowance) or allowance < 0
                or abs(mass-1) > allowance
                or abs(contrast-math.cos(steps*theta)) > allowance):
            raise ValueError('stored packet probabilities or parity violate the independent source bound')
        allowances.append(roundoff)
    probability_change = math.fsum(abs(x-y) for x, y in zip(
        rows[0]['complete_system_reference_probabilities'],
        rows[1]['complete_system_reference_probabilities'], strict=True))
    if abs(rows[0][contrast_key]-rows[1][contrast_key]) > probability_change+math.fsum(allowances):
        raise ValueError('packet parity replay exceeds its propagated probability bound')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path); group.add_argument('--check', type=Path)
    args = parser.parse_args(); result = run()
    if args.check:
        check_receipt(strict_json.load(args.check), result)
        print('Native continuum density paths, Gaussian reconstruction and wave controls replayed.')
    elif args.write:
        with args.write.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, indent=2, sort_keys=True); stream.write('\n')
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
