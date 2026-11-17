"""One native history gives jointly coupled continuum observation paths.

The time parameter is rescaled cumulative observation count. The field is a
long-observation Gaussian limit, not instantaneous matter or a wave amplitude.
The proof and access contract are in OCCUPATION_FIELD.md.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import sympy as sp

from . import occupation_field as occupation
from . import continuum_field as continuum
from . import local_continuum as spatial
from ..observer import native
from ..observer.scalar_response.experiment import canonical, compare_receipt
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*occupation.SOURCE_PATHS, *continuum.SOURCE_PATHS,
    *spatial.SOURCE_PATHS,
    'model/physical_space_time/observation_process.py'))))


def _integer(value, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError('an original integer in the admitted range is required')
    return value


def _positive(value):
    if (type(value) is not F or type(value.numerator) is not int
            or type(value.denominator) is not int or value.numerator <= 0
            or value.denominator <= 0 or math.gcd(value.numerator, value.denominator) != 1):
        raise ValueError('a positive canonical Fraction is required')
    return value


def clock_difference(level):
    """Exact covariance coefficient of the normalized same-prefix difference."""
    scale = 4**_integer(level)
    m, seams = 33, 630*scale
    probability = F(m, seams)
    return (1-probability)/(4*m*scale)


def path_bounds(level, *, horizon=F(1), count=None, degree=0):
    """All-level path bounds, requiring no allocation of a refined graph.

    A count supplies the exact finite-prefix reverse-martingale bound. Without
    it, the joint fixed-graph functional CLT has already been taken.
    """
    level, degree = _integer(level), _integer(degree)
    horizon = _positive(horizon)
    difference = clock_difference(level)
    result = {'joint_Gaussian_clock_sup_squared_upper': 8*horizon*difference,
        'Gaussian_harmonic_sup_tail_squared_upper': 4*horizon*continuum.harmonic_tail(degree),
        'Gaussian_fourth_increment_coefficient_upper': 3*continuum.harmonic_tail(0)**2,
        'clock_difference_covariance_coefficient': difference}
    if count is not None:
        _integer(count, 1)
        endpoint = (count*horizon.numerator+horizon.denominator-1)//horizon.denominator
        result['finite_prefix_clock_sup_squared_upper'] = 18*F(endpoint, count)*difference
    return result


@lru_cache(maxsize=1)
def marked_source():
    """Every native endpoint branch, with stationary marking verified exactly."""
    states, remote, incident, z, checked = occupation.native_counts()
    seams, _ = occupation.passage.recovered_carrier()
    index = {state: i for i, state in enumerate(states)}
    targets = np.empty((66, 60), dtype=np.int64)
    marks = np.array([int(0 in edge) for edge in seams for _ in (0, 1)], dtype=np.int64)
    for i, state in enumerate(states):
        for e, edge in enumerate(seams):
            for coin in (0, 1):
                final, (_, current) = native.native(state, edge, coin)
                decoded = tuple(state[p]+current*(int(p == edge[1])-int(p == edge[0]))
                                for p in range(12))
                if decoded != final:
                    raise ArithmeticError('the complete current tape lost field readback')
                targets[i, 2*e+coin] = index[final]
        # The instrument on each seam is exactly half identity plus half swap.
        for e, (a, b) in enumerate(seams):
            swapped = list(state); swapped[a], swapped[b] = swapped[b], swapped[a]
            if sorted(targets[i, 2*e:2*e+2]) != sorted((i, index[tuple(swapped)])):
                raise ArithmeticError('a native seam no longer preserves the uniform sector')
    uniform = sp.ones(66, 1)
    signed = 5*incident-remote
    if signed*uniform != sp.zeros(66, 1) or sp.ones(1, 66)*signed != sp.zeros(1, 66):
        raise ArithmeticError('a stationary marked output has nonzero conditional mean')
    score = np.array([[6*x-1 for x in state] for state in states], dtype=np.int64)
    if sp.Matrix(score.T) != 6*z:
        raise ArithmeticError('the integer current-decoded score changed')
    for array in (targets, marks, score):
        array.flags.writeable = False
    return states, targets, marks, score, checked


def finite_covariances(count):
    """Exact covariances of attempt sum A and clock-difference sum R.

    R=sum_j (I_j/p-1) eta_j uses post-repair eta_j; p=1/6.
    These formulas are valid for every count; this exact executable is bounded.
    """
    _integer(count, 1, 16)
    source = occupation.field('attempt')
    p = source['kernel']
    _, _, incident, _, _ = occupation.native_counts()
    z = sp.Matrix(marked_source()[3].T)/6
    # Density-mark matrix B/p-P is reconstructed from the complete state map.
    local = source['laplacian']-occupation.field('incident')['remote_laplacian']
    b = sp.eye(12)/6-local/60
    r = 6*b-p
    if z*(incident/10-(occupation.native_counts()[1]+incident)/60) != r*z:
        raise ArithmeticError('marked density closure disagrees with full native configurations')
    cross, power = sp.zeros(12), sp.eye(12)
    for lag in range(1, count):
        cross += (count-lag)*power*r
        power = power*p
    return tuple(sp.ImmutableMatrix(x) for x in (
        count*occupation.covariance_at_count('attempt', count),
        count*5*occupation.CHI*occupation.P0,
        occupation.CHI*cross))


def direct_two_steps(*, sticky=False, stationary=True):
    """Independent full endpoint-path moments, with no density-transfer use.

    sticky=True reuses the first seam at the second attempt. Each marginal
    remains stationary, but the IID-clock theorem need not survive.
    """
    if type(sticky) is not bool or type(stationary) is not bool:
        raise ValueError('original Boolean countercontrol switches are required')
    states, targets, marks, score, _ = marked_source()
    weights = 6*marks-1
    aa = np.zeros((12, 12), dtype=np.int64)
    rr = np.zeros_like(aa); ar = np.zeros_like(aa); lag = np.zeros_like(aa)
    mean = np.zeros(12, dtype=np.int64)
    terminal = np.zeros(66, dtype=np.int64)
    origins = range(66) if stationary else (0,)
    branches = 0
    for origin in origins:
        first = targets[origin]
        if sticky:
            columns = np.array([[2*(i//2), 2*(i//2)+1] for i in range(60)])
            final = targets[first[:, None], columns]
            second_weight = weights[columns]
        else:
            final = targets[first]
            second_weight = np.broadcast_to(weights, final.shape)
        one = np.broadcast_to(score[first, None, :], (*final.shape, 12))
        two = score[final]
        r1 = weights[:, None, None]*one
        r2 = second_weight[:, :, None]*two
        attempted = (one+two).reshape(-1, 12)
        difference = (r1+r2).reshape(-1, 12)
        aa += attempted.T@attempted
        rr += difference.T@difference
        ar += attempted.T@difference
        lag += r2.reshape(-1, 12).T@r1.reshape(-1, 12)
        mean += difference.sum(axis=0)
        terminal += np.bincount(final.ravel(), minlength=66)
        branches += final.size
    denominator = 36*branches
    return {'branches': branches, 'attempt_second_moment': sp.ImmutableMatrix(aa)/denominator,
        'difference_second_moment': sp.ImmutableMatrix(rr)/denominator,
        'cross_second_moment': sp.ImmutableMatrix(ar)/denominator,
        'difference_lag_one': sp.ImmutableMatrix(lag)/denominator,
        'difference_mean': sp.ImmutableMatrix(mean)/(6*branches),
        'terminal_counts': tuple(map(int, terminal)), 'initial_state': states[0]}


def joint_clt():
    """Cross-covariance derives the incident covariance on the same history."""
    attempt = occupation.field('attempt')
    incident = occupation.field('incident')
    local = attempt['laplacian']-incident['remote_laplacian']
    r = attempt['laplacian']/60-local/10
    cross = 60*occupation.CHI*attempt['green']*r
    delta = 5*occupation.CHI*occupation.P0
    if attempt['covariance']+cross+cross.T+delta != 6*incident['covariance']:
        raise ArithmeticError('joint clock covariance does not recover the exact stopped law')
    return sp.ImmutableMatrix(cross), sp.ImmutableMatrix(delta)


def martingale_certificate():
    """Exact forward innovations for the JOINT reward functional CLT.

    This is independent of using reverse orthogonality to obtain the clock
    difference covariance. It checks each conditional mean on all 66 states.
    """
    _, targets, marks, score, _ = marked_source()
    attempt = occupation.field('attempt')
    local = attempt['laplacian']-occupation.field('incident')['remote_laplacian']
    r = attempt['laplacian']/60-local/10
    h = 60*attempt['green']
    qa, qr = h-occupation.P0, r*h
    denominator = int(sp.ilcm(*(x.q for matrix in (qa, qr) for x in matrix)))
    qa = np.array(qa*denominator, dtype=object)
    qr = np.array(qr*denominator, dtype=object)
    after = score[targets].astype(object)
    jump = after-score[:, None, :]
    a = denominator*after+jump@qa.T
    z = denominator*(6*marks-1)[None, :, None]*after+jump@qr.T
    if np.any(a.sum(axis=1)) or np.any(z.sum(axis=1)):
        raise ArithmeticError('a joint reward innovation is not conditionally centered')
    a, z = a.reshape(-1, 12), z.reshape(-1, 12)
    normalizer = 36*denominator**2*len(a)
    blocks = tuple(sp.ImmutableMatrix(x)/normalizer for x in (a.T@a, z.T@z, a.T@z))
    cross, delta = joint_clt()
    if blocks != (attempt['covariance'], delta, cross):
        raise ArithmeticError('native joint martingale covariance disagrees with the Green formulas')
    return {'native_innovations': len(a), 'conditional_states': len(score),
            'all_conditional_means_zero': True, 'joint_asymptotic_covariance_checked': True}


def cylinder_action(times=(F(1, 2), F(1), F(2))):
    """Gaussian increment action AFTER the native functional CLT."""
    if type(times) is not tuple or not times or len(times) > 8:
        raise ValueError('one through eight increasing exact observation times required')
    for time in times:
        _positive(time)
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError('observation times must strictly increase')
    size = len(times)
    time_cov = sp.Matrix(size, size, lambda i, j: min(times[i], times[j]))
    increments = tuple(b-a for a, b in zip((F(), *times[:-1]), times, strict=True))
    difference = sp.eye(size)
    for i in range(1, size):
        difference[i, i-1] = -1
    time_precision = difference.T*sp.diag(*(1/x for x in increments))*difference
    spatial = occupation.field('attempt')['covariance'][:11, :11]/(120*occupation.CHI)
    precision = spatial.inv()
    covariance = sp.kronecker_product(time_cov, spatial)
    inverse = sp.kronecker_product(time_precision, precision)
    if covariance*inverse != sp.eye(11*size):
        raise ArithmeticError('the native Gaussian path action is not its joint precision')
    return {'times': times, 'increments': increments,
        'time_covariance': time_cov, 'time_precision': time_precision,
        'spatial_covariance': spatial, 'spatial_precision': precision,
        'dimension': 11*size}


def _rows(matrix):
    return [[str(x) for x in row] for row in matrix.tolist()]


def _radial_bridge_step(eigenvalues, start, stop, length):
    """Exact scalar formulas for a Dirichlet bridge; floats evaluate the formulas."""
    if not 0 <= start < stop < length:
        raise ValueError('strictly ordered interior radial coordinates required')
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    if np.any(eigenvalues < 0) or not np.all(np.isfinite(eigenvalues)):
        raise ValueError('nonnegative finite generator eigenvalues required')
    root = np.sqrt(eigenvalues)
    positive = root > 0
    transfer = np.full(len(root), (length-stop)/(length-start))
    innovation = np.full(len(root), (stop-start)*(length-stop)/(length-start))
    h = root[positive]
    denominator = np.sinh(h*(length-start))
    transfer[positive] = np.sinh(h*(length-stop))/denominator
    innovation[positive] = (
        np.sinh(h*(stop-start))*np.sinh(h*(length-stop))/(h*denominator))
    return transfer, innovation


def _radial_dirichlet_covariance(eigenvalues, first, second, length):
    """Green kernel of -d_r^2+lambda, including the Brownian-bridge zero mode."""
    lower, upper = sorted((first, second))
    if not 0 < lower <= upper < length:
        raise ValueError('interior radial coordinates required')
    root = np.sqrt(eigenvalues)
    positive = root > 0
    result = np.full(len(root), lower*(length-upper)/length)
    h = root[positive]
    result[positive] = np.sinh(h*lower)*np.sinh(h*(length-upper))/(h*np.sinh(h*length))
    return result


def radial_product_bounds(energy=F(12), sobolev_order=F(2)):
    """All-level product spectral and radial-path bounds, not physical units."""
    energy, sobolev_order = _positive(energy), _positive(sobolev_order)
    if energy < 3 or sobolev_order <= 1:
        raise ValueError('energy at least three and Sobolev order greater than one required')
    # On (0,pi) with Dirichlet ends, radial eigenvalues are k^2, k>=1.
    lower_k = math.isqrt((energy/2).numerator//(energy/2).denominator)
    lower_n = (energy/3).numerator//(energy/3).denominator
    upper_k = math.isqrt(energy.numerator//energy.denominator)
    upper_n = int(49920*energy)
    for index in (1, lower_n, upper_n):
        enclosure = spatial.continuum_eigenvalue_enclosure(index)
        if enclosure != {'lower': F(index, 49920), 'upper': F(3*index, 2)}:
            raise ArithmeticError('source spectral counting premise changed')
    harmonic_constant = 1+6*(1+1/(2*sobolev_order-2))
    return {
        'energy': str(energy), 'radial_domain': 'declared interval (0,pi), Dirichlet ends',
        'certified_positive_angular_mode_count_lower': lower_k*lower_n,
        'certified_including_zero_angular_mode_count_upper': upper_k*(upper_n+1),
        'spectral_count_power': '3/2',
        'all_energy_count_enclosure': 'Lambda^(3/2)/(12*sqrt(2)) <= N(Lambda) <= 49921*Lambda^(3/2), Lambda>=3',
        'angular_sobolev_order': str(sobolev_order),
        'radial_increment_second_moment_per_unit_distance_upper': str(harmonic_constant),
        'radial_increment_fourth_moment_per_squared_distance_upper': str(3*harmonic_constant**2),
        'continuum_radial_path_topology': f'C([0,L],H^-{sobolev_order}(S2))',
        'product_comparison_distance_squared': '|r-s|^2+2*d_screen^2',
        'product_wave_speed_in_this_comparison_distance': '1',
        'interval_boundary_length_or_metric_selected_by_repair': False,
    }


def product_mode_dictionary(native_eigenvalue, radial_mode=1, radial_points=3,
                            length=F(2), angular_scale=F(1)):
    """Exact radial-grid/native-mode dictionary; no clock or mode is selected."""
    radial_points = _integer(radial_points, 1, 64)
    radial_mode = _integer(radial_mode, 1, radial_points)
    length, angular_scale = _positive(length), _positive(angular_scale)
    if type(native_eigenvalue) not in (int, F) and not isinstance(native_eigenvalue, sp.Expr):
        raise ValueError('an exact nonnegative native eigenvalue is required')
    if type(native_eigenvalue) is F and (
            native_eigenvalue.denominator <= 0
            or math.gcd(native_eigenvalue.numerator, native_eigenvalue.denominator) != 1):
        raise ValueError('a canonical exact native eigenvalue is required')
    eigenvalue = sp.sympify(native_eigenvalue)
    if (eigenvalue.free_symbols or eigenvalue.atoms(sp.Float) or eigenvalue.is_real is not True
            or eigenvalue.is_nonnegative is not True or eigenvalue.is_finite is not True):
        raise ValueError('an exact nonnegative native eigenvalue is required')
    spacing = sp.Rational(length)/(radial_points+1)
    radial = sp.simplify(4*sp.sin(sp.pi*radial_mode/(2*(radial_points+1)))**2/spacing**2)
    angular = sp.Rational(angular_scale)*eigenvalue
    return {
        'native_eigenvalue': eigenvalue, 'angular_scale': sp.Rational(angular_scale),
        'scaled_native_eigenvalue': angular, 'radial_eigenvalue': radial,
        'product_eigenvalue': angular+radial,
        'continuum_radial_eigenvalue': (sp.pi*radial_mode/sp.Rational(length))**2,
        'radial_spacing': spacing, 'radial_mode': radial_mode,
        'radial_points': radial_points, 'interval_length': sp.Rational(length),
        'radial_boundary_condition': 'Dirichlet',
        'native_eigenvector_membership_requires_a_separate_certificate': True,
    }


def product_thermal_bounds(inverse_temperature=F(1), length=F(2),
                           radial_points=3, radial_band=1, horizon=F(1)):
    """Rational grid, thermal-state and H^-3 path bounds for the declared product."""
    beta, length, horizon = map(_positive, (inverse_temperature, length, horizon))
    points = _integer(radial_points, 1, 64)
    band = _integer(radial_band, 1, points)
    spacing = length/(points+1)
    # sin x >= x-x^3/6 and pi < 22/7 bound the squared-frequency defect.
    radial_defect = spacing**2*(F(22, 7)*band/length)**4/12
    minimum_frequency = 2/length
    covariance_norm = length/4+length**2/(4*beta)
    derivative_bound = length**3/32+length**4/(16*beta)
    response_derivative_bound = horizon**2/(4*minimum_frequency)
    ordered_derivative_bound = (derivative_bound
        +covariance_norm*horizon/(2*minimum_frequency)+response_derivative_bound/2)
    # N_B(Lambda) <= C Lambda^(3/2), using pi>3 and Lambda>=(pi/L)^2.
    count_constant = length*(49920+length**2/9)/3
    velocity_trace = 36*(1+F(11, 7)/length+1/beta)
    return {
        'inverse_temperature': beta, 'interval_length': length,
        'radial_points': points, 'radial_band': band, 'parameter_horizon': horizon,
        'radial_squared_frequency_error_upper': radial_defect,
        'frequency_gap_lower_including_angular_constant': minimum_frequency,
        'thermal_equal_time_covariance_norm_upper': covariance_norm,
        'thermal_covariance_grid_error_upper': radial_defect*derivative_bound,
        'sine_response_grid_error_upper': radial_defect*response_derivative_bound,
        'ordered_thermal_covariance_grid_error_upper': radial_defect*ordered_derivative_bound,
        'one_particle_heat_trace_upper': 6*count_constant/beta**3,
        'normal_ordered_Fock_log_partition_upper': 8*count_constant/beta**3,
        'fixed_product_H_minus_3_initial_second_moment_upper': 32*covariance_norm,
        'fixed_product_H_minus_3_velocity_second_moment_upper': velocity_trace,
        'fixed_product_H_minus_3_increment_fourth_moment_per_delta_fourth_upper': 3*velocity_trace**2,
        'full_angular_constant_included_without_a_zero_frequency': True,
        'radial_grid_errors_are_not_a_spatial_subsequence_rate': True,
        'classical_covariance_processing_establishes_CCR': False,
        'thermal_state_or_radial_domain_selected_by_repair': False,
    }


def native_block_covariance_bound(level, block_length, blocks, filter_norm_squared=F(1)):
    """A covariance error, not a finite-count Gaussian distribution guarantee.

    The internal gap is at least two, the cross-seam Laplacian norm at most
    two, and the existing Schur Green bound gives gap(L)>=1/(4096*4^l+1).
    """
    level, count, blocks = _integer(level), _integer(block_length, 1), _integer(blocks, 1)
    gain = _positive(filter_norm_squared)
    scale, seams = 4**level, 630*4**level
    contraction = 1-F(1, 2*seams*(4096*scale+1))
    cross = contraction/(4*seams*scale*count*(1-contraction)**2)
    return {'level': level, 'native_attempts_per_block': count, 'blocks': blocks,
        'filter_norm_squared_upper': gain,
        'native_density_contraction_upper': contraction,
        'one_block_covariance_defect_upper': 2*cross,
        'adjacent_block_cross_covariance_upper': cross,
        'filtered_joint_covariance_defect_upper': 2*blocks*cross*gain,
        'bound_applies_before_Gaussian_limit': True,
        'finite_count_Gaussian_total_variation_bound': False}


def product_thermal_control(inverse_temperature=F(1)):
    """Execute the declared radial-grid/thermal filter on the native records.

    The classical symmetric covariance is measured here. Noncommuting ordered
    readouts and the active product clock are certified by observer_time_fields.
    """
    beta = _positive(inverse_temperature)
    if not F(1, 16) <= beta <= 4:
        raise ValueError('the numerical thermal control admits inverse temperature in [1/16,4]')
    points, cells, length, spacing = 3, 20, 2., .5
    _, angular, _ = spatial._green_orthonormal_bundle(0)
    values, vectors = np.linalg.eigh(angular)
    values[0] = 0.
    grid = np.arange(1, points+1)
    sine = np.sqrt(2/(points+1))*np.sin(np.pi*np.outer(grid, grid)/(points+1))
    radial = np.array([float(product_mode_dictionary(
        0, k, points)['radial_eigenvalue']) for k in grid.tolist()])
    transform = np.kron(sine, vectors)
    spectrum = (radial[:, None]+values[None, :]).reshape(-1)
    frequency = np.sqrt(spectrum)
    thermal = 1/(2*frequency*np.tanh(float(beta)*frequency/2))
    product = (np.diag(np.full(points, 2.))+np.diag(np.full(points-1, -1.), 1)
               +np.diag(np.full(points-1, -1.), -1))/spacing**2
    product = np.kron(product, np.eye(cells))+np.kron(np.eye(points), angular)
    product_error = float(np.linalg.norm(
        product-(transform*spectrum)@transform.T, 2))
    positions = spacing*grid
    bridge = spacing*np.block([[(vectors*_radial_dirichlet_covariance(
        values, r, s, length))@vectors.T for s in positions] for r in positions])
    # The precision of sampled continuous Dirichlet bridges is tridiagonal.
    # Its DST eigenvalues are not rho_h+lambda at a finite radial mesh.
    bridge_spectrum = np.empty((points, cells))
    bridge_spectrum[:, 0] = 1/radial
    for j, lam in enumerate(values[1:], 1):
        root = math.sqrt(lam)
        bridge_spectrum[:, j] = spacing*np.sinh(spacing*root)/(
            2*root*(np.cosh(spacing*root)-np.cos(np.pi*grid/(points+1))))
    bridge_spectrum = bridge_spectrum.reshape(-1)
    bridge_error = float(np.linalg.norm(
        bridge-(transform*bridge_spectrum)@transform.T, 2))
    classical = np.linalg.inv(product)
    target = (transform*thermal)@transform.T
    gain = np.sqrt(thermal/bridge_spectrum)
    transfer = (transform*gain)@transform.T
    covariance_error = float(np.linalg.norm(transfer@bridge@transfer.T-target, 2))
    raw_bridge_discrepancy = float(np.linalg.norm(bridge-classical, 2))
    classical_thermal_discrepancy = float(np.linalg.norm(classical-target, 2))
    execution = radial_record_gluing()['finite_native_record_execution']
    record = math.sqrt(spacing)*np.asarray(execution['all_radial_slice_readouts']).reshape(-1)
    output = transfer@record
    # Ordered thermal covariance has a state-independent imaginary part.
    # The real record filter supplies only its symmetric, commuting part.
    time = .375
    ordered = thermal*np.cos(time*frequency)-.5j*np.sin(time*frequency)/frequency
    ordered_minus = thermal*np.cos(time*frequency)+.5j*np.sin(time*frequency)/frequency
    commutator_error = float(np.max(abs(
        ordered-ordered_minus+1j*np.sin(time*frequency)/frequency)))
    argument = (time-1j*float(beta))*frequency
    bose_occupation = 1/np.expm1(float(beta)*frequency)
    # Exponential weights avoid subtracting two O(exp(beta*omega)) terms
    # in the complex cosine/sine representation at the KMS boundary.
    shifted = ((1+bose_occupation)*np.exp(-1j*argument)
               +bose_occupation*np.exp(1j*argument))/(2*frequency)
    kms_error = float(np.max(abs(
        shifted-ordered_minus)))
    omitted = np.kron(np.eye(points), np.eye(cells)-np.ones((cells, cells))/cells)
    omitted_locality_defect = abs(float(omitted[0, 1]))
    # Actual full-port modes and cell Schur modes are not the same finite
    # operator. A shifted resolvent checks their common-limit source mapping.
    graph = spatial.source.source.federation(0)
    full = np.zeros((graph.ports, graph.ports))
    for a, b in graph.seams:
        full[a, a] += 1; full[b, b] += 1
        full[a, b] -= 1; full[b, a] -= 1
    lift = np.kron(np.eye(cells), np.ones((12, 1))/math.sqrt(12))
    resolvent_errors = []
    for rho in radial:
        shift = 1+rho
        compressed = lift.T@np.linalg.solve(full+shift*np.eye(graph.ports), lift)
        resolvent_errors.append(float(np.linalg.norm(
            np.linalg.inv(angular+shift*np.eye(cells))-compressed, 2)))
    if (max(product_error, bridge_error, covariance_error, commutator_error, kms_error)>1e-10
            or min(raw_bridge_discrepancy, classical_thermal_discrepancy)<1e-4
            or max(resolvent_errors) >= .125 or abs(omitted_locality_defect-.05)>1e-14):
        raise ArithmeticError('product thermal covariance or source-operator control failed')
    return {
        'inverse_temperature': str(beta), 'radial_points': points, 'angular_cells': cells,
        'interval_length': '2', 'radial_boundary_condition': 'Dirichlet',
        'radial_eigenvalues': [str(product_mode_dictionary(
            0, k, points)['radial_eigenvalue']) for k in grid.tolist()],
        'cell_generator_is_native_Schur_not_full_port_generator': True,
        'direct_product_matrix_diagonalization_error': product_error,
        'sampled_bridge_DST_covariance_error': bridge_error,
        'sampled_bridge_vs_grid_inverse_discrepancy': raw_bridge_discrepancy,
        'classical_inverse_vs_quantum_thermal_covariance_discrepancy': classical_thermal_discrepancy,
        'derived_record_filter_covariance_error': covariance_error,
        'record_filter_squared_operator_norm': float(np.max(gain**2)),
        'thermal_covariance_is_coth_beta_sqrt_B_over_two_divided_by_two_sqrt_B': True,
        'ordered_covariance_commutator_identity_error': commutator_error,
        'ordered_covariance_KMS_identity_error': kms_error,
        'smallest_mode_uncertainty_product': float(np.min(spectrum*thermal**2)),
        'omitted_angular_constant_disjoint_readout_response_derivative': -omitted_locality_defect,
        'full_product_response_derivative_is_identity': True,
        'source_shifted_product_resolvent_errors': resolvent_errors,
        'source_shifted_product_resolvent_error_upper': '1/8',
        'all_thermal_radial_readouts_from_complete_96_event_tape': output.reshape(points, cells).tolist(),
        'complete_signed_currents_retained_in_radial_record_innovation_controls': True,
        'finite_tape_certifies_Gaussian_thermal_ensemble': False,
        'classical_record_filter_establishes_quantum_CCR': False,
        'state_temperature_radial_domain_and_filter_are_observer_contracts': True,
        'physical_clock_radial_distance_or_interacting_QFT_identified': False}


def radial_record_gluing(level=0, slices=3, length=F(2)):
    """Native block innovations synthesize a nondegenerate local product field.

    Native count tends to infinity before this Gaussian encoding is applied.
    Radial labels, boundary conditions and filters are observer choices.
    """
    level, slices = _integer(level, 0, 1), _integer(slices, 3, 6)
    length = _positive(length)
    if not F(1, 16) <= length <= 4:
        raise ValueError('the finite numerical control admits length between 1/16 and four')
    green, generator, projector = spatial._green_orthonormal_bundle(level)
    cells, scale = len(green), 4**level
    correction = F(1, 2520*scale*scale)
    native_covariance = green-float(correction)*projector
    constant = np.ones((cells, cells))/cells
    # A second independent native block supplies one normalized nonzero
    # smearing. Encoding it in the constant angular direction does not make
    # the conserved source load fluctuate.
    contrast = np.zeros(cells); contrast[0], contrast[1] = 1/math.sqrt(2), -1/math.sqrt(2)
    scalar_variance = float(contrast@native_covariance@contrast)
    eigenvalues, vectors = np.linalg.eigh(generator)
    if abs(eigenvalues[0]) > 1e-12 or eigenvalues[1] <= 0 or scalar_variance <= 0:
        raise ArithmeticError('native covariance lost the radial innovation source')
    eigenvalues[0] = 0
    source_variances = np.concatenate(([1.], 1/eigenvalues[1:]-float(correction)))
    if np.min(source_variances) <= 0:
        raise ArithmeticError('native finite-clock covariance is not positive')
    source_covariance = native_covariance+constant
    times = [float(length*i/(slices+1)) for i in range(1, slices+1)]
    transfers, noises, filters = [], [], []
    for start, stop in zip([0., *times[:-1]], times, strict=True):
        transfer, variance = _radial_bridge_step(eigenvalues, start, stop, float(length))
        transfers.append((vectors*transfer)@vectors.T)
        noises.append((vectors*variance)@vectors.T)
        filters.append((vectors*np.sqrt(variance/source_variances))@vectors.T)
    filter_error = max(float(np.linalg.norm(f@source_covariance@f.T-q, 2))
                       for f, q in zip(filters, noises, strict=True))
    covariance = np.empty((slices, slices, cells, cells))
    for i in range(slices):
        covariance[i, i] = noises[i] if i == 0 else (
            transfers[i]@covariance[i-1, i-1]@transfers[i].T+noises[i])
        for j in range(i):
            covariance[i, j] = transfers[i]@covariance[i-1, j]
            covariance[j, i] = covariance[i, j].T
    expected = [[(vectors*_radial_dirichlet_covariance(
        eigenvalues, r, s, float(length)))@vectors.T for s in times] for r in times]
    joint = np.block([[covariance[i, j] for j in range(slices)] for i in range(slices)])
    direct = np.block(expected)
    covariance_error = float(np.linalg.norm(joint-direct, 2))
    difference = np.eye(slices*cells)
    innovation_precision = np.zeros_like(joint)
    for i in range(slices):
        block = slice(i*cells, (i+1)*cells)
        innovation_precision[block, block] = np.linalg.inv(noises[i])
        if i:
            difference[block, slice((i-1)*cells, i*cells)] = -transfers[i]
    precision = difference.T@innovation_precision@difference
    precision_error = float(np.linalg.norm(precision@joint-np.eye(len(joint)), 2))
    nonneighbor_precision = max(
        float(np.linalg.norm(precision[i*cells:(i+1)*cells, j*cells:(j+1)*cells], 2))
        for i in range(slices) for j in range(slices) if abs(i-j)>1)
    conditional = covariance[0, 2]-covariance[0, 1]@np.linalg.solve(
        covariance[1, 1], covariance[1, 2])
    conditional_error = float(np.linalg.norm(conditional, 2))
    # Gluing an inserted radial layer must leave the original joint law unchanged.
    ra, qa = _radial_bridge_step(eigenvalues, 0., float(length)/5, float(length))
    rb, qb = _radial_bridge_step(eigenvalues, float(length)/5, 3*float(length)/5, float(length))
    rc, qc = _radial_bridge_step(eigenvalues, 0., 3*float(length)/5, float(length))
    gluing_error = float(max(np.max(abs(rc-rb*ra)), np.max(abs(qc-(qb+rb*rb*qa)))))
    # Reusing one native block correlates innovations and changes the law.
    reused_maps = []
    for i in range(slices):
        reused_maps.append(filters[i] if i == 0 else transfers[i]@reused_maps[-1]+filters[i])
    reused = np.block([[f@source_covariance@g.T for g in reused_maps] for f in reused_maps])
    reused_error = float(np.linalg.norm(reused-joint, 2))
    # Execute a complete finite tape to check the actual record-to-field map.
    # This particular tape has positive probability under the admitted IID law;
    # it is a decoder control, not a finite-count Gaussian sample certificate.
    graph = spatial.source.source.federation(level)
    occupied = tuple(range(0, graph.ports, 3))
    decoded = np.zeros(graph.ports, dtype=int); decoded[list(occupied)] = 1
    block_length, block_sums, currents = 16, [], []
    for block in range(2*slices):
        total = np.zeros(cells, dtype=int)
        for offset in range(block_length):
            event = block*block_length+offset
            edge = graph.seams[(37*event+11) % len(graph.seams)]
            occupied, current = spatial.source.source.binary_step(occupied, edge, (event//2) % 2)
            decoded[edge[0]] -= current; decoded[edge[1]] += current
            if tuple(np.flatnonzero(decoded)) != occupied or np.any((decoded < 0) | (decoded > 1)):
                raise ArithmeticError('radial native current tape lost complete field readback')
            currents.append(current)
            total += decoded.reshape(cells, 12).sum(axis=1)-4
        block_sums.append(total)
    q = graph.ports//3
    chi = F(q*(graph.ports-q), graph.ports*(graph.ports-1))
    normalization = math.sqrt(12*4*float(chi)*len(graph.seams)*scale*block_length)
    fields = np.asarray(block_sums)/normalization
    finite_readouts, current_field = [], np.zeros(cells)
    unit_constant = np.ones(cells)/math.sqrt(cells)
    for i in range(slices):
        innovation = fields[2*i]+unit_constant*float(
            contrast@fields[2*i+1])/math.sqrt(scalar_variance)
        current_field = transfers[i]@current_field+filters[i]@innovation
        finite_readouts.append(current_field.tolist())
    if (max(filter_error, covariance_error, precision_error, conditional_error, gluing_error)>1e-10
            or nonneighbor_precision != 0 or reused_error < 0.01):
        raise ArithmeticError('native radial innovation, gluing or precision certificate failed')
    return {
        'level': level, 'angular_cells': cells, 'interior_radial_slices': slices,
        'interval_length': str(length), 'radial_coordinates': [str(length*i/(slices+1)) for i in range(1, slices+1)],
        'native_long_block_covariance_correction': str(correction),
        'independent_CLT_blocks_of_one_history_required': 2*slices,
        'constant_mode_native_smearing_variance_diagnostic': scalar_variance,
        'native_filter_covariance_error': filter_error,
        'joint_cross_radius_covariance_error': covariance_error,
        'joint_precision_identity_error': precision_error,
        'nonneighbor_radial_precision_norm': nonneighbor_precision,
        'conditional_past_future_covariance_norm': conditional_error,
        'inserted_layer_gluing_error': gluing_error,
        'reused_block_covariance_discrepancy': reused_error,
        'smallest_conditional_innovation_eigenvalue': float(min(np.linalg.eigvalsh(q)[0] for q in noises)),
        'joint_covariance_rank': int(np.linalg.matrix_rank(joint)),
        'expected_joint_covariance_rank': cells*slices,
        'boundary_only_extension_rank_at_most': cells-1,
        'finite_native_record_execution': {
            'initial_occupied_ports': '0,3,6,...,V-3',
            'seam_index_at_event_j': '(37*j+11) modulo M',
            'coin_at_event_j': '(j//2) modulo 2',
            'native_attempts_per_block': block_length,
            'complete_signed_current_records_including_zero': currents,
            'all_radial_slice_readouts': finite_readouts,
            'decoded_field_checked_after_every_native_event': True,
            'finite_tape_is_a_Gaussian_ensemble_certificate': False,
        },
        'zero_angular_mode_is_independent_encoded_readout': True,
        'native_conserved_total_load_fluctuates': False,
        'radial_Gaussian_markov_field_synthesized_from_native_blocks': True,
        'Gaussian_independence_holds_at_finite_native_count': False,
        'radial_domain_filters_and_boundary_conditions_are_declared': True,
        'physical_radial_transport_or_Lorentzian_clock_identified': False,
    }


def finite_radial_block_control():
    """Two adjacent native blocks are correlated before the functional CLT."""
    _, _, _, z, checked = occupation.native_counts()
    full = occupation.configuration_kernel('attempt')
    reduced = occupation.field('attempt')['kernel']
    direct = z*(full+2*full**2+full**3)*z.T/(66*240*occupation.CHI)
    expected = (reduced+2*reduced**2+reduced**3)*occupation.P0/240
    if direct != expected or sp.trace(direct) <= 0:
        raise ArithmeticError('complete native adjacent-block covariance control failed')
    return {
        'native_branch_binding': checked,
        'native_attempts_per_adjacent_block': 2,
        'normalized_adjacent_block_cross_covariance_trace': str(sp.trace(direct)),
        'complete_66_state_and_density_oracles_agree_exactly': True,
        'adjacent_blocks_independent_at_finite_count': False,
        'fixed_carrier_block_length_n_cross_covariance_norm_upper': '29/(4n)',
        'bound_is_uniform_over_source_refinement': False,
    }


def run():
    _, _, _, _, branches = marked_source()
    actual = direct_two_steps()
    predicted = finite_covariances(2)
    for key, expected in zip(('attempt_second_moment', 'difference_second_moment',
                              'cross_second_moment'), predicted, strict=True):
        if actual[key] != expected:
            raise ArithmeticError('the complete native two-step oracle differs from the theorem')
    if actual['difference_lag_one'] != sp.zeros(12) or actual['difference_mean'] != sp.zeros(12, 1):
        raise ArithmeticError('the reverse martingale covariance check failed')
    sticky, biased = direct_two_steps(sticky=True), direct_two_steps(stationary=False)
    if (sticky['difference_lag_one'] == sp.zeros(12)
            or len(set(sticky['terminal_counts'])) != 1
            or biased['difference_mean'] == sp.zeros(12, 1)):
        raise ArithmeticError('a hypothesis countercontrol stopped distinguishing the laws')
    cross, delta = joint_clt()
    action = cylinder_action()
    return {'schema': 'oph.native_joint_observation_process.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'native_exact_controls': {'one_step_branches': branches,
            'two_step_complete_branches': actual['branches'],
            'post_repair_mark_mean_zero_given_output_checked': True,
            'same_history_clock_difference_lag_one_zero': True,
            'difference_covariance_per_attempt': _rows(delta),
            'attempt_difference_asymptotic_cross_covariance': _rows(cross),
            'joint_covariance_recovers_exact_incident_clock': True},
        'joint_functional_CLT_innovations': martingale_certificate(),
        'countercontrols': {
            'sticky_seam_preserves_uniform_one_time_state': True,
            'sticky_seam_two_step_branches': sticky['branches'],
            'sticky_seam_difference_lag_one_trace': str(sp.trace(sticky['difference_lag_one'])),
            'nonstationary_preparation_difference_mean': _rows(biased['difference_mean']),
            'finite_native_adjacent_occupation_covariance': str(occupation.CHI*(
                occupation.field('attempt')['kernel']*occupation.P0)[0, 0]),
            'uncentered_constant_readout_random_time_change_noise_variance_per_attempt': '5'},
        'Gaussian_cylinder_action': {'dimension': action['dimension'],
            'times': list(map(str, action['times'])),
            'time_covariance': _rows(action['time_covariance']),
            'time_precision': _rows(action['time_precision']),
            'joint_precision_equals_time_precision_tensor_native_spatial_precision': True,
            'applies_after_functional_CLT_not_at_finite_native_count': True},
        'radial_record_innovation_controls': [radial_record_gluing(level) for level in (0, 1)],
        'radial_product_field_bounds': radial_product_bounds(),
        'thermal_product_covariance_control': product_thermal_control(),
        'thermal_product_continuum_bounds': {key: str(value) if type(value) is F else value
            for key, value in product_thermal_bounds().items()},
        'native_block_covariance_budget': {key: str(value) if type(value) is F else value
            for key, value in native_block_covariance_bound(0, 1000000, 6).items()},
        'finite_radial_block_countercontrol': finite_radial_block_control(),
        'all_level_path_bounds': {str(level): {k: str(v) for k, v in path_bounds(
            level, count=100, degree=20).items()} for level in (0, 1, 3, 10, 30)},
        'theorem': {'joint_limit_in_C_of_time_into_H_minus_2_of_S2': True,
            'time_covariance': 'min(t,u) times the native spatial subsequence covariance',
            'normalized_native_clocks_are_jointly_coupled_to_the_same_path': True,
            'finite_count_same_prefix_clock_error_bound': True,
            'incident_index_time_change_at_each_fixed_source_level': True,
            'limit_order': 'native observation count to infinity, then common spatial subsequence',
            'full_field_requires_initial_load_and_complete_current_read_contract': True,
            'instantaneous_density_process_identified_with_Brownian_path': False,
            'uniform_joint_count_refinement_CLT_error': False,
            'unique_spatial_limit_or_physical_radial_coordinate': False,
            'Lorentzian_propagator_interacting_QFT_or_modular_flow_identity': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--check', type=Path)
    modes.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        compare_receipt(strict_json.load(args.check), result)
        print('Joint native clock coupling, continuum observation paths and Gaussian action verified.')
    elif args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(canonical(result))
    else:
        print(canonical(result), end='')


if __name__ == '__main__':
    main()
