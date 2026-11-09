"""Exact binomial regulator and weighted local interaction certificates.

The smooth action is a declared alternative to discontinuous overlap sign.
Its own polynomial approximants share state, electric term, action and
clock with that smooth target. No zero-smoothing or continuum limit is used.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

import sympy as s

from model.particle_local_gauge_state import state
from model.particle_local_history_family.history import ceil_log2, ceil_fraction, cutoff_plan
from model.particle_local_wilson_density import density

SOURCE_INPUTS = (
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_field_causality/PROOF.md',
    'model/particle_field_causality/bounds.py',
    'model/particle_local_gauge_state/state.py',
    'model/particle_local_gauge_state/PROOF.md',
    'model/particle_local_history_family/history.py',
    'model/particle_local_history_family/PROOF.md',
    'model/particle_local_wilson_density/density.py',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/interacting_extension/chiral_lattice/PROOF.md',
)


def settings(delta='1', spacing='1', plaquette='0'):
    d = state.rational(delta, 'positive smoothing parameter', True)
    a = state.rational(spacing, 'positive spacing', True)
    k = state.rational(plaquette, 'plaquette coefficient')
    D = 64+d*d
    q = F(64)/D
    u = 1-q
    weight = 1+u/8
    rate = q*weight**4
    if not 0 < rate <= 1-u/2 < 1:
        raise ArithmeticError('weighted local series does not contract')
    # sqrt(1+delta²)<=1+delta and8/sqrt(64+delta²)<=1.
    root = 60*(1+d)/a
    kinetic = root*weight**2*4096*720/(1-rate)**7
    magnetic = 48*k*weight**2
    return {
        'delta': d, 'spacing': a, 'plaquette_coefficient': k,
        'normalization_squared': 1+d*d,
        'series_denominator_squared': D,
        'binomial_operator_ratio': q, 'one_minus_ratio': u,
        'spatial_weight_base': weight, 'weighted_series_ratio': rate,
        'rooted_nth_term_norm_prefactor': root,
        'weighted_kinetic_interaction_norm_upper': kinetic,
        'weighted_full_interaction_norm_upper': kinetic+magnetic,
        'lattice_weight_sum': ((weight+1)/(weight-1))**3,
        'one_body_H_norm_rational_upper': (2+d)/a,
        'full_Fock_kinetic_density_rational_upper': 60*(2+d)/a,
        'global_Wilson_gap_required': False,
        'same_as_discontinuous_sign_action': False,
        'normalization_uses_measured_particle_mass': False,
    }


def binomial_polynomial(kernel, degree, delta='1', normalized=True):
    """The declared smooth action's exact polynomial, not the old Legendre q_m."""
    m = state.integer(degree, 'binomial degree')
    d = state.rational(delta, 'positive smoothing parameter', True)
    h = s.ImmutableMatrix(kernel)
    if not h.rows or h.rows != h.cols or h != h.H:
        raise ValueError('a nonempty Hermitian source kernel is required')
    D = s.Rational(64)+s.Rational(d.numerator, d.denominator)**2
    z = (64*s.eye(h.rows)-h*h)/D
    power = s.eye(h.rows)
    value = s.zeros(h.rows)
    for n in range(m+1):
        value += s.Rational(comb(2*n, n), 4**n)*power
        power = power*z
    factor = s.sqrt(1+s.Rational(d.numerator, d.denominator)**2) if normalized else 1
    return s.ImmutableMatrix(factor*h*value/s.sqrt(D)).applyfunc(s.simplify)


def local_record_plan(config, duration='1', accuracy_bits=20,
                      observed_cells=1, intermediate_instruments=0):
    """Uniform complete local-history error for the smooth target and its truncation."""
    T = state.rational(duration, 'duration', True)
    bits = state.integer(accuracy_bits, 'accuracy bits', 1)
    cells = state.integer(observed_cells, 'observed cells', 1)
    k = state.integer(intermediate_instruments, 'intermediate instruments')
    W = config['weighted_full_interaction_norm_upper']
    A = config['weighted_kinetic_interaction_norm_upper']
    prefactor = (2*k+1)*2*cells*T*config['lattice_weight_sum']*A
    exponent0 = ceil_log2(prefactor)+ceil_fraction(4*W*T)
    u = config['one_minus_ratio']
    N = 1
    attempts = 0
    while True:
        exponent = exponent0+6*ceil_log2(F(N+1))-(u*N/2).__floor__()
        attempts += 1
        if exponent <= -bits:
            break
        N *= 2
    return {
        'duration': T, 'binomial_polynomial_degree': N-1,
        'maximum_coordinate_word_degree': 2*N-1,
        'first_omitted_series_term': N,
        'error_binary_exponent_upper': exponent,
        'complete_local_record_TV_upper': F(1, 2**bits),
        'accuracy_bits': bits, 'observed_cells': cells,
        'intermediate_local_instruments': k,
        'degree_search_steps': attempts,
        'shared_state_electric_plaquette_clock_and_full_Fock': True,
        'untouched_references_and_all_outcomes_retained': True,
        'depends_on_total_volume': False,
        'link_derivative_or_Wilson_gap_bound_required': False,
        'large_polynomial_materialized': False,
        'zero_delta_or_zero_spacing_limit_claimed': False,
    }


def finite_tail_bound(config, first):
    """Small exact control for the weighted-tail majorant, not huge default degrees."""
    n = state.integer(first, 'first omitted term')
    if n > 10000:
        raise ValueError('use the binary planner for large factored resource degrees')
    r = config['weighted_series_ratio']
    return config['weighted_kinetic_interaction_norm_upper']*(n+1)**6*r**n


def finite_volume_plan(config, approximation, accuracy_bits=20, observed_radius=0):
    """Open-box approximation of the same retained binomial polynomial.

    Use the smooth-series interaction bound, uniform in truncation degree,
    with the actual polynomial's finite boundary collar.
    """
    bits = state.integer(accuracy_bits, 'volume accuracy bits', 1)
    observed = state.integer(observed_radius, 'observed radius')
    d = approximation['maximum_coordinate_word_degree']
    collar = max(d, 2)
    R = observed+2*collar+1
    T = approximation['duration']
    k = approximation['intermediate_local_instruments']
    cells = approximation['observed_cells']
    if cells > (2*observed+1)**3:
        raise ValueError('observed cell count exceeds its declared observation cube')
    u = config['one_minus_ratio']
    J = config['full_Fock_kinetic_density_rational_upper']
    exp_bits = ceil_fraction(4*config['weighted_full_interaction_norm_upper']*T)
    attempts = 0
    while True:
        D = R-observed-2*collar
        C = 2*J*(2*(R+d)+1)**3+3*config['plaquette_coefficient']*(2*(R+2)+1)**3
        # log(lambda)>=u/(8+u)>=u/9.
        exponent = ceil_log2((2*k+1)*2*cells*C*T)+exp_bits-(u*D/9).__floor__()
        attempts += 1
        if exponent <= -bits:
            break
        R *= 2
    return {'box_radius': R, 'duration': T, 'intermediate_instruments': k,
            'observed_radius': observed, 'observed_cells': cells,
            'periodic_side_sufficient': 2*(R+2*collar+d)+1,
            'same_binomial_polynomial_degree': approximation['binomial_polynomial_degree'],
            'weighted_rate_uniform_in_polynomial_degree': True,
            'boundary_separation': D, 'boundary_perturbation_norm_upper': C,
            'error_binary_exponent_upper': exponent,
            'complete_record_TV_upper': F(1, 2**bits), 'radius_search_steps': attempts}


def free_symbol_squared(cosines, delta='1', spacing='1'):
    """Exact free energy squared for specified cos(k_i), keeping all branches."""
    if len(cosines) != 3:
        raise ValueError('three spatial cosine values are required')
    cs = [s.sympify(x) for x in cosines]
    if any(x.is_real is not True or x < -1 or x > 1 for x in cs):
        raise ValueError('exact real cosine data in[-1,1] are required')
    d = state.rational(delta, 'positive smoothing parameter', True)
    a = state.rational(spacing, 'positive spacing', True)
    B = sum(1-x for x in cs)-1
    sin_square = sum(1-x*x for x in cs)
    R2 = s.expand(B*B+sin_square)
    t = s.sqrt(1+s.Rational(d.numerator, d.denominator)**2)/s.sqrt(
        R2+s.Rational(d.numerator, d.denominator)**2)
    return s.simplify(((1+t*B)**2+t*t*sin_square)/s.Rational(a.numerator, a.denominator)**2)


@lru_cache(maxsize=1)
def exact_native_control():
    """Same native sixty-column zero-mode control, evaluated by the smooth action."""
    first = density.native_weak_link(0, s.Rational(1, 2), s.sqrt(3)/2)
    second = density.native_weak_link(1, s.Rational(1, 2), s.sqrt(3)/2)
    data = density.decomposition((first, second, s.eye(15)))
    h = data['kernel']; square = h*h; I = s.eye(60)
    p0 = (square-I)*(square-3*I)/3
    p1 = -square*(square-3*I)/2
    p3 = square*(square-I)/6
    ranks = [int(s.trace(p)) for p in (p0, p1, p3)]
    if ranks != [16, 28, 16]:
        raise ArithmeticError('the full native spectral control changed')
    # delta=1, normalized bysqrt2; exact spectral evaluation at0,1,3.
    target = s.ImmutableMatrix(h*(p1+p3/s.sqrt(2))).applyfunc(s.simplify)
    polynomial = binomial_polynomial(h, 2, '1')
    if polynomial != polynomial.H:
        raise ArithmeticError('the actual source polynomial is not Hermitian')
    c = settings()
    scalar_tail = s.sqrt(2)*c['binomial_operator_ratio']**3/(1-c['binomial_operator_ratio'])
    # Each spectral value is checked separately; no rounded eigenvalue admission.
    scalar_errors = []
    for eigen in (0, 1, s.sqrt(3)):
        approx = binomial_polynomial(s.Matrix([[eigen]]), 2, '1')[0]
        exact = s.sqrt(2)*eigen/s.sqrt(eigen*eigen+1)
        scalar_errors.append(s.simplify((exact-approx)**2))
        if s.simplify(scalar_tail**2-(exact-approx)**2).is_nonnegative is not True:
            raise ArithmeticError('exact smooth approximation exceeds its global tail')
    return {'native_one_body_dimension': 60, 'HW_squared_eigenvalues': [0, 1, 3],
            'full_spectral_ranks': ranks,
            'normalized_smooth_squared_eigenvalues': [0, 1, s.Rational(3, 2)],
            'smooth_target_annihilates_all16_Wilson_zero_modes': bool((target*p0).applyfunc(s.simplify) == s.zeros(60)),
            'degree2_scalar_error_squares': scalar_errors,
            'uniform_tail_bound_squared': s.simplify(scalar_tail**2),
            'point_background_claimed_a_quantum_preparation': False}


def certificate():
    source = state.density_plan(degree=0)
    config = settings('1', '1', source['declared_plaquette_coefficient'])
    plan = local_record_plan(config)
    volume = finite_volume_plan(config, plan)
    # The inherited cutoff planner depends only on maximum coordinate
    # degree and bounded norm. Its legacy degree key is a filtration
    # parameter here; no old Legendre coefficient is imported.
    filtration = {
        'Legendre_degree': plan['binomial_polynomial_degree'],
        'single_link': source['single_link'],
        'full_Fock_kinetic_norm_per_cell_upper': config['full_Fock_kinetic_density_rational_upper'],
        'declared_plaquette_coefficient': config['plaquette_coefficient']}
    finite = cutoff_plan(filtration, volume, accuracy_bits=20)
    # Same initial mixed state and fixed kappa, but the action and its
    # kinetic density allowance are explicitly changed to the smooth one.
    rough = 48000*(source['initial_magnetic_cost_per_cell_upper']+(
        source['initial_electric_energy_per_cell_upper']
        +2*config['full_Fock_kinetic_density_rational_upper'])/config['plaquette_coefficient'])
    return state.serial({
        'schema': 'oph.smooth_sign_field.v1', 'source_inputs': list(SOURCE_INPUTS),
        'smooth_action': config, 'fixed_time_polynomial_comparison': plan,
        'same_polynomial_finite_volume_plan': volume,
        'same_polynomial_finite_word_plan': finite,
        'total_complete_local_record_error_upper': plan['complete_local_record_TV_upper']+volume['complete_record_TV_upper']+finite['complete_record_TV_upper'],
        'initial_state': 'same referenced site-twirled product from particle_local_gauge_state, L fixed before observation',
        'per_link_preparation_degree': source['single_link']['per_link_word_degree'],
        'smooth_action_recomputed_roughness_bound': min(F(1), rough),
        'roughness_used_as_premise_of_record_comparison': False,
        'exact_native_control': exact_native_control(),
        'flat_zero_energy_squared': free_symbol_squared((1, 1, 1)),
        'free_first_corner_energy_squared': free_symbol_squared((-1, 1, 1)),
        'free_linear_Dirac_coefficient': '1 in the declared a,time convention',
        'infinite_smoothing_action': '(beta+H_W)/a: finite-range Wilson, not a physical chiral claim',
        'fixed_delta_complete_local_thermodynamic_limit': True,
        'same_as_original_Legendre_or_discontinuous_action': False,
        'finite_word_cutoff_applied_after_polynomial_approximation': True,
        'physical_GeV_mass_chiral_or_Lorentz_limit_claimed': False})
