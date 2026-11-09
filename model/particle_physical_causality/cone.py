"""Exact refinement-uniform causal majorants, without a continuum claim.

The clock is Hamiltonian duration and physical distance is a times the
cubic graph distance. The kinetic and magnetic actions scale at most 1/a;
arbitrary onsite electric generators are removed in the interaction picture.
"""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.particle_smooth_sign_field import smooth
from model.particle_field_causality import bounds
from model.particle_local_history_family.history import ceil_log2


SOURCE_INPUTS = (
    'model/particle_smooth_sign_field/smooth.py',
    'model/particle_smooth_sign_field/PROOF.md',
    'model/particle_field_causality/bounds.py',
    'model/particle_field_causality/PROOF.md',
    'model/particle_local_history_family/history.py',
    'model/particle_local_history_family/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/interacting_extension/chiral_lattice/source_data.py',
    'model/common_source/overlap_regulator.py',
    'model/full_fock_field/README.md',
    'model/source_family_compatibility/LOCAL_MATTER_READOUT.md',
    'Lean/QFT/FiniteCausalObserverNet.lean',
    'Lean/Geometry/RefiningLatticeCausalCone.lean',
)


def family(kind='smooth', delta='1', magnetic_scale='1/100', weight='2'):
    """W_a <= W_star/a for kappa(a) <= magnetic_scale/a.

    The smoothing parameter is fixed across spacings. For ordinary Wilson,
    weight is any rational greater than one; delta has no role in its action.
    Neither family inherits a vacuum, Lorentz boost or chiral measure.
    """
    k0 = bounds.rational(magnetic_scale, 'nonnegative magnetic scale')
    if kind == 'smooth':
        source = smooth.settings(delta, '1', k0)
        lam = source['spatial_weight_base']
        W = source['weighted_full_interaction_norm_upper']
        action = 'H_E(a) + dGamma((beta+c_delta s_delta(H_W))/a) - kappa(a) sum ReTr(U_p)/15'
        smoothing = source['delta']
    elif kind == 'wilson':
        lam = bounds.rational(weight, 'spatial weight', True)
        if lam <= 1:
            raise ValueError('spatial weight must exceed one')
        W = 360*lam+48*k0*lam*lam
        action = 'H_E(a) + dGamma((beta+H_W)/a) - kappa(a) sum ReTr(U_p)/15'
        smoothing = None
    else:
        raise ValueError('choose the smooth or wilson declared action')
    # Integral_1^lambda 1/x dx >= (lambda-1)/lambda.
    mu = (lam-1)/lam
    return {'kind': kind, 'action': action, 'fixed_delta': smoothing,
            'magnetic_scale_upper': k0, 'spatial_weight': lam,
            'log_weight_rational_lower': mu,
            'weighted_rate_times_spacing_upper': W,
            'physical_speed_rational_upper': 2*W/mu,
            'physical_speed_formula': '2 W_star / log(lambda)',
            'magnetic_scaling_premise': '0 <= kappa(a) <= magnetic_scale_upper/a',
            'onsite_electric_coefficient_restricted_by_bound': False,
            'measured_light_speed_or_mass_input': False,
            'fixed_delta_source_is_discontinuous_sign': False}


def support_count(radius, spacing):
    """Alignment-safe cell count for a cube of physical half-width radius."""
    R = bounds.rational(radius, 'physical support radius')
    a = bounds.rational(spacing, 'positive spacing', True)
    # A length 2R interval meets at most floor(2R/a)+2 integer centers;
    # one extra cell leaves the outgoing-link ownership convention explicit.
    return (int(2*R/a)+3)**3


def validate_family(config):
    """Reject a mutated interaction rate or logarithm lower bound."""
    try:
        expected = family(config['kind'], config['fixed_delta'] or '1',
                          config['magnetic_scale_upper'], config['spatial_weight'])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('a declared field-family configuration is required') from exc
    if config != expected:
        raise ValueError('field-family constants differ from their source derivation')
    return config


def spacing_member(config, spacing, magnetic=None, electric_coefficient='1'):
    """Instantiate a member and verify its actual magnetic scaling premise."""
    validate_family(config)
    a = bounds.rational(spacing, 'positive spacing', True)
    upper = config['magnetic_scale_upper']/a
    k = upper if magnetic is None else bounds.rational(magnetic, 'magnetic coefficient')
    e = bounds.rational(electric_coefficient, 'onsite electric coefficient')
    if k > upper:
        raise ValueError('magnetic coefficient violates the declared refinement envelope')
    lam = config['spatial_weight']
    if config['kind'] == 'smooth':
        rate = smooth.settings(config['fixed_delta'], a, k)['weighted_full_interaction_norm_upper']
    else:
        rate = 360*lam/a+48*k*lam*lam
    if a*rate > config['weighted_rate_times_spacing_upper']:
        raise ArithmeticError('source interaction exceeds its physical-rate envelope')
    return {'spacing': a, 'magnetic_coefficient': k, 'electric_coefficient': e,
            'weighted_interaction_rate_upper': rate,
            'rate_times_spacing_upper': a*rate,
            'electric_coefficient_enters_causal_rate': False}


def history_bound(config, spacing, duration, separation, physical_radius='0',
                  intervention_scale=0, intervention_power=0):
    """All local records versus one remote unitary or retained dilation.

    k(a)<=K*a^-p complete local dilations are allowed. The remote private
    apparatus is discarded when comparing the observer output. Every local
    outcome and arbitrary untouched reference is retained. A remote record
    handed directly to the observer is outside this locality experiment.
    """
    validate_family(config)
    a = bounds.rational(spacing, 'spacing', True)
    T = bounds.rational(duration, 'duration')
    d = bounds.rational(separation, 'physical support separation', True)
    R = bounds.rational(physical_radius, 'physical support radius')
    K = bounds.integer(intervention_scale, 'intervention scale')
    p = bounds.integer(intervention_power, 'intervention power')
    if a > 1:
        raise ValueError('the polynomial support envelope requires 0<a<=1')
    W = config['weighted_rate_times_spacing_upper']
    mu = config['log_weight_rational_lower']
    gap = mu*d-2*W*T
    degree = p+3
    C = 2*(2*K+1)*(2*R+3)**3
    # exp(-gap/a)<=2^-floor(gap/a) for gap>=0, since e>2.
    exponent = None if gap < 0 else ceil_log2(C/a**degree)-int(gap/a)
    certified_bits = max(0, -exponent) if exponent is not None else 0
    return {'spacing': a, 'duration': T, 'physical_separation_lower': d,
            'physical_support_cube_half_width': R,
            'alignment_safe_observed_cell_count_upper': support_count(R, a),
            'intermediate_local_dilation_count_upper': K/a**p,
            'support_and_apparatus_prefactor': C,
            'inverse_spacing_prefactor_power': degree,
            'physical_exponential_margin_lower': gap,
            'error_binary_exponent_upper': exponent,
            'certified_error_bits': certified_bits,
            'complete_local_history_half_trace_upper':
                {'base': 2, 'exponent': -certified_bits},
            'full_Fock_all_outcomes_and_untouched_references': True,
            'remote_support_cardinality_in_bound': False,
            'same_state_and_clock_within_each_comparison': True,
            'includes_boundary_representation_or_implementation_error': False}


def spacing_plan(config, duration, separation, physical_radius='0', accuracy_bits=20,
                 intervention_scale=0, intervention_power=0):
    """Dyadic a0 with the bound valid for EVERY real 0<a<=a0.

    The monotonic envelope C*a^-q*exp(-g/a) increases on a<=g/q.
    Requiring this domain certifies more than a list of sampled spacings.
    """
    bits = bounds.integer(accuracy_bits, 'accuracy bits', 1)
    a = F(1)
    result = history_bound(config, a, duration, separation, physical_radius,
                           intervention_scale, intervention_power)
    g = result['physical_exponential_margin_lower']
    if g <= 0:
        raise ValueError('positive physical cone margin is required by this certificate')
    q = result['inverse_spacing_prefactor_power']
    steps = 0
    while a > g/q or result['certified_error_bits'] < bits:
        a /= 2
        steps += 1
        result = history_bound(config, a, duration, separation, physical_radius,
                               intervention_scale, intervention_power)
    return {**result, 'all_smaller_spacing_upper': a,
            'dyadic_refinements': steps, 'requested_error_bits': bits,
            'monotonic_envelope_condition': a <= g/q,
            'uniform_for_every_real_spacing_below_cut': True,
            'physical_duration_held_fixed': True,
            'continuum_state_or_observable_convergence_proved': False}


def error_composition(parent_bound, left_errors=(), right_errors=()):
    """Separate complete-experiment errors on both sides add by triangle inequality."""
    values = [bounds.rational(parent_bound, 'parent half-trace error')]
    values += [bounds.rational(x, 'left complete-experiment error') for x in left_errors]
    values += [bounds.rational(x, 'right complete-experiment error') for x in right_errors]
    if any(x > 1 for x in values):
        raise ValueError('each complete-experiment half-trace bound must be at most one')
    return min(F(1), sum(values))


@lru_cache(maxsize=1)
def native_wilson_bonds():
    """Actual source spin matrices and all fifteen internal columns.

    Unitary gauge-coordinate factors preserve the bond singular values.
    The complete Fock norm is the sum of positive one-body eigenvalues.
    """
    data = bounds.field_route.source_data.load()
    gamma = [s.Matrix(g).applyfunc(lambda z: s.nsimplify(z, rational=True))
             for g in data['gamma'][0]]
    beta = gamma[3]
    rows = []
    for mu in range(3):
        P = (s.eye(4)-gamma[mu])/2
        B = -beta*P
        H = s.zeros(8)
        H[:4, 4:] = B
        H[4:, :4] = B.H
        spectrum = H.eigenvals()
        if P != P.H or P*P != P or B.H*B != P:
            raise ArithmeticError('the native Wilson projector convention differs')
        if spectrum != {-1: 2, 0: 4, 1: 2}:
            raise ArithmeticError('the complete two-site source spectrum differs')
        internal = len(data['charges'])
        if internal != 15:
            raise ArithmeticError('native internal dimension differs')
        rows.append({'axis': mu, 'spin_projector_rank': int(P.trace()),
                     'one_body_spectrum_per_internal_column': {'-1': 2, '0': 4, '1': 2},
                     'all_internal_columns': internal,
                     'two_site_one_body_dimension': 8*internal,
                     'full_Fock_bond_norm_times_spacing': 2*internal})
    return rows


def serial(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [serial(v) for v in x]
    return x


def certificate():
    fixed = family('smooth')
    wilson = family('wilson')
    # Choose a displayed physical separation with margin exactly one.
    # This exposes the conservative smooth speed instead of shrinking time.
    smooth_d = (2*fixed['weighted_rate_times_spacing_upper']+1)/fixed['log_weight_rational_lower']
    smooth_plan = spacing_plan(fixed, '1', smooth_d, '1', 20, 2)
    wilson_plan = spacing_plan(wilson, '1/1000', '4', '1/10', 20, 2)
    return serial({
        'schema': 'oph.physical_causality.v1', 'source_inputs': list(SOURCE_INPUTS),
        'fixed_delta_smooth_family': fixed, 'smooth_fixed_time_example': smooth_plan,
        'ordinary_Wilson_family': wilson, 'Wilson_fixed_time_example': wilson_plan,
        'native_Wilson_bond_census': native_wilson_bonds(),
        'complete_detector': 'all local occupation outcomes and two retained local apparatus records; remote Wilson-plaquette phase versus identity',
        'example_input': 'same full source/reference state in both experiments; no translation or vacuum assumption',
        'actual_large_field_history_allocated': False,
        'uniform_physical_cone_for_stated_scaling': True,
        'physical_support_growth_included': True,
        'strict_zero_tail_at_finite_spacing': False,
        'whole_degree_compression_assumed_local': False,
        'Lorentz_covariance_continuum_vacuum_or_chirality_proved': False,
        'measured_light_speed_or_particle_mass_inferred': False,
        'remaining_refinement_premise': 'compatible renormalized states, observables and interacting records must converge across spacings; a causal bound alone does not prove this',
    })
