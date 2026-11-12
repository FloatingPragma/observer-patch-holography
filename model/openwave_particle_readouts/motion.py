"""Atom-resolved charged arrivals after a retained position intervention."""
from functools import lru_cache

import sympy as s

from model.coherent_defect_wave.chiral_source.model import PAULI
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import response, polynomial
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import add, mul, divide, evaluate, outward, sqrt_interval
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.local import positive_response
from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.algebra import coefficient, routes
from model.openwave_energy_resolved_transport import response as energy_resolved
from model.openwave_native_charge_flux import algebra as native_flux


TIME = s.Integer(3)
HORIZONS = (0, 1, 3, 6)


@lru_cache(maxsize=1)
def spin_effects():
    """Contract the actual Pauli hopping with a selected local spectral atom."""
    w, c = s.symbols('w c', real=True)

    def tensor(i, j, scalar, vector):
        return (scalar*int(i == j)*s.eye(2)
                + vector*sum((a[j, i]*a/2 for a in PAULI), s.zeros(2)))

    effects = []
    for axis, sigma in enumerate(PAULI):
        effect = s.zeros(4)
        for spin in range(2):
            for alpha, beta in sigma.todok():
                for alphap, betap in sigma.todok():
                    effect += (s.conjugate(sigma[alphap, betap])*sigma[alpha, beta]
                               * s.kronecker_product(
                                   tensor(alphap, spin, w, c).H*tensor(alpha, spin, w, c),
                                   tensor(betap, beta, w, c))/32)
        exchange = sum(((1 if j == axis else -1)*s.kronecker_product(a/2, a/2)
                        for j, a in enumerate(PAULI)), s.zeros(4))
        expected = ((w**3/16+3*w*c*c/64)*s.eye(4)
                    + c*c*(2*w-c)*exchange/16)
        if (effect-expected).applyfunc(s.expand) != s.zeros(4):
            raise ArithmeticError('selected-atom native spin contraction failed')
        effects.append(effect.applyfunc(s.expand))
    return (w, c), tuple(effects)


def detuned_arrival_bounds(moment, coupling, time, vector_error,
                           edge_norm=15, hopping_norm=360):
    """Complete Dyson bound with the actual unit atomic flux detuning."""
    values = tuple(map(s.sympify, (*moment, coupling, time, vector_error, edge_norm, hopping_norm)))
    if any(value.is_Rational is not True for value in values):
        raise ValueError('exact rational source bounds are required')
    lo, hi, lam, t, epsilon, edge, tau = values
    if not (0 < lo <= hi <= edge**2 and lam > 0 and 0 <= t <= 6
            and epsilon >= 0 and 0 < edge <= tau):
        raise ValueError('positive source controls and a time between zero and six required')
    # On 0<=t/2<=3 the omitted tail decreases. Its first term is positive.
    argument = t/2
    sine_low = 2*sum((-1)**k*argument**(2*k+1)/s.factorial(2*k+1) for k in range(12))
    sine_high = sine_low+2*argument**25/s.factorial(25)
    root_low, root_high = sqrt_interval((lo, hi), 100)
    error = edge*t*epsilon+lam*tau*tau*t*t/2
    low = max(0, sine_low*root_low-error)
    high = sine_high*root_high+error
    return {'unit_detuning_amplitude_interval': (sine_low, sine_high),
            'scaled_amplitude_error': error,
            'strict_scaled_amplitude_margin': sine_low*root_low-error,
            'probability_interval': (lam*lam*low*low, lam*lam*high*high)}


def residue_interval(numerator, denominator, atom, ground):
    """Outward evaluation of an algebraic-field residue at a simple real root."""
    value = (s.Integer(0), s.Integer(0))
    for power in range(numerator.degree(), -1, -1):
        value = add(mul(value, atom), evaluate(polynomial(numerator.rep.nth(power)), ground))
    return divide(value, evaluate(denominator.diff(), atom))


def sampling_budget(early, late):
    """Original-input repetition counts for independent complete preparations."""
    values = tuple(map(s.sympify, (*early, *late)))
    if any(value.is_Rational is not True for value in values):
        raise ValueError('exact rational event-probability intervals are required')
    low, high, later_low, later_high = values
    if not (0 < low <= high <= 1 and 0 <= later_low <= later_high < low):
        raise ValueError('positive separated event-probability intervals are required')
    contrast = low-later_high
    detection = s.ceiling(3/low)
    comparison = s.ceiling(20*(high+later_high)/contrast**2)
    exponential_partial_sum = sum(s.Rational(3)**k/s.factorial(k) for k in range(9))
    if exponential_partial_sum <= 20:
        raise ArithmeticError('the exact exponential control must imply exp(-3)<1/20')
    return {'detection_trials': detection, 'comparison_trials_per_setting': comparison,
            'contrast_lower': contrast,
            'comparison_error_upper': (high+later_high)/(comparison*contrast**2),
            'exp_three_partial_sum_lower': exponential_partial_sum}


def validate_ec_pole_source(factor, weight, blocks):
    """Regenerate the shared atomic result; transfer no preparation bound."""
    actual = energy_resolved.native_pole_residue()
    native_factor = s.Poly(s.sympify(actual['selected_atomic_pole_polynomial']), factor.gen)
    if (native_factor != factor
            or actual['pump_pole_weight_interval'] != outward(weight, 100)
            or actual['native_removed_spin_target_dimensions'] != [block['H'].rows for block in blocks]
            or actual['actual_selected_target_spin'] != 1
            or not actual['exact_residue_difference_divisible_by_selected_polynomial']):
        raise ArithmeticError('the regenerated native ec atom does not match the shared spectral source')
    if not energy_resolved.tagged_spin_effects()['all_three_axes_exact']:
        raise ArithmeticError('the shared ec tagged-hop contraction failed')
    return {'provider': 'model/openwave_energy_resolved_transport/response.py',
            'complete_native_pole_and_weight_regenerated_and_equal': True,
            'all_three_ec_tagged_hop_contractions_regenerated': True,
            'preparation_or_timeline_probability_bound_transferred': False}


@lru_cache(maxsize=1)
def ec_flux_attachment():
    """Apply the all-state native Gauss identity to each atom-resolved subeffect."""
    endpoints = native_flux.endpoint_certificate()
    census = native_flux.gauss_cut_census()
    if (not endpoints['hypercharge_is_central_on_matter_scalar_and_both_link_endpoints']
            or endpoints['complete_119_state_interior_hypercharge_endpoint_sum'] != '0'
            or census['independent_cell_Gauss_comparisons'] != 24840):
        raise ArithmeticError('the complete native Gauss-cut operator identity failed')
    provider = native_flux.route_certificate()
    available = {row['edge']: row for row in provider['actual_routes']}
    actual_routes = routes()
    if set(available) != {row['edge'] for row in actual_routes}:
        raise ArithmeticError('charge flux and motion require the same actual six routes')
    rows = []
    for route in actual_routes:
        attached = available[route['edge']]
        if any(attached[key] != value for key, value in route.items()):
            raise ArithmeticError('a charge-flux route belongs to a different native source')
        links = [0]*24
        links[route['edge']] = attached['forced_full_link_basis_index']
        for mask in range(256):
            expected = (mask & 1)-((mask >> route['arrival']) & 1)
            if (native_flux.cut_flux_six(mask, links) != 6*expected
                    or s.Rational(attached['all_256_cut_fluxes'][mask]) != expected):
                raise ArithmeticError('the full native cut flux disagrees with the attached charge')
        if (attached['arrival_source_matter_plus_scalar_Y_and_Q'] != '-1'
                or attached['origin_reference_Y_and_Q'] != '1'
                or not attached['all_twelve_current_constraints_checked_at_both_endpoints']):
            raise ArithmeticError('the arrival requires the complete native origin/reference constraint')
        rows.append({'edge': route['edge'], 'arrival': route['arrival'],
                     'arrival_source_matter_plus_scalar_Y_and_Q': '-1',
                     'reference_Y_and_Q': '1',
                     'three_half_cube_fluxes': attached['three_half_cube_fluxes'],
                     'all_256_cut_fluxes': f'1_(0 in S)-1_({route["arrival"]} in S)',
                     'charge_flux_mismatch_probability': '0'})
    return {'provider': 'model/openwave_native_charge_flux/algebra.py',
            'regenerated_complete_endpoint_Ward_identities': endpoints['full_endpoint_Ward_identities'],
            'regenerated_additive_factor_basis_Gauss_comparisons': census['independent_cell_Gauss_comparisons'],
            'regenerated_route_cut_identities': 6*256,
            'all_state_domain': 'every joint-physical vector in each native ec arrival effect P_e; the atom-resolved A_e is a subeffect',
            'routes': rows,
            'complete_refinement': 'for every coarse record r, jointly measure the commuting central-Y cut-flux and total-charge projectors, retaining all outcomes; sum_(k,l) p(r,k,l)=p(r)',
            'joint_probability': 'p(r,k,l)=Tr[R_(S,l) Pi_(S,k) K_r rho_in K_r^dagger Pi_(S,k) R_(S,l)]',
            'all_branch_charge_flux_mismatch_probability': '0 for k != l',
            'all_preparation_failures_and_complete_nonarrival_refined': True,
            'other_preparation_arrival_probability_transferred': False,
            'arbitrary_nonarrival_Q_measurement_assumed': False,
            'unconditioned_mean_flux_or_origin_return_inferred': False}


@lru_cache(maxsize=2)
def atom_first_hop(family):
    """Derive both local residues from the complete native charged resolvents."""
    if family not in ('ec', 'L'):
        raise ValueError('the actual ec/L removal channels are required')
    positive = positive_response(family)
    signed = response(family)
    if positive['q'] != signed['q']:
        raise ArithmeticError('positive and spin residues must retain the same target support')
    atom = positive['rows'][0]['energy_interval']
    weight = positive['rows'][0]['weight_interval']
    spin = residue_interval(signed['numerator'], signed['q'], atom, positive['e_interval'])
    # Select the actual irreducible target factor by its isolated root.
    selected = []
    for factor, multiplicity in s.factor_list(positive['q'])[1]:
        if multiplicity != 1:
            raise ArithmeticError('the squarefree target annihilator is required')
        for (lo, hi), root_multiplicity in factor.intervals(eps=s.Rational(1, 10**70)):
            if max(lo, atom[0]) <= min(hi, atom[1]):
                if root_multiplicity != 1:
                    raise ArithmeticError('the selected spectral root must be simple')
                selected.append(factor)
    if len(selected) != 1:
        raise ArithmeticError('the complete target support must identify one selected factor')
    factor = selected[0]
    z = signed['numerator'].gen
    field = signed['numerator'].domain
    field_factor = s.Poly(factor.as_expr().subs(factor.gen, z), z, domain=field)
    # The positive numerator has twice the instrument weight at each atom.
    ratio = s.Rational(-2, 3) if family == 'ec' else s.Integer(2)
    if not (signed['numerator']-ratio*positive['numerator']/2).rem(field_factor).is_zero:
        raise ArithmeticError('selected spin residue is not the claimed native atom residue')
    expected_spin = mul(weight, (ratio, ratio))
    if max(spin[0], expected_spin[0]) > min(spin[1], expected_spin[1]):
        raise ArithmeticError('independent residue intervals disagree')
    native = coefficient(family)
    x = native['cube_4H_ground_interval']
    weight_cubed = mul(mul(weight, weight), weight)
    if family == 'ec':
        correction = mul(x, (s.Rational(1, 648),)*2)
        moment = mul(weight_cubed, add((s.Rational(1, 12),)*2, correction))
    else:
        moment = mul(weight_cubed, (s.Rational(1, 4),)*2)
    if not 0 < moment[0] <= moment[1] < 225:
        raise ArithmeticError('the full selected-route moment must be positive')
    common_source = validate_ec_pole_source(factor, weight, positive['blocks']) if family == 'ec' else None
    return {'family': family, 'atom': atom, 'weight': weight, 'spin_residue': spin,
            'spin_to_weight_ratio': ratio, 'selected_factor': factor,
            'cube_energy': x, 'moment': moment, 'shared_ec_atomic_source': common_source}


@lru_cache(maxsize=2)
def family_certificate(family):
    # This delayed import allows the owning readout certificate to consume motion.
    from .readouts import inputs, prepared_line_bounds

    data = inputs(family)
    source = data['source']
    local = atom_first_hop(family)
    epsilon = data['error']
    initial = prepared_line_bounds(data['weight'], epsilon, source['v'], 0)
    success_low, success_high = initial['probability']
    bounds = detuned_arrival_bounds(local['moment'], source['neutral']['coupling'], TIME,
                                   epsilon, 15, source['neutral']['tau'])
    lower, upper = bounds['probability_interval']
    q = data['preparation_probability_lower']
    if bounds['strict_scaled_amplitude_margin'] <= 0 or success_low <= 6*upper:
        raise ArithmeticError('prepared atom requires positive arrivals and retained nonarrival')
    energy_loss = min(s.Integer(1), epsilon**2/success_low)
    times = []
    time_bounds = {}
    for time in HORIZONS:
        controlled = detuned_arrival_bounds(local['moment'], source['neutral']['coupling'], time,
                                             epsilon, 15, source['neutral']['tau'])
        lo, hi = controlled['probability_interval']
        time_bounds[time] = (q*lo, hi)
        times.append({'source_time': str(time),
                      'unit_detuning_amplitude_interval': list(map(str, controlled['unit_detuning_amplitude_interval'])),
                      'per_route_probability_from_original_input': outward((q*lo, hi), 192),
                      'six_route_probability_from_original_input': outward((6*q*lo, 6*hi), 192),
                      'six_route_probability_from_original_input_decimal': [str(s.N(v, 16)) for v in (6*q*lo, 6*hi)],
                      'strict_scaled_amplitude_margin': str(controlled['strict_scaled_amplitude_margin'])})
    if not time_bounds[3][0] > time_bounds[6][1] or not time_bounds[1][0] > 0:
        raise ArithmeticError('the same-H excursion and return contrast is not certified')
    resources = sampling_budget(tuple(6*v for v in time_bounds[3]),
                                tuple(6*v for v in time_bounds[6]))
    rows = [
        {'record': 'physical_preparation_failure', 'exact_probability': '1-q'},
        {'record': 'native_pump_failure', 'exact_probability': 'q p_failure'},
        {'record': 'complete_energy_bin_complement', 'exact_probability': 'q p_outside'},
        {'record': 'line_position_preparation_failure', 'exact_probability': 'q(p_line-s)',
         'probability_interval': outward((0, epsilon**2), 160)},
    ]
    for route in routes():
        rows.append({'record': 'line_prepared_charged_arrival', **route,
                     'exact_probability': f'q r_edge_{route["edge"]}',
                     'probability_interval': outward((q*lower, upper), 192),
                     'static_reference_cell': 0})
    rows.append({'record': 'position_prepared_complete_nonarrival',
                 'exact_probability': 'q(s-sum_edge r_edge)',
                 'probability_interval': outward((q*(success_low-6*upper), success_high), 160)})
    flux = ec_flux_attachment() if family == 'ec' else None
    if flux is not None:
        flux_routes = {row['edge']: row for row in flux['routes']}
        for row in rows:
            if row['record'] == 'line_prepared_charged_arrival':
                row['joint_native_charge_flux_atom'] = flux_routes[row['edge']]
                row['joint_atom_probability_equals_actual_arrival_probability'] = True
            else:
                row['joint_native_flux_refinement'] = 'all central-Y cut-flux and charge outcomes retained; physical mismatches have zero probability'
        for row in times:
            row['joint_native_charge_flux_six_route_probability_from_original_input'] = row['six_route_probability_from_original_input']
    return {
        'family': family,
        'selected_atomic_energy_interval': outward(local['atom'], 80),
        'atom_weight_interval': outward(local['weight'], 100),
        'selected_spin_residue_interval': outward(local['spin_residue'], 100),
        'exact_spin_residue_divided_by_atom_weight': str(local['spin_to_weight_ratio']),
        'selected_complete_target_factor': str(local['selected_factor'].as_expr()),
        'exact_atom_route_moment': 'w^3(1/12+x/648)' if family == 'ec' else 'w^3/4',
        'atom_route_moment_interval': outward(local['moment'], 100),
        'atom_route_moment_decimal': str(s.N(local['moment'][0], 16)),
        'shared_ec_atomic_source': local['shared_ec_atomic_source'],
        'joint_native_charge_flux': flux,
        'prepared_line_to_atomic_atom_vector_error_upper': str(epsilon),
        'position_success_probability_conditional_on_physical_preparation': outward((success_low, success_high), 100),
        'position_failure_probability_before_line_renormalization_upper': str(epsilon**2),
        'energy_bin_loss_after_position_conditioned_on_position_success_upper': str(energy_loss),
        'energy_bin_loss_upper_decimal': str(s.N(energy_loss, 16)),
        'energy_bin_loss_bound_preserved_at_every_source_clock_time': True,
        'exact_energy_bin_support_after_position_asserted': False,
        'unit_detuning_amplitude_interval': list(map(str, bounds['unit_detuning_amplitude_interval'])),
        'first_order_scaled_amplitude_error_upper': str(bounds['scaled_amplitude_error']),
        'strict_scaled_arrival_amplitude_margin': str(bounds['strict_scaled_amplitude_margin']),
        'per_route_probability_conditional_only_on_physical_preparation': outward((lower, upper), 192),
        'per_route_probability_from_original_input': outward((q*lower, upper), 192),
        'six_route_probability_from_original_input': outward((6*q*lower, 6*upper), 192),
        'six_route_probability_from_original_input_decimal': [str(s.N(v, 16)) for v in (6*q*lower, 6*upper)],
        'arrival_probability_after_position_success': outward((6*lower/success_high, min(1, 6*upper/success_low)), 192),
        'clock_records': times,
        'time_three_arrival_lower_exceeds_time_six_upper': True,
        'per_route_excursion_return_probability_contrast_lower': str(time_bounds[3][0]-time_bounds[6][1]),
        'independent_complete_preparation_sampling': {
            'event': 'union of all six retained arrival routes',
            'trial_unit': 'one independently prepared original physical input, including every preparation and measurement failure',
            'independent_source_apparatus_repetitions_explicitly_required': True,
            'time_three_detection_trials': str(resources['detection_trials']),
            'detection_no_arrival_probability_upper': 'exp(-3)<1/20',
            'exact_exp_three_partial_sum_lower': str(resources['exp_three_partial_sum_lower']),
            'time_three_versus_six_trials_per_setting': str(resources['comparison_trials_per_setting']),
            'comparison_total_original_input_trials': str(2*resources['comparison_trials_per_setting']),
            'wrong_order_probability_upper': '1/20',
            'method': 'variance bound for independent Bernoulli sample means and Chebyshev inequality',
            'physical_preparation_rate_or_feasible_runtime_inferred': False,
        },
        'complete_records': {
            'records': rows, 'coarse_record_count': len(rows),
            'definitions': 'q=physical preparation success; p_line,p_outside,p_failure are conditional on q; s=joint line/position success; r_edge=joint line/position/arrival, also conditional only on q',
            'exact_completeness': '(1-q)+q(p_failure+p_outside+p_line-s)+q sum_edge r_edge+q(s-sum_edge r_edge)=1',
            'all_source_copies_reference_Fock_space_and_preparation_environments_retained': True,
            'all_joint_source_reference_Gauss_constraints_preserved': True,
        },
    }


def certificate():
    spin_effects()
    return {
        'source_time': str(TIME),
        'source_clock_settings': list(map(str, HORIZONS)),
        'source': 'same complete N2 H0+2^-34 T+4 N_ref with fixed reference at cell zero',
        'position_preparation': 'S0=all links blank and true local ground at every cell except zero; identity on full origin CAR/scalar and reference',
        'complete_position_effects': ['S0', 'I-S0'],
        'selected_arrival_effects': 'one incident link excited, all other links blank; true local ground at every nonarrival cell, complete local spectral projector for the selected charged atom at arrival',
        'complete_arrival_complement_retained': True,
        'exact_zero_hopping_kernel': 'P_edge exp(-it H0) S0=0 for every input and all real t',
        'exact_no_operation_arrival_probability': '0',
        'energy_filter_measurement_order': {
            'this_preparation': 'S0 P_I(H) B_r Psi, with S0 after the full-H energy filter',
            'energy_resolved_transport_preparation': 'P_I(H) B_ec P0 Psi, whose energy filter creates a positive time-zero tagged tail',
            'same_atomic_ec_result_different_explicit_preparation': True,
            'time_zero_dressing_probability_or_preparation_error_bound_transferred': False,
        },
        'full_hopping_norm_upper': '360',
        'selected_edge_norm_upper': '15',
        'actual_atomic_arrival_detuning': '1, the unchanged nonblank-link energy',
        'complete_Dyson_second_order_norm_upper': '(lambda*tau*t)^2/2',
        'three_axis_native_spin_identities': 3,
        'exact_two_spin_effect': '(w^3/16+3w*c^2/64)I+c^2*(2w-c) S0 dot R_axis Sy/16',
        'exact_cube_route_moment': 'w^3/16+3w*c^2/64+c^2*(2w-c)*x/768',
        'prepared_motion': [family_certificate(family) for family in ('ec', 'L')],
        'scope': {
            'strict_positive_same_H_arrival_after_energy_then_position_preparation': True,
            'position_intervention_and_energy_disturbance_explicit': True,
            'unchanged_source_clock_and_fixed_quantum_compensator': True,
            'arrival_claim_for_unmeasured_original_line': False,
            'arrival_claim_conditioned_on_a_final_energy_measurement': False,
            'off_resonant_flux_excursion_and_finite_return_contrast': True,
            'charged_band_dispersion_velocity_or_continuum_limit': False,
        },
    }
