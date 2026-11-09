"""Factored finite resources for fixed-time complete local polynomial records."""
from fractions import Fraction as F
from math import comb

import sympy as s

from model.particle_local_gauge_state import state
from model.particle_field_causality import bounds
from model.particle_prepared_field_readout import experiment

SOURCE_INPUTS = (
    'model/particle_local_gauge_state/state.py',
    'model/particle_local_gauge_state/PROOF.md',
    'model/particle_field_causality/bounds.py',
    'model/particle_field_causality/PROOF.md',
    'model/particle_observer_outputs/field_route.py',
    'model/particle_observer_outputs/FIELD_ROUTE.md',
    'model/particle_prepared_field_readout/experiment.py',
)


def ceil_log2(value):
    """Smallest integer b with value<=2^b, using exact positive rationals."""
    x = state.rational(value, 'positive logarithm argument', True)
    b = x.numerator.bit_length()-x.denominator.bit_length()
    power = F(2**b) if b >= 0 else F(1, 2**(-b))
    return b+int(x > power)


def ceil_fraction(value):
    return (value.numerator+value.denominator-1)//value.denominator


def local_marginal_patch(cells):
    """Exact incidence data determining the full retained-cell marginal."""
    cells = tuple(sorted(set(tuple(x) for x in cells)))
    if not cells or any(len(x) != 3 or any(type(v) is not int for v in x) for x in cells):
        raise ValueError('a nonempty set of integer cubic cells is required')
    vertices = set(cells)
    links = []
    for x in cells:
        for axis in range(3):
            y = tuple(v+int(j == axis) for j, v in enumerate(x))
            links.append((x, y)); vertices.add(y)
    return {'cells': cells, 'retained_outgoing_links': tuple(links),
            'independent_Haar_vertices': tuple(sorted(vertices)),
            'fermion_density': 'I on all retained60-mode cell factors / their full dimension',
            'outside_gauge_variables_and_factors_cancel_by_partial_trace': True}


def volume_plan(source, duration='1', accuracy_bits=24, observed_radius=0,
                observed_cells=1, intermediate_instruments=0):
    """Fixed-time locality bound without allocating exponentials or 2^R."""
    t = state.rational(duration, 'duration', True)
    bits = state.integer(accuracy_bits, 'accuracy bits', 1)
    r = state.integer(observed_radius, 'observed radius')
    count = state.integer(observed_cells, 'observed cell count', 1)
    k = state.integer(intermediate_instruments, 'intermediate instruments')
    if count > (2*r+1)**3:
        raise ValueError('observed cell count exceeds its cube')
    constants = bounds.interaction_constants(source['Legendre_degree'], source['spacing'],
                                            source['declared_plaquette_coefficient'])
    d, b = constants['root_support_radius'], max(constants['root_support_radius'], 2)
    rate = constants['weighted_local_interaction_upper']
    # e^(2 rate t)-1 <=2^ceil(4 rate t), since e<4.
    exp_bits = ceil_fraction(4*rate*t)
    R = r+2*b+1
    attempts = 0
    while True:
        D = R-r-2*b
        crossing = (2*constants['rooted_kinetic_norm_upper']*(2*(R+d)+1)**3
                    +3*constants['plaquette_norm']*(2*(R+2)+1)**3)
        prefactor = (2*k+1)*2*count*crossing*t
        exponent = ceil_log2(prefactor)+exp_bits-D
        attempts += 1
        if exponent <= -bits:
            break
        R *= 2
    return {'duration': t, 'accuracy_bits': bits, 'box_radius': R,
            'observed_radius': r, 'observed_cells': count,
            'intermediate_instruments': k, 'locality_constants': constants,
            'boundary_separation_lower': D, 'boundary_norm_upper': crossing,
            'exponential_binary_exponent_upper': exp_bits,
            'local_record_error_binary_exponent_upper': exponent,
            'complete_record_TV_upper': F(1, 2**bits),
            'periodic_side_sufficient': 2*(R+2*b+d)+1,
            'radius_search_steps': attempts,
            'duration_shrunk_with_box_or_volume': False,
            'large_exponential_or_2_to_radius_materialized': False,
            'locality_uses_same_fixed_source_parameters': True}


def cutoff_plan(source, volume, accuracy_bits=24):
    """Factorial-tail enclosure with integer resources, no giant factorials."""
    bits = state.integer(accuracy_bits, 'cutoff accuracy bits', 1)
    R, t = volume['box_radius'], volume['duration']
    cells, links = (2*R+1)**3, 3*(2*R+1)**3
    L = source['single_link']['per_link_word_degree']
    k0 = links*L
    B = cells*source['full_Fock_kinetic_norm_per_cell_upper']+links*source['declared_plaquette_coefficient']
    x = B*t
    intervals = volume['intermediate_instruments']+1
    first = max(1, ceil_fraction(4*x), 3*(bits+2))
    step = max(4, 2*source['Legendre_degree']+1)
    K = k0+step*(first-1)
    exponent = 2-first//3
    if exponent > -bits:
        raise ArithmeticError('the factorial-tail exponent exceeds its allocation')
    return {'initial_total_word_degree_upper': k0, 'total_word_cutoff': K,
            'word_step': step, 'first_omitted_Dyson_order': first,
            'integrated_bounded_action_norm_upper': x,
            'bounded_action_norm_upper': B, 'evolution_intervals': intervals,
            'whole_history_error_binary_exponent_upper': exponent,
            'complete_record_TV_upper': F(1, 2**bits),
            'gauge_dimension_upper': {'base': 1+450*links, 'exponent': K},
            'full_Fock_dimension': {'base': 2, 'exponent': 60*cells},
            'apparatus_menu': 'complete local occupation PVMs; coordinate-degree preserving',
            'all_original_outcomes_and_quantum_references_retained': True,
            'degree_preserving_pointers_are_interleaved_before_Dyson_count': True,
            'factorial_or_full_matrix_materialized': False}


def source_current_control():
    """Exact local two-time curvature on the full maximally mixed Fock input."""
    # One independent pair of CAR modes, before applying the actual native
    # hopping singular-value census. It checks the trace factor1/2.
    c1 = s.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
    c2 = s.Matrix([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, -1], [0, 0, 0, 0]])
    current = s.I*(c1.H*c2-c2.H*c1)
    if s.trace(current)/4 != 0 or s.trace(current*current)/4 != s.Rational(1, 2):
        raise ArithmeticError('the complete CAR trace factor differs')
    native = experiment.native_escape_coefficient()['squared_escape_coefficient']
    variance = native/2
    return {'Legendre_degree': 0, 'spacing': F(1),
            'single_directed_bond_Tr_BdaggerB': F(480, 289),
            'single_bond_current_mean': F(0), 'single_bond_current_variance': F(240, 289),
            'origin_number_derivative_mean_square': variance,
            'two_time_squared_number_increment_second_derivative_at_zero': 2*variance,
            'initial_complete_number_probabilities': [F(comb(60, n), 2**60) for n in range(61)],
            'two_time_number_menu_labels': 61**2,
            'all_initial_and_final_occupation_outcomes_retained': True,
            'finite_code_requires_K_at_least_K0_plus1': True,
            'finite_nonzero_time_probability_enclosure_claimed': False,
            'continuum_scattering_or_mass_claimed': False}


def certificate():
    source = state.density_plan(degree=0)
    volume = volume_plan(source, '1', accuracy_bits=24)
    cutoff = cutoff_plan(source, volume, accuracy_bits=24)
    return state.serial({
        'schema': 'oph.local_history_family.v1', 'source_inputs': list(SOURCE_INPUTS),
        'fixed_source_density_parameters': source,
        'local_marginal_control': local_marginal_patch([(0, 0, 0), (1, 0, 0)]),
        'fixed_time_volume_plan': volume, 'complete_finite_cutoff_plan': cutoff,
        'total_record_error_upper': volume['complete_record_TV_upper']+cutoff['complete_record_TV_upper'],
        'same_source_local_current_control': source_current_control(),
        'local_polynomial_record_limit_constructed_at_fixed_spacing': True,
        'full_canonical_purifiers_identified_across_volumes': False,
        'full_dynamical_exterior_output_compared': False,
        'external_reference_in_default_comparison': 'one independent untouched reference; correlated extensions require matching whole-box/reference marginal',
        'roughness_bound_reused_after_local_interventions': False,
        'exact_sign_or_relativistic_limit_claimed': False,
        'large_finite_simulation_executed': False})
