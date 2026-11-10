"""Outward complete records and capacity for the source-derived free member."""
from fractions import Fraction as F
import math
import mpmath as mp
import sympy as s

from .source import admit, preparation_parameters
from .algebra import complete_two_tests
from .controller import resource_bounds


def endpoint(value, upper=False):
    raw = value._mpi_[int(upper)]
    sign, mantissa, exponent, _ = raw
    return (-1 if sign else 1)*F(mantissa)*(F(2)**exponent)


def ivf(value):
    value = F(value)
    return mp.iv.mpf(value.numerator)/value.denominator


def outward(value, upper=False, places=30):
    value = F(value); scale = 10**places
    n = -(-value.numerator*scale//value.denominator) if upper else value.numerator*scale//value.denominator
    return str(F(n, scale))


def bounds(value):
    return [outward(endpoint(value)), outward(endpoint(value, True), True)]


def expression(value):
    value = s.sympify(value)
    if value.is_Rational:
        return ivf(F(int(value.p), int(value.q)))
    if value.is_Add:
        return sum(map(expression, value.args), mp.iv.mpf(0))
    if value.is_Mul:
        return math.prod(map(expression, value.args))
    if value.is_Pow and value.exp == s.Rational(1, 2):
        return mp.iv.sqrt(expression(value.base))
    if value.is_Pow and value.exp.is_Integer:
        return expression(value.base)**int(value.exp)
    raise ValueError('unsupported exact probability expression: '+str(value))


def finite_member(side=32, preparation_depth=3):
    if not isinstance(side, int) or side < 3 or not isinstance(preparation_depth, int) or preparation_depth < 0:
        raise ValueError('integer side>=3 and nonnegative preparation depth required')
    mp.iv.dps = 90
    parent, _ = admit()
    r0, contraction, length = preparation_parameters()
    V = side**3; blocks = 4*V; trials = blocks*2**preparation_depth
    filters = blocks*(2**preparation_depth-1)
    error_squared = 2*blocks*r0*contraction**preparation_depth
    error = mp.iv.sqrt(ivf(error_squared))
    error_upper = min(F(1), endpoint(error, True))
    log_success = -sum((blocks*2**(preparation_depth-j)*mp.iv.ln(1+ivf(r0*contraction**j))
                       for j in range(preparation_depth)), mp.iv.mpf(0))
    success_lower = endpoint(mp.iv.exp(log_success))
    assert 0 < success_lower <= 1 and 0 < error_upper <= 1
    actual = complete_two_tests()
    control = complete_two_tests(False)
    probabilities = []
    for row in actual['all_sixteen_records']:
        probability = expression(row['ideal_actual_source_probability'])
        p_lo, p_hi = endpoint(probability), endpoint(probability, True)
        # The unselected apparatus spin partners remain exactly empty in the
        # reference-relative frame for every source input, including its complement.
        spectator_forbidden = row['first_down'] or row['second_down']
        lo = F(0) if spectator_forbidden else max(F(0), p_lo-error_upper)
        hi = F(0) if spectator_forbidden else min(F(1), p_hi+error_upper)
        probabilities.append({**row, 'ideal_probability_interval': bounds(probability),
            'conditional_native_preparation_interval': [outward(lo), outward(hi, True)],
            'unconditional_native_preparation_interval': [outward(success_lower*lo), outward(hi, True)]})
    # Same normalized p/m=3/4 mode and length8pi/(3m), with side-dependent a.
    spacing_mass = 8*mp.iv.pi/(3*side)
    p_over_m = ivf(F(3, 4))
    omega_over_m = mp.iv.sqrt(1+4*mp.iv.sin(mp.iv.pi/side)**2/spacing_mass**2)
    duration_mass = mp.iv.pi/(2*omega_over_m)
    omega_error_over_m = spacing_mass**2*p_over_m**4/24
    spacing_record_error = mp.iv.pi/4*duration_mass*omega_error_over_m
    # Every nontrivial time dependence in the literal circuit comes from
    # interference of the vacuum and the a†b† pair in the first no-click branch.
    quarter_hit = expression(actual['active_four_records'][1]['ideal_actual_source_probability'])
    zero_hit = expression(control['active_four_records'][1]['ideal_actual_source_probability'])
    comparison_phase = ivf(F(5, 4))*duration_mass
    comparison_hit = zero_hit+(quarter_hit-zero_hit)*mp.iv.sin(comparison_phase)**2
    contrast = quarter_hit-zero_hit
    assert endpoint(contrast) > 0
    reference_dimension = 48*V
    # Four retained fermionic apparatus modes, two complete doublets.
    detector_log_capacity = mp.iv.ln(reference_dimension)+4*mp.iv.ln(2)
    log_native = 240*mp.iv.ln(2)+8*mp.iv.ln(5)+24*mp.iv.ln(119)
    log_constituent = log_native+mp.iv.ln(2**31+1)
    main_log_capacity = blocks*log_constituent+detector_log_capacity
    # Each deterministic native trial retains8 complete four-state spin
    # references. Every agreement-filter node retains its full M-state record.
    prepared_log_capacity = trials*(log_constituent+16*mp.iv.ln(2))+filters*mp.iv.ln(length)+detector_log_capacity
    attempts = 6
    amplified_failure = (1-success_lower)**attempts
    amplified_capacity = attempts*(prepared_log_capacity-detector_log_capacity)+detector_log_capacity+mp.iv.ln(attempts+1)
    return {'macro_side': side, 'macro_sites': V, 'main_native_constituents': blocks,
        'main_native_cells': 8*blocks, 'selected_canonical_Fock_modes': 4*V,
        'selected_Fock_dimension_expression': '2^'+str(4*V),
        'mass_parameter': 'm=E_star, the selected actual native eigenenergy',
        'mass_source_interval': parent['intervals']['native_selected_E_interval'],
        'source_box_side': '8 pi/(3 E_star)',
        'spacing_times_E_star_interval': bounds(spacing_mass),
        'selected_momentum_over_E_star': ['0', '0', '3/4'],
        'actual_mode_frequency_over_E_star_interval': bounds(omega_over_m),
        'actual_observation_time_times_E_star_interval': bounds(duration_mass),
        'finite_spacing_complete_record_error_upper': outward(endpoint(spacing_record_error, True), True),
        'same_window_continuum_second_only_probability_interval': bounds(comparison_hit),
        'actual_second_only_probability_interval': bounds(quarter_hit),
        'zero_elapsed_time_second_only_probability_exact': '1/20',
        'actual_time_response_contrast_interval': bounds(contrast),
        'prepared_conditional_time_response_contrast_lower': outward(endpoint(contrast)-2*error_upper),
        'complete_sixteen_detector_records': probabilities,
        'native_preparation': {'deterministic_initial_trial': True, 'base_trial_success': '1',
            'filter_depth_counts_all_rounds': True, 'binary_filter_depth': preparation_depth,
            'initial_excited_odds_upper': str(r0), 'geometric_odds_contraction': str(contraction),
            'retained_purification_vector_error_upper': outward(error_upper, True),
            'retained_purification_error_squared_upper_exact': str(error_squared),
            'source_trials_retained': trials, 'complete_filter_records_retained': filters,
            'each_filter_pointer_dimension': length,
            'all_preparation_success_probability_lower': outward(success_lower),
            'all_preparation_failure_probability_upper': outward(1-success_lower, True),
            'success_lower_formula': 'product_(j=0..d-1)(1+r0*c^j)^(-4V*2^(d-j))',
            'unconditional_complete_law': '(1-Q) preparation-failure; Q times each of16 complete detector probabilities',
            'unused_source_partners_and_spin_references_retained': True},
        'optional_complete_preparation_amplification': {
            'independent_full_attempts_retained': attempts,
            'first_success_index_dimension': attempts+1,
            'all_attempts_failed_probability_upper': outward(amplified_failure, True),
            'complete_record_TV_to_exact_ground_example_upper': outward(amplified_failure+error_upper, True),
            'failed_attempts_are_never_discarded': True,
            'source_trial_count': attempts*trials,
            'complete_source_preparation_reference_and_pointer_capacity_nats_interval': bounds(amplified_capacity),
            'capacity_formula': 'R[4V*2^d(logD_constituent+16ln2)+4V*(2^d-1)lnM]+ln(48V)+4ln2+ln(R+1)',
            'actual_source_clock_work_and_controller': resource_bounds(side, preparation_depth, attempts)},
        'actual_source_clock_work_and_controller': resource_bounds(side, preparation_depth),
        'capacity': {'native_source_and_measurement_nats_interval': bounds(main_log_capacity),
            'all_preparation_sources_references_records_and_measurement_nats_interval': bounds(prepared_log_capacity),
            'constituent_formula': '240 ln2+8 ln5+24 ln119+ln(2^31+1)',
            'source_and_measurement_formula': '4V logD_constituent+ln(48V)+4ln2',
            'all_preparation_formula': '4V*2^d*(logD_constituent+16ln2)+4V*(2^d-1)lnM+ln(48V)+4ln2',
            'regular_spatial_reference_dimension': reference_dimension,
            'retained_detector_Fock_dimension': 16,
            'gate_synthesis_and_additional_workspace_included': False,
            'capacity_is_full_kinematic_log_dimension_not_state_entropy': True,
            'memory_budget_derives_lattice_or_preparation_parameters': False},
        'execution': {'full_native_Hamiltonian_numerically_executed': False,
            'literal_active_mode_circuit_dimension': 16,
            'full_source_four_mode_and_two_reference_doublet_dimension': 256,
            'spectator_vacuum_reduction_is_exact': True,
            'selected_E_star_eigenvectors_numerically_executed': False}}
