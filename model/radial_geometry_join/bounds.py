"""Exact charged-sector to finite-geometry contrast on one Fock source."""
from fractions import Fraction as F

from model.radial_interacting_core.core import core_bounds


def joined_bounds():
    core = core_bounds(coupling=F(1, 65536), radius=1, spacing=F(1),
                       velocity=F(1), mass=F(64), hopping_bound=F(2))
    alpha, angle = core['origin_weight_lower'], core['projection_angle_bound']
    error = (1-alpha**2)/2+angle
    contrast = F(3, 40)-F(192, 5)*error
    if core['sites'] != 27 or core['exact_low_Kramers_ground_multiplicity'] != 2:
        raise ArithmeticError('the actual bounded radial core theorem does not match the join')
    if error >= F(1, 1000) or contrast <= F(59, 1000):
        raise ArithmeticError('the same-source charged geometry contrast is not certified')
    return {
        'sites': core['sites'], 'full_field_Fock_dimension': 2**54,
        'mass_texture': F(64), 'hopping_velocity': F(1), 'lattice_spacing': F(1),
        'contact_coupling': F(1, 65536), 'origin_zero_weight_lower': alpha,
        'interacting_ground_projection_angle_upper': angle,
        'actual_charged_gap_lower': core['actual_charged_gap_lower'],
        'actual_charged_gap_upper': core['actual_charged_gap_upper'],
        'ground_D0_expectation_interval': [F(0), error],
        'each_charged_D0_expectation_interval': [1-error, F(1)],
        'charged_minus_ground_D0_expectation_lower': 1-2*error,
        'charged_minus_ground_record1_probability_lower': F(24, 25)*(1-2*error),
        'charged_minus_ground_local_contact_energy_lower': F(1, 131072)*(1-2*error),
        'Gram_G01_Lipschitz_constant_in_D0_expectation': F(96, 5),
        'ground_Gram_G01_minus_each_charged_Gram_G01_strict_lower': contrast,
        'simple_geometry_contrast_strict_lower': F(59, 1000),
        'same_actual_ground_and_charged_Fock_projectors_used': True,
        'charged_signs_have_equal_parity_and_geometry_by_particle_hole_symmetry': True,
        'local_contact_energy_equals_total_excitation_energy_assumed': False,
        'fresh_independent_ground_or_charged_inputs_required_for_mean_iteration': True,
        'autonomous_reciprocal_geometry_dynamics_derived': False,
        'Lorentzian_metric_or_Einstein_equation_inferred': False,
    }


def observation_bound():
    # Hoeffding on two independent endpoint means: 4 exp(-2 N delta^2).
    from math import factorial
    delta, count = F(1, 4000), 48000000
    exponent = 2*count*delta**2
    lower_exp = sum((exponent**j/factorial(j) for j in range(21)), F(0))
    if exponent != 6 or lower_exp <= 400:
        raise ArithmeticError('the complete finite observation confidence bound changed')
    gram_error = 20*(delta+F(1, 10**9))
    return {'fresh_attempts_per_preparation': count,
        'record_probability_sampling_error': delta,
        'complete_control_trace_distance_strict_upper': F(1, 10**9),
        'inferred_ideal_Gram_error_per_preparation_upper': gram_error,
        'inferred_ideal_Gram_contrast_strict_lower': F(59, 1000)-2*gram_error,
        'joint_failure_probability_strict_upper': F(1, 100),
        'all_failure_records_retained_with_declared_binary_score_zero': True,
        'large_observation_or_spectral_state_preparation_executed': False}
