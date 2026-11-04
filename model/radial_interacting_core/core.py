"""Exact rational bounds for the same bounded radial interacting core."""
from fractions import Fraction as F

from ..native_mass_channels.algebra import rational
from ..radial_defect_gap.bounds import spectral_window


def core_bounds(coupling=F(1, 65536), radius=1, spacing=1, velocity=1,
                mass=64, hopping_bound=2, dimensions=3):
    if (type(radius) is not int or not 1 <= radius <= 7
            or type(dimensions) is not int or dimensions not in (1, 3)):
        raise ValueError('radius 1..7 and dimension one or three required')
    coupling = rational(coupling)
    if coupling <= 0:
        raise ValueError('positive exact centered contact required')
    source = spectral_window(spacing, velocity, mass, hopping_bound)
    volume = (2*radius+1)**dimensions
    gap = source['one_body_gap_lower']
    alpha = source['origin_weight_lower']
    ipr = source['ipr_lower']
    norm = coupling*volume/4
    if 4*norm >= gap:
        raise ValueError('complete-Fock contact requires ||W|| < d/4')
    epsilon = norm*norm/(gap-2*norm)
    angle = norm/(gap-2*norm)
    lower, upper = coupling*ipr/2-2*epsilon, coupling/2+2*epsilon
    certified = lower > 0
    if not certified:
        raise ValueError('the actual singly occupied ground is not resolved by this bound')
    effect_error = angle
    return {
        'radius': radius, 'dimensions': dimensions, 'sites': volume,
        'spacing': rational(spacing), 'velocity': rational(velocity),
        'mass': rational(mass), 'hopping_bound': rational(hopping_bound),
        'full_Fock_dimension': 2**(2*volume), 'complete_one_body_dimension': 2*volume,
        'one_body_gap_lower': gap,
        'origin_weight_lower': alpha, 'origin_weight_upper': F(1),
        'ipr_lower': ipr, 'ipr_upper': F(1),
        'centered_contact_coupling': coupling, 'full_contact_norm': norm,
        'projected_charged_gap_lower': coupling*ipr/2,
        'projected_charged_gap_upper': coupling/2,
        'ordered_low_energy_error': epsilon, 'projection_angle_bound': angle,
        'actual_charged_gap_lower': lower, 'actual_charged_gap_upper': upper,
        'charged_gap_certified_positive': certified,
        'exact_low_Kramers_ground_multiplicity': 2 if certified else None,
        'ground_particle_number': volume if certified else None,
        'charged_particle_numbers': [volume-1, volume+1],
        'other_charged_excitations_threshold_lower': gap-2*norm,
        'ground_unhalved_trace_error': 2*angle,
        'charged_unhalved_trace_error': 2*angle,
        'each_pole_matrix_operator_error': 3*angle,
        'each_pole_matrix_rank': 2 if 3*angle < F(1, 2) else None,
        'pole_rank_two_certified': 3*angle < F(1, 2),
        'each_zero_mode_low_pole_weight_lower': max(F(0), F(1, 2)-3*angle),
        'each_zero_mode_low_pole_weight_upper': F(1, 2),
        'each_zero_mode_high_energy_weight_upper': 3*angle,
        'each_origin_spin_low_pole_weight_lower': alpha/2-3*angle,
        'each_origin_spin_low_pole_weight_upper': F(1, 2)+3*angle,
        'center_parity_single_core_ideal_upper': (1-alpha*alpha)/2,
        'center_parity_empty_or_full_core_ideal_lower': (1+alpha*alpha)/2,
        'center_parity_actual_ground_lower': F(0),
        'center_parity_actual_ground_upper': min(F(1), (1-alpha*alpha)/2+effect_error),
        'center_parity_actual_charged_lower': max(F(0), (1+alpha*alpha)/2-effect_error),
        'center_parity_actual_charged_upper': F(1),
        'center_parity_actual_charged_minus_ground_lower': alpha*alpha-2*effect_error,
        'center_parity_effect_error': effect_error,
        'center_contact_charged_minus_ground_lower': coupling*(alpha*alpha-2*effect_error)/2,
        'entire_complement_retained': True,
        'supplied_bounded_radial_mass_required': True,
        'volume_uniform_interacting_gap_proved': False,
    }


def faithful_reference(bits=20, **parameters):
    if type(bits) is not int or not 1 <= bits <= 1024:
        raise ValueError('integer thermal accuracy bits 1..1024 required')
    row = core_bounds(**parameters)
    if not row['charged_gap_certified_positive']:
        raise ValueError('a positive complete-Fock ground gap is required')
    if row['actual_charged_gap_lower'] >= row['other_charged_excitations_threshold_lower']:
        raise ArithmeticError('the declared full gap ordering failed')
    beta = F(2*row['sites']+bits)/row['actual_charged_gap_lower']
    return {'inverse_temperature': beta, 'accuracy_bits': bits,
            'ground_mixture_unhalved_trace_distance_upper': F(1, 2**bits),
            'any_local_effect_expectation_error_upper': F(1, 2**(bits+1)),
            'actual_finite_interacting_Gibbs_state_faithful': True,
            'stationary_under_same_Hamiltonian': True,
            'modular_to_Hamiltonian_time_factor': -beta,
            'stationary_under_native_repair_inferred': False,
            'Gibbs_preparation_or_equilibration_executed': False,
            'geometric_modular_flow_identified': False}


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def certificate():
    retained = core_bounds()
    if (retained['exact_low_Kramers_ground_multiplicity'] != 2
            or retained['origin_weight_lower'] != F(7688, 7691)
            or retained['center_parity_actual_ground_upper'] >= F(1, 1000)
            or retained['center_parity_actual_charged_lower'] <= F(999, 1000)
            or retained['each_origin_spin_low_pole_weight_lower'] <= F(49, 100)):
        raise ArithmeticError('same-source radial ground, charged pole or local contrast failed')
    return {'same_bounded_radial_source_examples': encode([core_bounds(mass=4), retained]),
            'retained_27_site_mass64_source': encode(retained),
            'exact_center_parity_ideal_single': '(1-alpha^2)/2',
            'exact_center_parity_ideal_empty_or_full': '(1+alpha^2)/2',
            'exact_center_contact': 'W_center=(U/2)*(n_center-1)^2-U/4',
            'faithful_thermal_reference': encode(faithful_reference()),
            'interaction_background_same_as_exact_bounded_radial_zero_modes': True,
            'finite_linear_texture_spectrum_substituted': False,
            'actual_center_parity_is_an_exact_zero_mode_occupation': False,
            'two_to_54_Fock_matrix_executed': False,
            'numerical_eigenvalues_used_for_admission': False}
