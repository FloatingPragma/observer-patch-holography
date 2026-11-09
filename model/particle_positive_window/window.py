"""Separate surviving positive energy from a zero atom or escaping weight."""
from fractions import Fraction as F

from model.particle_common_contract import contract
from model.particle_joint_reconstruction import positive_energy
from model.particle_spatial_ground_state import ground

SOURCE_INPUTS = (
    'model/particle_positive_window/PROOF.md',
    'model/particle_positive_window/window.py',
    'model/particle_joint_reconstruction/POSITIVE_ENERGY.md',
    'model/source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md',
    'model/particle_yang_mills_gap_contract/PROOF.md',
)


def positive_window(contrast_estimate, contrast_error, time_abs_upper,
                    mass_upper, first_absolute_moment_upper,
                    negative_tail_upper, lower_energy, upper_energy):
    """Conditional on certified original-weight record and moment bounds.

    The estimate is F(0)-Re F(t). Error covers both scores and all apparatus
    errors. No division by a success probability is admitted here.
    """
    estimate = ground._signed(contrast_estimate, 'signed contrast estimate')
    error = ground._q(contrast_error, 'contrast error')
    time = ground._q(time_abs_upper, 'absolute source time upper bound')
    mass = ground._q(mass_upper, 'original spectral mass upper bound')
    moment = ground._q(first_absolute_moment_upper, 'first absolute moment upper bound')
    negative = ground._q(negative_tail_upper, 'negative spectral tail upper bound')
    delta = ground._q(lower_energy, 'positive lower energy', True)
    radius = ground._q(upper_energy, 'upper energy', True)
    if radius <= delta:
        raise ValueError('the upper energy must exceed the positive lower energy')
    lower = estimate-error
    if lower > 2*mass:
        raise ValueError('the contrast interval contradicts the spectral mass upper bound')
    infrared = mass*(delta*time)**2/2
    ultraviolet = 2*moment/radius
    numerator = lower-infrared-2*negative-ultraviolet
    return {'original_positive_band': [delta, radius],
            'contrast_lower': lower,
            'low_frequency_contrast_allowance': infrared,
            'negative_frequency_contrast_allowance': 2*negative,
            'high_positive_frequency_contrast_allowance': ultraviolet,
            'original_positive_band_weight_lower': max(F(0), numerator/2),
            'contrast_and_moment_hypotheses_supplied': True,
            'success_normalization_applied': False,
            'isolated_particle_or_mass_shell_proved': False}


def selected_requirement(member, contrast_estimate, contrast_error,
                         time_abs_upper='1', temporal_half_width='1',
                         lower_energy='1/4', upper_energy='16',
                         operator_norm_squared='1'):
    """Supply source-derived allowances, leaving the actual record contrast explicit."""
    contract.verify_member(member)
    filtered = positive_energy.triangular_filter_bounds(
        member, temporal_half_width, operator_norm_squared)
    tail = positive_energy.selected_tail_bound(member, lower_energy, operator_norm_squared)
    result = positive_window(contrast_estimate, contrast_error, time_abs_upper,
        filtered['spectral_total_mass_upper'],
        filtered['spectral_first_absolute_moment_upper'],
        tail['original_negative_spectral_weight_upper'], lower_energy, upper_energy)
    return contract.serial({
        'schema': 'oph.primary_positive_window_requirement.v1',
        'member_id': member['member_id'], **member['ids'],
        'temporal_filter': filtered, 'negative_tail': tail,
        'conditional_window_bound': result,
        'contrast_is_supplied_not_computed_for_primary': True,
        'primary_positive_energy_weight_derived': False,
        'physical_particle_or_mass_assignment': None,
    })


def exact_phase_controls():
    """Independent finite measures at t=pi/2; |t|<11/7.

    Integer frequencies give exact cosine values. These controls use their
    own finite diagonal Hamiltonians, not the selected spatial source.
    """
    controls = {
        'positive_excitation': [(1, F(1))],
        'zero_atom': [(0, F(1))],
        'ultraviolet_escape': [(101, F(1))],
        'small_negative_tail': [(-1, F(1, 100)), (1, F(99, 100))],
        'mixed_complete_response': [(0, F(1, 4)), (1, F(1, 2)), (3, F(1, 4))],
    }
    cosine = (1, 0, -1, 0)
    output = []
    delta, radius = F(1, 4), F(16)
    for name, measure in controls.items():
        mass = sum(weight for _, weight in measure)
        moment = sum(abs(energy)*weight for energy, weight in measure)
        negative = sum(weight for energy, weight in measure if energy <= -delta)
        contrast = sum((1-cosine[energy % 4])*weight for energy, weight in measure)
        result = positive_window(contrast, 0, F(11, 7), mass, moment, negative, delta, radius)
        actual = sum(weight for energy, weight in measure if delta <= energy <= radius)
        if not 0 <= result['original_positive_band_weight_lower'] <= actual:
            raise ArithmeticError('the exact finite spectral control violates the window bound')
        output.append({'name': name, 'original_measure': measure,
            'source_time': 'pi/2', 'absolute_time_upper': F(11, 7),
            'actual_original_band_weight': actual, 'bound': result,
            'is_selected_spatial_source': False})
    return output


def certificate():
    member = contract.member(2)
    return contract.serial({
        'schema': 'oph.positive_window_certificate.v1',
        'selected_requirements_without_a_positive_contrast_witness':
            selected_requirement(member, 0, 0),
        'exact_finite_phase_controls': exact_phase_controls(),
        'primary_noncollapse_or_species_identification_proved': False,
    })
