"""Attach the residue theorem to two independently proved full sources."""
from fractions import Fraction as F

from ..topological_defect_spectrum import core as textured
from ..native_odd_box_core import core as original
from . import bounds


def source_result(profile, coupling=F(1, 65536)):
    if profile == 'linear_texture_cube':
        parent = textured.core_bounds(coupling, 3)
        detector = F(8, 27)
    elif profile == 'unmodified_odd_cube':
        parent = original.core_bounds(coupling=coupling)
        detector = F(1, 8)
    else:
        raise ValueError('an explicitly proved source profile is required')
    if (parent['sites'] != 27 or parent['full_Fock_dimension'] != 2**54
            or parent['exact_low_Kramers_ground_multiplicity'] != 2):
        raise ArithmeticError('the parent source does not supply the stated ground')
    value = bounds.residue_bounds(
        sites=parent['sites'], coupling=coupling,
        gap_lower=parent['one_body_and_Fock_excitation_gap_lower'],
        inverse_participation=parent['zero_orbital_inverse_participation'],
        detector_probability=detector)
    if (value['full_contact_norm'] != parent['full_contact_norm']
            or value['charged_gap_lower'] != parent['actual_charged_gap_lower']
            or value['charged_gap_upper'] != parent['actual_charged_gap_upper']
            or not value['detector_pole_certified_visible']):
        raise ArithmeticError('pole energies or weights do not match the source theorem')
    value['imposed_mass_texture_required'] = profile == 'linear_texture_cube'
    value['boundary_standing_wave_source'] = profile == 'unmodified_odd_cube'
    return value


def certificate():
    sources = {name: source_result(name) for name in (
        'linear_texture_cube', 'unmodified_odd_cube')}
    archives = {
        name: bounds.poisson_archive(
            edges=54, rate=F(1), time=1/value['charged_gap_lower'],
            ground_angle=value['projection_angle_bound'])
        for name, value in sources.items()}
    continuing = {
        name: bounds.poisson_archive(
            edges=54, rate=value['charged_gap_lower']/54000,
            time=1/value['charged_gap_lower'],
            ground_angle=value['projection_angle_bound'])
        for name, value in sources.items()}
    for value in continuing.values():
        if value['negative_log_no_event_probability'] != F(1, 1000) or value['rate_per_seam'] <= 0:
            raise ArithmeticError('the declared continuing repair window changed')
    return bounds.encode({
        'complete_source_poles': sources,
        'unit_rate_time_one_over_gap_lower_archive_controls': archives,
        'declared_positive_rate_continuing_windows': continuing,
        'continuing_window_exact_ground_output_distance_upper': F(1, 1000),
        'same_gap_in_addition_and_removal': True,
        'all_virtual_Fock_states_retained': True,
        'positive_matrix_measures_not_only_eigenvalues': True,
        'full_residue_rank_is_not_a_physical_species_count': True,
        'full_residue_eigenvalue_equals_zero_mode_compression_assumed': False,
        'finite_box_gap_identified_as_thermodynamic_mass': False,
    })
