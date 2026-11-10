"""Exact channel labels and explicit limits of particle identification."""
from .candidate_data import candidate_data, certificate_data


def test_exact_candidate_channels_and_no_species_overclaim():
    result = candidate_data()
    charged = result['charged_channel_checks']
    assert charged['exact_gauge_generator_weights'] == ['0'] * 11 + ['1']
    assert charged['CAR_anticommutator_column_checks'] == 2 * 4 * 4 * 4096
    assert charged['reflected_conjugate_operator_columns'] == 4 * 4096
    assert charged['spin_reversal_and_parity_columns'] == 4096
    assert len(result['charged_probes']) == 4
    assert all(row['physical_species_assignment'] is None for row in result['charged_probes'])
    assert not charged['generation_count_or_four_distinct_species_inferred']
    assert not charged['physical_EM_current_identified']
    neutral = result['neutral_channels']
    assert neutral['actual_operator_current_commutators'] == 32
    assert neutral['scalar_occupation']['bounded_projector']
    assert neutral['Yukawa_composite']['safe_operator_norm_upper_bound'] == '8'
    assert not neutral['source_spin_scalar_implies_relativistic_scalar_particle']
    assert result['separate_periodic_control']['dimension'] == 3686400
    assert not result['separate_periodic_control']['open_ground_or_energy_spectrum_transferred']
    assert not result['separate_periodic_control']['larger_even_ring_menu_requirements']['menu_counts_suffice_for_a_particle_dispersion']
    assert not result['empirical_fit_or_confidence_assigned']
    family = result['native_active_family']
    assert [row['hypercharge'] for row in family['rows']] == ['1', '-1/2', '-1/2']
    assert [row['weak_T3'] for row in family['rows']] == ['0', '1/2', '-1/2']
    assert [row['algebraic_T3_plus_Y'] for row in family['rows']] == ['1', '0', '-1']
    assert not family['quark_particle_spectra_calculated_here']
    assert certificate_data() == result
