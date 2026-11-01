from fractions import Fraction as F

import pytest

from . import energy


def test_certified_energies_gaps_weights_and_distinct_energy_origins():
    result = energy.certificate()
    rows = result['exact_positive_probe_lines']
    assert len(rows) == 8
    assert sum(row['signed_gap_to_R_approximation'] < 0 for row in rows) == 5
    assert all(F(row['rigorous_gap_to_N6_ground_interval'][0]) > 0 for row in rows)
    assert rows[0]['gap_to_minimum_probe_energy_interval'] == ['0', '0']
    assert result['minimum_probe_energy_is_full_N6_ground'] is False
    assert sum(row['normalized_probe_weight_approximation'] for row in rows) == pytest.approx(1, abs=5e-12)
    for row in rows:
        low, high = map(F, row['H_energy_interval'])
        a, b = map(F, row['signed_gap_to_R_interval'])
        assert (low-3, high-3) == (a, b)


def test_ground_estimate_is_separate_from_exact_variational_bracket():
    result = energy.ground_certificate()
    low, high = map(F, result['rigorous_ground_interval'])
    assert float(low) < -.64 < result['numerical_lowest_eigenvalue_estimate'] < -.3 < float(high)
    assert result['numerical_lowest_eigenvalue_estimate'] == pytest.approx(-.3669041542038966, abs=1e-12)
    assert result['sparse_estimate_is_an_interval_certified_ground'] is False
    assert result['full_Fock_ground_across_all_number_sectors_determined'] is False


def test_mass_map_retains_hbar_length_and_time_units():
    result = energy.mass_certificate()
    assert result['conditional_inertial_mass'] == 'M=2 hbar^2/(E_unit L_unit^2)'
    assert result['required_time_map'] == 't_phys=T_unit t, E_unit T_unit=hbar'
    assert result['negative_binding_offset_is_rest_mass'] is False
    assert result['SI_or_GeV_mass_identified'] is False


def test_old_calibrations_are_not_promoted_to_current_predictions():
    result = energy.calibration_audit()
    assert len(result['sources_sha256']) == 6
    assert result['old_scale_inherited_by_cube_or_moving_pair'] is False
    assert result['measured_mass_used_to_fit_new_candidates'] is False
    assert result['old_dimer_scale']['role'].endswith('not an absolute mass prediction')
