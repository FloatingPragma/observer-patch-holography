from fractions import Fraction as F
import copy
import pytest
from . import calibrated_triplet as c


def test_complete_three_mass_calibration_and_exact_circle():
    out = c.certificate()
    x, y = F(out['phase_real']), F(out['phase_imaginary'])
    assert x*x+y*y == 1
    assert 0 < F(out['pass_amplitude']) < 1
    assert out['complete_retained_flags'] == 16
    assert F(out['loss_probability']) > 0
    assert out['phase_in_balanced_positive_chamber']
    assert out['all_three_masses_used_as_inputs']
    assert out['held_out_mass_predictions'] == 0
    assert not out['source_selected_controller']
    assert not out['physical_mass_matching_proved']
    for mass, box, error in zip(out['input_mass_energies_GeV'],
                               out['controller_output_mass_energies_GeV'],
                               out['absolute_mass_error_bounds_GeV']):
        lo, hi = map(F, box)
        assert max(abs(lo-F(mass)), abs(hi-F(mass))) <= F(error) < F(1, 10**36)
    assert F(out['target_Koide_quotient_interval'][1]) < F(2, 3)
    assert F(out['controller_Koide_quotient']) < F(2, 3)


def test_tau_is_a_fit_input_unlike_the_one_ratio_sanity_check():
    a = c.fit(['0.00051099895069', '0.1056583755', '1.77693'])
    b = c.fit(['0.00051099895069', '0.1056583755', '1.77692'])
    assert a['pass_amplitude'] != b['pass_amplitude']
    assert a['phase_half_angle_tangent'] != b['phase_half_angle_tangent']
    assert a['common_mass_scale_GeV'] != b['common_mass_scale_GeV']


def test_scale_changes_only_overall_mass_unit_to_certified_precision():
    a = c.fit(['0.00051099895069', '0.1056583755', '1.77693'])
    b = c.fit([F(x)*4 for x in a['input_mass_energies_GeV']])
    assert a['pass_amplitude'] == b['pass_amplitude']
    assert a['phase_half_angle_tangent'] == b['phase_half_angle_tangent']
    assert abs(F(b['common_mass_scale_GeV'])-4*F(a['common_mass_scale_GeV'])) < F(1, 10**38)


@pytest.mark.parametrize('masses', [[1, 1, 2], [0, 1, 2], [1, 1, 100],
                                  [1, 2, 1000], [True, 2, 3], [1.0, 2, 3], [1, 2]])
def test_invalid_or_unattainable_targets_are_not_promoted(masses):
    with pytest.raises(ValueError):
        c.fit(masses)


def test_mutating_export_cannot_select_a_source_or_change_a_later_certificate():
    before = c.certificate()
    altered = copy.deepcopy(before)
    altered['source_selected_controller'] = True
    altered['input_mass_energies_GeV'][2] = '2'
    assert c.certificate() == before
