from fractions import Fraction as F

import pytest

from . import algebra, bounds


def test_exact_mass_four_and_mass_sixty_four_interfaces():
    first = bounds.spectral_window()
    assert first['one_body_gap_lower'] == 2
    assert first['origin_coupling_norm_squared'] == F(3, 2)
    assert first['origin_weight_lower'] == F(8, 11)
    assert first['ipr_lower'] == F(64, 121)
    following = bounds.spectral_window(mass=64)
    assert following['one_body_gap_lower'] == 62
    assert following['origin_weight_lower'] == F(7688, 7691)
    assert following['origin_particle_number_deviation_contrast_lower'] == F(7688, 7691)**2
    assert following['hopping_ratio_upper'] == F(1, 32)
    assert not following['fixed_mass_velocity_continuum_refinement_stays_in_domain']


def test_rational_hopping_budget_and_neumann_bounds():
    sharper = bounds.spectral_window(hopping_bound=F(7, 4))
    assert sharper['one_body_gap_lower'] == F(9, 4)
    assert sharper['origin_weight_lower'] == F(27, 35)
    for count in range(8):
        expected = F(3, 8)*F(1, 2)**(2*count)
        assert bounds.neumann_error_squared(count) == expected
        assert bounds.tail_probability_upper(count+1) == expected
        if count:
            assert bounds.finite_volume_projector_error_squared(count) == min(F(1), 4*expected)
    # This admitted parameter point lies below twice the hopping norm. The
    # Schur result needs m>k, without the stronger m>2k separation assumption.
    assert bounds.spectral_window(mass=3)['one_body_gap_lower'] == 1


def test_source_quaternion_and_origin_density_identities():
    value = algebra.certificate()
    assert value['Kramers_square'] == -1
    assert value['origin_coupling_Gram_in_units_v2_over_a2'] == '3/2 * I2'
    assert value['local_zero_frame_determinant_squared'] == 'rho_n^2'
    assert value['origin_D0_charged_minus_single_contrast'] == 'rho_0^2'
    assert value['origin_contrast_is_noninteracting_reference_identity'] is True


@pytest.mark.parametrize('kwargs', [
    {'spacing': 0}, {'spacing': 1.}, {'velocity': True}, {'mass': 2},
    {'mass': 1}, {'hopping_bound': F(3, 2)}, {'hopping_bound': 0},
    {'spacing': F(1, 4), 'mass': 4, 'hopping_bound': 7}])
def test_uncontrolled_or_nonexact_domains_are_rejected(kwargs):
    with pytest.raises(ValueError):
        bounds.spectral_window(**kwargs)


def test_invalid_distance_and_term_counts_are_rejected():
    for value in (-1, True, 1.):
        with pytest.raises(ValueError):
            bounds.neumann_error_squared(value)
    with pytest.raises(ValueError):
        bounds.tail_probability_upper(0)
    with pytest.raises(ValueError):
        bounds.finite_volume_projector_error_squared(0)
