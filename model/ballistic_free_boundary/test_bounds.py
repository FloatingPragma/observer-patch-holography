from copy import deepcopy
from fractions import Fraction as F
import pytest
from . import bounds, packet


def test_constants_keep_spin_and_dimension_factors():
    k = bounds.constants()
    assert 2*k['weighted_row_square']*k['lattice_weight_sum'] < k['schur_coefficient']**2
    assert k['smoothed_error_coefficient'] == 2**17


@pytest.mark.parametrize('alpha', [F(-1), F(0), F(1, 5), F(1, 4), F(2, 7)-F(1, 10000)])
def test_strict_power_class_controls_all_three_errors(alpha):
    exponents, admitted = bounds.power_exponents(alpha)
    assert admitted and all(p < 0 for p in exponents.values())


@pytest.mark.parametrize('alpha', [F(2, 7), F(1, 3), F(1)])
def test_boundary_does_not_masquerade_as_interacting_or_free(alpha):
    exponents, admitted = bounds.power_exponents(alpha)
    assert not admitted
    assert exponents['weighted_composite'] >= 0


def test_zero_coupling_has_exact_free_dynamics():
    assert bounds.field_error(F(2), F(0), F(1), F(1, 2)) == 0


def test_bound_caps_weight_without_using_invalid_exponential_majorant():
    row = bounds.spectral_bound(F(1), F(1), F(100), F(1, 2))
    assert row['exponential_argument'] >= 1
    assert row['weight_upper'] == 1


def test_growing_power_examples_control_the_weighted_force():
    rows = [bounds.example(b) for b in (2, 8, 16, 64)]
    weights = [F(row['coupling_squared_weight_upper']) for row in rows]
    assert weights[0] > 1
    assert weights[1] < F(1, 100)
    assert weights[2] < F(1, 30_000)
    assert weights[3] < F(1, 2_000_000_000)
    assert all(a > b for a, b in zip(weights, weights[1:]))
    assert F(rows[-1]['unit_time_unit_window_field_error_upper']) < F(1, 90_000)


@pytest.mark.parametrize('value', [F(0), F(1, 3), F(2), F(10**20+1, 10**30)])
def test_square_root_is_an_outward_bound(value):
    result = bounds.sqrt_upper(value)
    assert result**2 >= value
    assert result == 0 or (result-F(1, 2**80))**2 <= value


@pytest.mark.parametrize('args', [(0, 1, 1, F(1, 2)), (1, -1, 1, F(1, 2)),
                                  (1, 1, -1, F(1, 2)), (1, 1, 1, 2),
                                  (True, 1, 1, F(1, 2))])
def test_invalid_spectral_parameters_rejected(args):
    with pytest.raises(ValueError):
        bounds.spectral_bound(*args)


def test_scope_and_arithmetic_cannot_be_promoted():
    value = packet.build()
    packet.validate(value)
    for field in ('physical_species_identified', 'interactions_above_boundary_proved',
                  'sharp_boundary_at_power_2_over_7_claimed', 'protected_spin_limit_excluded'):
        altered = deepcopy(value)
        altered['theorem'][field] = True
        with pytest.raises(ValueError):
            packet.validate(altered)
    altered = deepcopy(value)
    altered['constants']['weighted_row_square'] = 8
    with pytest.raises(ValueError):
        packet.validate(altered)


def test_boolean_alias_cannot_forge_scope():
    value = packet.build()
    value['theorem']['physical_species_identified'] = 0
    with pytest.raises(ValueError):
        packet.validate(value)
