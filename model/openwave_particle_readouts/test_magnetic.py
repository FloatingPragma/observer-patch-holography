"""Full quantum-link controls for the same-source static response theorem."""
import pytest
import sympy as s

from model.coherent_defect_wave.chiral_joint.periodic_family.model import (
    Configuration, kinetic_terms,
)
from model.coherent_defect_wave.chiral_joint.oriented_links.model import empty_key
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.model import (
    hopping_column,
)
from .magnetic import (
    COUPLING, PROFILE, background_hopping_column, coordinate_checks,
    fixed_preparation_control, full_source_checks, generator_value, phase,
    plaquette_checks, source_column, validate_profile,
)


def test_complete_PW_coordinates_include_excited_link_transitions():
    row = coordinate_checks()
    assert row['native_coordinate_matrices'] == 59
    assert row['nonzero_coordinate_entries'] == 1156
    assert row['occupied_to_occupied_entries'] == 1038
    assert row['central_endpoint_commutators'] == 24
    assert row['noncentral_weak_direction_transfer_rejected']


def test_actual_interacting_columns_and_both_current_derivatives():
    row = full_source_checks()
    assert row['complete_240_mode_basis_columns'] == 15
    assert row['nonzero_current_entries'] > 0
    assert row['nonzero_diamagnetic_entries'] > 0
    assert row['actual_hopping_coefficient'] == str(COUPLING)


@pytest.mark.parametrize('side', [2, 4])
def test_independent_path_phase_and_current_on_full_source_carrier(side):
    config = Configuration(side=side)
    base = empty_key(config)
    # At N=4, every axial spectral path length occurs. All native families
    # are occupied in this literal column; no free link value is substituted.
    mask = sum(1 << (30*(side**3-1)+a) for a in (0, 2, 9, 10, 12))
    key = mask, base[1], base[2]
    profile = tuple((edge % 3)-1 for edge in range(3*side**3))
    ordinary = hopping_column(key, kinetic_terms(config))
    changed = background_hopping_column(key, profile, config=config)
    expected = {target: phase(target, profile)*s.conjugate(phase(key, profile))*value
                for target, value in ordinary.items()}
    assert changed == expected
    derivative = background_hopping_column(key, profile, derivative=1,
                                           at_test_angle=False, config=config)
    expected_derivative = {target: s.I*(generator_value(target, profile)
                            - generator_value(key, profile))*value
                           for target, value in ordinary.items()
                           if generator_value(target, profile) != generator_value(key, profile)}
    assert derivative == expected_derivative


def test_nonzero_flux_cannot_be_dismissed_as_site_gauge_gradient():
    row = plaquette_checks()
    assert row['signed_profile_circulation'] == 1
    assert row['ec_Wilson_phase'] == '-1'
    assert row['added_plaquette_conjugation_difference_on_blank_norm_squared'] == '8'
    assert not row['bare_plaquette_term_present_in_actual_H']


def test_fixed_physical_preparation_is_distinct_from_paired_apparatus():
    row = fixed_preparation_control()
    assert s.Rational(row['fixed_state_energy_mean_change']) != 0
    assert row['paired_state_energy_mean_change'] == '0'
    assert not row['all_fixed_preparation_magnetic_responses_inferred_zero']


def test_zero_angle_recovers_actual_source():
    base = empty_key()
    key = 1 << (120+15+9), base[1], base[2]
    ordinary = hopping_column(key, kinetic_terms())
    assert background_hopping_column(key, at_test_angle=False) == ordinary
    assert source_column(key, 0) == source_column(key, 0, background=True)


def test_wrong_family_charge_breaks_full_link_current_covariance(monkeypatch):
    from . import magnetic
    source_representation = magnetic.representation()
    base = empty_key()
    key = 1 << (120+15+10), base[1], base[2]
    ordinary = hopping_column(key, kinetic_terms())
    expected = {target: s.I*(generator_value(target)-generator_value(key))*value
                for target, value in ordinary.items()
                if generator_value(target) != generator_value(key)}
    wrong_currents = (*source_representation[2][:-1], s.eye(15))
    wrong_representation = (*source_representation[:2], wrong_currents,
                            *source_representation[3:])
    monkeypatch.setattr(magnetic, 'representation', lambda: wrong_representation)
    wrong = background_hopping_column(key, derivative=1, at_test_angle=False)
    assert expected and wrong != expected


def test_invalid_backgrounds_are_rejected():
    with pytest.raises(ValueError):
        validate_profile(PROFILE[:-1])
    with pytest.raises(ValueError):
        validate_profile((s.Rational(1, 2),)+PROFILE[1:])
    with pytest.raises(ValueError):
        background_hopping_column(empty_key(), derivative=3)
    with pytest.raises(ValueError):
        source_column(empty_key(), s.Float(.1))


def test_degenerate_spectral_current_and_complete_static_cancellation():
    # An independent exact carrier checks the degenerate Ward identity and
    # exposes the false splitting obtained by dropping the diamagnetic term.
    h = s.diag(0, 0, 1, 3, 4)
    k = s.Matrix([[0, 1, 2, 0, 1], [1, 2, 1, 3, 0],
                  [2, 1, 0, 1, 2], [0, 3, 1, 1, 1], [1, 0, 2, 1, 0]])
    current = s.I*(k*h-h*k)
    double = -(k*(k*h-h*k)-(k*h-h*k)*k)
    for energy in (0, 1, 3, 4):
        selected = s.diag(*(int(h[i, i] == energy) for i in range(5)))
        resolvent = s.diag(*(0 if h[i, i] == energy else 1/(energy-h[i, i])
                              for i in range(5)))
        assert selected*current*selected == s.zeros(5)
        paramagnetic = 2*selected*current*resolvent*current*selected
        assert paramagnetic != s.zeros(5)
        assert selected*double*selected+paramagnetic == s.zeros(5)
