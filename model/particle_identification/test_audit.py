"""Independent controls for the exact identification exclusions."""
from fractions import Fraction as F

import pytest
import sympy as s

from . import audit


def test_actual_cube_density_has_even_grading_and_zero_number_change():
    result = audit.grading_certificate()
    assert result['complete_number_sector_dimension'] == 8008
    assert result['source_parity'] == 1
    assert result['transition_number_change'] == 0
    assert result['reference_occupation_columns'] == 56


def test_odd_operator_cannot_create_same_parity_state_but_even_one_can():
    parity = s.diag(1, 1, -1)
    source, desired = s.eye(3)[:, 0], s.eye(3)[:, 1]
    a, b, c, d = s.symbols('a b c d')
    odd = s.Matrix([[0, 0, a], [0, 0, b], [c, d, 0]])
    assert parity*odd == -odd*parity
    assert desired.dot(odd*source) == 0
    even = desired*source.T
    assert parity*even == even*parity
    assert even*source == desired


def test_quartic_obstruction_does_not_depend_on_energy_origin_or_mass():
    result = audit.shell_certificate()
    assert result['Lorentz_shell_squared_residual_quartic'] == 'alpha**2/16'
    assert result['nonrelativistic_approximation_excluded'] is False


@pytest.mark.parametrize('alpha,beta,gamma,m,c,l', [
    (1, 0, 1, 1, 1, 1), (-3, 7, 2, 4, 2, 5), (F(1, 3), -9, 4, 0, 7, 2)])
def test_fourth_derivative_retains_obstruction_for_changed_dictionary(alpha,beta,gamma,m,c,l):
    r = s.Symbol('r')
    residual = (alpha*(r*r/4-gamma*gamma)+beta)**2-c*c*l*l*r*r-m*m*c**4
    assert s.diff(residual, r, 4) == s.Rational(3, 2)*alpha*alpha


@pytest.mark.parametrize('r', [F(), F(1, 10), F(1), F(2)])
def test_rational_nonrelativistic_gap_encloses_exact_radical(r):
    result = audit.low_momentum_gap(3, F(1, 4), r, 2)
    gap = 3+s.Rational(r*r/4)-s.sqrt(9+s.Rational(3*r*r/2))
    assert gap >= s.Rational(result['lower'])
    assert gap <= s.Rational(result['upper'])
    if r:
        assert result['lower'] > 0


def test_nonrelativistic_gap_is_quartic_at_small_momentum():
    r = s.Symbol('r')
    assert s.series(3+r*r/4-s.sqrt(9+3*r*r/2), r, 0, 6).removeO() == r**4/96


@pytest.mark.parametrize('args', [(0,1,1,1),(1,0,1,1),(1,1,-1,1),(1,1,2,1),(1.,1,1,1),(True,1,1,1)])
def test_bad_or_nonexact_gap_inputs_fail(args):
    with pytest.raises(ValueError):
        audit.low_momentum_gap(*args)


def test_same_number_gaps_ignore_scalar_and_chemical_potential():
    result = audit.neutral_number_shift_certificate()
    assert result['transformed_gap'] == 'E - R'
    assert result['scalar_or_common_number_chemical_potential_repairs_negative_gaps'] is False


def test_prior_calibration_and_exposed_targets_are_not_new_predictions():
    from .packet import comparison_ancestry
    result = comparison_ancestry()
    assert result['physical_species_assignments'] == 0
    assert result['previous_dimer_excluded_assignments'] == 672
    assert result['previous_dimer_exclusion_applied_to_new_cube'] is False


def test_malformed_supplied_receipt_is_rejected(monkeypatch, tmp_path):
    from . import packet
    path = tmp_path/'bad.json'
    path.write_text('{"all_pass": false}')
    monkeypatch.setattr(packet, 'run', lambda: {'all_pass': True})
    monkeypatch.setattr('sys.argv', ['particle_identification', '--check', str(path)])
    with pytest.raises((ValueError, AssertionError, SystemExit)):
        packet.main()
