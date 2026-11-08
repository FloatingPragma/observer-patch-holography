"""Source force controls and original-probability spectral enclosures."""
from fractions import Fraction as F

import numpy as np
import pytest
import sympy as s

from . import charge_escape as e
from . import composite


def test_exact_four_mode_force_has_sharp_full_space_bounds():
    row = e.star_algebra()
    assert row['force_anticommutator_eigenvalues'] == {'1': 10, '5': 4, '9': 2}
    assert row['force_squared_eigenvalues'] == {'0': 8, '1': 4, '5': 4}
    assert row['full_contact_equation_exact']


def test_literal_original_weyl_axes_and_both_flavours():
    rows = e.literal_weyl_controls()
    assert [row['axis'] for row in rows] == [0, 1, 2]
    assert all(row['neighbor_orbital_norm_squared'] == F(1, 2) for row in rows)
    assert all(row['both_original_flavours_checked'] for row in rows)


@pytest.mark.parametrize('alpha2,beta2,energy', [
    (16, 9, F(1, 2)), (16, F(27, 2), 1), (16, 18, F(1, 1024)),
    (1, 0, F(1, 3)), (F(1, 9), 2, F(1, 4)), (16, 9, 0), (16, 9, 4)])
def test_rational_enclosure_matches_independent_analytic_bound(alpha2, beta2, energy):
    upper = e.low_window_upper(alpha2, beta2, energy, bits=40)
    a2, b2, value = map(s.Rational, (alpha2, beta2, energy))
    exact = (s.sqrt(a2)*value+s.sqrt(b2)*s.sqrt(a2+b2-value**2))**2/(a2+b2)**2
    if value**2 >= a2:
        exact = s.Integer(1)
    assert float(s.N(s.Rational(upper)-exact, 35)) >= -1e-30
    assert float(s.N(s.Rational(upper)-exact, 35)) <= 2**-40+1e-30
    assert e.upper_admissible(alpha2, beta2, energy, upper)
    if upper and energy and energy*energy < alpha2:
        assert not e.upper_admissible(alpha2, beta2, energy, upper-F(1, 2**40))


def test_orthogonality_is_essential_and_the_zero_window_bound_is_sharp():
    row = e.response_countercontrol()
    assert row['original_probabilities'][1] == e.low_window_upper(16, 9, 0)
    assert not row['actual_Hubbard_source_realization_claimed']
    # Without <v,Y>=0, L=0 and Y=-alpha*v would violate the bound.
    assert e.low_window_upper(16, 16, 0) == F(1, 2) < 1


@pytest.mark.parametrize('threshold', [0., 1., 3.8, 4., 4.8])
def test_actual_default_ratio_dimer_preserves_original_charge_probabilities(threshold):
    data = composite.source((2, 1, 1), 1, 8)
    h = data['H']/data['scale']
    spectrum, basis = np.linalg.eigh(h)
    omega = basis[:, 0]
    a, b = data['a'][0], data['b'][0]
    hopping = h-np.diag(np.diag(h))
    force = b@hopping-hopping@b
    # Retarded ground vector inner products use both original charge records.
    force_norm2 = np.vdot(force@omega, force@omega)+np.vdot(force.T@omega, force.T@omega)
    cross = np.vdot(a@omega, force@omega)+np.vdot(force.T@omega, a.T@omega)
    assert abs(cross) < 2e-14
    assert force_norm2 <= 9/4+2e-14
    addition = abs(basis.T@a.T@omega)**2
    removal = abs(basis.T@a@omega)**2
    np.testing.assert_allclose([addition.sum(), removal.sum()], [.5, .5], atol=2e-14, rtol=0)
    chosen = spectrum-spectrum[0] <= threshold+1e-13
    observed = float((addition+removal)[chosen].sum())
    budget = e.local_budget(1, 8, 1, F(str(threshold)))
    assert observed <= float(budget['original_signed_low_window_probability_upper'])+2e-14
    assert abs(float(addition[chosen].sum()-removal[chosen].sum())) < 2e-14


def test_default_local_escape_is_uniform_and_keeps_the_half_weight_factor():
    limiting = e.local_budget(1, 8, 6, 0)
    assert limiting['original_signed_low_window_probability_upper'] == F(27, 59)
    assert limiting['original_signed_outside_probability_lower'] == F(32, 59)
    assert limiting['each_local_PH_charge_low_window_probability_upper'] == F(27, 118)
    assert limiting['each_local_PH_charge_outside_probability_lower'] == F(16, 59)
    values = [e.local_budget(n, 8*n, 6, 1)['original_signed_low_window_probability_upper']
              for n in (16, 64, 256)]
    assert F(27, 59) < values[2] < values[1] < values[0] < F(1, 2)


def test_isolated_site_force_and_contact_control():
    below = e.local_budget(1, 8, 0, 1)
    assert below['force_retarded_norm_squared_upper'] == 0
    assert below['original_signed_low_window_probability_upper'] == F(1, 16)
    assert e.local_budget(1, 8, 0, 0)['original_signed_low_window_probability_upper'] == 0


def test_cubic_symbol_packet_bound_keeps_its_state_premise_and_full_mass():
    row = e.cubic_packet_budget(1, 8, 0)
    assert row['all_packet_force_squared_upper'] == 18
    assert row['original_signed_low_window_probability_upper'] == F(9, 17)
    assert row['original_signed_outside_probability_lower'] == F(8, 17)
    assert row['cubic_ground_symmetry_and_range_two_force_symbol_required']
    assert not row['individual_smeared_charge_has_half_weight']
    for g0, g2 in ((F(27, 2), -F(3, 4)), (F(27, 2), F(1, 4)), (6, -F(1, 2))):
        assert g0 >= 6*abs(g2)
        for cosine_sum in (-3, -2, 0, 1, 3):
            assert 0 <= g0+2*g2*cosine_sum <= row['all_packet_force_squared_upper']


def test_complete_time_average_lower_bound_and_nonzero_small_time_limit():
    for delta in (F(1, 4), F(1), F(4)):
        rows = [e.time_average_budget(c, 8*c, delta) for c in (64, 256, 1024)]
        values = [row['one_minus_real_time_average_lower'] for row in rows]
        assert 0 < values[0] < values[1] < values[2] < F(4, 17)
        assert all(row['complete_signed_response_mass'] == 1 for row in rows)
        assert all(not row['filtered_response_excluded'] for row in rows)
    # Independently integrate a complete symmetric measure, including a zero atom.
    energies = np.array([-17., 0., 17.])
    probabilities = np.array([4/17, 9/17, 4/17])
    delta = .25
    deficit = 1-np.dot(probabilities, np.sinc(delta*energies/np.pi))
    tail = probabilities[abs(energies) > 2/delta].sum()
    assert deficit >= tail/2
    assert abs(probabilities.sum()-1) < 1e-15


@pytest.mark.parametrize('args', [(0, 1, 6, 1), (1, 0, 6, 1), (1, 8, -1, 1),
                                 (1, 8, True, 1), (1., 8, 6, 1), (1, 8, 6, -1)])
def test_invalid_source_parameters_are_rejected(args):
    with pytest.raises(ValueError):
        e.local_budget(*args)


def test_invalid_projection_parameters_are_rejected():
    for call in (lambda: e.low_window_upper(0, 1, 1),
                 lambda: e.low_window_upper(1, -1, 0),
                 lambda: e.low_window_upper(1, 1, 0, bits=True),
                 lambda: e.upper_admissible(1, 1, 0, 2)):
        with pytest.raises(ValueError):
            call()
