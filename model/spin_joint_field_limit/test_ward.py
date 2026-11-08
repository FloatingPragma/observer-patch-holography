"""Source signs, Ward saturation, polarization mismatch and complete records."""
from fractions import Fraction as F

import pytest
import sympy as s

from . import ward
from .force import _quadratic_interval


def test_actual_joint_source_uniform_ward_closes_both_finite_equations():
    result = ward.finite_source()
    row = result['profiles'][0]
    u, b = map(s.sympify, row['source_coefficients_u_b'])
    energy = s.sympify(row['one_retained_positive_energy'])
    assert result['field'] == '1/32' and result['p'] == '1/4'
    assert result['uniform_spin_Ward_susceptibility_exact']
    assert row['positive_form_slacks'] == ['0', '0']
    assert b == s.Rational(1, 2)
    assert s.simplify(u*b-energy**2) == 0
    assert not result['continuum_particle_or_mass_inferred']


def test_actual_uniform_variance_uses_the_complete_unscaled_source_measure():
    row = ward.finite_source()['uniform_source_weighted_energy_variance']
    h, chi, chi3 = map(s.sympify, (row['first_moment'], row['inverse_moment'],
                                  row['third_inverse_moment']))
    energy = (16+s.sqrt(257))/32
    weight = s.sympify(row['original_preparation_weight'])
    assert s.simplify(h-energy*weight) == 0
    assert s.simplify(chi-weight/energy) == 0
    assert s.simplify(chi3-weight/energy**3) == 0
    assert s.simplify(chi/chi3-energy**2) == 0
    assert row['Q_force_form_slack'] == row['P_force_form_slack'] == '0'
    assert row['weighted_energy_squared_variance'] == '0'
    assert not row['auxiliary_weight_changes_original_event_probabilities']


def test_uniform_ward_equation_allows_a_genuine_two_energy_measure():
    k = s.diag(1, 2)
    q = s.Matrix([1, 1])/s.sqrt(2)
    b = s.Rational(1, 3)
    a = s.I*b*k.inv()*q
    assert k*a-s.I*b*q == s.zeros(2, 1)
    chi = (q.H*k.inv()*q)[0]
    chi3 = (q.H*k.inv()**3*q)[0]
    h = (q.H*k*q)[0]
    exact = lambda x: F(int(x.p), int(x.q))
    row = ward.uniform_variance_slack(exact(h), exact(chi), exact(chi3), exact(b))
    assert row['Q_force_form_slack'] == F(1, 2)
    assert row['P_force_form_slack'] == 0
    assert row['weighted_energy_squared_variance'] == F(8, 9)
    assert row['weighted_energy_squared_mean'] == F(4, 3)
    assert row['source_cross_moment'] == F(1, 4)
    u = s.Rational(row['source_coefficients_u_b'][0])
    # Independent inverse-form residual squared, without a chosen square root.
    residual = s.I*k*q-u*a
    assert s.simplify((residual.H*k.inv()*residual)[0]-s.Rational(1, 2)) == 0
    assert not row['continuum_variance_vanishing_inferred']
    smaller = ward.uniform_variance_slack(exact(h/4), exact(chi/4), exact(chi3/4), exact(b))
    assert smaller['Q_force_form_slack'] == row['Q_force_form_slack']/4
    assert smaller['source_cross_moment'] == row['source_cross_moment']/4
    assert smaller['weighted_energy_squared_variance'] == row['weighted_energy_squared_variance']
    with pytest.raises(ValueError, match='variance'):
        ward.uniform_variance_slack(1, 2, 1, 1)
    with pytest.raises(ValueError):
        ward.uniform_variance_slack(1, 1, 0, 1)


def test_positive_force_slack_does_not_exclude_an_actual_sharp_line():
    row = ward.finite_source()['profiles'][1]
    assert row['profile'] == [1, 2]
    assert all(_quadratic_interval(s.sympify(value))[0] > F(1, 100)
               for value in row['positive_form_slacks'])
    assert row['all_original_energy_weight_on_that_line']
    assert not row['positive_pair_slack_excludes_a_sharp_component']


def test_actual_charge_filter_resolves_both_positive_weight_force_pairs():
    row = ward.finite_source()['profiles'][1]
    branches = row['actual_charge_resolved_force_pairs']
    assert [x['internal_source_charge'] for x in branches] == [-1, 1]
    for branch in branches:
        assert branch['force_form_slacks'] == ['0', '0']
        assert _quadratic_interval(s.sympify(branch['original_Q_branch_probability']))[0] > F(1, 25)
        assert _quadratic_interval(s.sympify(branch['original_P_branch_probability']))[0] > F(1, 25)
        u, b = map(s.sympify, branch['source_coefficients_u_b'])
        energy = s.sympify(branch['source_excitation_energy'])
        assert s.simplify(u*b-energy**2) == 0
    assert _quadratic_interval(s.sympify(branches[0]['filtered_cross_moment']))[1] < 0
    assert _quadratic_interval(s.sympify(branches[1]['filtered_cross_moment']))[0] > 0
    assert s.simplify(sum(s.sympify(x['filtered_cross_moment']) for x in branches)
                      -s.sympify(row['order_cross_moment'])) == 0
    records = row['complete_preparation_charge_records']['records']
    measured = row['complete_preparation_charge_records']['original_cross_moment_by_retained_outcome']
    for branch in branches:
        assert measured[str(branch['internal_source_charge'])] == branch['filtered_cross_moment']
    assert len(records) == 16
    assert sum(s.sympify(x['original_probability']) for x in records) == 1
    assert sum(x['outcome'] == 'charge_complement' for x in records) == 4


def test_complete_source_settings_keep_every_failure_and_energy_label():
    for row in ward.finite_source()['profiles']:
        record = row['complete_preparation_energy_records']
        assert record['complete_record_count'] == 16
        assert len(record['records']) == 16
        assert sum(s.sympify(x['original_probability']) for x in record['records']) == 1
        assert sum(x['outcome'] == 'preparation_failure' for x in record['records']) == 4
        assert any(s.sympify(x['original_probability']) == 0 for x in record['records'])
        assert record['original_interference_cross_moment'] == row['order_cross_moment']
        assert not record['preparation_success_renormalized']
        assert not record['additional_native_tape_executed']


def test_source_matrix_mutation_is_rejected(monkeypatch):
    original = ward.source.matrix
    def changed(*args, **kwargs):
        matrix, denominator = original(*args, **kwargs)
        matrix = matrix.copy()
        matrix[0, 0] += 1
        return matrix, denominator
    ward.finite_source.cache_clear()
    try:
        monkeypatch.setattr(ward.source, 'matrix', changed)
        with pytest.raises(ArithmeticError, match='source frame'):
            ward.finite_source()
    finally:
        ward.finite_source.cache_clear()


def test_slack_intervals_keep_source_response_floor_and_signed_cross_moment():
    result = ward.interval_slack((F(3, 2), F(151, 100)),
        (F(49, 100), F(1, 2)), (F(1, 5), F(21, 100)))
    assert result['force_form_slack_interval'][0] == F(1, 4)
    assert result['source_coefficient_interval'] == [F(7, 3), F(5, 2)]
    neg = ward.interval_slack((2, 2), (-1, -1), (1, 1))
    assert neg['source_coefficient_interval'] == [-1, -1]
    assert not result['empirical_coefficient_fit_used']


@pytest.mark.parametrize('susceptibility', [(0, 1), (-1, 1), (2, 1)])
def test_susceptibility_inversion_rejects_uncertified_floor(susceptibility):
    with pytest.raises(ValueError):
        ward.interval_slack((1, 2), (0, 1), susceptibility)


def test_inconsistent_positive_slack_data_are_rejected():
    with pytest.raises(ValueError, match='contradict'):
        ward.interval_slack((1, 1), (2, 2), (1, 1))


def test_tube_bound_retains_original_preparation_probability():
    exact = ward.tube_weight(3, 2, 2, 0, 0, F(1, 2), F(1, 4), 2)
    assert exact['original_tube_probability_lower'] == F(1, 16)
    assert not exact['success_probability_divided_out']
    bounded = ward.tube_weight(3, 2, 2, F(1, 10000), F(1, 10000),
                               F(1, 2), F(1, 4), 1)
    assert bounded['original_tube_probability_lower'] == F(1211, 5000)
    negative = ward.tube_weight(3, -2, -2, 0, 0, F(1, 2), F(1, 4), 2)
    assert negative['original_tube_probability_lower'] == F(1, 16)


def test_tube_bound_against_independent_two_energy_vector_calculation():
    # A nearly monochromatic pair with a small genuine second spectral weight.
    k = s.diag(1, 2)
    q = s.Matrix([s.Rational(1, 2), s.Rational(1, 100)])
    a = s.I*q
    nq = (q.H*q)[0]
    c = s.im((q.H*a)[0])
    chi = (q.H*k.inv()*q)[0]
    h = (q.H*k*q)[0]
    u = c/chi
    delta = h-c*c/chi
    as_fraction = lambda x: F(int(x.p), int(x.q))
    result = ward.tube_weight(2, as_fraction(u), as_fraction(u),
        as_fraction(delta), as_fraction(delta), F(1, 4), as_fraction(nq), 1)
    assert abs(float(u)-1) < 0.001
    outside = F(1, 10000)
    assert result['outside_tube_source_norm_squared_upper'] >= outside
    assert 0 < result['original_tube_probability_lower'] <= F(1, 4)


def test_complete_score_error_charges_setting_and_operator_normalization():
    assert ward.moment_error(2, 5, F(1, 10000), F(1, 4)) == F(2, 125)
    with pytest.raises(ValueError):
        ward.moment_error(1, 1, F(1, 100), 2)
    with pytest.raises(TypeError):
        ward.moment_error(1.0, 1, F(1, 100))


def test_owner_replay_rejects_changed_ward_value_and_scope():
    from copy import deepcopy
    from . import packet
    receipt = packet.build()
    changed = deepcopy(receipt)
    changed['source_Ward_susceptibilities_and_force_equations']['actual_source'][
        'profiles'][0]['source_coefficients_u_b'][1] = '1'
    with pytest.raises(ValueError, match='differs'):
        packet.validate(changed)
    changed = deepcopy(receipt)
    changed['scope']['positive_full_pair_slack_excludes_sharp_subbranches'] = True
    with pytest.raises(ValueError, match='differs'):
        packet.validate(changed)
