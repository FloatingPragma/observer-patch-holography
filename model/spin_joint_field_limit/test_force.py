from fractions import Fraction as F

import pytest
import sympy as s

from model.spin_joint_field_limit import force


def test_current_filter_has_natural_scale_without_linear_closure():
    rows = [force.filtered_current(n, F(5), F(1), F(2)) for n in (2, 4, 8, 16)]
    for n, row in zip((2, 4, 8, 16), rows):
        assert row['filtered_current_norm_upper'] == F(1, n**4)*(5+F(2, n))
        assert row['selection_torque_norm_upper'] == F(2, n)
        assert row['p_minus_two_filtered_current_norm_upper'] == 5+F(2, n)
        assert row['same_H_B_ground_and_filter_required']
        assert not row['linear_force_residual_vanishes']
        assert not row['full_history_force_closure_claimed']
    assert all(a['filtered_current_norm_upper'] > b['filtered_current_norm_upper']
               for a, b in zip(rows, rows[1:]))


def test_actual_finite_ground_filter_and_complete_spectrum():
    row = force.finite_source()
    assert row['all_source_rows'] == 4
    assert row['selection_field'] == '1/32'
    assert sum(row['complete_source_spectrum'].values()) == 4
    assert row['complete_source_spectrum']['0'] == 2
    assert row['filter_rank'] == 2
    assert row['current_drift_identity_exact']
    assert row['selection_torque_nonzero']
    assert row['filtered_ground_current_identity_exact']
    assert not row['finite_level_identified_as_continuum_mass_shell']


def test_original_current_probability_and_complement_both_retained():
    row = force.finite_source()
    current = row['complete_current_and_comparison_records'][0]
    success, complement = map(s.sympify, current['complete_effect_probabilities'])
    c = s.sympify(current['source_norm_normalization'])
    original = s.sympify(row['original_filtered_current_squared_norm'])
    assert force._simple(success+complement) == 1
    assert force._quadratic_interval(success)[0] > 0
    assert force._quadratic_interval(complement)[0] > 0
    assert force._simple(c*c*success-original) == 0
    assert force._simple(c*c-original) != 0  # Renormalizing success loses source weight.


def test_all_pauli_records_reconstruct_both_current_and_residual_effects():
    for row in force.finite_source()['complete_current_and_comparison_records']:
        readout = row['complete_Pauli_score_readout']
        assert len(readout['records']) == 64
        assert len({x['setting'] for x in readout['records']}) == 16
        probabilities = [s.sympify(x['original_probability']) for x in readout['records']]
        scores = [s.sympify(x['score']) for x in readout['records']]
        assert any(p == 0 for p in probabilities)
        assert force._simple(sum(probabilities)) == 1
        assert force._simple(sum(p*x for p, x in zip(probabilities, scores))-
                             s.sympify(row['complete_effect_probabilities'][0])) == 0
        assert all(abs(score) <= s.sympify(readout['score_absolute_bound']) for score in scores)
        assert readout['zero_probability_records_retained']
        assert not readout['conditional_success_distribution_used']
        assert not row['linear_force_equation_claimed']


def test_wrong_pauli_score_and_deleted_complement_are_detectable():
    row = force.finite_source()['complete_current_and_comparison_records'][0]
    data = row['complete_Pauli_score_readout']['records']
    index = next(i for i, event in enumerate(data)
                 if s.sympify(event['original_probability']) != 0 and s.sympify(event['score']) != 0)
    selected = data[index]
    change = s.sympify(selected['original_probability'])*s.sympify(selected['score'])
    assert force._simple(change) != 0
    success = s.sympify(row['complete_effect_probabilities'][0])
    assert force._simple(success-1) != 0


def test_score_error_keeps_original_normalization_and_scaled_accuracy():
    row = force.score_error(F(1, 32), F(4), F(1, 1000), F(1, 4))
    assert row['original_squared_norm_error_upper'] == F(1, 128000)
    assert row['p_minus_four_squared_norm_error_upper'] == F(1, 500)
    assert not row['success_conditioning_used']
    with pytest.raises(ValueError):
        force.score_error(F(), F(4), F(1, 1000))
    with pytest.raises(ValueError):
        force.filtered_current(1, F(5), F(1), F(1))
    with pytest.raises(ValueError):
        force.filtered_current(2, F(-1), F(1), F(1))


def test_quadratic_enclosures_are_exact_and_outward():
    for value in (s.sqrt(257), 1/s.sqrt(257), s.Rational(1, 2)-8/s.sqrt(257)):
        low, high = force._quadratic_interval(value)
        assert low <= high
        assert high-low < F(1, 10**12)
    assert force._nonnegative(s.Rational(1, 2)-8/s.sqrt(257))
    assert not force._nonnegative(s.Rational(1, 2)-9/s.sqrt(257))
