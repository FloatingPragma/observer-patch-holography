"""Independent coherent-history, full-native and promised-input controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from math import comb

import pytest
import sympy as s

from ..radial_geometry_join import source
from . import controls, reuse


def test_exact_ordered_kraus_products_retain_quantum_cross_terms():
    operators = (s.diag(7, 1)/s.sqrt(50), s.I*s.diag(1, 7)/s.sqrt(50))
    rho = s.Matrix([[s.Rational(1, 3), s.I/5], [-s.I/5, s.Rational(2, 3)]])
    assert rho.det() > 0
    for rounds in range(5):
        unread = s.zeros(2)
        for history in product((0, 1), repeat=rounds):
            operator = s.eye(2)
            for bit in history:
                operator = operators[bit]*operator
            result = operator*rho*operator.H
            row = reuse.history_coefficients(rounds, sum(history))
            assert result[0, 0] == row['parity_zero']/3
            assert result[1, 1] == 2*row['parity_one']/3
            assert result[0, 1] == s.I*row['cross_parity']/5
            assert result[0, 1] != 0
            unread += result
        assert unread[0, 0] == rho[0, 0] and unread[1, 1] == rho[1, 1]
        assert unread[0, 1] == rho[0, 1]*s.Rational(7, 25)**rounds


@pytest.mark.parametrize('arguments', [(-1, 0), (True, 0), (2.0, 0), (1025, 0), (2, 3), (2, True)])
def test_invalid_complete_history_domain_rejected(arguments):
    with pytest.raises(ValueError):
        reuse.history_coefficients(*arguments)


def test_full_reference_conditional_states_are_not_initially_dephased():
    row = controls.quantum_histories()
    assert row['ordered_histories_checked'] == 62
    assert row['conditional_cross_parity_block_is_nonzero']
    assert not row['input_dephased_before_readout']


def test_all_native_states_and_exact_correlated_reuse_countercontrol():
    result = reuse.native_transition_control()
    assert result['complete_one_round_native_histories_checked'] == 79872
    assert result['full_78_state_two_round_mean_checked']
    assert result['two_round_reused_minus_fresh_matrix_00'] > 0
    q = F(3, 7)
    mixture = reuse.reused_mean(q, 2)
    wrong = reuse.matrix_power(source.mean_matrix(q), 2)
    delta = [[b-a for a, b in zip(row0, row1, strict=True)]
             for row0, row1 in zip(source.mean_matrix(F(0)), source.mean_matrix(F(1)), strict=True)]
    square = reuse.matmul(delta, delta)
    assert all(mixture[i][j]-wrong[i][j] == q*(1-q)*square[i][j]
               for i in range(12) for j in range(12))
    assert reuse.reused_mean(q, 1) == source.mean_matrix(q)


def test_deleted_native_current_record_is_rejected(monkeypatch):
    original = reuse.native.native
    def omitted_current(state, edge, coin):
        output, (_, _) = original(state, edge, coin)
        return output, (coin, 0)
    monkeypatch.setattr(reuse.native, 'native', omitted_current)
    with pytest.raises(ArithmeticError, match='reversible'):
        reuse.native_transition_control()


def test_classifier_tie_and_complete_failure_scoring_without_record_loss():
    assert reuse.classify([1]*16+[0]*16) == 'charged'
    assert reuse.classify([1]*15+[0]*17) == 'ground'
    frames = [{'kind': 'pointer', 'value': 1} for _ in range(20)]
    frames += [{'kind': 'failure', 'native_label': 'round20/pulse4/complement2'}]
    frames += [{'kind': 'halted'} for _ in range(11)]
    retained = deepcopy(frames)
    scores = reuse.score_records(frames)
    assert frames == retained and len(frames) == len(scores) == 32
    assert scores == (1,)*20+(0,)*12
    assert reuse.classify(scores) == 'charged'
    for reserved in (2, 3, 4, 5):
        frames[20] = {'kind': 'pointer', 'value': reserved}
        assert reuse.score_records(frames) == scores
    with pytest.raises(ValueError, match='32'):
        reuse.score_records(frames[:21])
    frames[21] = {'kind': 'pointer', 'value': 1}
    with pytest.raises(ValueError, match='halted'):
        reuse.score_records(frames)


def test_exact_binomial_tails_and_single_preparation_error():
    row = reuse.classification_bound()
    independent_ground = sum((F(comb(32, k)*49**(32-k), 50**32) for k in range(16, 33)), F())
    independent_charged = sum((F(comb(32, k)*49**k, 50**32) for k in range(16)), F())
    assert independent_ground == row['exact_ground_parity_conditional_error']
    assert independent_charged == row['exact_charged_parity_conditional_error']
    assert independent_charged < independent_ground < F(1, 10**6)
    assert row['each_spectral_input_misclassification_strict_upper'] < F(1, 2500)
    assert row['field_preparations'] == 1
    assert row['full_32_round_half_trace_error_upper'] < F(32, 10**9)


def test_one_field_record_average_keeps_nonvanishing_ensemble_variance():
    q, rounds = F(1, 3), 32
    p0, p1 = F(1, 50), F(49, 50)
    mean = (1-q)*p0+q*p1
    second = F()
    for k in range(rounds+1):
        weight = comb(rounds, k)*reuse.history_probability(q, rounds, k)
        second += weight*F(k, rounds)**2
    assert second-mean**2 == p0*(1-p0)/rounds+q*(1-q)*(p1-p0)**2
    assert second-mean**2 > q*(1-q)*(p1-p0)**2


def test_complete_fock_work_injection_includes_leaving_ground_support():
    row = controls.full_fock_work()
    assert row['complete_Fock_rows'] == 64
    assert row['actual_postmeasurement_ground_support_probability_below_one']
    assert not row['ground_sector_protection_assumed']
    assert reuse.unread_coefficients(0)['origin_hopping_energy_injection_coefficient'] == 0


def test_finite_resources_include_fresh_apparatus_but_only_one_field():
    row = reuse.resources()
    assert row['native_parity_pulses_per_round'] == 40
    assert row['native_pulses_total_upper'] == 2432
    assert row['native_integer_repair_attempts_upper'] == 32
    assert row['fresh_apparatus_logical_bits_without_recycling'] == 736
    assert row['field_preparations'] == 1
    assert row['retained_classical_record_frames'] == 32
    assert not row['complete_32_round_native_word_executed']
