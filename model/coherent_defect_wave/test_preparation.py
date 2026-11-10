"""Same-clock preparation and every recorded detector branch."""
from fractions import Fraction as F

import sympy as s

from .model import JointObserver, combine
from .preparation import certificate, executed_history


def same_vector(left, right):
    assert left.norm_denominator == right.norm_denominator
    a, b = left.as_dict(), right.as_dict()
    return all(s.simplify(a.get(k, 0)-b.get(k, 0)) == 0 for k in set(a) | set(b))


def test_filter_retains_both_branches_for_different_owned_inputs():
    observer = JointObserver()
    for port, expected_success in ((0, F(5, 9)), (6, F(17, 135))):
        initial = observer.basis_state(observer.initial_word, port)
        branches = observer.preparation_filter(initial)
        assert tuple(record['outcome'] for record, _ in branches) == ('success', 'failure')
        probabilities = tuple(state.norm_squared() for _, state in branches)
        assert probabilities == (expected_success, 1-expected_success)
        for _, state in branches:
            assert same_vector(observer.inverse_step(observer.step(state)), state)


def test_reference_numerator_commutes_with_actual_joint_tick():
    observer = JointObserver()
    initial = observer.basis_state(observer.initial_word, 0)

    def numerator(state):
        forwards = observer.step(state)
        backwards = observer.inverse_step(state)
        return combine(state, combine(forwards, backwards, F(1, 2), F(1, 2)), F(5, 4), 1)

    assert same_vector(observer.step(numerator(initial)), numerator(observer.step(initial)))


def test_complete_preparation_and_detector_history_has_no_dropped_failure():
    history = executed_history()
    assert history['input_is_uniform_ensemble'] is False
    rows = history['outcomes']
    assert len(rows) == 4
    assert [row['joint_probability'] for row in rows] == [
        '2267953/14580000', '5832047/14580000',
        '1132603/7290000', '2107397/7290000']
    assert sum(F(row['joint_probability']) for row in rows) == 1
    assert sum(F(row['joint_probability']) for row in rows
               if row['preparation'] == 'failure') == F(4, 9)


def test_trace_normalization_uses_all_link_words_and_readout_domains():
    result = certificate()
    assert result['joint_dimension'] == 72000
    assert result['joint_unitary_trace'] == '-47850'
    assert result['success_probability_from_uniform_basis_ensemble'] == '281/1080'
    assert result['reference_eigenvalue_lower'] == '1/168600'
    assert result['scope']['joint_archive_stationarity'] is False
    assert result['scope']['arbitrary_effect_tomography'] is False
