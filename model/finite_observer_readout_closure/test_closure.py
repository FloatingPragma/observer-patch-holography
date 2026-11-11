from fractions import Fraction as F
import sympy as s
from . import closure as c


def test_full_original_periodic_record_law_retains_failures_and_timeout():
    result = c.periodic_y()
    assert result['full_detector_Kraus_identity_exact']
    for law in result['cases']:
        for experiment in law['finite_complete_records']:
            assert sum(F(row['probability']) for row in experiment['records']) == 1
            assert experiment['records'][0]['probability'] == '28/33'
            gap = F(experiment['original_infinite_event_plus_probability'])-F(experiment['complete_original_plus_probability'])
            assert 0 <= gap <= F(experiment['original_timeout_probability'])
        assert not law['harmonic']['native_repair_generates_modular_parameter']


def test_joint_clock_instrument_keeps_full_original_reference_output():
    result = c.joint_clocks()
    assert [r['complete_three_attempt_runtime']['complete_field_plus_probability'] for r in result['cases']] == ['35/54', '14/27']
    for row in result['cases']:
        assert row['complete_reference_feedback_runtime']['reference_dimension'] == 2
        assert row['complete_three_attempt_runtime']['complete_record_count'] == 250
    cycle = result['complete_relative_count_law']['renewal_cycles_without_feedback']
    assert cycle['mean_B_ticks_per_cycle'] == '2' and cycle['variance_B_ticks_per_cycle'] == '16'


def test_full_source_cyclic_witness_proves_more_than_a_finite_window_fit():
    result = c.cyclic_transfer()
    assert result['cyclic_dimension'] == 7
    assert result['full_source_intertwiner_KC_equals_CR']
    assert len(result['full_source_inclusion']) == 66
    assert [r['lag'] for r in result['independent_unused_lags']] == list(range(15, 33))
    for row in result['independent_unused_lags']:
        assert row['full_source_centered_correlation'] == row['record_only_prediction']
    gram = s.Matrix(result['record_Gram'])
    transfer = s.Matrix(result['record_transfer'])
    assert gram*transfer == transfer.T*gram
    assert not result['sampled_data_certifies_exact_rank']
