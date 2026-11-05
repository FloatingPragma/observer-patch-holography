"""Independent fine-record sums, complete detector laws and exact resource controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
from scipy import linalg as la

from ..common_source import native_family_record_memory as memory
from ..common_source import native_repair_mass as repair
from ..compact_interacting_response.algebra import compact_memory_bound
from . import bounds as b, experiment as e, packet


@pytest.fixture(scope='module')
def receipt():
    return packet.produce()


def test_common_source_budget_matches_independent_parent_formula():
    for steps in (1, 17, 2**30):
        expected = compact_memory_bound(1, F(1, 128), F(1, 12), 1, steps)
        actual = b.segment(b.COMPACT_NORM, 1, F(1, 12), steps)
        assert F(expected['core_Fock_norm_upper']) == b.COMPACT_NORM
        assert F(expected['unread_arbitrary_initial_CQ_error_upper']) == min(F(2), actual['unread_trace_norm_upper'])
        assert min(F(4), F(expected['complete_history_arbitrary_CQ_error_squared_upper'])) == actual['complete_record_trace_norm_upper_squared']


def test_exact_schedule_allocates_complete_history_error_without_restarting_memory():
    result = b.schedule((F(1, 3),)*3)
    assert result['total_native_steps'] == 46370150349
    for row in result['intervals']:
        assert F(row['complete_record_error_squared_upper']) <= F(1, 300)**2
        assert row['stationary_product_used'] is False
        previous = b.segment(b.COMPACT_NORM, F(1, 3), F(1, 12), row['native_steps']-1)
        assert previous['complete_record_trace_norm_upper_squared'] > F(1, 300)**2
    prepared = b.schedule((F(1, 3),)*3, first_stationary_product=True)
    assert [row['stationary_product_used'] for row in prepared['intervals']] == [True, False, False]
    assert prepared['intervals'][0]['native_steps'] < result['intervals'][0]['native_steps']
    assert prepared['intervals'][1:] == result['intervals'][1:]


def test_equal_grid_formula_and_growing_intervention_boundary():
    m, t, p, count, n = F(3), F(2, 7), F(1, 9), 5, 1000000
    direct = 8*m*t*count**2/(n*(1-p))+2*(8/(1-p)+4)*m*m*t*t*count/n
    assert b.equal_grid_squared(m, t, p, count, n) == direct < 4
    assert b.equal_grid_squared(m, t, p, 100, 100) == 4


def test_all_native_fine_labels_and_exact_classical_refinement():
    data = e.law()
    for old, rows in enumerate(data['classes']):
        assert sorted(label for row in rows for label, _ in row['labels']) == list(range(512))
        assert sum(weight == 0 for row in rows for _, weight in row['labels']) == 416
        assert sum(row['weight'] for row in rows) == 1
        for row in rows:
            assert sum(weight for _, weight in row['labels']) == row['weight']
            if row['weight']:
                assert sum(weight/row['weight'] for _, weight in row['labels']) == 1
    assert e.law_certificate()['conditional_reset_port_next_step_vs_stationary_TV'] == '56025/67055104'


def test_grouped_quantum_step_equals_literal_512_branch_sum():
    operators = e.probe()
    initial, _ = e.initial_vectors()
    initial = initial[:1]
    duration = F(1, 64)
    actual, _, labels, _, _ = e.advance(initial, initial, np.array([2]), np.array([False]),
                                       np.array([0]), float(duration), operators)
    compressed = [np.zeros((32, 32), complex) for _ in range(12)]
    for vector, label in zip(actual, labels, strict=True):
        compressed[label] += np.outer(vector, vector.conj())
    literal = [np.zeros((32, 32), complex) for _ in range(12)]
    _, (gc, ge) = memory.geometry()
    off = la.expm(-1j*float(duration)*operators['fixed'])
    on = [la.expm(-1j*float(duration)*h) for h in operators['pulses']]
    count = 0
    for cc, cr, ec, er, _, raw in repair.native_records():
        for left, right in product((False, True), repeat=2):
            active = cc == 0 and cr == 1 and left
            event_e = ec == 0 and er in (4, 5) and right
            new = memory._destination(2, active, event_e, 2, gc, ge)
            weight = float(raw)*(3/8 if left else 5/8)*(3/8 if right else 5/8)
            vector = e.apply(on[new] if active else off, initial)[0]
            literal[new] += weight*np.outer(vector, vector.conj())
            count += 1
    assert count == 512
    assert max(la.norm(a-b) for a, b in zip(literal, compressed, strict=True)) < 1e-14


def test_complete_instruments_preserve_memory_and_all_full_carrier_inputs():
    operators = e.probe()
    assert la.norm(operators['middle'][0]@operators['final'][0]
                   -operators['final'][0]@operators['middle'][0]) > .1
    initial = np.eye(32, dtype=complex)
    labels = np.arange(32, dtype=np.int64) % 12
    active = np.arange(32) % 2 == 0
    outcome = np.zeros(32, np.int64)
    for final in (False, True):
        actual, ideal, actual_memory, _, bits = e.detector(
            initial, initial.copy(), labels, active, outcome, operators, final=final)
        assert np.array_equal(actual_memory, np.r_[labels, labels])
        assert np.array_equal(bits, np.r_[np.zeros(32), np.ones(32)])
        assert np.array_equal(actual, ideal)
        # The output operator sum is trace preserving on all32 input columns.
        assert la.norm(actual[:32].conj()@actual[:32].T
                       +actual[32:].conj()@actual[32:].T-np.eye(32)) < 1e-13


def test_complete_adaptive_execution_has_all_histories_and_nonvacuous_bound(receipt):
    rows = receipt['complete_record_execution']['executions']
    assert [row['retained_class_history_rows'] for row in rows] == [72, 648, 52488]
    assert F(rows[-1]['analytic_trace_bound_squared']) < 4
    intervals = [list(map(F, row['complete_record_trace_distance_observed_interval'])) for row in rows]
    assert intervals[0][0] > intervals[1][1] > intervals[2][1] > 0
    assert all(row['all_four_detector_outcomes_positive'] for row in rows)


def test_native_record_feedback_is_operational(monkeypatch):
    correct = e.execute(1)
    original = e.advance
    def erase_native_flag(*args, **kwargs):
        a, i, labels, active, outcomes = original(*args, **kwargs)
        return a, i, labels, np.zeros_like(active), outcomes
    monkeypatch.setattr(e, 'advance', erase_native_flag)
    wrong = e.execute(1)
    assert max(abs(a-b) for a, b in zip(correct['detector_probabilities'],
                                       wrong['detector_probabilities'], strict=True)) > 1e-5


def test_unread_rate_does_not_control_the_complete_archive():
    value = e.archive_countercontrol()
    unread = F(value['unread_trace_distance'])
    full = F(value['complete_record_trace_distance'])
    assert full > unread and full*full == 2*unread


@pytest.mark.parametrize('change', [
    lambda r: r['scope'].__setitem__('physical_species_or_masses_derived', True),
    lambda r: r['scope'].__setitem__('actual_memory_continues_across_interventions', 1),
    lambda r: r['scope'].__setitem__('native_compilation_of_added_detectors', True),
    lambda r: r['scope'].__setitem__('all_microscopic_adaptive_policies_converge', True),
    lambda r: r['exact_resource_certificate']['three_interval_one_percent_schedule'].__setitem__('total_native_steps', 1),
    lambda r: r['complete_record_execution']['law'].__setitem__('fine_labels_per_input', 96),
    lambda r: r['source_bindings'].pop('model/interacting_observation_process/PROOF.md'),
])
def test_receipt_rejects_changed_science_scope_and_custody(receipt, change):
    altered = deepcopy(receipt)
    change(altered)
    with pytest.raises(ValueError, match='reconstruction'):
        packet.validate(altered, receipt)


@pytest.mark.parametrize('durations,error', [([], F(1, 10)), ([F(1)], 0),
                                           ([F(-1)], F(1, 10)), ([0.5], F(1, 10)),
                                           ([F(1)], True)])
def test_invalid_exact_schedule_inputs_are_rejected(durations, error):
    with pytest.raises(ValueError):
        b.schedule(durations, error)


def test_retained_receipt_string_path_cold_replay():
    assert packet.verify(str(packet.OUTPUT)) is True
