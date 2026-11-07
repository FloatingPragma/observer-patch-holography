"""Actual native-program histories and adversarial provenance controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import pytest

from ..native_program import interpreter as program
from ..observer import native
from . import dataflow as d


def loads(*occupied):
    return tuple(int(i in occupied) for i in range(12))


@pytest.fixture(scope='module')
def pair():
    value = d.feedback_program()
    return [d.execute(value, loads(*occupied)) for occupied in ((0, 2), (1, 2))]


def test_exact_complete_interpreter_distributions_show_controller_influence(pair):
    assert [d.receiver_probability(result, 3) for result in pair] == [F(1, 4), F()]
    for result in pair:
        assert result['complete_raw_branch_count'] == 625
        assert result['positive_branch_count'] == 4
        assert sum(row['probability'] for row in result['branches']) == 1
        assert sum(row['events'] is None for row in result['branches']) == 621
        for row in result['branches']:
            assert len(row['records']) == 2
            if row['events']:
                first, second = row['records']
                assert first.seam == (0, 1)
                assert second.seam == ((2, 3) if first.current > 0 else (2, 10))


def test_full_footprints_authenticate_controller_and_storage_beyond_endpoints(pair):
    for result in pair:
        for row in result['branches']:
            if row['events'] is None:
                continue
            first, second = row['events']
            assert set(first['record']['seam']).isdisjoint(second['record']['seam'])
            assert second['endpoint_parents'] == []
            assert second['controller_parents'] == second['storage_parents'] == second['parents'] == [0]
            assert second['ancestors'] == [0] and second['height'] == 1
            assert second['before_writers']['fifo'] == 0
            assert second['before_values']['fifo'] == first['after_values']['fifo']
            assert not d.independent(first, second)
            endpoint_only = [{key: [x for x in event[key] if x.startswith('port:')]
                              for key in ('write_set', 'causal_support')} for event in (first, second)]
            assert d.independent(*endpoint_only)
            assert d.verify_trace(result['program'], result['loads'], row['records'], row['events'])


@pytest.mark.parametrize('field', ['before_writers', 'before_values', 'read_set', 'controller_parents'])
def test_forged_or_omitted_control_provenance_is_rejected(pair, field):
    result = pair[0]
    row = next(row for row in result['branches'] if row['events'])
    altered = deepcopy(row['events'])
    if field == 'before_writers':
        altered[1][field]['fifo'] = None
    elif field == 'before_values':
        altered[1][field]['fifo'] = ()
    elif field == 'read_set':
        altered[1][field].remove('fifo')
    else:
        altered[1][field] = []
    with pytest.raises(ValueError):
        d.verify_trace(result['program'], result['loads'], row['records'], altered)


def test_zero_probability_raw_histories_are_retained_but_not_authenticated_as_executions(pair):
    result = pair[0]
    row = next(row for row in result['branches'] if not row['probability'])
    assert row['events'] is None
    with pytest.raises(ValueError):
        d.replay_trace(result['program'], result['loads'], row['records'])


def test_same_fixed_remote_seam_removes_the_controller_signal():
    value = d.feedback_program(adaptive=False)
    results = [d.execute(value, loads(*occupied)) for occupied in ((0, 2), (1, 2))]
    assert [d.receiver_probability(result, 3) for result in results] == [F(1, 2), F(1, 2)]
    for result in results:
        for row in result['branches']:
            if row['events']:
                assert row['events'][1]['controller_parents'] == []
                assert row['events'][1]['storage_parents'] == [0]


def test_disjoint_endpoint_updates_commute_but_shared_record_storage_does_not():
    # Independent check of the native update component, over its full carrier.
    for state in native.states(binary=False):
        for first_coin, second_coin in product((0, 1), repeat=2):
            left = native.native(native.native(state, (0, 1), first_coin)[0], (2, 3), second_coin)[0]
            right = native.native(native.native(state, (2, 3), second_coin)[0], (0, 1), first_coin)[0]
            assert left == right
    result = d.execute(d.feedback_program(adaptive=False), loads(0, 2))
    for row in result['branches']:
        if row['events']:
            assert not d.independent(*row['events'])


def test_zero_current_commits_and_fifo_capacity_keep_the_actual_storage_edges(pair):
    result = d.execute(d.feedback_program(), loads(4, 5))
    for row in result['branches']:
        if row['events']:
            first, second = row['events']
            assert first['record']['current'] == second['record']['current'] == 0
            assert first['before_values']['port:0'] == first['after_values']['port:0'] == 0
            assert first['after_writers']['port:0'] == 0
            assert second['controller_parents'] == second['storage_parents'] == [0]
    larger = d.execute(d.feedback_program(capacity=2), loads(0, 2))
    assert [r['probability'] for r in larger['branches']] == [r['probability'] for r in pair[0]['branches']]
    for row in larger['branches']:
        if row['events']:
            assert 'archive' not in row['events'][1]['write_set']
            assert row['events'][1]['storage_parents'] == [0]


def test_holds_write_storage_and_complete_execution_boundary_is_explicit():
    value = program.program(biased=False, commands=[{'op': 'hold', 'repeat': 1},
        {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]}, 'law': 'C', 'repeat': 1}])
    result = d.execute(value, loads(0, 2))
    for row in result['branches']:
        if row['events']:
            first, second = row['events']
            assert first['read_set'] == ['archive', 'fifo']
            assert first['write_set'] == ['fifo']
            assert second['endpoint_parents'] == [] and second['storage_parents'] == [0]
    halted = deepcopy(value)
    halted['memory']['overflow'] = 'halt'
    halted['declared_resources'] = program._native_cost(halted)
    with pytest.raises(ValueError, match='complete bounded'):
        d.execute(halted, loads(0, 2))
    too_long = program.program(commands=[{'op': 'hold', 'repeat': 65}])
    with pytest.raises(ValueError, match='bounded fully'):
        d.replay_trace(too_long, loads(0, 2), ())


def test_longer_actual_history_has_authenticated_skipping_edges_and_canonical_height():
    fixed = {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]}, 'law': 'C', 'repeat': 1}
    value = program.program(biased=False, commands=[fixed, {'op': 'hold', 'repeat': 1},
                                                    fixed, {'op': 'hold', 'repeat': 1}])
    result = d.execute(value, loads(0, 2))
    for row in result['branches']:
        if row['events']:
            assert [event['height'] for event in row['events']] == [0, 1, 2, 3]
            assert [event['ancestors'] for event in row['events']] == [[], [0], [0, 1], [0, 1, 2]]
            assert row['events'][2]['endpoint_parents'] == [0]
            assert row['events'][2]['parents'] == [0, 1]
            assert all(parent < event['event_id'] for event in row['events'] for parent in event['parents'])


def test_collision_containing_classical_carrier_uses_its_actual_native_payload():
    value = program.program(shell='H78', biased=False, commands=[
        {'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]}, 'law': 'E', 'repeat': 1}])
    result = d.execute(value, (2,)+(0,)*11)
    assert result['positive_branch_count'] == 2
    for row in result['branches']:
        if row['events']:
            event, = row['events']
            assert event['record']['label'] == 2 and event['record']['current'] == 1
            assert event['after_values']['port:0'] == event['after_values']['port:1'] == 1


def test_certificate_and_domain_controls():
    result = d.certificate()
    assert json.loads(json.dumps(result)) == result
    assert result['receiver_probabilities'] == ['1/4', '0']
    assert result['complete_raw_branch_counts'] == [625, 625]
    assert result['full_footprint_is_operational_dependency_overapproximation']
    assert result['adaptive_receiver_probability_witnesses_influence']
    assert not result['storage_chain_alone_excludes_physical_Lorentzian_limit']
    assert not result['physical_spacetime_order_or_Lorentz_cone_derived']
    with pytest.raises(ValueError):
        d.execute(d.feedback_program(), (0,)*12)
    with pytest.raises(ValueError):
        d.execute(program.program(cells=2), loads(0, 2))
