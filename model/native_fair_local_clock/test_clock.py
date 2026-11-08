"""Clock, full raw-map, source admission and adversarial custody controls."""
from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as F
from math import comb
import json

import pytest

from . import clock as c, history as h, packet


EDGE = c.native.local.SEAMS[0]


def tiny_plan():
    return c.compile_plan(0, 1, (F(0), F(2, 630)),
                           (((EDGE, 'E'),),), ('E',))


def test_literal_preparation_prefix_and_one_counted_timeline():
    from model.native_fair_local_preparation import preparation
    plan = c.prepared_readout_plan(0, 80)
    assert plan['total_attempts'] == 7560
    assert plan['fixed_attempts'] == 1931
    assert plan['uniform_attempts'] == 5629
    assert plan['program']['commands'][:1931] == preparation.commands(32)
    assert plan['program']['declared_resources']['AST_nodes'] == 1933
    assert c.attempt_label(plan, 1931)['role'] == 'fixed'
    assert c.attempt_label(plan, 1932)['role'] == 'uniform'
    assert c.attempt_label(plan, 2520)['bin'] == 0
    assert c.attempt_label(plan, 2521)['bin'] == 1
    assert c.attempt_label(plan, 7560)['time'] == 12
    assert plan['epochs'][0]['background_exchange_time'] == F(589, 1260)
    assert plan['epochs'][1]['background_exchange_time'] == 4


@pytest.mark.parametrize('level', (0, 2, 9))
def test_growing_shared_duration_respects_original_admission(level):
    plan = c.prepared_readout_plan(level, 240*4**level//3)
    assert plan['total_attempts'] == 12*630*4**level
    assert c.native.admit(plan['program']) == plan['program']
    assert plan['program']['memory'] == {'capacity': 1, 'overflow': 'export_oldest'}


def test_capacity_does_not_override_literal_attempt_budget():
    assert 8*630*4**10 > c.native.MAX_REPEAT
    assert 12*630*4**10 > c.native.MAX_ATTEMPTS
    with pytest.raises(ValueError, match='resource/domain boundary'):
        c.prepared_readout_plan(10, 240*4**10//3)


def test_exact_clock_reservation_uses_uniform_count_not_total_count():
    plan = c.compile_plan(0, 1, (F(0), F(5, 630)),
                           (((EDGE, 'E'), (EDGE, 'C')),), ('mixture',))
    result = c.budget(plan, 33, 12, 6, 2)
    assert result['uniform_attempt_cone_upper'] == 12*comb(3, 3)*F(6, 630)**3
    assert result['complete_window_categorical_Poisson_upper'] == F(3*33**2, 630**2)
    assert result['fixed_attempt_clock_reservation_upper'] == F(2*33, 630)
    assert result['native_total_attempts'] == 5
    assert result['uniform_exchange_times'] == (F(3, 1260),)


def test_all_fixed_bin_has_a_nonzero_clock_reservation_cost():
    plan = c.compile_plan(0, 1, (F(0), F(1, 630)), (((EDGE, 'E'),),), ('E',))
    result = c.budget(plan, 33, 12, 6, 2)
    assert result['uniform_attempt_cone_upper'] == 0
    assert result['complete_window_categorical_Poisson_upper'] == 0
    assert result['fixed_attempt_clock_reservation_upper'] == F(33, 630)
    lo, hi = c.limit.exp_minus_interval(F(33, 630))
    assert 0 < 1-hi <= 1-lo <= result['fixed_attempt_clock_reservation_upper']


def test_coupon_clock_uses_only_actual_uniform_E_attempts():
    plan = c.prepared_readout_plan(0, 80)
    bound = c.readout_coverage_budget(plan, 1, 11, F(1, 10))
    assert bound['uniform_attempts'] == 5040
    assert bound['whole_seam_times'] == 8
    assert bound['coverage_failure_upper'] == F(11, 256)
    assert bound['coverage_or_erasure_upper'] == F(11, 256)+F(245, 2560)
    with pytest.raises(ValueError, match='uniform E'):
        c.readout_coverage_budget(plan, 0, 11, F(1, 10))


@pytest.mark.parametrize('partition,words', [
    ((F(0), F(1, 64)), ((),)),
    ((F(0), F(1, 630)), (((EDGE, 'E'), (EDGE, 'E')),)),
    ((F(0), F(1, 630)), ()),
])
def test_unaligned_or_overbooked_clock_rejected(partition, words):
    with pytest.raises(ValueError):
        c.compile_plan(0, 1, partition, words, ('E',))


@pytest.mark.parametrize('change', ['counter', 'law', 'bool', 'resource', 'command', 'FIFO', 'extra'])
def test_plan_mutations_rejected(change):
    value = tiny_plan()
    if change == 'counter': value['total_attempts'] += 1
    if change == 'law': value['program']['commands'][0]['law'] = 'C'
    if change == 'bool': value['level'] = False
    if change == 'resource': value['program']['declared_resources']['attempts_upper'] = 1
    if change == 'command': value['program']['commands'][0]['selector']['seam'][0] = False
    if change == 'FIFO': value['program']['memory']['overflow'] = 'halt'
    if change == 'extra': value['autonomous_unitarity'] = True
    with pytest.raises(ValueError): c.authenticate_plan(value)


def test_no_nested_loop_or_extra_command_budget():
    word = tuple((EDGE, 'E') for _ in range(4097))
    with pytest.raises(ValueError, match='bounded literal command'):
        c.compile_plan(0, 1, (F(0), F(4097, 630)), (word,), ('E',))


@pytest.fixture(scope='module')
def routed_pair():
    return tuple(c.routing.routed_window(level, 2) for level in (4, 6))


def attached(window, law='E', inner=0):
    plan = c.window_plan(window, 240*4**window['level']//3,
        (F(0), F(1, 630)), (((EDGE, law),),), ('E',))
    observed = tuple(e for e in window['internal_seams'] if set(e) <= set(window['root_ports']))
    row = c.attach_window(plan, window, inner, observed, window['root_ports'],
                         initial_error=c.canonical_error(plan, window, F(1, 3)),
                         initial_reference='Bernoulli(1/3) tensor the same independent R')
    return plan, window, row


def test_complete_rooted_transport_and_joint_error(routed_pair):
    left, right = map(attached, routed_pair)
    error = c.compare_attachments(*left, *right)
    assert error == min(F(1), sum(x[2]['budget']['native_to_common_finite_window_upper'] for x in (left, right)))
    assert right[2]['budget']['fixed_attempt_clock_reservation_upper'] < left[2]['budget']['fixed_attempt_clock_reservation_upper']


def test_mismatched_fixed_instrument_cannot_share_target(routed_pair):
    with pytest.raises(ValueError, match='fixed intervention'):
        c.compare_attachments(*attached(routed_pair[0]), *attached(routed_pair[1], law='C'))


@pytest.mark.parametrize('reference', [None, 'another joint reference'])
def test_initial_bounds_cannot_be_combined_without_one_joint_reference(routed_pair, reference):
    left = attached(routed_pair[0]); plan, window, declared = attached(routed_pair[1])
    changed = c.attach_window(plan, window, declared['root_cell_radius'],
        declared['observed_uniform_seams'], declared['final_ports'],
        initial_error=declared['budget']['initial_joint_half_trace_upper'], initial_reference=reference)
    with pytest.raises(ValueError, match='joint reference'):
        c.compare_attachments(*left, plan, window, changed)


def test_incident_record_requires_its_outside_endpoint_in_root(routed_pair):
    window = routed_pair[0]
    plan, _, _ = attached(window)
    root = set(window['root_ports'])
    crossing = next(e for e in window['internal_seams'] if bool(e[0] in root) != bool(e[1] in root))
    with pytest.raises(ValueError, match='both endpoints'):
        c.attach_window(plan, window, 0, (crossing,), window['root_ports'])
    allowed = c.attach_window(plan, window, 1, (crossing,), window['root_ports'])
    assert allowed['graph_radius_margin'] == 1
    assert len(allowed['root_ports']) == 48


def test_typed_window_and_observation_mutation_rejected(routed_pair):
    window = deepcopy(routed_pair[0]); window['orientation_preserved'] = 1
    with pytest.raises(ValueError): attached(window)
    left = attached(routed_pair[0]); right = attached(routed_pair[1])
    right[2]['root_cell_radius'] = False
    with pytest.raises(ValueError): c.compare_attachments(*left, *right)


def test_full_625_word_reference_control():
    control = h.complete_fixed_control()
    assert control['complete_raw_words'] == 625
    assert control['zero_words'] == 617
    assert control['reference_marginal_preserved_exactly'] is True


@pytest.fixture(scope='module')
def raw_query():
    plan = tiny_plan(); x, y = (1 << EDGE[0], 1 << EDGE[1])
    rho = {(x, 0, x, 0): F(1, 2), (y, 1, y, 1): F(1, 2),
           (x, 0, y, 1): F(1, 2), (y, 1, x, 0): F(1, 2)}
    first = c.native._outcomes(plan['program'], plan['program']['commands'][0], rho)
    s = next(i for i, out in enumerate(first) if out.density)
    out = first[s]
    second = c.native._outcomes(plan['program'], plan['program']['commands'][1], out.density,
                                out.memory, out.archive)
    t = next(i for i, out in enumerate(second) if out.density and out.record.seam == EDGE)
    return plan, rho, h.query_word(plan, (s, t), rho)


def test_actual_fixed_then_uniform_raw_query_retains_fifo(raw_query):
    plan, rho, query = raw_query
    assert h.authenticate_query(plan, query, rho)
    assert query['records'] == query['archive']+query['memory']
    assert len(query['archive']) == len(query['memory']) == 1
    assert query['attempts'] == 2
    assert query['records'][0].selector == 'fixed'
    assert query['records'][1].selector == 'uniform'


@pytest.mark.parametrize('change', ['density', 'coin', 'raw', 'FIFO', 'program', 'attemptbool', 'slotbool'])
def test_forged_query_rejected(raw_query, change):
    plan, rho, original = raw_query; query = deepcopy(original)
    if change == 'density': query['density'] = {}
    if change == 'coin': query['records'] = (replace(query['records'][0], coin=True), *query['records'][1:])
    if change == 'raw': query['records'] = (replace(query['records'][0], raw=100), *query['records'][1:])
    if change == 'FIFO': query['archive'], query['memory'] = query['memory'], query['archive']
    if change == 'program': query['program_digest'] = '0'*64
    if change == 'attemptbool': query['attempts'] = True
    if change == 'slotbool': query['slots'] = (True, query['slots'][1])
    with pytest.raises(ValueError): h.authenticate_query(plan, query, rho)


def test_zero_raw_word_is_preserved_without_success_normalization():
    plan = tiny_plan(); rho = {(1, 0, 1, 0): F(1)}
    query = h.query_word(plan, (0, 0), rho)
    assert query['density'] == {}
    assert tuple(record.kind for record in query['records']) == ('hold', 'hold')
    assert h.authenticate_query(plan, query, rho)


def test_actual_local_observation_transports_records_and_reference():
    window = c.routing.routed_window(2, 1)
    plan = c.window_plan(window, 1, (F(0), F(2, 10080)),
                          (((EDGE, 'E'), (EDGE, 'C')),), ('E',))
    attached = c.attach_window(plan, window, 0, (EDGE,), EDGE)
    x, y = (1 << EDGE[0], 1 << EDGE[1])
    rho = {(x, 0, x, 0): F(1, 2), (y, 1, y, 1): F(1, 2)}
    first = c.native._outcomes(plan['program'], plan['program']['commands'][0], rho)
    s = next(i for i, out in enumerate(first) if out.density)
    second = c.native._outcomes(plan['program'], plan['program']['commands'][1], first[s].density,
                                first[s].memory, first[s].archive)
    t = next(i for i, out in enumerate(second) if out.density)
    query = h.query_word(plan, (s, t), rho)
    observed = h.observe_query(plan, window, attached, query, rho)
    assert observed['global_tick_addresses_observed'] is False
    assert tuple(event[:3] for event in observed['events']) == ((0, 'fixed', 0), (0, 'fixed', 1))
    assert all(event[3].seam == EDGE for event in observed['events'])
    assert observed['density'] == query['density']
    zero = h.query_word(plan, (0, 0), rho)
    output = h.observe_query(plan, window, attached, zero, rho)
    assert output['density'] == {}
    assert len(output['events']) == 2


def test_literal_source_and_service_domain_controls():
    result = h.source_domain_controls()
    assert result['attempts_before_host_stop'] == 0
    assert result['unevaluated_fixed_intercell_command_rejected'] is True
    assert result['source_word_operation_rejected_by_native_grammar'] is True
    assert result['growing_federation_coherent_load_map_produced'] is False


def test_high_level_query_cannot_bypass_the_owning_source_constructor():
    plan = c.compile_plan(3, 1, (F(0), F(1, 630*4**3)), ((),), ('E',))
    with pytest.raises(ValueError, match='cannot construct'):
        h.query_word(plan, (0,), {(1, 0, 1, 0): F(1)})


def test_receipt_type_mutation_rejected(tmp_path, monkeypatch):
    original = {'schema': 'test', 'scope': {'source_preserved': True}, 'bound': '1/7'}
    monkeypatch.setattr(packet, 'produce', lambda: original)
    altered = deepcopy(original); altered['scope']['source_preserved'] = 1
    target = tmp_path/'bad.json'; target.write_text(json.dumps(altered))
    with pytest.raises(ValueError): packet.verify(target)


def test_current_deterministic_receipt():
    assert packet.verify()
