"""Complete actual-source -> native-record counter -> native-pointer readout.

The source is scheduled first and idles during terminal pointer service.
The detector reports the retained source prefix's relative charge. The
controller never receives the terminal loads or a precomputed degree.
"""
from fractions import Fraction as F
import hashlib
from itertools import cycle
import json
from pathlib import Path

import native_cycle_meter as meter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]


def source_prefix(initial_G, schedule, coins):
    """Advance the actual source and its finite record-only observer."""
    initial_G = tuple(initial_G)
    if len(initial_G) != 12 or any(type(x) is not int or x not in (0, 1) for x in initial_G) or sum(initial_G) != 1:
        return {'admitted': False, 'reason': 'source_outside_one_load_zero_tape_domain',
                'source_initial': list(initial_G), 'records': [], 'counter': None}
    schedule = tuple(tuple(edge) for edge in schedule)
    coins = tuple(coins)
    if len(coins) != len(schedule) or any(type(c) is not int or c not in (0, 1) for c in coins):
        raise ValueError('one literal coin for every scheduled source attempt')
    if any(edge not in meter.EDGES for edge in schedule):
        raise ValueError('actual oriented source seams required')
    loads = initial_G
    counter = [0]*12
    records = []
    for event, (edge, coin) in enumerate(zip(schedule, coins)):
        loads, (_, current) = meter.native(loads, edge, coin)
        # This is the entire observer update. It reads no occupation value.
        counter[edge[0]] = (counter[edge[0]]+current) % 6
        counter[edge[1]] = (counter[edge[1]]-current) % 6
        records.append({'stage': 'source', 'event': event, 'seam': list(edge),
                        'coin': coin, 'current': current, 'counter_after': list(counter)})
    return {'admitted': True, 'initial_G': list(initial_G), 'counter': counter,
            'records': records, 'source_final_loads_for_audit_only': list(loads),
            'source_history_probability': str(F(1, 2**len(records)))}


def execute_joined(initial_G, schedule, source_coins, region_mask,
                   pointer_initial_port, pointer_coins, retries=16):
    """Single joined policy, with source and pointer tapes in one output."""
    if type(region_mask) is not int or not 0 <= region_mask < 1 << 12:
        raise ValueError('declared twelve-port region mask required')
    source = source_prefix(initial_G, schedule, source_coins)
    if not source['admitted']:
        return {'success': False, 'reason': source['reason'], 'source': source,
                'pointer_executed': False, 'relative_charge': None, 'spatial_degree': None,
                'combined_native_tape': [], 'combined_history_probability': '1'}
    # The terminal source loads and stored expected degrees are unavailable
    # to this policy. Only its accumulated current-derived counter is read.
    regional_class = sum(source['counter'][v] for v in range(12) if region_mask >> v & 1) % 6
    pointer = meter.execute_class(pointer_initial_port, regional_class, pointer_coins, retries)
    combined = list(source['records'])
    for offset, row in enumerate(pointer['records'], start=len(combined)):
        combined.append({**row, 'stage': 'pointer_'+row['stage'], 'event': offset})
    probability = F(source['source_history_probability'])*F(pointer.get('fine_history_probability', '1'))
    return {'success': pointer['success'], 'reason': pointer['reason'],
            'region_mask': region_mask, 'source': source,
            'observer_regional_class': regional_class,
            'pointer_executed': True, 'pointer': pointer,
            'spatial_degree': pointer['winding'],
            'pointer_winding': pointer['winding'],
            'relative_charge': -pointer['winding'] if pointer['success'] else None,
            'combined_native_tape': combined,
            'combined_history_probability': str(probability),
            'clock': {'source_prefix_attempts': len(source['records']),
                      'pointer_attempts': len(pointer['records']),
                      'total_native_attempts': len(combined),
                      'source_during_pointer': 'idle under the declared addressed scheduler',
                      'reported_charge_time': 'end of the retained source prefix; unchanged through terminal meter'}}


def discovery_hold_coins(initial_port, retries):
    for a, b in meter.EDGES:
        for _ in range(retries):
            yield 1 if initial_port == a else 0


def run():
    source_path = HERE.parent/'spatial/NATIVE_DUAL_WINDING.json'
    paths = [Path(__file__), HERE/'native_cycle_meter.py', source_path,
             ROOT/'reverse-engineering-reality/model/observer/native.py',
             ROOT/'openwave-tuning/candidate/oph_openwave/geometry.py']
    source_pins_before = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    source_bytes = source_path.read_bytes()
    source = json.loads(source_bytes)
    source_hash = hashlib.sha256(source_bytes).hexdigest()
    retries = 16
    hold = F(1, 2**retries)
    reset = [meter.reset_success(v, retries) for v in range(12)]
    executed_source_histories = []
    budgets = []
    source_hold_events = 0
    source_move_events = 0
    counter_checks = 0
    complete_budget_checks = 0
    for case_index, seed in enumerate(source['seed_records']):
        initial_G = seed['initial_G']
        schedule = [tuple(row['seam']) for row in seed['actual_history']]
        source_law = {v: {q: F() for q in (-1, 0, 1)} for v in range(12)}
        for history_index, history in enumerate(seed['complete_prescribed_schedule_histories']):
            coins = [row['coin'] for row in history['events']]
            actual = source_prefix(initial_G, schedule, coins)
            assert actual['admitted']
            assert actual['source_final_loads_for_audit_only'] == history['final_loads']
            p = F(actual['source_history_probability'])
            assert p == F(history['probability'])
            decoded = {0: 0, 1: 1, 5: -1}
            for v in range(12):
                degree = decoded[actual['counter'][v]]
                assert degree == initial_G[v]-history['final_loads'][v]
                assert abs(degree-history['measured_phase_degrees'][v]) < 1e-12
                source_law[v][degree] += p
                counter_checks += 1
            source_hold_events += sum(row['current'] == 0 for row in actual['records'])
            source_move_events += sum(row['current'] != 0 for row in actual['records'])
            executed_source_histories.append({'source_case': case_index,
                                             'source_history': history_index, **actual})
        for vertex in range(12):
            assert sum(source_law[vertex].values(), F()) == 1
            for pointer_initial in range(12):
                conclusive = {q: source_law[vertex][q]*reset[pointer_initial]*(1-hold)**(6*abs(q))
                              for q in (-1, 0, 1)}
                timeout = sum(source_law[vertex][q]*(1-reset[pointer_initial]*(1-hold)**(6*abs(q)))
                              for q in (-1, 0, 1))
                assert sum(conclusive.values(), F())+timeout == 1
                complete_budget_checks += 1
            budgets.append({'source_case': case_index, 'region_mask': 1 << vertex,
                            'complete_source_degree_distribution': {str(q): str(p) for q, p in source_law[vertex].items()},
                            'joint_conclusive_degree_probability': 'source_degree_probability(q)*p_reset(pointer_initial)*(1-2^-R)^(6|q|)',
                            'joint_timeout_probability': 'sum_q source_degree_probability(q)*(1-p_reset(pointer_initial)*(1-2^-R)^(6|q|))'})
    assert len(executed_source_histories) == 456
    assert counter_checks == 456*12 and complete_budget_checks == 132*12*12
    initial_G = [1]+[0]*11
    edge = (0, 1)
    examples = []
    for name, coins, region, pointer_initial, expected in (
        ('positive_degree_actual_move', (0,), 1, 5, 1),
        ('negative_degree_actual_move', (0,), 2, 11, -1),
        ('zero_degree_actual_hold', (1,), 1, 5, 0),
        ('zero_degree_empty_prefix', (), 1, 11, 0),
        ('zero_degree_empty_region_after_move', (0,), 0, 5, 0),
        ('zero_degree_full_region_after_move', (0,), (1 << 12)-1, 11, 0),
    ):
        schedule = (edge,) if coins else ()
        result = execute_joined(initial_G, schedule, coins, region, pointer_initial,
                                meter.favorable_coins(pointer_initial, (expected,), retries), retries)
        assert result['success'] and result['spatial_degree'] == expected
        assert result['combined_history_probability'] == str(F(1, 2**len(result['combined_native_tape'])))
        examples.append({'case': name, **result})
    for name, pointer_coins, expected_reason in (
        ('discovery_failure', discovery_hold_coins(0, retries), 'discovery_timeout'),
        ('routing_failure', cycle((0,)), 'routing_timeout'),
    ):
        result = execute_joined(initial_G, (edge,), (0,), 1, 0, pointer_coins, retries)
        assert not result['success'] and result['reason'] == expected_reason and result['relative_charge'] is None
        examples.append({'case': name, **result})
    result = execute_joined(initial_G, (edge,), (0,), 1, 1, cycle((1,)), retries)
    assert not result['success'] and result['reason'] == 'measurement_timeout'
    examples.append({'case': 'measurement_failure', **result})
    result = execute_joined(initial_G, (edge,), (0,), (1 << 12)-1, 0,
                            discovery_hold_coins(0, retries), retries)
    assert not result['success'] and result['reason'] == 'discovery_timeout'
    assert result['observer_regional_class'] == 0 and result['relative_charge'] is None
    examples.append({'case': 'full_region_zero_class_retains_reset_failure', **result})
    enlarged = [7]+[0]*11
    result = execute_joined(enlarged, (), (), 1, 0, cycle((0, 1)), retries)
    assert not result['success'] and not result['pointer_executed']
    examples.append({'case': 'charge_seven_alias_admission_refused', **result})
    # A deliberately omitted current would leave class zero after this move.
    mutation = {'actual_source_prefix': source_prefix(initial_G, (edge,), (0,)),
                'mutated_counter_if_current_omitted': [0]*12,
                'rejected_by_actual_relative_degree': True}
    assert mutation['actual_source_prefix']['counter'][:2] == [1, 5]
    assert source_hash == hashlib.sha256(source_path.read_bytes()).hexdigest()
    source_pins_after = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    assert source_pins_before == source_pins_after
    return {'schema': 'oph.row1.joined_literal_native_readout.v1',
            'checks': {'one_joined_actual_source_and_pointer_executor': True,
                       'observer_counter_reads_only_actual_native_currents': True,
                       'all456_frozen_source_histories_reexecuted': True,
                       'all_19008_source_pointer_input_probability_budgets_exact': True,
                       'zero_move_hold_and_all_timeout_stages_retained': True,
                       'source_idle_terminal_clock_explicit': True,
                       'no_terminal_charge_or_degree_supplied_to_policy': True,
                       'larger_charge_alias_admission_refused': True,
                       'source_pins_unchanged': True},
            'source_history_count': len(executed_source_histories),
            'source_hold_events_retained': source_hold_events,
            'source_moving_events_retained': source_move_events,
            'counter_phase_degree_comparisons': counter_checks,
            'complete_probability_budgets_checked': complete_budget_checks,
            'pointer_retries': retries,
            'pointer_reset_success_by_initial_port': list(map(str, reset)),
            'maximum_source_prefix_attempts': max(len(x['records']) for x in executed_source_histories),
            'maximum_joined_native_attempts': 39*retries+max(len(x['records']) for x in executed_source_histories),
            'complete_source_histories': executed_source_histories,
            'complete_joint_probability_laws': budgets,
            'joined_execution_examples': examples,
            'record_omission_mutation': mutation,
            'protocol': 'known actual one-load initialG, zero observer counters; prescribed source seam schedule with all faircoin histories; primary source then idles; unknown one-load pointer uses native current-only reset and bounded cycle meter; source and apparatus tapes remain joined',
            'unconditional_outcomes': 'relative_charge=-spatial_degree on conclusive native readout, or explicit admission/discovery/routing/measurement/capacity flag with no charge assigned',
            'setting_probability_scope': 'Each source case, region and unknown pointer initial state specifies a separate protocol input. No random prior over those inputs is claimed. The fair-coin histories, reset and meter outcomes are unconditional and normalized within each configured protocol.',
            'full_tree_semantics': 'Every executed attempt contributes one literal coin/current letter and factor1/2. execute_joined and native_cycle_meter.execute enumerate the entire finite stopping tree; the exact probability law retains its grouped success and all timeout complements.',
            'scope': 'Complete operational readout of the primary source spatial class, under its declared native observer and scheduler. The spatial cohomology interpretation comes from the independently checked dual-cell adapter; no physical electromagnetic unit or automatic platform award is inferred.',
            'source_sha256': source_pins_after,
            'official_rows_changed': []}


if __name__ == '__main__':
    result = run()
    (HERE/'JOINED_READOUT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['checks', 'source_history_count',
                                           'source_hold_events_retained', 'source_moving_events_retained',
                                           'complete_probability_budgets_checked', 'maximum_joined_native_attempts']}))
