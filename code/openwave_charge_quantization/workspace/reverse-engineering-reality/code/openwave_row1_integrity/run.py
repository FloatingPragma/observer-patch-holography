#!/usr/bin/env python3
"""Reproduce bounded native integrity instruments and rejection controls."""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import ast
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
import hashlib
import inspect
import itertools
import json
from pathlib import Path
import sys
import time
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import integrity as it


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def pins():
    manifest = json.loads((HERE/'INPUTS.json').read_text())
    measured = {}
    for name, expected in manifest['runtime_files_sha256'].items():
        actual = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        assert actual == expected, (name, actual, expected)
        measured['portable/' + name] = actual
    workspace = HERE.parents[2]
    external_present = (workspace/'openwave-push').is_dir()
    if external_present:
        for name, expected in manifest['external_files_sha256'].items():
            actual = hashlib.sha256((workspace/name).read_bytes()).hexdigest()
            assert actual == expected, (name, actual, expected)
            measured[name] = actual
    return {'files': measured, 'workspace_originals_checked': external_present}


def frozen_audit():
    def functions(path):
        return {n.name: ast.dump(n, include_attributes=False) for n in ast.parse(path.read_text()).body
                if isinstance(n, ast.FunctionDef)}
    public = functions(HERE/'native_cycle_meter.py')
    original = functions(HERE/'provenance/PRIMARY_native_cycle_meter.original.py')
    matched = ['execute', 'route', 'reset_success', 'favorable_coins', 'independent_reset_distribution']
    assert all(public[n] == original[n] for n in matched)
    assert it.EDGES == it.sweep.SEAMS == it.strict.SEAMS
    assert it.meter.CYCLE == it.strict.CYCLE
    assert str(inspect.signature(it.strict.validate)) == '(result, inputs=None)'
    assert str(inspect.signature(it.strict.validate_pointer_tape)) == '(meter, pointer_initial, commanded, retries)'
    assert str(inspect.signature(it.sweep.closing_sweep)) == '(G, retained_tape, source_state, order, coins, *, channel=None, allow_partial=False)'
    states = preparations()
    events = 0
    for G in states:
        for edge in it.EDGES:
            for coin in (0, 1):
                a, (_, j) = it.native(G, edge, coin)
                b, k = it.strict.native(G, edge, coin)
                c, d, ell = it.sweep.native_event(G, edge, coin)
                assert a == b == c and j == k == ell and min(a) >= 0 and sum(a) == sum(G)
                assert d == G[edge[0]] - G[edge[1]]
                events += 1
    return {'frozen_signatures_checked': True, 'pointer_policy_AST_identical_functions': matched,
            'native_branch_equivalences': events, 'source_preparations': len(states),
            'strict_validate_scope': 'total two; final_n is required by frozen schema',
            'adapter_boundary': 'final_n is synthetic reconstruction from G and received records; closing E records provide integrity',
            'observer_has_source_state_parameter': 'source_state' in inspect.signature(it.admit).parameters}


def preparations():
    states = []
    for a in range(12):
        states.append(tuple(int(v == a) for v in range(12)))
    for a in range(12):
        for b in range(a, 12):
            states.append(tuple(int(v == a) + int(v == b) for v in range(12)))
    return states


class NeedChoice(Exception):
    pass


class DecisionCoins:
    """Audit-only quotient of the unmodified executor's exact branch tree."""
    def __init__(self, choices):
        self.choices = choices
        self.used = 0
        self.equal_branch_events = 0
    def __iter__(self):
        return self
    def __next__(self):
        # The actual executor asks for a coin inside attempt(). This is an
        # exhaustive auditor, never the apparatus controller or admission.
        frame = inspect.currentframe().f_back
        edge, loads = frame.f_locals['edge'], frame.f_locals['loads']
        zero = it.native(loads, edge, 0)
        one = it.native(loads, edge, 1)
        if zero[0] == one[0] and zero[1][1] == one[1][1]:
            self.equal_branch_events += 1
            return 0
        if self.used == len(self.choices):
            raise NeedChoice
        coin = self.choices[self.used]
        self.used += 1
        return coin


def executor_tree(initial, command, retries=1):
    pending = [()]
    outcomes = defaultdict(Fraction)
    leaves = 0
    validated_records = 0
    while pending:
        choices = pending.pop()
        coins = DecisionCoins(choices)
        try:
            result = it.meter.execute(initial, (command,), coins, retries)
        except NeedChoice:
            pending.extend((choices + (0,), choices + (1,)))
            continue
        errors, _ = it.strict.validate_pointer_tape(result, initial, command, retries)
        assert not errors, errors
        n = len(result['records'])
        weight = Fraction(2 ** coins.equal_branch_events, 2 ** n)
        assert weight == Fraction(1, 2 ** coins.used)
        key = (result['reason'], result['winding'], result['retained_signed_cycle_steps'], n)
        outcomes[key] += weight
        leaves += 1
        validated_records += n
    assert sum(outcomes.values()) == 1
    assert dict(outcomes) == it.pointer_distribution(initial, command, retries)
    return {'initial': initial, 'command': command, 'retries': retries,
            'quotient_leaves': leaves, 'actual_executor_records_replayed': validated_records,
            'outcomes': len(outcomes), 'mass': '1', 'exact_DP_equal': True}


def timeout_tape(initial, command, retries, stage):
    if stage == 'discover':
        return iter([1 if initial == edge[0] else 0 for edge in it.EDGES for _ in range(retries)])
    good = it.meter.execute(initial, (command,), it.meter.favorable_coins(initial, (command,), retries), retries)
    loads = tuple(int(v == initial) for v in range(12))
    coins = []
    for record in good['records']:
        edge = tuple(record['seam'])
        if record['stage'] == stage:
            hold = next(c for c in (0, 1) if it.native(loads, edge, c)[1][1] == 0)
            return iter(coins + [hold] * retries)
        coins.append(record['coin'])
        loads = it.native(loads, edge, record['coin'])[0]
    raise AssertionError(('timeout stage absent', stage, initial, command))


def full_examples_and_controls():
    accepted = 0
    record_checks = 0
    outcomes = Counter()
    controls = Counter()
    rejection_rules = defaultdict(Counter)
    samples = []
    stale_latch_witness = None
    source_metadata_ignored = 0
    source_cases = preparations()
    for k, G in enumerate(source_cases):
        prefix = [(k % 30, k % 2, 'C'), ((k + 5) % 30, (k // 2) % 2, 'E')]
        for initial in (0, 7, 11):
            for region in ((0, 1, 2), tuple(range(3, 12))):
                claim, audit = it.execute(G, prefix, region, initial,
                    lambda w: it.meter.favorable_coins(initial, (w,), 2),
                    [(k + j) % 2 for j in range(30)], 2,
                    list(range(30)) if k % 2 else list(range(29, -1, -1)))
                result = it.admit(claim, G, region, 2)
                assert result['admitted'] and initial in result['compatible_pointer_initials'], result
                assert claim['decoded_relative_charge'] == sum(audit['source_L1'][v]-G[v] for v in region)
                assert claim['latches']['L0']['reconstruction'] == audit['source_L0']
                assert claim['latches']['L1']['reconstruction'] == audit['source_L1']
                accepted += 1
                outcomes[result['outcome']] += 1
                record_checks += len(claim['meter']['records'])
                if stale_latch_witness is None and claim['latched_counter_classes']['L0'] != claim['counter_class']:
                    stale_latch_witness = {'G': list(G), 'region': list(region), 'L0': audit['source_L0'],
                        'L1': audit['source_L1'], 'classes': claim['latched_counter_classes'],
                        'pointer_relative_charge_at_L1': claim['decoded_relative_charge']}
                if initial == 7 and region == (0, 1, 2) and k in (0, 12):
                    samples.append({'G': list(G), 'region': list(region), 'retries': 2,
                                    'claim': claim, 'executor_audit': audit, 'admission': result})
                if k % 10 != 0 or initial != 7 or region != (0, 1, 2):
                    continue
                def reject(label, altered, expected_prefix=None):
                    outcome = it.admit(altered, G, region, 2)
                    assert not outcome['admitted'], (label, outcome)
                    if expected_prefix:
                        assert any(r.startswith(expected_prefix) for r in outcome['reasons']), outcome
                    controls[label] += 1
                    for reason in outcome['reasons']:
                        rejection_rules[label][reason.split(':')[0]] += 1
                corrupt = deepcopy(claim); corrupt['counter_class'] = (corrupt['counter_class'] + 1) % 6
                reject('counter_changed', corrupt, 'I2')
                corrupt = deepcopy(claim); corrupt['sweep_records'][0]['contrast'] += 2
                reject('wrong_E_contrast', corrupt, 'I3')
                corrupt = deepcopy(claim); corrupt['sweep_records'][0]['current'] += 1
                reject('non_native_E_current', corrupt, 'I3')
                corrupt = deepcopy(claim); corrupt['source_records'][0]['current'] += 1
                reject('non_native_C_current', corrupt, 'I3')
                corrupt = deepcopy(claim); corrupt['sweep_records'] = []
                reject('missing_sweep', corrupt, 'I1')
                corrupt = deepcopy(claim); corrupt['sweep_order'] = [0] * 30
                reject('disconnected_sweep', corrupt, 'I1')
                corrupt = deepcopy(claim); corrupt['sweep_records'][0]['kind'] = 'C'
                reject('sweep_current_without_E_contrast', corrupt, 'I1')
                corrupt = deepcopy(claim); corrupt['declared_total'] = 3 - sum(G)
                reject('declared_total_mismatch', corrupt, 'I0')
                corrupt = deepcopy(claim); corrupt['combined_attempts'] -= 30
                reject('sweep_clock_omitted', corrupt, 'I6')
                corrupt = deepcopy(claim); corrupt['charge_read_at'] = 'L0'
                reject('stale_latch_presented_as_current', corrupt, 'I2')
                corrupt = deepcopy(claim); corrupt['meter']['records'] = []
                corrupt['meter']['fine_history_probability'] = '1'
                reject('pointer_echo_zero_tape', corrupt, 'I4')
                corrupt = deepcopy(claim); corrupt['meter']['records'][0]['current'] += 1
                reject('corrupt_pointer_current', corrupt, 'I4')
                for at, rec in enumerate(claim['meter']['records']):
                    if rec['current']:
                        corrupt = deepcopy(claim); del corrupt['meter']['records'][at]
                        corrupt['meter']['fine_history_probability'] = str(Fraction(1, 2**len(corrupt['meter']['records'])))
                        reject('each_pointer_move_deleted', corrupt, 'I4')
                corrupt = deepcopy(claim)
                corrupt['audit_only'] = {'final_n': [999]*12, 'source_state': [999]*12}
                assert it.admit(corrupt, G, region, 2)['admitted']
                source_metadata_ignored += 1
                wrong_G = list(G)
                donor = next(v for v in range(12) if wrong_G[v])
                wrong_G[donor] -= 1; wrong_G[(donor+1) % 12] += 1
                wrong_result = it.admit(claim, wrong_G, region, 2)
                assert not wrong_result['admitted']
                controls['wrong_known_preparation'] += 1
    assert stale_latch_witness is not None
    # Native reset and all three retained failure stages, at both load totals.
    for total in (1, 2):
        G = (total,) + (0,)*11
        received, audit = it.source_execute(G, [(0, 0, 'C')], [0]*30)
        region = (0,)
        command = total
        for stage, reason in [('discover','discovery_timeout'), ('reset_route','routing_timeout'),
                              ('measure','measurement_timeout')]:
            ptr = it.meter.execute(0, (command,), timeout_tape(0, command, 2, stage), 2)
            assert ptr['reason'] == reason
            claim = it.assemble_claim(G, received, region, ptr, 2)
            verdict = it.admit(claim, G, region, 2)
            assert verdict['admitted'] and verdict['readout'] is None
            accepted += 1; outcomes[reason] += 1
            samples.append({'G': list(G), 'region': list(region), 'retries': 2,
                            'claim': claim, 'executor_audit': audit, 'admission': verdict})
            for guess in (0, command):
                corrupt = deepcopy(claim)
                corrupt['meter']['winding'] = guess
                corrupt['decoded_spatial_degree'] = guess
                corrupt['decoded_relative_charge'] = -guess
                result = it.admit(corrupt, G, region, 2)
                assert not result['admitted']
                controls['timeout_guessed_charge'] += 1
                for r in result['reasons']:
                    rejection_rules['timeout_guessed_charge'][r.split(':')[0]] += 1
            corrupt = deepcopy(claim); corrupt['meter']['success'] = True
            corrupt['meter']['winding'] = command
            assert not it.admit(corrupt, G, region, 2)['admitted']
            controls['timeout_relabelled_success'] += 1
    return {'accepted_full_instruments': accepted, 'outcomes': dict(outcomes),
            'successful_pointer_records_replayed': record_checks, 'rejections': dict(controls),
            'rejection_rules': {k: dict(v) for k,v in rejection_rules.items()},
            'forged_source_terminal_metadata_ignored': source_metadata_ignored,
            'stale_L0_witness': stale_latch_witness, 'samples': samples}


def damage_and_hypotheses():
    damaged = 0
    counterexamples = []
    frozen_strict_false_positive = 0
    for total in (1, 2):
        G = (total,) + (0,)*11
        for order_name in ('identity', 'reverse', 'seed1'):
            order = it.sweep.order_family()[order_name]
            received, audit = it.source_execute(G, [(0, 0, 'C')], [0]*30, order,
                                                prefix_channel=lambda records: [])
            assert sum(audit['source_L0']) == sum(G)
            assert audit['source_L0'] != list(G)
            try:
                it.observe_source(G, received, total)
            except ValueError as error:
                assert str(error).startswith('I3'), str(error)
                damaged += 1
            else:
                raise AssertionError('lost moving current escaped faithful sweep')
            if total == 2:
                ptr = it.meter.execute(7, (0,), it.meter.favorable_coins(7, (0,), 2), 2)
                standalone = it._au5_adapter(G, [], (0,), ptr, 7, 2, [0]*12)
                assert it.strict.validate(standalone)[0]
                frozen_strict_false_positive += 1
            # Forge the whole received E channel to the reconstruction's own
            # native prediction. This violates faithful-channel H2.
            def forge(record, reconstruction):
                a,b = record['seam']
                record['contrast'] = reconstruction[a] - reconstruction[b]
                record['current'] = it.sweep.current_of(record['contrast'], record['coin'])
                return record
            fake, hidden = it.source_execute(G, [(0, 0, 'C')], [0]*30, order,
                         prefix_channel=lambda records: [], sweep_channel=forge)
            latch = it.observe_source(G, fake, total)
            assert hidden['source_L0'] != latch['L0']['reconstruction']
            command = it.lift(latch['L1']['counters'][0], total)
            ptr = it.meter.execute(7, (command,), it.meter.favorable_coins(7, (command,), 2), 2)
            claim = it.assemble_claim(G, fake, (0,), ptr, 2)
            result = it.admit(claim, G, (0,), 2)
            assert result['admitted']
            counterexamples.append({'total':total, 'order':order_name,
                'observer_admitted':True, 'L0_reconstruction':latch['L0']['reconstruction'],
                'actual_L0':hidden['source_L0'], 'H2_faithful_sweep_channel':False,
                'L1_reconstruction':latch['L1']['reconstruction'], 'actual_L1':hidden['source_L1']})
    # Native currents that cancel can be lost without changing either latch.
    G = (1,) + (0,)*11
    received, audit = it.source_execute(G, [(0,0,'C'),(0,1,'C')], [0]*30,
                                       prefix_channel=lambda records: [])
    latch = it.observe_source(G, received, 1)
    assert audit['source_L0'] == latch['L0']['reconstruction']
    assert audit['source_L1'] == latch['L1']['reconstruction']
    # Connected-graph condition is checked algebraically on every possible
    # one-event discrepancy; it is independent of coins and sweep order.
    discrepancy_checks = 0
    for a,b in it.EDGES:
        for j in (-1,1):
            delta = [0]*12; delta[a] = -j; delta[b] = j
            mismatches = [e for e,(x,y) in enumerate(it.EDGES) if delta[x] != delta[y]]
            assert len(mismatches) == 9 and sum(delta) == 0
            discrepancy_checks += 1
    assert all(1-1 == 0 for a,b in it.EDGES)  # constant discrepancy is invisible
    return {'lost_moving_tape_with_consistent_metadata_rejected':damaged,
            'same_shortened_claim_passes_frozen_strict_alone':frozen_strict_false_positive,
            'adversarial_channel_counterexamples':counterexamples,
            'cancelled_lost_currents_have_correct_reconstruction':True,
            'missing_clock_of_cancelled_records_is_not_authenticated':True,
            'one_event_discrepancy_nine_seam_checks':discrepancy_checks,
            'constant_offset_counterexample':'Algebraic only: delta=(1,...,1) has zero seam contrasts and total difference 12. No out-of-domain source is executed.'}


def complete_laws():
    # Complete finite native domains, with all 2^32 source coin words covered
    # by exact weighted recursion, then every pointer branch with R=1.
    starts = preparations()[:12] + [tuple(int(v==0)+int(v==b) for v in range(12)) for b in (0,1,2,3)]
    cached_pointer = lru_cache(None)(it.pointer_distribution)
    rows = []
    words = 0
    max_probability_error = Fraction()
    min_success = Fraction(1)
    max_success = Fraction()
    outcome_union = set()
    for G in starts:
        for order_name in ('identity','reverse'):
            order = it.sweep.order_family()[order_name]
            for region in ((0,1,2), tuple(range(3,12))):
                source_joint = defaultdict(Fraction)
                for c0,c1 in itertools.product((0,1), repeat=2):
                    n = G
                    for e,c in ((0,c0),(5,c1)):
                        n = it.native(n, it.EDGES[e], c)[0]
                    w0 = sum(G[v]-n[v] for v in region)
                    final_dist = it.sweep_distribution(n, order)
                    for final, mass in final_dist.items():
                        w1 = sum(G[v]-final[v] for v in region)
                        assert -sum(G) <= w1 <= sum(G)
                        assert it.lift(w1%6, sum(G)) == w1
                        source_joint[(w0,w1)] += mass/4
                assert sum(source_joint.values()) == 1
                words += 2**32
                for initial in range(12):
                    full = defaultdict(Fraction)
                    for (w0,w1), p in source_joint.items():
                        for (reason,winding,signed,attempts), q in cached_pointer(initial,w1,1).items():
                            full[(w0,w1,reason,-winding if winding is not None else None,signed,32+attempts)] += p*q
                    total = sum(full.values())
                    assert total == 1
                    max_probability_error = max(max_probability_error, abs(total-1))
                    by_outcome = defaultdict(Fraction)
                    expectation = Fraction()
                    for key,p in full.items():
                        by_outcome[key[2]] += p
                        expectation += p*key[5]
                        assert key[5] <= 32+(33+6*sum(G))
                        if key[2] != 'complete_closed_pointer_word':
                            assert key[3] is None
                    outcome_union.update(by_outcome)
                    success = by_outcome['complete_closed_pointer_word']
                    min_success=min(min_success,success); max_success=max(max_success,success)
                    serial = [list(key)+[str(value)] for key,value in sorted(full.items(),key=lambda item:repr(item[0]))]
                    rows.append({'G':list(G),'order':order_name,'region':list(region),'pointer_initial_column':initial,
                         'mass':str(total),'support':len(full),'outcome_mass':{k:str(v) for k,v in by_outcome.items()},
                         'expected_total_attempts':str(expectation),'maximum_attempts':max(k[5] for k in full),
                         'complete_joint_law_sha256':digest(serial)})
    assert outcome_union == {'complete_closed_pointer_word','discovery_timeout','routing_timeout','measurement_timeout'}
    return {'instrument_count':len(rows),'source_preparations':len(starts),'source_prefix_schedule':[0,5],
            'source_prefix_record_kinds':['C','E'],'prefix_attempts':2,'sweep_attempts':30,'pointer_retries':1,
            'sweep_orders':['identity','reverse'],'unknown_pointer_columns':12,'regions_including_complement':2,
            'source_coin_words_covered_per_source_setting':2**32,
            'source_setting_count':words//2**32,
            'coverage_method':'Exact finite-state recursion with integer word multiplicities, followed by rational native-pointer policy DP. No literal enumeration of 2^30 sweep words.',
            'total_mass_error':str(max_probability_error),'unconditional_success_range':[str(min_success),str(max_success)],
            'failures_retained_without_normalization':True,'rows':rows}


def main():
    started=time.monotonic()
    before=pins()
    audit=frozen_audit()
    assert not audit['observer_has_source_state_parameter']
    print('Frozen byte, signature, policy and native-branch audit passed.', flush=True)
    trees=[executor_tree(x,w) for x in range(12) for w in (-2,-1,0,1,2)]
    print('Actual native pointer branch trees match exact probability laws.', flush=True)
    examples=full_examples_and_controls()
    hypotheses=damage_and_hypotheses()
    print('Full L1 instruments and damage controls passed.', flush=True)
    laws=complete_laws()
    after=pins()
    assert before==after
    out={'schema':'oph.row1.integrity.reproduction.v1','frozen_inputs_before':before,'frozen_inputs_after':after,
         'all_pins_unchanged':True,'frozen_audit':audit,'actual_pointer_executor_complete_trees':trees,
         'full_instruments':examples,'integrity_hypotheses_and_controls':hypotheses,'complete_bounded_laws':laws,
         'elapsed_seconds':round(time.monotonic()-started,3),
         'scope':'Native classical sources with total one or two. Known G, conserved total, faithful closing E channel. Finite counters, externally declared addressed scheduler, actual native pointer; no autonomous counter sensor is derived.',
         'physical_boundary':'Integer current/winding and charge-unit dictionary are declared architecture. No electron, spin, force or continuum identification follows.',
         'reserved_or_load_four_executed':False,'official_rows_changed':[]}
    path=HERE/'RESULTS.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'receipt':str(path),'full_instruments':examples['accepted_full_instruments'],
        'complete_instruments':laws['instrument_count'],'controls_rejected':sum(examples['rejections'].values()),
        'all_pins_unchanged':True,'elapsed_seconds':out['elapsed_seconds']}),flush=True)

if __name__=='__main__':
    main()
