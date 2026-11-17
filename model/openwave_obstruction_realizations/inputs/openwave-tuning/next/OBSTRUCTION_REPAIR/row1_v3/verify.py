"""New fixed-policy evidence; no writes to the reviewed v2 package."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import random
import sys

sys.dont_write_bytecode = True
import observer

HERE = Path(__file__).resolve().parent
V2 = HERE.parent/'row1_v2'
EDGES = observer.v2.EDGES


def states():
    out = [(0,)*12]
    for a in range(12):
        for b in range(12):
            if a != b:
                n = [0]*12
                n[a], n[b] = 1, -1
                out.append(tuple(n))
    return out


def step(n, e, coin):
    a, b = EDGES[e]
    total = n[a]+n[b]
    m = list(n)
    m[a] = total//2+coin*(total % 2)
    m[b] = total-m[a]
    return tuple(m), n[a]-m[a]


def sweep(initial, coins, hidden=None):
    n, records = initial, []
    for k, (e, c) in enumerate(zip(observer.ORDER, coins)):
        if hidden is not None and k == hidden[0]:
            n, _ = step(n, hidden[1], hidden[2])
        a, b = EDGES[e]
        contrast = n[a]-n[b]
        n, j = step(n, e, c)
        records.append({'kind': 'E', 'seam': e, 'coin': c, 'current': j, 'contrast': contrast})
    return records, n


def coin_stream(rng):
    while True:
        yield rng.randrange(2)


def main():
    rng = random.Random(411721)
    coverage, outcomes = Counter(), Counter()
    examples = []
    later, middle = Counter(), Counter()
    counterexamples = []
    for initial in states():
        for repeat in range(4):
            coins = [rng.randrange(2) for _ in observer.ORDER]
            records, final = sweep(initial, coins)
            coverage['registered_sweeps'] += 1
            # Regions are selected by a declared independent diagnostic policy:
            # every singleton, its complement, and the complete carrier.
            regions = [[v] for v in range(12)]+[[u for u in range(12) if u != v] for v in range(12)]+[list(range(12))]
            for region in regions:
                got = observer.observe(records, region)
                charge = sum(final[v] for v in region)
                assert got['charge'] == charge and got['final_relative_state'] == list(final)
                assert not got['channel_completeness_authenticated']
                coverage['regional_reads'] += 1
                coverage['regional_charge_'+str(charge)] += 1
                x = rng.randrange(12)
                pointer = observer.v2.original.meter.execute(x, (-charge,), coin_stream(rng), 2)
                verdict = observer.admit(records, region, pointer)
                assert verdict['admitted'] and x in verdict['compatible_pointer_initials']
                assert verdict['readout'] == (charge if pointer['success'] else None)
                outcomes[verdict['reason']] += 1
                coverage['effect_charge_'+str(charge)] += 1
                if pointer['success']:
                    coverage['complete_readout_'+str(charge)] += 1
                if charge and pointer['success'] and len(examples) < 6:
                    examples.append({'records': records, 'region': region, 'pointer': pointer, 'verdict': verdict})
            bad = records[1:]+records[:1]
            assert not observer.admit(bad, [0], pointer)['admitted']
            # An actual extra native event, not an invented endpoint array.
            for e in range(30):
                for c in (0, 1):
                    altered, current = step(final, e, c)
                    if not current:
                        continue
                    region = [EDGES[e][0]]
                    before = observer.observe(records, region)
                    after = observer.observe(records, region)
                    assert before == after and before['charge'] != altered[region[0]]
                    later['actual_changed_native_events'] += 1
                    later['undetected_changed_charge'] += 1
                    if len(counterexamples) < 2:
                        counterexamples.append({'kind': 'after', 'records': records, 'hidden_seam': e,
                            'hidden_coin': c, 'region': region, 'read_charge': after['charge'],
                            'true_charge_after_hidden_event': altered[region[0]]})
            # A concealed native event during a genuine otherwise faithful sweep.
            for _ in range(12):
                hidden = (rng.randrange(1, 30), rng.randrange(30), rng.randrange(2))
                changed_records, changed_final = sweep(initial, coins, hidden)
                for v in range(12):
                    try:
                        got = observer.observe(changed_records, [v])
                    except ValueError:
                        middle['rejected'] += 1
                    else:
                        right = got['charge'] == changed_final[v]
                        middle['admitted_correct' if right else 'admitted_wrong'] += 1
                        if not right and not any(x['kind'] == 'during' for x in counterexamples):
                            counterexamples.append({'kind': 'during', 'initial': initial, 'coins': coins,
                                'hidden_event': hidden, 'records': changed_records, 'region': [v],
                                'read_charge': got['charge'], 'true_charge': changed_final[v]})
    laws = json.loads((V2/'EXACT_LAWS.json').read_text())
    ptr = {(r['initial_pointer'], r['degree']): {k: F(v) for k, v in r['probabilities'].items()}
           for r in laws['pointer_laws']}
    normalizations, survival = [], []
    for row in laws['complete_source_sweep_laws']:
        initial = row['initial_relative_state']
        dist = [(n, F(count, row['denominator'])) for n, count in row['final_counts']]
        assert sum(p for n, p in dist) == 1
        if any(initial):
            a, b = initial.index(1), initial.index(-1)
            adjacent = (a, b) in EDGES or (b, a) in EDGES
            survival.append({'initial': initial, 'adjacent': adjacent,
                'survival_probability': str(sum(p for n, p in dist if any(n)))})
        for region in ([0], list(range(1, 12)), list(range(12))):
            for x in range(12):
                totals = Counter()
                for n, p in dist:
                    degree = -sum(n[v] for v in region)
                    for outcome, q in ptr[x, degree].items():
                        totals[outcome] += p*q
                total = sum(totals.values())
                assert total == 1
                normalizations.append({'initial': initial, 'region': region, 'pointer': x,
                    'outcomes': {k: str(v) for k, v in sorted(totals.items())}, 'total': str(total)})
    survival_bounds = {}
    for adjacent in (True, False):
        values = [F(r['survival_probability']) for r in survival if r['adjacent'] == adjacent]
        survival_bounds[str(adjacent)] = {'count': len(values), 'minimum': str(min(values)), 'maximum': str(max(values))}
    assert coverage['complete_readout_1'] and coverage['complete_readout_-1']
    assert later['undetected_changed_charge'] and middle['admitted_wrong']
    result = {'schema': 'oph.row1.fixed_policy.v3', 'registered_order': observer.ORDER,
        'coverage': dict(coverage), 'pointer_outcomes': dict(outcomes),
        'computed_joint_normalizations': len(normalizations),
        'all_computed_normalizations_one': all(r['total'] == '1' for r in normalizations),
        'survival_bounds': survival_bounds, 'post_sweep_hidden_native_events': dict(later),
        'during_sweep_hidden_native_events': dict(middle),
        'channel_completeness_authenticated': False, 'physical_charge_identification_solved': False,
        'official_rows_changed': [],
        'consumed_laws_sha256': hashlib.sha256((V2/'EXACT_LAWS.json').read_bytes()).hexdigest()}
    for name, data in [('RESULT.json', result), ('JOINT_LAWS.json', normalizations),
                       ('SURVIVAL.json', survival), ('EXAMPLES.json', examples),
                       ('CHANNEL_COUNTEREXAMPLES.json', counterexamples)]:
        (HERE/name).write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
