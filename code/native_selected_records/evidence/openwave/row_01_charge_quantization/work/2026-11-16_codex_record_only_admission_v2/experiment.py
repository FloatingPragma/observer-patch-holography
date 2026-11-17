"""Registered metadata isolation and unknown-count admission; retain every tape."""
from pathlib import Path
from fractions import Fraction
import gzip, hashlib, json, random, time
import numpy as np
from native_source import graph, prepare, source, CALLS
from observer import admit, decode
from legacy_decoder import observe
from legacy_selectors import selector
from continuation import columns, prefix
HERE = Path(__file__).resolve().parent
SEEDS = (72361000, 72362000)
PROGRESS = []

def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))

def digest(x):
    return hashlib.sha256(canonical(x).encode()).hexdigest()

def first_cell(bound, total):
    sign = 1 if total >= 0 else -1
    left = abs(total)
    q = []
    for _ in range(12):
        amount = min(bound, left)
        q.append(sign * amount)
        left -= amount
    empty = [i for i, x in enumerate(q) if x == 0]
    for a, b in zip(empty[::2], empty[1::2]):
        q[a] = bound
        q[b] = -bound
    require(sum(q) == total and max(map(abs, q)) <= bound, 'scientific check failed')
    return q

def legacy_blinding():
    cases = json.loads(gzip.decompress((HERE / 'LEGACY_RECORDS.json.gz').read_bytes()))
    baseline = {}
    trials = 0
    for index, c in enumerate(cases):
        r = c['record']
        answer = observe(c['n'], c['edges'], r['body'], r['probe'], 32)
        require(answer == r['reading'], 'scientific check failed')
        baseline[index] = answer
    for permutation in (list(range(len(cases))), list(reversed(range(len(cases)))), random.Random(9001).sample(range(len(cases)), len(cases))):
        for index in permutation:
            c = cases[index]
            r = dict(c['record'])
            r.update(seed=-100000 - index, total_input_audit=999999, initial_audit_only=['wrong'], filename=f'Q{index % 7}.json', run_order=-index)
            answer = observe(c['n'], c['edges'], r['body'], r['probe'], 32)
            require(answer == baseline[index], 'scientific check failed')
            trials += 1
    selectors = []
    for m in (4, 8):
        g = graph(2, m)
        for count in (2, 3, 5, 2 * m):
            a = selector(count, 8, g)
            altered = dict(g, seed=-777, total=123, filename='wrong', run_order=0)
            b = selector(count, 8, altered)
            require(a == b, 'scientific check failed')
            selectors.append({'M': m, 'count': count, 'output_sha256': digest(a), 'metadata_invariant': True})
    return {'records': len(cases), 'metadata_permutation_trials': trials, 'selectors': selectors, 'legacy_total_coded_seed_channel_retained': True, 'finding': 'The exact old observer does not consume seed/total/filename/order. Its harness nevertheless encoded total in seed; the new source uses both common public streams for every hidden input.'}

def run():
    for law in CALLS: CALLS[law] = 0
    PROGRESS.clear()
    started = time.perf_counter()
    legacy = legacy_blinding()
    rows = []
    prefixes = {}
    count = 0
    control_challenges = 0
    for bound in (2, 3, 4):
        for m in (4, 8):
            g = graph(bound, m)
            for seed in SEEDS:
                for total in range(-12 * bound, 12 * bound + 1):
                    initial = prepare(first_cell(bound, total), g)
                    for control in (False, True):
                        request, audit = source(initial, seed, control)
                        out = admit(request, g)
                        audit = audit()
                        count += len(audit['events'])
                        require(out.get('initial_core') == initial[:g['core_N']], 'scientific check failed')
                        require(out.get('observed_final') == audit['final'], 'scientific check failed')
                        result = {'B': bound, 'M': m, 'seed_audit_only': seed, 'total_audit_only': total, 'control': control, 'initial_audit_only': initial, 'controller': out, 'audit': audit}
                        PROGRESS.append(result)
                        if not control:
                            require(out['outcome'] == ('no_unit' if total == 0 else 'signed_target_unit'), 'scientific check failed')
                            if total:
                                sign = out['sign']
                                require(sign == (1 if total > 0 else -1), 'scientific check failed')
                                require(out['retained_core_total'] == total - sign, 'scientific check failed')
                                fine, R, cols = columns(audit['final'], g, sign)
                                count += len(cols)
                                result['consumer_global_columns'] = cols
                                key = f'{bound}:{m}:{sign}'
                                if key not in prefixes:
                                    prefixes[key] = prefix(audit['final'], g, sign)
                                    count += 40
                                require(R.tolist() == prefixes[key]['R'], 'scientific check failed')
                        else:
                            nonunit = max(out['post_probe_core']) - min(out['post_probe_core']) > 1
                            if nonunit:
                                control_challenges += 1
                                require(out['outcome'] == 'repair_failed', 'scientific check failed')
                            else:
                                require(out['outcome'] in ('no_unit', 'signed_target_unit'), 'scientific check failed')
                        rows.append(result)
    require(control_challenges >= 10, 'scientific control did not challenge enough unselected states')
    g = graph(2, 4)
    initial = prepare(first_cell(2, 1), g)
    request, audit = source(initial, SEEDS[0], forced_coin=0)
    timeout = admit(request, g)
    timeout_audit = audit()
    require(timeout['outcome'] == 'probe_timeout', 'scientific check failed')
    replay_order = random.Random(317).sample(range(len(rows)), len(rows))
    blind = 0
    for index in replay_order:
        row = rows[index]
        events = iter(row['audit']['events'])

        def replay(edge, stage):
            saved_stage, r = next(events)
            require(stage == saved_stage and list(edge) == r[:2], 'scientific check failed')
            return list(r)
        answer = admit(replay, graph(row['B'], row['M']))
        require(answer == row['controller'], 'scientific check failed')
        require(next(events, None) is None, 'scientific check failed')
        blind += 1
    conjugacy = []
    for bound in (2, 3, 4):
        for m in (4, 8):
            plus = np.array(prefixes[f'{bound}:{m}:1']['R'])
            minus = np.array(prefixes[f'{bound}:{m}:-1']['R'])
            permutation = [11 * p + (h ^ 1 if h < 10 else h) for p in range(3) for h in range(11)]
            require(np.array_equal(minus, plus[np.ix_(permutation, permutation)]), 'scientific check failed')
            conjugacy.append({'B': bound, 'M': m, 'exact_R_coin_conjugacy': True})
    bounds = {str(b): str(Fraction(2 * (12 * b - 1) + 6 * b * b * (12 * b - 1) + 12 * b, 2 ** 64)) for b in (2, 3, 4)}
    return {'candidate': 'A-native-record-only-admission-v2', 'passed': True, 'legacy_blinding': legacy, 'rows': rows, 'rows_count': len(rows), 'blind_replays': blind, 'complete_prefixes': prefixes, 'sign_conjugacy': conjugacy, 'control_unselected_challenges': control_challenges, 'literal_source_calls': dict(CALLS), 'timeout_control': {'controller': timeout, 'audit': timeout_audit}, 'inconclusive_upper_bounds': bounds, 'seconds': time.perf_counter() - started, 'external_G6': False, 'official_awards_added': 0}

def require(condition, message='scientific check failed'):
    if not condition:
        raise AssertionError(message)

def main():
    import argparse,sys,traceback
    parser=argparse.ArgumentParser();parser.add_argument('--record-to',type=Path)
    args=parser.parse_args();destination=args.record_to
    if destination is not None:destination.mkdir(parents=True,exist_ok=False)
    try:
        require(__debug__,'Optimized Python is unsupported: legacy source assertions must remain active')
        registry=HERE/('AMENDMENT_PREREG.json' if (HERE/'AMENDMENT_PREREG.json').exists() else 'PREREG.json')
        prereg=json.loads(registry.read_text())
        for n,h in prereg['files'].items():require(hashlib.sha256((HERE/n).read_bytes()).hexdigest()==h,n)
        if destination is None:
            receipt=json.loads((HERE/'records/RECEIPT.json').read_text())
            require(receipt['passed'] is True,'saved receipt reports scientific failure')
            require(hashlib.sha256((HERE/'records/EVIDENCE.json.gz').read_bytes()).hexdigest()==receipt['evidence_sha256'],'saved evidence custody mismatch')
        result=run();require(result['passed'] is True,'failed scientific verdict')
        if destination is None:
            old=json.loads(gzip.decompress((HERE/'records/EVIDENCE.json.gz').read_bytes()))
            require({k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in old.items() if k!='seconds'},'saved evidence mismatch')
            print(json.dumps({'passed':True,'read_only_replay':True,'rows':result['rows_count'],'seconds':result['seconds']}))
        else:
            data=gzip.compress((canonical(result)+'\n').encode(),mtime=0)
            with (destination/'EVIDENCE.json.gz').open('xb') as f:f.write(data)
            receipt={k:v for k,v in result.items() if k not in ('rows','complete_prefixes','timeout_control')}
            receipt['evidence_sha256']=hashlib.sha256(data).hexdigest()
            with (destination/'RECEIPT.json').open('x') as f:f.write(json.dumps(receipt,indent=2)+'\n')
            print(json.dumps(receipt,indent=2))
        return 0
    except Exception as error:
        failure={'passed':False,'error':repr(error),'traceback':traceback.format_exc(),'completed_or_current_rows':PROGRESS}
        if destination is not None:
            with (destination/'FAILURE.json').open('x') as f:f.write(json.dumps(failure,indent=2)+'\n')
        print(json.dumps({'passed':False,'error':repr(error)}));return 1

if __name__=='__main__':raise SystemExit(main())
