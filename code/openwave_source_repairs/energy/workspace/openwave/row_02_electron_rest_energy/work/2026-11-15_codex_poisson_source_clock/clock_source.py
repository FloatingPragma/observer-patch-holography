"""Native-count Poisson clock and full original-source compact generators."""
from collections import Counter, defaultdict
from fractions import Fraction as F
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
from native_event import native

HERE = Path(__file__).resolve().parent
GRAPH = json.loads((HERE / 'graph.json').read_text())
EDGES = [tuple(e) for e in GRAPH['seams']]
K = 2 * len(EDGES)
CALLS = Counter()


def call(x, edge, c, kind):
    CALLS[kind] += 1
    return native(x, edge, c)


def exp_interval(z, terms=160):
    z = F(z)
    assert 0 <= z <= 1 and terms % 2 == 0
    value, term, before = F(1), F(1), None
    for n in range(1, terms + 1):
        before = value
        term *= -z / n
        value += term
    assert before <= value
    return before, value


ELO, EHI = exp_interval(1)
CDF = []
partial = F()
for n in range(33):
    partial += F(1, math.factorial(n))
    CDF.append((ELO * partial, min(F(1), EHI * partial)))
TAIL = F(34, 33 * math.factorial(33))
CLOCK_FAILURE = F(68, 2**128) + 68 * (EHI - ELO) + TAIL


def classify(prefix, bits):
    low, high = F(prefix, 1 << bits), F(prefix + 1, 1 << bits)
    previous_upper = F(0)
    for n, (clo, chi) in enumerate(CDF):
        if low >= previous_upper and high <= clo:
            return n
        previous_upper = chi
    return -2 if low > CDF[-1][1] else -1


def coin(aux, rng, mode='native'):
    c = rng.randrange(2)
    if mode == 'native':
        out, (_, j) = call(aux, (0, 1), c, 'clock_auxiliary')
    else:
        out, j = aux, 0
    return out, abs(j), [c, j, out.index(1)]


def sample_count(aux, rng, mode='native'):
    prefix, tape = 0, []
    for bits in range(1, 129):
        aux, bit, row = coin(aux, rng, mode)
        tape.append(row)
        prefix = 2 * prefix + bit
        value = classify(prefix, bits)
        if value >= 0 or value == -2:
            return aux, value, tape
    return aux, -1, tape


def address(aux, rng):
    tape = []
    for _ in range(64):
        label, trial = 0, []
        for _ in range(10):
            aux, bit, row = coin(aux, rng)
            trial.append(row)
            label = 2 * label + bit
        tape.append(trial)
        if label < len(EDGES):
            return aux, label, tape
    return aux, None, tape


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n').encode()


def source_modes():
    bases = json.loads((HERE / 'DEVELOPMENT_BASIS.json').read_text())['basis_by_cell']
    assert [len(b) for b in bases] == [2, 2]
    singles = []
    for cell, base in enumerate(bases):
        for mode, v in enumerate(base):
            singles.append({'name': f'S{cell}{mode}', 'kind': 'single',
                            'inputs': [[cell * 12 + i, value] for i, value in enumerate(v) if value]})
    pairs = []
    for i in range(2):
        for j in range(2):
            pairs.append({'name': f'P{i}{j}', 'kind': 'pair',
                'inputs': [[a * 256 + b, av * bv]
                           for a, av in singles[i]['inputs'] for b, bv in singles[2 + j]['inputs']]})
    result = []
    for mode in singles + pairs:
        counts, records = defaultdict(int), []
        for code, amplitude in mode['inputs']:
            initial = [0] * 240
            if mode['kind'] == 'single':
                initial[code] = 1
            else:
                p, n = divmod(code, 256)
                initial[p], initial[n] = 1, -1
            initial = tuple(initial)
            for e, edge in enumerate(EDGES):
                for c in (0, 1):
                    out, (_, current) = call(initial, edge, c, 'source_mode_columns')
                    occupied = [(i, v) for i, v in enumerate(out) if v]
                    if mode['kind'] == 'single':
                        assert len(occupied) == 1 and occupied[0][1] == 1
                        dest = occupied[0][0]
                    elif not occupied:
                        dest = -1
                    else:
                        assert sorted(v for _, v in occupied) == [-1, 1]
                        p = next(i for i, v in occupied if v == 1)
                        n = next(i for i, v in occupied if v == -1)
                        dest = p * 256 + n
                    counts[dest] += amplitude
                    records.append([code, e, c, current, dest])
        counts = {q: v for q, v in counts.items() if v}
        eigen_numerator = K - (6 if mode['kind'] == 'single' else 12)
        wanted = {q: eigen_numerator * v for q, v in mode['inputs']}
        assert counts == wanted
        result.append({**mode, 'branches': records,
                       'mean_numerator': sorted(counts.items()), 'eigen_numerator': eigen_numerator,
                       'hold_generator_eigenvalue': 0})
    return result


def clock_tree():
    level, rows = {0: (1, 0)}, []
    for depth in range(16):
        children = {}
        for prefix, aux in sorted(level.items()):
            emitted = set()
            for c in (0, 1):
                out, (_, j) = call(aux, (0, 1), c, 'clock_prefix_columns')
                bit = abs(j)
                emitted.add(bit)
                child = 2 * prefix + bit
                assert child not in children
                children[child] = out
                rows.append([depth, prefix, aux.index(1), c, j, child, out.index(1)])
            assert emitted == {0, 1}
        level = children
    leaves = [classify(p, 16) for p in range(1 << 16)]
    counts = Counter(leaves)
    fail = F(sum(v for k, v in counts.items() if k < 0), 1 << 16)
    assert fail <= F(34, 1 << 16)
    for n in range(33):
        prob_lo, prob_hi = ELO / math.factorial(n), EHI / math.factorial(n)
        resolved = F(counts[n], 1 << 16)
        assert resolved <= prob_hi and resolved + fail + TAIL >= prob_lo
    factors = []
    for a, target_exp in ((F(209, 210), F(1, 210)), (F(104, 105), F(1, 105))):
        lower = sum(F(count, 1 << 16) * a**n for n, count in counts.items() if n >= 0)
        expected_lo, expected_hi = exp_interval(target_exp)
        assert lower <= expected_lo <= expected_hi <= lower + fail
        factors.append({'source_eigenvalue': str(a), 'retained_factor': str(lower),
                        'unresolved_mass': str(fail), 'ideal_factor_interval': list(map(str, (expected_lo, expected_hi)))})
    return {'native_prefix_branches': rows, 'leaf_outcomes': leaves,
            'leaf_counts': {str(k): v for k, v in sorted(counts.items())},
            'failure_mass': str(fail), 'transfer_factors': factors}


def actual_histories():
    result = []
    for index, (kind, port) in enumerate((kind, port) for kind in ('single', 'pair') for port in range(4, 8)):
        x = [0] * 240
        x[port] = 1
        if kind == 'pair':
            x[16] = -1
        initial = x[:]
        aux, rng, ticks = (1, 0), random.Random(72331000 + index), []
        for tick in range(64):
            aux, n, clock_records = sample_count(aux, rng)
            item = {'count': n, 'clock_records': clock_records, 'events': []}
            ticks.append(item)
            if n < 0:
                break
            for _ in range(n):
                aux, e, address_records = address(aux, rng)
                if e is None:
                    item['events'].append({'address_records': address_records, 'source': None})
                    break
                c = rng.randrange(2)
                x, (_, j) = call(x, EDGES[e], c, 'clocked_source_events')
                item['events'].append({'address_records': address_records, 'source': [e, c, j]})
        complete = len(ticks) == 64 and all(t['count'] >= 0 and len(t['events']) == t['count']
                        and all(e['source'] is not None for e in t['events']) for t in ticks)
        result.append({'kind': kind, 'port': port, 'seed': 72331000 + index,
                       'initial_audit_only': initial, 'ticks': ticks, 'complete': complete,
                       'final_source_audit_only': x, 'final_auxiliary': aux})
    assert all(r['complete'] for r in result)
    seen = {t['count'] for r in result for t in r['ticks']}
    assert 0 in seen and any(n >= 2 for n in seen)
    aux, rng, hold = (1, 0), random.Random(72331999), []
    for _ in range(16):
        aux, n, tape = sample_count(aux, rng, 'hold')
        hold.append({'count': n, 'records': tape})
        assert n == 0
    return result, hold


def log_one_minus(z):
    z, value = F(z), F()
    for n in range(1, 33):
        value += z**n / n
    return value, value + z**33 / (33 * (1 - z))


def rare_event():
    result = []
    for denominator in (1, 16, 256, 4096):
        epsilon = F(1, denominator)
        e1 = [K / epsilon * v for v in log_one_minus(6 * epsilon / K)]
        e2 = [K / epsilon * v for v in log_one_minus(12 * epsilon / K)]
        excess = [e2[0] - 2 * e1[1], e2[1] - 2 * e1[0]]
        upper = 36 * epsilon / (K - 12 * epsilon)
        assert 0 < excess[0] <= excess[1] < upper
        result.append({'epsilon': str(epsilon), 'single': list(map(str, e1)),
                       'pair': list(map(str, e2)), 'excess': list(map(str, excess)),
                       'excess_upper': str(upper)})
    return result


def calculate():
    CALLS.clear()
    assert GRAPH['ports'] == 240 and K == 1260
    modes = source_modes()
    tree = clock_tree()
    histories, hold = actual_histories()
    rare = rare_event()
    fixed_gap = F(209, 210)**2 - F(104, 105)
    assert fixed_gap == F(1, 44100)
    assert CLOCK_FAILURE < F(1, 10**30)
    evidence = {'source_modes': modes, 'clock_tree': tree, 'clock_histories': histories,
                'hold_auxiliary': hold, 'rare_event': rare,
                'exp_minus_one_interval': [str(ELO), str(EHI)],
                'poisson_count_cap_tail_upper': str(TAIL)}
    gates = {'source_single_and_four_pair_residuals_zero': len(modes) == 8,
             'complete_16bit_native_prefix_tree': len(tree['native_prefix_branches']) == 131070,
             'prefix_probabilities_and_transfer_enclosed': True,
             'all_eight_clock_histories_complete': all(r['complete'] for r in histories),
             'counts_zero_and_multiple_observed': True,
             'hold_clock_counts_zero': all(r['count'] == 0 for r in hold),
             'fixed_count_nonadditivity_retained': fixed_gap == F(1, 44100),
             'rare_event_positive_excess_and_limit_bound': len(rare) == 4,
             'finite_128bit_clock_error_bounded': CLOCK_FAILURE < F(1, 10**30)}
    receipt = {'candidate': 'A-native-poisson-event-clock-v1', 'gates': gates,
               'all_gates_pass': all(gates.values()), 'literal_native_calls_by_kind': dict(CALLS),
               'literal_native_calls': sum(CALLS.values()), 'native_single_generator_values': [6] * 4,
               'native_pair_generator_values': [12] * 4, 'poisson_pair_equals_sum_exact_ideal': True,
               'clock_tick_duration': '1/1260', 'clock_mean_count': 1,
               'finite_prefix_failure_mass': tree['failure_mass'],
               'per128bit_clock_sampler_failure_upper': str(CLOCK_FAILURE),
               'per_source_address_failure_upper': str(F(197, 512)**64),
               'fixed_one_attempt_factor_gap': str(fixed_gap),
               'rare_event_finite_excess_is_zero': False,
               'physical_rest_energy_identified': False, 'coherent_UL_execution_in_this_packet': False,
               'full_row_solved': False, 'source_pin_changed': False}
    return receipt, evidence


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--record', action='store_true')
    args = parser.parse_args()
    freeze = json.loads((HERE / 'PREREG.json').read_text())
    for name, digest in freeze['sha256'].items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    receipt, evidence = calculate()
    data = gzip.compress(encoded(evidence), mtime=0)
    receipt['evidence_sha256'] = hashlib.sha256(data).hexdigest()
    receipt['prereg_sha256'] = hashlib.sha256((HERE / 'PREREG.json').read_bytes()).hexdigest()
    if args.record:
        for name, raw in [('RECEIPT.json', encoded(receipt)), ('RECORDS.json.gz', data)]:
            with (HERE / name).open('xb') as f:
                f.write(raw)
    else:
        assert (HERE / 'RECORDS.json.gz').read_bytes() == data
        assert json.loads((HERE / 'RECEIPT.json').read_text()) == receipt
    print(json.dumps(receipt, indent=2))
    if not receipt['all_gates_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
