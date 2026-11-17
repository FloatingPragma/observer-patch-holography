#!/usr/bin/env python3
"""Exhaustive native-current compatibility witness, not a spacetime simulation."""
from fractions import Fraction as F
from pathlib import Path
import gzip
import hashlib
import importlib.util
import json
import math
import sys

P = Path(__file__).resolve().parent
sys.dont_write_bytecode = True


def population(weights):
    # The coordinate law is read from complete native records below.
    levels = []
    distribution = [1]
    for n in range(1, 129):
        following = [0] * (len(distribution)+2)
        for i, count in enumerate(distribution):
            for j, weight in enumerate(weights):
                following[i+j] += count*weight
        distribution = following
        if n not in (1, 2, 4, 8, 16, 32, 64, 128):
            continue
        norm = second = fourth = 0
        for i, wi in enumerate(distribution):
            zi = i-n
            for j, wj in enumerate(distribution):
                k = 3*n-i-j
                if 0 <= k <= 2*n:
                    count = wi*wj*distribution[k]
                    norm += count
                    second += zi*zi*count
                    fourth += zi**4*count
        var = F(8*second, n*norm)
        moment4 = F(64*fourth, n*n*norm)
        levels.append({'n': n, 'normalizer_integer': str(norm),
                       'acceptance': str(F(norm, sum(weights)**(3*n))),
                       'variance': str(var), 'fourth_moment': str(moment4)})
    return levels


def main():
    pin = P/'inputs/native.py'
    assert hashlib.sha256(pin.read_bytes()).hexdigest() == 'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
    spec = importlib.util.spec_from_file_location('compatible_native_pin', pin)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    cases = []
    native_calls = 0
    with (P/'RECORDS.jsonl.gz').open('wb') as stream:
        with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as compressed:
            for mode in ('native', 'hold_mutation'):
                accepted = nonzero = norm2 = 0
                marginals = [[0, 0, 0] for _ in range(3)]
                moments = [[0]*3 for _ in range(3)]
                for start in (0, 1):
                    for word in range(4096):
                        loads = tuple(int(i == start) for i in range(12))
                        currents = []
                        for bit in range(12):
                            coin = (word >> bit) & 1
                            if mode == 'native':
                                loads, record = module.native(loads, (0, 1), coin)
                                native_calls += 1
                                assert record[0] == coin
                                current = record[1]
                            else:
                                current = 0
                            currents.append(current)
                        values = []
                        for i in range(3):
                            gate = all(abs(j) == 1 for j in currents[4*i:4*i+3])
                            value = (1 if abs(currents[4*i+3]) else -1) if gate else 0
                            values.append(value)
                            marginals[i][value+1] += 1
                        success = sum(values) == 0
                        if success:
                            accepted += 1
                            nonzero += any(values)
                            norm2 += sum(x*x for x in values)
                            for i in range(3):
                                for j in range(3):
                                    moments[i][j] += values[i]*values[j]
                        line = [mode, start, word, currents, values, success]
                        compressed.write((json.dumps(line, separators=(',', ':'))+'\n').encode())
                cases.append({'mode': mode, 'complete_words': 8192, 'preparations': 8192,
                              'event_calls': 98304, 'accepted': accepted,
                              'accepted_nonzero': nonzero, 'accepted_square_sum': norm2,
                              'coordinate_counts': marginals,
                              'conditional_F_covariance': [[str(F(8*x, accepted)) for x in row] for row in moments]})
    assert native_calls == 98304
    native, control = cases
    gcd = math.gcd(*native['coordinate_counts'][0])
    weights = [v//gcd for v in native['coordinate_counts'][0]]
    exact = population(weights)
    final = exact[-1]
    checks = {
        'all_marginals_from_native_currents': all(x == [512, 7168, 512] for x in native['coordinate_counts']),
        'exact_native_compatibility_probability': F(native['accepted'], 8192) == F(707, 1024),
        'all_nonzero_accepted_words_retained': native['accepted_nonzero'] == 168,
        'finite_native_conditional_covariance': native['conditional_F_covariance'] == [[str(F(112 if i == j else -56, 707)) for j in range(3)] for i in range(3)],
        'endpoint_hold_kills_flux_variance': control['accepted'] == 8192 and control['accepted_nonzero'] == 0,
        'record_derived_ternary_weights': weights == [1, 14, 1],
        'compressed_n128_variance_error_below_002': abs(F(final['variance'])-F(2, 3)) < F(1, 50),
        'compressed_n128_fourth_error_below_01': abs(F(final['fourth_moment'])-F(4, 3)) < F(1, 10),
    }
    result = {'candidate': 'A-native-compatible-flux-v1', 'all_pass': all(checks.values()),
              'checks': checks, 'cases': cases, 'literal_native_calls': native_calls,
              'changed_law_event_calls': 98304, 'actual_new_preparations': 16384,
              'rank_one_integer_constraint': [1, 1, 1], 'compressed_populations': exact,
              'population_is_actual_long_tape': False, 'full_Maxwell_execution': False,
              'G6_pass': False, 'records_sha256': hashlib.sha256((P/'RECORDS.jsonl.gz').read_bytes()).hexdigest()}
    (P/'RECEIPT.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'all_pass': result['all_pass'], 'native_calls': native_calls,
                      'native_accepts': native['accepted'], 'native_nonzero_accepts': native['accepted_nonzero'],
                      'n128_variance': float(F(final['variance'])), 'n128_fourth': float(F(final['fourth_moment']))}))
    return 0 if result['all_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
