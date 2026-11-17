#!/usr/bin/env python3
"""Independent current-law and whole-word reader; no producer imports."""
from fractions import Fraction
from pathlib import Path
import gzip
import hashlib
import json

P = Path(__file__).resolve().parent


def main():
    receipt = json.loads((P/'RECEIPT.json').read_text())
    checked = 0
    with gzip.open(P/'RECORDS.jsonl.gz', 'rt') as tape:
        for case in receipt['cases']:
            native = case['mode'] == 'native'
            accepted = nonzero = 0
            counts = [[0]*3 for _ in range(3)]
            covariance = [[0]*3 for _ in range(3)]
            for initial in range(2):
                for word in range(1 << 12):
                    left = 1-initial
                    js = []
                    for bit in range(12):
                        c = word // (1 << bit) % 2
                        j = left-c if native else 0
                        left -= j
                        assert left in (0, 1)
                        js.append(j)
                    xs = []
                    for slot in range(3):
                        segment = js[slot*4:slot*4+4]
                        xs.append((2*abs(segment[3])-1)*abs(segment[0]*segment[1]*segment[2]))
                    ok = sum(xs) == 0
                    expected = [case['mode'], initial, word, js, xs, ok]
                    assert json.loads(next(tape)) == expected
                    checked += 12
                    for i, x in enumerate(xs):
                        counts[i][x+1] += 1
                    if ok:
                        accepted += 1
                        nonzero += sum(x*x for x in xs) != 0
                        for i in range(3):
                            for j in range(3):
                                covariance[i][j] += xs[i]*xs[j]
            assert counts == case['coordinate_counts']
            assert accepted == case['accepted'] and nonzero == case['accepted_nonzero']
            assert [[str(Fraction(8*x, accepted)) for x in row] for row in covariance] == case['conditional_F_covariance']
        assert tape.read() == ''
    # Independently recompute each population directly from multinomial coefficients.
    from math import factorial
    for level in receipt['compressed_populations']:
        n = level['n']
        law = []
        for z in range(-n, n+1):
            coefficient = 0
            for minus in range(n+1):
                plus = minus+z
                zero = n-minus-plus
                if plus >= 0 and zero >= 0:
                    coefficient += factorial(n)//(factorial(minus)*factorial(plus)*factorial(zero))*14**zero
            law.append(coefficient)
        normalizer = m2 = m4 = 0
        for i, x in enumerate(law):
            for j, y in enumerate(law):
                k = 3*n-i-j
                if 0 <= k < len(law):
                    value = x*y*law[k]
                    normalizer += value
                    m2 += (i-n)**2*value
                    m4 += (i-n)**4*value
        assert str(normalizer) == level['normalizer_integer']
        assert Fraction(level['acceptance']) == Fraction(normalizer, 16**(3*n))
        assert Fraction(level['variance']) == Fraction(8*m2, n*normalizer)
        assert Fraction(level['fourth_moment']) == Fraction(64*m4, n*n*normalizer)
    result = {'pass': True, 'native_and_changed_law_events_checked': checked,
              'full_current_reader': True, 'independent_multinomial_population_check': True,
              'same_group': True, 'G6_pass': False,
              'receipt_sha256': hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest()}
    (P/'VERIFICATION.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
