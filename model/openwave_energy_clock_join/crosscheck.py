"""Internal opposing calculation: spectral projector and cosine bisection.

Imports neither check.py nor bounds.py. It shares the source equations and
the frozen graph, so this is an implementation crosscheck, not an external
review. Optional --compare reads primary output only after computing.
"""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import argparse
import json

import sympy as s

ROOT = Path(__file__).resolve().parents[2]
GRAPH = ROOT/'code/openwave_static_signed_source/evidence/openwave/shared/localized_rest_join/evidence/graph.json'


def oppose():
    data = json.loads(GRAPH.read_text())
    edges = data['seams']
    A = s.zeros(12)
    for a, b in edges:
        if a < 12 and b < 12:
            A[a, b] = A[b, a] = 1
    L = 5*s.eye(12)-A
    r = s.sqrt(5)
    high, low = 5+r, 5-r
    projector = (L*(L-6*s.eye(12))*(L-low*s.eye(12))/(high*(high-6)*(high-low))).applyfunc(s.simplify)
    assert (projector*projector-projector).applyfunc(s.simplify) == s.zeros(12)
    assert s.trace(projector) == 3
    external = [(a, b) for a, b in edges if a//12 != b//12]
    profiles, v = {}, []
    for cell in range(20):
        ports = tuple(sorted({p % 12 for pair in external for p in pair if p//12 == cell}))
        if ports not in profiles:
            block = projector[list(ports), list(ports)]
            profiles[ports] = (projector[:, list(ports)]*block.inv()*s.ones(3, 1)).applyfunc(s.simplify)
        v.extend(profiles[ports])
    # Edgewise action, rather than a full Laplacian matrix multiplication.
    out = [s.Integer(0)]*240
    for a, b in edges:
        delta = v[a]-v[b]
        out[a] += delta
        out[b] -= delta
    assert all(s.simplify(a-high*b) == 0 for a, b in zip(out, v))
    assert all(v[a] == v[b] == 1 for a, b in external)
    assert v[0] == -1 and s.simplify(v[5]-(r-2)) == 0

    grid = 10**40
    j = isqrt(5*grid*grid)
    rlo, rhi = F(j, grid), F(j+1, grid)
    t = ((F(209, 210), F(209, 210)),
         (1-(5+rhi)/1260, 1-(5+rlo)/1260))

    def cosine(x):
        up = sum(((-1)**n*x**(2*n)/factorial(2*n) for n in range(19)), F())
        return up-x**38/factorial(38), up

    def root(target):
        a, b = F(0), F(1)
        for _ in range(92):
            m = (a+b)/2
            lo, hi = cosine(m)
            if lo > target:
                a = m
            elif hi < target:
                b = m
            else:
                break
        return a, b

    def rate(x):
        d = 1-x
        lo = sum((d**n/n for n in range(1, 18)), F())
        return 1260*lo, 1260*(lo+d**18/(18*(1-d)))

    e = [(rate(hi)[0], rate(lo)[1]) for lo, hi in t]
    theta = [(root(hi)[0], root(lo)[1]) for lo, hi in t]
    k = (e[0][0]/theta[0][1], e[0][1]/theta[0][0])
    predicted = (k[0]*theta[1][0], k[1]*theta[1][1])
    gap = (e[1][0]-predicted[1], e[1][1]-predicted[0])
    assert gap[0] > F(13, 20)
    return {'method': 'cell spectral projector, edgewise source equation, log-power series, cosine bisection',
            'energy_intervals': [list(map(str, row)) for row in e],
            'phase_intervals': [list(map(str, row)) for row in theta],
            'second_energy_residual_interval': list(map(str, gap)),
            'second_mode_norm_squared': str(s.simplify(sum(x*x for x in v))),
            'positive_lower_residual_display': float(gap[0]),
            'independent_reviewer': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--compare', type=Path)
    args = parser.parse_args()
    result = oppose()
    if args.compare:
        primary = json.loads(args.compare.read_text())
        for key in ('energy_intervals', 'phase_intervals'):
            for a, b in zip(result[key], primary[key]):
                a, b = tuple(map(F, a)), tuple(map(F, b))
                assert max(a[0], b[0]) <= min(a[1], b[1])
        a, b = tuple(map(F, result['second_energy_residual_interval'])), tuple(map(F, primary['second_energy_residual_interval']))
        assert max(a[0], b[0]) <= min(a[1], b[1])
        result['primary_intervals_overlap_after_independent_computation'] = True
    print(json.dumps(result, indent=2))
