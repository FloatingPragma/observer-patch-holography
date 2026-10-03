"""Candidate producer: use the committed operator, execute every local update."""
import json
import importlib.util
import sys
from fractions import Fraction as F
from pathlib import Path
from .format import digest
from .observables import exp_bounds

ROOT = Path(__file__).resolve().parents[2]
# Legacy callers also import source_scalar_execution.py as a top-level module.
# Load this exact producer dependency without changing their module namespace.
_spec = importlib.util.spec_from_file_location('_paired_gravity_producer_algebra',
    ROOT/'code/source_scalar_execution/scalar_execution_algebra.py')
_algebra = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _algebra
_spec.loader.exec_module(_algebra)
Q, parse, rational_bounds = _algebra.Q, _algebra.parse, _algebra.rational_bounds
DEN = 2**36
SWEEPS = 32
GAMMAS = (-1, 0, 1, 2)


def data():
    packet = json.loads((ROOT/'code/source_scalar_execution/source_scalar_execution_receipt.json').read_text())
    rows = [dict((j, parse(a)) for j, a in row) for row in packet['action_rows_Qphi']]
    for i, row in enumerate(rows):
        row[i] -= 1
    return rows, [parse(a) for a in packet['mass_Qphi']]


def floor_field(x):
    # Floating arithmetic proposes a nearby integer; exact signs certify it.
    from math import sqrt, floor
    k = floor((float(x.a)+float(x.b)*(1+sqrt(5))/2)*DEN)
    while (x-F(k, DEN)).sign() < 0:
        k -= 1
    while (x-F(k+1, DEN)).sign() >= 0:
        k += 1
    return F(k, DEN)


def nonlinear_floor(diagonal, off, strength):
    from math import sqrt, exp, floor
    number = lambda z: float(z.a)+float(z.b)*(1+sqrt(5))/2
    a, b = number(diagonal), number(off)
    x = max(0., (float(strength)-b)/a)
    for _ in range(8):
        x -= (a*x+b-float(strength)*exp(-x))/(a+float(strength)*exp(-x))
    def sign(n):
        lo, hi = exp_bounds(-F(n, DEN))
        f0, f1 = diagonal*F(n, DEN)+off-strength*hi, diagonal*F(n, DEN)+off-strength*lo
        if f0.sign() > 0:
            return 1
        if f1.sign() < 0:
            return -1
        raise ValueError('unresolved nonlinear rounding sign')
    n = floor(x*DEN)
    while sign(n) > 0:
        n -= 1
    while sign(n+1) < 0:
        n += 1
    return F(n, DEN)


def execute(strength, reverse, seed, nonlinear):
    rows, mass = data()
    order = list(range(64))[::(-1 if reverse else 1)]
    u, w = [F(0)]*64, [F(0)]*64
    w[42] = strength/16 if seed else F(0)
    versions = [0]*64
    chain = {g: '0'*64 for g in GAMMAS}
    projection = '0'*64
    events = reads = 0
    def record(i, inputs, new_u, new_w):
        nonlocal projection, events, reads
        event = dict(site=i, version=versions[i]+1, reads=inputs,
                     u=str(new_u), w=str(new_w))
        projection = digest(dict(parent=projection, **event))
        for g in GAMMAS:
            v_reads = [str(g*F(a[2])+F(a[3])) for a in inputs]
            chain[g] = digest(dict(parent=chain[g], gamma=g, v=str(g*new_u+new_w), v_reads=v_reads, **event))
        versions[i] += 1
        events += 1
        reads += len(inputs)
    for i in order:
        record(i, [], u[i], w[i])
    for sweep in range(SWEEPS):
        for i in order:
            inputs = [[j, versions[j], str(u[j]), str(w[j])] for j in sorted(rows[i])]
            rhs = Q(strength if i == 21 else 0)
            off = sum((a*u[j] for j, a in rows[i].items() if j != i), Q())
            du = (rhs-off)/rows[i][i]
            dw = -sum((a*w[j] for j, a in rows[i].items() if j != i), Q())/rows[i][i]
            nu = nonlinear_floor(rows[i][i], off, strength) if nonlinear and i == 21 else floor_field(du)
            nw = floor_field(dw)
            record(i, inputs, nu, nw)
            u[i], w[i] = nu, nw
    residuals = []
    for values, source in ((u, True), (w, False)):
        maximum, maximum_lower = Q(), Q()
        for i, row in enumerate(rows):
            raw = sum((a*values[j] for j, a in row.items()), Q())
            forcing = (exp_bounds(-values[i]) if nonlinear else (F(1), F(1))) if source and i == 21 else (F(0), F(0))
            left, right = raw-strength*forcing[1], raw-strength*forcing[0]
            absolutes = [-a if a.sign() < 0 else a for a in (left, right)]
            low, high = absolutes if (absolutes[0]-absolutes[1]).sign() <= 0 else absolutes[::-1]
            if left.sign() <= 0 <= right.sign():
                low = Q()
            if (high-maximum).sign() > 0:
                maximum = high
            if (low-maximum_lower).sign() > 0:
                maximum_lower = low
        residuals.append([rational_bounds(maximum_lower/8, 10**16)[0], rational_bounds(maximum/8, 10**16)[1]])
    return dict(strength=str(strength), reverse=reverse, seed=seed, nonlinear=nonlinear, sweeps=SWEEPS, denominator=DEN,
                u=list(map(str, u)), w=list(map(str, w)), projection=projection,
                chains={str(g): h for g, h in chain.items()}, events=events, reads=reads,
                error_intervals=residuals)


def build():
    return [execute(s, reverse, seed, nonlinear) for s in (F(1, 32), F(1, 16))
            for reverse in (False, True) for seed in (False, True) for nonlinear in (False, True)]
