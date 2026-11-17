"""Finite confidence from one continuing local signed-current history.

The stopped twelve-class source is a proof coordinate, never observer data.
For each visible transition reward h solve (I-Q)g = Qh-mu exactly. Then
sum(h-mu) = g(X_0)-g(X_n)+sum D_k, with conditionally centered D_k and
conditional range width at most W. If R=osc(g), bounded-difference iteration
gives P(|mean(h)-mu| > (R+W*t)/n) <= 2 exp(-2*t*t/n).
Using integer t=ceil(sqrt(n*(s+1)/2)) bounds this by 2**(-s), without floats.

The horizon n is fixed in advance in OWN INCIDENT TICKS. No stationary
preparation, reset, independent histories, initial-bit disclosure, or global
attempt counter is required. C/E complete current laws are population laws
also for coherent inputs and untouched references. Bit moments are decoded
retrospectively; an all-zero history returns [0,1] for those moments. Coverage
is unconditional, not conditioned on resolution and not valid for arbitrary
optional stopping. Uniform IID seam selection and fair placements remain.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import random

import sympy as sp

from . import bit_readout as bits
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('trajectory_results.json')
NAMES = ('occupation', 'adjacent_product', 'movement')


def _integer(value, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError('an original integer in the stated range is required')
    return value


@lru_cache(maxsize=12)
def source(port=0):
    """Exact stationary means and martingale certificates on actual source."""
    native, blocks, _, marked = bits.source(port)
    d = len(blocks)
    q = sum(marked.values(), sp.zeros(d))
    pi = sp.Matrix([[sp.Rational(len(b), len(native.states)) for b in blocks]])
    one = sp.ones(d, 1)
    if q*one != one or pi*q != pi or (sp.eye(d)-q).rank() != d-1:
        raise ArithmeticError('stopped native chain lost its unique stationary law')
    f = [native.states[b[0]][port] for b in blocks]
    rewards = {
        'occupation': sp.Matrix(d, d, lambda i, j: f[i]),
        'adjacent_product': sp.Matrix(d, d, lambda i, j: f[i]*f[j]),
        'movement': sp.Matrix(d, d, lambda i, j: abs(f[j]-f[i])),
    }
    inverse = (sp.eye(d)-q+one*pi).inv(method='DM')
    certificates = {}
    for name, h in rewards.items():
        r = sp.Matrix([sum(q[i, j]*h[i, j] for j in range(d)) for i in range(d)])
        mu = (pi*r)[0]
        g = inverse*(r-mu*one)
        if (sp.eye(d)-q)*g != r-mu*one or (pi*g)[0] != 0:
            raise ArithmeticError('native Poisson identity failed')
        widths = []
        for i in range(d):
            values = [h[i, j]+g[j] for j in range(d) if q[i, j]]
            widths.append(max(values)-min(values))
            if sum(q[i, j]*(h[i, j]-mu+g[j]-g[i]) for j in range(d)):
                raise ArithmeticError('native martingale difference is not centered')
        certificates[name] = {
            'stationary_mean': F(mu), 'boundary_span': F(max(g)-min(g)),
            'increment_range': F(max(widths)), 'poisson_solution': tuple(map(F, g))}
    return q, pi, rewards, certificates


def radius(local_ticks, failure_bits=12, *, observable='occupation'):
    """Rational confidence radius; per-observable failure at most 2**(-s)."""
    n, s = _integer(local_ticks, 1), _integer(failure_bits, 1)
    if observable not in NAMES:
        raise ValueError('a supported local observable is required')
    # ceil(sqrt(n*(s+1)/2)), with only integer operations.
    target = (n*(s+1)+1)//2
    t = isqrt(target)
    if t*t < target:
        t += 1
    c = source()[3][observable]
    return (c['boundary_span']+c['increment_range']*t)/n


def estimate(currents, failure_bits=12):
    """Consume only one local tape; never discard unresolved histories."""
    decoded = bits.decode(currents)
    n = _integer(len(currents), 1)
    s = _integer(failure_bits, 1)
    observed = {'occupation': None, 'adjacent_product': None,
                'movement': F(sum(abs(j) for j in currents), n)}
    if decoded['initial_bit'] is not None:
        bit = decoded['initial_bit']
        occupation = product = 0
        for current in currents:
            after = bit+current
            occupation += bit
            product += bit*after
            bit = after
        observed.update(occupation=F(occupation, n), adjacent_product=F(product, n))
    entries = {}
    for name in NAMES:
        value = observed[name]
        error = radius(n, s, observable=name)
        box = (F(0), F(1)) if value is None else (max(F(0), value-error), min(F(1), value+error))
        mean = source()[3][name]['stationary_mean']
        entries[name] = {'observed_mean': None if value is None else str(value),
                         'confidence_interval': list(map(str, box)),
                         'stationary_prediction': str(mean),
                         'prediction_in_interval': box[0] <= mean <= box[1]}
    return {'local_incident_ticks': n, 'resolved': decoded['initial_bit'] is not None,
            'observables': entries, 'simultaneous_failure_upper': str(min(F(1), F(3, 2**s))),
            'coverage': 'unconditional at a fixed local horizon; no independence of moments',
            'all_predictions_retained': all(row['prediction_in_interval'] for row in entries.values())}


def native_tape(local_ticks, *, seed=731, port=0):
    """Source execution for a reproducible diagnostic, not a sampling axiom."""
    _integer(local_ticks, 1)
    bits._port(port)
    kernel = bits.transfer.native_counts()
    rng = random.Random(seed)
    state = kernel.states[0]  # Deliberately nonstationary, with no disclosed bit.
    tape = []
    while len(tape) < local_ticks:
        seam = kernel.seams[rng.randrange(len(kernel.seams))]
        out, (_, current) = bits.transfer.native.native(state, seam, rng.randrange(2))
        if port in seam:
            tape.append(-current if seam[0] == port else current)
        state = out
    return tape


def run():
    _, _, _, c = source()
    # Port relabeling preserves the exact proof constants, not a hypothesis
    # of independent local observers.
    for port in range(12):
        for name in NAMES:
            for key in ('stationary_mean', 'boundary_span', 'increment_range'):
                if source(port)[3][name][key] != c[name][key]:
                    raise ArithmeticError('port-independent readout certificate failed')
    n = 32768
    baseline = estimate(native_tape(n))
    # Actual counterprotocol: choose an incident seam whose two endpoints
    # are empty and repeatedly repair only that seam. Every current is zero.
    kernel = bits.transfer.native_counts()
    edge = next(e for e in kernel.seams if 0 in e)
    state = next(s for s in kernel.states if not s[edge[0]] and not s[edge[1]])
    if any(bits.transfer.native.native(state, edge, coin)[0] != state for coin in (0, 1)):
        raise ArithmeticError('scheduler control is not a native fixed state')
    control = estimate([0]*n)
    if not baseline['all_predictions_retained'] or control['all_predictions_retained']:
        raise ArithmeticError('retained source/control diagnostic changed')
    paths = sorted(set(bits.transfer.SOURCE_PATHS+(
        'model/native_clock/bit_readout.py', 'model/native_clock/trajectory_readout.py',
        'model/observer/quantum_mixing.py')))
    return {'schema': 'oph.native_local_trajectory_readout.v1', 'all_pass': True,
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}),
            'certificates': {name: {k: list(map(str, v)) if isinstance(v, tuple) else str(v)
                                    for k, v in row.items()} for name, row in c.items()},
            'ports_checked': 12, 'baseline': baseline, 'fixed_empty_seam_control': control,
            'contract': {'laws': 'C and E complete local-current probabilities',
                         'initial_state': 'arbitrary binary total-two, including coherent/reference inputs',
                         'schedule': 'IID uniform seams; independent fair placements',
                         'readout': 'own signed incident tape, including zeros',
                         'horizon': 'fixed local incident count, before inspecting the tape',
                         'stationary_preparation_or_independent_trials_required': False,
                         'physical_observable_identified': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', nargs='?', const=str(OUTPUT))
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, sort_keys=True)+'\n'
    if args.check:
        if Path(args.check).read_text(encoding='utf-8') != payload:
            raise SystemExit('native trajectory readout receipt differs')
        print('native trajectory readout receipt verified')
    else:
        args.output.write_text(payload, encoding='utf-8', newline='\n')
        print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
