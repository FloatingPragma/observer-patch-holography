"""Recover a binary local field from its own signed current records.

Before its first nonzero current the local bit is unknown. A +1 identifies
the outgoing bit as one, and a -1 identifies it as zero; subsequent currents
update that bit. The earlier bit can then be reconstructed retrospectively.
For coherent inputs this is an instrument-domain statement, not a hidden
classical value assigned before the read. Every native C/E Kraus branch is a
partial injection, remote/zero-current branches preserve the local-number
projectors, and a signed branch has its specified input/output support. These
identities survive arbitrary untouched references and branch composition.

Uniform independent seams and fair placement supply the waiting bound. A
checked twelve-class strong lumping of the actual 66-state source preserves
all local signed-current words. Summing remote waits gives Q_j=(I-A)^-1 B_j.
Q_0^16 1 <= (3/8)1, so unresolved mass after m LOCAL INCIDENT TICKS is at most
(3/8)^floor(m/16), uniformly over initial states. All failure mass is retained.
There is no stationary preparation, initial-bit disclosure or hidden attempt
counter in the decoder. Coherent interventions and record-adaptive schedules
are outside this passive observation contract. For a window starting at a
declared checkpoint, let r_b be the unconditional mass of resolved initial
label b and u the unresolved mass. Then r_1<=Tr(rho N_port)<=r_1+u. This
interval retains the unknown outcomes; its endpoints are ensemble laws,
not confidence limits inferred from a finite sample. The complete native
instrument supports the inequality also for coherent, reference-entangled
inputs. Retrodiction supplies no control access before resolution.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as sp

from ..primitive_observer import transfer
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('bit_readout_results.json')
MARKS = (None, -1, 0, 1)  # remote, then signed local changes


def decode(currents):
    """Only signed incident records enter; an all-zero prefix stays unknown."""
    if (type(currents) not in (tuple, list)
            or any(type(j) is not int or j not in (-1, 0, 1) for j in currents)):
        raise ValueError('a finite list of signed binary incident currents is required')
    initial, bit, resolved_at = None, None, None
    for tick, current in enumerate(currents, 1):
        if bit is None and current:
            bit = int(current == 1)
            initial = 1-bit
            resolved_at = tick
        elif bit is not None:
            bit += current
            if bit not in (0, 1):
                raise ValueError('current records violate binary local continuity')
    return {'initial_bit': initial, 'current_bit': bit,
            'resolved_at_local_tick': resolved_at, 'local_ticks': len(currents)}


def _port(port):
    if type(port) is not int or port not in range(12):
        raise ValueError('a native port from zero through eleven is required')


def occupation_interval(resolved_zero, resolved_one, unresolved):
    """Field-mean interval from a complete resolved/unknown ensemble law."""
    values = (resolved_zero, resolved_one, unresolved)
    if (any(type(x) not in (int, F) or x < 0 for x in values)
            or sum(values) != 1):
        raise ValueError('three exact nonnegative masses summing to one are required')
    return F(resolved_one), F(resolved_one)+F(unresolved)


@lru_cache(maxsize=12)
def source(port=0):
    """Exact marked strong lumping; classes are proof coordinates, not reads."""
    _port(port)
    kernel = transfer.native_counts()
    index = {state: i for i, state in enumerate(kernel.states)}
    raw = {mark: [[0]*66 for _ in range(66)] for mark in MARKS}
    for i, state in enumerate(kernel.states):
        for seam in kernel.seams:
            for coin in (0, 1):
                out, _ = transfer.native.native(state, seam, coin)
                change = out[port]-state[port]
                if port not in seam and change:
                    raise ArithmeticError('a remote branch changes the local field')
                mark = change if port in seam else None
                raw[mark][i][index[out]] += 1
    blocks = [tuple(i for i, s in enumerate(kernel.states) if s[port] == bit)
              for bit in (0, 1)]
    while True:
        refined = []
        for block in blocks:
            groups = {}
            for i in block:
                key = tuple(sum(raw[m][i][j] for j in dest)
                            for m in MARKS for dest in blocks)
                groups.setdefault(key, []).append(i)
            refined.extend(tuple(group) for group in groups.values())
        if refined == blocks:
            break
        blocks = refined
    matrices = {m: sp.ImmutableMatrix([
        [sp.Rational(sum(raw[m][block[0]][j] for j in dest), 60)
         for dest in blocks] for block in blocks]) for m in MARKS}
    size = len(blocks)
    one = sp.ones(size, 1)
    if size != 12 or matrices[None]*one != sp.Rational(5, 6)*one:
        raise ArithmeticError('native local-event partition changed')
    inverse = (sp.eye(size)-matrices[None]).inv(method='DM')
    stopped = {m: sp.ImmutableMatrix(inverse*matrices[m]) for m in (-1, 0, 1)}
    if sum(stopped.values(), sp.zeros(size))*one != one:
        raise ArithmeticError('summing remote waits lost probability')
    if any(x < 0 for matrix in stopped.values() for x in matrix):
        raise ArithmeticError('stopped current maps are not positive')
    return kernel, tuple(blocks), matrices, stopped


def quantum_support_audit(port=0):
    """Exact source Kraus-domain identities for both C and E, all columns."""
    _port(port)
    kernel = transfer.native_counts()
    index = {state: i for i, state in enumerate(kernel.states)}
    checked = {}
    for law in ('C', 'E'):
        groups, marks = defaultdict(dict), {}
        for edge, seam in enumerate(kernel.seams):
            for coin in (0, 1):
                for i, state in enumerate(kernel.states):
                    out, (_, current) = transfer.native.native(state, seam, coin)
                    label = current if law == 'C' else state[seam[0]]-state[seam[1]]
                    key = edge, coin, label
                    mark = out[port]-state[port] if port in seam else None
                    if key in marks and marks[key] != mark:
                        raise ArithmeticError('one Kraus branch has conflicting visible currents')
                    marks[key] = mark
                    groups[key][i] = index[out]
        coverage = [0]*66
        for key, mapping in groups.items():
            if len(set(mapping.values())) != len(mapping):
                raise ArithmeticError('native branch is not a partial injection')
            mark = marks[key]
            for i, j in mapping.items():
                before, after = kernel.states[i][port], kernel.states[j][port]
                if mark in (None, 0):
                    valid = before == after
                else:
                    valid = (before, after) == ((0, 1) if mark == 1 else (1, 0))
                if not valid:
                    raise ArithmeticError('local number-projector intertwining failed')
                coverage[i] += 1
        if set(coverage) != {60}:
            raise ArithmeticError('complete native Kraus coverage failed')
        checked[law] = {'branches': len(groups), 'source_columns': sum(coverage)}
    return checked


def timeout_bound(local_ticks):
    if type(local_ticks) is not int or local_ticks < 0:
        raise ValueError('a nonnegative local incident count is required')
    return F(3, 8)**(local_ticks//16)


def certificate(port=0):
    kernel, blocks, _, stopped = source(port)
    one = sp.ones(len(blocks), 1)
    survival = one
    for _ in range(16):
        survival = stopped[0]*survival
    if max(survival) > F(3, 8):
        raise ArithmeticError('the uniform sixteen-tick resolution bound failed')
    expected = (sp.eye(len(blocks))-stopped[0]).inv(method='DM')*one
    arrival = (stopped[-1]+stopped[1])*one
    class_bits = [kernel.states[block[0]][port] for block in blocks]
    for mark, matrix in stopped.items():
        for i in range(len(blocks)):
            for j in range(len(blocks)):
                if matrix[i, j] and class_bits[j]-class_bits[i] != mark:
                    raise ArithmeticError('stopped readout lost number-projector support')
    return {'port': port, 'class_sizes': list(map(len, blocks)),
            'first_incident_success_min': str(min(arrival)),
            'first_incident_success_max': str(max(arrival)),
            'sixteen_tick_timeout_max': str(max(survival)),
            'expected_local_ticks_min': str(min(expected)),
            'expected_local_ticks_max': str(max(expected)),
            'timeout_bound_128_local_ticks': str(timeout_bound(128)),
            'twelve_port_union_bound_160_ticks_each': str(12*timeout_bound(160)),
            'quantum_support': quantum_support_audit(port)}


def run():
    paths = sorted(set(transfer.SOURCE_PATHS+('model/native_clock/bit_readout.py',
                                           'model/observer/quantum_mixing.py')))
    certs = [certificate(port) for port in range(12)]
    for field in ('first_incident_success_min', 'first_incident_success_max',
                  'sixteen_tick_timeout_max', 'expected_local_ticks_min',
                  'expected_local_ticks_max'):
        if len({row[field] for row in certs}) != 1:
            raise ArithmeticError('native port symmetry failed for the readout law')
    return {'schema': 'oph.native_local_bit_readout.v1', 'all_pass': True,
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                              for p in paths}),
            'reference_port': certs[0], 'ports_checked': len(certs),
            'waiting_law_port_invariance_verified': True,
            'contract': {
                'input': 'arbitrary binary total-two quantum state with untouched references',
                'law': 'C or E; uniform independent seams and fair placements',
                'access': 'only own signed incident currents, including zero',
                'initial_bit_disclosure_required': False,
                'timeout': 'unresolved bit remains unknown and retains all probability mass',
                'clock': 'local incident ticks; no global attempt counter or seconds',
                'initial_read': 'retrospective on resolved histories, not earlier causal access',
                'field_mean_interval': 'resolved_one_mass <= mean_initial_local_bit <= resolved_one_mass+unresolved_mass',
                'sampling': 'ensemble-law interval; finite-sample confidence requires a separate sampling contract',
                'postselection': 'retaining resolution can bias later field statistics; use complete law',
                'all_port_bound': 'union bound, no independence; joint access is an additional service',
                'physical_clock_or_field_identified': False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', nargs='?', const=str(OUTPUT))
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, sort_keys=True)+'\n'
    if args.check:
        if Path(args.check).read_text(encoding='utf-8') != payload:
            raise SystemExit('native local-bit readout receipt differs')
        print('native local-bit readout receipt verified')
    else:
        args.output.write_text(payload, encoding='utf-8', newline='\n')
        print(f'wrote {args.output}')


if __name__ == '__main__':
    main()
