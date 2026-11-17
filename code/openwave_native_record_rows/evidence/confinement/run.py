#!/usr/bin/env python3
"""Preregistered native centre instrument and complete finite flux audit."""
from collections import Counter, deque
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import random
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
NATIVE_PATH = ROOT / 'reverse-engineering-reality/model/observer/native.py'
NATIVE_HASH = 'dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d'
assert hashlib.sha256(NATIVE_PATH.read_bytes()).hexdigest() == NATIVE_HASH
spec = importlib.util.spec_from_file_location('row22_triality_native', NATIVE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def dump(path, obj):
    with path.open('x') as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write('\n')


def primitive(loads, coin, mutated=False):
    if mutated:
        return loads, (coin, 0)
    return module.native(loads, (0, 1), coin)


def local_instrument(mutated=False):
    rows = []
    matrices = []
    for initial_position in (0, 1):
        matrix = [[Fraction(0) for _ in range(3)] for _ in range(3)]
        for initial_counter in range(3):
            for word in itertools.product((0, 1), repeat=2):
                loads = tuple(int(i == initial_position) for i in range(12))
                records = []
                for coin in word:
                    loads, record = primitive(loads, coin, mutated)
                    records.append(list(record))
                a, b = (abs(r[1]) for r in records)
                delta = 0 if a == 0 else 1 if b == 0 else -1
                final_counter = (initial_counter + delta) % 3
                matrix[final_counter][initial_counter] += Fraction(1, 4)
                rows.append(dict(initial_position=initial_position,
                                 initial_counter=initial_counter,
                                 records=records, final_loads=list(loads),
                                 delta=delta, final_counter=final_counter,
                                 probability='1/4'))
        matrices.append(matrix)
    expected = [[Fraction(int(i == j), 1) if mutated else
                 Fraction(1, 2) if i == j else Fraction(1, 4)
                 for j in range(3)] for i in range(3)]
    assert all(matrix == expected for matrix in matrices)
    return dict(rows=rows, matrix=[[str(x) for x in row] for row in expected],
                calls=48, eigenvalues=['1', '1', '1'] if mutated else ['1', '1/4', '1/4'])


def all_fluxes(name, vertices, edges):
    m = len(edges)
    count = 3 ** m
    indices = np.arange(count, dtype=np.int64)
    q = np.zeros((count, vertices), dtype=np.int8)
    support = np.zeros(count, dtype=np.int8)
    for e, (a, b) in enumerate(edges):
        t = ((indices // (3 ** e)) % 3).astype(np.int8)
        q[:, a] = (q[:, a] + t) % 3
        q[:, b] = (q[:, b] - t) % 3
        support += (t != 0)
    qcode = q @ np.array([3 ** i for i in range(vertices)], dtype=np.int64)
    histogram = np.bincount(qcode * (m + 1) + support,
                           minlength=3 ** vertices * (m + 1)).reshape(-1, m + 1)
    assert int(histogram.sum()) == count
    assert np.all(q.sum(axis=1) % 3 == 0)
    adjacency = [[] for _ in range(vertices)]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    comparisons = []
    for a in range(vertices):
        distances = {a: 0}
        queue = deque([a])
        while queue:
            v = queue.popleft()
            for w in adjacency[v]:
                if w not in distances:
                    distances[w] = distances[v] + 1
                    queue.append(w)
        for b in range(a + 1, vertices):
            population = histogram[3 ** a + 2 * 3 ** b]
            minimum = int(np.flatnonzero(population)[0])
            assert minimum == distances[b]
            comparisons.append(dict(endpoints=[a, b], distance=distances[b],
                                    minimum_support=minimum,
                                    count_energy=str(Fraction(3 * minimum, 4)),
                                    multiplicity_at_minimum=int(population[minimum]),
                                    total_charged_flux_states=int(population.sum())))
    complete_sectors = [{'charge_base3': i, 'support_multiplicities': row.tolist()}
                        for i, row in enumerate(histogram) if row.sum()]
    return dict(name=name, vertices=vertices, edges=edges, full_flux_states=count,
                complete_sectors=complete_sectors, endpoint_comparisons=comparisons,
                vacuum_minimum_support=int(np.flatnonzero(histogram[0])[0]))


class NativeBits:
    def __init__(self, seed, mutated=False):
        self.rng = random.Random(seed)
        self.loads = (1,) + (0,) * 11
        self.raw = bytearray()
        self.mutated = mutated

    def bit(self):
        coin = self.rng.getrandbits(1)
        self.loads, (coin, current) = primitive(self.loads, coin, self.mutated)
        assert current in (-1, 0, 1)
        self.raw.append(coin + 2 * (current + 1))
        return abs(current)

    def uniform(self, bits):
        k = 0
        for _ in range(bits):
            k = 2 * k + self.bit()
        return (k + 0.5) / 2 ** bits


def colour_draw(source, bits):
    reals = []
    for _ in range(9):
        u, v = source.uniform(bits), source.uniform(bits)
        radius = math.sqrt(-2 * math.log(u))
        reals += [radius * math.cos(2 * math.pi * v),
                  radius * math.sin(2 * math.pi * v)]
    z = (np.array(reals[::2]) + 1j * np.array(reals[1::2])).reshape(3, 3) / math.sqrt(2)
    q, r = np.linalg.qr(z)
    diagonal = np.diag(r)
    threshold = 2 ** (-bits / 4)
    failed = float(np.min(np.abs(diagonal))) <= threshold
    if failed:
        u = np.eye(3, dtype=complex)
    else:
        q = q @ np.diag(diagonal / np.abs(diagonal))
        q[:, -1] *= np.conjugate(np.linalg.det(q))
        u = q
    residual = max(float(np.linalg.norm(u.conj().T @ u - np.eye(3))),
                   float(abs(np.linalg.det(u) - 1)))
    return u, failed, residual


def preparation(bits, samples, mutated=False):
    source = NativeBits(202611150022 + bits, mutated)
    mean = np.zeros(9, dtype=complex)
    second = np.zeros((9, 9), dtype=complex)
    failures = []
    residual = 0.0
    for i in range(samples):
        u, failed, defect = colour_draw(source, bits)
        vector = u.reshape(9)
        mean += vector
        second += np.outer(vector.conj(), vector)
        residual = max(residual, defect)
        if failed:
            failures.append(i)
    mean /= samples
    second /= samples
    name = f'preparation_B{bits}' + ('_hold' if mutated else '')
    with (HERE / f'{name}.tape').open('xb') as f:
        f.write(source.raw)
    mean_error = float(np.max(np.abs(mean)))
    second_error = float(np.max(np.abs(second - np.eye(9) / 3)))
    result = dict(bits=bits, samples=samples, source_calls=len(source.raw),
                  raw_record_encoding='byte=coin+2*(current+1); initial loads=(1,0,...,0); seam=(0,1)',
                  initial_centre_preparations=[0, 1, 2],
                  centre_averaging='analytic equal mixture; no extra native calls claimed',
                  failure_indices=failures, failure_count=len(failures),
                  maximum_unitarity_determinant_residual=residual,
                  mean_real=mean.real.tolist(), mean_imag=mean.imag.tolist(),
                  second_real=second.real.tolist(), second_imag=second.imag.tolist(),
                  max_uncentred_mean_error=mean_error, max_second_moment_error=second_error,
                  diagnostic_gate=mean_error < .12 and second_error < .12,
                  precision_gate=residual < 2e-12)
    if mutated:
        assert len(failures) == samples and not result['diagnostic_gate']
    return result


def main():
    assert not (HERE / 'RECEIPT.json').exists(), 'Outputs are immutable; use a new candidate for reruns.'
    native = local_instrument()
    hold = local_instrument(True)
    dump(HERE / 'LOCAL_RECORDS.json', dict(native=native, changed_endpoint_hold=hold))
    graphs = [all_fluxes('chain_7', 7, [(i, i + 1) for i in range(6)]),
              all_fluxes('cycle_6', 6, [(i, (i + 1) % 6) for i in range(6)]),
              all_fluxes('ladder_2_by_3', 6,
                         [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]),
              all_fluxes('cube', 8, [(i, i ^ (1 << k)) for i in range(8)
                                    for k in range(3) if i < (i ^ (1 << k))])]
    dump(HERE / 'ALL_FLUX_SECTORS.json', graphs)
    preparations = [preparation(bits, 1024) for bits in (8, 12, 16)]
    control = preparation(16, 1, True)
    dump(HERE / 'COLOUR_PREPARATION.json', dict(native=preparations, changed_endpoint_hold=control))
    # Character expectation is obtained exactly from the measured local matrix.
    energy = [dict(distance=d, exact_one_layer_character=str(Fraction(1, 4 ** d)),
                   count_energy=str(Fraction(3 * d, 4)),
                   layer_energy=d * math.log(4), hold_character='1', hold_energy='0')
              for d in range(1, 7)]
    receipt = dict(candidate='A-native-triality-string-v1', row=22, class_='N+ with conditional colour interface',
                   native_sha256=NATIVE_HASH, native_calls=48 + sum(p['source_calls'] for p in preparations),
                   mutated_hold_calls=48 + control['source_calls'],
                   exact_local_histories_per_law=24,
                   full_flux_states=sum(g['full_flux_states'] for g in graphs),
                   endpoint_pairs=sum(len(g['endpoint_comparisons']) for g in graphs),
                   exact_flux_ground_energy_gate=True, exact_local_source_gate=True,
                   native_preparation_diagnostic_gates=[p['diagnostic_gate'] for p in preparations],
                   native_preparation_precision_gates=[p['precision_gate'] for p in preparations],
                   changed_law_control_rejected=True, all_failures_retained=True,
                   exact_character_slope=energy,
                   clocks={'count_per_address_coefficient':'3/4', 'whole_layer_coefficient':'log(4)',
                           'physical_clock':'unidentified'},
                   full_row_criterion='OPEN: physical endpoints, shared clock/ruler and ~1 GeV/fm tension',
                   no_source_execution_for_flux_enumeration=True)
    dump(HERE / 'RECEIPT.json', receipt)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
