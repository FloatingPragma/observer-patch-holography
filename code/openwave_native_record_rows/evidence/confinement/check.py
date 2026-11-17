#!/usr/bin/env python3
"""Producer self-check using copied source and separate scalar flux enumeration."""
from collections import Counter
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import platform
import random
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
for line in (HERE / 'PREREG_MANIFEST.sha256').read_text().splitlines():
    digest, name = line.split('  ', 1)
    assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
spec = importlib.util.spec_from_file_location('custody_native_triality', HERE / 'inputs/native.py')
native = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = native
spec.loader.exec_module(native)

branches = json.loads((HERE / 'LOCAL_RECORDS.json').read_text())
local_calls = 0
for row in branches['native']['rows']:
    loads = tuple(int(i == row['initial_position']) for i in range(12))
    for coin, current in row['records']:
        loads, got = native.native(loads, (0, 1), coin)
        assert got == (coin, current)
        local_calls += 1
    assert list(loads) == row['final_loads']
    first, second = [abs(r[1]) for r in row['records']]
    delta = first * (1 - 2 * second)
    assert delta == row['delta']
    assert (row['initial_counter'] + delta) % 3 == row['final_counter']
for row in branches['changed_endpoint_hold']['rows']:
    assert all(current == 0 for coin, current in row['records'])
    assert row['initial_counter'] == row['final_counter']

audit = json.loads((HERE / 'COLOUR_PREPARATION.json').read_text())
native_tape_calls = 0
preparation_checks = []
for saved in audit['native'] + [audit['changed_endpoint_hold']]:
    bits, samples = saved['bits'], saved['samples']
    mutated = samples == 1
    suffix = '_hold' if mutated else ''
    tape = (HERE / f'preparation_B{bits}{suffix}.tape').read_bytes()
    assert len(tape) == samples * 18 * bits == saved['source_calls']
    rng = random.Random(202611150022 + bits)
    loads = (1,) + (0,) * 11
    movement = []
    for byte in tape:
        coin, current = byte % 2, byte // 2 - 1
        assert coin == rng.getrandbits(1)
        if mutated:
            assert current == 0
        else:
            loads, record = native.native(loads, (0, 1), coin)
            assert record == (coin, current)
            native_tape_calls += 1
        movement.append(abs(current))
    # Reconstruct every matrix from its complete raw bit window.
    matrices = []
    failed = []
    for sample in range(samples):
        start = sample * 18 * bits
        uniforms = []
        for k in range(18):
            value = sum(movement[start + k * bits + j] * 2 ** (bits - 1 - j)
                        for j in range(bits))
            uniforms.append((value + .5) / 2 ** bits)
        numbers = []
        for k in range(0, 18, 2):
            magnitude = math.sqrt(-2 * math.log(uniforms[k])) / math.sqrt(2)
            angle = 2 * math.pi * uniforms[k + 1]
            numbers.append(magnitude * complex(math.cos(angle), math.sin(angle)))
        q, r = np.linalg.qr(np.array(numbers).reshape(3, 3))
        diagonal = np.diag(r)
        if min(abs(diagonal)) <= 2 ** (-bits / 4):
            q = np.eye(3, dtype=complex)
            failed.append(sample)
        else:
            for k in range(3):
                q[:, k] *= diagonal[k] / abs(diagonal[k])
            q[:, 2] *= np.conj(np.linalg.det(q))
        matrices.append(q.reshape(9))
    values = np.array(matrices)
    mean = values.mean(axis=0)
    second = values.conj().T @ values / samples
    recorded_mean = np.array(saved['mean_real']) + 1j * np.array(saved['mean_imag'])
    recorded_second = np.array(saved['second_real']) + 1j * np.array(saved['second_imag'])
    assert failed == saved['failure_indices']
    assert np.max(abs(mean - recorded_mean)) < 2e-14
    assert np.max(abs(second - recorded_second)) < 2e-14
    preparation_checks.append(dict(bits=bits, mutated=mutated, failures=len(failed),
                                   records=len(tape), all_raw_letters_replayed=True))

flux_file = json.loads((HERE / 'ALL_FLUX_SECTORS.json').read_text())
flux_checks = []
for graph in flux_file:
    vertices, edges = graph['vertices'], graph['edges']
    table = Counter()
    for word in itertools.product(range(3), repeat=len(edges)):
        charge = [0] * vertices
        for (a, b), t in zip(edges, word):
            charge[a] += t
            charge[b] -= t
        code = sum((x % 3) * 3 ** i for i, x in enumerate(charge))
        table[code, sum(t != 0 for t in word)] += 1
    supplied = {(row['charge_base3'], k): count
                for row in graph['complete_sectors']
                for k, count in enumerate(row['support_multiplicities']) if count}
    assert dict(table) == supplied
    assert sum(table.values()) == graph['full_flux_states']
    for pair in graph['endpoint_comparisons']:
        a, b = pair['endpoints']
        code = 3 ** a + 2 * 3 ** b
        levels = [k for (q, k), count in table.items() if q == code]
        assert min(levels) == pair['distance'] == pair['minimum_support']
    flux_checks.append(dict(graph=graph['name'], flux_states=sum(table.values()),
                            full_multiplicity_table_match=True))

result = dict(all_pass=True, role='producer self-check; not independent G6 review',
              preregistration_leaves_unchanged=True,
              copied_native_replay_calls=local_calls + native_tape_calls,
              preparation_checks=preparation_checks, complete_flux_checks=flux_checks,
              runtime={'python':sys.version,'numpy':np.__version__,'platform':platform.platform()},
              scope='All raw source letters and all finite flux-sector multiplicities replayed. Matrix comparison is floating-point diagnostic, not certified exact-decoder forward error.')
with (HERE / 'SELF_CHECK.json').open('x') as f:
    json.dump(result, f, indent=2, sort_keys=True)
    f.write('\n')
print(json.dumps(result, sort_keys=True))
