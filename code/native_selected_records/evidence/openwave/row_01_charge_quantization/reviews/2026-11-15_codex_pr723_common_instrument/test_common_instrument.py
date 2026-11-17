"""Independent same-group test of the shared positive source interface.

The source under review is imported only after its content hash is checked.
Expected endpoints use separate floor/ceiling arithmetic. The inverse
Householder check retains all failure source kets, including higher currents.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import random

HERE = Path(__file__).resolve().parent
OW = HERE.parents[2]
SOURCE = HERE.parent / '2026-11-15_codex_pr723_integration'


def main():
    bindings = json.loads((HERE / 'SOURCE_BINDINGS.json').read_text())
    for relative, digest in bindings['files'].items():
        assert hashlib.sha256((OW / relative).read_bytes()).hexdigest() == digest, relative
    path = SOURCE / 'common_source_v2.py'
    spec = importlib.util.spec_from_file_location('reviewed_common_v2', path)
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    original = json.loads((OW / 'row_01_charge_quantization/work/'
        '2026-11-15_codex_pr723_unit_selection/seams.json').read_text())
    assert source.ALLOWED == frozenset(map(tuple, original)) | {(12, 13)}

    rng = random.Random(72314001)
    states = [(a, b) + (0,) * 12 for a in range(-8, 9) for b in range(-8, 9)]
    amplitudes = {q: complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for q in states}
    norm = math.sqrt(sum(abs(z) ** 2 for z in amplitudes.values()))
    amplitudes = {q: z / norm for q, z in amplitudes.items()}
    results = []
    for tickets in [[(0, 1)], [(0, 1), (1, 0), (12, 13)]]:
        observed = source.observe(amplitudes, tickets, 'optical_erase')
        fine = observed['fine']
        expected = {}
        gain = 1 / math.sqrt(2 * len(tickets))
        for q, value in amplitudes.items():
            for ticket, (a, b) in enumerate(tickets):
                for coin in (0, 1):
                    total = q[a] + q[b]
                    low, high = total // 2, total - total // 2
                    out = list(q)
                    out[a], out[b] = (low, high) if coin == 0 else (high, low)
                    current = q[a] - out[a]
                    branch = expected.setdefault((ticket, coin, current), {})
                    out = tuple(out)
                    branch[out] = branch.get(out, 0) + gain * value
        assert fine.keys() == expected.keys()
        assert all(fine[r].keys() == branch.keys() for r, branch in expected.items())
        endpoint_error = max(abs(z - fine[r][q]) for r, br in expected.items()
                             for q, z in br.items())
        assert endpoint_error < 1e-14
        fine_norm = sum(abs(z) ** 2 for v in fine.values() for z in v.values())
        full_norm = observed['success_probability'] + observed['all_failure_probability']
        returned_norm = sum(abs(z) ** 2 for z in observed['success'].values()) + sum(
            abs(z) ** 2 for v in observed['failure_states'].values() for z in v.values())
        assert abs(returned_norm - full_norm) < 1e-14
        assert abs(fine_norm - 1) < 1e-12 and abs(full_norm - 1) < 1e-12

        support = [(t, c, j) for t in range(len(tickets))
                   for c in (0, 1) for j in (-1, 0, 1)]
        w = 1 / math.sqrt(len(support))
        u = {r: (w - int(i == 0)) / math.sqrt(2 * (1 - w))
             for i, r in enumerate(support)}
        outcomes = dict(observed['failure_states'])
        outcomes[support[0]] = observed['success']
        contraction = {}
        for flag, coefficient in u.items():
            for q, value in outcomes[flag].items():
                contraction[q] = contraction.get(q, 0) + coefficient * value
        inverse_error = 0
        for flag in set(outcomes) | set(fine):
            target, actual = fine.get(flag, {}), outcomes.get(flag, {})
            for q in set(target) | set(actual) | set(contraction):
                inverse_error = max(inverse_error, abs(actual.get(q, 0)
                    - 2 * u.get(flag, 0) * contraction.get(q, 0) - target.get(q, 0)))
        assert inverse_error < 1e-13
        results.append({'tickets': len(tickets), 'input_coherent_terms': len(amplitudes),
            'independent_endpoint_error': endpoint_error, 'fine_norm': fine_norm,
            'complete_outcome_norm': full_norm, 'householder_inverse_error': inverse_error,
            'failure_flags_retained': len(observed['failure_states'])})
    rejected = []
    for edge in ((0, 2), (0, 12)):
        try:
            source.event((0,) * 14, edge, 0)
        except ValueError:
            rejected.append(edge)
        else:
            raise AssertionError('Invalid nonedge accepted')
    return {'candidate': 'A-pr723-common-positive-instrument-v2',
        'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'same_group': True, 'checks': results, 'allowed_graph_matches': True,
        'rejected_nonedges': rejected, 'passed': True,
        'complete_common_controller_established': False,
        'native_coherent_apparatus_established': False, 'physical_row_awards': False}


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
