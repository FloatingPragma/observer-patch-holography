#!/usr/bin/env python3
"""Independent route and full-superoperator check; no producer imports."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
import numpy as np


def main(packet):
    started = time.perf_counter()
    read = lambda name: json.loads((packet / name).read_text())
    route = read('ROUTE_RECORDS.json')
    calls = 0
    for trial in route['trials']:
        state = list(trial['input'])
        assert len(state) == 240 and sum(state) == 1
        for leg in trial['legs']:
            for a, b, coin, current in leg['records']:
                total = state[a] + state[b]
                # Ceiling/floor placement is independent of the producer's
                # native call and record-only subtraction implementation.
                high = -((-total) // 2)
                low = total // 2
                left, right = (high, low) if coin else (low, high)
                assert current == state[a] - left
                state[a], state[b] = left, right
                calls += 1
            assert state[leg['target']] == 1 and leg['success']
        assert state == trial['observed_final_source'] == [1] + [0] * 239
        assert trial['frame_native_event'] == [0, 1, 0, 1]
        h = 1 - Fraction(1, 2**32)
        assert Fraction(trial['success_probability']) == h ** len(trial['legs'])
        assert sum(map(Fraction, trial['first_failed_leg_probabilities'])) + h ** len(trial['legs']) == 1
    assert calls + 32 == route['fresh_native_calls'] == 293
    assert Fraction(route['route_failure_upper']) == 1 - h**5

    data = np.load(packet / 'inputs/first_cell_maps.npz')
    Q = data['certified_Q_proposal']
    D = data['dark_integer_basis']
    graph = read('inputs/graph.json')
    L = np.zeros((240, 240), dtype=int)
    for a, b in graph['seams']:
        L[a, a] += 1
        L[b, b] += 1
        L[a, b] -= 1
        L[b, a] -= 1
    assert np.array_equal(L @ D, 6 * D)
    projector = D[:12] @ np.linalg.inv(D.T @ D) @ D[:12].T
    seed = np.zeros((12, 12), complex)
    seed[0, 0] = 1
    assert np.linalg.norm(projector @ seed) == 0
    uniform = np.eye(12) / 12
    effect = np.eye(12) - Q.conj().T @ Q
    # Column-major vectorization yields a distinct 144-dimensional linear
    # representation of the complete escape/reinjection map.
    channel = np.kron(Q.conj(), Q) + np.outer(
        uniform.reshape(-1, order='F'), effect.T.reshape(-1, order='F'))
    final = (np.linalg.matrix_power(channel, 512) @ seed.reshape(-1, order='F')).reshape((12, 12), order='F')
    # Committing coordinate records gives a separate classical kernel.
    kernel = abs(Q)**2 + np.diag(effect).real[np.newaxis, :] / 12
    committed = np.linalg.matrix_power(kernel, 512)[:, 0]
    records = np.load(packet / 'FIELD_A_RECORDS.npz')
    assert np.max(abs(final - records['native_density_history'][-1])) < 2e-11
    assert np.max(abs(committed - np.diag(records['dephasing_density_history'][-1]).real)) < 2e-11
    result = read('FIELD_A_RESULT.json')
    ry = result['exact_rotations']['Ry']
    assert np.max(abs(projector[np.ix_(ry, ry)] - projector)) < 1e-14
    native_tv = np.sum(abs(np.diag(final).real - np.diag(final)[ry].real)) / 2
    control_tv = np.sum(abs(committed - committed[ry])) / 2
    assert abs(native_tv - result['native_coordinate_halfturn_TV']) < 2e-11
    assert abs(control_tv - result['dephasing_coordinate_halfturn_TV']) < 2e-11
    assert control_tv > .03 and native_tv < 1e-8
    transfer = read('TRANSFER_RESULT.json')
    native_upper = float(Fraction(transfer['full_position_halfturn_TV_upper']))
    control_lower = float(Fraction(transfer['active_dephasing_halfturn_TV_lower']))
    assert native_upper < .0005 and control_lower > .035
    print(json.dumps(dict(
        verified=True, producer_imports=[], continuing_routes=len(route['trials']),
        independently_reconstructed_route_events=calls,
        independent_source_superoperator_dimension=144,
        independent_committed_kernel_dimension=12,
        native_position_TV=native_tv, committed_position_TV=control_tv,
        certified_native_TV_upper=native_upper, certified_control_TV_lower=control_lower,
        seconds=time.perf_counter()-started), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('packet', type=Path)
    main(parser.parse_args().packet.resolve())
