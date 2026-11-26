import json
from pathlib import Path

import numpy as np
import pytest

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation
from . import field, scatter
from .readback import intervals


def receipt():
    return json.loads(Path(__file__).with_name('results.json').read_text())


def test_passive_source_and_complete_written_blocks():
    h0, h1 = law.local_term_parts()
    h = h0 + 8 * h1
    np.testing.assert_allclose(h[9:, :9], 0, atol=2e-14)
    for word in [(0, 2), (2, 0), (1, 1)]:
        source = law.echo()[:, 3*word[0] + word[1]]
        np.testing.assert_allclose(h @ source, 0, atol=2e-14)
        np.testing.assert_allclose((abs(source.reshape(-1, 9))**2).sum(axis=0), np.eye(9)[4], atol=2e-15)
    for sign in (1, -1):
        blank, left, right = scatter.central_blocks(h, sign, [-2, 2])
        np.testing.assert_allclose(blank, 0, atol=2e-14)
        np.testing.assert_allclose(left, [[32.5, -.5], [-.5, 32.5]], atol=2e-13)
        np.testing.assert_allclose(right, left)
        np.testing.assert_allclose(scatter.central_blocks(h, sign, [0]), 0, atol=2e-14)


def test_direct_native_work_on_all_cuts_and_reference_routes():
    board = federation(0)
    term = law.seam_term(8)
    ends = np.array(board.seams)
    def cost(word):
        ids = 3 * word[ends[:, 0]] + word[ends[:, 1]]
        return np.sum(term[ids, ids]) + 8 * np.dot(6 - board.degree, word**2)
    for seam in (board.external[0], board.external[-1]):
        cuts, _ = field.mount_geometry(board, seam)
        for s in (1, -1):
            stock = np.ones(board.n, int)
            stock[seam[0]] += s
            host = np.ones(board.n, int)
            for radius, cut in cuts.items():
                total = 0.
                for a, b in cut:
                    path = field.private_route(board, b, {a, *seam})
                    assert board.degree[path[-1]] == 5 and not set(path) & {a, *seam}
                    assert all(y in board.neighbours[x] for x, y in zip(path, path[1:]))
                    for p in (1, -1):
                        u = np.eye(1, board.n, a, dtype=int)[0] * p
                        v = np.eye(1, board.n, b, dtype=int)[0] * p
                        total += p * (cost(stock + v) - cost(stock + u) - cost(host + v) + cost(host + u)) / 2
                assert total == (162.5 * s if radius == 1 else 0)


def test_counts_and_field_identification():
    data = receipt()
    computed = intervals(data['boundary_work'], data['mobile_probe'])
    expected = data['count_only_readback'].copy()
    # libm phase reductions and inverse beta quantiles differ in their last bits.
    # Preserve exact equality for the counts, selected appointment and coverage.
    tolerance = dict(rel=2e-14, abs=0)
    for name in ('maximum_abs_charged_outer_ratio', 'maximum_abs_written_extreme_outer_ratio',
                 'simultaneous_upper_on_abs_charged_outer_ratio'):
        expected[name] = pytest.approx(expected[name], **tolerance)
    bands = 'simultaneous_selected_endpoint_intervals'
    expected[bands] = {source: pytest.approx(bounds, **tolerance)
                       for source, bounds in expected[bands].items()}
    assert computed == expected
    assert computed['simultaneous_upper_on_abs_charged_outer_ratio'] < .1
    a = computed['simultaneous_selected_endpoint_intervals']
    assert a['write02'][1] < a['write11'][0]
    for bank in data['mobile_probe']['banks']:
        for source in bank['receiver_records'].values():
            assert sum(source['owner_counts']) == 2**20
    for diag in data['mobile_probe']['diagnostics']:
        assert max(diag[key] for key in ('exact_barrier_twin_max_port_error', 'cut_max_port_error',
                                         'reversal_max_port_error')) < 5e-11
