"""Cross-check the archive reader and the uneven-region prediction on the native service."""

import numpy as np

from model.executable_model import federation
from model.candidate_laws.interval_echo.openwave.row_28_acceleration.experiment import (
    history,
)
from model.experiments.interval_echo_consensus.asynchronous import priority
from settled_load import measure_shard


def test_four_observer_readout():
    # Each pair shares a different gateway. Only one port starts double.
    neighbors = np.array(
        [[j for j in range(4) if j != i] for i in range(4)], dtype=np.uint32
    )
    raw = np.zeros((4, 48), np.uint8)
    raw[:, 26:32] = 255
    raw[:, 44:47] = [0, 1, 2]
    raw[0, 0] = 1
    for i in (0, 1):
        raw[i, 2] = 3
        raw[i, 26] = 0xFC
    values, pairs, regions = measure_shard(raw, 0, neighbors)
    initial = np.array([2, 0, 0, 0])
    settled = np.array([1, 1, 0, 0])
    expected = np.array([initial - 12, settled - 12, settled - initial]).T
    np.testing.assert_array_equal(values, expected)
    for field in range(3):
        v = expected[:, field]
        products = sum(v[i] * v[j] for i in range(4) for j in range(i))
        assert pairs[0, field].tolist() == [6, 3 * v.sum(), 3 * (v * v).sum(), products]
    assert regions[0, 6] == 12
    assert not pairs[1].any()


def test_polarized_region_in_native_history():
    for level in (0, 1, 2):
        board = federation(level)
        n_region = 4**level
        initial = np.zeros(board.n, dtype=np.int8)
        initial[: 12 * n_region] = 2
        edges = np.array(board.seams)
        ranks = np.array([priority(e) for e in board.seams])
        ticks, _, _ = history(
            board.n, edges, ranks, initial, 0.25, np.random.default_rng(20261126)
        )
        inside = edges < 12 * n_region
        crossing = inside[:, 0] != inside[:, 1]
        np.testing.assert_array_equal(ticks > 0, crossing)
        assert crossing.sum() == 3 * 2**level
        final = initial.copy()
        final[edges[ticks > 0].ravel()] = 1
        assert final[: 12 * n_region].sum() == 24 * n_region - crossing.sum()
