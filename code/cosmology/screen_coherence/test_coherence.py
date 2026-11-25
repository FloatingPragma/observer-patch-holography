from itertools import product
from unittest.mock import patch

import numpy as np

from measure import apply_mean, near_half, paired_moves, federation, observations
from survey import counts
from model.candidate_laws.interval_echo.openwave.row_28_acceleration import (
    experiment as repair,
)


def test_routed_ports_match_seeded_half_board_intervention():
    tables = observations.instruments()
    for level in (0, 1, 2):
        board = federation(level)
        groups = [(g[0], g[1]) for g in observations.schedules(board)[0]]
        near, _ = near_half(board.cells, groups)
        initial = np.zeros(board.cells, dtype=np.uint8)
        initial[near] = 1
        seeds = np.random.SeedSequence([6182, level]).spawn(3)
        prep, coin, move = map(np.random.default_rng, seeds)

        class SharpPreparation:
            def random(self, size):
                if size == board.cells:
                    return 1 - initial
                return prep.random(size)

        with patch.object(
            np.random, "default_rng", side_effect=[SharpPreparation(), coin, move]
        ):
            _, raw = observations.measure(level, 24, 6182, tables)
        expected = raw["singleton_history"][:, :, 1]
        state = initial.copy()
        repeat = np.random.default_rng(seeds[2])
        np.testing.assert_array_equal(state, expected[0])
        for t in range(24):
            x, y = groups[t % len(groups)]
            paired_moves(state, x, y, repeat.random(len(x)) < tables[1][0, 2])
            np.testing.assert_array_equal(state, expected[t + 1])


def test_unequal_time_covariance_by_enumerating_all_histories():
    p = 0.16
    mean_kernel = np.eye(3)
    edges = [(np.array([0]), np.array([1])), (np.array([1]), np.array([2]))]
    for x, y in edges:
        apply_mean(mean_kernel, x, y, 0.5)
    same, lagged = np.zeros((3, 3)), np.zeros((3, 3))
    for bits in product((0, 1), repeat=3):
        initial = np.array(bits)
        prob = p ** sum(bits) * (1 - p) ** (3 - sum(bits))
        for flags in product((False, True), repeat=2):
            state = initial.copy()
            for (x, y), take in zip(edges, flags):
                if take:
                    state[x], state[y] = state[y].copy(), state[x].copy()
            same += prob / 4 * np.outer(state - p, state - p)
            lagged += prob / 4 * np.outer(state - p, initial - p)
    np.testing.assert_allclose(same, np.eye(3) * p * (1 - p), atol=1e-14)
    np.testing.assert_allclose(lagged, p * (1 - p) * mean_kernel, atol=1e-14)


def test_reflector_counts_agree_with_independent_roundtrip_distances():
    n, edges, ranks, *_ = repair.layout(0)
    rng = np.random.default_rng(101)
    commit, _, _ = repair.history(
        n, edges, ranks, 2 * rng.integers(2, size=n, dtype=np.int8), 0.25, rng
    )
    for open_edges in (edges, edges[commit == 0]):
        d = repair.geometry(n, open_edges, np.arange(0, n, 12))
        np.testing.assert_array_equal(
            counts(n, open_edges),
            [((d >= 0) & (2 * d <= ticks)).sum(axis=1) for ticks in (4, 6)],
        )


def test_bounded_hand_stationarity_at_every_contact():
    transition = np.zeros((10, 10))
    for source, hand in product(range(2), range(5)):
        for flip in (0, 1):
            new = source ^ flip
            transition[source * 5 + hand, new * 5 + (hand + new) % 5] += 0.5
    np.testing.assert_array_equal(transition.sum(axis=0), np.ones(10))
    np.testing.assert_array_equal(transition.sum(axis=1), np.ones(10))
