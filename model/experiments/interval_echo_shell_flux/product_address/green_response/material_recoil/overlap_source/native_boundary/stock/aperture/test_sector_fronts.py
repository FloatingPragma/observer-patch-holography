"""Independent checks of geodesic counting and full-register covariance."""

from itertools import combinations

import numpy as np
from scipy.sparse.csgraph import shortest_path

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    peer_population,
    sector_fronts,
    sector_fronts_counts,
    sector_fronts_lab,
)


def test_prompt_counts_and_ports_are_separate_observations():
    adj = np.array([[0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 0, 1], [0, 1, 1, 0]])
    d = shortest_path(adj).astype(int)
    routes = list(sector_fronts_lab.enumerate_routes(adj, d, 0, 3))
    assert routes == [(0, 1, 3), (0, 2, 3)]
    occupied = np.array([False, True, False, False])
    n, port = sector_fronts.prompt(adj, d, occupied)
    u, v = np.nonzero(adj)
    assert n[0, 3] == 1
    assert u[(v == 3) & (port[0] > 0)].tolist() == [2]


def test_collapsed_register_second_moments_against_all_six_allocations():
    l, r = np.array(peer_population.fibers.BLANK_PAIRS).T
    slots = np.array([0, 0, 1, 2])
    hist = []
    for selection in combinations(range(4), 2):
        v = np.zeros(16, dtype=int)
        v[r] = np.bincount(slots[list(selection)], minlength=3)
        v[l] = [2, 1, 1] - v[r]
        hist.append(v)
    mean, cov = sector_fronts_counts.moments([(6, hist[0])])
    np.testing.assert_allclose(np.array(mean, float), np.mean(hist, 0))
    np.testing.assert_allclose(
        np.array(cov, float), np.cov(np.array(hist).T, bias=True), atol=1e-15
    )
