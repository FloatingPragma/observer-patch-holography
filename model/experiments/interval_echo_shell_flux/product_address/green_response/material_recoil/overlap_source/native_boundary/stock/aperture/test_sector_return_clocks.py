"""Exact small-walk checks of the clock identification and observed pair estimator."""

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_return_clock_lab as lab,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_return_clocks as clocks,
)


def test_killed_return_mean_and_variance_match_green_diagonal():
    A = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
    edges = np.array([[0, 1], [1, 2]])
    for times in ([2, 2], [3, 4], [4, 2]):
        K, G, R, eps = clocks.kernels(A, edges, np.array(times))
        for i in range(len(K)):
            mean, var = clocks.killed_moments(K, i)
            assert abs(mean - 3) < 1e-12
            assert abs(var - (108 * G[i, i] - 6)) < 1e-11
        assert abs(R[0, 2] - sum(times) / 2) < eps


def test_native_pair_counts_exclude_pairing_a_member_with_itself():
    # Each register contains one of each word. Distinct members must
    # differ and two independent registers have chance one half to match.
    h = np.array([[1, 1], [1, 1], [1, 1]])
    d, p = lab.pair_law(3, h.sum(0), h.T @ h, 2)
    np.testing.assert_allclose(d, [[-0.25, 0.25], [0.25, -0.25]])
    np.testing.assert_allclose(p, [0.5, 0.5])
