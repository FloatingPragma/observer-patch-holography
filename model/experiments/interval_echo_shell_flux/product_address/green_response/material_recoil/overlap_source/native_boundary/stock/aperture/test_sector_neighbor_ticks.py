"""Independent elementary networks with exactly known travel and source costs."""

from fractions import Fraction

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_neighbor_ticks as ticks,
)


def test_triangle_slowed_edge():
    adj = np.ones((3, 3), int) - np.eye(3, dtype=int)
    uv = np.array([[0, 1], [0, 2], [1, 2]])
    u, v = uv.T
    _, G, blank, _ = ticks.clocks.kernels(adj, uv, np.array([2, 2, 2]))
    K, _, slow, _ = ticks.clocks.kernels(adj, uv, np.array([4, 2, 2]))
    # A resistance-two edge is in parallel with a two-edge path of cost two.
    expected = [1, float(Fraction(3, 4)), float(Fraction(3, 4))]
    np.testing.assert_allclose(slow[u, v], expected)
    np.testing.assert_allclose((slow - blank)[u, v].sum(), 0.5)
    np.testing.assert_allclose(
        ticks.predict_edges(ticks.transfer(G, u, v, [0]), [0], np.array([2, 1, 1])),
        (slow - blank)[u, v],
    )
    # Solve each directed passage as a killed chain without using the Green matrix.
    times = []
    for target, start in ((0, 1), (1, 0)):
        keep = np.array([i for i in range(3) if i != target])
        H = np.linalg.solve(np.eye(2) - K[np.ix_(keep, keep)] / 36, np.ones(2))
        times.append(H[np.flatnonzero(keep == start).item()])
    np.testing.assert_allclose(sum(times), 18)


def test_separate_articulation_blocks():
    uv = np.array([[0, 1], [1, 2]])
    u, v = uv.T
    adj = np.zeros((3, 3), int)
    adj[u, v] = adj[v, u] = 1
    _, G, _, _ = ticks.clocks.kernels(adj, uv, np.array([2, 2]))
    changed = ticks.predict_edges(ticks.transfer(G, u, v, [0]), [0], np.array([2, 1]))
    np.testing.assert_allclose(changed, [1, 0], atol=1e-13)
