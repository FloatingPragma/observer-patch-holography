"""Small exact reference checks for the joint register instruments."""

from itertools import product

import numpy as np
from scipy.sparse import csr_matrix

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_register_balls as balls,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_register_escape as escape,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    sector_register_price as price,
)


def test_exit_count_in_a_one_marker_domain():
    W = csr_matrix([[30, 6], [6, 30]])
    actual = escape.first_exits(W, np.array([True, False]))
    expected = [1 - (5 / 6) ** n for n in escape.DEADLINES]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-14)


def test_mixed_price_matches_enumerated_joint_distribution():
    cells = {(0, 0): 0.1, (0, 1): 0.2, (1, 0): 0.3, (1, 1): 0.4}
    n = 4
    prob = np.zeros((n + 1, n + 1))
    for terms in product(cells, repeat=n):
        x, y = map(sum, zip(*terms, strict=True))
        prob[x, y] += np.prod([cells[t] for t in terms])
    K = -np.log(prob)
    e = np.log(cells[0, 0] * cells[1, 1] / (cells[1, 0] * cells[0, 1]))
    for s, t in [(0, 0), (1, 2), (3, 3)]:
        exact = K[s, t] + K[s + 1, t + 1] - K[s + 1, t] - K[s, t + 1]
        val, (lo, hi) = price.mixed(n, s, t, [e, e])
        assert lo <= exact <= hi and abs(val - exact) < 1e-12


def test_pooling_changes_price_even_with_identical_marginals():
    p = np.array([0.5, 0.5])
    bands = np.column_stack([p, p])
    delta = np.array([[0.25, -0.25], [-0.25, 0.25]])
    for labels, expected in [([0, 1, 2, 3], 0.0), ([0, 1, 1, 2], np.log(2) / 2)]:
        value, bounds = balls.tangent(*balls.prices(labels, p, bands), delta, 0)
        assert abs(value - expected) < 1e-13 and bounds[0] <= expected <= bounds[1]


def test_connected_score_agrees_with_distinct_register_repairing():
    xy = np.array([[0, 1], [1, 2], [1, 2], [1, 0], [0, 0]])
    K = np.array([[0, 1, 3], [2, -1, 4]], float)
    hist = np.zeros(K.shape, int)
    np.add.at(hist, tuple(xy.T), 1)
    raw = np.mean([K[x, y] for x, y in xy]) - np.mean(
        [K[a[0], b[1]] for i, a in enumerate(xy) for j, b in enumerate(xy) if i != j]
    )
    assert abs(raw - price.connected(K, hist)["connected_modular_nats"]) < 1e-13
