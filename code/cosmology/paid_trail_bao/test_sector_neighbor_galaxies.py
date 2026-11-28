"""Analytic disk ring, finite support, and the shared-gain shape check."""

import numpy as np

from . import sector_neighbor_galaxies as sky


def test_ring_and_finite_support():
    r = np.array([0.0, 1.0])
    np.testing.assert_allclose(sky.within(1.0, r, 1.0), [0, 1 / 3])
    np.testing.assert_allclose(sky.within(0.0, r, 0.5), [1, 0])
    result = sky.convolve(
        np.array([0.0, 1.0, 2.0]), np.array([0.0]), 1.0, [3.0, 2.0], np.array([[2.0]])
    )
    np.testing.assert_allclose(result, [[6.0, 4.0, 0.0]])


def test_shared_positive_gain():
    v = np.array([1.0, -1.0])
    e = np.full(2, 0.1)
    assert sky.feasible(v, e, [0.95, -1.1], [1.05, -0.9])
    assert not sky.feasible(v, e, [1, 1], [2, 2])
