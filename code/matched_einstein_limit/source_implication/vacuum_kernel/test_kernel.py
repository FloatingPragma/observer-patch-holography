from math import pi

import numpy as np

from .compiled import read
from .measure import SPEED, averaged_read, ball_sites, covariances, modular_read


def test_native_words_preserve_probability_and_change_local_rays():
    panels = read()
    np.testing.assert_allclose(panels["nominal"], [1/8, 0], atol=1e-13, rtol=0)
    np.testing.assert_allclose(panels["shift_only"], [1/8, 1], atol=1e-13, rtol=0)


def test_zero_mode_and_ball_measurement_by_independent_real_space_sum():
    sites, h = 24, 0.25
    points = ball_sites(1/h)
    xmat, pmat = covariances(sites, h)
    np.testing.assert_allclose([xmat.sum(), pmat.sum()], 0, atol=1e-13)
    coord = tuple((points[:, None] - points).transpose(2, 0, 1))
    variance = xmat[coord].mean()/h**3
    var_k, relative_entropy, _ = averaged_read(sites, h)
    np.testing.assert_allclose(variance, var_k, atol=1e-13, rtol=0)
    np.testing.assert_allclose(relative_entropy, 1/(2*variance), atol=1e-12, rtol=0)


def test_regional_modular_read_dominates_single_quadrature():
    for m in (1, 2, 3):
        h = 1/m
        full = modular_read(48, h)
        _, coarse, _ = averaged_read(48, h)
        assert full["field"] >= coarse
        assert full["symplectic_gap"] > 1e-12
        assert full["momentum"] > 0
    # The finite ratios compare independently with the analytic normalization.
    assert 1.08 < full["field_ratio"] < 1.10
    assert 1.12 < full["momentum_ratio"] < 1.15


def test_quadrature_variance_tracks_continuum_value():
    variance, _, _ = averaged_read(96, 0.125)
    expected = 9/(16*pi**2*SPEED)
    assert abs(variance / expected - 1) < 0.013
