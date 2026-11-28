import numpy as np

from .renewal_disks import integral_at


def test_log_radius_integral_with_sign_change_and_intermediate_radii():
    radius = np.exp([0.0, 0.4, 1.2, 2.0])
    v = 2 * np.log(radius) - 1
    log_x = np.array([0.0, 0.2, 0.4, 0.9, 1.7, 2.0])
    np.testing.assert_allclose(
        integral_at(radius, v, np.exp(log_x)), 2 - log_x**2 + log_x, atol=1e-14
    )
