import numpy as np
from scipy.integrate import dblquad

from .dispersive_disks import LENGTH, disk_derivative


def test_annular_force_matches_independent_potential_difference():
    rings = np.array([3.3, 7.4, 12.5])
    profile = np.array([5.0, 3.0, 2.0, 0.0])

    # Integrate the actual potential of the middle annulus. This uses
    # ordinary polar area integration, with no derivative of the kernel.
    def potential(radius):
        def integrand(angle, rho):
            dist = np.sqrt(
                (radius - rho) ** 2 + 4 * radius * rho * np.sin(angle / 2) ** 2
            )
            return -rho * np.interp(dist / LENGTH, np.arange(4), profile)

        value, _ = dblquad(
            integrand, rings[0], rings[1], 0, np.pi, epsabs=4e-8, epsrel=2e-9
        )
        return 2 * value / (np.pi * (rings[1] ** 2 - rings[0] ** 2))

    r = rings[-1]
    width = 0.002
    independent = r * (potential(r + width) - potential(r - width)) / (2 * width)
    reported = disk_derivative(rings, 3, 32)[-1, 1] @ (-np.diff(profile))
    np.testing.assert_allclose(reported, independent, atol=2e-5, rtol=0)
