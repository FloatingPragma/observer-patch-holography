"""Independent identities for the released shear stacks and stress transition."""

import json
import unittest

import numpy as np
from scipy.integrate import quad

from measure import FROZEN, HERE, analyse, check_values, joint_edge_distance


class PressureTransferTest(unittest.TestCase):
    def test_archived_result_and_permutation(self):
        profiles = json.loads((HERE / "profiles.json").read_text())
        frozen = json.loads(FROZEN.read_text())
        results = analyse(profiles, frozen)
        check_values(results, json.loads((HERE / "results.json").read_text()))
        # Check the two-stack statistic with the covariance blocks explicitly.
        for name, data in profiles["splits"].items():
            joint = np.array(data["covariance"])
            for a, b in ((0, 6), (6, 0)):
                d = np.diff(np.array(data["tangential_esd"]), axis=0)[0]
                sign = 1 if a == 0 else -1
                i, j = slice(a, a + 6), slice(b, b + 6)
                cov = joint[i, i] + joint[j, j] - joint[i, j] - joint[j, i]
                template = (
                    np.sqrt(
                        frozen["fitted_sparc_a0_m_s2"]
                        * np.array(data["bin_centers_m_s2"])
                    )
                    / frozen["four_G_m_s2_per_Msun_pc2"]
                )
                white = np.linalg.cholesky(cov)
                x = np.linalg.solve(white, template[:, None])
                q = np.linalg.lstsq(x, np.linalg.solve(white, sign * d), rcond=None)[0][
                    0
                ]
                self.assertAlmostEqual(
                    q,
                    sign
                    * results[name]["difference_at_centers"]["envelope_high_minus_low"],
                    places=10,
                )

    def test_extreme_pair_selection(self):
        edges = 2.0 ** np.arange(4)
        k = np.sqrt(2.0)
        a0, four_g = 2.0, 3.0
        coeff = 4 / (np.pi * four_g)
        for q in (0.0, 1.0, 7.0):

            def mean(b):
                return coeff * b + q * np.sqrt(a0 * b) / four_g

            lo, hi = edges[:-1], edges[1:]
            np.testing.assert_allclose(mean(hi) - k * mean(lo), coeff * (hi - k * lo))
            for high_mix, low_mix in ((0, 1), (1, 0), (0.9, 0.2)):
                early = high_mix * mean(hi) + (1 - high_mix) * mean(lo)
                late = low_mix * mean(hi) + (1 - low_mix) * mean(lo)
                self.assertTrue(
                    np.all(early - k * late <= coeff * (hi - k * lo) + 1e-12)
                )

        y = np.r_[mean(1.1 * lo), mean(1.9 * lo)]
        result = joint_edge_distance(y, np.eye(len(y)), edges, a0, four_g)
        self.assertLess(result["any_nonnegative_envelope"]["distance_upper"], 1e-10)
        self.assertGreater(result["frozen"]["distance_lower"], 1.0)

    def test_conserved_tension_taper(self):
        # Continuous radial pressure with compact support and m = h*r up to Rt.
        h, u, r0, rt = 0.01, 0.4, 1.0, 2.0

        def rho(r):
            return h / (4 * np.pi * r**2)

        def radial(r):
            return -u * rho(r) * ((rt - r) / (rt - r0)) ** 2

        def phi_prime(r):
            return (h * r + 4 * np.pi * r**3 * radial(r)) / (r * r * (1 - 2 * h))

        def tangential(r):
            derivative = 2 * u * h * (rt - r) / (4 * np.pi * (rt - r0) ** 2)
            return (derivative + r * r * (rho(r) + radial(r)) * phi_prime(r)) / (2 * r)

        for r in np.linspace(r0, rt, 50):
            self.assertLessEqual(abs(tangential(r)), rho(r))
            self.assertGreaterEqual(rho(r) + radial(r), 0.0)
            self.assertGreaterEqual(phi_prime(r), 0.0)
        stress = quad(lambda r: 2 * r * tangential(r), r0, rt)[0]
        gravity = quad(lambda r: r * r * (rho(r) + radial(r)) * phi_prime(r), r0, rt)[0]
        self.assertAlmostEqual(stress - gravity, -r0 * r0 * radial(r0), places=13)
        mass = quad(lambda r: 4 * np.pi * r * r * rho(r), r0, rt)[0]
        self.assertGreaterEqual(mass, 2 * np.pi * r0**3 * abs(radial(r0)))
        self.assertGreaterEqual(rt / r0, np.exp(u / 2))


if __name__ == "__main__":
    unittest.main()
