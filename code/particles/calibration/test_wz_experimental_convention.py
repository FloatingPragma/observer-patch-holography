#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import unittest

from wz_experimental_convention import (
    convention_jacobian,
    make_target_receipt,
    propagate_covariance,
    running_width_to_energy_pole,
)


ROOT = Path(__file__).resolve().parent


class ExperimentalConventionTests(unittest.TestCase):
    def test_exact_complex_root(self) -> None:
        mass, width = 91.1879, 2.4955
        got = running_width_to_energy_pole(mass, width)
        direct = mass * mass / complex(1.0, width / mass)
        readout = complex(got["mass_GeV"], -0.5 * got["width_GeV"]) ** 2
        self.assertLess(abs(direct - readout), 2e-12)
        self.assertGreater(got["mass_GeV"], got["sqrt_Re_s_pole_mass_GeV"])

    def test_jacobian_against_central_difference(self) -> None:
        mass, width = 80.3625, 2.14
        jacobian = convention_jacobian(mass, width)
        for column, (delta_mass, delta_width) in enumerate(((1e-5, 0.0), (0.0, 1e-5))):
            plus = running_width_to_energy_pole(mass + delta_mass, width + delta_width)
            minus = running_width_to_energy_pole(mass - delta_mass, width - delta_width)
            numerical = [
                (plus["mass_GeV"] - minus["mass_GeV"]) / 2e-5,
                (plus["width_GeV"] - minus["width_GeV"]) / 2e-5,
            ]
            for row in range(2):
                self.assertAlmostEqual(jacobian[row][column], numerical[row], places=8)

    def test_covariance_rejects_indefinite_input(self) -> None:
        with self.assertRaises(ValueError):
            propagate_covariance(
                [[1.0, 0.0], [0.0, 1.0]],
                [[1.0, 2.0], [2.0, 1.0]],
            )

    def test_frozen_target_receipt(self) -> None:
        packet = json.loads((ROOT / "wz_pdg_2026_target_fixture.json").read_text())
        receipt = make_target_receipt(packet)
        self.assertEqual(receipt["output_convention"], "s_pole=(M-i*Gamma/2)^2")
        self.assertFalse(receipt["joint_chi2_licensed"])
        self.assertTrue(all(value > 0.0 for value in receipt["output_sigma_GeV"]))
        self.assertAlmostEqual(receipt["W"]["sqrt_Re_s_pole_mass_GeV"], 80.3340217557, places=9)
        self.assertAlmostEqual(receipt["Z"]["sqrt_Re_s_pole_mass_GeV"], 91.1537725310, places=9)


if __name__ == "__main__":
    unittest.main()
