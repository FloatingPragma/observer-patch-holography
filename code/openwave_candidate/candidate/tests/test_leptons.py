import unittest
from fractions import Fraction as F

import mpmath as mp
import numpy as np

from oph_openwave import leptons


class LeptonTests(unittest.TestCase):
    def test_full_tick_and_loss_control(self):
        tick = leptons.register_tick()
        leptons.verify_tick(tick)
        tick[0] = 0
        with self.assertRaisesRegex(ArithmeticError, 'not unitary'):
            leptons.verify_tick(tick)

    def test_orientation_mutation_reaches_compression(self):
        leptons.verify_compression(leptons.register_tick())
        broken = leptons.register_tick(oriented=False)
        leptons.verify_tick(broken)
        with self.assertRaisesRegex(ArithmeticError, 'balanced circulant'):
            leptons.verify_compression(broken)

    def test_complete_read_probabilities(self):
        self.assertEqual(leptons.exact_read_counts(leptons.register_tick()),
                         [F(1, 12), F(1, 6), F(1, 72)])

    def test_matrix_diagonalization_independent_of_character_formula(self):
        theta = 2/9
        shift = np.roll(np.eye(3), 1, axis=0)
        matrix = np.eye(3) + (np.exp(1j*theta)*shift + np.exp(-1j*theta)*shift.T)/np.sqrt(2)
        actual = np.linalg.eigvalsh(matrix)
        with mp.workdps(60):
            character = sorted(float(v) for v in leptons.spectrum())
        np.testing.assert_allclose(actual, character, atol=2e-15, rtol=0)

    def test_frozen_reference_rejects_exact_muon_match(self):
        result = leptons.comparison()
        self.assertGreater(float(result['muon_reference_only_standardized_discrepancy']), 400)
        self.assertLess(float(result['tau_reference_only_discrepancy_high']), 1)
        self.assertAlmostEqual(float(result['koide_ratio']), 2/3, places=15)
        self.assertAlmostEqual(float(result['electron_calibrated_muon_GeV']) /
                               float(result['electron_calibrated_electron_GeV']),
                               float(result['muon_electron_ratio']), places=10)


if __name__ == '__main__':
    unittest.main()
