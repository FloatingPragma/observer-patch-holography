"""Independent full-spectrum checks of the exact scalar-cost derivatives."""
import unittest
from fractions import Fraction

import numpy as np
from scipy.linalg import eigh

from .response import (CENTRE_CLASS, calculate, channel, ground_pair, lepton,
                       quark, scalar_projector)


def numeric(matrix):
    return np.asarray(matrix.tolist(), dtype=float)


def complete_spectral_moment(family, scalar_cost):
    """Use every generalized target eigenvector, including zero-weight rows."""
    source = ground_pair()['block']
    source_gram = numeric(source['gram'])
    source_h = numeric(source['H'])+(scalar_cost-1)*numeric(ground_pair()['n'])
    source_energies, source_vectors = eigh(source_gram@source_h, source_gram)
    e, g = source_energies[0], source_vectors[:, 0]
    if family in quark.CHANNELS:
        blocks = [(quark.transfer(family), quark.target(family))]
    else:
        blocks = [(lepton.transfer(family, spin), lepton.target(family, 2*spin))
                  for spin in (0, 1)]
    moment, total = 0.0, 0.0
    for block, target in blocks:
        gram = numeric(block['gram'])
        h = numeric(block['H'])+(scalar_cost-1)*numeric(scalar_projector(target))
        energies, vectors = eigh(gram@h, gram)
        amplitudes = vectors.T@gram@numeric(block['map'])@g
        probabilities = np.abs(amplitudes)**2/len(blocks)
        total += float(probabilities.sum())
        moment += float(probabilities@((energies-e)**2))
    return total, moment


class ScalarCostResponse(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = calculate()

    def test_all_complete_spectra_and_ground_derivatives(self):
        step = 0.001
        for family in CENTRE_CLASS:
            with self.subTest(family=family):
                expected = channel(family)['complete_unconditioned_second_energy_transfer_moment']
                values = []
                for cost in (1-step, 1.0, 1+step):
                    weight, moment = complete_spectral_moment(family, cost)
                    self.assertAlmostEqual(weight, 0.5, places=10)
                    values.append(moment)
                center = sum(float(Fraction(x)) for x in expected['value_interval'])/2
                derivative = sum(float(Fraction(x)) for x in expected['derivative_interval'])/2
                self.assertAlmostEqual(values[1], center, places=9)
                self.assertAlmostEqual((values[2]-values[0])/(2*step), derivative, delta=3e-8)

    def test_class_labels_do_not_replace_matrix_response(self):
        rows = self.receipt['channels']
        self.assertEqual(rows['uc']['complete_unconditioned_second_energy_transfer_moment'],
                         rows['dc']['complete_unconditioned_second_energy_transfer_moment'])
        self.assertEqual(rows['uc']['centre_class_of_created_hole'], 4)
        self.assertEqual(rows['Q']['centre_class_of_created_hole'], 5)
        self.assertEqual(rows['ec']['centre_class_of_created_hole'], 0)
        self.assertEqual(rows['L']['centre_class_of_created_hole'], 3)


if __name__ == '__main__':
    unittest.main()
