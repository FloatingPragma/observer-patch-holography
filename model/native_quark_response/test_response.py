"""Independent literal CAR normalization and complete positive-law controls."""
from fractions import Fraction as F
import unittest

import sympy as s

from .quark import copy_alignment, elementary_matrix, exterior_data, occupations, target_module
from .local import CHANNELS, target, transfer
from .response import positive_response
from .full import windows


class QuarkResponse(unittest.TestCase):
    def test_literal_eight_mode_split(self):
        _, _, _, two, vector = exterior_data()
        v = s.zeros(16, 5)
        for i, mask in enumerate(two):
            for j in range(5):
                v[mask, j] = vector[i, j]
        keys = occupations(2, 5)
        source = s.zeros(256, 15)
        for col, key in enumerate(keys):
            for i in range(5):
                for j in range(5):
                    if tuple(int(k == i)+int(k == j) for k in range(5)) != key:
                        continue
                    coefficient = s.prod(s.factorial(n) for n in key)/s.sqrt(2)
                    for a in range(16):
                        for b in range(16):
                            source[a+16*b, col] += coefficient*v[a, i]*v[b, j]
        for create in (False, True):
            odd = tuple(mask for mask in range(16) if mask.bit_count() == (3 if create else 1))
            embed = s.zeros(256, 20)
            for i in range(5):
                for j, a in enumerate(odd):
                    for b in range(16):
                        embed[a+16*b, 4*i+j] = v[b, i]/s.sqrt(2)
            for mode in range(4):
                # This matrix uses all eight CAR bits directly.
                entries = {}
                for mask in range(256):
                    if bool(mask >> mode & 1) != create:
                        entries[mask ^ (1 << mode), mask] = (-1)**((mask & ((1 << mode)-1)).bit_count())
                literal = s.ImmutableSparseMatrix(256, 256, entries)
                local = elementary_matrix(mode, create).extract(odd, two)*vector
                derivative = s.zeros(20, 15)
                for col, key in enumerate(keys):
                    for i in range(5):
                        if not key[i]:
                            continue
                        remaining = list(key)
                        remaining[i] -= 1
                        position = remaining.index(1)
                        for j in range(4):
                            derivative[4*position+j, col] += key[i]*local[j, i]
                self.assertEqual(literal*source, embed*derivative)
                self.assertNotEqual(literal*source, s.sqrt(6)*embed*derivative)

    def test_actual_six_copy_frames_and_full_images(self):
        frames = copy_alignment()
        self.assertEqual(frames[0][1], s.eye(4))
        self.assertEqual(len(frames), 6)
        for create in (False, True):
            data = target_module(create)
            self.assertEqual(data['image_dimension'], 224)
            self.assertEqual(data['Casimir'], s.Rational(95, 2))
            self.assertEqual(data['number'], 13 if create else 11)

    def test_complete_scalar_lepton_kernels(self):
        for family, dimension in (('uc', 26), ('dc', 26), ('Q', 52)):
            block = transfer(family)
            self.assertEqual(block['H'].rows, dimension)
            self.assertEqual(block['full_weight_dimension'], 236)
            self.assertEqual(target(family)['H'], block['H'])

    def test_exact_positive_measures_and_computed_degeneracy(self):
        uc = positive_response('uc')
        dc = positive_response('dc')
        self.assertEqual(uc['moments'], dc['moments'])
        self.assertEqual(uc['numerator'], dc['numerator'])
        self.assertEqual(uc['denominator'], dc['denominator'])
        for family, count in (('uc', 13), ('dc', 13), ('Q', 25)):
            data = positive_response(family)
            self.assertEqual(len(data['records']), count)
            self.assertEqual(data['visible_support'].degree(), 13)
            self.assertEqual(sum(row['exactly_zero_residue'] for row in data['records']),
                             12 if family == 'Q' else 0)
            self.assertTrue(all(row['weight_interval'][0] > 0 or row['exactly_zero_residue']
                                for row in data['records']))
            self.assertTrue(all(row['difference_interval'][0] > 0 for row in data['records']))

    def test_same_H_complete_prepared_records(self):
        for family in CHANNELS:
            data = windows(family)
            self.assertEqual(data['success_reference_energy'], '4')
            self.assertEqual(data['exact_ground_pump_failure'], '1/2')
            self.assertEqual(len(data['complete_law']), 4)
            for record in data['records']:
                lo, hi = map(F, record['prepared_conditional_probability'])
                self.assertTrue(0 <= lo <= hi <= 1)
            self.assertGreater(max(F(r['original_attempt_probability_lower']) for r in data['records']), F(1, 10))
        self.assertGreater(F(windows('uc')['records'][0]['original_attempt_probability_lower']), F(49, 100))


if __name__ == '__main__':
    unittest.main()
