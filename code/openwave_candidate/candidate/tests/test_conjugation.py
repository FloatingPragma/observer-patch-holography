import unittest
from fractions import Fraction as F

from oph_openwave import conjugation as c


class ConjugationTests(unittest.TestCase):
    def test_exhaustive_full_branches_and_class_symmetry(self):
        self.assertEqual(c.branch_mismatches(), [])
        self.assertEqual(c.kernel_defect(F(1, 2))[0], F(0))
        self.assertTrue(all(c.paired_odd_mean(s, F(1, 2)) == 0 for s in c.STATES))

    def test_mutations_change_actual_targets_and_kernel(self):
        for failures in (c.branch_mismatches(cost=(0, 5, 8, 10, 8, 6)),
                         c.branch_mismatches(complement_coin=False)):
            self.assertTrue(any(row['target_mismatch'] for row in failures))
        self.assertEqual(c.kernel_defect(F(3, 4))[0], F(1, 2))
        self.assertTrue(any(c.paired_odd_mean(s, F(3, 4)) != 0 for s in c.STATES))

    def test_integer_load_conservation_and_valid_input(self):
        for s in c.STATES:
            for coin in (0, 1):
                target, (_, current, applied) = c.step(s, coin)
                self.assertEqual(sum(target[:2]), sum(s[:2]))
                self.assertEqual((target[2]+target[3]) % 6, (s[2]+s[3]) % 6)
                self.assertEqual(current, s[0]-target[0])
                self.assertTrue(applied or target[2:] == s[2:])
        with self.assertRaises(ValueError):
            c.step((0, 1, 0, 0), True)
        with self.assertRaises(ValueError):
            c.kernel_row((0, 1, 0, 0), F(3, 2))


if __name__ == '__main__':
    unittest.main()
