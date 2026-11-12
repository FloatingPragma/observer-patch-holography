import unittest
from fractions import Fraction

from oph_openwave import charge


class ChargeTests(unittest.TestCase):
    def test_complete_exact_weights_and_counterexample(self):
        result = charge.run()
        self.assertTrue(result["computational_pass"])
        self.assertFalse(result["full_criterion_met"])
        self.assertEqual(result["centre_multiplicities"], [4, 6, 6, 4, 6, 6])
        self.assertEqual(result["equal_class_halves"]["0"], [0, 3])
        self.assertEqual(sum(Fraction(r["Q"]) for r in result["algebraically_allowed_coloured_daughters"]), -1)

    def test_bad_normalization_is_detected(self):
        with self.assertRaisesRegex(ValueError, "normalization"):
            charge.verify_charge_lattice(charge.enumerate_code(electric_scale=Fraction(1, 2)))

    def test_deleted_outcome_is_detected(self):
        with self.assertRaisesRegex(ValueError, "complete"):
            charge.verify_charge_lattice(charge.enumerate_code()[:-1])


if __name__ == "__main__":
    unittest.main()
