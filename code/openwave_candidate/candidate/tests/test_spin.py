import unittest

import sympy as s

from oph_openwave import spin


class SpinTests(unittest.TestCase):
    def test_exact_kernel_and_complete_exchange(self):
        result = spin.run()
        self.assertTrue(result["computational_pass"])
        self.assertFalse(result["full_criterion_met"])
        self.assertEqual(result["Casimir_multiplicities"], {"3/4": 2, "15/4": 4})
        self.assertEqual(result["exchange_basis_comparisons"], 147456)

    def test_central_phase_cannot_be_quotiented_out(self):
        q = spin.exterior_rotation((s.eye(2)-s.I*spin.PAULI[2])/s.sqrt(2))
        with self.assertRaisesRegex(ValueError, "360-degree"):
            spin.verify_turns((1+s.I)/s.sqrt(2)*q)

    def test_commuting_statistics_fail_the_same_exchange_check(self):
        with self.assertRaisesRegex(ValueError, "exchange sign"):
            spin.verify_exchange(fermionic=False, full_domain=False)


if __name__ == "__main__":
    unittest.main()
