from __future__ import annotations

import sys
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import collar_poisson_witness as witness  # noqa: E402
import validate_collar_poisson_certificate as validator  # noqa: E402


class CollarPoissonWitnessTests(unittest.TestCase):
    def test_bernoulli_sum_law_is_exact_binomial(self) -> None:
        law = witness.bernoulli_sum_law([Fraction(1, 3)] * 4)
        self.assertEqual(sum(law), 1)
        self.assertEqual(law[0], Fraction(16, 81))
        self.assertEqual(law[4], Fraction(1, 81))

    def test_one_dependent_law_normalizes_with_correct_mean(self) -> None:
        m, q = 12, 6
        law = witness.one_dependent_law(m, q)
        self.assertEqual(sum(law), 1)
        mean = sum(n * mass for n, mass in enumerate(law))
        self.assertEqual(mean, Fraction(m, q * q))

    def test_total_variation_zero_case(self) -> None:
        # Bernoulli(p) sum with one cell vs Poisson(p): TV = (1 - p) - e^-p + tail terms.
        tv = witness.total_variation([Fraction(1, 2), Fraction(1, 2)], Fraction(1, 2))
        self.assertGreater(tv, 0)
        self.assertLess(tv, Decimal("0.25"))

    def test_all_declared_families_sit_within_their_bounds(self) -> None:
        payload = witness.compute()
        self.assertTrue(payload["all_families_within_bound"])
        self.assertEqual(payload["issue"], 320)
        for family in payload["families"]:
            self.assertTrue(family["tv_le_bound"], family["family_id"])
            tv = Decimal(family["exact_tv"])
            bound = Decimal(family["bound_decimal"])
            self.assertLessEqual(tv, bound, family["family_id"])

    def test_refinement_shrinks_tv_within_each_family_line(self) -> None:
        payload = witness.compute()
        by_id = {f["family_id"]: Decimal(f["exact_tv"]) for f in payload["families"]}
        pairs = [
            (
                "independent_uniform_unit_sqrtx_1_1_m8",
                "independent_uniform_unit_sqrtx_1_1_m64",
            ),
            (
                "independent_nonuniform_z6_presence_m8",
                "independent_nonuniform_z6_presence_m64",
            ),
            ("mean_drift_unit_m8", "mean_drift_unit_m64"),
        ]
        for shallow, deep in pairs:
            self.assertLess(by_id[deep], by_id[shallow])

    def test_tracked_certificate_validates(self) -> None:
        cert_path = (
            Path(__file__).resolve().parents[1]
            / "certificates"
            / "collar_poisson_witness_certificate.json"
        )
        self.assertTrue(cert_path.is_file())
        self.assertEqual(validator.main(str(cert_path)), 0)


if __name__ == "__main__":
    unittest.main()
