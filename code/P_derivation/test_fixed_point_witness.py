#!/usr/bin/env python3
"""Smoke tests for the P/alpha fixed-point witness surface."""

from __future__ import annotations

from decimal import Decimal
import unittest

from paper_math import CODATA_2022_ALPHA_INV, PaperMathContext, build_fixed_point_witness


class FixedPointWitnessTests(unittest.TestCase):
    def test_codata_alpha_maps_to_expected_pixel_ratio(self) -> None:
        ctx = PaperMathContext(precision=24, su2_cutoff=4, su3_cutoff=4)
        p_value = ctx.observed_p_from_alpha_inv(CODATA_2022_ALPHA_INV)

        self.assertAlmostEqual(float(p_value), 1.6309682094039595, places=15)

    def test_witness_keeps_codata_as_compare_only_metadata(self) -> None:
        witness = build_fixed_point_witness(
            precision=10,
            mode="mz_anchor",
            su2_cutoff=6,
            su3_cutoff=4,
            scan_points=8,
            max_iterations=3,
            derivative_step="0.0001",
            sample_points=1,
        )

        self.assertEqual(witness["claim_status"], "numerical_witness_not_interval_certificate")
        self.assertEqual(witness["codata_2022_compare_only"]["alpha_inv"], str(CODATA_2022_ALPHA_INV))
        self.assertGreater(Decimal(witness["finite_difference"]["max_abs_sample_slope"]), Decimal("0"))
        self.assertIn("fixed_point_minus_codata_alpha_inv", witness["codata_2022_compare_only"])


if __name__ == "__main__":
    unittest.main()
