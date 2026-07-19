#!/usr/bin/env python3
from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import unittest

from operational_readback_contract import (
    compare_with_public_capacity,
    discrimination_error,
    example_receipt,
    finest_stable_scale,
    public_error,
    scalarize_rho_fiber,
)


HERE = Path(__file__).resolve().parent


class OperationalRhoEstimatorTests(unittest.TestCase):
    def test_discrimination_error_extremes(self) -> None:
        self.assertEqual(
            discrimination_error(
                [Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]
            ),
            0,
        )
        self.assertEqual(
            discrimination_error(
                [Fraction(1, 2), Fraction(1, 2)],
                [Fraction(1, 2), Fraction(1, 2)],
            ),
            Fraction(1, 2),
        )

    def test_finest_scale_requires_all_coarser_scales(self) -> None:
        errors = {
            Fraction(1, 10): Fraction(1, 2),
            Fraction(1, 5): Fraction(1, 100),
            Fraction(2, 5): Fraction(0),
        }
        self.assertEqual(
            finest_stable_scale(errors, Fraction(1, 10)), Fraction(1, 5)
        )

    def test_public_error_union_bound_is_capped(self) -> None:
        rows = [
            {
                "p0": [Fraction(1, 2), Fraction(1, 2)],
                "p1": [Fraction(1, 2), Fraction(1, 2)],
                "distribution_stage": "pre_checkpoint",
                "checkpoint_error": Fraction(3, 4),
                "reread_error": Fraction(3, 4),
            }
        ]
        self.assertEqual(public_error(rows), 1)

    def test_post_checkpoint_errors_cannot_be_double_counted(self) -> None:
        rows = [
            {
                "p0": [Fraction(1), Fraction(0)],
                "p1": [Fraction(0), Fraction(1)],
                "distribution_stage": "post_checkpoint",
                "checkpoint_error": Fraction(1, 100),
                "reread_error": Fraction(0),
            }
        ]
        with self.assertRaisesRegex(ValueError, "zero separately added"):
            public_error(rows)

    def test_rho_fiber_requires_completeness_and_agreement(self) -> None:
        ambiguous = scalarize_rho_fiber(
            [Fraction(1, 5), Fraction(2, 5)], fiber_manifest_complete=True
        )
        incomplete = scalarize_rho_fiber(
            [Fraction(1, 5)], fiber_manifest_complete=False
        )
        self.assertEqual(ambiguous["status"], "AMBIGUOUS_RHO_ESTIMATE")
        self.assertEqual(incomplete["status"], "INCOMPLETE_TERMINAL_FIBER")

    def test_rho_never_defines_F(self) -> None:
        comparison = compare_with_public_capacity(
            4, Fraction(1, 5), rho_derived_from_capacity=False
        )
        tautology = compare_with_public_capacity(
            4, Fraction(1, 5), rho_derived_from_capacity=True
        )
        self.assertEqual(comparison["status"], "PASS")
        self.assertFalse(comparison["defines_F"])
        self.assertEqual(tautology["status"], "RHO_DERIVED_FROM_CAPACITY")

    def test_example_receipt_is_current_and_nonproducing(self) -> None:
        receipt = example_receipt()
        runtime = json.loads(
            (
                HERE
                / "runtime"
                / "operational_capacity_readback_contract_example.json"
            ).read_text()
        )
        self.assertEqual(runtime, receipt)
        self.assertFalse(receipt["defines_F"])
        self.assertFalse(receipt["moves_cl7"])
        self.assertFalse(receipt["capacity_input_accessed_by_estimator"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
