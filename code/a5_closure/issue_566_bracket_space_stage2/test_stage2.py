#!/usr/bin/env python3
"""Regression and boundary tests for the isolated issue #566 stage-two packet."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
PRODUCER = HERE / "certify.py"
VERIFIER = HERE / "verify.py"
SYSTEM = HERE / "a5_jacobi_system_reduction.json"
RECEIPT = HERE / "a5_jacobi_stage2.receipt.json"


class JacobiStage2Tests(unittest.TestCase):
    def test_generated_artifacts_are_current(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(PRODUCER), "--check"],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertEqual(summary["raw_jacobi_coordinates"], 2640)
        self.assertEqual(summary["equivariant_orbits"], 44)
        self.assertEqual(summary["independent_quadrics"], 38)
        self.assertEqual(summary["rank_split"], [11, 27])
        self.assertTrue(summary["full_classification_open"])

    def test_independent_verifier(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(VERIFIER)],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        self.assertTrue(summary["verified"])
        self.assertEqual(summary["jacobi_rank"], 38)
        self.assertEqual(summary["fixed_line_rank"], 11)
        self.assertEqual(summary["residual_rank"], 27)
        self.assertTrue(all(summary["independent_mutations"].values()))
        required = {
            "rehashed_false_full_classification_caught",
            "rehashed_false_preferred_algebra_caught",
            "rehashed_false_source_selection_caught",
            "rehashed_false_compactness_caught",
            "rehashed_firewall_change_caught",
            "rehashed_system_rank_change_caught",
            "rehashed_system_split_change_caught",
            "rehashed_system_field_change_caught",
            "rehashed_minimal_polynomial_change_caught",
            "rehashed_transform_determinant_change_caught",
            "rehashed_sector_dimension_change_caught",
            "rehashed_receipt_extra_key_caught",
            "rehashed_system_extra_key_caught",
        }
        self.assertTrue(required <= set(summary["independent_mutations"]))

    def test_system_reduction_counts(self) -> None:
        system = json.loads(SYSTEM.read_text(encoding="utf-8"))
        jacobi = system["jacobi_tensor"]
        reduction = system["fixed_line_reduction"]
        self.assertEqual(len(jacobi["representative_equations"]), 44)
        self.assertEqual(len(jacobi["reduced_integer_equations"]), 38)
        self.assertEqual(len(reduction["contracted_integer_equations_in_x"]), 11)
        self.assertEqual(len(reduction["residual_integer_equations_in_x"]), 27)
        self.assertEqual(reduction["full_rank_split"], [11, 27, 38])
        self.assertEqual(reduction["derivation_weight_arrangement"]["flat_count"], 28)
        self.assertEqual(
            reduction["derivation_only_subspace"]["scope"],
            "Jacobi vanishing only; compactness and source selection are not established",
        )
        self.assertEqual(
            reduction["derivation_weight_arrangement"]["scope"],
            "product-sector forced-zero stratification only; the rank-27 residual equations still apply on every flat",
        )
        self.assertIn("coefficient-row rank 38", jacobi["zero_locus_equivalence"])

    def test_claim_boundary_is_fail_closed(self) -> None:
        receipt = json.loads(RECEIPT.read_text(encoding="utf-8"))
        self.assertEqual(
            receipt["status"],
            "EXACT_JACOBI_SYSTEM_AND_FIXED_LINE_REDUCTION__FULL_CLASSIFICATION_OPEN",
        )
        self.assertTrue(receipt["target_firewall"]["enabled"])
        self.assertEqual(receipt["semantic_inputs"]["count"], 2)
        imports = receipt["target_firewall"]["audited_import_roots_by_file"]
        self.assertNotIn("subprocess", imports["certify.py"])
        self.assertNotIn("subprocess", imports["verify.py"])
        self.assertIn("subprocess", imports["test_stage2.py"])
        self.assertTrue(
            any(row["name"] == "promote_verifier_runtime_import_boundary" and row["passed"] for row in receipt["mutation_tests"])
        )
        producer_mutations = {row["name"] for row in receipt["mutation_tests"]}
        self.assertTrue(
            {
                "rehashed_false_full_classification_claim",
                "rehashed_false_preferred_algebra_claim",
                "rehashed_false_source_selection_claim",
                "rehashed_false_compactness_claim",
                "rehashed_semantic_path_change",
                "rehashed_firewall_change",
                "rehashed_system_rank_change",
                "rehashed_system_split_change",
                "rehashed_system_field_change",
                "rehashed_minimal_polynomial_change",
                "rehashed_transform_determinant_change",
                "rehashed_sector_dimension_change",
                "rehashed_receipt_extra_key",
                "rehashed_system_extra_key",
            }
            <= producer_mutations
        )
        self.assertFalse(any(receipt["not_attained"].values()))
        self.assertTrue(all(row["passed"] for row in receipt["mutation_tests"]))


if __name__ == "__main__":
    unittest.main()
