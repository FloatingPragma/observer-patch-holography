#!/usr/bin/env python3
"""Regression and adversarial tests for the kinetic-form selection receipt."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import kinetic_form_selection_certificate as cert  # noqa: E402

pcc = cert.pcc


class KineticFormSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.receipt = cert.build_receipt()
        cls.current = cert.rebuild_current()
        cls.gram = cert.hs_gram(cls.current["generators"])
        cls.table = cert.structure_constants(
            cls.current["generators"], cls.current["flat"]
        )

    def test_runtime_receipt_is_exactly_recomputable(self) -> None:
        stored = json.loads(cert.RECEIPT_PATH.read_text(encoding="ascii"))
        self.assertEqual(stored, self.receipt)
        cert.verify_runtime()

    def test_ad_invariance_verified_on_all_triples(self) -> None:
        self.assertEqual(
            self.receipt["ad_invariance"]["verified_basis_triples"], 12 ** 3
        )

    def test_killing_relative_coefficients(self) -> None:
        row = self.receipt["killing_relative_coefficients"]
        self.assertEqual(row["su2"], "1")
        self.assertEqual(row["su3"], "1/6")
        self.assertEqual(row["ratio_su2_over_su3"], "6")

    def test_dichotomy_is_exact(self) -> None:
        row = self.receipt["dichotomy"]
        self.assertEqual(row["port_response_ratio_su2_over_su3"], "6")
        self.assertEqual(row["matter_trace_ratio_su2_over_su3"], "3/2")
        self.assertTrue(row["branches_exactly_distinct"])

    def test_port_metric_band_data_retained_at_corrected_typing(self) -> None:
        row = self.receipt["port_metric_band_data"]
        self.assertEqual(row["band_coefficients"]["unit_band"], "1/4")
        self.assertEqual(row["band_coefficients"]["frame_band"], "5 + 1*sqrt(5)")
        self.assertEqual(row["band_coefficients"]["kernel_band"], "5 + -1*sqrt(5)")
        self.assertEqual(row["band_coefficients"]["quintet_band"], "3 + 1*sqrt(5)")
        self.assertEqual(row["su3_dimension_weighted_average"], "15/4 + 1*sqrt(5)")
        self.assertIn("not a canonical invariant coefficient", row["typing"])

    def test_matter_branch_cofactors_frozen(self) -> None:
        row = self.receipt["matter_trace_branch"]["frozen_rg_statistic"]
        self.assertEqual(row["exact_cofactors"], ["-23/3", "37", "-218/9"])
        self.assertEqual(row["integer_zero_locus"], "69 x1 - 333 x2 + 218 x3 = 0")
        self.assertEqual(row["kinetic_column_k"], ["10/3", "2", "2"])
        self.assertEqual(row["beta_column_b"], ["41/6", "-19/6", "-7"])

    def test_comparison_boundary_sealed(self) -> None:
        row = self.receipt["comparison_boundary"]
        self.assertFalse(row["public_measurement_read"])
        self.assertFalse(row["comparison_permitted"])

    def test_superseded_typing_names_the_withdrawn_claim(self) -> None:
        row = self.receipt["superseded_typing"]
        self.assertIn("not an ad-invariant kinetic form on su(3)", row["withdrawn"])
        self.assertIn("Killing-form multiple", row["replacement"])

    def test_no_float_values_anywhere(self) -> None:
        def walk(node: object) -> None:
            if isinstance(node, float):
                self.fail("receipt contains a float value")
            elif isinstance(node, dict):
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(self.receipt)

    def test_mutated_gram_fails_ad_invariance(self) -> None:
        mutated = [row[:] for row in self.gram]
        mutated[2][7] = mutated[2][7] + pcc.ONE
        mutated[7][2] = mutated[2][7]
        with self.assertRaises(cert.KineticSelectionError):
            cert.check_ad_invariance(mutated, self.table)

    def test_mutated_table_fails_ad_invariance_or_closure(self) -> None:
        mutated = [
            [[entry for entry in row] for row in plane] for plane in self.table
        ]
        mutated[0][1][5] = mutated[0][1][5] + pcc.ONE
        mutated[1][0][5] = -mutated[0][1][5]
        with self.assertRaises(cert.KineticSelectionError):
            cert.check_ad_invariance(self.gram, mutated)

    def test_killing_proportionality_guard(self) -> None:
        bases = cert.band_field_bases(self.current["frame"])
        ideals = cert.verify_ideal_structure(self.table, bases)
        su3_b = cert.form_gram_on(self.gram, ideals["su3"])
        su3_k = cert.killing_gram_on(self.table, ideals["su3"])
        mutated = [row[:] for row in su3_b]
        mutated[0][0] = mutated[0][0] + pcc.ONE
        with self.assertRaises(cert.KineticSelectionError):
            cert.killing_relative_coefficient(mutated, su3_k, "su3")

    def test_matter_branch_cofactor_arithmetic_is_forced(self) -> None:
        k = (Fraction(10, 3), Fraction(2), Fraction(2))
        b = (Fraction(41, 6), Fraction(-19, 6), Fraction(-7))
        det = (
            k[1] * b[2] - k[2] * b[1],
            -(k[0] * b[2] - k[2] * b[0]),
            k[0] * b[1] - k[1] * b[0],
        )
        self.assertEqual(det, (Fraction(-23, 3), Fraction(37), Fraction(-218, 9)))
        x = (Fraction(3), Fraction(5), Fraction(7))
        direct = sum(c * v for c, v in zip(det, x, strict=True))
        integer_form = 69 * x[0] - 333 * x[1] + 218 * x[2]
        self.assertEqual(direct * Fraction(-9), integer_form)

    def test_parent_pins_match_disk(self) -> None:
        for pin in self.receipt["parent_pins"]:
            path = cert.REPO_ROOT / pin["path"]
            data = path.read_bytes()
            self.assertEqual(len(data), pin["bytes"])
            self.assertEqual(cert.tagged_sha256(data), pin["sha256"])


if __name__ == "__main__":
    unittest.main()
