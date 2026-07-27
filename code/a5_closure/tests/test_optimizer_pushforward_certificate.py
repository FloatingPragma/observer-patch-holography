#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #619 pushforward certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import optimizer_pushforward_certificate as cert  # noqa: E402


class OptimizerPushforwardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = (
            MODULE_DIR / "manifests" / "optimizer_pushforward_reference.json"
        )
        cls.receipt_path = (
            MODULE_DIR / "receipts" / "optimizer_pushforward_reference.receipt.json"
        )
        cls.negative_path = (
            MODULE_DIR / "negative_controls" / "issue_619_negative_controls.json"
        )
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_payload_is_deterministic(self) -> None:
        again = cert.certificate_payload(self.manifest)
        self.assertEqual(again, self.expected)
        self.assertEqual(cert.sha256_json(again), cert.sha256_json(self.expected))

    def test_schema_and_verdict_fields(self) -> None:
        self.assertEqual(self.manifest["schema"], cert.SCHEMA)
        self.assertEqual(self.expected["schema"], cert.RECEIPT_SCHEMA)
        self.assertEqual(self.expected["issue"], 619)
        verdict = self.expected["verdict"]
        self.assertEqual(verdict["positive_lane"], "exact_sufficiency_family")
        self.assertEqual(verdict["negative_lane"], "exact_nonclosure_defect_5_66")
        self.assertEqual(verdict["axiom_status"], "conditional_open_interface")
        self.assertEqual(verdict["bounded_exit"], "conditional_open_interface")
        self.assertEqual(len(verdict["sufficiency_conditions"]), 2)
        self.assertIn("not an instance of", verdict["gravity_ladder_composition"])

    def test_channel_pairing_is_sourced_from_565_receipt(self) -> None:
        channel = self.expected["refinement_channel"]
        self.assertEqual(
            channel["axis_pairs"],
            [[0, 3], [1, 2], [4, 7], [5, 6], [8, 11], [9, 10]],
        )
        self.assertEqual(
            self.expected["axis_receipt_sha256"],
            self.manifest["axis_receipt_sha256"],
        )
        references = self.expected["references"]
        self.assertTrue(references["coarse_equals_pushforward_of_fine"])
        self.assertEqual(references["fine"], ["1/12"] * 12)
        self.assertEqual(references["coarse"], ["1/6"] * 6)

    def test_positive_family_battery(self) -> None:
        family = self.expected["positive_sufficiency_family"]
        self.assertGreaterEqual(family["battery_size"], 5)
        self.assertGreaterEqual(family["fully_asymmetric_targets"], 2)
        self.assertEqual(family["free_rational_parameters"], 5)
        for row in family["battery"]:
            self.assertTrue(row["commutes"])
            self.assertEqual(
                row["pushforward_of_fine_projection"], row["coarse_target"]
            )
            self.assertEqual(row["coarse_projection"], row["coarse_target"])
            target = [Fraction(x) for x in row["coarse_target"]]
            projection = [Fraction(x) for x in row["fine_projection"]]
            pairs = self.expected["refinement_channel"]["axis_pairs"]
            for axis, (i, j) in enumerate(pairs):
                self.assertEqual(projection[i], target[axis] / 2)
                self.assertEqual(projection[j], target[axis] / 2)

    def test_nonclosure_defect_is_exactly_5_66(self) -> None:
        countermodel = self.expected["nonclosure_countermodel"]
        self.assertEqual(countermodel["nonclosure_defect"], "5/66")
        self.assertEqual(countermodel["pushed_axis_weight"], "8/33")
        self.assertEqual(countermodel["free_face_value"], "5/66")
        self.assertEqual(
            Fraction(countermodel["nonclosure_defect"]),
            Fraction(8, 33) - Fraction(1, 6),
        )
        projection = [Fraction(x) for x in countermodel["fine_projection"]]
        self.assertEqual(projection[0], Fraction(1, 6))
        self.assertEqual(projection[1:], [Fraction(5, 66)] * 11)
        self.assertEqual(sum(projection), 1)
        pushed = [Fraction(x) for x in countermodel["pushforward_of_fine_projection"]]
        self.assertEqual(pushed[0], Fraction(8, 33))
        self.assertEqual(pushed[1:], [Fraction(5, 33)] * 5)
        self.assertEqual(
            countermodel["coarse_projection"], ["1/6"] * 6
        )

    def test_countermodel_image_and_non_expressibility(self) -> None:
        countermodel = self.expected["nonclosure_countermodel"]
        battery = countermodel["image_surjectivity_battery"]
        self.assertGreaterEqual(len(battery), 5)
        pairs = self.expected["refinement_channel"]["axis_pairs"]
        for row in battery:
            q = [Fraction(x) for x in row["coarse_state"]]
            preimage = [Fraction(x) for x in row["preimage"]]
            self.assertEqual(preimage[0], Fraction(1, 6))
            self.assertGreaterEqual(q[0], Fraction(1, 6))
            for axis, (i, j) in enumerate(pairs):
                self.assertEqual(preimage[i] + preimage[j], q[axis])
        witness = countermodel["no_preimage_witness"]
        self.assertLess(Fraction(witness["coarse_state"][0]), Fraction(1, 6))
        expressibility = countermodel["non_expressibility_witness"]
        feasible = [Fraction(x) for x in expressibility["feasible_state"]]
        infeasible = [
            Fraction(x)
            for x in expressibility["infeasible_state_with_equal_pushforward"]
        ]
        for i, j in pairs:
            self.assertEqual(feasible[i] + feasible[j], infeasible[i] + infeasible[j])
        self.assertEqual(feasible[0], Fraction(1, 6))
        self.assertNotEqual(infeasible[0], Fraction(1, 6))

    def test_control_battery_fails_closed_with_exact_witnesses(self) -> None:
        controls = self.expected["control_battery"]
        wrong_reference = controls["wrong_reference"]
        self.assertTrue(wrong_reference["fails_closed"])
        self.assertEqual(wrong_reference["defect_on_axis_1"], "2/15")
        self.assertNotEqual(
            wrong_reference["fine_then_push"], wrong_reference["push_then_optimize"]
        )
        changed_weight = controls["changed_weight"]
        self.assertTrue(changed_weight["fails_closed"])
        invisible = changed_weight["t_invisible_tilt"]
        self.assertTrue(invisible["commutation_survives"])
        self.assertTrue(invisible["equal_split_formula_breaks"])
        projection = [Fraction(x) for x in invisible["fine_projection"]]
        self.assertEqual(projection[0], Fraction(2, 9))
        self.assertEqual(projection[3], Fraction(1, 9))
        mismatched = changed_weight["mismatched_pair"]
        self.assertFalse(mismatched["coarse_reference_is_pushforward"])
        self.assertEqual(mismatched["defect_on_axis_1"], "2/15")
        incomplete = controls["incomplete_constraint"]
        self.assertTrue(incomplete["fails_closed"])
        self.assertEqual(incomplete["defect_on_axis_5"], "5/42")
        self.assertEqual(
            Fraction(incomplete["defect_on_axis_5"]),
            Fraction(6, 21) - Fraction(1, 6),
        )
        channel = controls["non_sufficient_channel"]
        self.assertTrue(channel["fails_closed"])
        self.assertEqual(channel["defect"], "5/66")
        self.assertEqual(channel["cross_reference"], "nonclosure_countermodel")

    def test_no_floats_anywhere_in_the_receipt(self) -> None:
        def walk(value: object) -> None:
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for key, item in value.items():
                    walk(key)
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)

        walk(self.expected)

    def test_negative_controls_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        stored = cert.load_json(self.negative_path)
        self.assertEqual(payload, stored)
        names = {row["name"] for row in payload["finite_controls"]}
        self.assertEqual(
            names,
            {
                "wrong_schema",
                "float_arithmetic_declaration",
                "wrong_axis_receipt_pin",
                "wrong_axis_receipt_target",
                "float_tolerance_injection",
                "approximate_projection_injection",
                "quantum_state_space_injection",
                "declared_defect_injection",
            },
        )
        self.assertTrue(all(row["passed"] for row in payload["finite_controls"]))

    def test_manifest_mutations_are_rejected(self) -> None:
        for name, mutant, expected_code in cert.negative_control_cases(self.manifest):
            with self.assertRaises(cert.CertificateError, msg=name) as caught:
                cert.certificate_payload(mutant)
            self.assertEqual(caught.exception.code, expected_code, name)

    def test_tampered_receipt_is_rejected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["nonclosure_countermodel"]["nonclosure_defect"] = "0"
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(caught.exception.code, "RECEIPT_MISMATCH")

    def test_tampered_pairing_is_rejected(self) -> None:
        mutant = copy.deepcopy(dict(self.manifest))
        mutant["axis_receipt_path"] = "manifests/optimizer_pushforward_reference.json"
        mutant["axis_receipt_sha256"] = cert.sha256_json(mutant)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertIn(caught.exception.code, {"UPSTREAM_HASH", "UPSTREAM_RECEIPT"})

    def test_projection_helpers_are_exact(self) -> None:
        reference = cert.uniform(4)
        blocks = ([0, 1], [2, 3])
        masses = (Fraction(1, 3), Fraction(2, 3))
        projection = cert.block_mass_projection(reference, blocks, masses)
        self.assertEqual(
            projection,
            (Fraction(1, 6), Fraction(1, 6), Fraction(1, 3), Fraction(1, 3)),
        )
        for block in blocks:
            self.assertTrue(
                cert.conditional_matches_reference(projection, reference, block)
            )
        tilted = (Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8))
        tilted_projection = cert.block_mass_projection(tilted, blocks, masses)
        self.assertEqual(
            tilted_projection,
            (Fraction(2, 9), Fraction(1, 9), Fraction(1, 3), Fraction(1, 3)),
        )
        self.assertFalse(
            cert.conditional_matches_reference(projection, tilted, [0, 1])
        )


if __name__ == "__main__":
    unittest.main()
