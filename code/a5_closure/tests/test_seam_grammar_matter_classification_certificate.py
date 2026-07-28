#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #627 classification certificate."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import seam_grammar_matter_classification_certificate as cert  # noqa: E402


class SeamGrammarMatterClassificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()

    def test_reference_manifest_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(cert.MANIFEST_PATH)
        body = {k: v for k, v in stored.items() if k != "manifest_sha256"}
        self.assertEqual(body, self.payload)
        self.assertEqual(stored["manifest_sha256"], "sha256:" + cert.sha256_json(body))
        self.assertEqual(cert.verify_stored()["status"], "PASS")

    def test_rerun_writes_byte_identical_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_grammar_selection_is_measured_and_target_free(self) -> None:
        selection = self.payload["grammar_selection"]
        self.assertEqual(selection["selected_group"], "Z6")
        self.assertTrue(selection["target_free"])
        self.assertIn("Z7", selection["excluded"])

    def test_module_action_table_is_exact(self) -> None:
        actions = {
            row["subgroup_order"]: row
            for row in self.payload["module_action_classification"]["subgroup_actions"]
        }
        self.assertTrue(actions[6]["faithful_on_module"])
        self.assertEqual(actions[6]["fixed_subspace_dimension"], 1)
        self.assertEqual(actions[3]["fixed_subspace_dimension"], 3)
        self.assertEqual(actions[2]["fixed_subspace_dimension"], 7)
        self.assertEqual(actions[1]["fixed_subspace_dimension"], 15)

    def test_mechanism_rows_are_classified(self) -> None:
        rows = {r["module"]: r for r in self.payload["mechanism_classification"]["rows"]}
        for module in ("Z2 -> 1", "Z3 -> 1", "Z6 -> 1"):
            self.assertTrue(rows[module]["realized_module_action"]["faithful_on_module"])
        self.assertEqual(rows["id: S3 -> S3"]["realized_module_action"]["through"], "none")
        self.assertEqual(
            rows["Z2 -> Z4 (doubling)"]["realized_module_action"]["through"], "none"
        )

    def test_selection_interface_stays_open(self) -> None:
        interface = self.payload["matter_action_interface"]
        self.assertEqual(interface["id"], "physical_sector_mechanism_selection")
        self.assertEqual(interface["class"], "conditional_open_interface")
        self.assertEqual(
            self.payload["bounded_exit"],
            "classification_landed_with_named_selection_interface",
        )

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)

    def test_doctored_charges_are_refused(self) -> None:
        with self.assertRaises(cert.CertificateError):
            cert.module_action_table({"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 5})
        with self.assertRaises(cert.CertificateError):
            cert.module_action_table({"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/7": 6})


if __name__ == "__main__":
    unittest.main()
