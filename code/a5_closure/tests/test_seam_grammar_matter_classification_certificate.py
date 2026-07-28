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

    def test_grammar_branch_keeps_exhaustiveness_open(self) -> None:
        branch = self.payload["grammar_branch"]
        self.assertEqual(branch["named_branch"], "Z6")
        self.assertTrue(branch["target_free"])
        self.assertEqual(branch["general_grammar_classification"], "open")
        self.assertIn("does not prove", branch["z7_result_scope"])

    def test_hypercharge_character_menu_is_complete(self) -> None:
        menus = {
            row["group_order"]: row
            for row in self.payload["hypercharge_character_menu"]["groups"]
        }
        self.assertEqual(menus[2]["fixed_dimension_menu"], [7, 15])
        self.assertEqual(menus[3]["fixed_dimension_menu"], [3, 15])
        self.assertEqual(menus[6]["fixed_dimension_menu"], [1, 3, 7, 15])
        self.assertEqual(len(menus[6]["characters"]), 6)

    def test_diagonal_kernel_fixes_all_realized_matter(self) -> None:
        diagonal = self.payload["diagonal_kernel_action"]
        self.assertEqual(diagonal["generator_color_weak_hypercharge"], [1, 1, 1])
        self.assertFalse(diagonal["faithful_on_module"])
        self.assertEqual(diagonal["fixed_subspace_dimension"], 15)
        self.assertTrue(all(row["phase_sixths"] == 0 for row in diagonal["fields"]))

    def test_mechanism_rows_are_classified(self) -> None:
        rows = {r["module"]: r for r in self.payload["mechanism_classification"]["rows"]}
        for module in ("Z2 -> 1", "Z3 -> 1", "Z6 -> 1"):
            action = rows[module]["realized_module_action"]
            self.assertEqual(action["status"], "supplied_character_required")
            self.assertNotIn("selected_character_exponent", action)
        self.assertEqual(rows["id: S3 -> S3"]["realized_module_action"]["through"], "none")
        self.assertEqual(
            rows["Z2 -> Z4 (doubling)"]["realized_module_action"]["through"], "none"
        )

    def test_selection_interface_stays_open(self) -> None:
        interface = self.payload["matter_action_interface"]
        self.assertEqual(interface["id"], "physical_sector_mechanism_selection")
        self.assertEqual(interface["class"], "conditional_open_interface")
        self.assertEqual(interface["owner_issue"], 569)
        self.assertEqual(
            self.payload["bounded_exit"],
            "exact_named_character_and_diagonal_action_classification",
        )

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)

    def test_doctored_charges_are_refused(self) -> None:
        with self.assertRaises(cert.CertificateError):
            cert.hypercharge_character_menu(
                {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/6": 5}
            )
        with self.assertRaises(cert.CertificateError):
            cert.hypercharge_character_menu(
                {"-1/2": 2, "-2/3": 3, "1": 1, "1/3": 3, "1/7": 6}
            )

    def test_nonidentity_refinement_is_refused_even_when_rehashed(self) -> None:
        artifact = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.RESPONSE_ARTIFACT_NAME
        )
        row = artifact["physical_refinement_maps"]["port_persistence_maps"][0]
        row["port_map"][0], row["port_map"][1] = row["port_map"][1], row["port_map"][0]
        map_body = {key: value for key, value in row.items() if key != "map_hash"}
        row["map_hash"] = "sha256:" + cert.sha256_json(map_body)
        body = {
            key: value for key, value in artifact.items()
            if key != "artifact_sha256"
        }
        artifact["artifact_sha256"] = "sha256:" + cert.sha256_json(body)
        with self.assertRaises(cert.CertificateError):
            cert.validate_refinement_artifact(artifact)


if __name__ == "__main__":
    unittest.main()
