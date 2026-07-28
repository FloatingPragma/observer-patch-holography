#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #623 blindness certificate."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import scalar_sector_blindness_certificate as cert  # noqa: E402


class ScalarSectorBlindnessTests(unittest.TestCase):
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

    def test_bounded_exit_and_interface(self) -> None:
        self.assertEqual(
            self.payload["bounded_exit"],
            "conditional_open_interface_with_countermodels_retained",
        )
        interface = self.payload["discriminating_interface"]
        self.assertEqual(interface["interface_id"], "physical_scalar_attachment")
        self.assertEqual(interface["interface_class"], "conditional_open_interface")
        self.assertEqual(len(interface["channels"]), 3)

    def test_input_closure_covers_all_completions(self) -> None:
        closure = self.payload["input_closure"]
        self.assertEqual(
            closure["completions_covered"],
            ["empty", "one_doublet_declared", "duplicate_doublet", "inert_doublet"],
        )
        self.assertEqual(len(closure["chain"]), 2)
        for entry in closure["chain"]:
            self.assertIn("carrier manifest", entry["inputs"])

    def test_countermodels_are_retained(self) -> None:
        retained = self.payload["countermodels"]
        self.assertIn("not source-determined", retained["grammar_statement"])
        self.assertEqual(
            len(retained["declared_completion"]["yukawa_channels"]), 3
        )

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        self.assertEqual(
            sorted(self.payload["controls"]["wrong_multiplicity"]["countermodels"]),
            ["n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"],
        )

    def test_scalar_key_detection_is_live(self) -> None:
        keys: set[str] = set()
        cert.collect_keys({"outer": {"scalar_mass": 1}}, keys)
        offending = [k for k in keys if any(t in k for t in cert.SCALAR_TOKENS)]
        self.assertTrue(offending)

    def test_doctored_parent_pin_fails(self) -> None:
        response = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.RESPONSE_ARTIFACT_NAME
        )
        keys: set[str] = set()
        cert.collect_keys(response, keys)
        for token in cert.SCALAR_TOKENS:
            self.assertFalse([k for k in keys if token in k], token)


if __name__ == "__main__":
    unittest.main()
