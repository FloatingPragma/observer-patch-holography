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
            "limited_subchain_audit_physical_interface_open",
        )
        interface = self.payload["discriminating_interface"]
        self.assertEqual(interface["interface_id"], "physical_scalar_attachment")
        self.assertEqual(interface["interface_class"], "conditional_open_interface")
        self.assertEqual(len(interface["declared_yukawa_channels"]), 3)
        self.assertIn("scalar two-point pole and residue observables", interface["classes"])

    def test_selected_subchain_scope_is_explicit(self) -> None:
        closure = self.payload["input_closure"]
        self.assertEqual(
            closure["grammar_countermodel_labels_not_physical_completions"],
            ["empty", "one_doublet_declared", "duplicate_doublet", "inert_doublet"],
        )
        self.assertEqual(len(closure["chain"]), 2)
        for entry in closure["chain"]:
            self.assertIn("carrier manifest", entry["declared_inputs"])
        self.assertIn(
            "not a physical scalar-sector non-identifiability theorem",
            closure["scope_limits"],
        )
        self.assertFalse(closure["dynamics_hash_semantics"]["interchangeable"])

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
            self.payload["controls"]["target_injection"]["code"],
            "EXPLICIT_SCALAR_INPUT",
        )
        self.assertEqual(
            self.payload["controls"]["scalar_key_mutation"]["code"],
            "EXPLICIT_SCALAR_INPUT",
        )
        self.assertEqual(
            sorted(self.payload["controls"]["wrong_multiplicity"]["countermodels"]),
            ["n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"],
        )

    def test_doctored_parent_pin_fails(self) -> None:
        carrier = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.CARRIER_MANIFEST_NAME
        )
        response = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.RESPONSE_ARTIFACT_NAME
        )
        pole = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.POLE_RESIDUE_ARTIFACT_NAME
        )
        pole["carrier_binding"]["parent_artifact_sha256"] = "sha256:" + "0" * 64
        pole["artifact_sha256"] = cert.artifact_self_hash(pole)
        with self.assertRaises(cert.CertificateError):
            cert.validate_selected_subchain(carrier, response, pole)

    def test_explicit_scalar_input_fails_even_when_rehashed(self) -> None:
        carrier = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.CARRIER_MANIFEST_NAME
        )
        response = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.RESPONSE_ARTIFACT_NAME
        )
        pole = cert.load_json(
            cert.MODULE_DIR / "manifests" / cert.POLE_RESIDUE_ARTIFACT_NAME
        )
        response["inputs"] = {"neutral_name": "inert_doublet_scalar"}
        response["artifact_sha256"] = cert.artifact_self_hash(response)
        pole["carrier_binding"]["parent_artifact_sha256"] = response["artifact_sha256"]
        pole["artifact_sha256"] = cert.artifact_self_hash(pole)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.validate_selected_subchain(carrier, response, pole)
        self.assertEqual(caught.exception.code, "EXPLICIT_SCALAR_INPUT")


if __name__ == "__main__":
    unittest.main()
