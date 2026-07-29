#!/usr/bin/env python3
"""Regression and adversarial tests for the bounded #631 precursor."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import check_common_ew_order_unit_carrier as checker  # noqa: E402
import common_ew_order_unit_carrier_certificate as cert  # noqa: E402


def rehash(manifest: dict) -> dict:
    body = {key: value for key, value in manifest.items() if key != "manifest_sha256"}
    manifest["manifest_sha256"] = "sha256:" + cert.sha256_json(body)
    return manifest


class CommonEWOrderUnitCarrierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()
        cls.manifest = cert.build_manifest()
        cls.stored_manifest = cert.load_json(cert.MANIFEST_PATH)
        cls.stored_receipt = cert.load_json(checker.RECEIPT_PATH)

    def test_manifest_is_byte_exact_rebuild(self) -> None:
        self.assertEqual(self.stored_manifest, self.manifest)
        self.assertEqual(cert.verify_stored()["status"], "PASS")

    def test_independent_receipt_is_exact_rebuild(self) -> None:
        expected = checker.build_receipt(self.stored_manifest)
        self.assertEqual(self.stored_receipt, expected)
        checker.verify_receipt(self.stored_manifest, self.stored_receipt)
        replay = expected["independent_replay"]
        self.assertEqual(replay["screen_fixed_space_dimension"], 1)
        self.assertEqual(replay["weak_doublet_copies"], 4)
        self.assertEqual(replay["weak_multiplicity_center_dimension"], 1)

    def test_screen_line_is_unique_and_refinement_natural(self) -> None:
        screen = self.payload["screen_order_unit_line"]
        self.assertEqual(screen["basis"], [1] * 12)
        self.assertEqual(screen["fixed_space_dimension"], 1)
        self.assertEqual(screen["rotation_checks"], 60)
        self.assertEqual(screen["orbit_of_p00"], list(range(12)))
        refinement = screen["refinement"]
        self.assertEqual(refinement["declared_nonidentity_maps"], 3)
        self.assertEqual(refinement["checked_cocycle_triangles"], 1)
        self.assertTrue(refinement["order_unit_natural"])
        self.assertTrue(all(row["order_unit_fixed"] for row in refinement["maps"]))

    def test_weak_multiplicity_center_is_exact(self) -> None:
        electroweak = self.payload["electroweak_order_unit_line"]
        copies = electroweak["copy_derivation"]
        self.assertEqual(copies["Q_weak_doublet_copies"], 3)
        self.assertEqual(copies["L_weak_doublet_copies"], 1)
        self.assertEqual(copies["total_weak_doublet_copies"], 4)
        self.assertEqual(electroweak["multiplicity_algebra"], "End(C^4)")
        self.assertEqual(electroweak["multiplicity_algebra_dimension"], 16)
        self.assertEqual(electroweak["central_constraint_rank"], 15)
        self.assertEqual(electroweak["center_dimension"], 1)

    def test_intertwiner_is_exact_but_nonphysical(self) -> None:
        intertwiner = self.payload["local_order_unit_intertwiner"]
        self.assertEqual(intertwiner["formula"], "T(c * 1_12) = c * I_4")
        self.assertTrue(intertwiner["linear"])
        self.assertTrue(intertwiner["positive"])
        self.assertTrue(intertwiner["unital"])
        self.assertTrue(intertwiner["normalized_trace_preserving"])
        self.assertTrue(
            intertwiner["unique_among_unital_linear_maps_between_the_two_lines"]
        )
        self.assertIn("no physical common-load identity", intertwiner["scope"])

    def test_capacity_scalar_and_physical_promotions_stay_false(self) -> None:
        firewall = self.payload["source_firewall"]
        self.assertEqual(
            firewall["allowed_issue_inputs"], [565, 566, 314, 567, 616]
        )
        self.assertNotIn(505, firewall["allowed_issue_inputs"])
        self.assertFalse(firewall["issue_505_consumed"])
        self.assertFalse(firewall["cosmic_capacity_consumed"])
        self.assertFalse(firewall["physical_scalar_assumed"])
        self.assertFalse(firewall["physical_common_carrier_assumed"])
        self.assertTrue(
            all(value is False for value in self.payload["promotion"].values())
        )
        self.assertFalse(self.stored_receipt["promotion_allowed"])

    def test_scalar_nonselection_boundary_is_retained(self) -> None:
        boundary = self.payload["electroweak_order_unit_line"]["scalar_boundary"]
        self.assertEqual(boundary["scalar_existence"], "not_source_determined")
        self.assertEqual(boundary["scalar_multiplicity"], "independence_limited")
        self.assertEqual(boundary["one_doublet_completion"], "declared_only")
        self.assertFalse(boundary["physical_scalar_selected"])
        self.assertEqual(
            self.payload["open_gates"][0],
            {
                "gate": "source_selected_scalar_carrier",
                "status": "open",
                "owners": [636],
            },
        )

    def test_every_control_fails_closed(self) -> None:
        self.assertEqual(len(self.payload["controls"]), 6)
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        self.assertEqual(
            self.payload["controls"]["capacity_injection"]["code"],
            "CAPACITY_INPUT_FORBIDDEN",
        )
        self.assertEqual(
            self.payload["controls"]["physical_common_carrier_promotion"]["code"],
            "PHYSICAL_COMMON_CARRIER_OPEN",
        )

    def test_rehashed_physical_promotion_is_rejected(self) -> None:
        tampered = json.loads(json.dumps(self.stored_manifest))
        tampered["promotion"]["promotion_allowed"] = True
        rehash(tampered)
        with self.assertRaises(cert.CertificateError) as caught:
            checker.verify_manifest(tampered)
        self.assertEqual(caught.exception.code, "PROMOTION_GUARD")

    def test_rehashed_open_gate_closure_is_rejected(self) -> None:
        tampered = json.loads(json.dumps(self.stored_manifest))
        tampered["open_gates"][1]["status"] = "closed"
        rehash(tampered)
        with self.assertRaises(cert.CertificateError) as caught:
            checker.verify_manifest(tampered)
        self.assertEqual(caught.exception.code, "OPEN_GATES")

    def test_rehashed_capacity_injection_is_rejected(self) -> None:
        tampered = json.loads(json.dumps(self.stored_manifest))
        tampered["source_firewall"]["allowed_issue_inputs"].append(505)
        tampered["source_firewall"]["issue_505_consumed"] = True
        rehash(tampered)
        with self.assertRaises(cert.CertificateError) as caught:
            checker.verify_manifest(tampered)
        self.assertEqual(caught.exception.code, "SOURCE_FIREWALL")

    def test_rehashed_weak_copy_mutation_is_rejected(self) -> None:
        tampered = json.loads(json.dumps(self.stored_manifest))
        tampered["electroweak_order_unit_line"]["copy_derivation"][
            "total_weak_doublet_copies"
        ] = 5
        rehash(tampered)
        with self.assertRaises(cert.CertificateError) as caught:
            checker.verify_manifest(tampered)
        self.assertEqual(caught.exception.code, "EW_LINE")

    def test_rehashed_upstream_pin_mutation_is_rejected(self) -> None:
        tampered = json.loads(json.dumps(self.stored_manifest))
        tampered["upstream_pins"]["carrier_manifest"]["sha256"] = "0" * 64
        rehash(tampered)
        with self.assertRaises(cert.CertificateError) as caught:
            checker.verify_manifest(tampered)
        self.assertEqual(caught.exception.code, "UPSTREAM_PINS")

    def test_no_floating_point_enters_the_manifest(self) -> None:
        cert.require_no_floats(self.stored_manifest)

    def test_cli_outputs_are_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            receipt = Path(scratch) / "receipt.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(first.read_bytes(), cert.MANIFEST_PATH.read_bytes())
            self.assertEqual(
                checker.main(
                    [
                        "--manifest",
                        str(first),
                        "--receipt",
                        str(receipt),
                        "--write",
                    ]
                ),
                0,
            )
            self.assertEqual(
                receipt.read_bytes(), checker.RECEIPT_PATH.read_bytes()
            )
            self.assertEqual(
                checker.main(
                    ["--manifest", str(first), "--receipt", str(receipt)]
                ),
                0,
            )


if __name__ == "__main__":
    unittest.main()
