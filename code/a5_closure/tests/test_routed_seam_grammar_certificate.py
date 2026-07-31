#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #613 grammar certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import routed_seam_grammar_certificate as cert  # noqa: E402


class RoutedSeamGrammarTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "routed_seam_grammar_reference.json"
        cls.receipt_path = MODULE_DIR / "receipts" / "routed_seam_grammar_reference.receipt.json"
        cls.negative_path = MODULE_DIR / "negative_controls" / "issue_613_negative_controls.json"
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    # -- determinism and reference receipts ---------------------------------

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_payload_is_deterministic(self) -> None:
        again = cert.certificate_payload(self.manifest)
        self.assertEqual(again, self.expected)
        self.assertEqual(cert.sha256_json(again), cert.sha256_json(self.expected))

    def test_reference_negative_controls_recompute(self) -> None:
        stored = cert.load_json(self.negative_path)
        recomputed = cert.negative_control_payload(self.manifest)
        self.assertEqual(stored, recomputed)
        self.assertTrue(all(row["passed"] for row in recomputed["finite_controls"]))

    # -- schema and verdict --------------------------------------------------

    def test_receipt_schema_and_issue(self) -> None:
        self.assertEqual(self.expected["schema"], cert.RECEIPT_SCHEMA)
        self.assertEqual(self.expected["issue"], 613)
        self.assertEqual(
            self.expected["manifest_sha256"], cert.sha256_json(self.manifest)
        )

    def test_verdict_fields(self) -> None:
        verdict = self.expected["verdict"]
        self.assertEqual(verdict["central_column_completeness"], "exact")
        self.assertEqual(
            verdict["central_class_group_smith_invariants"], [1, 1, 1, 1, 1, 6]
        )
        self.assertEqual(verdict["composition_laws"], "closed")
        self.assertEqual(verdict["general_grammar"], "conditional_open_interface")
        self.assertEqual(
            verdict["four_dimensional_instanton_theta"], "separate_gates"
        )
        witness = verdict["general_grammar_witness"]
        self.assertEqual(witness["transport_group_order"], 120)
        self.assertEqual(witness["centre_order"], 2)
        self.assertEqual(witness["noncommuting_lift_pairs"], 3)
        self.assertTrue(witness["artifact_sha256"].startswith("sha256:"))
        self.assertTrue(all(verdict["controls"].values()))
        self.assertEqual(
            self.expected["claim_boundary"]["bounded_exit"],
            "conditional_open_interface",
        )

    def test_consumed_upstream_pins(self) -> None:
        upstream = self.expected["consumed_upstream"]
        descent = cert.load_json(
            MODULE_DIR / "manifests" / "axis_center_descent_reference.json"
        )
        self.assertEqual(
            upstream["descent_manifest_sha256"], cert.sha256_json(descent)
        )
        self.assertEqual(
            upstream["global_form_artifact_sha256"],
            descent["global_form_artifact_sha256"],
        )
        self.assertEqual(
            upstream["spin_statistics_artifact_sha256"],
            self.manifest["spin_statistics_artifact_sha256"],
        )
        scope = upstream["physical_scope"]
        self.assertTrue(scope["declared_axis_coefficient_system_exact"])
        self.assertFalse(scope["axis_relation_lattice_source_selected"])
        self.assertFalse(scope["physical_global_form_selected"])
        self.assertFalse(scope["same_source_loop_to_tensor_kernel_identification"])

    # -- central-column completeness -----------------------------------------

    def test_central_column_completeness_is_exact(self) -> None:
        central = self.expected["central_column_grammar"]
        self.assertEqual(central["status"], "exact")
        self.assertEqual(
            central["axis_class_lattice"]["smith_invariants"], [1, 1, 1, 1, 1, 6]
        )
        self.assertEqual(central["boundary_smith_invariants"], [1] * 19)
        self.assertEqual(central["realized_class_menu"], [0, 1, 2, 3, 4, 5])
        self.assertEqual(central["generator_tubes_constructed"], 19)
        self.assertEqual(central["pair_realization_checks"], 19 * 18 * 5)
        self.assertTrue(central["puncture_faces_antipodal"])
        self.assertTrue(central["subgroup_obstruction_menu_recomputed"])
        self.assertEqual(len(central["pinned_witnesses_recomputed"]), 6)

    def test_composition_laws_close(self) -> None:
        laws = self.expected["composition_laws"]
        self.assertEqual(laws["status"], "closed")
        self.assertEqual(laws["pairwise_composition_checks"], 36)
        self.assertEqual(laws["seam_scaling_checks"], 180)
        self.assertEqual(laws["triple_overlap_coherence_entries"], 240)
        self.assertEqual(laws["face_loop_holonomy_checks"], 120)
        self.assertEqual(laws["homotopy_move_checks"], 120)
        self.assertEqual(laws["loop_associativity_checks"], 8000)
        self.assertTrue(laws["coboundary_lands_in_class_group"])

    def test_frame_and_conjugate_invariance(self) -> None:
        frame = self.expected["frame_reconstruction"]
        self.assertTrue(frame["antipode_is_exact_negation"])
        self.assertEqual(frame["seams_from_squared_distance_four"], 30)
        self.assertEqual(frame["oriented_face_determinants_positive"], 20)
        self.assertTrue(frame["conjugate_frame_same_seams"])
        self.assertTrue(frame["class_menu_invariant_under_conjugation"])

    # -- nonabelian scope boundary -------------------------------------------

    def test_nonabelian_boundary_witness(self) -> None:
        boundary = self.expected["nonabelian_scope_boundary"]
        self.assertEqual(boundary["transport_group_order"], 120)
        self.assertEqual(boundary["centre_order"], 2)
        self.assertEqual(boundary["unique_nontrivial_involution"], "-1")
        self.assertEqual(boundary["klein_four_lift_squares_minus_one"], 3)
        self.assertEqual(boundary["noncommuting_lift_pairs"], 3)
        self.assertTrue(boundary["double_cover_non_split"])

    def test_quaternion_arithmetic_is_exact(self) -> None:
        i = (Fraction(0), Fraction(1), Fraction(0), Fraction(0))
        j = (Fraction(0), Fraction(0), Fraction(1), Fraction(0))
        k = (Fraction(0), Fraction(0), Fraction(0), Fraction(1))
        self.assertEqual(cert.quaternion_multiply(i, j), k)
        self.assertNotEqual(
            cert.quaternion_multiply(i, j), cert.quaternion_multiply(j, i)
        )
        self.assertEqual(
            cert.quaternion_multiply(i, i), cert.QUATERNION_MINUS_ONE
        )

    # -- structural fail-closed controls -------------------------------------

    def test_structural_controls_fail_closed(self) -> None:
        controls = {row["name"]: row for row in self.expected["structural_controls"]}
        self.assertEqual(
            set(controls),
            {
                "wrong_routing_one_seam_tampered",
                "missing_triple_overlap_deleted_face",
                "conjugate_frame_extra_sector_claim",
                "extra_flux_seventh_sector_claim",
            },
        )
        self.assertEqual(
            controls["wrong_routing_one_seam_tampered"]["actual_error"], "SEAM_ROUTING"
        )
        self.assertEqual(
            controls["missing_triple_overlap_deleted_face"]["actual_error"],
            "NERVE_INCIDENCE",
        )
        self.assertEqual(
            controls["conjugate_frame_extra_sector_claim"]["actual_error"],
            "CONJUGATE_FRAME",
        )
        self.assertEqual(
            controls["extra_flux_seventh_sector_claim"]["actual_error"], "EXTRA_FLUX"
        )

    def test_extra_flux_rejection_is_direct(self) -> None:
        self.assertEqual(cert.admit_sector_classes([0, 1, 2, 3, 4, 5], 6), 6)
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.admit_sector_classes([0, 1, 2, 3, 4, 5, 6], 6)
        self.assertEqual(ctx.exception.code, "EXTRA_FLUX")

    def test_frame_dependent_sector_rejection_is_direct(self) -> None:
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.reject_frame_dependent_sector(9, [0, 1, 2, 3, 4, 5], 6)
        self.assertEqual(ctx.exception.code, "CONJUGATE_FRAME")

    # -- mutation gates -------------------------------------------------------

    def test_descent_manifest_pin_tamper_is_detected(self) -> None:
        tampered = copy.deepcopy(dict(self.manifest))
        tampered["descent_manifest_sha256"] = "0" * 64
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.certificate_payload(tampered)
        self.assertEqual(ctx.exception.code, "UPSTREAM_HASH")

    def test_spin_artifact_pin_tamper_is_detected(self) -> None:
        tampered = copy.deepcopy(dict(self.manifest))
        tampered["spin_statistics_artifact_sha256"] = "sha256:" + "0" * 64
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.certificate_payload(tampered)
        self.assertEqual(ctx.exception.code, "UPSTREAM_HASH")

    def test_carrier_swap_is_detected(self) -> None:
        tampered = copy.deepcopy(dict(self.manifest))
        tampered["carrier_manifest_path"] = tampered["spin_statistics_artifact_path"]
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.certificate_payload(tampered)
        self.assertEqual(ctx.exception.code, "CARRIER_BINDING")

    def test_forbidden_promotion_keys_are_rejected(self) -> None:
        for key in cert.FORBIDDEN_MANIFEST_KEYS:
            tampered = copy.deepcopy(dict(self.manifest))
            tampered[key] = True
            with self.assertRaises(cert.CertificateError) as ctx:
                cert.certificate_payload(tampered)
            self.assertEqual(ctx.exception.code, "FORBIDDEN_DEPENDENCY")

    def test_undeclared_manifest_key_is_rejected(self) -> None:
        tampered = copy.deepcopy(dict(self.manifest))
        tampered["extra_input"] = 1
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.certificate_payload(tampered)
        self.assertEqual(ctx.exception.code, "SCHEMA_FIELDS")

    def test_receipt_tamper_is_detected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["verdict"]["general_grammar"] = "exact_named_realization"
        with self.assertRaises(cert.CertificateError) as ctx:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(ctx.exception.code, "RECEIPT_MISMATCH")

    # -- command line ---------------------------------------------------------

    def test_cli_verify_exits_zero(self) -> None:
        exit_code = cert.main(
            [
                "verify",
                "--manifest",
                str(self.manifest_path),
                "--receipt",
                str(self.receipt_path),
            ]
        )
        self.assertEqual(exit_code, 0)


if __name__ == "__main__":
    unittest.main()
