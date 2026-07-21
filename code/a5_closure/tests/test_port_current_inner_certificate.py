#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #566 certificate."""

from __future__ import annotations

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402
import port_current_inner_certificate as cert  # noqa: E402


class PortCurrentInnerCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "port_current_response_reference.json"
        cls.receipt_path = MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json"
        cls.negative_path = MODULE_DIR / "negative_controls" / "issue_566_negative_controls.json"
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_gate_passes_with_all_rows_true(self) -> None:
        gate = self.expected["physical_current_gate"]
        self.assertTrue(gate["passed"])
        self.assertTrue(all(gate.values()))

    def test_port_to_generator_map_shape(self) -> None:
        row = self.expected["port_to_generator_map"]
        self.assertTrue(row["injective"])
        self.assertEqual(row["image_real_dimension"], 12)
        self.assertTrue(row["skew_adjoint"])
        self.assertEqual(row["block_dimensions"], {"even_block_u3": 9, "kernel_block_so3": 3})
        self.assertTrue(row["block_dimensions_verified"])
        self.assertTrue(row["kernel_block_real_verified"])
        self.assertEqual(row["compact_lie_type"], "u(3) (+) so(3) = u(1) (+) su(3) (+) su(2)")

    def test_closure_center_and_adjoint_rank(self) -> None:
        row = self.expected["closure"]
        self.assertTrue(row["commutator_closed"])
        self.assertEqual(row["derived_dimension"], 11)
        self.assertEqual(
            row["derived_block_dimensions"],
            {"even_block_su3": 8, "kernel_block_so3": 3},
        )
        self.assertEqual(row["center_dimension"], 1)
        self.assertTrue(row["center_is_constant_even_port_line"])
        self.assertEqual(row["adjoint_rank"], 11)
        self.assertEqual(len(row["structure_constants"]), 66)

    def test_derivation_chain_and_acceptance_status(self) -> None:
        chain = self.expected["derivation_chain"]
        self.assertEqual(len(chain), 13)
        self.assertEqual([row["step"] for row in chain], list(range(1, 14)))
        for row in chain:
            self.assertTrue(row["premise"])
            self.assertTrue(row["uses"])
            self.assertTrue(row["source_artifact"])
            self.assertTrue(row["conclusion"])
        status = self.expected["acceptance_criteria_status"]
        self.assertEqual(len(status), 5)
        self.assertEqual(set(status.values()), {True})
        self.assertIn("branch_scope", self.expected)
        self.assertIn("factor_origins", self.expected)
        self.assertIn("dependency_acyclicity_note", self.expected)
        self.assertIn("verify", self.expected["verifier_command"])

    def test_response_provenance_is_typed_not_asserted(self) -> None:
        definedness = self.expected["source_definedness"]
        self.assertEqual(definedness["response_packet_provenance"], "declared_source_data")
        self.assertFalse(definedness["measurement_artifact_pinned"])
        self.assertTrue(definedness["all_source_defined"])
        closure_condition = self.expected["issue_closure_condition"]
        self.assertIs(closure_condition["met_locally"], True)
        self.assertIn("open source data", closure_condition["source_data"])
        self.assertIn("do not require it", closure_condition["optional_strengthening"])

    def test_pinned_measurement_artifact_upgrades_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            artifact_path = Path(tmp) / "response_measurement.json"
            artifact = {"kind": "response_measurement_stub_for_contract_test"}
            cert.write_json(artifact_path, artifact)
            manifest = copy.deepcopy(self.manifest)
            manifest["response_measurement_contract"]["measurement_status"] = "measured"
            manifest["response_measurement_contract"]["measurement_artifact"] = {
                "path": str(artifact_path),
                "sha256": cert.sha256_json(artifact),
            }
            receipt = cert.certificate_payload(manifest)
        self.assertIs(receipt["issue_closure_condition"]["met_locally"], True)
        self.assertEqual(receipt["issue_closure_condition"]["response_field_provenance"], "measured")
        self.assertIsNone(receipt["issue_closure_condition"]["optional_strengthening"])
        self.assertTrue(receipt["source_definedness"]["measurement_artifact_pinned"])
        wrong_hash = copy.deepcopy(self.manifest)
        wrong_hash["response_measurement_contract"]["measurement_status"] = "measured"
        wrong_hash["response_measurement_contract"]["measurement_artifact"] = {
            "path": "manifests/port_current_response_reference.json",
            "sha256": "0" * 64,
        }
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(wrong_hash)
        self.assertEqual(caught.exception.code, "RESPONSE_ARTIFACT")

    def test_register_relabeling_no_go(self) -> None:
        row = self.expected["response_versus_register_relabeling"]
        self.assertEqual(row["element_orders"], [1, 2, 3, 5])
        self.assertEqual(row["order_five_elements_with_irrational_sector_characters"], 24)
        self.assertEqual(row["sector_character_norms"], {"even_block": "1", "kernel_block": "1"})
        self.assertIn("cannot generate these currents", row["conclusion"])

    def test_hilbert_schmidt_band_coefficients_are_galois_paired(self) -> None:
        bands = self.expected["compactness"]["hilbert_schmidt_pullback_band_coefficients"]
        self.assertEqual(bands["unit_band"], "1/4")
        self.assertEqual(bands["quintet_band"], "3 + 1*sqrt(5)")
        self.assertEqual(bands["frame_band"], "5 + 1*sqrt(5)")
        self.assertEqual(bands["kernel_band"], "5 + -1*sqrt(5)")
        self.assertTrue(self.expected["compactness"]["positive_definite"])
        self.assertTrue(self.expected["compactness"]["a5_invariant"])

    def test_covariance_homomorphism_and_innerness_witness_counts(self) -> None:
        intertwiner = self.expected["icosahedral_intertwiner"]
        self.assertEqual(intertwiner["covariance_checks"], 60 * 12)
        self.assertEqual(intertwiner["implementer_homomorphism_pairs"], 60 * 60)
        self.assertTrue(intertwiner["implementers_faithful"])
        inner = self.expected["inner_action"]
        self.assertTrue(inner["block_skew_pairs_in_image"])
        self.assertEqual(inner["witness_count"], 60)
        orders = sorted({row["element_order"] for row in inner["witnesses"]})
        self.assertEqual(orders, [1, 2, 3, 5])

    def test_order_five_witness_cosines_are_galois_conjugate(self) -> None:
        witnesses = [row for row in self.expected["inner_action"]["witnesses"] if row["element_order"] == 5]
        self.assertEqual(len(witnesses), 24)
        pairs = {(row["even_block"]["cosine"], row["kernel_block"]["cosine"]) for row in witnesses}
        self.assertEqual(
            pairs,
            {
                ("-1/4 + -1/4*sqrt(5)", "-1/4 + 1/4*sqrt(5)"),
                ("-1/4 + 1/4*sqrt(5)", "-1/4 + -1/4*sqrt(5)"),
            },
        )

    def test_schur_rigidity_and_response_moduli(self) -> None:
        rigidity = self.expected["icosahedral_intertwiner"]["kernel_band_schur_rigidity"]
        self.assertEqual(rigidity["multiplicity_of_kernel_band_in_even_block"], "0")
        moduli = self.expected["response_moduli"]
        self.assertEqual(moduli["equivariant_lift_dimension"], 4)
        self.assertEqual(len(moduli["open_source_data"]), 5)

    def test_refinement_naturality(self) -> None:
        row = self.expected["refinement"]
        self.assertTrue(row["natural"])
        self.assertEqual(len(row["naturality"]), 3)
        self.assertTrue(all(item["intertwined"] for item in row["naturality"]))
        self.assertTrue(row["carrier_tower"]["all_maps_in_A5"])

    def test_repair_response_distinction(self) -> None:
        row = self.expected["repair_response_distinction"]
        self.assertTrue(row["reversible_response_automorphisms_define_currents"])
        self.assertTrue(row["strict_descent_repairs_typed_irreversible"])
        self.assertTrue(row["disjoint"])

    def test_classification_realization_distinction_recorded(self) -> None:
        row = self.expected["classification_vs_realization"]
        self.assertTrue(row["distinguished"])
        self.assertIn("commute", row["coefficient_layer"])
        self.assertIn("CENTER_NOT_ONE_DIMENSIONAL", row["separating_witness"])

    def test_all_finite_negative_controls_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        stored = cert.load_json(self.negative_path)
        self.assertEqual(payload, stored)
        self.assertGreaterEqual(len(payload["finite_controls"]), 10)
        self.assertTrue(all(row["passed"] for row in payload["finite_controls"]))
        by_name = {row["name"]: row["actual_error"] for row in payload["finite_controls"]}
        self.assertEqual(by_name["abelian_record_model"], "CENTER_NOT_ONE_DIMENSIONAL")
        self.assertEqual(by_name["rank_deficient_kernel_band"], "IMAGE_RANK_DEFICIENT")
        self.assertEqual(by_name["register_relabeling_conflated_as_response"], "REGISTER_RELABELING_CONFLATION")
        self.assertEqual(by_name["phantom_measurement_claim"], "RESPONSE_ARTIFACT")

    def test_forbidden_downstream_data_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.manifest)
        mutant["target_hint"] = {"standard_model_current": "adjoint"}
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "FORBIDDEN_DEPENDENCY")

    def test_carrier_hash_pin_is_enforced(self) -> None:
        mutant = copy.deepcopy(self.manifest)
        mutant["carrier_manifest_sha256"] = "0" * 64
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "CARRIER_HASH")

    def test_tampered_receipt_is_rejected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["closure"]["derived_dimension"] = 12
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(caught.exception.code, "RECEIPT_MISMATCH")

    def test_arbitrary_port_relabeling_preserves_invariants(self) -> None:
        carrier_manifest = cert.load_json(MODULE_DIR / "manifests" / "echosahedral_federation_reference.json")
        sigma = (7, 2, 10, 4, 0, 11, 5, 9, 1, 6, 3, 8)
        relabeled = e565.relabel_manifest(carrier_manifest, sigma)
        with tempfile.TemporaryDirectory() as tmp:
            relabeled_path = Path(tmp) / "relabeled_carrier.json"
            cert.write_json(relabeled_path, relabeled)
            relabeled_carrier = cert.load_json(relabeled_path)
            manifest = copy.deepcopy(self.manifest)
            manifest["carrier_manifest_path"] = str(relabeled_path)
            manifest["carrier_manifest_sha256"] = cert.sha256_json(relabeled_carrier)
            receipt = cert.certificate_payload(manifest)
        self.assertEqual(
            receipt["compactness"]["hilbert_schmidt_pullback_band_coefficients"],
            self.expected["compactness"]["hilbert_schmidt_pullback_band_coefficients"],
        )
        self.assertEqual(receipt["closure"]["derived_dimension"], 11)
        self.assertEqual(receipt["closure"]["center_dimension"], 1)
        self.assertTrue(receipt["physical_current_gate"]["passed"])

    def test_galois_field_arithmetic(self) -> None:
        phi = cert.PHI
        self.assertEqual(phi * phi, phi + cert.ONE)
        norm = cert.VERTEX_NORM
        sqrt5 = cert.F5(0, 1)
        self.assertEqual(norm, phi * sqrt5)
        self.assertEqual(phi.conj() * phi, -cert.ONE)
        self.assertTrue((phi.inv() * phi - cert.ONE).is_zero())

    def test_cli_certify_verify_and_negative_controls(self) -> None:
        script = MODULE_DIR / "port_current_inner_certificate.py"
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            receipt = tmpdir / "receipt.json"
            controls = tmpdir / "negative.json"
            subprocess.run(
                [sys.executable, str(script), "certify", "--manifest", str(self.manifest_path), "--output", str(receipt)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sys.executable, str(script), "verify", "--manifest", str(self.manifest_path), "--receipt", str(receipt)],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sys.executable, str(script), "negative-controls", "--manifest", str(self.manifest_path), "--output", str(controls)],
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertEqual(cert.load_json(receipt), self.expected)
            self.assertEqual(cert.load_json(controls), cert.negative_control_payload(self.manifest))


if __name__ == "__main__":
    unittest.main()
