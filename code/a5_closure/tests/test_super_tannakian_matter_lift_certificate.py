#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #314 matter-lift certificate."""

from __future__ import annotations

import copy
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import port_current_inner_certificate as p566  # noqa: E402
import super_tannakian_matter_lift_certificate as cert  # noqa: E402


class SuperTannakianMatterLiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "super_tannakian_matter_reference.json"
        cls.receipt_path = MODULE_DIR / "receipts" / "super_tannakian_matter_reference.receipt.json"
        cls.negative_path = MODULE_DIR / "negative_controls" / "issue_314_negative_controls.json"
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_conditional_gate_passes_and_physical_gate_is_open(self) -> None:
        gate = self.expected["conditional_algebraic_gate"]
        self.assertEqual(set(gate.values()), {True})
        self.assertEqual(len(gate), 13)
        self.assertTrue(gate["passed"])
        physical_gate = self.expected["physical_source_gate"]
        self.assertEqual(set(physical_gate.values()), {False})
        self.assertFalse(physical_gate["passed"])

    def test_upstream_pins_match_stored_packets(self) -> None:
        upstream = self.expected["upstream"]
        current_manifest = cert.load_json(MODULE_DIR / "manifests" / "port_current_response_reference.json")
        current_receipt = cert.load_json(MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json")
        self.assertEqual(upstream["current_manifest_sha256"], cert.sha256_json(current_manifest))
        self.assertEqual(upstream["current_receipt_sha256"], cert.sha256_json(current_receipt))
        self.assertEqual(upstream["carrier_manifest_sha256"], current_receipt["carrier_manifest_sha256"])

    def test_port_spin_lift_is_nonsplit_binary_icosahedral(self) -> None:
        row = self.expected["port_spin_lift"]
        self.assertEqual(row["witness_count"], 60)
        self.assertEqual(row["lift_group_order"], 120)
        self.assertTrue(row["unique_involution"])
        self.assertEqual(
            row["order_profile"],
            {"1": 1, "2": 1, "3": 20, "4": 30, "5": 24, "6": 20, "10": 24},
        )
        self.assertEqual(row["involution_lift_order"], 4)
        self.assertEqual(row["irrational_order_five_spinor_traces"], 24)
        self.assertIn("non-split", row["conclusion"])

    def test_exact_spinor_lift_solver(self) -> None:
        # The half-angle data of an order-five icosahedral rotation lives in
        # Q(sqrt5): sqrt((3+sqrt5)/8) = (1+sqrt5)/4.
        value = cert.F5(Fraction(3, 8), Fraction(1, 8))
        root = cert.sqrt_f5(value)
        self.assertIsNotNone(root)
        self.assertEqual(root, cert.F5(Fraction(1, 4), Fraction(1, 4)))
        self.assertIsNone(cert.sqrt_f5(cert.F5(2)))
        self.assertEqual(cert.sqrt_f5(cert.F5(Fraction(9, 4))), cert.F5(Fraction(3, 2)))
        self.assertEqual(cert.sqrt_f5(cert.F5(5)), cert.F5(0, 1))

    def test_matter_transport_is_faithful_homomorphism(self) -> None:
        row = self.expected["matter_carrier"]
        self.assertEqual(row["faithful_rank_on_carrier"], 12)
        self.assertEqual(row["transport_homomorphism_bracket_checks"], 66)
        self.assertEqual(row["transport_covariance_checks"], 720)
        self.assertIn("BLOCK-DETERMINANT-BALANCE", row["central_charge_provenance"])

    def test_car_fock_structure(self) -> None:
        row = self.expected["auxiliary_car_fock"]
        self.assertEqual(row["dimension"], 32)
        self.assertEqual(row["car_relation_checks"], 50)
        self.assertEqual(row["vacuum_cyclic_rank"], 32)
        self.assertEqual(row["second_quantization_derivation_checks"], 60)
        self.assertIn("(-1)^N", row["fermionic_parity"])

    def test_selection_projector_rank_fifteen(self) -> None:
        row = self.expected["selection"]
        self.assertEqual(row["rule"], "parity_even_minus_derived_invariants")
        self.assertEqual(row["projector_rank"], 15)
        self.assertEqual(row["equivariance_checks"], 12)
        self.assertTrue(row["commutes_with_parity"])
        self.assertIn("not representation arithmetic", row["realization"])

    def test_realized_package_spectrum_and_fields(self) -> None:
        row = self.expected["realized_package"]
        self.assertEqual(row["dimension"], 15)
        self.assertEqual(row["faithful_rank_on_matter"], 12)
        self.assertEqual(row["integrality_normalization"], 6)
        self.assertEqual(
            row["charge_spectrum"],
            {"1/6": 6, "-2/3": 3, "1": 1, "1/3": 3, "-1/2": 2},
        )
        self.assertTrue(row["multiplicity_free"])
        self.assertTrue(row["contains_no_invariant_line"])
        self.assertEqual(row["irreducible_block_commutants"], 5)
        dimensions = {name: data["dimension"] for name, data in row["fields"].items()}
        self.assertEqual(dimensions, {"Q": 6, "u_c": 3, "e_c": 1, "d_c": 3, "L": 2})
        for data in row["fields"].values():
            self.assertEqual(data["commutant_dimension"], 1)

    def test_chirality_and_conjugation(self) -> None:
        chirality = self.expected["chirality"]
        self.assertTrue(chirality["matter_spectrum_disjoint_from_dual"])
        self.assertEqual(chirality["hom_dimension_with_dual"], 0)
        conjugation = self.expected["conjugation"]
        self.assertEqual(conjugation["invariance_checks"], 12)
        self.assertTrue(conjugation["nondegenerate"])
        self.assertIn("top line", conjugation["pairing"])

    def test_realized_anomalies_and_witten_parity(self) -> None:
        anomalies = self.expected["anomalies"]
        self.assertEqual(
            anomalies["traces"],
            {
                "gravity_squared_U1": "0",
                "U1_cubed": "0",
                "SU3_squared_U1": "0",
                "SU2_squared_U1": "0",
                "SU3_cubed": "0",
            },
        )
        self.assertEqual(anomalies["su3_d_symbol_checks"], 120)
        self.assertEqual(anomalies["witten_parity"]["weak_doublets"], 4)
        self.assertTrue(anomalies["witten_parity"]["even"])

    def test_yukawa_invariant_sector(self) -> None:
        row = self.expected["yukawa_sector"]
        self.assertEqual(row["invariant_sector_dimension"], 3)
        self.assertEqual(
            [item["channel"] for item in row["channels"]],
            [["Q", "S", "u_c"], ["Q", "Sbar", "d_c"], ["L", "Sbar", "e_c"]],
        )
        self.assertTrue(all(item["invariant_dimension"] == 1 for item in row["channels"]))
        self.assertEqual(row["forbidden_channel_control"]["invariant_dimension"], 0)

    def test_kernel_is_emitted_on_the_cover_not_assumed(self) -> None:
        row = self.expected["kernel_emission"]
        # On the simply connected cover R x SU(3) x SU(2) the kernel is
        # infinite: the unit deck translation (one full central turn) acts
        # trivially on every integral weight without being the identity.
        self.assertIn("infinite cyclic", row["kernel_group_on_cover"])
        self.assertIn("non-compact R", row["central_factor"])
        self.assertEqual(row["integrality_normalization"], 6)
        self.assertEqual(
            row["kernel_generator"],
            {"u1_phase_turns": "1/6", "su3_center_power": 1, "su2_center_power": 1},
        )
        self.assertIn("generator^6", row["deck_relation"])
        self.assertIn("not the identity on the cover", row["deck_relation"])
        self.assertEqual(row["residual_order_modulo_deck_translations"], 6)
        self.assertEqual(len(row["kernel_residues_modulo_deck_translations"]), 6)
        self.assertFalse(row["global_quotient_assumed"])
        self.assertIn("AXIS-CENTER-DESCENT", row["downstream_consumer"])

    def test_kernel_generator_sixth_power_is_deck_translation_not_identity(self) -> None:
        # Exact recomputation of the inconsistency fix: compose the generator
        # six times in R x Z3 x Z2 WITHOUT reducing the R coordinate mod 1.
        row = self.expected["kernel_emission"]
        gen = row["kernel_generator"]
        r = Fraction(gen["u1_phase_turns"])
        a = gen["su3_center_power"]
        b = gen["su2_center_power"]
        sixth = (6 * r, (6 * a) % 3, (6 * b) % 2)
        self.assertEqual(sixth, (Fraction(1), 0, 0))
        self.assertNotEqual(sixth, (Fraction(0), 0, 0))

    def test_refinement_descent(self) -> None:
        row = self.expected["refinement"]
        self.assertTrue(row["natural"])
        self.assertEqual(len(row["maps"]), 3)
        self.assertTrue(all(item["intertwined"] for item in row["maps"]))

    def test_mar_class_nonempty_without_uniqueness(self) -> None:
        row = self.expected["mar_class"]
        self.assertTrue(row["nonempty"])
        self.assertFalse(row["uniqueness_promoted"])
        self.assertIn("nonemptiness precondition", row["note"])

    def test_acceptance_criteria_chain_and_scope(self) -> None:
        status = self.expected["acceptance_criteria_status"]
        self.assertEqual(len(status), 9)
        # The source-derivation criterion is honestly False: the structures
        # are derived exactly, but from declared branch premises that are not
        # yet source-bound (mirroring the #566 conditional scoping).
        self.assertIs(
            status["fermionic_parity_spin_lift_chirality_conjugation_tensor_product_source_derived"],
            False,
        )
        algebraic_rows = {k: v for k, v in status.items() if not k.endswith("source_derived")}
        self.assertEqual(set(algebraic_rows.values()), {True})
        chain = self.expected["derivation_chain"]
        self.assertEqual([row["step"] for row in chain], list(range(1, 17)))
        for row in chain:
            self.assertTrue(row["premise"])
            self.assertTrue(row["uses"])
            self.assertTrue(row["source_artifact"])
            self.assertTrue(row["conclusion"])
        self.assertIn("branch_scope", self.expected)
        self.assertIn("factor_origins", self.expected)
        self.assertIn("dependency_acyclicity_note", self.expected)
        self.assertIn("verify", self.expected["verifier_command"])

    def test_matter_lift_is_conditional_not_source_bound(self) -> None:
        closure_condition = self.expected["issue_closure_condition"]
        self.assertTrue(closure_condition["conditional_algebraic_gate_passed"])
        self.assertFalse(closure_condition["physical_source_realization_gate_passed"])
        self.assertFalse(closure_condition["met_locally"])
        self.assertIn("source binding", closure_condition["remaining_producer"])
        boundary = self.expected["claim_boundary"]
        self.assertEqual(boundary["status"], "proved_conditional_on_declared_matter_contracts")
        self.assertIn("PORT-SPIN-LIFT as a physical source-bound receipt", boundary["does_not_close"][0])
        self.assertTrue(any("BLOCK-DETERMINANT-BALANCE" in row for row in boundary["does_not_close"]))
        self.assertTrue(any("AXIS-CENTER-DESCENT" in row for row in boundary["does_not_close"]))

    def test_control_contracts_are_not_valid_production_manifests(self) -> None:
        for path, value, code in (
            (("category_contract", "typing"), "vec", "CATEGORY_TYPING"),
            (("category_contract", "typing"), "svec", "CATEGORY_TYPING"),
            (("category_contract", "selection_rule"), "lambda2_only", "SELECTION_RULE"),
            (("statistics_contract", "matter_statistics"), "bosonic_even", "STATISTICS_TYPING"),
        ):
            mutant = copy.deepcopy(self.manifest)
            mutant[path[0]][path[1]] = value
            with self.assertRaises(cert.CertificateError) as caught:
                cert.validate_manifest(mutant)
            self.assertEqual(caught.exception.code, code)
        split = copy.deepcopy(self.manifest)
        split["category_contract"]["spin_lift"]["double_cover"] = False
        with self.assertRaises(cert.CertificateError) as caught:
            cert.validate_manifest(split)
        self.assertEqual(caught.exception.code, "CATEGORY_TYPING")

    def test_redundant_trace_balance_flag_is_not_required(self) -> None:
        # Trace balance is an arithmetic consequence of the declared charge
        # pair; the manifest carries no redundant declared balance flag.
        self.assertNotIn("trace_balanced", self.manifest["exterior_matter_contract"])
        unbalanced = copy.deepcopy(self.manifest)
        unbalanced["exterior_matter_contract"]["block_trace_charges"]["weak_block"] = "1/3"
        with self.assertRaises(cert.CertificateError) as caught:
            cert.validate_manifest(unbalanced)
        self.assertEqual(caught.exception.code, "TRACE_BALANCE")

    def test_all_finite_negative_controls_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        stored = cert.load_json(self.negative_path)
        self.assertEqual(payload, stored)
        self.assertEqual(len(payload["finite_controls"]), 15)
        self.assertTrue(all(row["passed"] for row in payload["finite_controls"]))
        by_name = {row["name"]: row["actual_error"] for row in payload["finite_controls"]}
        self.assertEqual(by_name["vec_typing"], "VEC_TYPING")
        self.assertEqual(by_name["svec_split_spin"], "SPIN_LIFT_SPLIT")
        self.assertEqual(by_name["opposite_weyl_selection"], "YUKAWA_CHANNEL_EMPTY")
        self.assertEqual(by_name["truncated_lambda2_selection"], "WITTEN_PARITY")
        self.assertEqual(by_name["full_even_clifford_module"], "TRIVIAL_LINE_IN_MATTER")
        self.assertEqual(by_name["kernel_killing_extra_scalar"], "KERNEL_TRIVIAL")
        self.assertEqual(by_name["charge_dead_package"], "CURRENT_ACTION_NOT_FAITHFUL")

    def test_forbidden_matter_targets_are_rejected(self) -> None:
        for hint in (
            {"attachment": "family attachment rank"},
            {"scalar_sector": "scalar potential"},
            {"mass": "pole mass input"},
        ):
            mutant = copy.deepcopy(self.manifest)
            mutant["downstream_hint"] = hint
            with self.assertRaises(cert.CertificateError) as caught:
                cert.certificate_payload(mutant)
            self.assertEqual(caught.exception.code, "FORBIDDEN_DEPENDENCY")

    def test_upstream_hash_pins_are_enforced(self) -> None:
        mutant = copy.deepcopy(self.manifest)
        mutant["current_manifest_sha256"] = "0" * 64
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "UPSTREAM_HASH")
        mutant = copy.deepcopy(self.manifest)
        mutant["current_receipt_sha256"] = "0" * 64
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "UPSTREAM_HASH")

    def test_tampered_upstream_receipt_is_rejected(self) -> None:
        receipt = cert.load_json(MODULE_DIR / "receipts" / "port_current_inner_reference.receipt.json")
        tampered = copy.deepcopy(receipt)
        tampered["conditional_algebraic_gate"]["passed"] = False
        with tempfile.TemporaryDirectory() as tmp:
            tampered_path = Path(tmp) / "tampered_receipt.json"
            cert.write_json(tampered_path, tampered)
            mutant = copy.deepcopy(self.manifest)
            mutant["current_receipt_path"] = str(tampered_path)
            mutant["current_receipt_sha256"] = cert.sha256_json(tampered)
            with self.assertRaises(cert.CertificateError) as caught:
                cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "UPSTREAM_RECEIPT")

    def test_tampered_receipt_is_rejected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["kernel_emission"]["residual_order_modulo_deck_translations"] = 12
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(caught.exception.code, "RECEIPT_MISMATCH")

    def test_selection_rule_typo_is_rejected(self) -> None:
        mutant = copy.deepcopy(self.manifest)
        mutant["category_contract"]["selection_rule"] = "hand_picked_states"
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "SELECTION_RULE")

    def test_spin_lift_solver_matches_all_kernel_rotations(self) -> None:
        upstream = cert.load_upstream(self.manifest, MODULE_DIR)
        algebra = cert.CurrentAlgebra(upstream["current_manifest"], MODULE_DIR)
        identity = cert.cidentity(2)
        for g in algebra.plus[:10]:
            lift = cert.spin_lift_of_rotation(algebra.kernel_rotations[g])
            self.assertTrue(cert.spin_lift_matches(lift, algebra.kernel_rotations[g]))
            product = cert.cmul(cert.cdagger(lift), lift)
            self.assertTrue(cert.c_is_zero(cert.csub(product, identity)))

    def test_cli_certify_and_verify(self) -> None:
        script = MODULE_DIR / "super_tannakian_matter_lift_certificate.py"
        with tempfile.TemporaryDirectory() as tmp:
            receipt = Path(tmp) / "receipt.json"
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
            self.assertEqual(cert.load_json(receipt), self.expected)

    def test_issue_566_receipt_schema_is_pinned(self) -> None:
        # The upstream gate check relies on the #566 schema constant; make the
        # cross-module coupling explicit so schema drift fails loudly here.
        self.assertEqual(p566.RECEIPT_SCHEMA, "oph.port_current_inner_receipt.v2")


if __name__ == "__main__":
    unittest.main()
