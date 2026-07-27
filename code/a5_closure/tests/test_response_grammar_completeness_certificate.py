#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #611 response-grammar
completeness certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import response_grammar_completeness_certificate as cert  # noqa: E402


class ResponseGrammarCompletenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = (
            MODULE_DIR / "manifests" / "response_grammar_completeness_reference.json"
        )
        cls.receipt_path = (
            MODULE_DIR
            / "receipts"
            / "response_grammar_completeness_reference.receipt.json"
        )
        cls.negative_path = (
            MODULE_DIR / "negative_controls" / "issue_611_negative_controls.json"
        )
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_commutant_dimension_carries_the_lean_reference(self) -> None:
        commutant = self.expected["commutant"]
        self.assertEqual(commutant["dimension"], 4)
        self.assertEqual(commutant["pair_orbit_count"], 4)
        self.assertEqual(commutant["lean_reference"], "Lean/Screen/A5Commutant.lean")
        self.assertIn("commutant_decomposition", commutant["lean_theorems"])
        self.assertIn("orbitals_independent", commutant["lean_theorems"])
        self.assertIn("equivariant_iff_commutes", commutant["lean_theorems"])
        self.assertEqual(
            set(commutant["pair_orbit_representatives"]),
            {"diagonal", "adjacent", "distance_two", "antipodal"},
        )
        self.assertEqual(
            commutant["orbital_in_projector_basis"]["antipode"],
            {"1": "1", "3": "-1", "3p": "-1", "5": "1"},
        )

    def test_projector_system_is_exact(self) -> None:
        projectors = self.expected["isotypic_projectors"]
        self.assertTrue(projectors["idempotent"])
        self.assertTrue(projectors["mutually_orthogonal"])
        self.assertTrue(projectors["partition_of_unity"])
        self.assertEqual(
            projectors["ranks_by_exact_trace"], {"1": 1, "3": 3, "3p": 3, "5": 5}
        )
        self.assertEqual(
            projectors["diagonal_entries"],
            {"1": "1/12", "3": "1/4", "3p": "1/4", "5": "5/12"},
        )
        self.assertEqual(
            self.expected["port_space"]["adjacency_eigenvalues"],
            {"1": "5", "3": "1*sqrt(5)", "3p": "-1*sqrt(5)", "5": "-1"},
        )

    def test_sixteen_patterns_and_fourteen_admissible(self) -> None:
        completeness = self.expected["involution_completeness"]
        self.assertEqual(completeness["pattern_count"], 16)
        self.assertEqual(completeness["trivial_patterns"], ["++++", "----"])
        self.assertEqual(completeness["admissible_count"], 14)
        self.assertEqual(len(completeness["admissible_patterns"]), 14)
        self.assertEqual(completeness["sector_order"], ["1", "3", "3p", "5"])

    def test_exactly_two_signed_automorphism_realizations(self) -> None:
        selection = self.expected["signed_automorphism_selection"]
        self.assertEqual(selection["realized_patterns"], ["+--+", "-++-"])
        self.assertEqual(selection["realized_count"], 2)
        self.assertEqual(selection["J_sign_pattern"], [1, -1, -1, 1])
        self.assertTrue(selection["J_is_antipode_permutation"])
        self.assertTrue(selection["identity_verified"])
        self.assertEqual(
            selection["central_involution_identity"],
            "J = (A^3 - 4*A^2 - 5*A + 10*I)/10",
        )
        self.assertEqual(len(selection["non_realized"]), 12)
        for row in selection["non_realized"]:
            self.assertNotIn(row["pattern"], ("+--+", "-++-", "++++", "----"))
            self.assertIn("value", row["offending_entry"])

    def test_selected_producer_response_is_minus_J(self) -> None:
        selected = self.expected["selected_producer_response"]
        self.assertEqual(selected["pattern"], [-1, 1, 1, -1])
        self.assertEqual(selected["pattern_name"], "-++-")
        self.assertTrue(selected["realized_as_signed_automorphism"])
        self.assertIn("#566", selected["convention"])
        reference = self.expected["producer_reference"]
        self.assertEqual(
            reference["manifest_sha256"],
            self.manifest["producer_manifest_sha256"],
        )

    def test_countermodel_pair_is_inequivalent_and_load_bearing(self) -> None:
        pair = self.expected["countermodel_pair"]
        self.assertEqual(pair["clause_dropped"], "signed_graph_automorphism_readback")
        self.assertTrue(pair["inequivalent"])
        self.assertTrue(
            pair["both_satisfy_equivariance_involution_signedness_nontriviality"]
        )
        model_a, model_b = pair["model_a"], pair["model_b"]
        self.assertEqual(model_a["pattern"], "+--+")
        self.assertTrue(model_a["signed_graph_automorphism"])
        self.assertEqual(model_a["trace"], "0")
        self.assertEqual(model_a["entry_values"], ["0", "1"])
        self.assertEqual(model_b["pattern"], "+++-")
        self.assertFalse(model_b["signed_graph_automorphism"])
        self.assertEqual(model_b["trace"], "2")
        self.assertEqual(
            model_b["offending_entry"], {"row": 0, "col": 0, "value": "1/6"}
        )
        self.assertNotEqual(model_a["entry_values"], model_b["entry_values"])

    def test_presentation_invariance_and_verdict(self) -> None:
        invariance = self.expected["presentation_invariance"]
        self.assertEqual(invariance["rotations_checked"], 60)
        self.assertTrue(invariance["group_closure_verified"])
        self.assertTrue(invariance["conjugation_fixes_every_admissible_response"])
        self.assertEqual(invariance["verdict"], "presentation_invariant")
        verdict = self.expected["verdict"]
        self.assertEqual(verdict["grammar_internal_uniqueness"], "up_to_sign")
        self.assertEqual(verdict["axiom_forcing"], "independence_limited")
        self.assertEqual(verdict["operational_clause"], "load_bearing")
        self.assertEqual(verdict["bounded_exit"], "independence_limited")
        self.assertEqual(
            self.expected["claim_boundary"]["status"], "exact_named_realization"
        )

    def test_projector_construction_recomputes(self) -> None:
        _, adjacency, _, antipode = cert.port_model()
        projectors = cert.isotypic_projectors(adjacency)
        j_mat = cert.central_involution(adjacency, antipode, projectors)
        self.assertTrue(
            cert.mat_eq(
                cert.sign_pattern_matrix(projectors, cert.PATTERN_J), j_mat
            )
        )
        minus_j = cert.sign_pattern_matrix(projectors, cert.PATTERN_MINUS_J)
        self.assertTrue(cert.mat_eq(minus_j, cert.smul(cert.MINUS_ONE, j_mat)))

    def test_tampered_adjacency_fails_closed(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.tampered_adjacency_control()
        self.assertEqual(caught.exception.code, "PROJECTOR_SYSTEM")

    def test_fifteenth_admissible_pattern_is_rejected(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.fifteenth_pattern_control()
        self.assertEqual(caught.exception.code, "ADMISSIBLE_SET")

    def test_non_automorphism_promotion_is_rejected(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.non_automorphism_promotion_control()
        self.assertEqual(caught.exception.code, "SIGNED_AUTOMORPHISM")

    def test_negative_controls_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        stored = cert.load_json(self.negative_path)
        self.assertEqual(payload, stored)
        names = {row["name"] for row in payload["finite_controls"]}
        self.assertEqual(
            names,
            {
                "wrong_schema",
                "involution_clause_dropped",
                "operational_clause_dropped",
                "plus_J_promotion",
                "producer_pin_drift",
                "mass_target_injection",
                "source_firewall_token",
                "tampered_adjacency_breaks_projectors",
                "fifteenth_admissible_pattern_rejected",
                "non_automorphism_promotion_rejected",
            },
        )
        self.assertTrue(all(row["passed"] for row in payload["finite_controls"]))

    def test_tampered_receipt_is_rejected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["signed_automorphism_selection"]["realized_count"] = 3
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(caught.exception.code, "RECEIPT_MISMATCH")


if __name__ == "__main__":
    unittest.main()
