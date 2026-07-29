#!/usr/bin/env python3
"""Regression tests for the finite topological packet associated with #311."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import flux_defect_criterion_certificate as cert  # noqa: E402


def polynomial_product(*factors: list[int]) -> list[int]:
    result = [1]
    for factor in factors:
        out = [0] * (len(result) + len(factor) - 1)
        for i, a in enumerate(result):
            for j, b in enumerate(factor):
                out[i + j] += a * b
        result = out
    return result


class FluxDefectCriterionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "flux_defect_criterion_reference.json"
        cls.receipt_path = (
            MODULE_DIR / "receipts" / "flux_defect_criterion_reference.receipt.json"
        )
        cls.negative_path = (
            MODULE_DIR / "negative_controls" / "issue_311_negative_controls.json"
        )
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_untwisted_spectrum_is_the_icosahedral_anchor(self) -> None:
        # (x - 5)(x + 1)^5 (x^2 - 5)^3, expanded by exact integer arithmetic.
        anchor = polynomial_product(
            [1, -5],
            [1, 1],
            [1, 1],
            [1, 1],
            [1, 1],
            [1, 1],
            [1, 0, -5],
            [1, 0, -5],
            [1, 0, -5],
        )
        per_class = self.expected["twisted_spectra"]["per_class"]
        self.assertEqual(per_class["0"]["characteristic_polynomial"], anchor)
        self.assertEqual(per_class["0"]["delta_from_untwisted"], [0] * 13)

    def test_spectral_criterion_distinctness_and_pairing(self) -> None:
        verdict = self.expected["spectral_criterion"]
        self.assertTrue(verdict["class_0_differs_from_class_3"])
        self.assertTrue(verdict["class_0_differs_from_class_1"])
        self.assertTrue(verdict["class_1_differs_from_class_3"])
        self.assertFalse(verdict["all_six_pairwise_distinct"])
        self.assertEqual(
            verdict["coincidence_partition"], [[0], [1, 5], [2, 4], [3]]
        )
        self.assertEqual(verdict["coincident_unordered_pairs"], [[1, 5], [2, 4]])
        self.assertEqual(len(verdict["distinct_unordered_pairs"]), 13)
        per_class = self.expected["twisted_spectra"]["per_class"]
        self.assertEqual(
            per_class["1"]["characteristic_polynomial"],
            per_class["5"]["characteristic_polynomial"],
        )
        self.assertNotEqual(
            per_class["0"]["characteristic_polynomial"],
            per_class["3"]["characteristic_polynomial"],
        )

    def test_hilbert_representation_and_gauge_invariance(self) -> None:
        hilbert = self.expected["hilbert_representation"]
        self.assertTrue(hilbert["self_adjoint"])
        self.assertEqual(hilbert["hermiticity_entry_checks"], 6 * 60)
        gauge = self.expected["gauge_invariance"]
        self.assertTrue(gauge["characteristic_polynomials_invariant"])
        self.assertEqual(len(gauge["sample_gauges"]), 4)
        self.assertEqual(gauge["diagonal_conjugation_entry_checks"], 4 * 6 * 60)

    def test_local_coboundary_control_keeps_physical_boundary_open(self) -> None:
        control = self.expected["local_coboundary_control"]
        self.assertEqual(control["scaled_basis_coboundaries_with_zero_holonomy"], 60)
        self.assertTrue(control["marked_vertex_leaves_operator_family_unchanged"])
        self.assertTrue(
            control["local_coboundary_regauge_preserves_class_and_spectra"]
        )
        self.assertEqual(control["ordered_class_pairs_with_no_connecting_regauge"], 30)
        self.assertFalse(control["excludes_every_classical_localized_record"])
        self.assertFalse(control["physical_particle_discrimination"])
        self.assertIn(
            "two-component harmonic vector-spring",
            control["same_support_classical_completion"],
        )

    def test_exact_support_vector_spring_realification(self) -> None:
        control = self.expected["exact_support_classical_realification"]
        self.assertEqual(
            control["domain"],
            {
                "name": "certified_icosahedral_support",
                "vertex_count": 12,
                "seam_count": 30,
                "regular_degree": 5,
                "separate_from_local_domain_seam_complex": True,
            },
        )
        self.assertIn(
            "two-component harmonic vector-spring",
            control["classical_model"],
        )
        self.assertEqual(control["stiffness_rule"], "K_k = 5 I - A_k")
        self.assertEqual(control["phase_unit_checks"], 6 * 60)
        self.assertEqual(control["hessian_entry_identity_checks"], 6 * 72)
        self.assertEqual(control["phase_metric_isometry_checks"], 6 * 30)
        self.assertEqual(
            control["edge_hessian_identity_entry_checks"],
            6 * 24 * 24,
        )
        self.assertEqual(len(control["per_class"]), 6)
        self.assertTrue(
            control["declared_adjacency_spectral_family_recoverable"]
        )
        self.assertFalse(control["scalar_operator_or_gap_matched"])
        self.assertFalse(control["complete_interface_ontology_no_go"])
        self.assertFalse(control["extended_domain_non_identifiability"])
        per_class = self.expected["twisted_spectra"]["per_class"]
        for row in control["per_class"]:
            self.assertEqual(
                row["adjacency_characteristic_polynomial"],
                per_class[str(row["flux_class"])][
                    "characteristic_polynomial"
                ],
            )
            self.assertEqual(
                row["real_generalized_characteristic_polynomial"],
                cert.polynomial_square(
                    row["stiffness_characteristic_polynomial"]
                ),
            )
            hessian = row["classical_hessian_certificate"]
            self.assertEqual(hessian["coordinate_metric"], [[2, 1], [1, 2]])
            self.assertEqual(
                hessian["metric_leading_principal_minors"],
                [2, 3],
            )
            self.assertTrue(hessian["metric_positive_definite"])
            self.assertEqual(hessian["phase_metric_isometry_checks"], 30)
            self.assertTrue(
                hessian["energy_hessian_equals_metric_times_stiffness"]
            )
            self.assertTrue(hessian["energy_hessian_symmetric"])
            self.assertTrue(
                hessian["positive_semidefinite_by_edge_sum_of_squares"]
            )
        with self.assertRaises(cert.CertificateError) as caught:
            cert.assert_exact_support_spectra_exclude_classical(control)
        self.assertEqual(caught.exception.code, "SPECTRAL_PROJECTION_CONTROL")

    def test_local_domain_vector_spring_context_is_separate_and_bounded(self) -> None:
        control = self.expected["local_domain_classical_spectral_context"]
        self.assertEqual(
            control["verdict"],
            "CLASSICAL_REALIZATION_MATCHES_DECLARED_FINITE_SPECTRAL_INTERFACE",
        )
        self.assertTrue(control["separate_from_exact_flux_support"])
        self.assertFalse(control["exact_flux_identity_bridge"])
        interface = control["spectral_interface_identity"]
        self.assertEqual(
            interface["producer_schema"],
            "oph.local-domain-defect-sector-spectra.v1",
        )
        self.assertEqual(interface["main_domain"]["visible_node_count"], 8662)
        self.assertEqual(interface["main_domain"]["visible_edge_count"], 11816)
        self.assertEqual(
            interface["ladder_domain"]["source_carrier_count"],
            2048,
        )
        self.assertEqual(
            interface["ladder_domain"]["visible_node_count"],
            1052,
        )
        self.assertEqual(
            interface["ladder_domain"]["visible_edge_count"],
            1663,
        )
        self.assertFalse(
            interface["rer_exact_flux_12_42_vertex_identity_bridge"]
        )
        self.assertIn(
            "two-component classical harmonic network",
            control["classical_model"],
        )
        self.assertTrue(control["sector_payload_identity"])
        self.assertTrue(control["scalar_gap_payload_identity"])
        self.assertTrue(control["ladder_payload_identity"])
        readback = control["numerical_readback"]
        self.assertLess(
            readback["maximum_sector_gap_abs_residual"],
            control["numerical_gates"]["gap_match_abs_tolerance"],
        )
        self.assertLess(
            readback["scalar_gap_abs_residual"],
            control["numerical_gates"]["gap_match_abs_tolerance"],
        )
        self.assertLess(
            readback["maximum_ladder_gap_abs_residual"],
            control["numerical_gates"]["gap_match_abs_tolerance"],
        )
        self.assertTrue(readback["all_kernel_counts_match_exact_doubling"])
        self.assertFalse(control["complete_interface_ontology_no_go"])
        self.assertFalse(control["extended_domain_non_identifiability"])
        self.assertFalse(
            control["bounded_lexical_census"]["completeness_theorem"]
        )

    def test_charge_fusion_and_finite_chain_composition(self) -> None:
        fusion = self.expected["charge_and_fusion"]
        self.assertEqual(fusion["fusion_pairs_checked"], 36)
        self.assertEqual(fusion["multi_defect_pairs_checked"], 36)
        self.assertTrue(fusion["selected_seam_path_supports_disjoint"])
        first = set(self.expected["defect_sector_definition"]["witness_seam_support"])
        second = set(fusion["second_pair_seam_support"])
        self.assertFalse(first & second)
        self.assertFalse(fusion["puncture_pairs_are_asymptotically_separated"])
        self.assertTrue(fusion["puncture_face_shared_boundary_edges"])
        self.assertIn("no asymptotic state", fusion["composition_scope"])

    def test_one_step_refinement_check_does_not_claim_a_tower(self) -> None:
        refinement = self.expected["one_step_refinement_check"]
        self.assertEqual(
            refinement["refined_complex"], {"vertices": 42, "seams": 120, "faces": 80}
        )
        self.assertTrue(refinement["matches_pinned_refinement_counts"])
        self.assertEqual(refinement["realized_flux_menu"], [0, 1, 2, 3, 4, 5])
        self.assertEqual(
            refinement["coincidence_partition"], [[0], [1, 5], [2, 4], [3]]
        )
        self.assertTrue(refinement["one_step_partition_persists"])
        self.assertFalse(refinement["coarse_chain_transported_to_refined_chain"])
        self.assertFalse(refinement["coarse_to_refined_operator_intertwiner_proved"])
        self.assertFalse(refinement["all_depth_refinement_stability_proved"])

    def test_gate_rows_and_open_interfaces(self) -> None:
        gate = self.expected["finite_topological_sector_gate"]
        self.assertTrue(gate["passed"])
        self.assertTrue(gate["local_coboundary_regauge_cannot_change_flux_class"])
        self.assertFalse(gate["physical_particle_discrimination"])
        self.assertFalse(gate["classical_same_flux_countermodel_excluded"])
        self.assertFalse(gate["asymptotic_states_controlled"])
        self.assertFalse(gate["all_depth_refinement_stability"])
        self.assertFalse(gate["continuum_quantum_pole"])
        self.assertFalse(gate["scattering_amplitude_interface"])
        self.assertFalse(gate["laboratory_identification"])
        boundary = self.expected["claim_boundary"]
        self.assertEqual(
            boundary["status"],
            "finite_topological_flux_spectrum_with_exact_support_"
            "classical_realification_and_separate_local_domain_context",
        )
        self.assertGreaterEqual(len(boundary["does_not_close"]), 7)
        exact_boundary = self.expected["exact_support_classical_boundary"]
        self.assertEqual(
            exact_boundary["status"], "exact_classical_realification"
        )
        self.assertTrue(
            exact_boundary["declared_adjacency_spectral_family_recoverable"]
        )
        self.assertFalse(exact_boundary["scalar_operator_or_gap_matched"])
        self.assertFalse(
            exact_boundary["complete_interface_discriminators_covered"]
        )
        local_boundary = self.expected[
            "local_domain_spectral_context_boundary"
        ]
        self.assertEqual(
            local_boundary["status"],
            "bounded_local_domain_classical_match",
        )
        self.assertFalse(local_boundary["same_domain_as_exact_flux_support"])
        self.assertFalse(
            local_boundary["identity_bridge_to_exact_flux_support"]
        )
        self.assertFalse(local_boundary["extended_domain_no_go"])
        self.assertFalse(local_boundary["lexical_scan_completeness_theorem"])

        acceptance = self.expected["acceptance_criteria_status"]
        self.assertTrue(acceptance["defect_object_and_equivalence_source_defined"])
        self.assertTrue(acceptance["charge_invariant_and_target_independent"])
        self.assertFalse(acceptance["mass_invariant_and_target_independent"])
        self.assertFalse(
            acceptance["mass_and_charge_invariant_and_target_independent"]
        )
        self.assertFalse(
            acceptance["multi_defect_composition_and_asymptotic_states_controlled"]
        )
        self.assertFalse(
            acceptance[
                "quantum_pole_or_equivalent_physical_spectral_criterion_proved"
            ]
        )
        self.assertFalse(
            acceptance["classical_localization_alone_cannot_pass_physical_gate"]
        )
        self.assertFalse(acceptance["all_depth_refinement_stability_proved"])
        self.assertFalse(acceptance["all_issue_acceptance_criteria_satisfied"])

    def test_negative_controls_all_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        self.assertEqual(cert.load_json(self.negative_path), payload)
        self.assertTrue(all(row["passed"] for row in payload["manifest_controls"]))
        self.assertTrue(all(row["passed"] for row in payload["tamper_controls"]))
        self.assertEqual(
            [row["name"] for row in payload["tamper_controls"]],
            [
                "non_coboundary_presented_as_regauge",
                "equal_spectra_claim_classes_0_3",
                "equal_spectra_claim_classes_0_1",
                "non_hermitian_seam_tamper",
                "gauge_dependence_tamper",
                "declared_spectra_exclude_classical_claim",
                "regular_block_claimed_euclidean_orthogonal",
            ],
        )
        self.assertIn(
            "exact_support_classical_completion",
            payload["countermodel_witnesses"],
        )
        self.assertIn(
            "local_domain_context_boundary",
            payload["countermodel_witnesses"],
        )

    def test_manifest_pin_drift_is_rejected(self) -> None:
        mutant = copy.deepcopy(dict(self.manifest))
        mutant["carrier_manifest_sha256"] = "0" * 64
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "UPSTREAM_HASH")

    def test_forbidden_target_keys_are_rejected(self) -> None:
        for key in cert.FORBIDDEN_MANIFEST_KEYS:
            mutant = copy.deepcopy(dict(self.manifest))
            mutant[key] = True
            with self.assertRaises(cert.CertificateError) as caught:
                cert.certificate_payload(mutant)
            self.assertEqual(caught.exception.code, "FORBIDDEN_DEPENDENCY")

    def test_witness_tamper_is_rejected(self) -> None:
        support = cert.load_carrier_complex(self.manifest, MODULE_DIR)
        artifact = cert.load_global_form_artifact(
            self.manifest, MODULE_DIR, support["carrier_manifest_sha256"]
        )
        tampered = copy.deepcopy(artifact)
        witness = tampered["sector_menu"]["flux_tube_witnesses"][1]
        first_seam = sorted(witness["seam_values"])[0]
        witness["seam_values"][first_seam] = (witness["seam_values"][first_seam] + 1) % 6
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_witness_chains(support, tampered)
        self.assertEqual(caught.exception.code, "WITNESS")

    def test_eisenstein_arithmetic_identities(self) -> None:
        omega = cert.OMEGA_POWERS[1]
        self.assertEqual(cert.eis_mul(omega, omega), cert.OMEGA_POWERS[2])
        product = (1, 0)
        for _ in range(6):
            product = cert.eis_mul(product, omega)
        self.assertEqual(product, (1, 0))
        self.assertEqual(cert.eis_mul(omega, cert.eis_conj(omega)), (1, 0))
        for exponent in range(6):
            block = cert.omega_block_power(exponent)
            expected_a, expected_b = cert.OMEGA_POWERS[exponent]
            # a + b omega maps to a I + b W in the companion representation.
            identity = ((1, 0), (0, 1))
            omega_block = cert.OMEGA_BLOCK
            image = tuple(
                tuple(
                    expected_a * identity[i][j] + expected_b * omega_block[i][j]
                    for j in range(2)
                )
                for i in range(2)
            )
            self.assertEqual(block, image)


if __name__ == "__main__":
    unittest.main()
