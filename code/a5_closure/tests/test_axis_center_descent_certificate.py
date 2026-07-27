#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #567 descent certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import axis_center_descent_certificate as cert  # noqa: E402


class AxisCenterDescentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "axis_center_descent_reference.json"
        cls.receipt_path = MODULE_DIR / "receipts" / "axis_center_descent_reference.receipt.json"
        cls.negative_path = MODULE_DIR / "negative_controls" / "issue_567_negative_controls.json"
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_kernel_is_computed_on_realized_tensors(self) -> None:
        kernel = self.expected["kernel_on_realized_tensors"]
        self.assertEqual(kernel["candidates_enumerated"], 36)
        self.assertEqual(kernel["kernel_order"], 6)
        self.assertEqual(kernel["cyclic_generator"], [1, 1, 1])
        self.assertTrue(kernel["matches_emitted_kernel_data"])
        self.assertEqual(kernel["tensor_additivity_checked_pairs"], 49)

    def test_maximal_effective_image_and_dual_lattice(self) -> None:
        form = self.expected["maximal_effective_image"]
        self.assertEqual(form["group"], "(SU(3) x SU(2) x U(1)) / Z6")
        self.assertEqual(form["character_lattice_smith_invariants"], [1, 1, 6])
        self.assertEqual(form["character_residue_class_count"], 6)
        self.assertIn("not a selected physical", form["status"])
        lattice = self.expected["dual_cocharacter_lattice"]
        self.assertEqual(
            lattice["primitive_correlated_cocharacter"],
            {
                "color_coweight": "1/3",
                "weak_coweight": "1/2",
                "u1_charge_electron_dirac_units": "1/6",
            },
        )
        self.assertEqual(
            lattice["basis_coordinates_color_weak_u1"],
            [["1", "0", "0"], ["0", "1", "0"], ["1/3", "1/2", "1/6"]],
        )
        self.assertTrue(lattice["monopole_dynamics_not_inferred"])
        self.assertIn("not derived", lattice["theta_periodicity_status"])

    def test_four_global_forms_are_locally_indistinguishable(self) -> None:
        control = self.expected["global_form_nonidentifiability"]
        self.assertEqual(len(control["subgroup_menu"]), 4)
        self.assertEqual(
            [
                row["character_residue_class_count"]
                for row in control["quotient_candidates"]
            ],
            [36, 18, 12, 6],
        )
        self.assertTrue(
            all(
                row["all_declared_local_tensors_descend"]
                for row in control["quotient_candidates"]
            )
        )
        # Local tensors alone still cannot select; the selection is carried by
        # the measured sector menu, and the receipt records both facts.
        self.assertTrue(control["physical_global_form_selected"])
        self.assertEqual(control["selection_scope"], "finite source-model scope")
        self.assertIn("cannot select", control["selection_boundary"])
        self.assertEqual(control["adjoint_only_kernel_order"], 36)
        self.assertEqual(control["fractional_singlet_countermodel_kernel_order"], 1)
        self.assertIn("not every subgroup", control["subgroup_menu_scope"])

    def test_weak_center_relation_is_not_fermion_parity(self) -> None:
        relation = self.expected["algebraic_weak_center_u1_relation"]
        self.assertEqual(
            relation["h_cubed"],
            [0, 1, 3],
        )
        self.assertEqual(relation["phase_on_every_declared_weight_sixths"], 0)
        self.assertEqual(relation["universal_fermion_minus_one_candidates"], [])
        matter = {label: cert.REALIZED_WEIGHTS[label] for label in cert.MATTER_LABELS}
        self.assertEqual(cert.common_scalar_phase_elements(matter, 3), [])

    def test_weight_and_sector_refinement_are_proved(self) -> None:
        refinement = self.expected["weight_level_refinement_invariance"]
        self.assertEqual(refinement["carrier_rotations"], 60)
        self.assertEqual(refinement["artifact_persistence_maps"], 2)
        self.assertTrue(
            refinement["physical_loop_or_bundle_refinement_naturality_derived"]
        )

    def test_kernel_enumeration_matches_lean_module(self) -> None:
        kernel = cert.common_kernel(cert.REALIZED_WEIGHTS)
        self.assertEqual(set(kernel), cert.cyclic_generated((1, 1, 1)))
        self.assertEqual(len(cert.common_kernel(cert.ADJOINT_WEIGHTS)), 36)
        fractional = dict(cert.REALIZED_WEIGHTS)
        fractional["fractional_singlet"] = (1, 0, 0)
        self.assertEqual(cert.common_kernel(fractional), [(0, 0, 0)])

    def test_primitive_correlated_cocharacter_pairs_integrally(self) -> None:
        character_basis = [(1, 0, -2), (0, 1, -3), (0, 0, 6)]
        cocharacter = (Fraction(1, 3), Fraction(1, 2), Fraction(1, 6))
        pairings = [
            sum(Fraction(character[index]) * cocharacter[index] for index in range(3))
            for character in character_basis
        ]
        self.assertEqual(pairings, [0, 0, 1])

    def test_negative_controls_fail_closed(self) -> None:
        payload = cert.negative_control_payload(self.manifest)
        stored = cert.load_json(self.negative_path)
        self.assertEqual(payload, stored)
        names = {row["name"] for row in payload["finite_controls"]}
        self.assertEqual(
            names,
            {
                "wrong_hypercharge_convention",
                "unsupported_physical_promotion",
                "wrong_matter_receipt_pin",
                "monopole_dynamics_injection",
                "global_form_artifact_pin_drift",
                "wrong_global_form_artifact_schema",
                "declared_deck_loop_injection",
                "declared_spacetime_spin_injection",
                "declared_line_category_injection",
                "declared_uv_polarization_injection",
                "declared_instanton_sector_injection",
                "declared_theta_period_injection",
            },
        )
        self.assertTrue(all(row["passed"] for row in payload["finite_controls"]))

    def test_physical_global_form_gate_passes_with_named_open_lanes(self) -> None:
        self.assertTrue(self.expected["conditional_algebraic_gate"]["passed"])
        gate = self.expected["physical_global_form_gate"]
        self.assertTrue(gate["passed"])
        self.assertTrue(gate["source_derived_deck_loop_class"])
        self.assertTrue(gate["screen_spin_attachment_at_source_scope"])
        self.assertTrue(gate["genuine_line_category_selected"])
        self.assertTrue(gate["uv_mutual_locality_polarization_selected"])
        self.assertTrue(gate["screen_flux_sector_lattice_and_witnesses_measured"])
        # Continuum and laboratory rows stay false and are excluded from
        # "passed" through the explicit composition block.
        self.assertFalse(gate["four_dimensional_instanton_action_normalization"])
        self.assertFalse(gate["theta_periodicity_derived"])
        self.assertFalse(gate["laboratory_global_form_attachment"])
        composition = gate["composition"]
        for row in composition["passed_over"]:
            self.assertTrue(gate[row], row)
        self.assertIn("theta_periodicity_derived", composition["deferred"])
        self.assertIn("laboratory_global_form_attachment", composition["deferred"])
        self.assertTrue(self.expected["issue_closure_condition"]["met_locally"])

    def test_sector_transport_and_polarization_are_exhaustive(self) -> None:
        consistency = self.expected["sector_transport_consistency"]
        self.assertEqual(consistency["monodromy_checks"], 42)
        self.assertTrue(consistency["all_measured_sectors_carry_realized_matter"])
        self.assertEqual(
            consistency["fractional_singlet_obstructed_sectors"], [1, 2, 3, 4, 5]
        )
        self.assertEqual(consistency["unique_menu_matching_form"], "z6_quotient")
        polarization = self.expected["line_polarization"]
        self.assertEqual(polarization["maximal_mutually_local_lattices"], 12)
        self.assertEqual(polarization["lattices_containing_realized_electric"], 1)
        deck = self.expected["source_deck_loop_measurement"]
        self.assertEqual(deck["deck_group_order"], 120)
        self.assertEqual(deck["six_axis_class_group_order"], 6)
        self.assertEqual(deck["measured_source_admissible_menu"], [0, 1, 2, 3, 4, 5])
        self.assertEqual(deck["reference_federation_sector_class"], 0)

    def test_tampered_receipt_is_rejected(self) -> None:
        receipt = copy.deepcopy(self.expected)
        receipt["kernel_on_realized_tensors"]["kernel_order"] = 3
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(caught.exception.code, "RECEIPT_MISMATCH")

    def test_smith_normal_form(self) -> None:
        reduced = cert.smith_normal_form([[1, 0, -2], [0, 1, -3], [0, 0, 6]])
        self.assertEqual([reduced[i][i] for i in range(3)], [1, 1, 6])
        reduced_two = cert.smith_normal_form([[2, 0], [0, 3]])
        self.assertEqual([reduced_two[i][i] for i in range(2)], [1, 6])


if __name__ == "__main__":
    unittest.main()
