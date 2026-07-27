#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #625 load-fiber and
quadratic-readback certificate."""

from __future__ import annotations

import copy
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import load_fiber_readback_certificate as cert  # noqa: E402


class LoadFiberReadbackTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = (
            MODULE_DIR / "manifests" / "load_fiber_readback_reference.json"
        )
        cls.expected = cert.build_payload()

    def test_reference_manifest_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(self.manifest_path)
        self.assertEqual(stored, self.expected)
        body = {
            key: value
            for key, value in stored.items()
            if key != "payload_sha256"
        }
        self.assertEqual(
            stored["payload_sha256"], "sha256:" + cert.sha256_json(body)
        )

    def test_rerun_writes_byte_identical_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            first_bytes = first.read_bytes()
            self.assertEqual(first_bytes, second.read_bytes())
            self.assertEqual(first_bytes, self.manifest_path.read_bytes())

    def test_schema_and_identity_fields(self) -> None:
        self.assertEqual(
            self.expected["schema"], "oph.load_fiber_readback_certificate.v1"
        )
        self.assertEqual(self.expected["issue"], 625)
        self.assertEqual(
            self.expected["generated_by"],
            "code/a5_closure/load_fiber_readback_certificate.py",
        )
        self.assertEqual(
            set(self.expected),
            {
                "schema",
                "issue",
                "generated_by",
                "lean_spine",
                "integer_fiber",
                "readback_quadratic_form",
                "integer_fiber_lane",
                "quadratic_readback_lane",
                "consumer_ledger",
                "controls",
                "claim_boundary",
                "payload_sha256",
            },
        )
        spine = self.expected["lean_spine"]
        self.assertEqual(spine["commutant_module"], "Lean/Screen/A5Commutant.lean")
        self.assertEqual(spine["unit_split_module"], "Lean/Screen/UnitSplit12.lean")
        self.assertEqual(
            spine["unit_split_theorem"],
            "OPH.UnitSplit12.unit_split_of_positive_sum",
        )

    def test_conclusion_fields(self) -> None:
        self.assertEqual(
            self.expected["integer_fiber"],
            "axiom_forced_from_A1_finiteness_and_atomicity",
        )
        self.assertEqual(
            self.expected["readback_quadratic_form"],
            "forced_to_identity_class_by_A3_reference",
        )
        self.assertEqual(
            self.expected["integer_fiber_lane"]["conclusion"],
            "axiom_forced_from_A1_finiteness_and_atomicity",
        )
        self.assertEqual(
            self.expected["quadratic_readback_lane"]["conclusion"],
            "forced_to_identity_class_by_A3_reference",
        )

    def test_integer_fiber_lane(self) -> None:
        lane = self.expected["integer_fiber_lane"]
        atomicity = lane["atomicity"]
        self.assertEqual(atomicity["atom_count"], 12)
        self.assertEqual(atomicity["record_projection_count"], 4096)
        self.assertEqual(atomicity["atom_multiplicities"], [0, 1])
        self.assertTrue(atomicity["atoms_pairwise_orthogonal"])
        self.assertTrue(atomicity["atoms_sum_to_unit"])
        self.assertTrue(atomicity["every_atom_minimal"])
        grammar = lane["load_observable_grammar"]
        self.assertEqual(grammar["defect_domain"], "integer_port_charges")
        self.assertEqual(grammar["declared_total_charge"], 12)
        self.assertEqual(grammar["value_fiber"], "Z")
        self.assertEqual(len(grammar["observables"]), 13)
        for row in grammar["observables"]:
            self.assertTrue(row["integer_combination_of_atom_counters"])
            self.assertTrue(
                all(type(c) is int for c in row["atom_coefficients"])
            )
        total_row = grammar["observables"][-1]
        self.assertEqual(total_row["name"], "total_charge")
        self.assertEqual(total_row["atom_coefficients"], [1] * 12)
        carrier_manifest = cert.load_json(
            MODULE_DIR / "manifests" / "echosahedral_federation_reference.json"
        )
        self.assertEqual(
            lane["carrier_manifest_sha256"],
            cert.sha256_json(carrier_manifest),
        )

    def test_second_order_theorem_fisher_coefficient_is_six(self) -> None:
        theorem = self.expected["quadratic_readback_lane"][
            "second_order_theorem"
        ]
        self.assertEqual(theorem["first_order_coefficients"], ["1"] * 12)
        self.assertTrue(theorem["first_order_vanishes_on_tangent_space"])
        self.assertEqual(theorem["tangent_basis_evaluations"], [0] * 11)
        self.assertEqual(theorem["second_order_diagonal"], ["6"] * 12)
        self.assertEqual(theorem["fisher_coefficient_per_port"], "6")
        self.assertEqual(
            theorem["fisher_quadratic_form"], "Q = diag(1 / (2 tau_p)) = 6 I"
        )
        self.assertTrue(theorem["positive_definite"])
        self.assertFalse(theorem["numerics_used"])

    def test_port_divergence_jet_is_exact(self) -> None:
        jet = cert.port_divergence_jet(Fraction(1, 12))
        self.assertEqual(
            jet, {(1, 1): Fraction(1), (2, 2): Fraction(6)}
        )
        tilted = cert.port_divergence_jet(Fraction(2, 13))
        self.assertEqual(
            tilted, {(1, 1): Fraction(1), (2, 2): Fraction(13, 4)}
        )
        with self.assertRaises(cert.CertificateError) as caught:
            cert.port_divergence_jet(Fraction(0))
        self.assertEqual(caught.exception.code, "REFERENCE_SUPPORT")
        with self.assertRaises(cert.CertificateError) as caught:
            cert.divergence_expansion((Fraction(1, 2),) * 12)
        self.assertEqual(caught.exception.code, "REFERENCE_NORMALIZATION")

    def test_invariance_menu(self) -> None:
        menu = self.expected["quadratic_readback_lane"]["invariance_menu"]
        self.assertEqual(menu["menu_dimension"], 4)
        self.assertTrue(menu["orbitals_symmetric"])
        self.assertTrue(menu["orbitals_equivariant_under_sixty_rotations"])
        self.assertTrue(menu["orbital_supports_partition_entries"])
        fisher = menu["fisher_form"]
        self.assertEqual(fisher["coefficient_vector"], ["6", "6", "6", "6"])
        self.assertEqual(
            fisher["isotypic_coefficients"],
            {"1": "6", "3": "6", "3p": "6", "5": "6"},
        )
        self.assertTrue(fisher["identity_class"])
        countermodel = menu["retained_countermodel"]
        self.assertTrue(countermodel["passes_incidence_equivariance"])
        self.assertTrue(countermodel["positive_definite"])
        self.assertFalse(countermodel["identity_class"])
        self.assertEqual(
            countermodel["isotypic_coefficients"],
            {
                "1": "11",
                "3": "6 + 1*sqrt(5)",
                "3p": "6 + -1*sqrt(5)",
                "5": "5",
            },
        )

    def test_consumer_ledger_rows_are_typed(self) -> None:
        ledger = self.expected["consumer_ledger"]
        self.assertEqual(len(ledger), len(cert.CONSUMERS))
        paths = [row["consumer"] for row in ledger]
        self.assertIn(
            "code/a5_closure/echosahedral_selector_certificate.py", paths
        )
        self.assertIn("Lean/Screen/UnitSplit12.lean", paths)
        self.assertIn("code/capacity_readback/F_candidate_capK.py", paths)
        for row in ledger:
            self.assertTrue(
                row["integer_fiber_source"] is not None
                or row["quadratic_form_source"] is not None,
                row["consumer"],
            )
            if row["integer_fiber_source"] is not None:
                self.assertEqual(
                    row["integer_fiber_source"], cert.INTEGER_FIBER_SOURCE
                )
            if row["quadratic_form_source"] is not None:
                self.assertEqual(
                    row["quadratic_form_source"], cert.QUADRATIC_FORM_SOURCE
                )
            self.assertTrue((cert.REPO_ROOT / row["consumer"]).is_file())

    def test_missing_consumer_fails_closed(self) -> None:
        original = cert.CONSUMERS
        cert.CONSUMERS = original + (
            ("code/a5_closure/absent_consumer.py", "absent", "absent", True, False),
        )
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.consumer_ledger()
            self.assertEqual(caught.exception.code, "CONSUMER_MISSING")
        finally:
            cert.CONSUMERS = original

    def test_controls_record_required_failures(self) -> None:
        controls = self.expected["controls"]
        self.assertEqual(
            set(controls),
            {
                "non_atomic_load_model",
                "tilted_reference",
                "linear_readback",
                "adjacency_form_claimed_forced",
            },
        )
        for name, verdict in controls.items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        non_atomic = controls["non_atomic_load_model"]["witness"]
        self.assertEqual(non_atomic["rejection_code"], "A1_FINITENESS")
        tilted = controls["tilted_reference"]["witness"]
        self.assertEqual(tilted["tau_prime"][0], "2/13")
        self.assertEqual(
            tilted["distinct_diagonal_values"], ["13/2", "13/4"]
        )
        self.assertFalse(tilted["proportional_to_identity"])
        linear = controls["linear_readback"]["witness"]
        self.assertEqual(linear["tangent_basis_evaluations"], ["0"] * 11)
        adjacency = controls["adjacency_form_claimed_forced"]["witness"]
        self.assertTrue(adjacency["passes_incidence_equivariance"])
        self.assertFalse(adjacency["identity_class"])
        self.assertEqual(
            adjacency["fisher_coefficient_vector"], ["6", "6", "6", "6"]
        )

    def test_non_atomic_model_rejected_exactly(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.accept_load_model(cert.NON_ATOMIC_LOAD_MODEL)
        self.assertEqual(caught.exception.code, "A1_FINITENESS")
        accepted = cert.accept_load_model(
            {
                "record_spectrum": cert.ACCEPTED_SPECTRUM,
                "primitive_central_atom_count": 12,
                "load_observable": cert.ACCEPTED_LOAD_TYPE,
            }
        )
        self.assertTrue(accepted["accepted"])
        with self.assertRaises(cert.CertificateError) as caught:
            cert.accept_load_model(
                {
                    "record_spectrum": cert.ACCEPTED_SPECTRUM,
                    "primitive_central_atom_count": 12,
                    "load_observable": "lebesgue_density_readback",
                }
            )
        self.assertEqual(caught.exception.code, "A1_LOAD_TYPE")

    def test_forced_atomic_countermodel_fails_closed(self) -> None:
        original = copy.deepcopy(cert.NON_ATOMIC_LOAD_MODEL)
        cert.NON_ATOMIC_LOAD_MODEL.update(
            {
                "record_spectrum": cert.ACCEPTED_SPECTRUM,
                "primitive_central_atom_count": 12,
                "load_observable": cert.ACCEPTED_LOAD_TYPE,
            }
        )
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.build_payload()
            self.assertEqual(caught.exception.code, "CONTROL_NOT_FAILED")
            with tempfile.TemporaryDirectory() as scratch:
                output = Path(scratch) / "forced.json"
                self.assertEqual(cert.main(["--output", str(output)]), 1)
                self.assertFalse(output.exists())
        finally:
            cert.NON_ATOMIC_LOAD_MODEL.clear()
            cert.NON_ATOMIC_LOAD_MODEL.update(original)

    def test_forced_uniform_tilted_reference_fails_closed(self) -> None:
        original = cert.TILTED_REFERENCE_TAU
        cert.TILTED_REFERENCE_TAU = cert.UNIFORM_REFERENCE_TAU
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.build_payload()
            self.assertEqual(caught.exception.code, "CONTROL_NOT_FAILED")
        finally:
            cert.TILTED_REFERENCE_TAU = original

    def test_forced_identity_adjacency_form_fails_closed(self) -> None:
        original = cert.ADJACENCY_FORM_COEFFICIENTS
        cert.ADJACENCY_FORM_COEFFICIENTS = (6, 0)
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.build_payload()
            self.assertEqual(caught.exception.code, "CONTROL_NOT_FAILED")
        finally:
            cert.ADJACENCY_FORM_COEFFICIENTS = original

    def test_tilted_positive_lane_reference_fails_closed(self) -> None:
        original = cert.UNIFORM_REFERENCE_TAU
        cert.UNIFORM_REFERENCE_TAU = cert.TILTED_REFERENCE_TAU
        try:
            with self.assertRaises(cert.CertificateError) as caught:
                cert.build_payload()
            self.assertEqual(caught.exception.code, "FISHER_NOT_IDENTITY")
        finally:
            cert.UNIFORM_REFERENCE_TAU = original

    def test_tampered_atoms_break_atomicity(self) -> None:
        good = [
            tuple(1 if k == i else 0 for k in range(12)) for i in range(12)
        ]
        overlapping = list(good)
        overlapping[1] = good[0]
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_record_lattice(overlapping)
        self.assertEqual(caught.exception.code, "ATOM_ORTHOGONALITY")
        doubled = list(good)
        doubled[0] = tuple(2 if k == 0 else 0 for k in range(12))
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_record_lattice(doubled)
        self.assertEqual(caught.exception.code, "ATOM_MULTIPLICITY")
        absent = list(good)
        absent[0] = (0,) * 12
        with self.assertRaises(cert.CertificateError) as caught:
            cert.verify_record_lattice(absent)
        self.assertEqual(caught.exception.code, "ATOM_COMPLETENESS")
        missing_flag = [{"atom_id": "e_p00", "primitive": False}] + [
            {"atom_id": f"e_p{i:02d}", "primitive": True} for i in range(1, 12)
        ]
        with self.assertRaises(cert.CertificateError) as caught:
            cert.atom_indicator_vectors(missing_flag)
        self.assertEqual(caught.exception.code, "ATOM_PRIMITIVITY")

    def test_equivariance_checker_detects_non_invariant_form(self) -> None:
        verts, adjacency, distance, antipode = cert.r611.port_model()
        rotations = cert.r611.rotation_permutations(verts, adjacency)
        fisher = [
            [Fraction(6) if i == j else Fraction(0) for j in range(12)]
            for i in range(12)
        ]
        self.assertTrue(cert.verify_equivariance(fisher, rotations))
        tilted = copy.deepcopy(fisher)
        tilted[0][0] = Fraction(7)
        self.assertFalse(cert.verify_equivariance(tilted, rotations))

    def test_no_floats_anywhere_in_payload(self) -> None:
        cert.require_no_floats(self.expected)
        poisoned = copy.deepcopy(self.expected)
        poisoned["quadratic_readback_lane"]["poison"] = 1 / 6
        with self.assertRaises(cert.CertificateError) as caught:
            cert.require_no_floats(poisoned)
        self.assertEqual(caught.exception.code, "FLOAT_FORBIDDEN")


if __name__ == "__main__":
    unittest.main()
