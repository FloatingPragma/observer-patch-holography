#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #624 noncentral routed-seam
reduction certificate."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import noncentral_seam_reduction_certificate as cert  # noqa: E402


class NoncentralSeamReductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = (
            MODULE_DIR / "manifests" / "noncentral_seam_reduction_reference.json"
        )
        cls.receipt_path = (
            MODULE_DIR / "receipts" / "noncentral_seam_reduction_reference.receipt.json"
        )
        cls.negative_path = (
            MODULE_DIR / "negative_controls" / "issue_624_negative_controls.json"
        )
        cls.manifest = cert.load_json(cls.manifest_path)
        cls.expected = cert.certificate_payload(cls.manifest)
        cls.spin_artifact = cert.load_json(
            MODULE_DIR / "manifests" / "spin_statistics_semantic_artifact.json"
        )
        quaternions, _ = cert.pinned_lift_quaternions(cls.spin_artifact)
        cls.s3 = cert.build_symmetric_group_three()
        cls.q8 = cert.build_quaternion_group_eight(quaternions)

    def test_reference_receipt_is_exactly_recomputable(self) -> None:
        receipt = cert.load_json(self.receipt_path)
        cert.verify_receipt(self.manifest, receipt)
        self.assertEqual(receipt, self.expected)

    def test_negative_controls_file_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(self.negative_path)
        recomputed = cert.negative_control_payload(self.manifest)
        self.assertEqual(stored, recomputed)
        self.assertTrue(all(row["passed"] for row in stored["finite_controls"]))
        names = {row["name"] for row in stored["finite_controls"]}
        self.assertIn("spin_manifest_pin_tampered", names)
        self.assertIn("routed_manifest_pin_drift", names)

    def test_value_groups_are_reconstructed_from_source(self) -> None:
        groups = self.expected["value_groups"]
        two_i = groups["binary_icosahedral"]
        self.assertEqual(two_i["order"], 120)
        self.assertEqual(two_i["centre_order"], 2)
        self.assertTrue(two_i["order_profile_matches_pinned_lift_measurement"])
        self.assertTrue(two_i["contains_pinned_klein_four_lifts"])
        self.assertEqual(groups["s3"], {
            "order": 6,
            "centre_order": 1,
            "source": "the symmetric group on three letters as permutation tuples",
        })
        self.assertEqual(groups["q8"]["order"], 8)
        self.assertEqual(groups["q8"]["centre_order"], 2)

    def test_flat_part_is_gauge_trivial_exact(self) -> None:
        flat = self.expected["flat_part"]
        self.assertEqual(flat["status"], "gauge_trivial_exact")
        self.assertEqual(flat["spanning_tree_seams"], 11)
        self.assertEqual(len(flat["worklist_schedule"]), 19)
        self.assertTrue(flat["schedule_value_independent"])
        names = [row["group"] for row in flat["per_group"]]
        self.assertEqual(names, ["2I", "S3", "Q8"])
        for row in flat["per_group"]:
            self.assertEqual(row["bijection_section_checks"], cert.GAUGE_TRIALS)
            self.assertEqual(row["bijection_retraction_checks"], cert.GAUGE_TRIALS)
            self.assertTrue(row["identity_assignment_maps_to_trivial_gauge"])
            self.assertGreaterEqual(row["twisted_candidates_rejected"], 1)
        counts = {row["group"]: row["coherent_assignment_count"] for row in flat["per_group"]}
        self.assertEqual(counts["S3"], f"6^11 = {6 ** 11}")
        self.assertEqual(counts["Q8"], f"8^11 = {8 ** 11}")
        self.assertEqual(counts["2I"], f"120^11 = {120 ** 11}")

    def test_reduced_exhaustive_census_matches_the_count_identity(self) -> None:
        reduced = {row["group"]: row for row in self.expected["flat_part"]["reduced_exhaustive"]}
        self.assertEqual(set(reduced), {"S3", "Q8"})
        self.assertEqual(reduced["S3"]["assignments_enumerated"], 6 ** 5)
        self.assertEqual(reduced["S3"]["coherent_found"], 6 ** 3)
        self.assertEqual(reduced["Q8"]["assignments_enumerated"], 8 ** 5)
        self.assertEqual(reduced["Q8"]["coherent_found"], 8 ** 3)
        for row in reduced.values():
            self.assertTrue(row["every_coherent_assignment_is_a_gauge"])

    def test_central_obstruction_is_centre_valued_exact(self) -> None:
        central = self.expected["central_obstruction"]
        self.assertEqual(central["status"], "centre_valued_exact")
        smith = central["h2_smith"]
        self.assertEqual(smith["coboundary_smith_invariants"], [1] * 19)
        self.assertEqual(smith["integer_tube_generators"], 19)
        self.assertEqual(smith["h2_by_coefficient_group"], {"Z2": 2, "Z3": 3, "Z6": 6})
        self.assertEqual(
            smith["h2_by_value_group_centre"],
            {"Z(2I) = Z2": 2, "Z(Q8) = Z2": 2, "Z(S3) = 1": 1},
        )
        lifted = central["lifted_witness"]
        self.assertEqual(lifted["quotient_group_order"], 60)
        self.assertTrue(lifted["commutes_with_every_seam_value"])
        self.assertTrue(lifted["gauge_invariant"])
        self.assertGreater(lifted["minus_one_faces"], 0)
        self.assertEqual(lifted["minus_one_faces"] % 2, 0)
        twist = central["central_twist_flux"]
        self.assertEqual(twist["total_class"], "trivial")
        self.assertIn("nontrivial Z2 class", twist["nontrivial_class_witness"])

    def test_measured_transport_matches_the_nontrivial_z2_class(self) -> None:
        section = self.expected["central_obstruction"]["section_obstruction_match"]
        self.assertEqual(section["sign_assignments_tested"], 8)
        self.assertEqual(section["sections_found"], 0)
        self.assertTrue(section["matches_pinned_section_obstruction"])
        self.assertEqual(
            section["recovered_class"], "the nontrivial class of the order-two centre"
        )
        self.assertEqual(self.expected["verdict"]["measured_transport_class"], "nontrivial_Z2")

    def test_consequence_keeps_the_order_six_menu(self) -> None:
        consequence = self.expected["consequence"]
        self.assertEqual(consequence["measured_flux_menu"], [0, 1, 2, 3, 4, 5])
        self.assertTrue(consequence["flux_menu_exhausted_by_central_classes"])
        self.assertEqual(
            consequence["transport_centre_embedding"],
            {"centre_order": 2, "menu_subgroup": [0, 3]},
        )
        self.assertEqual(
            consequence["axis_class_lattice"]["smith_invariants"], [1, 1, 1, 1, 1, 6]
        )

    def test_out_of_class_control_records_the_flat_sectors(self) -> None:
        control = self.expected["consequence"]["out_of_class_control"]
        self.assertTrue(control["worklist_stalls"])
        self.assertTrue(control["flat_nontrivial_holonomy_witness_not_a_gauge"])
        self.assertEqual(
            control["flat_sector_count_equals_conjugacy_classes"],
            {"2I": 9, "Q8": 5, "S3": 3},
        )

    def test_verdict_closes_the_general_grammar_on_the_spherical_class(self) -> None:
        verdict = self.expected["verdict"]
        self.assertEqual(
            verdict["general_grammar"],
            "exact_reduction_to_central_obstructions_on_simply_connected_nerve",
        )
        self.assertTrue(verdict["measured_order_six_menu_stands"])
        self.assertEqual(
            verdict["central_class_group_orders"],
            {"binary_icosahedral": 2, "q8": 2, "s3": 1},
        )
        self.assertEqual(
            self.expected["claim_boundary"]["bounded_exit"], "exact_named_realization"
        )
        self.assertTrue(all(verdict["controls"].values()))

    # -- adversarial unit checks on the algorithm itself ---------------------

    def test_incoherent_assignment_fails_closed(self) -> None:
        cx = cert.make_complex(4, cert.REDUCED_FACES)
        gauge = [0, 1, 2, 3]
        assignment = cert.pure_gauge_assignment(cx, self.s3, gauge)
        twist = min(i for i in range(self.s3.order) if i != self.s3.identity)
        assignment[0] = self.s3.mul(assignment[0], twist)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.trivialize(cx, self.s3, assignment)
        self.assertEqual(caught.exception.code, "FACE_COHERENCE")

    def test_noncentral_obstruction_fails_the_commutation_check(self) -> None:
        cx = cert.make_complex(4, cert.REDUCED_FACES)
        gauge = [1, 3, 5, 7]
        assignment = cert.pure_gauge_assignment(cx, self.q8, gauge)
        noncentral = min(i for i in range(self.q8.order) if i not in self.q8.centre)
        assignment[0] = self.q8.mul(assignment[0], noncentral)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.central_discrepancy_certificate(cx, self.q8, assignment)
        self.assertEqual(caught.exception.code, "NONCENTRAL_OBSTRUCTION")

    def test_central_twist_passes_the_commutation_check_on_the_reduced_complex(self) -> None:
        cx = cert.make_complex(4, cert.REDUCED_FACES)
        assignment = cert.pure_gauge_assignment(cx, self.q8, [0, 2, 4, 6])
        minus_one = cert.minus_one_index(self.q8)
        assignment[0] = self.q8.mul(assignment[0], minus_one)
        certificate = cert.central_discrepancy_certificate(cx, self.q8, assignment)
        self.assertTrue(certificate["centre_membership"])
        self.assertTrue(certificate["gauge_invariant"])

    def test_s3_extra_flux_sector_is_rejected(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.admit_flux_sectors([0, 1], len(self.s3.centre))
        self.assertEqual(caught.exception.code, "EXTRA_FLUX")

    def test_faceless_cycle_stalls_the_worklist(self) -> None:
        cycle = cert.make_complex(3, [], extra_edges=[(0, 1), (0, 2), (1, 2)])
        with self.assertRaises(cert.CertificateError) as caught:
            cert.trivialize(cycle, self.s3, [self.s3.identity, self.s3.identity, 1])
        self.assertEqual(caught.exception.code, "FLAT_PROPAGATION")

    def test_tampered_spin_pin_fails_closed(self) -> None:
        mutant = copy.deepcopy(dict(self.manifest))
        mutant["spin_statistics_artifact_sha256"] = "sha256:" + "0" * 64
        with self.assertRaises(cert.CertificateError) as caught:
            cert.certificate_payload(mutant)
        self.assertEqual(caught.exception.code, "UPSTREAM_HASH")

    def test_forbidden_promotion_keys_fail_closed(self) -> None:
        for key in ("noncentral_flux_sector", "flat_sector_promotion", "instanton_sector"):
            mutant = copy.deepcopy(dict(self.manifest))
            mutant[key] = {"declared_without_source_receipt": True}
            with self.assertRaises(cert.CertificateError) as caught:
                cert.certificate_payload(mutant)
            self.assertEqual(caught.exception.code, "FORBIDDEN_DEPENDENCY")


if __name__ == "__main__":
    unittest.main()
