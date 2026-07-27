#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #614 scheduler certificate."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path
from unittest import mock

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import a3_scheduler_kernel_certificate as cert  # noqa: E402


class A3SchedulerKernelCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_path = MODULE_DIR / "manifests" / "a3_scheduler_kernel_reference.json"
        cls.stored = cert.load_json(cls.manifest_path)
        cls.payload = cert.certificate_payload()
        cls.vertices, cls.adjacency, cls.edges = cert.build_move_simplex()

    def test_stored_manifest_verifies_and_recomputes(self) -> None:
        cert.verify_manifest(self.stored)
        body = {k: v for k, v in self.stored.items() if k != "manifest_sha256"}
        self.assertEqual(body, self.payload)
        self.assertEqual(
            self.stored["manifest_sha256"], "sha256:" + cert.sha256_json(body)
        )

    def test_rerun_is_deterministic_and_byte_identical(self) -> None:
        self.assertEqual(cert.certificate_payload(), self.payload)
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "rerun.json"
            self.assertEqual(cert.main(["certify", "--output", str(out)]), 0)
            self.assertEqual(out.read_bytes(), self.manifest_path.read_bytes())
            self.assertEqual(cert.main(["verify", "--manifest", str(out)]), 0)

    def test_schema_and_verdict_fields(self) -> None:
        self.assertEqual(self.payload["schema"], "oph.a3_scheduler_kernel_certificate.v1")
        self.assertEqual(self.payload["issue"], 614)
        verdict = self.payload["verdict"]
        self.assertEqual(verdict["kernel"], "uniform_1_over_30")
        self.assertEqual(verdict["kernel_status"], "A3_selected_exact")
        self.assertEqual(verdict["dynamics_tier"], "exact_named_realization")
        self.assertEqual(
            verdict["general_model_fairness_liveness"], "conditional_open_interface"
        )
        self.assertEqual(
            set(verdict["controls"]),
            {"wrong_reference", "reducible", "periodic", "presentation_dependent"},
        )
        for name, record in self.payload["controls"].items():
            self.assertEqual(record["verdict"], "fails_closed", name)
            self.assertIn("witness", record)
            self.assertIn("conclusion", record["witness"])

    def test_move_simplex_and_kernel_are_exact(self) -> None:
        simplex = self.payload["move_simplex"]
        self.assertEqual(simplex["ports"], 12)
        self.assertEqual(simplex["seams"], 30)
        self.assertEqual(simplex["degree"], 5)
        self.assertEqual(len(simplex["seams_by_port_pair"]), 30)
        selection = self.payload["a3_selection"]
        self.assertEqual(selection["selected_kernel_probability_per_seam"], "1/30")
        self.assertTrue(selection["kernel_sums_to_one"])
        deck = selection["deck_invariance"]
        self.assertEqual(deck["carrier_automorphism_group_order"], 120)
        self.assertEqual(deck["seam_orbit_size"], 30)
        walk = self.payload["induced_port_walk"]
        self.assertEqual(walk["per_step_seam_traversal_probability"], "1/30")
        self.assertTrue(walk["matches_selected_kernel"])

    def test_spectral_block_is_exact(self) -> None:
        contraction = self.payload["proved_dynamics"]["contraction"]
        self.assertEqual(contraction["trace_powers_of_A"], [12, 0, 60, 120])
        self.assertEqual(
            contraction["multiplicities_from_exact_trace_system"], [1, 3, 3, 5]
        )
        self.assertEqual(
            contraction["spectral_gap"],
            {
                "display": "1 - (1/5)*sqrt(5)",
                "rational_part": "1",
                "sqrt5_coefficient": "-1/5",
            },
        )
        self.assertEqual(contraction["second_largest_modulus"], "sqrt(5)/5")

    def test_fairness_bound_matches_exact_fraction(self) -> None:
        fairness = self.payload["proved_dynamics"]["fairness"]
        self.assertEqual(fairness["expected_first_visit_steps"], "30")
        self.assertEqual(fairness["horizon_steps"], 120)
        exact = Fraction(29, 30) ** 120
        self.assertEqual(fairness["unvisited_probability_at_horizon"], str(exact))
        self.assertEqual(
            fairness["unvisited_probability_numerator_digits"],
            len(str(exact.numerator)),
        )
        self.assertEqual(
            fairness["unvisited_probability_denominator_digits"],
            len(str(exact.denominator)),
        )

    def test_tampered_adjacency_breaks_the_spectral_identity(self) -> None:
        matrix = cert.adjacency_int_matrix(self.adjacency)
        u, v = self.edges[0]
        matrix[u][v] = 0
        matrix[v][u] = 0
        with self.assertRaises(cert.CertificateError) as caught:
            cert.spectral_certificate(matrix)
        self.assertEqual(caught.exception.code, "SPECTRAL")

    def test_wrong_reference_control_forced_to_pass_fails_closed(self) -> None:
        uniform = {edge: Fraction(1, 30) for edge in self.edges}
        with self.assertRaises(cert.CertificateError) as caught:
            cert.control_wrong_reference(self.edges, weights=uniform)
        self.assertEqual(caught.exception.code, "CONTROL_WRONG_REFERENCE")

    def test_reducible_control_forced_to_pass_fails_closed(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.control_reducible(self.edges, self.adjacency, subset=self.edges)
        self.assertEqual(caught.exception.code, "CONTROL_REDUCIBLE")

    def test_periodic_control_forced_to_pass_fails_closed(self) -> None:
        triangle = cert.enumerate_triangles(self.adjacency)[0]
        with self.assertRaises(cert.CertificateError) as caught:
            cert.control_periodic(self.adjacency, cycle=list(triangle))
        self.assertEqual(caught.exception.code, "CONTROL_PERIODIC")

    def test_presentation_control_forced_to_pass_fails_closed(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.control_presentation(
                self.vertices,
                self.adjacency,
                self.edges,
                rotations=[cert.IDENTITY_ROTATION],
            )
        self.assertEqual(caught.exception.code, "CONTROL_PRESENTATION")

    def test_forced_control_makes_certificate_exit_nonzero(self) -> None:
        def forced(edges, adjacency):
            return list(edges), {"first_triangle": None, "second_triangle": None}

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "forced.json"
            with mock.patch.object(cert, "default_reducible_edges", forced):
                self.assertEqual(cert.main(["certify", "--output", str(out)]), 1)
            self.assertFalse(out.exists())

    def test_verify_rejects_tampered_manifest(self) -> None:
        tampered = copy.deepcopy(self.stored)
        tampered["verdict"]["kernel"] = "uniform_1_over_29"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tampered.json"
            path.write_text(json.dumps(tampered), encoding="utf-8")
            self.assertEqual(cert.main(["verify", "--manifest", str(path)]), 1)

    def test_periodic_witness_is_a_bipartite_six_cycle(self) -> None:
        witness = self.payload["controls"]["periodic"]["witness"]
        self.assertEqual(witness["cycle_length"], 6)
        self.assertEqual(witness["period"], 2)
        even = set(witness["even_class_ports"])
        odd = set(witness["odd_class_ports"])
        self.assertEqual(even & odd, set())
        self.assertEqual(even | odd, set(witness["cycle_ports"]))

    def test_presentation_witness_uses_a_nonidentity_rotation(self) -> None:
        witness = self.payload["controls"]["presentation_dependent"]["witness"]
        self.assertFalse(witness["rotation_is_identity"])
        self.assertNotEqual(witness["induced_port_permutation"], list(range(12)))
        self.assertNotEqual(
            witness["first_move_reference_presentation"],
            witness["first_move_relabeled_presentation_physical_seam"],
        )
        agreement = witness["quotient_visible_agreement"]
        self.assertEqual(agreement["seam_count"], 30)
        self.assertEqual(agreement["degree_sequence"], [5] * 12)
        self.assertTrue(agreement["relabeled_seam_multiset_equals_reference"])


if __name__ == "__main__":
    unittest.main()
