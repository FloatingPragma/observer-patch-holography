#!/usr/bin/env python3
"""Regression and adversarial tests for the directed-seam repair bridge."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import directed_seam_repair_certificate as cert  # noqa: E402


class DirectedSeamRepairCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.orbit = cert.incidence_and_orbit()
        cls.payload = cert.build_payload()

    def test_orientation_is_constructed_from_unoriented_incidence(self) -> None:
        carrier = self.orbit["carrier"]
        faces = cert.coherent_orientation(carrier.adjacency)
        self.assertEqual(len(faces), 20)
        self.assertEqual(
            {frozenset(face) for face in faces},
            {frozenset(face) for face in cert.unoriented_triangles(carrier.adjacency)},
        )
        directed = [
            seam for a, b, c in faces for seam in ((a, b), (b, c), (c, a))
        ]
        for u, v in carrier.edges:
            self.assertEqual(directed.count((u, v)), 1)
            self.assertEqual(directed.count((v, u)), 1)
        self.assertFalse(self.payload["scope"]["coordinate_geometry_used"])

    def test_global_orientation_reversal_keeps_the_proper_subgroup(self) -> None:
        faces = self.orbit["faces"]
        reversed_faces = tuple(cert.canonical_cycle((a, c, b)) for a, b, c in faces)
        proper_reversed = {
            permutation
            for permutation in self.orbit["automorphisms"]
            if cert.orientation_sign(permutation, reversed_faces) == 1
        }
        self.assertEqual(proper_reversed, set(self.orbit["proper"]))

    def test_proper_group_action_is_simply_transitive(self) -> None:
        seams = set(self.orbit["seams"])
        proper = self.orbit["proper"]
        self.assertEqual(len(seams), 60)
        self.assertEqual(len(proper), 60)
        for source in seams:
            images = {
                (permutation[source[0]], permutation[source[1]])
                for permutation in proper
            }
            self.assertEqual(images, seams)
            self.assertEqual(len(images), len(proper))
        block = self.payload["directed_seam_primitive_orbit"]
        self.assertTrue(block["simply_transitive"])
        self.assertEqual(block["representative_stabilizer_order"], 1)
        self.assertEqual(block["unique_transporter_checks"], 3600)

    def test_628_microhistory_reaches_exact_balanced_shell(self) -> None:
        for left in range(-17, 18):
            for right in range(-17, 18):
                terminal, history = cert.settle_pair_via_628(left, right)
                lower, upper = cert.nearest_balanced_pair(left, right)
                self.assertIn(terminal, ((lower, upper), (upper, lower)))
                self.assertEqual(len(history), abs(left - right) // 2)
                self.assertEqual(sum(terminal), left + right)

    def test_directed_expectation_is_pair_average(self) -> None:
        for left in range(-17, 18):
            for right in range(-17, 18):
                self.assertEqual(
                    cert.directed_pair_expectation(left, right),
                    cert.pair_average(left, right),
                )

    def test_directed_s1_primitive_is_not_preaveraged(self) -> None:
        u, v = self.orbit["seams"][0]
        source_basis = tuple(
            1 if index == u else 0 for index in range(cert.PORT_COUNT)
        )
        forward = cert.directed_balanced_update(source_basis, (u, v))
        reverse = cert.directed_balanced_update(source_basis, (v, u))
        self.assertNotEqual(forward, reverse)
        expected = cert.expected_directed_update(source_basis, (u, v))
        self.assertEqual(
            tuple(
                (Fraction(a) + Fraction(b)) / 2
                for a, b in zip(forward, reverse)
            ),
            expected,
        )

    def test_odd_total_boundary_is_not_hidden(self) -> None:
        forward = cert.directed_balanced_pair(1, 0)
        reverse = cert.opposite_directed_balanced_pair(1, 0)
        self.assertEqual({forward, reverse}, {(0, 1), (1, 0)})
        self.assertNotEqual(forward[0], forward[1])
        self.assertNotEqual(reverse[0], reverse[1])
        self.assertEqual(
            cert.directed_pair_expectation(1, 0),
            (Fraction(1, 2), Fraction(1, 2)),
        )
        self.assertEqual(
            self.payload["verdict"]["pathwise_odd_total_agreement"],
            "not attained",
        )
        self.assertEqual(
            self.payload["verdict"]["opposite_shell_placement_as_628_descent"],
            "not generally attained",
        )

    def test_uniform_s1_channel_is_exactly_i_minus_laplacian_over_60(self) -> None:
        block = self.payload["uniform_s1_channel"]
        channel = tuple(
            tuple(Fraction(entry) for entry in row) for row in block["exact_matrix"]
        )
        laplacian = tuple(
            tuple(Fraction(entry) for entry in row)
            for row in block["exact_laplacian"]
        )
        expected = cert.matrix_add(
            cert.identity_matrix(12),
            cert.matrix_scale(Fraction(-1, 60), laplacian),
        )
        self.assertEqual(channel, expected)
        self.assertEqual(block["channel_identity"], "R = I - Delta_ico/60")
        self.assertEqual(block["projected_undirected_probability"], "1/30")
        self.assertEqual(block["primitive_basis_checks"], 720)
        self.assertEqual(block["opposite_direction_pair_average_checks"], 60)

    def test_incomplete_directed_orbit_fails_closed(self) -> None:
        with self.assertRaises(cert.CertificateError) as caught:
            cert.s1_channel_certificate(
                self.orbit["carrier"].edges,
                self.orbit["seams"][:-1],
            )
        self.assertEqual(caught.exception.code, "DIRECTED_ORBIT_SUPPORT")

    def test_nonmanifold_incidence_fails_orientation(self) -> None:
        adjacency = [set(row) for row in self.orbit["carrier"].adjacency]
        u, v = self.orbit["carrier"].edges[0]
        adjacency[u].remove(v)
        adjacency[v].remove(u)
        with self.assertRaises(cert.CertificateError):
            cert.coherent_orientation(adjacency)

    def test_manifest_round_trip_is_byte_identical_and_fail_closed(self) -> None:
        stored = cert.load_json(cert.MANIFEST_PATH)
        verification = cert.verify_manifest(stored)
        self.assertEqual(verification["status"], "PASS")
        self.assertTrue(verification["producer_replay"])
        self.assertFalse(verification["independent_implementation"])
        self.assertEqual(
            stored,
            {**self.payload, "manifest_sha256": "sha256:" + cert.sha256_json(self.payload)},
        )
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["certify", "--output", str(first)]), 0)
            self.assertEqual(cert.main(["certify", "--output", str(second)]), 0)
            self.assertEqual(first.read_bytes(), cert.MANIFEST_PATH.read_bytes())
            self.assertEqual(first.read_bytes(), second.read_bytes())
            manifest = cert.load_json(first)
            self.assertEqual(cert.verify_manifest(manifest)["status"], "PASS")

            tampered = copy.deepcopy(manifest)
            tampered["uniform_s1_channel"]["channel_identity"] = "I"
            bad = Path(scratch) / "tampered.json"
            bad.write_text(json.dumps(tampered), encoding="utf-8")
            self.assertEqual(cert.main(["verify", "--manifest", str(bad)]), 1)

    def test_full_universe_and_physical_claims_remain_open(self) -> None:
        verdict = self.payload["verdict"]
        self.assertEqual(verdict["refinement_semigroup_naturality"], "open")
        self.assertEqual(
            verdict["physical_clock_and_laboratory_attachment"],
            "open",
        )
        self.assertEqual(
            verdict["directed_orbit_schedule_as_physical_selection"],
            "open",
        )
        self.assertEqual(
            verdict["full_self_readback_and_universe_selection"],
            "open",
        )


if __name__ == "__main__":
    unittest.main()
