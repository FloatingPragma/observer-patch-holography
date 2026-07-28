#!/usr/bin/env python3
"""Regression and adversarial tests for the GitHub #628 mechanism certificate."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(MODULE_DIR))

import record_counting_mechanism_certificate as cert  # noqa: E402


class RecordCountingMechanismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()
        cls.carrier, cls.rotations, cls.edges, _ = cert.load_carrier()

    def test_reference_manifest_is_exactly_recomputable(self) -> None:
        stored = cert.load_json(cert.MANIFEST_PATH)
        body = {k: v for k, v in stored.items() if k != "manifest_sha256"}
        self.assertEqual(body, self.payload)
        self.assertEqual(stored["manifest_sha256"], "sha256:" + cert.sha256_json(body))
        self.assertEqual(cert.verify_stored()["status"], "PASS")

    def test_rerun_writes_byte_identical_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as scratch:
            first = Path(scratch) / "first.json"
            second = Path(scratch) / "second.json"
            self.assertEqual(cert.main(["--output", str(first)]), 0)
            self.assertEqual(cert.main(["--output", str(second)]), 0)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_counting_map_is_atomic_and_integer(self) -> None:
        counts = cert.counting_map([[1, 1, -1], [], [1]] + [[]] * 9)
        self.assertEqual(counts[:3], (1, 0, 1))
        with self.assertRaises(cert.CertificateError):
            cert.counting_map([[1], [2]] + [[]] * 10)
        with self.assertRaises(cert.CertificateError) as caught:
            cert.counting_map([[]] * 11)
        self.assertEqual(caught.exception.code, "REGISTER_ARITY")

    def test_separation_witnesses_are_exact(self) -> None:
        separation = self.payload["cost_separation"]
        self.assertEqual(
            separation["adjacent_dipole"],
            {
                "settle_cost": 1,
                "seam_quadratic": 12,
                "load_square": 2,
                "a3_hessian_12I_value": 24,
            },
        )
        self.assertEqual(
            separation["distance_two_dipole"],
            {
                "settle_cost": 0,
                "seam_quadratic": 10,
                "load_square": 2,
                "a3_hessian_12I_value": 24,
            },
        )
        self.assertEqual(
            separation["quadratic_candidate_disposition"],
            "retained_as_distinct_typed_implementation",
        )

    def test_candidate_dispositions(self) -> None:
        dispositions = self.payload["candidate_dispositions"]
        self.assertEqual(
            dispositions["half_unit_readback"],
            "retained_as_same_complete_mechanism_under_unit_rescaling",
        )
        lift = self.payload["half_unit_complete_lift"]
        self.assertEqual(lift["scale"], "1/2")
        self.assertEqual(lift["displayed_threshold"], "1/1")
        self.assertGreater(lift["naturality_basis_checks"], 0)
        self.assertGreater(lift["minimum_cost_checks"], 0)
        self.assertEqual(self.payload["bounded_exit"], "exact_named_realization")

    def test_consensus_witness_replays(self) -> None:
        rows = self.payload["theorems"]["fiber_obstruction"]["rows"]
        witness = next(r for r in rows if r.get("consensus"))
        self.assertEqual(witness["witness_schedule_moves"], 18)
        self.assertTrue(witness["seam_quadratic_increase_steps"])
        state = tuple([12] + [0] * 11)
        for move in witness["witness_schedule"]:
            move = tuple(move)
            self.assertIn(move, cert.admissible_moves(state, self.edges))
            state = cert.apply_move(state, move)
        self.assertEqual(max(state), min(state))
        for row in rows:
            if not row.get("consensus"):
                self.assertNotEqual(row["total"] % 12, 0)

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)
        self.assertEqual(
            self.payload["controls"]["wrong_reference"]["code"],
            "A3_REFERENCE",
        )
        self.assertEqual(
            self.payload["controls"]["alternative_weight"]["code"],
            "A3_WEIGHT",
        )
        self.assertEqual(
            self.payload["controls"]["higher_order_objective"]["code"],
            "A3_OBJECTIVE",
        )

    def test_cost_naturality_is_explicit(self) -> None:
        clauses = self.payload["a2_clauses"]
        rotation = clauses["rotation_naturality"]["cost_move_graph_isomorphism"]
        refinement = clauses["refinement_naturality"][
            "cost_move_graph_isomorphism"
        ]
        self.assertEqual(rotation["graph_isomorphisms"], 60)
        self.assertEqual(refinement["graph_isomorphisms"], 3)
        self.assertEqual(rotation["shortest_path_witnesses"], 60)
        self.assertEqual(refinement["shortest_path_witnesses"], 3)

    def test_strict_descent_on_a_fresh_probe(self) -> None:
        probe = tuple([4, -4] + [0] * 10)
        v0 = cert.descent_potential(probe)
        for move in cert.admissible_moves(probe, self.edges):
            moved = cert.apply_move(probe, move)
            v1 = cert.descent_potential(moved)
            difference = probe[move[0]] - probe[move[1]]
            self.assertEqual(v1 - v0, -2 * (difference - 1))
            self.assertLessEqual(v1, v0 - 2)

    def test_seam_quadratic_is_not_the_lyapunov_function(self) -> None:
        probe = (-2, -2, -2, -2, -2, -2, -2, 0, -2, -2, 1, 2)
        move = (7, 5)
        self.assertIn(move, cert.admissible_moves(probe, self.edges))
        moved = cert.apply_move(probe, move)
        self.assertGreater(
            cert.quadratic(moved, self.edges),
            cert.quadratic(probe, self.edges),
        )
        self.assertLess(
            cert.descent_potential(moved),
            cert.descent_potential(probe),
        )

    def test_record_repair_commutes_with_counting(self) -> None:
        registers = cert.full_pile_registers()
        before = cert.counting_map(registers)
        move = cert.admissible_moves(before, self.edges)[0]
        after_registers = cert.apply_record_move(registers, move, self.edges)
        self.assertEqual(
            cert.counting_map(after_registers),
            cert.apply_move(before, move),
        )

    def test_settled_characterization(self) -> None:
        lipschitz = tuple([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1])
        self.assertTrue(cert.is_settled(lipschitz, self.edges))
        self.assertEqual(cert.admissible_moves(lipschitz, self.edges), [])

    def test_doctored_settle_search_bound_fails_closed(self) -> None:
        counts = tuple([6, -6] + [0] * 10)
        with self.assertRaises(cert.CertificateError):
            cert.settle_cost_exact(counts, self.edges, limit=3)


if __name__ == "__main__":
    unittest.main()
