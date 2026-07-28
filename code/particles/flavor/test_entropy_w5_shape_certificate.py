#!/usr/bin/env python3
"""Regression and adversarial tests for the entropy W5 shape certificate."""

from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import entropy_w5_shape_certificate as cert  # noqa: E402

Q5 = cert.Q5


class EntropyW5ShapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()

    def test_stored_certificate_matches_rebuild(self) -> None:
        stored = json.loads(cert.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, self.payload)

    def test_verdict_is_excluded_with_declared_mismatch(self) -> None:
        comparison = self.payload["comparison"]
        self.assertEqual(comparison["verdict"], "EXCLUDED")
        self.assertEqual(comparison["closest_pairing_mismatch"], comparison["relative_mismatch_flipped"])
        self.assertGreater(float(comparison["closest_pairing_mismatch"]), 0.1)

    def test_golden_orbit_identities_are_exact(self) -> None:
        golden = cert.golden_orbit()
        triple = sorted(
            (cert.invert_profile((Q5(2), Q5(-1), Q5(-1)))),
            key=lambda v: v.to_float(),
        )
        ratio = (triple[2] - triple[1]) / (triple[1] - triple[0])
        self.assertEqual(ratio, cert.PHI)
        self.assertEqual(golden["sorted_gap_ratio"], "phi = (1 + sqrt5)/2, exactly")
        self.assertFalse(golden["double_eigenvalue"])
        self.assertEqual(
            golden["eigenvalues"],
            ["-5/2 + 0*sqrt5", "5/4 + -3/4*sqrt5", "5/4 + 3/4*sqrt5"],
        )

    def test_prolate_orbit_is_exactly_degenerate(self) -> None:
        prolate = cert.prolate_orbit()
        self.assertTrue(prolate["double_eigenvalue"])
        self.assertEqual(prolate["port_values"]["pole"], "2/3 + 0*sqrt5")
        self.assertEqual(prolate["port_values"]["ring"], "-2/15 + 0*sqrt5")

    def test_cyclic_inversion_roundtrip(self) -> None:
        profile = (Q5(Fraction(3, 7)), Q5(Fraction(-1, 7)), Q5(Fraction(-2, 7)))
        q = cert.invert_profile(profile)
        for got, want in zip(cert.port_profile(q), profile):
            self.assertEqual(got, want)
        self.assertTrue((q[0] + q[1] + q[2]).is_zero())

    def test_mutated_profile_breaks_the_golden_identity(self) -> None:
        mutated = (Q5(2), Q5(Fraction(-6, 7)), Q5(Fraction(-8, 7)))
        q = cert.invert_profile(mutated)
        triple = sorted(q, key=lambda v: v.to_float())
        ratio = (triple[2] - triple[1]) / (triple[1] - triple[0])
        self.assertNotEqual(ratio, cert.PHI)

    def test_every_control_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)

    def test_crossing_uses_the_exact_invariants(self) -> None:
        crossing = self.payload["crossing"]
        self.assertIn("s3n2_prolate", crossing)
        self.assertIn("s4n_golden", crossing)
        r_c = float(crossing["r_c_display"])
        self.assertGreater(r_c, 2.0)
        self.assertLess(r_c, 3.0)

    def test_no_lepton_value_outside_the_comparison_block(self) -> None:
        payload = {k: v for k, v in self.payload.items() if k != "comparison"}
        text = json.dumps(payload)
        for token in ("0.000510", "0.105658", "1.77693", "1.8890", "0.5293"):
            self.assertNotIn(token, text)


if __name__ == "__main__":
    unittest.main()
