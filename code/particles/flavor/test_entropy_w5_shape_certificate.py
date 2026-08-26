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

    def test_global_verdict_keeps_boundary_reentry_viable(self) -> None:
        comparison = self.payload["comparison"]
        self.assertEqual(
            comparison["golden_branch_compare"]["verdict"],
            "GOLDEN_BRANCH_EXCLUDED_COMPARE_ONLY",
        )
        self.assertEqual(
            comparison["global_packet_verdict"],
            "EXACTLY_CLASSIFIED_NOT_GLOBALLY_EXCLUDED",
        )
        self.assertEqual(
            comparison["boundary_reentry_compare"]["verdict"],
            "BOUNDARY_BRANCH_CAN_MATCH_CENTRAL_SHAPE_COMPARE_ONLY",
        )
        self.assertAlmostEqual(
            float(comparison["boundary_reentry_compare"]["target_attached_inferred_r"]),
            5.0613260829,
            places=9,
        )

    def test_global_minimizers_are_classified_but_packet_no_go_is_false(self) -> None:
        boundary = self.payload["exhaustiveness_boundary"]
        self.assertEqual(self.payload["schema"], "oph.entropy_w5_shape_certificate.v3")
        self.assertFalse(boundary["full_critical_orbit_classification_proved"])
        self.assertTrue(boundary["global_minimizer_classification_proved"])
        self.assertFalse(boundary["quartic_packet_globally_excluded"])
        self.assertIn(
            "zero-weight closed-simplex boundary states",
            boundary["viable_routes_not_excluded"],
        )
        self.assertFalse(self.payload["promotion_allowed"])

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

    def test_exact_two_value_table_and_crossing_identities(self) -> None:
        theorem = self.payload["exact_global_quartic_certificate"]
        rows = theorem["two_value_stationary_table"]
        self.assertEqual([row["positive_root_multiplicity"] for row in rows], [1, 2, 3, 4, 5])
        self.assertEqual(rows[0]["S3"], "2*sqrt(15)/15")
        self.assertEqual(rows[1]["S4"], "1/8")
        self.assertTrue(all(all(row["exact_checks"].values()) for row in rows))
        crossing = theorem["crossing"]
        self.assertEqual(crossing["r_c_exact"], "(32*sqrt(15)-20*sqrt(6))/27")
        self.assertEqual(crossing["comparisons"]["G2_minus_G3"], "(r - 4*sqrt(6))/48")
        r_c = float(crossing["r_c_display"])
        self.assertGreater(r_c, 2.0)
        self.assertLess(r_c, 3.0)

    def test_three_root_saddles_are_exact_not_seeded(self) -> None:
        theorem = self.payload["exact_global_quartic_certificate"]
        saddle = theorem["stationarity_and_saddles"]
        self.assertEqual(
            saddle["middle_multiplicity_at_least_two"]["second_variation_factor"],
            "-8*A*B*r",
        )
        self.assertEqual(
            [row["bracket"] for row in saddle["middle_multiplicity_one"]["outer_multiplicity_sign_checks"]],
            ["-3*A", "-4*A - 3*B", "-3*A - 4*B", "-3*B"],
        )
        self.assertTrue(all(saddle["checks"].values()))
        self.assertFalse(self.payload["finite_seed_or_numerical_optimization_used_as_proof"])

    def test_closed_simplex_factors_and_boundary_energies_are_exact(self) -> None:
        closed = self.payload["exact_global_quartic_certificate"]["closed_probability_simplex"]
        factors = [
            row["H_k_minus_H_2_factor"]
            for row in closed["high_amplitude_face_comparison"]["factor_table"]
        ]
        self.assertEqual(
            factors,
            [
                "4*(t_3 - 6)**2*(t_3 + 3)/9",
                "3*(t_4 - 6)**2*(t_4 + 6)**2/64",
                "12*(t_5 - 6)**2*(t_5**2 + 6*t_5 + 18)/125",
                "5*(t_6 - 6)**2*(t_6**2 + 4*t_6 + 12)/36",
            ],
        )
        boundary = closed["boundary_minimizer"]
        self.assertEqual(boundary["G_boundary"], "(r**4 - 480)/(8*r**3)")
        self.assertEqual(boundary["original_energy_E_boundary"], "(r**4 - 480)/48")
        self.assertEqual(boundary["b_values"]["four_values"], "0")
        self.assertEqual(boundary["b_values"]["b_minus"], "3-sqrt(r^2-24)/2")
        self.assertEqual(boundary["x_values"]["four_values"], "-1")
        self.assertEqual(boundary["x_values"]["x_minus"], "2-sqrt(r^2-24)/2")
        self.assertIn("(r - 2*sqrt(15))**2", boundary["G_boundary_minus_G_m1_factor"])
        nonnegative = closed["independent_nonnegative_global_certificate"]
        self.assertEqual(
            nonnegative["exact_identity"],
            "J=sum_{i<j<k} p_i*p_j*p_k*(p_i+p_j+p_k)",
        )
        self.assertTrue(nonnegative["nonnegative"])
        face_checks = closed["face_minimum_structure"]["exact_sign_checks"]
        self.assertTrue(all(face_checks["checks"].values()))
        self.assertTrue(
            all(row["strictly_negative"] for row in face_checks["group_tangent_cases"])
        )
        self.assertTrue(all(closed["checks"].values()))

    def test_boundary_quadrupole_ratio_and_inverse_are_exact(self) -> None:
        quadrupole = self.payload["exact_global_quartic_certificate"]["closed_probability_simplex"]["boundary_quadrupole"]
        self.assertEqual(quadrupole["spectrum_up_to_positive_scale"], ["-2", "1-d", "1+d"])
        self.assertEqual(quadrupole["sorted_gap_ratio_R"], "2*d/(3-d)")
        self.assertEqual(quadrupole["inverse_r_squared"], "15+45*R^2/(R+2)^2")
        self.assertTrue(all(quadrupole["exact_checks"].values()))

    def test_no_lepton_value_outside_the_comparison_block(self) -> None:
        payload = {k: v for k, v in self.payload.items() if k != "comparison"}
        text = json.dumps(payload)
        for token in ("0.000510", "0.105658", "1.77693", "1.8890", "5.0613", "0.5293"):
            self.assertNotIn(token, text)


if __name__ == "__main__":
    unittest.main()
