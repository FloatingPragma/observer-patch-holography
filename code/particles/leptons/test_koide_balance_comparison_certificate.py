#!/usr/bin/env python3
"""Regression and adversarial tests for the Koide balance certificates."""

from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import koide_arithmon_cross_control as cross  # noqa: E402
import koide_balance_comparison_certificate as cert  # noqa: E402


class KoideBalanceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cert.build_payload()

    def test_stored_certificate_matches_rebuild(self) -> None:
        stored = json.loads(cert.OUT_PATH.read_text(encoding="utf-8"))
        self.assertEqual(stored, self.payload)

    def test_conditional_tau_enclosure_and_distance(self) -> None:
        row = self.payload["conditional_tau"]
        lo = Decimal(row["tau_enclosure_mev_outward"][0])
        hi = Decimal(row["tau_enclosure_mev_outward"][1])
        self.assertLess(lo, hi)
        self.assertLess(hi - lo, Decimal("0.001"))
        self.assertLess(Decimal(row["distance_sigma"]), Decimal("1"))
        self.assertTrue(row["measured_central_inside_one_sigma_of_output"])
        self.assertEqual(
            row["row_class"], "conditional_postdiction_with_measured_premise_ancestry"
        )
        self.assertLess(Decimal(row["spurious_root_mev"]), Decimal("105"))

    def test_balance_comparison_is_consistent(self) -> None:
        row = self.payload["balance_comparison"]
        self.assertTrue(row["target_inside_enclosure"])
        self.assertLess(Decimal(row["distance_in_enclosure_half_widths"]), Decimal("1"))

    def test_controls_failed_closed(self) -> None:
        for name, verdict in self.payload["controls"].items():
            self.assertTrue(verdict["expected_failure"], name)
            self.assertTrue(verdict["failed"], name)

    def test_tau_solution_solves_the_balance_exactly(self) -> None:
        with localcontext() as context:
            context.prec = cert.DECIMAL_PRECISION
            m_e = cert.PDG_MEV["electron"][0]
            m_mu = cert.PDG_MEV["muon"][0]
            tau_plus, tau_minus = cert.conditional_tau_roots(m_e, m_mu)
            target = Decimal(2) / Decimal(3)
            for tau in (tau_plus, tau_minus):
                q = cert.koide_q(m_e, m_mu, tau)
                self.assertLess(abs(q - target), Decimal("1e-90"))

    def test_ancestry_is_declared(self) -> None:
        boundary = self.payload["claim_boundary"]
        self.assertIn("measured-target ancestry", boundary)
        self.assertIn("neither lane is a prospective prediction", boundary.lower())


class KoideArithmonCrossTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = cross.build_payload()

    def test_stored_cross_certificate_matches_rebuild(self) -> None:
        stored = json.loads(cross.OUT.read_text(encoding="utf-8"))
        self.assertEqual(stored, self.payload)

    def test_cross_fixture_is_version_pinned_and_target_free(self) -> None:
        fixture = cross.load_fixture()
        self.assertFalse(fixture["experimental_data_consumed"])
        self.assertFalse(fixture["prediction_freeze_claimed"])
        self.assertEqual(
            fixture["freeze"]["arithmon_k7_lean"]["theorem"],
            cross.EXPECTED["lean_theorem"],
        )
        self.assertEqual(
            fixture["version_pinned_arithmon_mass_pair"]["m_mu_over_m_e"]["formula"],
            "27^phi",
        )
        self.assertEqual(
            fixture["version_pinned_arithmon_mass_pair"]["m_tau_over_m_e"]["formula"],
            "3477",
        )

    def test_local_recompute_has_the_pinned_strict_sign(self) -> None:
        row = self.payload["local_recompute"]
        self.assertTrue(row["strictly_below_two_thirds"])
        self.assertLess(Decimal(row["Q_minus_two_thirds"]), Decimal(0))
        self.assertLess(Decimal(row["Q"]), Decimal(2) / Decimal(3))

    def test_joint_exact_realization_is_excluded_only_at_declared_scope(self) -> None:
        verdict = self.payload["logical_verdict"]
        self.assertFalse(verdict["joint_exact_realization_possible"])
        self.assertIn("balanced positive-chamber", verdict["scope"])
        self.assertFalse(self.payload["formal_composition_in_one_kernel"])
        self.assertFalse(self.payload["public_physical_promotion_allowed"])
        self.assertFalse(self.payload["prospective_prediction_freeze"])

    def test_balance_root_mutations_close_the_gap(self) -> None:
        for name, control in self.payload["controls"].items():
            self.assertTrue(control["restores_exact_balance_numerically"], name)
            self.assertLess(
                Decimal(control["abs_Q_minus_two_thirds_after_replacement"]),
                Decimal("1e-110"),
                name,
            )

    def test_cross_boundary_does_not_adjudicate_physics(self) -> None:
        boundary = self.payload["claim_boundary"].lower()
        self.assertIn("compare-only", boundary)
        self.assertIn("not composed in one lean environment", boundary)
        self.assertIn("physical adjudication", boundary)
        self.assertFalse(self.payload["experimental_data_consumed"])


if __name__ == "__main__":
    unittest.main()
