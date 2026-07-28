#!/usr/bin/env python3
"""Regression and adversarial tests for the Koide balance certificate."""

from __future__ import annotations

import json
import sys
import unittest
from decimal import Decimal, localcontext
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

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


if __name__ == "__main__":
    unittest.main()
