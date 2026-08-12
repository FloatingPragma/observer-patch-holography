#!/usr/bin/env python3
"""Focused tests for the three-axiom drift guard."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import check_axiom_consistency as guard  # noqa: E402


class AxiomConsistencyGuardTests(unittest.TestCase):
    def test_registry_declares_exactly_three_axioms(self) -> None:
        errors: list[str] = []
        registry = guard.registry_checks(errors)
        self.assertEqual(errors, [])
        self.assertEqual(registry["core_axiom_count"], 3)
        self.assertEqual([a["id"] for a in registry["axioms"]], ["A1", "A2", "A3"])
        a1, a2, _ = registry["axioms"]
        self.assertIn("commutator-closed", a1["formal_concise"])
        self.assertIn(
            "complete for the declared quotient-visible infinitesimal port response",
            a1["formal_concise"],
        )
        self.assertIn("one projective implementer", a2["formal_concise"])
        self.assertIn("PU(H)", a2["formal_concise"])
        self.assertIn("endogenous", a2["formal_concise"])

    def test_active_tree_is_clean(self) -> None:
        errors: list[str] = []
        guard.scan_surfaces(errors)
        guard.entry_surface_checks(errors)
        self.assertEqual(errors, [], "\n".join(errors))

    def test_injected_stale_five_axiom_block_fails(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.md",
            "The theory uses five axioms.\nAxiom 5 selects the sector.\n",
            errors,
        )
        self.assertEqual(len(errors), 2)
        self.assertIn("five-axiom count", errors[0])
        self.assertIn("Axiom 5 as principle", errors[1])

    def test_cross_line_five_axiom_count_fails(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.tex",
            "The theory uses five\naxioms in its basis.\n",
            errors,
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("cross-line stale token (five-axiom count)", errors[0])

    def test_cross_line_retired_selector_name_fails(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.tex",
            "The proof invokes Minimal\nAdmissible Realization.\n",
            errors,
        )
        self.assertEqual(len(errors), 1)
        self.assertIn("cross-line stale token (retired selector name)", errors[0])

    def test_injected_oph5_and_range_fail(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.tex",
            "Under OPH5 the axioms A1--A5 give recovery.\n",
            errors,
        )
        self.assertEqual(len(errors), 2)

    def test_retired_recovery_and_refinement_clauses_fail(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.tex",
            "The Recoverable Generalized Entropy axiom supplies focusing.\n"
            "The refinement-closure clause of Axiom 3 supplies the RG map.\n"
            "Minimal admissibility selects the light sector.\n",
            errors,
        )
        self.assertEqual(len(errors), 3)

    def test_exact_migration_table_row_is_allowed(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "handoff.md",
            '| Five axioms (also "OPH5") | Three axioms | Retired terms are search aliases. |\n',
            errors,
        )
        self.assertEqual(errors, [])

    def test_mathematical_a5_group_reference_passes(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "group.md",
            "The proper rotation group is $A_5$, the alternating group on "
            "five letters, realized by alternatingGroup (Fin 5) in Lean and "
            "by the a5_closure certificates.\n"
            "The sixty rotations of $A_5$ act on the twelve ports.\n",
            errors,
        )
        self.assertEqual(errors, [], "\n".join(errors))

    def test_inventory_writer_reports_surfaces(self) -> None:
        guard.write_inventory()
        payload = (guard.ROOT / "claims" / "active_surface_inventory.json").read_text(
            encoding="utf-8"
        )
        self.assertIn('"oph.active_surface_inventory.v1"', payload)

    def test_committed_inventory_matches_live_scan_without_writing(self) -> None:
        errors: list[str] = []
        guard.inventory_check(errors)
        self.assertEqual(errors, [], "\n".join(errors))


if __name__ == "__main__":
    unittest.main()
