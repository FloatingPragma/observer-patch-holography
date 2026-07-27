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

    def test_injected_oph5_and_range_fail(self) -> None:
        errors: list[str] = []
        guard.scan_text(
            "injected.tex",
            "Under OPH5 the axioms A1--A5 give recovery.\n",
            errors,
        )
        self.assertEqual(len(errors), 2)

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


if __name__ == "__main__":
    unittest.main()
