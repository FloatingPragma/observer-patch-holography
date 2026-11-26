"""Premise-edge validation at the transition to three-digit canonical identifiers."""
import unittest

from tools.check_claim_registry import check_premise_dependencies


class ClaimPremiseTests(unittest.TestCase):
    def setUp(self):
        self.known = {"PR-01", "PR-99", "PR-100", "PR-146", "PR-180"}

    def check_edges(self, consumed, boundary):
        check_premise_dependencies(
            {
                "claim_id": "fixture",
                "premise_dependencies": {
                    "classification": "explicit_edges",
                    "consumed": consumed,
                    "open": [],
                    "boundary": boundary,
                },
            },
            self.known,
        )

    def test_canonical_boundary_and_antecedent(self):
        self.check_edges(["PR-146"], ["PR-180"])
        self.check_edges(["PR-01"], ["PR-99"])
        self.check_edges([], ["PR-100"])

    def test_aliases_unknown_and_overlapping_roles(self):
        for alias in ("PR-1", "PR-099", "PR-0180", "PR-999"):
            with self.subTest(alias=alias), self.assertRaises(SystemExit):
                self.check_edges([], [alias])
        with self.assertRaises(SystemExit):
            self.check_edges(["PR-180"], ["PR-180"])
