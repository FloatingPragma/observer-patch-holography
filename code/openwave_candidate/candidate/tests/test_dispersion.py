"""Exact structural and direct-symbol checks of the definition-only recompute."""
import importlib.util
from pathlib import Path
import unittest

MODULE = Path(__file__).resolve().parents[1] / "oph_openwave/dispersion.py"
SPEC = importlib.util.spec_from_file_location("definition_only_dispersion", MODULE)
dispersion = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dispersion)


class DispersionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = dispersion.run()

    def test_exact_identities_and_direct_bounds(self):
        self.assertTrue(self.payload["checks"])
        self.assertTrue(all(value is True for value in self.payload["checks"].values()))

    def test_broken_controls_rejected(self):
        self.assertEqual(len(self.payload["mutations"]), 2)
        self.assertTrue(all(value is True for value in self.payload["mutations"].values()))

    def test_direct_check_coverage_and_precision(self):
        numeric = self.payload["outputs"]["direct_symbol"]
        self.assertGreaterEqual(numeric["decimal_precision"], 60)
        self.assertEqual(len(numeric["rows"]), 40)
        self.assertEqual(len({row["direction"] for row in numeric["rows"]}), 4)
        self.assertEqual(len({row["k"] for row in numeric["rows"]}), 5)

    def test_scope_is_definition_only(self):
        self.assertIs(self.payload["full_criterion_met"], False)
        self.assertIs(self.payload["target_values_not_provided"], True)
        self.assertEqual(len(self.payload["scientific_input_files_read"]), 1)
        self.assertIn("No physical identification", self.payload["scope"])


if __name__ == "__main__":
    unittest.main()
