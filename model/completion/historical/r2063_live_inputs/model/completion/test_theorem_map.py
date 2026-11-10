"""Custody checks for the focused source-to-theorem review, not proof tests."""
import hashlib
import json
from pathlib import Path
import unittest

from .family import binding


ROOT = Path(__file__).resolve().parents[2]


class TheoremMapCustodyTests(unittest.TestCase):
    def test_reviewed_sources_and_canonical_navigation(self):
        packet = json.loads((ROOT / "model/completion/theorem_map.json").read_text())
        self.assertEqual(packet["campaign"], binding())
        self.assertFalse(packet["full_gate_accepted"])
        ids = {row["id"] for row in json.loads((ROOT / packet["canonical_registry"]).read_text())["rows"]}
        for source in packet["sources"].values():
            path = (ROOT / source["path"]).resolve()
            self.assertTrue(path.is_relative_to(ROOT))
            data = path.read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), source["sha256"])
            for anchor in source["reviewed_anchors"]:
                self.assertIn(anchor, data.decode())
        for row in packet["mappings"]:
            self.assertTrue(set(row["sources"]) <= set(packet["sources"]))
            self.assertTrue(set(row["canonical_premises"]) <= ids)
            self.assertTrue(row["remaining"])
            for path in row.get("evidence", []):
                self.assertTrue((ROOT / path).is_file())


if __name__ == "__main__":
    unittest.main()
