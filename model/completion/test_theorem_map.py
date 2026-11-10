"""Custody checks for the focused source-to-theorem review, not proof tests."""
import hashlib
import json
from pathlib import Path
import unittest

from .family import binding
from . import current_theorem_navigation as navigation


ROOT = Path(__file__).resolve().parents[2]


class TheoremMapCustodyTests(unittest.TestCase):
    def test_reviewed_sources_and_canonical_navigation(self):
        packet = navigation.validate()
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


def test_current_projection_cannot_promote_science_or_rewrite_history():
    import copy
    import pytest
    original = navigation.build()
    for damage in ('acceptance', 'classification', 'scope', 'retained', 'source'):
        packet = copy.deepcopy(original)
        if damage == 'acceptance':
            packet['reviewed_map_with_current_navigation']['full_gate_accepted'] = True
        elif damage == 'classification':
            packet['reviewed_map_with_current_navigation']['mappings'][0]['classification'] = 'proved_physical'
        elif damage == 'scope':
            packet['reviewed_map_with_current_navigation']['mappings'][0]['remaining'] = []
        elif damage == 'retained':
            packet['retained_review']['sha256'] = 'f' * 64
        else:
            packet['reviewed_map_with_current_navigation']['sources']['collision']['sha256'] = 'f' * 64
        with pytest.raises(ValueError):
            navigation.validate(packet)
