"""A changed scientific leaf or original manifest must fail custody."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('record_replay', Path(__file__).with_name('replay.py'))
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class CustodyTest(unittest.TestCase):
    def test_changed_evidence_and_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            base = root / 'evidence/example'
            base.mkdir(parents=True)
            leaf = base / 'PROOF.md'
            leaf.write_text('fixed source statement\n')
            manifest = base / 'MANIFEST.sha256'
            manifest.write_text(replay.sha(leaf) + '  PROOF.md\n')
            pins = {str(p.relative_to(root)): replay.sha(p) for p in (leaf, manifest)}
            (root / 'SNAPSHOT_MANIFEST.json').write_text(json.dumps(pins))
            (root / 'ORIGINS.json').write_text(json.dumps(
                {'example': {'manifest_sha256': replay.sha(manifest)}}))
            self.assertEqual(replay.verify(root), 2)
            original = leaf.read_bytes()
            leaf.write_text('changed scientific claim\n')
            with self.assertRaisesRegex(ValueError, 'Changed native-record evidence'):
                replay.verify(root)
            leaf.write_bytes(original)
            # Updating the outer snapshot cannot replace the bound origin seal.
            manifest.write_text(hashlib.sha256(b'other').hexdigest() + '  PROOF.md\n')
            pins[str(manifest.relative_to(root))] = replay.sha(manifest)
            (root / 'SNAPSHOT_MANIFEST.json').write_text(json.dumps(pins))
            with self.assertRaisesRegex(ValueError, 'Changed original manifest'):
                replay.verify(root)


if __name__ == '__main__':
    unittest.main()
