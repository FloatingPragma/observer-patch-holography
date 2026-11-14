"""Adversarial custody tests use disposable minimal package fixtures."""
from pathlib import Path
import hashlib
import importlib.util
import json
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('public_replay', HERE / 'replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


class Custody(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.original = replay.HERE
        replay.HERE = self.root
        (self.root / 'workspace/model').mkdir(parents=True)
        (self.root / 'workspace/model/source.py').write_bytes(b'exact source\n')
        snapshot = {'files_sha256': {'model/source.py': sha(b'exact source\n')}}
        (self.root / 'SNAPSHOT_MANIFEST.json').write_text(json.dumps(snapshot))
        public = {'files_sha256': {'SNAPSHOT_MANIFEST.json': sha((self.root / 'SNAPSHOT_MANIFEST.json').read_bytes())}}
        (self.root / 'PUBLIC_MANIFEST.json').write_text(json.dumps(public))

    def tearDown(self):
        replay.HERE = self.original
        self.temporary.cleanup()

    def test_exact_fixture_passes(self):
        replay.verify()

    def test_one_byte_source_change_rejected(self):
        (self.root / 'workspace/model/source.py').write_bytes(b'Exact source\n')
        with self.assertRaises(ValueError):
            replay.verify()

    def test_manifest_repin_rejected(self):
        (self.root / 'SNAPSHOT_MANIFEST.json').write_text('{}')
        with self.assertRaises(ValueError):
            replay.verify()

    def test_unmanifested_source_rejected(self):
        (self.root / 'workspace/model/injected.py').write_text('pass')
        with self.assertRaises(ValueError):
            replay.verify()

    def test_exact_byte_symlink_rejected(self):
        source = self.root / 'workspace/model/source.py'
        external = self.root / 'external.py'
        external.write_bytes(source.read_bytes())
        source.unlink()
        source.symlink_to(external)
        with self.assertRaises(ValueError):
            replay.verify()


if __name__ == '__main__':
    unittest.main()
