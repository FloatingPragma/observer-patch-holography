"""Mutation checks for the publication transport, independent of science."""
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "openwave_wave_reproduce", Path(__file__).with_name("reproduce.py")
)
r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(r)


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / 'record.txt').write_bytes(b'sealed\n')
        self.write_manifest(hashlib.sha256(b'sealed\n').hexdigest() + '  record.txt\n')

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self, text):
        path = self.root / 'MANIFEST.sha256'
        path.write_bytes(text.encode())
        self.anchor = {'manifest_sha256': r.digest(path), 'entries': 1}

    def test_exact_source_passes(self):
        self.assertEqual(r.sealed_directory(self.root, self.anchor), 1)

    def test_changed_source_fails(self):
        (self.root / 'record.txt').write_bytes(b'sealed\r\n')
        with self.assertRaisesRegex(ValueError, 'Evidence digest mismatch'):
            r.sealed_directory(self.root, self.anchor)

    def test_extra_source_fails(self):
        (self.root / 'extra.py').write_bytes(b'pass\n')
        with self.assertRaisesRegex(ValueError, 'inventory mismatch'):
            r.sealed_directory(self.root, self.anchor)

    def test_changed_manifest_fails(self):
        (self.root / 'MANIFEST.sha256').write_bytes(b'wrong\n')
        with self.assertRaisesRegex(ValueError, 'Manifest digest mismatch'):
            r.sealed_directory(self.root, self.anchor)

    def test_escape_fails(self):
        self.write_manifest('0' * 64 + '  ../outside.txt\n')
        with self.assertRaisesRegex(ValueError, 'escapes evidence'):
            r.sealed_directory(self.root, self.anchor)

    def test_optimized_public_entry_fails(self):
        env = os.environ.copy()
        env['PYTHONOPTIMIZE'] = '1'
        result = subprocess.run([sys.executable, '-B', str(r.HERE / 'reproduce.py'), 'custody'],
                                env=env, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b'Optimized Python is unsupported', result.stderr)


if __name__ == '__main__':
    unittest.main()
