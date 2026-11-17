"""A Python build change must not permit changed scientific evidence."""
from pathlib import Path
import copy
import importlib.util
import json
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('runtime_wave_replay', HERE / 'replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class RuntimeMetadata(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='oph-wave-runtime-')
        self.addCleanup(self.temporary.cleanup)
        self.actual = Path(self.temporary.name) / 'IMPORTS.json'
        self.relative = replay.RUNTIME_IMPORT_RECEIPT
        self.original = HERE / 'workspace' / self.relative
        self.expected = replay.sha(self.original)
        self.raw = self.original.read_bytes()
        self.data = json.loads(self.raw)
        self.runtime = '3.13.0 (portable comparison control) [Clang test]'
        self.portable = self.raw.replace(json.dumps(self.data['python']).encode(),
                                         json.dumps(self.runtime).encode(), 1)

    def compare(self, data, *, relative=None, generated=None, runtime=None):
        self.actual.write_bytes(data)
        return replay.compare_regenerated(relative or self.relative, self.actual,
            self.expected, [self.relative] if generated is None else generated,
            runtime or self.runtime)

    def test_exact_bytes_and_declared_runtime_only(self):
        self.assertIsNone(self.compare(self.raw))
        result = self.compare(self.portable)
        self.assertEqual(result['frozen_value'], self.data['python'])
        self.assertEqual(result['runtime_value'], self.runtime)
        self.assertTrue(result['every_other_byte_matches'])
        self.assertNotEqual(result['frozen_sha256'], result['regenerated_sha256'])

    def test_import_hash_and_dependency_mutations_rejected(self):
        for field in ('files', 'numpy', 'sympy', 'numerical_threads'):
            data = copy.deepcopy(self.data)
            data['python'] = self.runtime
            if field == 'files':
                data[field][next(iter(data[field]))] = '0' * 64
            else:
                data[field] += '-changed'
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.compare((json.dumps(data, indent=2) + '\n').encode())

    def test_unlisted_outputs_and_sources_rejected(self):
        for relative, generated in ((self.relative, []),
                ('openwave-tuning/other/IMPORTS.json', ['openwave-tuning/other/IMPORTS.json']),
                ('openwave-tuning/source.py', [])):
            with self.subTest(path=relative, generated=generated), self.assertRaises(ValueError):
                self.compare(self.portable, relative=relative, generated=generated)

    def test_false_runtime_and_unsupported_version_rejected(self):
        with self.assertRaises(ValueError):
            self.compare(self.portable, runtime='3.13.99 (not this runtime)')
        changed = self.portable.replace(b'3.13.0', b'3.12.0')
        with self.assertRaises(ValueError):
            self.compare(changed, runtime=self.runtime.replace('3.13.0', '3.12.0'))

    def test_missing_or_extra_metadata_rejected(self):
        for mutation in ('missing', 'extra', 'whitespace'):
            data = copy.deepcopy(self.data)
            data['python'] = self.runtime
            if mutation == 'missing':
                del data['python']
            elif mutation == 'extra':
                data['forged_scientific_result'] = True
            raw = (json.dumps(data, indent=4 if mutation == 'whitespace' else 2) + '\n').encode()
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                self.compare(raw)

    def test_symlink_cannot_bypass_exact_hash(self):
        self.actual.symlink_to(self.original)
        with self.assertRaises(ValueError):
            replay.compare_regenerated(self.relative, self.actual, self.expected,
                                       [self.relative], self.runtime)


if __name__ == '__main__':
    unittest.main()
