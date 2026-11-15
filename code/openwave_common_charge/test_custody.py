"""The public source admission must reject altered or incomplete packages."""
from pathlib import Path
import importlib.util
import shutil
import tempfile
import unittest
from copy import deepcopy

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('charge_replay', HERE / 'replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class Custody(unittest.TestCase):
    def test_identical_source_path_controls(self):
        value = {'source_files_read': {'copied/source.py': 'abc'}, 'passed': True}
        mapping = {'copied/source.py': 'archive/source.py'}
        pins = {'copied/source.py': 'abc', 'archive/source.py': 'abc'}
        normalized, changes = replay.normalize_source_paths(value, mapping, pins)
        self.assertEqual(normalized['source_files_read'], {'archive/source.py': 'abc'})
        self.assertEqual(len(changes), 1)
        for bad in ({'copied/source.py': 'def', 'archive/source.py': 'abc'},
                    {'copied/source.py': 'abc'}, {}):
            with self.assertRaises(ValueError):
                replay.normalize_source_paths(value, mapping, bad)

    def test_only_declared_roundoff_is_admitted(self):
        expected = {'Gram_residual': 1e-14, 'mass': .23,
                    'bound_exact': '1/1000', 'count': 36, 'passed': True}
        actual = deepcopy(expected)
        actual['Gram_residual'] = 2e-14
        self.assertEqual(len(replay.compare_diagnostics(expected, actual, {'Gram_residual'})), 1)
        for key, value in [('mass', .23000000000001), ('bound_exact', '1/999'),
                           ('count', True), ('passed', False),
                           ('Gram_residual', float('nan')), ('Gram_residual', 1e-10)]:
            actual = deepcopy(expected)
            actual[key] = value
            with self.assertRaises(ValueError):
                replay.compare_diagnostics(expected, actual, {'Gram_residual'})

    def test_custody_controls(self):
        replay.verify()
        with tempfile.TemporaryDirectory(prefix='oph-charge-controls-') as name:
            destination = Path(name) / 'package'
            shutil.copytree(HERE, destination)
            original = replay.HERE
            replay.HERE = destination
            try:
                replay.verify()
                source = destination / 'workspace/openwave-tuning/next/CHAIN5/common_charge_topology/bridge.py'
                data = source.read_bytes()
                source.write_bytes(data + b'\n# altered input\n')
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                source.write_bytes(data)
                source.unlink()
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                source.write_bytes(data)
                extra = destination / 'workspace/unpinned.py'
                extra.write_text('raise RuntimeError()\n')
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                extra.unlink()
                target = destination / 'outside.py'
                target.write_bytes(data)
                source.unlink()
                source.symlink_to(target)
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                source.unlink()
                source.write_bytes(data)
                manifest = destination / 'PUBLIC_MANIFEST.json'
                raw = manifest.read_bytes()
                manifest.unlink()
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                manifest.write_bytes(raw)
                replay.verify()
            finally:
                replay.HERE = original


if __name__ == '__main__':
    unittest.main()
