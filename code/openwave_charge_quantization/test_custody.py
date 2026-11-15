"""The public source admission must reject altered or incomplete packages."""
from pathlib import Path
import importlib.util
import shutil
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('charge_replay', HERE / 'replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class Custody(unittest.TestCase):
    def test_custody_controls(self):
        replay.verify()
        with tempfile.TemporaryDirectory(prefix='oph-charge-controls-') as name:
            destination = Path(name) / 'package'
            shutil.copytree(HERE, destination)
            original = replay.HERE
            replay.HERE = destination
            try:
                replay.verify()
                source = destination / 'workspace/openwave-tuning/row_completion/ROW_01/topology.py'
                data = source.read_bytes()
                source.write_bytes(data + b'\n# altered input\n')
                with self.assertRaises(ValueError): replay.verify()
                source.write_bytes(data)
                source.unlink()
                with self.assertRaises((FileNotFoundError, ValueError)): replay.verify()
                source.write_bytes(data)
                extra = destination / 'workspace/unpinned.py'
                extra.write_text('raise RuntimeError()\n')
                with self.assertRaises(ValueError): replay.verify()
                extra.unlink()
                target = destination / 'outside.py'
                target.write_bytes(data)
                source.unlink()
                source.symlink_to(target)
                with self.assertRaises(ValueError): replay.verify()
                source.unlink()
                source.write_bytes(data)
                manifest = destination / 'PUBLIC_MANIFEST.json'
                raw = manifest.read_bytes()
                manifest.unlink()
                with self.assertRaises(ValueError): replay.verify()
                manifest.write_bytes(raw)
                replay.verify()
            finally:
                replay.HERE = original


if __name__ == '__main__':
    unittest.main()
