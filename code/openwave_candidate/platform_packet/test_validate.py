"""Mutation checks for the preparation validator, including optimized Python."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('packet_validate', ROOT / 'validate.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class PacketValidationTests(unittest.TestCase):
    def setUp(self):
        self.packet = json.loads((ROOT / 'packet.json').read_text())

    def test_valid_and_optimized(self):
        self.assertTrue(module.validate(self.packet)['preparation_valid'])
        subprocess.run([sys.executable, '-O', str(ROOT / 'validate.py')], check=True,
                       capture_output=True)

    def test_mutations_rejected(self):
        for key, value in [('id', True), ('id', 2), ('id', 1.0),
                           ('literal_gate', 'different'), ('criterion', 'different'),
                           ('candidate_id', 'other'), ('status', 'validated'),
                           ('official_status', 'validated'), ('result', {})]:
            packet = copy.deepcopy(self.packet)
            packet['rows'][0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                module.validate(packet)
        for key, value in [('candidate_id', 'other'), ('official_status_change', 0),
                           ('platform_commit', 'other')]:
            packet = copy.deepcopy(self.packet)
            packet[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                module.validate(packet)
        with self.assertRaises(ValueError):
            module.validate(self.packet, b'{}')

    def test_optimized_rejects_promotion(self):
        code = ("import sys,json;sys.path.insert(0,sys.argv[1]);import validate;"
                "p=json.load(open(sys.argv[1]+'/packet.json'));"
                "p['rows'][0]['id']=True;validate.validate(p)")
        result = subprocess.run([sys.executable, '-O', '-c', code, str(ROOT)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('ValueError', result.stderr)


if __name__ == '__main__':
    unittest.main()
