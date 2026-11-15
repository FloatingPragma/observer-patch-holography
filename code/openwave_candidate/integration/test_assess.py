"""Persisted integrity and input controls for all-row local assessment."""
import copy
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import json

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('assessment_dispatch', HERE/'assess.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)

class Integrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = a.assess()

    def rejected(self, receipt):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'receipt.json'; path.write_text(json.dumps(receipt))
            result=subprocess.run([sys.executable, '-O', str(HERE/'assess.py'), '--verify', str(path)],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)

    def test_evidence_references_must_resolve(self):
        a.validate_evidence_links(self.receipt)
        r=copy.deepcopy(self.receipt)
        r['rows'][3]['measurement']['evidence_key']='missing.scientific.output'
        with self.assertRaises(ValueError): a.validate_evidence_links(r)

    def test_missing_sources(self):
        r=copy.deepcopy(self.receipt);r['sources']={};self.rejected(r)

    def test_missing_row(self):
        r=copy.deepcopy(self.receipt);r['rows'].pop();self.rejected(r)

    def test_duplicate_row(self):
        r=copy.deepcopy(self.receipt);r['rows'][1]=r['rows'][0];self.rejected(r)

    def test_source_invalidation(self):
        r=copy.deepcopy(self.receipt);r['sources'][next(iter(r['sources']))]='0'*64;self.rejected(r)

    def test_measurement_tampering(self):
        r=copy.deepcopy(self.receipt);r['rows'][0]['measurement']['gauss_error']=1;self.rejected(r)

    def test_boolean_and_nan_receipts(self):
        for value in (True,float('nan')):
            r=copy.deepcopy(self.receipt);r['rows'][0]['id']=value;self.rejected(r)

if __name__=='__main__':unittest.main()
