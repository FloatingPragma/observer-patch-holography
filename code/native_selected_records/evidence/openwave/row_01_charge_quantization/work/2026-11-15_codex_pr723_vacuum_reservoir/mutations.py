"""Integrity tests only; the scientific competing law is in the main packet."""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
from verify import check

HERE = Path(__file__).resolve().parent


def main():
    original = json.loads(gzip.decompress((HERE/'EVIDENCE.json.gz').read_bytes()))
    original_receipt = json.loads((HERE/'RECEIPT.json').read_text())
    names = ['census_current','body_current','probe_current','vacuum_reference',
             'declared_preparation','reader_final','missing_control','calibrated_summary']
    results = []
    for name in names:
        data, receipt = copy.deepcopy(original), copy.deepcopy(original_receipt)
        if name == 'census_current':
            data['censuses'][0]['sectors'][0]['branches'][0][3] += 1
        elif name == 'body_current':
            data['records'][0]['body'][0][3] += 1
        elif name == 'probe_current':
            data['records'][0]['probe'][0][3] += 1
        elif name == 'vacuum_reference':
            data['records'][0]['initial_audit_only'][12] = 1
        elif name == 'declared_preparation':
            data['records'][0]['initial_audit_only'][0] -= 1
        elif name == 'reader_final':
            data['records'][0]['reading']['final_from_records'][0] += 1
        elif name == 'missing_control':
            data['records'].pop()
        elif name == 'calibrated_summary':
            receipt['calibrated_unit'] *= 2
        payload = gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0)
        receipt['evidence_sha256'] = hashlib.sha256(payload).hexdigest()
        with tempfile.TemporaryDirectory(prefix='row01-reservoir-mutant-') as tmp:
            output = Path(tmp)
            (output/'EVIDENCE.json.gz').write_bytes(payload)
            (output/'RECEIPT.json').write_text(json.dumps(receipt))
            try:
                check(output)
                rejected = False
            except (AssertionError, ValueError, KeyError, IndexError):
                rejected = True
        results.append({'mutation':name,'rejected':rejected})
    assert all(r['rejected'] for r in results)
    return {'mutations':results,'passed':True,'same_group':True,
            'scientific_emergence_credit':False}


if __name__ == '__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
