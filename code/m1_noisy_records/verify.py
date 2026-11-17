"""Strict custody plus independent channel, fault and resource replay."""

import argparse
import hashlib
import json
from pathlib import Path

from m1_fixed_noise import verify as parent
from .independent import keys, need, verify_evidence
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'circuits.py', 'independent.py', 'archive.py', 'archive_check.py',
       'computation_code.py', 'computation_code_check.py',
       'instruments.py', 'instrument_check.py',
       'build.py', 'verify.py', 'test_records.py', 'README.md', 'CONTRACT.md')
SOURCES = sorted(set(parent.SOURCES+['code/m1_fixed_noise/receipt.json']+
                    ['code/m1_noisy_records/'+p for p in OWN]+
                    ['docs/research/NOISY_SOURCE_RECORDS.md', '.github/workflows/m1-noisy-records.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-NOISY-CONTROL-HISTORIES', 'OPH-SOURCE-LIVE-RECORD-PROTECTION')


def load(path):
    def object_pairs(pairs):
        result = {}
        for key, value in pairs:
            need(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError('non-finite JSON constant '+value)
    return json.loads(Path(path).read_text(encoding='utf-8'),
                      object_pairs_hook=object_pairs, parse_constant=nonfinite)


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claims():
    rows = load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        found = [r for r in rows if r['claim_id'] == name]
        need(len(found) == 1, 'unique registered claim')
        raw = json.dumps(found[0], sort_keys=True, separators=(',', ':'), ensure_ascii=True)
        result[name] = hashlib.sha256(raw.encode('ascii')).hexdigest()
    return result


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources claims evidence')
    need(packet['schema'] == 'oph-noisy-source-records-v1', 'schema')
    need(pin_scope.bound_pins(packet['sources']) == pins(), 'source custody')
    need(packet['claims'] == claims(), 'registered claim custody')
    parent.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified complete adaptive histories, live record faults and charged refinement')
