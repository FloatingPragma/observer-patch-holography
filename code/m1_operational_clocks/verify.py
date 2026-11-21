"""Fail-closed receipt verification with independent reconstruction and source custody."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_operational_clocks'

from .check import need, same, verify_evidence
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCES = [f'code/m1_operational_clocks/{p}' for p in
           ('__init__.py', 'model.py', 'check.py', 'verify.py', 'build.py', 'test_clocks.py', 'README.md', 'CONTRACT.md')]
SOURCES += ['code/m1_quantum_transport/tetra_model.py', 'code/m1_quantum_transport/tetra_check.py',
            'code/m1_quantum_transport/DERIVATION.md', 'requirements.txt', '.gitattributes']
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-M1-MINIMAL-TETRAHEDRAL-TRANSPORT',)


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def unique(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, f'duplicate JSON key: {key}')
        result[key] = value
    return result


def nonfinite(value):
    raise ValueError(f'nonfinite JSON constant: {value}')


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique, parse_constant=nonfinite)


def claim_pins():
    claims = load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        rows = [row for row in claims if row['claim_id'] == name]
        need(len(rows) == 1, f'claim must occur exactly once: {name}')
        raw = json.dumps(rows[0], sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode('ascii')
        result[name] = hashlib.sha256(raw).hexdigest()
    return result


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    need(type(packet) is dict and packet.keys() == {'schema', 'sources', 'claims', 'evidence'}, 'invalid envelope')
    same(packet['schema'], 'oph-m1-operational-clocks-v1')
    same(pin_scope.bound_pins(packet['sources']), pins(), 'sources')
    same(packet['claims'], claim_pins(), 'claims')
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    args = parser.parse_args()
    verify(args.path)
    print('Verified full massive spectra, all-channel dynamics, native records and finite clock readouts')
