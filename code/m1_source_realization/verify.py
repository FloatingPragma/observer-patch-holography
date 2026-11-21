"""Fail-closed envelope, source custody and independent evidence replay."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_source_realization'

import pin_scope

from .check import need, same, verify_evidence

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SOURCES = [f'code/m1_source_realization/{p}' for p in
           ('__init__.py', 'model.py', 'circuits.py', 'topology.py', 'check.py', 'verify.py', 'build.py',
            'test_source.py', 'README.md', 'CONTRACT.md')]
SOURCES += ['claims/axiom_registry.yaml', 'docs/AXIOM_REFERENCE.md',
            'code/source_selection_model/response.json', 'code/source_selection_model/response.py',
            'code/source_selection_model/verify_response.py', 'code/source_selection_model/geometry.py',
            'code/source_selection_model/DERIVATION.md', 'code/source_selection_model/RECORD_GLUING.md',
            'code/a5_closure/port_current_inner_certificate.py',
            'code/m1_operational_clocks/check.py', 'requirements.txt', '.gitattributes']
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ()


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def unique(pairs):
    result = {}
    for k, v in pairs:
        need(k not in result, f'duplicate JSON key: {k}')
        result[k] = v
    return result


def invalid(value):
    raise ValueError(f'nonfinite JSON: {value}')


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique, parse_constant=invalid)


def claim_pins():
    rows = load(ROOT/'claims/claim_registry.yaml')['claims']
    out = {}
    for name in CLAIMS:
        matches = [r for r in rows if r['claim_id'] == name]
        need(len(matches) == 1, f'unique claim: {name}')
        out[name] = hashlib.sha256(json.dumps(matches[0], sort_keys=True, separators=(',', ':'),
                                            ensure_ascii=True).encode('ascii')).hexdigest()
    return out


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    need(type(packet) is dict and packet.keys() == {'schema', 'sources', 'claims', 'evidence'}, 'envelope')
    same(packet['schema'], 'oph-source-realization-v1')
    same(pin_scope.bound_pins(packet['sources']), pins(), 'sources')
    same(packet['claims'], claim_pins(), 'claims')
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified source controls, code selection, support, clock compiler and charged error budget')
