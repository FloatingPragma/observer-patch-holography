"""Scoped custody plus producer-free graph, fault, decoder and channel replay."""

import argparse
import hashlib
import json
from pathlib import Path
from m1_noisy_records import verify as parent
from .format import keys, need
import pin_scope
from . import graph_check, circuit_check, decoder_check, bounds_check, export_check

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'format.py', 'graphs.py', 'graph_check.py', 'circuit.py',
       'circuit_check.py', 'decoder.py', 'decoder_check.py', 'bounds.py', 'bounds_check.py',
       'export.py', 'export_check.py', 'verify.py', 'build.py', 'test_archive.py',
       'README.md', 'CONTRACT.md')
SOURCES = sorted(set(parent.SOURCES+['code/m1_noisy_records/receipt.json']+
                     ['code/m1_expander_archive/'+p for p in OWN]+
                     ['docs/research/FIXED_STRENGTH_PUBLIC_RECORDS.md', '.github/workflows/m1-expander-archive.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-EXPANDER-LIVE-ARCHIVE', 'OPH-SOURCE-FIXED-STRENGTH-PUBLIC-HISTORY')


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def nonfinite(token):
        raise ValueError('non-finite JSON token '+token)
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs, parse_constant=nonfinite)


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claims():
    rows = load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        selected = [row for row in rows if row['claim_id'] == name]
        need(len(selected) == 1, 'unique registered archive claim')
        raw = json.dumps(selected[0], sort_keys=True, separators=(',', ':'), ensure_ascii=True)
        result[name] = hashlib.sha256(raw.encode('ascii')).hexdigest()
    return result


def verify_evidence(row):
    keys(row, 'graphs circuit decoder bounds export')
    graph_check.verify(row['graphs'])
    circuit_check.verify(row['circuit'])
    decoder_check.verify(row['decoder'])
    bounds_check.verify(row['bounds'])
    export_check.verify(row['export'])


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources claims evidence')
    need(packet['schema'] == 'oph-fixed-strength-public-records-v1', 'receipt schema')
    need(pin_scope.bound_pins(packet['sources']) == pins(), 'source custody')
    need(packet['claims'] == claims(), 'scoped claim custody')
    parent.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', type=Path, nargs='?', default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified fixed-strength live archive, actual-code decoder and complete export')
