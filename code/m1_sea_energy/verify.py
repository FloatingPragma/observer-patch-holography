"""Source custody and producer-free semantic replay of the filled-source result."""

import argparse
import hashlib
import json
from pathlib import Path
from m1_expander_archive import verify as protected
from m1_operational_clocks import verify as clocks
from . import spectrum_check, preparation_check, clock_check, bounds_check
from .format import keys, need
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'format.py', 'spectrum.py', 'spectrum_check.py', 'preparation.py',
       'preparation_check.py', 'clock.py', 'clock_check.py', 'bounds.py', 'bounds_check.py',
       'verify.py', 'build.py', 'test_sea.py', 'README.md', 'CONTRACT.md')
SOURCES = sorted(set(protected.SOURCES+clocks.SOURCES+
                    ['code/m1_expander_archive/receipt.json', 'code/m1_operational_clocks/receipt.json']+
                    ['code/m1_sea_energy/'+name for name in OWN]+
                    ['docs/research/SOURCE_SEA_ENERGY.md', '.github/workflows/m1-sea-energy.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-SEA-GENERATOR', 'OPH-SOURCE-LOCAL-SEA-LIMIT', 'OPH-SOURCE-SEA-CLOCK')


def load(path):
    def pairs(items):
        out = {}
        for k, v in items:
            need(k not in out, 'duplicate JSON key')
            out[k] = v
        return out
    def nonfinite(value):
        raise ValueError('nonfinite JSON token '+value)
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs, parse_constant=nonfinite)


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claims():
    rows = load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        selected = [row for row in rows if row['claim_id'] == name]
        need(len(selected) == 1, 'unique registered sea claim')
        raw = json.dumps(selected[0], sort_keys=True, separators=(',', ':'), ensure_ascii=True)
        result[name] = hashlib.sha256(raw.encode('ascii')).hexdigest()
    return result


def verify_evidence(evidence):
    keys(evidence, 'spectrum preparation clock bounds')
    spectrum_check.verify(evidence['spectrum'])
    preparation_check.verify(evidence['preparation'])
    clock_check.verify(evidence['clock'])
    bounds_check.verify(evidence['bounds'])


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources claims evidence')
    need(packet['schema'] == 'oph-source-sea-energy-v1', 'receipt schema')
    need(pin_scope.bound_pins(packet['sources']) == pins(), 'source custody')
    need(packet['claims'] == claims(), 'claim custody')
    protected.verify()
    clocks.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified filled source generator, local vacuum, complete reads and charged protection')
