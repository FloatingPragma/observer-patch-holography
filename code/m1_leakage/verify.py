"""Strict custody and producer-free native, noise, correction and budget replay."""

import argparse
import hashlib
import json
from pathlib import Path
from m1_noisy_records import verify as parent
from .format import keys, need
from .interface_check import verify as verify_interfaces
from .noise_check import verify as verify_noise
from .correction_check import verify as verify_correction
from .scaling_check import verify as verify_scaling
from .composition_check import verify as verify_composition
import pin_scope


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'format.py', 'interfaces.py', 'interface_check.py',
       'noise.py', 'noise_check.py', 'correction.py', 'correction_check.py',
       'composition.py', 'composition_check.py',
       'scaling.py', 'scaling_check.py', 'build.py', 'verify.py',
       'test_leakage.py', 'README.md', 'CONTRACT.md')
SOURCES = sorted(set(parent.SOURCES+['code/m1_noisy_records/receipt.json']
                    +['code/m1_leakage/'+p for p in OWN]
                    +['docs/research/LEAKAGE_SOURCE_REFINEMENT.md', '.github/workflows/m1-leakage.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-NATIVE-LEAKAGE-REDUCTION', 'OPH-SOURCE-GENERAL-NOISE-REFINEMENT')


def load(path):
    def pairs(items):
        result = {}
        for key, value in items:
            need(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    def reject(token):
        raise ValueError('non-finite JSON token '+token)
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=pairs, parse_constant=reject)


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claims():
    rows = load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        selected = [row for row in rows if row['claim_id'] == name]
        need(len(selected) == 1, 'unique registered leakage claim')
        data = json.dumps(selected[0], sort_keys=True, separators=(',', ':'), ensure_ascii=True)
        result[name] = hashlib.sha256(data.encode('ascii')).hexdigest()
    return result


def verify_evidence(evidence):
    keys(evidence, 'interfaces noise correction scaling composition')
    verify_interfaces(evidence['interfaces'])
    verify_noise(evidence['noise'])
    verify_correction(evidence['correction'])
    verify_scaling(evidence['scaling'])
    verify_composition(evidence['composition'])


def verify(path=HERE/'receipt.json'):
    packet = load(path)
    keys(packet, 'schema sources claims evidence')
    need(packet['schema'] == 'oph-native-leakage-v1', 'receipt schema')
    need(pin_scope.bound_pins(packet['sources']) == pins(), 'source custody')
    need(packet['claims'] == claims(), 'registered claim custody')
    parent.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified native leakage reduction, general noise, complete correction and paid refinement')
