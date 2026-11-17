"""Strict source custody and independent recovery/noise/operational replay."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_fixed_noise'

from m1_fermionic_source import verify as parent
from m1_fermionic_source.check import exact, keys, need
from .algebra_check import check_accounting, check_decoder, check_interfaces, check_memory, check_noise
from .recovery_check import verify_recovery
from .scaling_check import verify as verify_scaling
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN = ('__init__.py', 'algebra.py', 'algebra_check.py', 'recovery.py', 'recovery_check.py',
       'executor.py', 'instrument_check.py', 'noise.py', 'scaling.py', 'scaling_check.py', 'build.py', 'verify.py',
       'test_fixed_noise.py', 'README.md', 'CONTRACT.md')
SOURCES = sorted(set(parent.SOURCES+['code/m1_fermionic_source/receipt.json']
                    +['code/m1_fixed_noise/'+p for p in OWN]
                    +['docs/research/FIXED_RATE_SOURCE_READS.md', '.github/workflows/m1-fixed-noise.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-FIXED-RATE-PROTECTED-READS', 'OPH-SOURCE-RAW-DECODED-NOISE-SEPARATION')


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claims():
    registry = parent.parent.load(ROOT/'claims/claim_registry.yaml')['claims']
    result = {}
    for name in CLAIMS:
        found = [r for r in registry if r['claim_id'] == name]
        need(len(found) == 1, 'unique registered fixed-rate claim')
        result[name] = hashlib.sha256(json.dumps(found[0], sort_keys=True, separators=(',', ':'),
                                                 ensure_ascii=True).encode('ascii')).hexdigest()
    return result


def verify_evidence(row):
    keys(row, 'decoder interfaces recovery memory continuous_noise accounting scaling')
    code, maps = check_decoder(row['decoder'])
    check_interfaces(row['interfaces'])
    verify_recovery(row['recovery'])
    check_memory(row['memory'], code, maps)
    check_noise(row['continuous_noise'])
    check_accounting(row['accounting'], code, maps)
    verify_scaling(row['scaling'])


def verify(path=HERE/'receipt.json'):
    packet = parent.parent.load(path)
    keys(packet, 'schema sources claims evidence')
    exact(packet['schema'], 'oph-fixed-rate-source-v1')
    exact(pin_scope.bound_pins(packet['sources']), pins())
    exact(packet['claims'], claims())
    parent.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified complete decoder, bounded noisy recovery, native noise and fixed-rate operational bounds')
