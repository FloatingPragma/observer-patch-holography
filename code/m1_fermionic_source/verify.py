"""Strict custody envelope and independent fermionic-source replay."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = 'm1_fermionic_source'

from m1_source_realization import verify as parent
from .check import exact, keys, need, verify_gadgets, verify_small
from .experiment_check import (verify_conversions, verify_detector, verify_interaction,
                               verify_measurements, verify_noise, verify_resolved_reads, verify_weighted_measurements)
from .geometry_check import verify_budgets, verify_clock_witness, verify_graphs, verify_walks
from .preparation_check import verify_preparation
import pin_scope

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OWN_FILES = ('__init__.py', 'pauli.py', 'model.py', 'circuits.py', 'spatial.py', 'walk.py',
             'preparation.py', 'experiments.py', 'check.py', 'geometry_check.py',
             'preparation_check.py', 'experiment_check.py', 'verify.py', 'build.py',
             'test_fermionic.py', 'README.md', 'CONTRACT.md')
SOURCES = sorted(set(parent.SOURCES + ['code/m1_source_realization/receipt.json']
                    + ['code/m1_fermionic_source/'+p for p in OWN_FILES]
                    + ['docs/research/FERMIONIC_SOURCE_CLOCKS.md', '.github/workflows/m1-fermionic-source.yml']))
SOURCES = pin_scope.bound(SOURCES)
CLAIMS = ('OPH-SOURCE-LOCAL-FERMIONIC-REALIZATION', 'OPH-SOURCE-ENCODED-FERMION-CLOCK',
          'OPH-SOURCE-FERMIONIC-VACUUM-NOISE')


def pins():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def claim_pins():
    rows = parent.load(ROOT/'claims/claim_registry.yaml')['claims']
    out = {}
    for name in CLAIMS:
        matches = [r for r in rows if r['claim_id'] == name]
        need(len(matches) == 1, 'unique claim '+name)
        out[name] = hashlib.sha256(json.dumps(matches[0], sort_keys=True, separators=(',', ':'),
                                            ensure_ascii=True).encode('ascii')).hexdigest()
    return out


def verify_evidence(evidence):
    keys(evidence, 'small_graphs native_rotations native_measurements weighted_measurements geometry preparation '
                   'budgets walks conversions resolved_reads detector interaction noise clock_witness')
    cache = verify_small(evidence['small_graphs'])
    verify_gadgets(evidence['native_rotations'])
    verify_measurements(evidence['native_measurements'])
    verify_weighted_measurements(evidence['weighted_measurements'])
    verify_graphs(evidence['geometry'])
    verify_preparation(evidence['preparation'])
    verify_budgets(evidence['budgets'])
    verify_walks(evidence['walks'])
    verify_conversions(evidence['conversions'])
    verify_resolved_reads(evidence['resolved_reads'], cache)
    verify_detector(evidence['detector'], cache)
    verify_interaction(evidence['interaction'], cache)
    verify_noise(evidence['noise'], cache)
    verify_clock_witness(evidence['clock_witness'])


def verify(path=HERE/'receipt.json'):
    packet = parent.load(path)
    keys(packet, 'schema sources claims evidence')
    exact(packet['schema'], 'oph-fermionic-source-v1')
    exact(pin_scope.bound_pins(packet['sources']), pins())
    exact(packet['claims'], claim_pins())
    parent.verify()
    verify_evidence(packet['evidence'])
    return packet


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('path', nargs='?', type=Path, default=HERE/'receipt.json')
    verify(parser.parse_args().path)
    print('Verified both parity blocks, local preparation, native Fock clock, reads, interaction and noise')
