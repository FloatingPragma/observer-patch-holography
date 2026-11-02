"""Replay the complete onsite mass catalogue with the actual centered reference.

The centered packet has its own independent cold replay in model.verify. This
consumer checks that parent's declared conclusions and complete current source
manifest, then independently replays the mass algebra. It does not rerun the
parent's large integer Gram proof for each mass-catalogue validation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..centered_source_ground import packet as centered_packet
from . import algebra


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent/'results.json'
PARENT = 'model/centered_source_ground/results.json'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/native_mass_channels/README.md',
        'source_receipts': [PARENT], 'independent_checks': [],
    })
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(paths)}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def validate(value):
    required = {'schema', 'all_pass', 'catalogue', 'actual_centered_ground',
                'physical_comparison_count_added', 'integrated_gate_accepted',
                'source_bindings'}
    if (type(value) is not dict or set(value) != required
            or value.get('schema') != 'oph.native_mass_channels.receipt.v1'
            or value.get('all_pass') is not True):
        raise ValueError('invalid native mass catalogue receipt')
    # Canonical JSON distinguishes true from 1 and false from 0, unlike dict
    # equality. Every nested proof claim and scope flag participates.
    if canonical(value['catalogue']) != canonical(algebra.certificate()):
        raise ValueError('the complete mass/texture classification changed')
    actual = {'full_finite_gap_strict_lower': '1/128',
              'onsite_mass_order_expectations': ['0']*4,
              'finite_gap_identified_with_a_classical_mass_order': False}
    if canonical(value['actual_centered_ground']) != canonical(actual):
        raise ValueError('the actual centered reference was changed or overinterpreted')
    if (type(value['physical_comparison_count_added']) is not int
            or value['physical_comparison_count_added'] != 0
            or value['integrated_gate_accepted'] is not False):
        raise ValueError('a finite catalogue cannot accept a physical completion gate')
    if value['source_bindings'] != source_bindings():
        raise ValueError('mass catalogue source manifest is incomplete or stale')


def _produce():
    pins = source_bindings()
    parent = strict_json.load(ROOT/PARENT)
    if (type(parent) is not dict
            or parent.get('source_bindings') != centered_packet.source_bindings()):
        raise ValueError('centered reference source manifest is incomplete or stale')
    try:
        centered_packet.validate(parent)
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError('malformed centered reference proof') from error
    value = {
        'schema': 'oph.native_mass_channels.receipt.v1', 'all_pass': True,
        'catalogue': algebra.certificate(),
        'actual_centered_ground': {
            'full_finite_gap_strict_lower': parent['finite_gap']['full_centered_gap_strict_lower'],
            'onsite_mass_order_expectations': parent['charged_particles']['onsite_mass_order']['actual_unique_ground_order_expectations'],
            'finite_gap_identified_with_a_classical_mass_order': False,
        },
        'physical_comparison_count_added': 0, 'integrated_gate_accepted': False,
        'source_bindings': pins,
    }
    validate(value)
    if source_bindings() != pins:
        raise ValueError('mass catalogue sources changed during replay')
    return value


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-mass-catalogue-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, PYTHONUTF8='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(folder)/'pycache'),
                        '-m', 'model.native_mass_channels.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        value = strict_json.load(path)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('sources changed across cold mass replay')
    validate(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = _produce() if args.worker else report()
    if args.check:
        retained = strict_json.load(args.check)
        validate(retained)
        if canonical(retained) != canonical(value):
            raise ValueError('mass catalogue or source bindings changed')
    else:
        (args.worker or args.write).write_text(canonical(value), encoding='utf-8', newline='\n')
    if not args.worker:
        print(json.dumps({'all_pass': True, 'source_bindings': len(value['source_bindings']),
                          'physical_species_identified': False}))


if __name__ == '__main__':
    main()
