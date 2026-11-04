"""Source-bound replay of the nonlinear ground-texture energy theorem."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import energy

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENT = 'model/spin_canonical_energy/results.json'


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/spin_texture_energy/README.md',
        'source_receipts': [PARENT],
        'independent_checks': ['model/spin_texture_energy/PROOF.md'],
    })
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


_IMPORTED_BINDINGS = bindings()


def report():
    pins = bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; replay in a fresh interpreter')
    parent = strict_json.load(ROOT/PARENT)
    if parent.get('schema') != 'oph.spin_canonical_energy.v1':
        raise ValueError('specified actual-ground source required')
    for path, expected in parent['source_bindings'].items():
        if pins[path] != expected:
            raise ValueError('inherited energy source changed: '+path)
    value = {'schema': 'oph.spin_texture_energy.v1', 'all_pass': True,
             'source_bindings': pins, 'energy': energy.certificate(),
             'scope': 'Static nonlinear energy on the finite local-unitary orbit of the actual selected ordered spin ground; fixed-spacing and controlled long-wavelength texture families.',
             'parent_source_custody_checked': True,
             'inherited_expensive_producers_reexecuted': False}
    if bindings() != pins:
        raise ValueError('source changed during texture replay')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('texture energy proof, scope or source custody changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    if args.write:
        args.write.write_text(canonical(report()))
    else:
        validate(strict_json.load(args.check))
    print('actual-ground nonlinear texture energy and topology scope: PASS')


if __name__ == '__main__':
    main()
