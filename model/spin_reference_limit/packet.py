"""Source-bound replay of the finite and local reference estimates."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENT = 'model/bipartite_spin_limit/results.json'


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/spin_reference_limit/README.md',
        'source_receipts': [PARENT],
        'independent_checks': ['model/spin_reference_limit/PROOF.md'],
    })
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; use a fresh interpreter')
    parent = strict_json.load(ROOT/PARENT)
    if parent.get('all_pass') is not True or parent.get('schema') != 'oph.bipartite_spin_limit.v1':
        raise ValueError('verified local spin source required')
    for path, expected in parent['source_bindings'].items():
        if pins[path] != expected:
            raise ValueError('local spin source custody changed: '+path)
    result = {
        'schema': 'oph.spin_reference_limit.v1', 'all_pass': True,
        'issues': [1138, 1142], 'gates': ['CP-02', 'CP-06'],
        'bounds': bounds.certificate(), 'complete_carrier_controls': algebra.certificate(),
        'source_bindings': pins,
        'scope': 'Actual symmetry-averaged full-source grounds, faithful full-source Gibbs references, local ground condition and volume-uniform pressure bounds at fixed lattice spacing.',
        'state_selection_rule_is_a_specified_input': True,
        'parent_source_custody_checked': True,
        'parent_expensive_controls_reexecuted': False,
    }
    if pins != source_bindings():
        raise ValueError('source changed during reference replay')
    return result


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('source reference proof, bounds or custody changed')


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
    print('source ground and thermal spin references: PASS')


if __name__ == '__main__':
    main()
