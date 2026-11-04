"""Replay exact bipartite identities and bind their analytic proof."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENT = 'model/strong_coupling_refinement/results.json'


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/bipartite_spin_limit/README.md',
        'source_receipts': [PARENT], 'independent_checks': [],
    })
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; replay in a fresh interpreter')
    parent = strict_json.load(ROOT/PARENT)
    mismatches = sorted(path for path, expected in parent['source_bindings'].items()
                        if pins[path] != expected)
    if mismatches or parent.get('new_checks_pass') is not True:
        raise ValueError('parent source custody or finite checks failed: '+str(mismatches))
    result = {
        'schema': 'oph.bipartite_spin_limit.v1',
        'all_pass': True,
        'source_bindings': pins,
        'bounds': bounds.certificate(),
        'full_sector_controls': algebra.certificate(),
        'scope': {
            'source': 'same centered signed cubic Hubbard graph at half filling',
            'spectral_domain': 'finite simple bipartite graph, U>tau=c|E|, c>0',
            'canonical_dynamics_domain': 'same graph, U>=4tau',
            'ordered_energy_error': '2tau^4/[U^2(U-tau)]',
            'canonical_generator_error': '7tau^4/U^3',
            'complete_dynamics_error': 'min(2,4tau/U+7|t|tau^4/U^3)',
            'all_virtual_charge_states_retained': True,
            'infinite_volume_estimate': False,
            'interacting_relativistic_continuum': False,
            'physical_particle_assignment': False,
            'integrated_gate_accepted': False,
        },
        'parent_receipt_custody_checked': True,
        'parent_retained_aggregate_pass': parent.get('all_pass') is True,
        'parent_expensive_experiments_replayed': False,
    }
    if source_bindings() != pins:
        raise ValueError('source changed during replay')
    return result


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('bipartite spin proof, bounds, scope or source custody changed')


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
    print('bipartite spin identities, bounds and source custody: PASS')


if __name__ == '__main__':
    main()
