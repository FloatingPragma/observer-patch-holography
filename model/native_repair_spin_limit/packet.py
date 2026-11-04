"""Replay exact native repair controls, analytic bounds and source custody."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra, bounds

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')
PARENTS = {
    'model/bipartite_spin_limit/results.json': 'oph.bipartite_spin_limit.v1',
    'model/defect_field_bridge/full_domain_results.json': 'oph.full_spatial_CAR_repair.v1',
}


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/native_repair_spin_limit/README.md',
        'source_receipts': list(PARENTS),
        'independent_checks': ['model/native_repair_spin_limit/PROOF.md'],
    })
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed after import; replay in a fresh interpreter')
    for path, schema in PARENTS.items():
        parent = strict_json.load(ROOT/path)
        if parent.get('schema') != schema:
            raise ValueError('specified source parent required: '+path)
        if path.endswith('bipartite_spin_limit/results.json') and parent.get('all_pass') is not True:
            raise ValueError('verified local spin theorem required')
        for source, expected in parent.get('source_bindings', parent.get('source_sha256', {})).items():
            if pins[source] != expected:
                raise ValueError('inherited source custody changed: '+source)
    value = {
        'schema': 'oph.native_repair_spin_limit.v1', 'all_pass': True,
        'issues': [1138, 1142], 'gates': ['CP-02', 'CP-06'],
        'source_bindings': pins,
        'bounds': bounds.certificate(),
        'complete_source_controls': algebra.certificate(),
        'scope': 'Fixed positive homogeneous native C or E repair rate on every edge, bounded positive initial source energy density and period-two source marginal; complete local records and volume-independent interacting spin limit at fixed spacing.',
        'input_state_preparation_is_specified': True,
        'parent_source_custody_checked': True,
        'inherited_expensive_source_pulse_experiments_reexecuted': False,
        'scope_exclusions': {
            'arbitrary_feedback': True,
            'thermal_selection_by_repair': True,
            'infinite_archive_global_trace_norm': True,
            'relativistic_continuum_or_Einstein_completion': True,
        },
    }
    if source_bindings() != pins:
        raise ValueError('source changed during native repair replay')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('native repair spin proof, bounds, scope or custody changed')


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
    print('fixed positive native repair rate and local interacting spin limit: PASS')


if __name__ == '__main__':
    main()
