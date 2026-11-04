"""Cold replay of all rational spectral intervals and full-Fock consequences."""
import argparse
import hashlib
import json
from pathlib import Path

import model
from tools import strict_json
from .catalogue import certified_case

ROOT = Path(model.__file__).resolve().parent.parent
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/winding_sector_search/results.json', 'model/radial_interacting_core/results.json')
OWN_FILES = ('__init__.py', '__main__.py', 'intervals.py', 'spectrum.py', 'core.py',
             'catalogue.py', 'packet.py', 'test_spectrum.py', 'test_packet.py', 'PROOF.md', 'README.md')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        if parent.get('all_pass') is not True or type(parent.get('source_bindings')) is not dict:
            raise ValueError('parent has no complete passing manifest')
        for path, digest in parent['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('stale parent source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('parent source manifests disagree')
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS)
    paths.update('model/certified_winding_core/'+name for name in OWN_FILES)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def report():
    pins = source_bindings()
    search = strict_json.load(ROOT/PARENTS[0])['finite_search']
    targets = search['all_initial_cases']+search['volume_followup']['all_cases']
    if len(targets) != 28:
        raise ValueError('the complete retained search family changed')
    cases = [certified_case(t) for t in targets]
    value = {'schema': 'oph.certified_winding_core.v1', 'all_pass': True,
             'complete_cases': cases, 'case_count': len(cases),
             'initial_grid': search['declared_search_grid'],
             'adaptive_followup': {k: v for k, v in search['volume_followup'].items() if k != 'all_cases'},
             'exact_two_zero_kernel_cases': sum(c['exact_kernel_dimension'] == 2 for c in cases),
             'certified_invertible_cases': sum(c['exact_kernel_dimension'] == 0 for c in cases),
             'unresolved_kernel_cases': sum(c['exact_kernel_dimension'] is None for c in cases),
             'admitted_interacting_contact_cases': sum(b['bound_proves_ground_and_charged_pair']
                 for c in cases for b in c['full_Fock_contact_cases']),
             'scope': {'all_finite_source_eigenvalues_retained': True,
                       'exact_binary64_limb_accumulation_domain_checked': True,
                       'source_matrix_rounding_included': True,
                       'numerical_residual_supplies_the_proof': False,
                       'zero_mode_index_identifies_particle_generations': False,
                       'supplied_parameters_are_selected_by_native_repair': False,
                       'full_Fock_complement_is_retained': True,
                       'common_interacting_relativistic_limit': False,
                       'observed_particle_identification': False},
             'canonical_physical_premises_discharged': [], 'physical_comparison_count_added': 0,
             'integrated_gate_accepted': False, 'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during certification')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('receipt differs from the full fresh certification')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=str(OUTPUT))
    group.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        path = Path(args.write); path.write_text(canonical(report()))
        print('wrote '+str(path))
    else:
        validate(strict_json.load(Path(args.check)))
        print('all28 complete spectra and supplied-contact cases certified')


if __name__ == '__main__':
    main()
