"""Full replay of exact index algebra and every finite search candidate."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import model
from tools import strict_json
from . import algebra, search

ROOT = Path(model.__file__).resolve().parent.parent
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/native_mass_channels/results.json',
           'model/bounded_defect_index/results.json',
           'model/radial_defect_gap/results.json')
OWN_FILES = ('__init__.py', '__main__.py', 'algebra.py', 'search.py', 'packet.py',
             'test_search.py', 'test_packet.py', 'PROOF.md', 'README.md')
DIRECT_SOURCES = ('model/topological_defect_spectrum/algebra.py',
                  'model/topological_defect_spectrum/PROOF.md',
                  'model/radial_defect_gap/bounds.py',
                  'model/radial_defect_gap/PROOF.md')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        if parent.get('all_pass') is not True or type(parent.get('source_bindings')) is not dict:
            raise ValueError('a parent has no complete passing source manifest')
        for path, digest in parent['source_bindings'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
                raise ValueError('stale parent source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('parent source manifests disagree')
            inherited[path] = digest
    paths = set(inherited) | set(PARENTS) | set(DIRECT_SOURCES)
    paths.update('model/winding_sector_search/'+name for name in OWN_FILES)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def report():
    pins = source_bindings()
    value = {
        'schema': 'oph.winding_sector_search.v1', 'all_pass': True,
        'continuum_index': algebra.index_certificate(),
        'degree_family': [algebra.degree_certificate(d) for d in (-4, -3, -2, -1, 1, 2, 3, 4)],
        'strong_mass_examples': [algebra.strong_mass_certificate(d, mass=4)
                                 for d in (-5, -3, -1, 1, 3, 5)],
        'finite_search': search.search_report(),
        'scope': {
            'prescribed_texture': True,
            'all_source_spin_and_valley_rows': True,
            'infinite_lattice_strong_mass_odd_texture_kernel_theorem': True,
            'analytic_Callias_theorem_imported_with_hypotheses': True,
            'index_proved_by_numerical_diagonalization': False,
            'degree_selected_by_source_dynamics': False,
            'mass_couplings_selected': False,
            'same_current_particle_assignment': False,
            'native_repair_protection': False,
            'interacting_relativistic_particle_spectrum': False,
            'GeV_masses_derived': False,
        },
        'physical_comparison_count_added': 0, 'integrated_gate_accepted': False,
        'source_bindings': pins,
    }
    if source_bindings() != pins:
        raise ValueError('sources changed during the calculation')
    return value


def compare(actual, expected, path='receipt'):
    """Exact structure/types and exact algebra; tolerance only for float diagnostics."""
    if type(actual) is not type(expected):
        raise ValueError('fresh calculation differs at '+path)
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError('fresh calculation keys differ at '+path)
        for key in expected:
            compare(actual[key], expected[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError('fresh calculation length differs at '+path)
        for index, (x, y) in enumerate(zip(actual, expected, strict=True)):
            compare(x, y, path+f'[{index}]')
    elif type(expected) is float:
        field = path.rsplit('.', 1)[-1]
        diagnostic = field.startswith('all_eigenvalues[') or field in {
            'window_origin_projector_trace', 'window_boundary_projector_trace',
            'window_total_radius_squared_trace',
        }
        if not (math.isfinite(actual) and math.isfinite(expected)
                and (math.isclose(actual, expected, rel_tol=1e-11, abs_tol=1e-10)
                     if diagnostic else actual == expected)):
            raise ValueError('fresh calculation diagnostic differs at '+path)
    elif actual != expected:
        raise ValueError('fresh calculation differs at '+path)


def validate(value):
    compare(value, report())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument('--write', nargs='?', const=str(OUTPUT))
    operation.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        path = Path(args.write)
        path.write_text(canonical(report()))
        print('wrote '+str(path))
    else:
        validate(strict_json.load(Path(args.check)))
        print('all 28 winding-search cases and exact index algebra replayed')


if __name__ == '__main__':
    main()
