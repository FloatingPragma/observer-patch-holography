"""Source-bound replay of exact lattice degree and energy selection bounds."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import energy, repair, shells

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'
PARENTS = ('model/spin_texture_energy/results.json',
           'model/native_repair_spin_limit/results.json')
DIRECT = ('model/__init__.py', 'model/requirements.txt', 'tools/strict_json.py',
          'model/spin_texture_energy/PROOF.md', 'model/native_repair_spin_limit/PROOF.md',
          'model/defect_field_bridge/full_domain.py', 'model/native_repair_spin_limit/algebra.py')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    inherited = {}
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        if parent.get('all_pass') is not True or type(parent.get('source_bindings')) is not dict:
            raise ValueError('the exact parent source receipt must have complete passing custody')
        for path, digest in parent['source_bindings'].items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != digest:
                raise ValueError('stale or external source: '+path)
            if path in inherited and inherited[path] != digest:
                raise ValueError('parent source manifests disagree')
            inherited[path] = digest
    paths = set(inherited)|set(PARENTS)|set(DIRECT)
    paths.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir() if p.is_file() and p.suffix in ('.py', '.md'))
    if any(p.startswith(('tracking/', 'model/completion/')) or p == 'model/verify.py' for p in paths):
        raise ValueError('central integration files cannot enter the scientific source closure')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def report():
    pins = source_bindings()
    value = {'schema': 'oph.lattice_texture_selection.v1', 'all_pass': True,
             'exact_degree_and_edges': shells.certificate(),
             'same_source_energy': energy.certificate(),
             'complete_native_repair': repair.certificate(),
             'scope': {
                 'actual_selected_ordered_quantum_ground_rotation_orbit': True,
                 'ordered_ground_selection_is_an_inherited_input': True,
                 'nearest_neighbor_source_energy_used': True,
                 'finite_acute_interpolation_hypothesis_required': True,
                 'unscreened_nonzero_degree_on_infinitely_many_admissible_shells_has_infinite_energy': True,
                 'finite_energy_derives_eventual_shell_admissibility_and_zero_degree': True,
                 'exact_regular_target_degree_uses_all_triangles': True,
                 'controlled_finite_unwinding_path_constructed': True,
                 'absolute_crossing_floor_is_lifetime_or_activation_barrier': False,
                 'native_repair_selects_winding_or_spin_direction': False,
                 'interleaved_Hamiltonian_evolution_is_repair_identity': False,
                 'full_source_vacuum_or_texture_formation_rule_derived': False,
                 'spin_direction_equals_prescribed_Dirac_mass_texture': False,
                 'all_quantum_defects_or_particles_excluded': False,
                 'physical_mass_scale_or_species_selected': False,
                 'interacting_relativistic_QFT_or_Einstein_completed': False,
             },
             'canonical_physical_premises_discharged': [],
             'physical_comparison_count_added': 0, 'integrated_gate_accepted': False,
             'source_bindings': pins}
    if source_bindings() != pins:
        raise ValueError('source changed during the complete calculation')
    return value


def validate(value):
    if canonical(value) != canonical(report()):
        raise ValueError('lattice texture receipt differs from the exact fresh calculation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    operation = parser.add_mutually_exclusive_group(required=True)
    operation.add_argument('--write', nargs='?', const=str(OUTPUT))
    operation.add_argument('--check', nargs='?', const=str(OUTPUT))
    args = parser.parse_args()
    if args.write:
        Path(args.write).write_text(canonical(report()))
        print('wrote '+args.write)
    else:
        validate(strict_json.load(Path(args.check)))
        print('exact lattice degree, source energy and repair scope replay passed')


if __name__ == '__main__':
    main()
