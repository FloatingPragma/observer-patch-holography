"""Recompute sector-relative Weyl stress/modular identities and honest custody."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import algebra, entropy, geometry

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
OUTPUT = HERE / 'results.json'
PRIOR = 'model/continuum_free_limit/results.json'


def same_json(left, right):
    """Retain bool/integer distinctions in scientific and custody assertions."""
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def source_bindings():
    """Narrow mathematical dependency surface; never assert broad old custody."""
    paths = {p.relative_to(ROOT).as_posix() for p in HERE.rglob('*.py')}
    paths.update({
        'model/weyl_geometric_modular/README.md',
        'model/weyl_geometric_modular/imports.json',
        'paper/tex_fragments/WEYL_GEOMETRIC_MODULAR.tex',
        PRIOR,
        'model/continuum_free_limit/README.md',
        'model/continuum_free_limit/algebra.py',
        'model/continuum_free_limit/bounds.py',
        'model/continuum_field_net/README.md',
        'model/continuum_field_net/valleys.py',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean',
        'tools/strict_json.py',
    })
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def custody_diagnostic():
    prior = strict_json.load(ROOT/PRIOR)
    required = {
        'complete_centered_field_hierarchy_is_free': True,
        'massless_Weyl_field_count': 8,
        'Weyl_chiralities_each_sign': 4,
        'original_uncentered_ground_reidentified': False,
    }
    for key, expected in required.items():
        if type(prior['theorem'].get(key)) is not type(expected) or prior['theorem'][key] != expected:
            raise ArithmeticError('FA82 input theorem changed: '+key)
    mismatches = []
    for path, digest in sorted(prior['source_bindings'].items()):
        candidate = ROOT/path
        actual = hashlib.sha256(candidate.read_bytes()).hexdigest() if candidate.is_file() else None
        if actual != digest:
            mismatches.append({'path': path, 'recorded_sha256': digest, 'current_sha256': actual})
    return {
        'prior_receipt': PRIOR,
        'prior_binding_count': len(prior['source_bindings']),
        'current_prior_source_mismatch_count': len(mismatches),
        'current_prior_source_mismatches': mismatches,
        'prior_full_source_replay': False,
        'prior_expensive_experiments_reexecuted': False,
        'current_broad_source_map_matches': not mismatches,
        'historical_theorem_basis': strict_json.load(HERE/'imports.json')['historical_FA82'],
    }


def theorem_scope():
    return {
        'reference': 'FA82 centered bounded-U free eight-Weyl continuum vacuum',
        'source_clock': 'x0=v*t; model energy=v*P0; source hbar=1',
        'same_field_symmetric_conserved_flat_stress': True,
        'Ward_identity_independent_of_Einstein_equation': True,
        'same_vacuum_round_double_cone_modular_group': True,
        'modular_Wick_stress_identification_on_common_form_domain': True,
        'matter_first_law_on_explicit_normal_state_directions': True,
        'small_ball_matter_kernel_and_C2_remainder': True,
        'all_eight_valleys_retained': True,
        'absolute_type_III_entropy_assumed_finite': False,
        'specific_lattice_entropy_subtraction_convergence_proved': False,
        'native_event_count_geometry_identified': False,
        'common_lattice_collar_chart_limit_proved': False,
        'edge_area_density_or_Newton_constant_derived': False,
        'generalized_entropy_stationarity_proved': False,
        'interacting_or_massive_sector_supplied': False,
        'measured_particle_species_identified': False,
        'joint_EinsteinRegisterData_inhabited_by_this_field': False,
        'dynamical_Einstein_equation_derived': False,
        'full_G1_through_G6_closed': False,
    }


def validate(result):
    if result.get('schema') != 'oph.weyl_geometric_modular.v1':
        raise ArithmeticError('wrong Weyl modular schema')
    if not same_json(result.get('theorem'), theorem_scope()):
        raise ArithmeticError('unsupported theorem promotion or missing sector result')
    if result.get('finite_attachment') != 'FA87' or result.get('scope') != 'SAME_FREE_SECTOR_MATTER_MODULAR_THEOREM':
        raise ArithmeticError('wrong attachment or scope')
    if not same_json(result.get('imports'), strict_json.load(HERE/'imports.json')['imported_theorems']):
        raise ArithmeticError('imported theorem premises changed')
    expected = {'algebra': algebra.certificate(), 'geometry': geometry.certificate(),
                'entropy': entropy.certificate()}
    if not same_json(result.get('certificate'), expected):
        raise ArithmeticError('modular/stress identities or normalization changed')
    custody = result.get('custody', {})
    if not same_json(custody, custody_diagnostic()):
        raise ArithmeticError('prior custody diagnostic is stale or misrepresented')
    matched = custody['current_broad_source_map_matches']
    flags = {'new_checks_pass': True, 'dependency_custody_pass': matched,
             'all_pass': matched, 'integrated_gate_accepted': False,
             'physical_comparison_count_added': 0}
    for key, expected_value in flags.items():
        if type(result.get(key)) is not type(expected_value) or result[key] != expected_value:
            raise ArithmeticError('misleading status flag: '+key)
    if result.get('source_bindings') != source_bindings():
        raise ArithmeticError('current narrow mathematical source bindings changed')


def build():
    before = source_bindings()
    custody = custody_diagnostic()
    result = {
        'schema': 'oph.weyl_geometric_modular.v1',
        'finite_attachment': 'FA87',
        'scope': 'SAME_FREE_SECTOR_MATTER_MODULAR_THEOREM',
        'new_checks_pass': True,
        'dependency_custody_pass': custody['current_broad_source_map_matches'],
        'all_pass': custody['current_broad_source_map_matches'],
        'integrated_gate_accepted': False, 'physical_comparison_count_added': 0,
        'theorem': theorem_scope(),
        'certificate': {'algebra': algebra.certificate(), 'geometry': geometry.certificate(),
                        'entropy': entropy.certificate()},
        'imports': strict_json.load(HERE/'imports.json')['imported_theorems'],
        'custody': custody, 'source_bindings': before,
    }
    if source_bindings() != before:
        raise RuntimeError('source changed during Weyl modular replay')
    validate(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    group.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    result = build()
    if args.check:
        actual = strict_json.load(args.check)
        validate(actual)
        if not same_json(actual, result):
            raise ArithmeticError('Weyl modular receipt differs from cold replay')
    if args.write:
        args.write.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: result[k] for k in ('schema', 'new_checks_pass',
                     'dependency_custody_pass', 'all_pass')}))


if __name__ == '__main__':
    main()
