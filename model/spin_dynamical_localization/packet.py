"""Replay the source dynamic-position certificate and proof input custody."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import finite

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
INPUTS = (
    'model/spin_canonical_field/README.md',
    'model/spin_covariant_localization/CONTINUUM_POSITION.md',
    'model/spin_covariant_localization/PROOF.md',
    'model/native_ordered_ground/README.md',
    'model/native_ordered_ground/witness.npz',
    'model/requirements.txt',
)


def bindings():
    paths = set(INPUTS)
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    pending = list(paths)
    visited = set()

    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        stem = '/'.join(parts)
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in paths:
                paths.add(name)
                pending.append(name)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def produce():
    before = bindings()
    result = {
        'schema': 'oph.spin_dynamical_localization.v1',
        'source_bindings': before,
        'finite_source_readout': finite.report(),
        'scope': {
            'same_continuum_state_spatial_generator_and_time_clock': True,
            'time_invariant_complete_response_space_constructed_in_proof': True,
            'source_history_gram_constructs_covariant_position_povm': True,
            'restriction_to_original_field_position_pvm_exact': True,
            'dark_outcome_commutes_with_actual_time_generator': True,
            'uncertainty_for_complete_position_outcomes_when_moments_finite': True,
            'finite_moments_preserved_under_all_source_times': False,
            'source_fiber_frame_smoothness_proved': False,
            'finite_arithmetic_alone_proves_continuum_theorem': False,
            'native_position_instrument_tape_executed': False,
            'carrier_load_position_equals_this_detector': False,
            'local_continuum_detector_realization_proved': False,
            'sharp_particle_mass_shell_or_observed_species_identified': False,
            'all_inherited_owner_producers_replayed': False,
        },
    }
    if before != bindings():
        raise ArithmeticError('source changed during dynamic-position replay')
    return result


def validate(value, expected=None):
    if type(value) is not dict or value.get('schema') != 'oph.spin_dynamical_localization.v1':
        raise ValueError('invalid dynamic position receipt')
    expected = produce() if expected is None else expected
    if canonical(value) != canonical(expected):
        raise ValueError('dynamic position receipt does not reproduce')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--write', action='store_true')
    action.add_argument('--check', nargs='?', const=str(OUTPUT), metavar='RECEIPT')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()), encoding='utf-8')
    else:
        validate(strict_json.load(Path(args.check)))
    print('Complete source dynamic-position certificate verified.')


if __name__ == '__main__':
    main()
