"""Replay exact complete-clock tube bounds and the finite spin witness."""
import argparse
import ast
import importlib.util
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import bounds, finite, threshold

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
OUTPUT = HERE/'results.json'
INPUTS = (
    'model/native_ordered_ground/source.py',
    'model/native_ordered_ground/certificate.py',
    'model/native_ordered_ground/witness.npz',
    'model/native_ordered_ground/README.md',
    'model/native_joint_spin_spectroscopy/algebra.py',
    'model/native_joint_spin_spectroscopy/spectral.py',
    'model/native_joint_spin_spectroscopy/README.md',
    'model/native_spin_observer/precision.py',
    'model/centered_particle_spectrum/certificate.py',
    'model/spin_shell_diagnostic/README.md',
    'model/spin_packet_readout/README.md',
    'model/requirements.txt', 'tools/strict_json.py')


def encode(value):
    if type(value) is F:
        return str(value)
    if type(value) in (tuple, list):
        return [encode(x) for x in value]
    if type(value) is dict:
        return {k: encode(v) for k, v in value.items()}
    if value is None or type(value) in (str, int, bool):
        return value
    raise TypeError('exact typed scientific value required')


def canonical(value):
    return json.dumps(encode(value), sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    legacy='evidence/observer_dynamics_20260925/oph-physics-sim/'
    names=set(INPUTS)|{legacy+'oph_exact/carrier.py',legacy+'oph_exact/__init__.py'}
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py','.md'))
    pending=list(names);visited=set()
    def add(module):
        parts=module.split('.')
        prefix=legacy if parts[0] in ('oph_exact','oph_fpe') else ''
        candidates=[prefix+'/'.join(parts)+'.py',prefix+'/'.join(parts)+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1,len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name);pending.append(name)
    while pending:
        name=pending.pop()
        if name in visited or not name.endswith('.py'):continue
        visited.add(name)
        module=name.removeprefix(legacy)[:-3].replace('/','.')
        package=module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:add(alias.name)
            elif isinstance(node,ast.ImportFrom):
                target='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(target,package) if node.level else target
                add(target)
                for alias in node.names:add(target+'.'+alias.name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def produce():
    before = bindings()
    result = {
        'schema': 'oph.spin_shell_readout.v1', 'source_bindings': before,
        'complete_clock_example': bounds.clock_plan(F(1, 1000), F(1, 100)),
        'finite_source_readout': finite.certificate_report(),
        'infrared_threshold_readout': threshold.certificate(),
        'scope': {
            'complete_original_spectral_tube_weight_bounds': True,
            'simultaneous_predeclared_scan_sampling_bound': True,
            'same_selected_finite_spin_state_and_clock': True,
            'finite_source_full_integer_witness_rechecked': True,
            'full128bin_native_tape_executed': False,
            'infinite_ordered_state_preparation_modulus_derived': False,
            'fixed_history_convergence_controls_growing_spectrometer': False,
            'uniform_positive_shrinking_tube_fraction_proved': False,
            'source_selected_common_relativistic_speed_derived': False,
            'sharp_continuum_mass_shell_or_measured_mass_derived': False,
            'inherited_native_compiler_owner_replayed': False}}
    if before != bindings():
        raise ArithmeticError('consumed source changed during spectral-tube replay')
    return encode(result)


def validate(value, expected=None):
    expected = produce() if expected is None else expected
    if canonical(value) != canonical(expected):
        raise ArithmeticError('complete spectral-tube certificate does not replay')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', nargs='?', const=str(OUTPUT), metavar='RECEIPT')
    args = parser.parse_args()
    result = produce()
    if args.write:
        OUTPUT.write_text(canonical(result), encoding='utf-8')
    else:
        validate(strict_json.load(Path(args.check)), result)
    print('Complete spectral-tube source certificate verified.')


if __name__ == '__main__':
    main()
