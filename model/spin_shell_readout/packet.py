"""Replay exact complete-clock tube bounds and the finite spin witness."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import bounds, finite

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
    paths = set(INPUTS)
    paths.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def produce():
    before = bindings()
    result = {
        'schema': 'oph.spin_shell_readout.v1', 'source_bindings': before,
        'complete_clock_example': bounds.clock_plan(F(1, 1000), F(1, 100)),
        'finite_source_readout': finite.certificate_report(),
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
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = produce()
    if args.write:
        OUTPUT.write_text(canonical(result), encoding='utf-8')
    else:
        validate(strict_json.load(OUTPUT), result)
    print('Complete spectral-tube source certificate verified.')


if __name__ == '__main__':
    main()
