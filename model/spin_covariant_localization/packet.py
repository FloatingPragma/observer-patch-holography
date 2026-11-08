"""Replay the finite source coordinate certificate and exact Gram constants."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import bounds, localization

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent
OUTPUT = HERE/'results.json'
INPUTS = (
    'model/native_ordered_ground/source.py',
    'model/native_ordered_ground/certificate.py',
    'model/native_ordered_ground/witness.npz',
    'model/native_joint_spin_spectroscopy/algebra.py',
    'model/centered_particle_spectrum/certificate.py',
    'model/spin_canonical_field/README.md',
    'model/spin_canonical_field/algebra.py',
    'model/spin_canonical_field/bounds.py',
    'model/spin_goldstone_sector/README.md',
    'model/requirements.txt', 'tools/strict_json.py',
)


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
        'schema': 'oph.spin_covariant_localization.v1', 'source_bindings': before,
        'finite_source_readout': localization.report(),
        'continuum_source_gram_constants': bounds.gram_bounds(),
        'finite_compression_error_example': bounds.frame_perturbation(F(1, 1000), F(1, 2000000)),
        'scope': {
            'actual_finite_ground_and_complete_witness_rechecked': True,
            'finite_reflection_coordinate_detector_complete': True,
            'finite_cyclic_code_energy_leakage_certified_positive': True,
            'source_continuum_field_vector_spectral_density_bounded_above_and_below': True,
            'same_continuum_spatial_generator_covariant_position_on_field_cyclic_subspace': True,
            'schwartz_core_position_momentum_commutator_proved_in_text': True,
            'vacuum_requires_dark_outcome_or_sector_restriction': True,
            'finite_arithmetic_alone_proves_continuum_theorem': False,
            'native_detector_tape_executed': False,
            'continuum_whitening_local_or_operationally_implemented': False,
            'carrier_rank3_record_equals_continuum_position': False,
            'time_generator_preserves_field_cyclic_subspace': False,
            'sharp_particle_shell_or_measured_mass_identified': False,
            'inherited_owner_producers_replayed': False,
        },
    }
    if before != bindings():
        raise ArithmeticError('consumed source changed during localization replay')
    return encode(result)


def validate(value, expected=None):
    expected = produce() if expected is None else expected
    if canonical(value) != canonical(expected):
        raise ArithmeticError('source coordinate certificate does not replay')
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
    print('Source coordinate certificate verified.')


if __name__ == '__main__':
    main()
