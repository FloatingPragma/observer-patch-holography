"""Replay exact high-field source coefficients and complete local preparation."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

from tools import strict_json
from model.defect_localization_readout import packet as native_dependencies
from . import algebra, preparation, readout

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PARENTS = (
    'model/high_field_particle/PROOF.md',
    'model/native_ordered_ground/source.py',
    'model/native_ordered_ground/__init__.py',
    'model/spin_twist_closure/field_gap.py',
    'model/spin_twist_closure/__init__.py',
    'model/spin_twist_closure/FINITE_FIELD_GAP.md',
    'model/spin_joint_field_limit/PROOF.md',
    'model/spin_fixed_field_locality/PROOF.md',
    'model/particle_identification/README.md',
    'model/defect_localization_readout/PROOF.md',
    'model/native_spin_observer/circuits.py',
    'model/native_spin_observer/precision.py',
    'model/native_spin_observer/algebra.py',
    'model/requirements.txt',
    'tools/strict_json.py',
)


def _exact(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: _exact(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_exact(item) for item in value]
    return value


def canonical(value):
    return json.dumps(_exact(value), sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    # Reuse the transitive native compiler/executor dependency closure, including
    # the source's non-Python data inputs; bind the collector itself as well.
    names = set(native_dependencies.bindings()) | set(PARENTS)
    names.update(p.relative_to(ROOT).as_posix() for p in HERE.iterdir()
                 if p.suffix in ('.py', '.md'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(names)}


def build():
    before = bindings()
    value = {
        'schema': 'oph.high_field_particle.v1',
        'source_bindings': before,
        'exact_source_algebra': algebra.certificate(),
        'complete_local_charge_creation': readout.creation_certificate(),
        'finite_ground_preparation_plan': readout.finite_preparation_example(),
        'short_complete_native_preparation': preparation.execution_certificate(),
        'bounded_mass_speed_controls': [
            {'spacing': str(a), 'J_over_B_ceiling': '1/100',
             'physical_rest_energy_ceiling': '1',
             'physical_speed_upper': str(readout.speed_bound(a, F(1, 100), F(1)))}
            for a in (F(1, 8), F(1, 64), F(1, 512), F(1, 4096))],
        'external_theorem_contract': {
            'isolated_band': 'https://arxiv.org/abs/math-ph/0411042v1',
            'gap': 'https://arxiv.org/abs/2607.10765v1',
            'velocity_bound': 'https://arxiv.org/abs/1412.2970v2',
            'theorems_reproved_by_finite_matrix_tests': False,
            'explicit_high_field_admission_threshold': None,
            'band_existence': 'for sufficiently small positive J/B',
        },
        'scope': {
            'exact_finite_coefficients_through_second_order': True,
            'complete_local_charge_instrument_certified': True,
            'same_source_finite_ground_planner_reused': True,
            'native_preparation_tape_executed_in_this_packet': True,
            'executed_preparation_scope': 'finite eight-site J=1 B=400 product-residual certificate',
            'finite_coupling_dispersion_remainder_certified': False,
            'uniform_numeric_original_particle_residue_floor': None,
            'bounded_J_over_B_finite_mass_local_relativistic_refinement_excluded': True,
            'exclusion_uses_same_ground_relative_clock_and_spatial_dictionary': True,
            'all_source_regimes_excluded': False,
            'physical_particle_assignments': [],
            'measured_mass_inputs_used': False,
            'physical_mass_predictions': [],
        },
    }
    if bindings() != before:
        raise RuntimeError('source changed during the high-field replay')
    return _exact(value)


def verify(path=OUTPUT):
    if canonical(strict_json.load(Path(path))) != canonical(build()):
        raise ValueError('high-field receipt differs from exact source replay or scope')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(build()))
        print(OUTPUT)
    else:
        verify()
        print('high-field source coefficients, local instrument and scope match')


if __name__ == '__main__':
    main()
