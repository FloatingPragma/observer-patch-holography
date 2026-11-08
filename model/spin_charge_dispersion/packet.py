"""Source-bound certificate for bounded charged dispersive-band readouts."""
import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.spin_dynamical_memory import packet as memory_packet
from model.protected_defect_catalog import field
from . import bounds, dimer, finite, order

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PARENTS = ('Lean/Geometry/MassShellKinematics.lean',
           'Lean/Screen/LightSignalAdequacySurface.lean',
           'Lean/QFT/SourceHistoryGNSDynamics.lean',
           'Lean/Screen/FieldSectorEnergyInnerProduct.lean',
           'model/spin_shell_diagnostic/README.md',
           'model/spin_current_spectrum/README.md',
           'model/spin_goldstone_sector/README.md',
           'model/spin_goldstone_sector/infrared.py',
           'model/spin_goldstone_sector/__init__.py',
           'model/spin_operational_position/PROOF.md')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    result = memory_packet.bindings()
    names = set(PARENTS)
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir() if path.suffix in ('.py', '.md'))
    names.add('model/spin_operational_position/bounds.py')
    names.add('model/spin_operational_position/__init__.py')
    for name in names:
        result[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(sorted(result.items()))


def build():
    before = bindings()
    value = {'schema': 'oph.spin_charge_dispersion.v1', 'source_bindings': before,
        'actual_configuration_charge': field.configuration_charge(),
        'stronger_rigorous_source_order': order.certificate(),
        'source_dimer_energy_and_current_bound': dimer.certificate(),
        'bounded_creator_positive_joint_band': bounds.band(),
        'complete_event_error_example': bounds.event_budget(),
        'exact_finite_source_creator': finite.build(),
        'scope': {'same_source_state_and_clock': True,
            'bounded_exact_charge_plus_and_minus_preparations_constructed': True,
            'uniform_original_positive_joint_energy_momentum_band_weight': True,
            'all_preparation_and_band_complements_retained': True,
            'full_probability_and_state_limit_is_an_ordered_qualitative_limit': True,
            'effective_infinite_ground_preparation_modulus': False,
            'all_energy_momentum_weight_forced_onto_one_shell': False,
            'sharp_nonzero_residue_dispersion_or_Lorentz_covariance_proved': False,
            'physical_particle_or_mass_assignment': None}}
    if before != bindings():
        raise RuntimeError('charged-band inputs changed during replay')
    return value


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('bounded charged source readout or scope differs')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    mode.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    args = parser.parse_args()
    if args.write:
        args.write.write_text(canonical(build()))
    else:
        validate(strict_json.load(args.check))
    print('PASS: bounded charged source preparation and positive joint band')
