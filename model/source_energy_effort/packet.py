"""Same-source charge work, energy and signed action certificate."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from model.spin_modular_geometry.packet import compare as compare_native_records, rational
from . import energy, readout, clock

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent

# The inherited complete-record comparator admits this absolute error on each
# probability. A difference quotient must propagate that error through 1/time.
_RECORD_ERROR = 2e-11
_CLOCK_DERIVATIVE = 'unnormalized_mean_energy_from_sine_derivative'


def _clock_derivative(value, path):
    try:
        clock_value = value['native_clock_action']
        duration = clock_value['time']
        probabilities = clock_value['joint_charge_clock_success_probabilities']
        derivative = clock_value[_CLOCK_DERIVATIVE]
        if type(duration) is not str:
            raise ValueError('time must retain its exact rational string')
        time = float(F(duration))
        if not math.isfinite(time) or time <= 0:
            raise ValueError('time must be positive and finite')
        if (type(probabilities) is not list or len(probabilities) != 2
                or type(derivative) is not list or len(derivative) != 2):
            raise ValueError('both original charge and clock records are required')
        for row, score in zip(probabilities, derivative, strict=True):
            if (type(row) is not list or len(row) != 2
                    or any(type(p) is not float or not math.isfinite(p) or p < 0 for p in row)
                    or type(score) is not float or not math.isfinite(score) or score < 0):
                raise ValueError('invalid complete clock record')
            # This is the producer's exact binary64 operation on the retained
            # probabilities, not a comparison with a different platform's score.
            if score != (row[0]-row[1])/time:
                raise ValueError('sine derivative differs from its own original probabilities')
    except (KeyError, TypeError, ZeroDivisionError, OverflowError) as error:
        raise ValueError(path+': invalid clock derivative record') from error
    return time, probabilities, derivative


def compare(actual, expected, path='receipt'):
    """Compare complete records, propagating only the derived clock score error."""
    a_time, a_probabilities, a_scores = _clock_derivative(actual, path)
    e_time, e_probabilities, e_scores = _clock_derivative(expected, path)
    masked = []
    for value in (actual, expected):
        root = dict(value)
        root['native_clock_action'] = dict(value['native_clock_action'])
        root['native_clock_action'][_CLOCK_DERIVATIVE] = [0.0, 0.0]
        masked.append(root)
    # All raw successes, failures, rational bounds, types, source pins and
    # physical-scope flags retain the unchanged inherited comparison.
    compare_native_records(*masked, path)
    if a_time != e_time:
        raise ValueError(path+': source clock changed')
    for a_row, e_row, a_score, e_score in zip(
            a_probabilities, e_probabilities, a_scores, e_scores, strict=True):
        roundoff = 4*(sum(math.ulp(p) for p in a_row+e_row)/a_time
                      + math.ulp(a_score)+math.ulp(e_score))
        if abs(a_score-e_score) > 2*_RECORD_ERROR/a_time+roundoff:
            raise ValueError(path+': propagated clock derivative error exceeded')


def source_bindings():
    paths = set(HERE.glob('*.py'))
    queue = list(paths)
    legacy = ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim'

    def include(path):
        if path.is_file() and path not in paths:
            paths.add(path)
            queue.append(path)

    include(legacy/'oph_exact/carrier.py')
    include(legacy/'oph_exact/__init__.py')
    while queue:
        path = queue.pop()
        root = legacy if path.is_relative_to(legacy) else ROOT
        package = '.'.join(path.relative_to(root).with_suffix('').parts[:-1])
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                name = '.'*node.level+(node.module or '')
                base = importlib.util.resolve_name(name, package) if node.level else name
                names = [base]+[base+'.'+a.name for a in node.names]
            for name in names:
                root = legacy if name.split('.')[0] in ('oph_fpe', 'oph_exact') else ROOT
                segments = name.split('.')
                include(root.joinpath(*segments).with_suffix('.py'))
                for j in range(1, len(segments)+1):
                    include(root.joinpath(*segments[:j])/'__init__.py')
    paths.update(ROOT/p for p in (
        'code/source_selection_model/response.json',
        'evidence/native_geometric_source_20260925/spec.json',
        'code/m1_source_realization/receipt.json',
        'model/native_spectral_detector/contract.cpp',
        'model/native_ordered_ground/witness.npz',
        'model/native_ordered_ground/README.md',
        'model/protected_defect_catalog/PROOF.md',
        'model/spin_modular_geometry/PROOF.md',
        'model/native_spin_clock_join/README.md',
        'model/native_spin_observer/README.md',
        'model/native_gravity_identification/PROOF.md',
        'model/native_regional_modular_clock/README.md',
        'Lean/ObserverPatchHolography/EinsteinBranch/Entropy.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/SmallBall.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/Composition.lean',
        'Lean/ObserverPatchHolography/EinsteinBranch/RegisterSurface.lean'))
    paths.update((HERE/'README.md', HERE/'PROOF.md'))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def claims():
    return {'same_selected_finite_H_B_source_and_charge_channels': True,
            'exact_record_measurement_sufficiency_and_repair_blindness': True,
            'source_local_commutator_work_observable_derived': True,
            'complete_native_blank_charge_creation_energy_readout_executed': True,
            'complete_native_relative_clock_action_readout_executed': True,
            'actual_twelve_port_signed_action_functional': True,
            'mean_energy_equals_relative_H_B_action_rate': True,
            'zero_system_power_under_H_B_evolution': True,
            'zero_bare_protected_repair_energy_input': True,
            'actual_source_modular_relative_entropy_without_Gibbs_reference': True,
            'bare_repair_counts_or_control_count_determine_energy': False,
            'native_clock_is_autonomous_unread_repair': False,
            'source_H_B_clock_has_zero_native_implementation_cost': False,
            'apparatus_energy_memory_reset_heat_or_temperature_derived': False,
            'successful_charge_sectors_postselected_from_complete_output': False,
            'Pauli_energy_estimator_is_spectral_projection': False,
            'finite_conditional_mean_energy_is_particle_rest_mass': False,
            'finite_local_channel_energies_are_collective_continuum_energies': False,
            'laboratory_mass_energy_units_or_particles_identified': False}


def table(value):
    rows = ['# Finite created-sector energy readout', '',
            'The selected eight-spin source has J=1 and B=1/4. These are conditional',
            'mean excitation energies, with the same certified H_B ground reference.',
            'Both charge channels and all failures belong to the complete experiment.',
            'No row is a particle identification, a sharp spectral energy or a rest mass.', '',
            '| Internal charge | Creation weight | Unconditioned excitation moment | Mean excitation energy |',
            '|---|---:|---:|---:|']
    for row in value['native_energy_measurement']['charge_energy_readouts']:
        rows.append(f"| {row['charge']:+d} | {row['unconditioned_success_probability']:.9f} | "
                    f"{row['unconditioned_excitation_moment']:.9f} | "
                    f"{row['conditional_mean_excitation_diagnostic']:.9f} |")
    measurement = value['native_energy_measurement']
    if (F(measurement['complete_protocol_TV_error_upper']) >= F(665, 10**15)
            or F(measurement['unconditioned_energy_moment_error_upper']) >= F(284, 10**13)):
        raise ArithmeticError('finite energy table error caption is no longer valid')
    rows += ['', 'These are native estimator central values. The probability error is below',
             '6.65e-13 and the unconditioned energy-moment error is below 2.84e-11.',
             'The complete Pauli instrument measures source energy first moments.',
             'An independent controlled-clock phase measurement reproduces the two',
             'unconditioned energy moments with error below 0.000359.',
             'Exact rational intervals and native errors are retained in `results.json`.', '',
             'Under H_B evolution the relative action over duration T equals T times',
             'the corresponding mean excitation energy. System power and bare protected',
             'repair energy input are zero. Native clock implementation cost, apparatus',
             'energy and memory-reset heat are separate quantities.', '']
    return '\n'.join(rows)


def produce():
    pins = source_bindings()
    value = {'schema': 'oph.source_energy_effort.v1', 'all_pass': True,
             'repair_record_sufficiency': energy.repair_blindness(),
             'exact_same_source_work': rational(energy.work()),
             'literal_full_port_action': energy.primitive_action(),
             'native_energy_measurement': readout.protocol(),
             'native_clock_action': clock.run(),
             'native_reversible_control_cycle': readout.control_cycle(),
             'maintenance_and_relative_action': rational(clock.maintenance()),
             'claims': claims(),
             'finite_thermodynamic_attachment': {
                 'complete_system_energy_balance': True,
                 'actual_faithful_regional_modular_identity': True,
                 'physical_H_B_equals_temperature_times_modular_generator': False,
                 'Landauer_temperature_and_reservoir_supplied': False},
             'source_bindings': pins}
    value['energy_table_markdown'] = table(value)
    if source_bindings() != pins:
        raise ValueError('source changed during energy/action replay')
    return value


def instruments(value):
    measurement = value['native_energy_measurement']
    return ([measurement['creation'], measurement['selector'], value['native_clock_action'],
             value['native_reversible_control_cycle']]+measurement['settings'])


def validate(value):
    expected = produce()
    compare(value, expected)
    for current, retained in zip(instruments(expected), instruments(value), strict=True):
        a, b = (np.asarray(row['native_first_failure_weights']) for row in (current, retained))
        if np.linalg.norm(np.sqrt(a)-np.sqrt(b)) > 2e-11:
            raise ValueError('complete global failure amplitude changed')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', type=Path, const=HERE/'results.json')
    group.add_argument('--check', nargs='?', type=Path, const=HERE/'results.json')
    args = parser.parse_args()
    if args.write:
        value = produce()
        args.write.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
        args.write.with_name('ENERGY_READOUT.md').write_text(value['energy_table_markdown'])
    else:
        value = validate(strict_json.load(args.check))
        if args.check.with_name('ENERGY_READOUT.md').read_text() != value['energy_table_markdown']:
            raise ValueError('finite energy table differs from the same-source readout')
    print('PASS: complete source energy, work and signed action readouts')


if __name__ == '__main__':
    main()

