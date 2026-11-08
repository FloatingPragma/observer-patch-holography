"""Exact common-state clock certificates and strict live source bindings."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import algebra, bounds, energy

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).parent


def source_bindings():
    paths = set(HERE.glob('*.py')) | {HERE/'README.md'}
    paths.update(ROOT/path for path in (
        'model/__init__.py', 'tools/strict_json.py',
        'model/native_spin_observer/README.md',
        'model/native_spin_observer/bounds.py',
        'model/native_spin_observer/precision.py',
        'model/native_ordered_ground/README.md',
        'model/native_ordered_ground/source.py',
        'model/spin_current_spectrum/README.md',
        'model/spin_current_spectrum/algebra.py',
    ))
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def produce():
    pins = source_bindings()
    finite = bounds.finite_history(F(1, 100000), 8, [F(1, 2), F(1, 2)],
                                  native_error=bounds.spin_port_error(10000, 60))
    local = bounds.local_clock_plan(support_radius=1, support_sites=2,
                                   coupling=F(1), duration=F(1, 4),
                                   interventions=3, error=F(1, 10000))
    large = bounds.local_clock_plan(support_radius=10**20, support_sites=1,
                                   coupling=F(1), duration=F(10**30),
                                   interventions=2, error=F(1, 10**80))
    result = bounds.encode({
        'schema': 'oph.native_spin_clock_join.v1',
        'all_pass': True,
        'complete_finite_history_example': finite,
        'local_common_state_clock_plan': local,
        'growing_time_clock_plan': large,
        'finite_ground_energy_bounds': bounds.ground_energy_bound(F(3, 8), 2, 1),
        'exact_two_spin_control': algebra.controls(),
        'exact_energy_continuity': energy.exact_checks(),
        'claims': {
            'same_prepared_state_complete_clock_comparison': True,
            'all_failure_and_zero_records_retained': True,
            'arbitrary_spectator_reference_allowed': True,
            'finite_local_energy_continuity_and_weighted_identity': True,
            'finite_controlled_phase_compensation_required_if_recentered': True,
            'local_bound_applies_to_global_energy_spectrometer': False,
            'finite_field_ground_is_zero_field_vacuum': False,
            'effective_infinite_ground_selection_modulus': False,
            'native_pulse_tape_executed_here': False,
            'sharp_mass_shell_or_relativistic_rest_mass': False,
            'interacting_relativistic_continuum': False,
        },
        'source_bindings': pins,
    })
    if pins != source_bindings():
        raise RuntimeError('source changed during exact clock verification')
    return result


def validate(value):
    expected = produce()
    # Canonical serialization also distinguishes bool from integer and rejects
    # forged extra fields, modified failure weights and weakened nonclaims.
    if json.dumps(value, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
        raise ValueError('clock receipt differs from exact replay or live sources')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check', action='store_true')
    group.add_argument('--write', action='store_true')
    args = parser.parse_args()
    path = HERE/'results.json'
    if args.write:
        path.write_text(json.dumps(produce(), indent=2, sort_keys=True, allow_nan=False)+'\n')
    else:
        validate(strict_json.load(path))
    print('PASS: complete common-state spin clock comparison')


if __name__ == '__main__':
    main()
