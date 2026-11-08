"""Owning exact receipt for the complete finite native-source CHSH score."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from tools import strict_json
from model.native_coherent_occupation_control import packet as controls
from . import chsh

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'


def source_bindings():
    paths = set(controls.source_bindings())
    paths.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md'))
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sorted(paths)}


def canonical(value):
    # Reuse the owning exact scalar discipline; binary floats, NaN and
    # implicit bool/int equivalence are not permitted in the receipt.
    return controls.canonical(value)


def produce():
    before = source_bindings()
    chsh.beam.source.specification.cache_clear()
    # Every full stored law is freshly evaluated, including all zero raw
    # slots. Platform-dependent rational masses are checked diagnostics;
    # only the stable analytic certificate is the immutable owner receipt.
    result = chsh.run()
    delta = F(result['one_setting_complete_halfdiamond_upper'])
    value = {
        'schema': 'oph.native_source_chsh.v1',
        'source_bindings': before,
        'exact_setting_angles': [['0', '25/64'], ['0', '-25/64'],
                                 ['25/32', '25/64'], ['25/32', '-25/64']],
        'same_pre_setting_herald_for_all_settings': True,
        'abort_outputs_independent_of_settings': [1, 1],
        'complete_unconditional_score_certified_above': '12/5',
        'exact_ideal_complete_score_interval': result['ideal_CHSH_interval'],
        'one_setting_complete_halfdiamond_upper': str(delta),
        'four_setting_score_error_upper': str(8 * delta),
        'actual_stored_score_and_gain_checked_without_normalization': True,
        'herald_native_program': result['herald_native_program'],
        'readout_native_program': result['readout_native_program'],
        'complete_records': {
            'herald_raw_slots': result['herald_raw_slots_retained'],
            'abort_slots': len(result['abort_herald_records']),
            'accepted_herald_slots': 2,
            'readout_raw_slots_per_accepted_herald': 625,
            'readout_raw_slots_per_setting': result['native_readout_raw_slots_per_setting'],
            'readout_raw_slots_all_settings': 5000,
            'all_zero_and_failure_records_retained': True,
        },
        'resources_per_complete_path': {
            'native_attempts_upper': result['maximum_native_attempts_per_path'],
            'coherent_source_pulses_upper': result['maximum_coherent_source_pulses_per_path'],
            'retained_native_records_upper': result['maximum_retained_native_records_per_path'],
            'proper_code_complement_outcomes_upper': result['source_complement_outputs_per_accepted_path'],
            'setting_input_bits': 2,
            'controller_herald_records': 1,
            'local_readout_FIFO_records': 2,
        },
        'scope': {
            'known_A0_B1_source_input_interface_retained': True,
            'declared_coherent_occupation_load_export_retained': True,
            'exact_existing_native_E_readout': True,
            'complete_native_failures_scored_and_retained': True,
            'conditional_acceptance_normalization_used': False,
            'stored_map_exact_trace_preservation_claimed': False,
            'physical_spacelike_separation_claimed': False,
            'laboratory_Bell_experiment_claimed': False,
            'continuum_field_or_particle_identification_claimed': False,
            'frozen_V2_native_word_grammar_changed': False,
        },
    }
    if source_bindings() != before:
        raise ArithmeticError('source bytes changed during complete CHSH replay')
    return controls.encode(value)


def check_payload(retained, fresh):
    if canonical(retained) != canonical(fresh):
        raise ArithmeticError('complete CHSH receipt differs from exact fresh replay')
    return True


def verify(path=OUTPUT):
    return check_payload(strict_json.load(Path(path)), produce())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--write', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        value = produce()
        OUTPUT.write_text(canonical(value), encoding='utf-8')
        check_payload(strict_json.load(OUTPUT), value)
        print('Complete native source CHSH receipt written.')
    else:
        verify()
        print('Complete native source CHSH receipt verified.')


if __name__ == '__main__':
    main()
