"""Owned cold replay of the complete native dark-state Ramsey experiment."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import exact, experiment


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = '9afe9c88ec238d141fd23496d3e27a85a1a5514c573fa219cd682df6bced67b8'
REPLAY_ATOL = 1e-10


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/dark_spectral_source/README.md',
        'source_receipts': ['model/repair_spectral_sectors/results.json',
                            'model/certified_reference/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/dark_spectral_source/contract.json',
                  'paper/tex_fragments/DARK_SPECTRAL_SOURCE.tex'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(f'{path}: receipt type changed')
    if type(actual) is dict:
        if expected.keys() != actual.keys():
            raise ValueError(f'{path}: receipt keys changed')
        for key in actual:
            compare(expected[key], actual[key], path+'.'+key)
    elif type(actual) is list:
        if len(expected) != len(actual):
            raise ValueError(f'{path}: receipt length changed')
        for index, (left, right) in enumerate(zip(expected, actual, strict=True)):
            compare(left, right, f'{path}[{index}]')
    elif type(actual) is float:
        field = path.rsplit('.', 1)[-1].split('[', 1)[0]
        if not math.isfinite(expected) or not math.isfinite(actual):
            raise ValueError(f'{path}: nonfinite numerical value')
        if field not in ('K_energy', 'nominal_ideal_characteristic_imag') and (expected < 0 or actual < 0):
            raise ValueError(f'{path}: negative magnitude or probability')
        tolerance = 0 if ('admission' in field or field == 'replay_float_absolute_tolerance') else REPLAY_ATOL
        if abs(expected-actual) > tolerance:
            raise ValueError(f'{path}: numerical replay changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def _probabilities(values):
    if any(type(x) is not float or not math.isfinite(x) or not 0 <= x <= 1 for x in values):
        raise ValueError('invalid complete-output probability')


def validate_source_record(record):
    table = record['all_first_export_failure_weights']
    pulses = record['pulses']
    pairs = record['pair_code_event_indices']
    if (type(pulses) is not int or pulses < 1 or len(table) != pulses
            or any(len(row) != 4 for row in table)
            or pairs != sorted(set(pairs))
            or any(type(j) is not int or not 0 <= j < pulses for j in pairs)
            or record['source_abort_slots'] != 4*pulses-2*len(pairs)
            or record['padding_slots_in_four_column_table'] != 2*len(pairs)):
        raise ValueError('complete source event/abort-slot count changed')
    for row in table:
        _probabilities(row)
    if any(table[j][2:] != [0., 0.] for j in pairs):
        raise ValueError('padding cannot be represented as a source outcome')
    mass = math.fsum(x for row in table for x in row)
    if abs(mass-record['first_export_failure_probability']) > 1e-24:
        raise ValueError('source abort mass differs from all retained slots')
    return mass


def validate_accounting(receipt):
    history = receipt['complete_native_Ramsey_history']
    prefix = history['common_prefix_source_records']
    if [row['stage'] for row in prefix] != ['preparation', 'field_first', 'field_second']:
        raise ValueError('complete source prefix changed')
    prefix_mass = sum(validate_source_record(row) for row in prefix)
    repairs = history['all40_native_repair_records']
    expected = [(law, list(seam), ['native', coin, raw])
                for law, seam in (('E', (0, 1)), ('C', (1, 2)), ('E', (0, 2)))
                for coin in (0, 1) for raw in range(8 if law == 'E' else 4)]
    if [(row['law'], row['seam'], row['record']) for row in repairs] != expected:
        raise ValueError('the complete native repair record family changed')
    _probabilities([row['probability'] for row in repairs])
    repair_mass = math.fsum(row['probability'] for row in repairs)
    if abs(repair_mass-history['repair_native_probability']) > 1e-24:
        raise ValueError('repair mass differs from retained native records')
    settings = history['pointer_settings']
    if [row['phase_over_pi'] for row in settings] != ['0', '1/2']:
        raise ValueError('both Ramsey quadratures are required')
    for setting in settings:
        pointer = setting['pointer_probabilities']
        if setting['pointer_records'] != ['zero', 'one', 'abort2', 'abort3', 'abort4', 'abort5'] or len(pointer) != 6:
            raise ValueError('complete pointer export changed')
        _probabilities(pointer)
        analyzer_mass = validate_source_record(setting['analyzer_source'])
        total = prefix_mass+repair_mass+analyzer_mass+math.fsum(pointer)
        if abs(total-setting['complete_probability']) > 1e-13 or abs(total-1) > experiment.execute.ADMISSION:
            raise ValueError('complete Ramsey output mass changed')
    spectrum = receipt['exact_source_accessible_spectrum']
    lines, groups = spectrum['all15_mode_pairs'], spectrum['all13_distinct_energy_groups']
    if len(lines) != 15 or len(groups) != 13 or any(row['weight'] <= 0 for row in lines):
        raise ValueError('a positive dark spectral response disappeared')
    _probabilities([row['weight'] for row in lines])
    if abs(math.fsum(row['weight'] for row in lines)-1) > 1e-12:
        raise ValueError('complete spectral response is not normalized')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    if pins['model/dark_spectral_source/contract.json'] != CONTRACT_SHA256:
        raise ValueError('dark source contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.dark_spectral_source.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_sparse_preparation': exact.preparation_certificate(),
             'exact_source_accessible_spectrum': exact.spectral_certificate(),
             'exact_source_and_field_error_certificate': experiment.json_rationals(experiment.source_certificate()),
             'complete_native_Ramsey_history': experiment.report(),
             'integrated_gate_accepted': False,
             'physical_comparison_count_added': 0,
             'physical_particle_identification': False,
             'primitive_or_repair_law_changed': False,
             'native_source_preparation_and_detector_executed': True,
             'physical_hardware_execution_or_accuracy_proved': False,
             'source_bindings': pins,
             'replay_float_absolute_tolerance': REPLAY_ATOL}
    if source_bindings() != pins:
        raise ValueError('source changed during dark source replay')
    validate_accounting(value)
    json.dumps(value, allow_nan=False)
    return value


def report():
    """Own a fresh interpreter and empty bytecode cache for each replay."""
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-dark-source-') as temporary:
        target = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = '1'
        environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.dark_spectral_source.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold dark-source replay failed: '+result.stderr[-3000:])
        value = strict_json.load(target)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source manifest changed across cold replay')
    validate_accounting(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    value = _produce() if args.worker else report()
    if args.check:
        retained = strict_json.load(args.check)
        compare(retained, value)
        validate_accounting(retained)
    else:
        with (args.worker or args.write).open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    if not args.worker:
        print('Verified sparse native dark preparation, all-line Ramsey detection and complete outputs.')


if __name__ == '__main__':
    main()
