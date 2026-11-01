"""Owned cold replay of exact active/dark algebra and complete repair response."""
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

import numpy as np

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import algebra, response


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = '9085cc3ac7a1eb22ad18eb9f99acc7f0f0ec8a3b77a4c4f72efb5e48a25f6afd'
REPLAY_ATOL = 1e-10
SIGNED_FIELDS = {'all90_input_and90_output_energies', 'mean_K_transfer_per_energy_input',
                 'numerical_display_energies', 'numerical_energies'}


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/repair_spectral_sectors/README.md',
        'source_receipts': ['model/common_spin_protection/results.json',
                            'model/charged_defect/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/repair_spectral_sectors/contract.json',
                  'paper/tex_fragments/REPAIR_SPECTRAL_SECTORS.tex'))
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
        if field not in SIGNED_FIELDS and (expected < 0 or actual < 0):
            raise ValueError(f'{path}: negative magnitude or probability')
        tolerance = 0 if ('admission' in field or field == 'replay_float_absolute_tolerance') else REPLAY_ATOL
        if abs(expected-actual) > tolerance:
            raise ValueError(f'{path}: numerical replay changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def validate_accounting(receipt):
    active = receipt['complete_active_energy_response']
    matrix = np.asarray(active['complete_energy_transition_matrix'])
    outside = np.asarray(active['outside_energy_transition_matrix'])
    energy = np.asarray(active['all90_input_and90_output_energies'])
    if matrix.shape != (90, 90) or outside.shape != matrix.shape or energy.shape != (90,):
        raise ValueError('the complete active response dimension changed')
    for probabilities in (matrix, outside):
        if (probabilities.dtype.kind != 'f' or not np.isfinite(probabilities).all()
                or np.any(probabilities < 0) or np.any(probabilities > 1)):
            raise ValueError('invalid retained active probabilities')
    if (np.any(outside > matrix+1e-14)
            or np.max(abs(matrix.sum(axis=0)-1)) > response.ADMISSION
            or np.max(abs(matrix.sum(axis=1)-1)) > response.ADMISSION):
        raise ValueError('active complete output accounting failed')
    if not np.isfinite(energy).all():
        raise ValueError('nonfinite energy chart')
    drift = energy@matrix-energy
    spread = np.sum(matrix*(energy[:, None]-energy[None, :])**2, axis=0)
    if (np.max(abs(drift-np.asarray(active['mean_K_transfer_per_energy_input']))) > response.ADMISSION
            or np.max(abs(spread-np.asarray(active['mean_squared_energy_change_per_input']))) > response.ADMISSION
            or abs(drift[0]-active['lowest_energy_input_mean_transfer']) > response.ADMISSION):
        raise ValueError('active energy moments differ from the complete response')
    dark = receipt['complete_dark_spectral_response']
    rows = dark['complete_stopping_history_outputs']
    expected = [(attempt, law, list(seam), ['native', coin, raw])
                for attempt, (law, seam) in enumerate((('E', (0, 2)), ('C', (0, 1)), ('E', (1, 2))))
                for coin in (0, 1) for raw in range(8 if law == 'E' else 4)]
    expected.append((3, 'stop', [], ['outside']))
    if (len(rows) != 41 or [(r['attempt'], r['law'], r['seam'], r['record']) for r in rows] != expected):
        raise ValueError('complete first-native/outside stopping records changed')
    weights = [r['probability'] for r in rows]
    if (any(type(x) is not float or not math.isfinite(x) or not 0 <= x <= 1 for x in weights)
            or abs(math.fsum(weights)-1) > response.ADMISSION
            or abs(math.fsum(weights[:-1])-dark['native_outcome_probability']) > 1e-24
            or weights[-1] != dark['final_outside_outcome_probability']):
        raise ValueError('complete stopping-history probability accounting failed')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    if pins['model/repair_spectral_sectors/contract.json'] != CONTRACT_SHA256:
        raise ValueError('repaired spectral sector contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.repair_spectral_sectors.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_active_cyclicity': algebra.cyclic_certificate(),
             'actual_complete_repair_corner': algebra.repair_corner_certificate(),
             'exact_dark_eta_spectrum': algebra.dark_certificate(),
             'complete_energy_transfer_identity': response.complete_work_identity(),
             'complete_active_energy_response': response.active_response(),
             'complete_dark_spectral_response': response.dark_response(),
             'integrated_gate_accepted': False,
             'physical_comparison_count_added': 0,
             'physical_particle_identification': False,
             'primitive_or_repair_law_changed': False,
             'native_pulse_source_executed': False,
             'all_Fock_common_algebra_classified': False,
             'controller_work_or_heat_identified': False,
             'source_bindings': pins,
             'replay_float_absolute_tolerance': REPLAY_ATOL}
    if source_bindings() != pins:
        raise ValueError('source changed during repaired spectral replay')
    validate_accounting(value)
    json.dumps(value, allow_nan=False)
    return value


def report():
    """Own a fresh interpreter and fresh bytecode location for every replay."""
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-repair-spectrum-') as temporary:
        target = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = '1'
        environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.repair_spectral_sectors.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold repaired spectral replay failed: '+result.stderr[-3000:])
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
        print('Verified active/dark repaired spectral sectors and complete energy transfer.')


if __name__ == '__main__':
    main()
