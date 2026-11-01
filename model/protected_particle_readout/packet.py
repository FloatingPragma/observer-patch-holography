"""Owned cold replay and lossless complete-record archive for native readout."""
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
from ..cube_eta_protection import repair
from . import bounds, experiment, identification


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = 'a547a3f3b99de0fd42f1a22ba8cf21db1398e3dff5e83c3d87cfe87544f6847d'
REPLAY_ATOL = 1e-12


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/protected_particle_readout/README.md',
        'source_receipts': ['model/cube_eta_protection/results.json',
                            'model/particle_identification/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/protected_particle_readout/contract.json',
                  'paper/tex_fragments/PROTECTED_PARTICLE_READOUT.tex'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual): raise ValueError(f'{path}: receipt type changed')
    if type(actual) is dict:
        if expected.keys() != actual.keys(): raise ValueError(f'{path}: receipt keys changed')
        for key in actual: compare(expected[key], actual[key], path+'.'+key)
    elif type(actual) is list:
        if len(expected) != len(actual): raise ValueError(f'{path}: receipt length changed')
        for index, (a, b) in enumerate(zip(expected, actual, strict=True)): compare(a, b, f'{path}[{index}]')
    elif type(actual) is float:
        if not math.isfinite(expected) or not math.isfinite(actual): raise ValueError(f'{path}: nonfinite evidence')
        field = path.rsplit('.', 1)[-1]
        exact = any(word in field for word in ('admission', 'bound', 'tolerance'))
        tolerance = 0 if exact else min(REPLAY_ATOL, max(1e-42, abs(actual)*1e-8))
        if abs(expected-actual) > tolerance: raise ValueError(f'{path}: numerical evidence changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def validate_accounting(value, arrays):
    result = value['executed_native_readout']
    records = result['source_prefix_records']
    if [(row['stage'], row['pulses']) for row in records] != [
            ('tree', 5379), ('copy', 24), ('pointer', 19),
            ('field_first', 337920), ('field_second', 337920)]:
        raise ValueError('native source stages or counts changed')
    settings = result['pointer_settings']
    if len(settings) != 2 or [s['phase_over_pi'] for s in settings] != ['0', '1/2']:
        raise ValueError('both pointer quadratures required')
    all_records = records+[s['analyzer_source_record'] for s in settings]
    keys = {r['failure_weight_array'] for r in all_records}
    keys.update(f'analyzer_{j}_full_success_state' for j in (0, 1))
    if set(arrays) != keys: raise ValueError('complete record archive key set changed')
    for row in all_records:
        table = arrays[row['failure_weight_array']]
        if table.shape != (row['pulses'], 4) or table.dtype != np.dtype('float64'):
            raise ValueError('source failure table shape/type changed')
        if not np.isfinite(table).all() or np.any(table < 0): raise ValueError('invalid source failure weights')
        if row['native_abort_slots']+row['padding_slots'] != 4*row['pulses']:
            raise ValueError('source failure/padding coverage changed')
        mass = float(table.sum())
        if abs(mass-row['first_failure_probability']) > row['failure_sum_admission']+1e-40:
            raise ValueError('source first-failure mass changed')
    native = result['all130_repair_records']
    expected = [(law, [repair.PORTS[j] for j in seam], ['native', coin, raw])
                for law in ('C', 'E') for seam in repair.SEAMS
                for coin in (0, 1) for raw in range(4 if law == 'C' else 8)]
    actual = [(r['law'], r['native_seam'], r['record']) for r in native if r['record'][0] == 'native']
    outside = [(r['law'], r['native_seam']) for r in native if r['record'] == ['outside']]
    if (len(native) != 130 or actual != expected or outside !=
            [(law, [repair.PORTS[j] for j in seam]) for law in ('C', 'E') for seam in repair.SEAMS]):
        raise ValueError('complete native/outside repair records changed')
    if any(type(r['probability']) is not float or not math.isfinite(r['probability']) or
           not 0 <= r['probability'] <= 1 for r in native): raise ValueError('invalid native repair weight')
    native_mass = math.fsum(r['probability'] for r in native if r['record'][0] == 'native')
    if abs(native_mass-result['terminal_native_repair_probability']) > 1e-40:
        raise ValueError('native terminal mass changed')
    prefix = math.fsum(r['first_failure_probability'] for r in records)+native_mass
    for index, setting in enumerate(settings):
        state = arrays[f'analyzer_{index}_full_success_state']
        if state.shape != (131072, 1) or state.dtype != np.dtype('complex128') or not np.isfinite(state).all():
            raise ValueError('full final field/pointer output changed')
        probabilities = setting['pointer_probabilities']
        if len(probabilities) != 6 or any(not math.isfinite(x) or not 0 <= x <= 1 for x in probabilities):
            raise ValueError('complete six-outcome pointer export changed')
        calculated = [float(np.vdot(out, out).real) for out in experiment.pointer_outputs(state)]
        if max(abs(a-b) for a, b in zip(calculated, probabilities, strict=True)) > REPLAY_ATOL:
            raise ValueError('pointer marginal disagrees with retained field outputs')
        total = prefix+setting['analyzer_source_record']['first_failure_probability']+math.fsum(probabilities)
        if abs(total-setting['complete_probability']) > REPLAY_ATOL or abs(total-1) > result['numerical_simulation_admission']:
            raise ValueError('complete unconditioned probability changed')
    if (result['resolved_all_eight_frequencies_from_these_two_settings'] is not False or
            result['absolute_energy_zero_determined_by_Ramsey'] is not False or
            result['observed_inelastic_cyclic_dimension'] != 8):
        raise ValueError('operational identification scope changed')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    if pins['model/protected_particle_readout/contract.json'] != CONTRACT_SHA256: raise ValueError('readout contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.protected_particle_readout.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_native_error_certificate': bounds.source_certificate(),
             'executed_native_readout': experiment.report(),
             'prospective_spectral_resolution': identification.spectral_resolution_certificate(),
             'prospective_sampling_budget': identification.prospective_sampling_budget(),
             'source_bindings': pins, 'replay_float_absolute_tolerance': REPLAY_ATOL,
             'integrated_gate_accepted': False, 'physical_particle_identification': False,
             'physical_hardware_execution_or_accuracy_proved': False,
             'primitive_or_field_or_repair_law_changed': False}
    arrays = experiment.arrays()
    if source_bindings() != pins: raise ValueError('source changed during native readout')
    validate_accounting(value, arrays)
    json.dumps(value, allow_nan=False)
    return value, arrays


def write_pair(target, value, arrays):
    archive = target.with_suffix('.npz')
    if target.exists() or archive.exists(): raise FileExistsError('retain existing receipt and archive')
    with archive.open('xb') as stream: np.savez_compressed(stream, **arrays)
    value = dict(value, complete_record_archive={'filename': archive.name,
                      'sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
                      'format': 'lossless NumPy NPZ; no object arrays or pickle'})
    with target.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def read_pair(target):
    value = strict_json.load(target); metadata = value.pop('complete_record_archive')
    if Path(metadata['filename']).name != metadata['filename']: raise ValueError('archive path must be a local filename')
    archive = target.parent/metadata['filename']
    if hashlib.sha256(archive.read_bytes()).hexdigest() != metadata['sha256']: raise ValueError('archive hash changed')
    with np.load(archive, allow_pickle=False) as stored:
        arrays = {name: stored[name] for name in stored.files}
    validate_accounting(value, arrays)
    return value, arrays


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-protected-readout-') as temporary:
        target = Path(temporary)/'worker.json'; environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.protected_particle_readout.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode: raise RuntimeError('cold native readout failed: '+result.stderr[-4000:])
        value, arrays = read_pair(target)
    if source_bindings() != pins or value['source_bindings'] != pins: raise ValueError('source changed across cold readout replay')
    return value, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(); value, arrays = _produce() if args.worker else report()
    if args.check:
        retained, stored = read_pair(args.check); compare(retained, value)
        if stored.keys() != arrays.keys(): raise ValueError('archive keys changed')
        for name, actual in arrays.items():
            atol, rtol = (1e-42, 1e-8) if name.endswith('failure_weights') else (REPLAY_ATOL, 0)
            if stored[name].shape != actual.shape or not np.allclose(stored[name], actual, atol=atol, rtol=rtol):
                raise ValueError('archive replay changed: '+name)
    else: write_pair(args.worker or args.write, value, arrays)
    if not args.worker: print('Verified full17-bit native protected response, complete records and separate prospective resolution.')


if __name__ == '__main__': main()
