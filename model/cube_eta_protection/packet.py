"""Owned cold replay of the finite three-dimensional protected eta sector."""
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
from . import cube, repair, spectrum, source


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = '78944a11a35486c34349075978de9c4f07119dcf6b39fd08809ee7b0f150369a'
REPLAY_ATOL = 1e-12


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/cube_eta_protection/README.md',
        'source_receipts': ['model/repair_spectral_sectors/results.json',
                            'model/certified_reference/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/cube_eta_protection/contract.json',
                  'paper/tex_fragments/CUBE_ETA_PROTECTION.tex'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual): raise ValueError(f'{path}: receipt type changed')
    if type(actual) is dict:
        if expected.keys() != actual.keys(): raise ValueError(f'{path}: receipt keys changed')
        for key in actual: compare(expected[key], actual[key], path+'.'+key)
    elif type(actual) is list:
        if len(expected) != len(actual): raise ValueError(f'{path}: receipt length changed')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)): compare(a, b, f'{path}[{i}]')
    elif type(actual) is float:
        if not math.isfinite(expected) or not math.isfinite(actual) or min(expected, actual) < 0:
            raise ValueError(f'{path}: invalid nonnegative finite numerical evidence')
        field = path.rsplit('.', 1)[-1]
        exact = any(word in field for word in ('admission', 'bound', 'tolerance'))
        if abs(expected-actual) > (0 if exact else REPLAY_ATOL): raise ValueError(f'{path}: numerical replay changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def validate_accounting(value):
    execution = value['native_reference_and_local_phase']
    records = execution['source_records']
    if [(r['stage'], r['pulses']) for r in records] != [('signed_tree', 5379), ('pair_copy', 24), ('local_phase', 3)]:
        raise ValueError('native source stage or pulse count changed')
    source_mass = 0.
    for row in records:
        table, pairs, count = row['all_first_export_failure_weights'], row['pair_code_event_indices'], row['pulses']
        if (len(table) != count or any(len(r) != 4 for r in table)
                or pairs != sorted(set(pairs)) or any(type(j) is not int or not 0 <= j < count for j in pairs)
                or row['source_abort_slots'] != 4*count-2*len(pairs)):
            raise ValueError('complete source abort slots changed')
        if any(type(x) is not float or not math.isfinite(x) or not 0 <= x <= 1 for r in table for x in r):
            raise ValueError('invalid source abort weight')
        if any(table[j][2:] != [0., 0.] for j in pairs): raise ValueError('source padding is not an outcome')
        mass = math.fsum(x for r in table for x in r)
        if abs(mass-row['first_export_failure_probability']) > 1e-24: raise ValueError('source abort mass changed')
        source_mass += mass
    repair_records = execution['all130_following_repair_records']
    expected = [(law, [repair.PORTS[j] for j in seam], ['native', coin, raw])
                for law in ('C', 'E') for seam in repair.SEAMS
                for coin in (0, 1) for raw in range(4 if law == 'C' else 8)]
    actual = [(r['law'], r['seam'], r['record']) for r in repair_records if r['record'][0] == 'native']
    outside = [(r['law'], r['seam']) for r in repair_records if r['record'] == ['outside']]
    if (len(repair_records) != 130 or actual != expected or outside !=
            [(law, [repair.PORTS[j] for j in seam]) for law in ('C', 'E') for seam in repair.SEAMS]):
        raise ValueError('complete finite four-block repair records changed')
    if any(type(r['probability']) is not float or not math.isfinite(r['probability']) or not 0 <= r['probability'] <= 1 for r in repair_records):
        raise ValueError('invalid repair probability')
    native_mass = math.fsum(r['probability'] for r in repair_records if r['record'][0] == 'native')
    if abs(native_mass-execution['following_repair_native_probability']) > 1e-24: raise ValueError('native repair mass changed')
    total = source_mass+native_mass+repair_records[-1]['probability']
    if abs(total-execution['complete_output_probability']) > 1e-13 or abs(total-1) > source.ADMISSION:
        raise ValueError('complete source and repair output mass changed')
    lines = value['exact_fixed_number_local_response']['eight_interacting_energy_lines']
    if len(lines) != 8 or any(r['normalized_probe_weight'] <= 0 for r in lines): raise ValueError('positive interacting response line removed')
    if abs(math.fsum(r['normalized_probe_weight'] for r in lines)-1) > 2e-12: raise ValueError('spectral response mass changed')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    if pins['model/cube_eta_protection/contract.json'] != CONTRACT_SHA256: raise ValueError('cube eta contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.cube_eta_protection.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_full_field_eta_protection': cube.certificate(),
             'exact_complete_native_four_block_repair': repair.certificate(),
             'exact_interacting_spectral_census': spectrum.polynomial_certificate(),
             'exact_fixed_number_local_response': spectrum.local_probe_certificate(),
             'exact_native_source_error_certificate': source.certificate(),
             'native_reference_and_local_phase': source.report(),
             'integrated_gate_accepted': False, 'physical_comparison_count_added': 0,
             'physical_particle_identification': False, 'primitive_or_repair_law_changed': False,
             'native_time_dependent_detector_executed': False,
             'physical_hardware_execution_or_accuracy_proved': False,
             'source_bindings': pins, 'replay_float_absolute_tolerance': REPLAY_ATOL}
    if source_bindings() != pins: raise ValueError('source changed during cube eta replay')
    validate_accounting(value)
    json.dumps(value, allow_nan=False)
    return value


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-cube-eta-') as temporary:
        target = Path(temporary)/'worker.json'; environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.cube_eta_protection.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode: raise RuntimeError('cold cube eta replay failed: '+result.stderr[-3000:])
        value = strict_json.load(target)
    if source_bindings() != pins or value['source_bindings'] != pins: raise ValueError('source manifest changed across cold replay')
    validate_accounting(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(); value = _produce() if args.worker else report()
    if args.check:
        retained = strict_json.load(args.check); compare(retained, value); validate_accounting(retained)
    else:
        with (args.worker or args.write).open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    if not args.worker: print('Verified complete3D eta protection, interacting spectrum and native local preparation.')


if __name__ == '__main__': main()
