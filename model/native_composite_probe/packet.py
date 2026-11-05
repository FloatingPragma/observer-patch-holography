"""Cold complete composite detector replay with an authenticated output archive."""
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
from . import circuits, detection

ROOT = detection.ROOT
HERE = Path(__file__).parent
OUTPUT = HERE/'results.json'
ATOL = 2e-11


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/native_composite_probe/PROOF.md',
        'source_receipts': ['model/centered_native_preparation/results.json',
                            'model/charged_composite_search/results.json'],
        'independent_checks': []})
    paths.update({'model/centered_native_preparation/results.npz',
                  'model/centered_particle_spectrum/witness.npz',
                  'model/charged_composite_search/witness.npz',
                  'model/native_spectral_detector/contract.cpp'})
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def validate(value, arrays):
    if value.get('schema') != 'oph.native_composite_probe.v1' or value.get('all_pass') is not True:
        raise ValueError('invalid native composite receipt')
    if value['control_certificate'] != circuits.certificate():
        raise ValueError('literal source-control certificate changed')
    shapes = {'complete_pointer_field_reference_amplitudes': (32, 256, 256),
              'first_export_failure_weights': (1008, 4), 'pointer_probabilities': (32,),
              'native_pointer_probabilities': (7776,)}
    if set(arrays) != set(shapes):
        raise ValueError('complete native archive keys changed')
    for name, shape in shapes.items():
        a = arrays[name]
        dtype = np.complex128 if name.endswith('amplitudes') else np.float64
        if a.shape != shape or a.dtype != dtype or not np.isfinite(a).all():
            raise ValueError('complete native archive carrier changed: '+name)
        if dtype == np.float64 and np.min(a) < 0:
            raise ValueError('negative record probability')
    logical = np.sum(abs(arrays['complete_pointer_field_reference_amplitudes'])**2, axis=(1, 2))
    if not np.allclose(logical, arrays['pointer_probabilities'], rtol=0, atol=ATOL):
        raise ValueError('pointer probabilities disagree with retained amplitudes')
    if not np.array_equal(arrays['native_pointer_probabilities'],
                          detection.native_pointer_probabilities(arrays['pointer_probabilities'])):
        raise ValueError('incomplete native six-outcome pointer table')
    pairs = np.array([event['kind'] == 'J' for event in circuits.word()])
    if np.any(arrays['first_export_failure_weights'][pairs, 2:] != 0):
        raise ValueError('padding is not a native outcome')
    run = value['execution']
    for key, expected in {'physical_species_or_mass_identification': False,
                          'native_energy_clock_or_resolved_138_lines_executed': False,
                          'inherited_preparation_reexecuted': False,
                          'inherited_source_custody_verified': True,
                          'success_postselection_or_renormalization': False,
                          'every_reference_column_is_coherent_not_a_classical_sample': True}.items():
        if run.get(key) is not expected:
            raise ValueError('detector scope changed: '+key)
    for key, expected in {'native_pointer_outcomes': 7776, 'logical_pointer_labels': 32,
                          'native_failure_slots': 3492, 'padding_slots': 540,
                          'complete_joint_amplitudes': 2097152}.items():
        if type(run.get(key)) is not int or run[key] != expected:
            raise ValueError('complete native outcome accounting changed')
    if abs(run['success_pointer_probability']-logical[0]) > ATOL:
        raise ValueError('success probability changed')
    if abs(run['first_export_failure_probability']-arrays['first_export_failure_weights'].sum()) > 1e-38:
        raise ValueError('first-failure mass changed')
    if run['preparation_and_detector_probability_residual'] > detection.execution.ADMISSION:
        raise ValueError('full preparation and detector lost probability')
    for key in ('native_to_ideal_complete_output_error_diagnostic',
                'complete_probability_residual', 'largest_chunk_probability_residual'):
        if not 0 <= run[key] <= detection.execution.ADMISSION:
            raise ValueError('native numerical diagnostic exceeds its admission')
    if not 0 <= run['first_export_failure_probability'] <= 1e-24:
        raise ValueError('native failure diagnostic exceeds its numerical admission')
    if (type(run['rounded_negative_slot_count']) is not int
            or not 0 <= run['rounded_negative_slot_count'] <= 8*1008*4
            or not 0 <= run['rounded_negative_slot_mass'] <= 1e-30):
        raise ValueError('native floating roundoff exceeds its admission')
    if value['source_bindings'] != bindings():
        raise ValueError('source dependency changed')


def produce():
    pins = bindings()
    result, arrays = detection.report()
    value = {'schema': 'oph.native_composite_probe.v1', 'all_pass': True,
             'issues': [1141, 1143], 'gates': ['CP-05', 'CP-07'],
             'all_pass_scope': 'finite native composite filter and complete records; no resolved energy clock',
             **result, 'source_bindings': pins}
    validate(value, arrays)
    return value, arrays


def write_pair(path, value, arrays):
    archive = path.with_suffix('.npz')
    np.savez_compressed(archive, **arrays)
    value = dict(value, archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
    path.write_text(canonical(value))


def read_pair(path):
    value = strict_json.load(path)
    archive = path.with_suffix('.npz')
    if value['archive_sha256'] != hashlib.sha256(archive.read_bytes()).hexdigest():
        raise ValueError('native output archive hash changed')
    with np.load(archive, allow_pickle=False) as stored:
        arrays = {key: stored[key] for key in stored.files}
    validate(value, arrays)
    return value, arrays


def compare(left, right, path='receipt'):
    if type(left) is not type(right):
        raise ValueError(path+': type changed')
    if isinstance(left, dict):
        if left.keys() != right.keys():
            raise ValueError(path+': keys changed')
        for key in left:
            # Both receipts are validated independently. The count of
            # tiny negative floating roundoff slots depends on host BLAS.
            if key not in ('archive_sha256', 'rounded_negative_slot_count'):
                compare(left[key], right[key], path+'.'+key)
    elif isinstance(left, list):
        if len(left) != len(right):
            raise ValueError(path+': length changed')
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            compare(a, b, path+str(i))
    elif type(left) is float:
        # Floating diagnostics use the declared replay precision. Rational
        # theorem bounds, tape hashes and carrier counts stay exact.
        if not math.isfinite(left) or not math.isfinite(right) or abs(left-right) > ATOL:
            raise ValueError(path+': numerical replay changed')
    elif left != right:
        raise ValueError(path+': value changed')


def report():
    pins = bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-composite-') as temp:
        output = Path(temp)/'results.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        worker = subprocess.run([sys.executable, '-W', 'error', '-X',
                                 'pycache_prefix='+str(Path(temp)/'cache'),
                                 '-m', 'model.native_composite_probe.packet', '--worker', str(output)],
                                cwd=ROOT, env=env, capture_output=True, text=True)
        if worker.returncode:
            raise RuntimeError('cold native composite replay failed:\n'+worker.stderr)
        value, arrays = read_pair(output)
    if pins != bindings():
        raise ValueError('source changed during cold native replay')
    return value, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    group.add_argument('--worker', type=Path)
    args = parser.parse_args()
    if args.worker:
        write_pair(args.worker, *produce())
        return
    value, arrays = report()
    if args.write:
        write_pair(OUTPUT, value, arrays)
    else:
        retained, old = read_pair(OUTPUT)
        compare(retained, value)
        for key in old:
            if not np.allclose(old[key], arrays[key], rtol=0, atol=ATOL):
                raise ValueError('full native array changed: '+key)
    print(json.dumps({'all_pass': True, 'native_pulses': 1008,
                      'raw_pointer_records': 7776, 'resolved_energy_clock': False}))


if __name__ == '__main__':
    main()
