"""Cold replay of the common interacting-field spin factor and full spectrum."""
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
from . import history, spectrum, spin


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = '9b3ff2d67929be8d43d1e0ad27514e0417c35f1d9be32d028714e5beb3dedf22'
REPLAY_ATOL = 1e-10


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/common_spin_protection/README.md',
        'source_receipts': ['model/charged_defect/results.json', 'model/field_symmetry/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/common_spin_protection/contract.json',
                  'paper/tex_fragments/COMMON_FIELD_SPIN_PROTECTION.tex'))
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
        if field not in ('energies', 'ground_energy') and (expected < 0 or actual < 0):
            raise ValueError(f'{path}: negative magnitude or probability')
        if field in ('probability', 'first_complete_probabilities', 'all_highest_weight_residues',
                     'lowest_multiplet_residue', 'probe_norm_squared') and expected > 1:
            raise ValueError(f'{path}: weight exceeds one')
        tolerance = 0 if ('admission' in field or field == 'replay_float_absolute_tolerance') else REPLAY_ATOL
        if abs(expected-actual) > tolerance:
            raise ValueError(f'{path}: numerical replay changed')
    elif type(actual) not in (str, int, bool, type(None)) or expected != actual:
        raise ValueError(f'{path}: receipt value changed')


def validate_joint_probabilities(receipt):
    value = receipt['complete_joint_history']
    first = value['first_complete_probabilities']
    joint = value['complete_joint_outcomes']
    if len(first) != 17 or len(joint) != 153:
        raise ValueError('complete history outcome count changed')
    for probabilities in (first, [row['probability'] for row in joint]):
        if (any(type(x) is not float or not math.isfinite(x) or not 0 <= x <= 1 for x in probabilities)
                or abs(math.fsum(probabilities)-1) > history.ADMISSION):
            raise ValueError('complete retained probabilities are invalid')
    if abs(math.fsum(row['probability'] for row in joint)-value['total_probability']) > 1e-13:
        raise ValueError('joint output probabilities disagree with retained total')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    if pins['model/common_spin_protection/contract.json'] != CONTRACT_SHA256:
        raise ValueError('common-spin contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.common_spin_protection.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_local_premises': spin.local_certificate(),
             'same_full_field_and_repair_checks': spin.full_certificate(),
             'complete_joint_history': history.report(),
             'complete_finite_spectrum': spectrum.report(),
             'general_block_alignment': spin.alignment_certificate(),
             'unchanged_cube_dimer_obstruction': spin.frame_obstruction(),
             'integrated_gate_accepted': False,
             'physical_comparison_count_added': 0,
             'physical_particle_identification': False,
             'primitive_or_repair_law_changed': False,
             'native_pulse_source_executed': False,
             'source_bindings': pins,
             'replay_float_absolute_tolerance': REPLAY_ATOL}
    if source_bindings() != pins:
        raise ValueError('source changed during common-spin replay')
    validate_joint_probabilities(value)
    json.dumps(value, allow_nan=False)
    return value


def report():
    """Own a fresh interpreter and empty bytecode-cache location for evidence."""
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS:
        raise ValueError('source changed since import; replay in a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-common-spin-') as temporary:
        target = Path(temporary)/'worker.json'
        environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = '1'
        environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.common_spin_protection.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode:
            raise RuntimeError('cold common-spin replay failed: '+result.stderr[-3000:])
        value = strict_json.load(target)
    if source_bindings() != pins or value['source_bindings'] != pins:
        raise ValueError('source manifest changed across cold replay')
    validate_joint_probabilities(value)
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
        validate_joint_probabilities(retained)
    else:
        path = args.worker or args.write
        with path.open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    if not args.worker:
        print('Verified common full-field spin factor, complete spectrum/history and geometry/energy boundaries.')


if __name__ == '__main__':
    main()
