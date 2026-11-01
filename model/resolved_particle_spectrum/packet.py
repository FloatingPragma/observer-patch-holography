"""Owned cold replay of certified eight-line spectral weight resolution."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from ..protected_particle_readout.packet import compare
from . import budgets, certified_power as cp, experiment, factorization, inversion


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
CONTRACT_SHA256 = '69ef2d5a06eaa2bec03f63e9c604bf4bc4fdcd1056a4cd6418550359397372b2'


def source_bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/resolved_particle_spectrum/README.md',
        'source_receipts': ['model/protected_particle_readout/results.json'],
        'independent_checks': [],
    })
    paths.update(('model/resolved_particle_spectrum/contract.json',
                  'paper/tex_fragments/RESOLVED_PARTICLE_SPECTRUM.tex'))
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


_IMPORTED_BINDINGS = source_bindings()


def validate(value):
    result = value['executed_accelerated_resolution']
    if (result['actual_raw_native_long_word_executed'] is not False or
            result['unknown_frequency_support_inferred'] is not False or
            result['physical_particle_species_assigned'] != 0):
        raise ValueError('accelerated resolution scope changed')
    settings, matrices = result['settings'], result['retained_dyadic_power_matrices']
    if len(settings) != 8 or len(matrices) != 8: raise ValueError('all eight times and powers required')
    zero, quadrature = [], []
    for index, (setting, matrix) in enumerate(zip(settings, matrices, strict=True)):
        if setting['time'] != 39*index or matrix['time'] != setting['time']:
            raise ValueError('sampling time grid changed')
        if len(setting['pointer_settings']) != 2: raise ValueError('both quadratures required')
        certificate = setting['exact_word_power']
        actual = cp.CertifiedUnitary(matrix['real_integer_numerators'], matrix['imag_integer_numerators'],
                                     matrix['precision_bits'], F(certificate['operator_error_upper']),
                                     certificate['rounded_matrix_multiplications'])
        if actual.dimension != 36 or actual.certificate() != certificate:
            raise ValueError('retained exact dyadic power changed')
        for phase, row in enumerate(setting['pointer_settings']):
            if row['phase_over_pi'] != ('0' if phase == 0 else '1/2') or row['trials'] != experiment.TRIALS:
                raise ValueError('pointer phase or sample count changed')
            count = row['pointer_zero_count']
            if (type(count) is not int or not 0 <= count <= experiment.TRIALS or
                    count+row['all_other_complete_outcomes_count'] != experiment.TRIALS):
                raise ValueError('complete unconditional sample count changed')
            if F(row['total_systematic_probability_bound']) >= experiment.SYSTEMATIC:
                raise ValueError('systematic probability budget exceeded')
            if F(row['sampling_probability_binary64']) != F(row['sampling_probability_exact_binary_rational']):
                raise ValueError('binary sampler probability changed')
            if abs(F(count, experiment.TRIALS)-F(row['sampling_probability_exact_binary_rational'])) > experiment.SAMPLING:
                raise ValueError('retained sample misses the stated simultaneous event')
            (zero if phase == 0 else quadrature).append(count)
    inferred = inversion.reconstruct(zero, quadrature, experiment.TRIALS,
                                     experiment.SYSTEMATIC, experiment.SAMPLING)
    if inferred != result['weight_reconstruction'] or not inferred['all_eight_weights_strictly_positive']:
        raise ValueError('eight positive count-derived weights changed')
    if result['total_simulated_reset_trials'] != 16*experiment.TRIALS:
        raise ValueError('total trial count changed')


def _produce():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    if pins['model/resolved_particle_spectrum/contract.json'] != CONTRACT_SHA256: raise ValueError('resolution contract changed')
    contract = strict_json.load(HERE/'contract.json')
    value = {'schema': 'oph.resolved_particle_spectrum.receipt.v1',
             'campaign_id': contract['campaign_id'], 'scope': contract['scope'],
             'contract_sha256': CONTRACT_SHA256,
             'exact_source_word_factorization': factorization.certificate(),
             'complete_native_error_and_resources': budgets.certificate(),
             'certified_rational_inverse': inversion.matrix_certificate(),
             'executed_accelerated_resolution': experiment.execution(),
             'integrated_gate_accepted': False,
             'physical_hardware_execution_or_accuracy_proved': False,
             'primitive_or_field_or_repair_law_changed': False,
             'source_bindings': pins}
    if source_bindings() != pins: raise ValueError('source changed during spectral resolution')
    validate(value); json.dumps(value, allow_nan=False)
    return value


def report():
    pins = source_bindings()
    if pins != _IMPORTED_BINDINGS: raise ValueError('source changed since import; use a fresh interpreter')
    with tempfile.TemporaryDirectory(prefix='oph-resolved-spectrum-') as temporary:
        target = Path(temporary)/'worker.json'; environment = dict(os.environ)
        environment['OPENBLAS_NUM_THREADS'] = environment['OMP_NUM_THREADS'] = '1'
        result = subprocess.run([sys.executable, '-W', 'error', '-X', 'pycache_prefix='+str(Path(temporary)/'pycache'),
                                 '-m', 'model.resolved_particle_spectrum.packet', '--worker', str(target)],
                                cwd=ROOT, env=environment, capture_output=True, text=True, check=False)
        if result.returncode: raise RuntimeError('cold spectral resolution failed: '+result.stderr[-4000:])
        value = strict_json.load(target)
    if source_bindings() != pins or value['source_bindings'] != pins: raise ValueError('source changed across cold replay')
    validate(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--write', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--check', type=Path, nargs='?', const=OUTPUT)
    actions.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(); value = _produce() if args.worker else report()
    if args.check:
        retained = strict_json.load(args.check); compare(retained, value); validate(retained)
    else:
        with (args.worker or args.write).open('x', encoding='utf-8', newline='\n') as stream:
            stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    if not args.worker: print('Verified eight positive count-derived weights with exact word acceleration and complete native error bounds.')


if __name__ == '__main__': main()
