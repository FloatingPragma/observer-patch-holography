"""Cold replay of the corrected regional modular source operation."""
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
from . import flow


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('results.json')


def encode(value):
    if type(value) is F:
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'


def bindings():
    paths = integration_source_paths(ROOT, {
        'domain_and_proof': 'model/native_regional_modular/PROOF.md',
        'source_receipts': [flow.entropy.PRIOR, 'model/prepared_regional_entropy/results.json'],
        'independent_checks': [],
    })
    paths.update(flow.entropy.SOURCE_PATHS)
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
            for path in sorted(paths)}


def produce():
    pins = bindings()
    prior = strict_json.load(ROOT/flow.entropy.PRIOR)
    exact = flow.exact_bounds(prior)
    if exact['complete_half_diamond_error_upper'] >= F(1, 10**8):
        raise ArithmeticError('the retained corrected regional clock exceeds its complete budget')
    execution = flow.execute()
    measured = flow.diagnostic(execution)
    for key in ('complete_operator_error_diagnostic', 'complete_Gram_error', 'failure_effect_norm'):
        if not 0 <= measured[key] < 1e-9:
            raise ArithmeticError('raw native execution exceeds its separate numerical admission: '+key)
    actual_tape = [{'kind': kind, 'angle': angle, 'targets': list(targets)}
                   for kind, angle, targets in execution.word]
    if flow.literal.word_digest(actual_tape) != exact['native_word_sha256']:
        raise ArithmeticError('executed source tape differs from the certified operation')
    value = encode({
        'schema': 'oph.native_regional_modular.v1', 'all_pass': True,
        'issues': [1138, 1145], 'gates': ['CP-02', 'CP-09'],
        'interface': {
            'state': 'literal certified interacting plaquette success marginal sigma, subnormalized with every preparation record and abort retained',
            'algebra': 'complete M16 on CAR modes 0,1,2,3; exterior four field bits and arbitrary untouched reference',
            'operation': 'complete native Pauli symmetric product for exp(+is P8), retaining every first export failure',
            'clock': 'dimensionless regional reference modular parameter s=1/2 at beta=1/32; source c=1/8 U=1',
            'refinement': 'one complete 256-state interacting plaquette; 32 actuator subdivisions, no spatial limit'},
        'exact_certificate': exact,
        'floating_execution': {
            'all16_regional_input_columns_executed': measured['all_input_columns_executed'] == 16,
            'complete_operator_residual_admission': F(1, 10**9),
            'complete_Gram_residual_admission': F(1, 10**9),
            'failure_effect_norm_admission': F(1, 10**9),
            'raw_arrays_polarized_or_normalized': False,
            'floating_matrix_exponential_accuracy_used_as_proof': False},
        'composition': {
            'covers_actual_prepared_reference_by_all_input_theorem': True,
            'all256_preparation_records_and_prior_aborts_retained': True,
            'all256_preparation_words_reexecuted_by_default_replay': False,
            'full_preparation_continuation_command': 'python -m model.native_regional_modular --prepared-check',
            'regional_reference_state_half_trace_error_upper': exact['complete_half_diamond_error_upper'],
            'untouched_exterior_and_reference_complete_channel_bound': exact['complete_half_diamond_error_upper'],
            'global_reference_stationarity_asserted': False},
        'scope': {
            'corrected_regional_modular_action_constructed': True,
            'bare_regional_Gibbs_action_substituted': False,
            'regional_action_identified_with_global_modular_restriction': False,
            'geometric_diamond_or_stress_charge_identified': False,
            'physical_time_or_energy_unit_calibrated': False,
            'Einstein_or_interacting_continuum_gate_accepted': False},
        'source_bindings': pins,
    })
    if bindings() != pins:
        raise ValueError('regional modular sources changed during replay')
    return value


def report():
    pins = bindings()
    with tempfile.TemporaryDirectory(prefix='oph-native-regional-') as folder:
        path = Path(folder)/'worker.json'
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
        subprocess.run([sys.executable, '-W', 'error', '-X',
                        'pycache_prefix='+str(Path(folder)/'pycache'),
                        '-m', 'model.native_regional_modular.packet', '--worker', str(path)],
                       cwd=ROOT, env=env, check=True, capture_output=True, text=True)
        result = strict_json.load(path)
    if bindings() != pins or result['source_bindings'] != pins:
        raise ValueError('regional modular sources changed across cold replay')
    return result


def validate(value, expected):
    if canonical(value) != canonical(expected):
        raise ValueError('regional modular receipt differs from exact replay')


def prepared_check():
    """Explicitly rerun all 256 preparation histories, then continue each."""
    prepared = flow.preparation.prepare()
    continued = flow.continue_prepared(prepared)
    result = flow.prepared_diagnostic(prepared, continued)
    for key in ('actual_reference_modular_action_error_up_to_scalar',
                'regional_reference_stationarity_error', 'complete_trace_error'):
        if not 0 <= result[key] < 2e-8:
            raise ArithmeticError('complete prepared continuation failed: '+key)
    if result['all256_classical_records_retained'] is not True or result['success_renormalized'] is not False:
        raise ArithmeticError('prepared source output was conditioned or truncated')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    group.add_argument('--worker', type=Path)
    group.add_argument('--prepared-check', action='store_true')
    args = parser.parse_args()
    if args.prepared_check:
        print(json.dumps(prepared_check(), indent=2))
        return
    value = produce() if args.worker else report()
    if args.worker or args.write:
        (args.worker or args.write).write_text(canonical(value))
    else:
        validate(strict_json.load(args.check), value)
    if not args.worker:
        print(json.dumps({'all_pass': True, 'native_pulses': 37568,
                          'complete_error_upper': '1/100000000',
                          'physical_clock_calibrated': False}))


if __name__ == '__main__':
    main()
