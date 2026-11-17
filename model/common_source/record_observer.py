"""Replay optimal recorded native memory, protected fields and information.

The retained evidence includes every complex field/reference output, including
zero outcomes and both selector timeouts. Native tensor outputs are represented
by the proved local factorization, not by a simulated 512**6 array. Source
execution diagnostics are numerical guards, not outward or hardware bounds.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from ..native_equilibration import common_variation
from . import open_observer, record_capacity, record_code, record_field, record_information


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = tuple(sorted(set((*common_variation.SOURCE_PATHS, *open_observer.SOURCE_PATHS,
    'model/common_source/record_observer.py', 'model/common_source/record_capacity.py',
    'model/common_source/record_code.py', 'model/common_source/record_field.py',
    'model/common_source/record_information.py', 'model/common_source/slow_geometry.py',
    'model/observer/record_recovery.py', 'model/repair_scheduler/coins.py',
    'requirements.txt',
    'model/common_source/test_record_observer.py', 'model/common_source/test_record_capacity.py',
    'model/common_source/test_record_code.py', 'model/common_source/test_record_field.py',
    'model/common_source/test_record_information.py',
    'model/common_source/ENCODED_MEMORY.md',
    'evidence/native_geometric_source_20260925/spec.json',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
))))
NUMERICAL_COMPARISON_TOLERANCE = 2e-9
FIELD_OUTPUT_OPERATOR_GATE = 1e-9
MACRO_GATES = {'raw_unitarity_defect': 1e-10, 'projection_change': 5.1e-11,
               'restored_unitarity_defect': 2e-13}


def _ideal_field_blocks():
    """Independent complete-output oracle, with no executed amplitudes as input."""
    target = record_field.second_quantization(record_field.permutation())
    permutation = list(range(record_field.MODES))
    permutation[2], permutation[3] = permutation[3], permutation[2]
    feedback = record_field.second_quantization(permutation)
    first_number = record_field.annihilator(2).T@record_field.annihilator(2)
    second_number = record_field.annihilator(1).T@record_field.annihilator(1)
    identity = np.eye(record_field.DIMENSION)
    first = (identity-first_number, first_number)
    second = (identity-second_number, second_number)
    initial = np.zeros((record_field.DIMENSION, 2), complex)
    initial[1, 0] = initial[2, 1] = 1/math.sqrt(2)
    expected = {'initial': initial}
    accepted = 15/16
    for bit in (0, 1):
        state = first[bit]@target@initial
        state = feedback@state if bit else state
        expected[f'after_first_{bit}'] = state
        for raw in (30, 31):
            expected[f'selector_timeout_second_{bit}_{raw}'] = state*math.sqrt(accepted/32)
        for last in (0, 1):
            terminal = second[last]@target@state
            expected[f'terminal_{bit}_{last}'] = terminal
            expected[f'weighted_terminal_{bit}_{last}'] = terminal*accepted
    for raw in (30, 31):
        expected[f'selector_timeout_first_{raw}'] = initial/math.sqrt(32)
    return expected


def _check_field_blocks(arrays):
    expected = _ideal_field_blocks()
    if arrays.keys() != expected.keys():
        raise ArithmeticError('field/reference output menu differs from complete oracle')
    for name, value in arrays.items():
        if (value.shape != expected[name].shape or not np.isfinite(value).all()
                or np.linalg.norm(value-expected[name], 2) > FIELD_OUTPUT_OPERATOR_GATE):
            raise ArithmeticError(f'field/reference output failed independent oracle: {name}')


def _field_blocks():
    """Check every source-executed quantum output before retaining it."""
    _, arrays = record_field.detector_experiment()
    _check_field_blocks(arrays)
    return {name: {'shape': list(value.shape), 'real': value.real.tolist(),
                   'imaginary': value.imag.tolist()} for name, value in sorted(arrays.items())}


def execution_diagnostics():
    """Fresh local observations, including backend-dependent numerical norms.

    These observations are not portable historical values authenticated by
    receipt replay. The canonical run retains their fixed admission intervals.
    """
    # Cheap exact/oracle checks precede the full native source execution.
    capacity = record_capacity.certificate()
    code = record_code.exact_certificate()
    outputs = _field_blocks()
    executed = record_code.source_execution_report()
    budget = record_code.source_composition_report()
    field = record_field.report(native_source_budget=budget)
    for law in ('C', 'E'):
        error = field['complete_source_composition']['laws'][law][
            'maximum_complete_two_cycle_operator_error_diagnostic']
        if not math.isfinite(error) or not 0 <= error <= 1e-6:
            raise ArithmeticError('complete field composition exceeds numerical guard')
    result = {
        'schema': 'oph.common_source.record_observer.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'optimal_native_code': code, 'rank15_recovery_obstruction': capacity,
        'native_source_execution': executed, 'all_seam_source_composition': budget,
        'protected_field': field, 'complete_field_reference_arrays': outputs,
        'information': record_information.report(),
        'scope': {'physical_theory_completed': False,
                  'native_record_law_selects_transport': False,
                  'coherent_code_preparation_storage_control_and_read_services_supplied': True,
                  'fixed_native_occupation_is_physical_electric_charge': False,
                  'all_public_measurements_reproduced': False},
    }
    # No NaN/Infinity can become a retained diagnostic or pass comparison.
    json.dumps(result, allow_nan=False)
    _compare_diagnostics(result, result)
    return result


def _numeric_interval(path):
    """The producer's named gates, applied before tolerant replay comparison."""
    if path.startswith('root.native_source_execution.'):
        for direction in ('forward', 'inverse'):
            prefix = f'root.native_source_execution.{direction}_full512_source_macro_guard.'
            if path.startswith(prefix):
                name = path[len(prefix):]
                if name in MACRO_GATES:
                    return 0., MACRO_GATES[name]
                if name == 'exact_polar_distance_formula_evaluated':
                    return 0., MACRO_GATES['projection_change']
        if (path == 'root.native_source_execution.encoder_operator_error'
                or path.startswith('root.native_source_execution.laws.') and
                path.rsplit('.', 1)[-1] in (
                    'reference_gram_error', 'maximum_branch_operator_error',
                    'complete_retained_output_operator_error', 'operator_error')):
            return 0., 1e-8
    if path.startswith('root.all_seam_source_composition.'):
        if '.mixing_angle_synthesis_diagnostics.' in path:
            return 0., 1e-10
        if path.endswith('_full512_macro_polar_change'):
            return 0., MACRO_GATES['projection_change']
        return 0., 1e-6
    if path.startswith('root.protected_field.complete_source_composition.'):
        if path in (
                'root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_complement_operator_norm',
                'root.protected_field.complete_source_composition.proper_code_export_diagnostic.maximum_success_norm_excess'):
            return 0., 2e-11
        if path.endswith('.maximum_complete_cq_diamond_error_diagnostic'):
            return 0., float(Fraction(1, 10**6)*(2+Fraction(1, 10**6)))
        if path.endswith('.maximum_complete_cq_half_trace_error_diagnostic'):
            return 0., float(Fraction(1, 10**6)*(1+Fraction(1, 2*10**6)))
        if path.endswith(('.walk_operator_error', '.feedback_operator_error')):
            # Both source words gate max entry error; ||M||_2 <= 64 max_ij |M_ij|.
            return 0., record_field.DIMENSION*2e-11
        return 0., 1e-6
    if (path.startswith('root.protected_field.dispersion.modes[')
            and path.endswith('.source_residual')):
        return 0., 2e-11
    if path in (
            'root.protected_field.source.source_walk_max_entry_residual',
            'root.protected_field.source.walk_resources.maximum_numerical_valid_code_leakage',
            'root.protected_field.source.feedback_resources.maximum_numerical_valid_code_leakage'):
        return 0., 2e-11
    if path == 'root.protected_field.dispersion.charged_retarded_anticommutator_response_residual':
        return 0., 3e-11
    clock = 'root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.'
    if path == clock+'complete_retained_stopped_operator_error':
        return 0., 5e-11
    if path == clock+'complete_stopped_gram_error':
        return 0., 1e-10
    if path == 'root.protected_field.two_read_experiment.weighted_terminal_source_probability_sum':
        return 1-4e-11, 1+4e-11
    if path.startswith('root.information.'):
        name = path.rsplit('.', 1)[-1]
        if name in ('reference_marginal_error', 'complete_reference_marginal_error'):
            return 0., record_information.CHECK_TOLERANCE
        if name == 'quantum_envelope_residual':
            return -record_information.CHECK_TOLERANCE, record_information.CHECK_TOLERANCE
        if name == 'information_loss':
            return -record_information.CHECK_TOLERANCE, math.inf
    return None


def _acceptance_gate(path):
    if path == 'root.native_source_execution.maximum_allowed_numerical_execution_error':
        return 1e-8
    for direction in ('forward', 'inverse'):
        prefix = f'root.native_source_execution.{direction}_full512_source_macro_guard.maximum_'
        if path.startswith(prefix):
            return MACRO_GATES.get(path[len(prefix):])
    return None


def _derived_equal(actual, expected, path):
    # Relative comparison on the residual's own scale; a nanoscopic diagnostic
    # cannot be replaced by zero under the receipt's ordinary absolute tolerance.
    if not math.isclose(actual, expected, rel_tol=2e-12, abs_tol=0.):
        raise ArithmeticError(f'inconsistent derived diagnostic at {path}')


def _check_relations(row, path):
    """Verify within-receipt formulas separately from cross-backend agreement."""
    if path == 'root.complete_field_reference_arrays':
        arrays = {name: np.array(value['real'])+1j*np.array(value['imaginary'])
                  for name, value in row.items()}
        if any(list(arrays[name].shape) != value['shape'] for name, value in row.items()):
            raise ArithmeticError('retained field array shape differs from its declaration')
        _check_field_blocks(arrays)
    elif path.endswith('_full512_source_macro_guard'):
        delta = row['raw_unitarity_defect']
        _derived_equal(row['exact_polar_distance_formula_evaluated'],
                       delta/(1+math.sqrt(1-delta)), path)
    elif path == 'root.all_seam_source_composition':
        for direction in ('forward', 'inverse'):
            _derived_equal(row[f'{direction}_full512_macro_operator_error'],
                           row[f'{direction}_full512_stabilized_macro_operator_error']
                           +row[f'{direction}_full512_macro_polar_change'], path)
        for action in ('encoder', 'decoder'):
            _derived_equal(row[f'{action}_operator_error_bound'],
                           math.expm1(row[f'{action}_log_error']), path)
        for local in row['laws'].values():
            log = local['maximum_corrected_attempt_log_error']
            _derived_equal(local['maximum_corrected_attempt_operator_error'], math.expm1(log), path)
            _derived_equal(local['maximum_encoded_corrected_attempt_operator_error'],
                           math.expm1(row['encoder_log_error']+log), path)
    elif path.startswith('root.all_seam_source_composition.mixing_angle_synthesis_diagnostics.'):
        _derived_equal(row['uniform_target_unitary_completion_error_upper'],
                       2*row['requested_word_uniform_column_error'], path)
    elif path.startswith('root.protected_field.complete_source_composition.laws.') and path.count('.') == 4:
        epsilon = row['maximum_complete_two_cycle_operator_error_diagnostic']
        _derived_equal(epsilon, row['maximum_successful_block_two_cycle_operator_error_diagnostic']
                       +row['maximum_proper_code_failure_amplitude_diagnostic'], path)
        _derived_equal(row['maximum_complete_cq_diamond_error_diagnostic'], epsilon*(2+epsilon), path)
        _derived_equal(row['maximum_complete_cq_half_trace_error_diagnostic'], epsilon*(1+epsilon/2), path)
    elif path.startswith('root.protected_field.two_read_experiment.complete_detector_outputs['):
        expected = .5 if row['first_read'] == row['second_read'] else 0.
        actual = row['source_probability_before_selector_weights']
        if actual < 0 or abs(actual-expected) > 3e-11:
            raise ArithmeticError(f'detector probability exceeds producer gate at {path}')
    elif path == 'root':
        native = row['all_seam_source_composition']
        field = row['protected_field']['complete_source_composition']
        for direction in ('forward', 'inverse'):
            _derived_equal(native[f'{direction}_full512_macro_polar_change'],
                           row['native_source_execution'][f'{direction}_full512_source_macro_guard']['projection_change'], path)
        for law in ('C', 'E'):
            log = (6*native['encoder_log_error']
                   +12*(native['laws'][law]['maximum_corrected_attempt_log_error']
                        +native['decoder_log_error']+native['encoder_log_error'])
                   +2*math.log1p(field['walk_operator_error'])+math.log1p(field['feedback_operator_error']))
            local = field['laws'][law]
            _derived_equal(local['maximum_successful_block_two_cycle_operator_error_diagnostic'],
                           math.expm1(log), path)
            export = field['proper_code_export_diagnostic']
            attempts = local['maximum_two_cycle_resources']['native_pulses']
            failure = (math.sqrt(attempts)*export['maximum_complement_operator_norm']
                       *math.exp((attempts-1)*math.log1p(export['maximum_success_norm_excess'])))
            _derived_equal(local['maximum_proper_code_failure_amplitude_diagnostic'], failure, path)


def _compare(retained, replayed, path='root', *, _diagnostics=False):
    """Compare canonical evidence; historical backend observations are absent."""
    if type(retained) is not type(replayed):
        raise ArithmeticError(f'receipt type mismatch at {path}')
    if type(replayed) is dict:
        if retained.keys() != replayed.keys():
            raise ArithmeticError(f'receipt fields mismatch at {path}')
        for key in replayed:
            _compare(retained[key], replayed[key], f'{path}.{key}', _diagnostics=_diagnostics)
        if (_diagnostics or path == 'root.complete_field_reference_arrays'
                or path.startswith('root.protected_field.two_read_experiment.complete_detector_outputs[')):
            _check_relations(retained, path)
            _check_relations(replayed, path)
    elif type(replayed) is list:
        if len(retained) != len(replayed):
            raise ArithmeticError(f'receipt length mismatch at {path}')
        for index, (left, right) in enumerate(zip(retained, replayed, strict=True)):
            _compare(left, right, f'{path}[{index}]', _diagnostics=_diagnostics)
    elif type(replayed) is float:
        if not _diagnostics and path in OBSERVATION_PATHS:
            raise ArithmeticError(f'backend observation is not canonical evidence at {path}')
        gate = _acceptance_gate(path)
        if gate is not None:
            if retained != gate or replayed != gate:
                raise ArithmeticError(f'acceptance gate changed at {path}')
            return
        interval = _numeric_interval(path)
        if interval is not None and any(not math.isfinite(value) or not interval[0] <= value <= interval[1]
                                        for value in (retained, replayed)):
            raise ArithmeticError(f'source diagnostic exceeds gate at {path}')
        # This branch is used only to validate fresh local observations before
        # omission. It never authenticates historical numerical measurements.
        if _diagnostics and (path.startswith('root.all_seam_source_composition.')
                or path.startswith('root.protected_field.complete_source_composition.')):
            return
        elif (not math.isfinite(retained) or not math.isfinite(replayed)
              or not math.isclose(retained, replayed, rel_tol=NUMERICAL_COMPARISON_TOLERANCE,
                                  abs_tol=NUMERICAL_COMPARISON_TOLERANCE)):
            raise ArithmeticError(f'numerical replay mismatch at {path}')
    elif retained != replayed:
        raise ArithmeticError(f'exact replay mismatch at {path}')


def _compare_diagnostics(left, right, path='root'):
    """Internal gate/formula checks, not a historical-observation verifier."""
    _compare(left, right, path, _diagnostics=True)


def _observation_paths():
    """Explicit inventory: new or missing numerical diagnostics fail closed."""
    paths = {'root.native_source_execution.encoder_operator_error'}
    for direction in ('forward', 'inverse'):
        paths.update(f'root.native_source_execution.{direction}_full512_source_macro_guard.{name}'
                     for name in (*MACRO_GATES, 'exact_polar_distance_formula_evaluated'))
    for law, branches in (('C', 8), ('E', 16)):
        paths.update(f'root.native_source_execution.laws.{law}.{name}' for name in (
            'reference_gram_error', 'maximum_branch_operator_error', 'complete_retained_output_operator_error'))
        paths.update(f'root.native_source_execution.laws.{law}.branches[{i}].operator_error'
                     for i in range(branches))
        paths.update(f'root.all_seam_source_composition.laws.{law}.{name}' for name in (
            'maximum_corrected_attempt_log_error', 'maximum_corrected_attempt_operator_error',
            'maximum_encoded_corrected_attempt_operator_error'))
        paths.update(f'root.all_seam_source_composition.laws.{law}.complete_local_native_dilation_operator_errors_by_coin[{i}]'
                     for i in (0, 1))
        paths.update(f'root.protected_field.complete_source_composition.laws.{law}.{name}' for name in (
            'maximum_successful_block_two_cycle_operator_error_diagnostic',
            'maximum_proper_code_failure_amplitude_diagnostic',
            'maximum_complete_two_cycle_operator_error_diagnostic',
            'maximum_complete_cq_diamond_error_diagnostic',
            'maximum_complete_cq_half_trace_error_diagnostic'))
    paths.update('root.all_seam_source_composition.'+name for name in (
        'encoder_log_error', 'encoder_operator_error_bound', 'decoder_log_error',
        'decoder_operator_error_bound', 'fair_coin_preparation_operator_error',
        'source_X_conjugation_operator_error',
        *(f'{direction}_full512_{name}' for direction in ('forward', 'inverse') for name in (
            'macro_operator_error', 'macro_polar_change', 'stabilized_macro_operator_error'))))
    paths.update(f'root.all_seam_source_composition.mixing_angle_synthesis_diagnostics.{count}.{name}'
                 for count in (1, 4, 5, 6) for name in (
                     'requested_inverse_adjoint_error', 'requested_inverse_uniform_to_zero_error',
                     'requested_word_uniform_column_error', 'requested_word_unitarity_residual',
                     'uniform_target_unitary_completion_error_upper'))
    paths.update('root.protected_field.complete_source_composition.'+name for name in (
        'walk_operator_error', 'feedback_operator_error',
        'proper_code_export_diagnostic.maximum_complement_operator_norm',
        'proper_code_export_diagnostic.maximum_success_norm_excess'))
    paths.update(f'root.protected_field.dispersion.modes[{i}].source_residual' for i in range(6))
    paths.update(('root.protected_field.dispersion.charged_retarded_anticommutator_response_residual',
        'root.protected_field.source.source_walk_max_entry_residual',
        'root.protected_field.source.walk_resources.maximum_numerical_valid_code_leakage',
        'root.protected_field.source.feedback_resources.maximum_numerical_valid_code_leakage',
        'root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.complete_retained_stopped_operator_error',
        'root.protected_field.clocks.joint_clock_histories.source_executed_stopped_channel.complete_stopped_gram_error'))
    paths.update(f'root.information.same_field_experiment_reference.{law}.complete_reference_marginal_error'
                 for law in ('C', 'E'))
    paths.update(f'root.information.same_native_population_reference.{sample}.{name}'
                 for sample in ('complete_C_with_timeouts', 'complete_E',
                                'destroyed_core_coherence', 'surviving_spectator_coherence')
                 for name in ('reference_marginal_error', 'quantum_envelope_residual'))
    return frozenset(paths)


OBSERVATION_PATHS = _observation_paths()
DIAGNOSTIC_SUBTREES = ('root.native_source_execution.', 'root.all_seam_source_composition.',
                      'root.protected_field.complete_source_composition.')


def _canonical_receipt(diagnostics):
    # Validate before omitting anything. This private transformer does not
    # establish provenance; run()/verify() obtain their own fresh execution.
    _compare_diagnostics(diagnostics, diagnostics)
    omitted, skip = set(), object()

    def retain(value, path):
        if path in OBSERVATION_PATHS:
            if type(value) is not float:
                raise ArithmeticError(f'original floating diagnostic required at {path}')
            omitted.add(path)
            return skip
        if type(value) is dict:
            result = {key: retain(item, f'{path}.{key}') for key, item in value.items()}
            return {key: item for key, item in result.items() if item is not skip}
        if type(value) is list:
            result = [retain(item, f'{path}[{i}]') for i, item in enumerate(value)]
            if result and all(item is skip for item in result):
                return skip
            if any(item is skip for item in result):
                raise ArithmeticError(f'partially registered observation list at {path}')
            return result
        if type(value) is float and path.startswith(DIAGNOSTIC_SUBTREES) and _acceptance_gate(path) is None:
            raise ArithmeticError(f'unregistered numerical diagnostic at {path}')
        return value

    result = retain(diagnostics, 'root')
    if omitted != OBSERVATION_PATHS:
        raise ArithmeticError('registered numerical diagnostics are missing from fresh execution')
    intervals = {}
    for path in sorted(OBSERVATION_PATHS):
        interval = _numeric_interval(path)
        if interval is None or any(not math.isfinite(value) for value in interval):
            raise ArithmeticError(f'a fixed finite diagnostic admission interval required at {path}')
        intervals[path] = [str(Fraction(str(value))) for value in interval]
    result['schema'] = 'oph.common_source.record_observer.v2'
    result['numerical_admission'] = {
        'intervals_exact': intervals,
        'complete_two_cycle_operator_ceiling_exact': '1/1000000',
        'complete_cq_diamond_ceiling_exact': '2000001/1000000000000',
        'complete_cq_half_trace_ceiling_exact': '2000001/2000000000000',
        'every_interval_checked_on_owned_fresh_execution': True,
        'backend_observations_retained': False,
        'outward_rounding_or_physical_error_certificate': False,
        'fresh_observation_interface': 'execution_diagnostics() or --diagnostics',
    }
    return result


def run():
    """Portable evidence with fixed admission ceilings and owned fresh replay."""
    return _canonical_receipt(execution_diagnostics())


def verify(retained):
    """Own the entire replay; a caller cannot supply a matching authority."""
    _compare(retained, run())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument('--write', type=Path)
    files.add_argument('--check', type=Path)
    files.add_argument('--diagnostics', action='store_true', help='print fresh local numerical observations')
    args = parser.parse_args()
    if args.diagnostics:
        print(json.dumps(execution_diagnostics(), sort_keys=True, allow_nan=False))
        return
    if args.check:
        verify(strict_json.loads(args.check.read_text(encoding='utf-8')))
    else:
        result = run()
        if args.write:
            with args.write.open('x', encoding='utf-8', newline='\n') as stream:
                # Machine evidence is compact; the proof remains in ENCODED_MEMORY.
                stream.write(json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n')
    print(json.dumps({'optimal_native_code': 6, 'protected_CAR_modes': 6,
                      'numerical_operator_admission_ceiling_exact': '1/1000000',
                      'complete_field_reference_outputs_verified': True}))


if __name__ == '__main__':
    main()
