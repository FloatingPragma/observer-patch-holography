"""Recompute source words and every retained exact quantum output."""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as F
import hashlib
from pathlib import Path

from ..observer.scalar_response.experiment import compare_receipt
from ..quantum_observer_correspondence import construction as ideal
from . import history, source
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
ERROR = F(1, 10**11)
COMPONENT_GATE = F(1, 10**12)


def bound(*, gamma=1, reference=False, cap=2):
    history._arguments(gamma, cap, reference)
    prep = source.certify('preparation', reference=reference)
    gate = source.certify('apparatus', gamma=gamma)
    read = source.certify('readout')
    for item in (prep, gate, read):
        if item['operator_error_upper'] > COMPONENT_GATE:
            raise ArithmeticError('emitted source array exceeds its operation budget')
    success = 1-F(3, 5)**cap
    gain = 1+success*(gate['channel_norm_upper']-1)
    result = read['channel_norm_upper']*(gain*prep['half_diamond_upper']
              +success*gate['half_diamond_upper'])+read['half_diamond_upper']
    if result > ERROR:
        raise ArithmeticError('complete source-history bound exceeds its budget')
    return result


def costs(*, reference=False, cap=2):
    history._arguments(0, cap, reference)
    prep = source.resources('preparation', reference=reference)
    gate = source.resources('apparatus')
    read = source.resources('readout')
    success = 1-F(3, 5)**cap
    maximum = prep['native_pulses']+gate['native_pulses']+read['native_pulses']
    return {'preparation_pulses': prep['native_pulses'],
            'conditional_intervention_pulses': gate['native_pulses'],
            'terminal_pointer_pulses': read['native_pulses'],
            'source_pulses_maximum': maximum,
            'proper_code_loads_maximum': maximum,
            'proper_code_exports_maximum': maximum,
            'scheduled_native_attempts_maximum': 2*cap,
            'feedback_native_attempts_maximum': 1,
            'scheduled_native_attempts_ideal_expectation': str(success*(1+success)/F(2, 5)),
            'feedback_native_attempts_ideal_expectation': str(success),
            'source_pulses_ideal_expectation': str(prep['native_pulses']+success*gate['native_pulses']+read['native_pulses']),
            'work_code_qubits': 9, 'reference_qubits': int(reference),
            'fresh_terminal_pointer_qubits': 1,
            'peak_code_qubits_including_pointer_and_reference': 10+int(reference),
            'additional_source_processor_dimension': 6,
            'blank_work_registers': 1, 'blank_reference_qubits': int(reference),
            'blank_terminal_pointer_qubits': 1,
            'retained_clock_current_marks_maximum': 3,
            'fixed_record_buffer_bits': 11+(2*cap+3).bit_length(),
            'reusable_attempt_cap_counter_bits': cap.bit_length(),
            'terminal_pointer_bits': 1, 'attempt_cap_service_required': True,
            'complement_checks_before_native_attempts_maximum': 2*cap+1,
            'padding_abort_is_native_tick': False,
            'ideal_expectations_are_not_normalized_raw_masses': True}


def joint(outputs):
    result = defaultdict(F)
    for (word, status, bit), state in outputs.items():
        a = sum(mark[0][0] for _, mark in word)
        b = sum(mark[0][1] for _, mark in word)
        result[status, a, b, bit] += history.trace(state)
    return dict(result)


def audit_output(actual, target, *, cap=2, reference=False):
    actual = history.validate_history(actual, cap=cap, reference=reference)
    target = history.validate_history(target, cap=cap, reference=reference, native_target=True)
    error = history.distance_upper(actual, target)
    mass = sum(map(history.trace, actual.values()), F(0))
    if error > ERROR or abs(mass-1) > 2*ERROR:
        raise ArithmeticError('exact complete emitted history misses the independent target')
    abort = sum((history.trace(state) for (_, status, _), state in actual.items()
                 if status.startswith('code_abort')), F(0))
    if not 0 <= abort <= 2*ERROR:
        raise ArithmeticError('code-complement mass is outside its complete allowance')
    return {'complete_output_records': len(actual),
            'inherited_nonabort_records': len(target),
            'explicit_code_complement_records': len(actual)-len(target),
            'exact_entrywise_half_nuclear_error_upper': str(ERROR),
            'exact_mass_defect_upper': str(2*ERROR),
            'code_complement_mass_upper': str(2*ERROR),
            'all_quantum_reference_entries_compared': True,
            'no_amplitude_projection_or_normalization': True}


def run():
    cases, joints = [], {}
    for reference in (False, True):
        for gamma in (0, 1):
            bound(gamma=gamma, reference=reference)
            checkpoint = None
            if reference and gamma == 1:
                checkpoint = history.checkpoint_json(history.first_stage(gamma=gamma, reference=True))
            actual = history.execute(gamma=gamma, reference=reference, checkpoint=checkpoint)
            target = history.ideal_history(gamma=gamma, reference=reference)
            row = audit_output(actual, target, reference=reference)
            row.update({'gamma': gamma, 'reference_dimension': 1+int(reference),
                        'reference_uniform_source_composition_error_upper': str(ERROR),
                        'ideal_timeout_probability': '369/625'})
            if checkpoint is not None:
                row['serialized_reference_prefix_exact_and_full_continuation_executed'] = True
            cases.append(row)
            if not reference:
                joints[gamma] = (joint(actual), joint(target))
    first, second = joints[0][1], joints[1][1]
    ideal_response = sum((abs(first.get(key, F(0))-second.get(key, F(0)))
                         for key in first.keys() | second.keys()), F(0))/2
    raw_first, raw_second = joints[0][0], joints[1][0]
    actual_response = sum((abs(raw_first.get(key, F(0))-raw_second.get(key, F(0)))
                          for key in raw_first.keys() | raw_second.keys()), F(0))/2
    if ideal_response <= 2*ERROR or abs(actual_response-ideal_response) > 2*ERROR:
        raise ArithmeticError('compiled clock/field intervention response was lost')
    paths = set(source.bindings()) | set(ideal.SOURCE_PATHS) | {
        'model/compiled_clock_history/history.py', 'model/compiled_clock_history/receipt.py',
        'model/compiled_clock_history/__main__.py', 'model/requirements.txt',
        'model/observer/scalar_response/experiment.py', 'tools/strict_json.py'}
    return {'schema': 'oph.compiled_clock_history.v1', 'all_pass': True,
            'source_sha256': bound_pins({path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                              for path in sorted(paths)}),
            'contract': {'native_law': 'C', 'binary_total': 2, 'physical_ports': 12,
                         'admitted_dimension': 66, 'full_encoded_dimension': 512,
                         'union_ports': list(source.UNION), 'cap_per_interval': 2,
                         'native_incident_probability': '2/5',
                         'raw_binary_arrays_propagated_exactly': True,
                         'code_complement_abort_retained': True,
                         'supplied_ideal_coherent_gates_in_executor': False},
            'executions': cases,
            'intervention_response': {'observable': 'complete (status, clock A, clock B, detector bit) law',
                'ideal_half_l1_response': str(ideal_response),
                'compiled_half_l1_response_lower': str(ideal_response-2*ERROR),
                'compiled_half_l1_response_upper': str(ideal_response+2*ERROR),
                'strictly_positive': True, 'native_clock_marginal_is_not_the_complete_output': True},
            'resources': [costs(reference=reference) for reference in (False, True)],
            'scope': {'source_controls_emitted_and_executed': True,
                      'all_later_density_arithmetic_exact': True,
                      'compiled_preparation_from_declared_blank': True,
                      'source_native_repair_instrument_is_primitive': True,
                      'encoding_and_code_complement_check_are_declared_services': True,
                      'physical_service_error_certified': False,
                      'physical_seconds_or_relativistic_field_recovery': False,
                      'original_eight_mode_benchmark_completed': False,
                      'raw_array_maps_exactly_trace_preserving_claimed': False,
                      'effective_correspondence_join': 'complete bound to the existing exact decoded ideal history'}}
