"""Data-current feedback driven by the same retained native clock as cooling.

All eighteen input currents pass through x <- (I-G/2)x. G is derived from
the biased native reference, not a target current menu. One real C record
triggers both bath replacement and one update. The finite run retains every
clock record, bath, timeout and native data-repair output. The limiting
twelve-current comparison is an output theorem, not a finite exact cutoff.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy.linalg import expm

from . import native_response_biased as biased, native_response_cooling as cooling
from . import native_response_dynamics as dynamics, native_response_reference as reference
from . import native_code_transport as transport
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
HORIZON = 512
CAP = 16
BIAS = F(16)
STEP = F(1, 2)
PARAMETER = F(1, 8)
SOURCE_ADMISSION = 2e-7
COEFFICIENT_ADMISSION = 2e-13


def _count(value, name, maximum=4096):
    if type(value) is not int or not 0 <= value <= maximum:
        raise ValueError(name+' must be an original bounded nonnegative integer')
    return value


def terminal_count(records, cap=CAP):
    """Literal controller state transition on retained raw native C records."""
    cap = _count(cap, 'cap')
    if type(records) not in (list, tuple) or len(records) > 4096:
        raise ValueError('a bounded complete sequence of raw clock records is required')
    labels = {(coin, raw) for coin, raw, _ in cooling.raw_records()}
    count = 0
    for record in records:
        if (type(record) is not tuple or len(record) != 2
                or any(type(x) is not int for x in record) or record not in labels):
            raise ValueError('an original native clock record pair is required')
        count = min(cap, count+int(record == cooling.TRIGGER))
    return count


def distribution(horizon=HORIZON, cap=CAP, *, weights=None):
    """All raw histories aggregated by count; no sampling or normalization."""
    horizon, cap = _count(horizon, 'horizon'), _count(cap, 'cap')
    records = cooling.raw_records()
    if weights is None:
        weights = [row[2] for row in records]
        one, zero = F(1), F(0)
    else:
        if (type(weights) is not np.ndarray or weights.shape != (8,)
                or weights.dtype != np.dtype(float) or not np.isfinite(weights).all()
                or np.any(weights < 0)):
            raise ValueError('eight finite nonnegative raw source probabilities required')
        one, zero = 1., 0.
    trigger = sum(w for row, w in zip(records, weights, strict=True) if row[:2] == cooling.TRIGGER)
    rest = sum(w for row, w in zip(records, weights, strict=True) if row[:2] != cooling.TRIGGER)
    result = [one]+[zero]*cap
    for _ in range(horizon):
        update = [zero]*(cap+1)
        for k, value in enumerate(result):
            update[k] += rest*value
            update[min(k+1, cap)] += trigger*value
        result = update
    return result


@lru_cache(maxsize=1)
def _source_data():
    p = biased.exact_data(BIAS)['first_triplet_weight']
    gram = np.asarray(dynamics.cost_gram(p).evalf(30), float)
    basis = np.asarray([np.asarray(x.evalf(30), complex) for x in dynamics.full_basis()])
    return p, gram, basis


def update(coefficients):
    """The full source-Gram recurrence, including all eighteen coefficients."""
    if (type(coefficients) is not np.ndarray or coefficients.shape not in ((18,), (18, 18))
            or coefficients.dtype != np.dtype(float) or not np.isfinite(coefficients).all()):
        raise ValueError('original finite eighteen-coefficient data required')
    return _update_matrix()@coefficients


@lru_cache(maxsize=1)
def _update_matrix():
    return np.eye(18)-float(STEP)*_source_data()[1]


def _coefficients():
    values = [np.eye(18)]
    for _ in range(CAP):
        values.append(update(values[-1]))
    return np.asarray(values)


def _generator(coefficients, mu):
    return np.einsum('j,jab->ab', coefficients[:, mu], _source_data()[2])


def _target(generator):
    target = np.eye(8, dtype=complex)
    target[:6, :6] = expm(float(PARAMETER)*generator)
    return target


def _word(generator):
    # Every full-frame primitive after feedback remains a real rotation or
    # an imaginary real-symmetric phase, including the coupled block centers.
    if np.max(abs(generator.imag)) < 1e-15:
        return compiler.compile_real_orthogonal(_target(generator).real)
    if np.max(abs(generator.real)) < 1e-15:
        return reference._phase_word(generator.imag, PARAMETER)
    raise ArithmeticError('the complete-frame primitive has no registered source compilation')


def _schema(arrays):
    shapes = {'coefficients': ((CAP+1, 18, 18), np.dtype(float)),
              'controls': ((CAP+1, 18, 8, 8), np.dtype(complex))}
    if type(arrays) is not dict or arrays.keys() != shapes.keys():
        raise ValueError('complete controller array menu required')
    for name, (shape, dtype) in shapes.items():
        x = arrays[name]
        if type(x) is not np.ndarray or x.shape != shape or x.dtype != dtype or not np.isfinite(x).all():
            raise ValueError('controller source array schema differs: '+name)


def work_output(arrays, reference_arrays, count, primitive, records, syndromes, incoming):
    """Actual two-carrier C/E output block, retaining all64 low rows.

    The two syndrome indices label coherent output tensor legs. The incoming
    eight-row columns may carry any untouched reference. The theorem compares
    the first six input rows; the complete padding columns remain available.
    """
    count = _count(count, 'count', CAP)
    primitive = _count(primitive, 'primitive', 17)
    if type(records) is not tuple or len(records) != 2 or type(syndromes) is not tuple or len(syndromes) != 2:
        raise ValueError('two native records and two coherent syndrome indices required')
    if (type(incoming) is not np.ndarray or incoming.ndim not in (1, 2) or incoming.shape[0] != 8
            or incoming.dtype.kind not in 'iufc' or not np.isfinite(incoming).all()):
        raise ValueError('finite eight-row input/reference columns required')
    blocks = []
    for law, record, high in zip(('C', 'E'), records, syndromes, strict=True):
        if type(record) is not tuple or len(record) != 2:
            raise ValueError('native coin/raw pairs required')
        coin = _count(record[0], 'coin', 1)
        raw = _count(record[1], 'raw', 3 if law == 'C' else 7)
        high = _count(high, 'coherent syndrome', 63)
        row = coin*(4 if law == 'C' else 8)+raw
        blocks.append(reference_arrays['decoded_'+law][row, 8*high:8*(high+1)])
    state = arrays['controls'][count, primitive]@np.asarray(incoming, complex).reshape(8, -1)
    state = (reference_arrays['geometric_preparation'][:, 8*np.arange(8)]@state).reshape(8, 8, -1)
    state = np.einsum('ai,bj,ijr->abr', *blocks, state, optimize=True).reshape(64, -1)
    state = reference_arrays['geometric_inverse']@state
    return state[:, 0] if incoming.ndim == 1 else state


def _binary_scale(value):
    """Exact power-of-two scale with an explicit lost-component refusal."""
    if not np.isfinite(value).all():
        raise ArithmeticError('history contraction exceeds finite arithmetic')
    peak = max(float(np.max(abs(value.real), initial=0.)),
               float(np.max(abs(value.imag), initial=0.)))
    if peak == 0:
        return value, 0
    exponent = math.frexp(peak)[1]
    scaled = np.ldexp(value.real, -exponent)+1j*np.ldexp(value.imag, -exponent)
    if np.any((value.real != 0) & (scaled.real == 0)) or np.any((value.imag != 0) & (scaled.imag == 0)):
        raise ArithmeticError('history components exceed the supported relative dynamic range')
    return scaled, exponent


def history_output_block(arrays, reference_arrays, cooling_arrays, *, clock_records,
                         counter_rows, bath_rows, primitive, work_records,
                         work_syndromes, incoming, scaled=False):
    """Sequential complete-history tensor slice, with the actual clock count.

    Incoming axes are reference64, work8, untouched-reference columns. Each
    bath/counter index is a coherent output leg, not an observed or accepted
    event. Varying the indices reconstructs all outputs. The same reference
    propagates through successive baths; it is never reinitialized in place.
    Short prefixes support direct independent contraction tests. The complete
    run uses512 stages, represented by these reusable factors without an
    exponentially large output array. With scaled=True the result is
    (mantissa, exponent) and the physical amplitude is mantissa*2**exponent.
    The exponent is retained, not a conditioning or probability normalization.
    The ordinary array export refuses overflow or a vanished nonzero component;
    use the scaled result for exponentially small complete-history slices.
    """
    if type(scaled) is not bool:
        raise ValueError('scaled must be an original Boolean')
    count = terminal_count(clock_records)
    length = len(clock_records)
    if length > HORIZON or type(counter_rows) is not tuple or type(bath_rows) is not tuple or len(counter_rows) != length or len(bath_rows) != length:
        raise ValueError('one coherent counter and bath output index per clock stage required')
    if (type(incoming) is not np.ndarray or incoming.ndim != 3 or incoming.shape[:2] != (64, 8)
            or incoming.dtype.kind not in 'fc' or not np.isfinite(incoming).all()):
        raise ValueError('finite reference64/work8/ancilla-column input required')
    state, exponent = _binary_scale(incoming.astype(complex, copy=True))
    for record, counter_row, bath_row in zip(clock_records, counter_rows, bath_rows, strict=True):
        counter_row = _count(counter_row, 'coherent counter output', 511)
        bath_row = _count(bath_row, 'coherent bath output', 63)
        factor = cooling_arrays['replaced'] if record == cooling.TRIGGER else cooling_arrays['unchanged']
        block = factor.reshape(64, 64, 64)[bath_row]
        amplitude, amplitude_exponent = _binary_scale(np.asarray(
            cooling_arrays['counter'][4*record[0]+record[1], counter_row]))
        state, step_exponent = _binary_scale(amplitude*np.einsum('ab,bir->air', block, state, optimize=True))
        exponent += amplitude_exponent+step_exponent
    columns = state.transpose(1, 0, 2).reshape(8, -1)
    output = work_output(arrays, reference_arrays, count, primitive, work_records, work_syndromes, columns)
    output, final_exponent = _binary_scale(output.reshape(64, 64, incoming.shape[2]).transpose(1, 0, 2))
    exponent += final_exponent
    if scaled:
        return output, exponent
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        ordinary = np.ldexp(output.real, exponent)+1j*np.ldexp(output.imag, exponent)
    if (not np.isfinite(ordinary).all() or np.any((output.real != 0) & (ordinary.real == 0))
            or np.any((output.imag != 0) & (ordinary.imag == 0))):
        raise ArithmeticError('history amplitude cannot be exported in binary64; request scaled=True')
    return ordinary


def metrics(arrays, reference_arrays, cooling_arrays):
    _schema(arrays)
    # Every parent factor used by history_output_block is part of this
    # admission, including unselected bath columns and native raw branches.
    reference.metrics(reference_arrays)
    cooling.metrics(cooling_arrays)
    p = _source_data()[0]
    # Independent exact spectral power checks the literal numerical recurrence.
    exact = np.asarray([np.asarray(dynamics.feedback_matrix(p, k, STEP).evalf(30), float)
                        for k in range(CAP+1)])
    coefficient_errors = np.linalg.norm(arrays['coefficients']-exact, axis=(1, 2))
    if max(coefficient_errors) > COEFFICIENT_ADMISSION:
        raise ValueError('source-Gram feedback differs from independent exact matrix powers')
    errors = np.zeros((CAP+1, 18))
    for k in range(CAP+1):
        for mu in range(18):
            errors[k, mu] = np.linalg.norm(arrays['controls'][k, mu]-_target(_generator(exact[k], mu)), 2)
    if np.max(errors) > SOURCE_ADMISSION:
        raise ValueError('complete data-current source words do not implement the feedback')
    weights = np.sum(abs(cooling_arrays['counter'])**2, axis=1).astype(float)
    observed = np.asarray(distribution(weights=weights))
    ideal = np.asarray([float(x) for x in distribution()])
    if max(abs(observed-ideal)) > HORIZON*6*SOURCE_ADMISSION:
        raise ValueError('retained native clock probabilities exceed whole-run admission')
    records, weight = [], 1.
    for law in ('C', 'E'):
        code = transport.PAIR_CODES[law]
        record = next((c, r) for c in (0, 1) for r in range(4 if law == 'C' else 8)
                      if code.branch_weight(law, transport.SEAMS[law], c, r))
        records.append(record)
        weight *= float(code.branch_weight(law, transport.SEAMS[law], *record))
    amplitude_errors = []
    for k in (0, 1, CAP):
        for mu in range(18):
            actual = work_output(arrays, reference_arrays, k, mu, tuple(records), (0, 0), np.eye(8, 6, dtype=complex))
            target = np.zeros((64, 6), complex)
            target[8*np.arange(8)] = math.sqrt(weight)*_target(_generator(exact[k], mu))[:, :6]
            amplitude_errors.append(float(np.linalg.norm(actual-target, 2)))
    if max(amplitude_errors) > 6*SOURCE_ADMISSION:
        raise ValueError('actual controller/native work output amplitudes disagree')
    # Exercise the full clock horizon through the history API. Independent
    # high-precision scalar powers check the retained scale without forming
    # exponentially small binary64 amplitudes. This supplements, rather than
    # replaces, the complete-factor operator bounds above.
    ground = cooling_arrays['bath_preparation'][:, 0]
    bath_row = int(np.argmax(abs(ground)))
    drive = np.zeros(8, complex); drive[[3, 5]] = 1/math.sqrt(2)
    initial = np.outer(ground, drive)[:, :, None]
    long_errors, long_exponents = [], []
    for triggers in (0, CAP):
        clock_records = (cooling.TRIGGER,)*triggers+((0, 0),)*(HORIZON-triggers)
        counter_rows = tuple(int(np.argmax(abs(cooling_arrays['counter'][4*c+r]))) for c, r in clock_records)
        block, exponent = history_output_block(arrays, reference_arrays, cooling_arrays,
            clock_records=clock_records, counter_rows=counter_rows, bath_rows=(bath_row,)*HORIZON,
            primitive=12, work_records=tuple(records), work_syndromes=(0, 0), incoming=initial, scaled=True)
        data = work_output(arrays, reference_arrays, triggers, 12, tuple(records), (0, 0), drive)
        base = np.outer(ground, data)
        if np.linalg.norm(base) < .01:
            raise ValueError('the long-history positive control has no resolved output')
        with mp.workdps(90):
            expected = mp.mpc(1)
            for c, r in ((0, 0), cooling.TRIGGER):
                row = int(np.argmax(abs(cooling_arrays['counter'][4*c+r])))
                scalar = mp.mpc(complex(cooling_arrays['counter'][4*c+r, row]))*mp.mpc(complex(ground[bath_row]))
                expected *= scalar**(triggers if (c, r) == cooling.TRIGGER else HORIZON-triggers)
            # Compare the complete retained block at its own scale, including
            # relative phases. A single unaffected matrix entry could hide a
            # disconnected clock controlling a different work unitary.
            expected_block = complex(expected*mp.power(2, -exponent))*base
            error = float(np.linalg.norm(block[:, :, 0]-expected_block)/np.linalg.norm(expected_block))
        if not math.isfinite(error) or error > SOURCE_ADMISSION:
            raise ValueError('complete scaled history loses its amplitude, scale or shared clock')
        long_errors.append(error); long_exponents.append(exponent)
    return {'coefficient_operator_errors': coefficient_errors.tolist(),
        'complete_control_operator_errors': errors.tolist(),
        'native_clock_terminal_probabilities': observed.tolist(),
        'native_clock_total_probability': float(sum(observed)),
        'source_timeout_probability': float(sum(observed[:-1])),
        'selected_complete_work_amplitude_errors': amplitude_errors,
        'long_history_scaled_probe_relative_errors': long_errors,
        'long_history_scaled_probe_binary_exponents': long_exponents,
        'long_history_scaled_probe_terminal_counts': [0, CAP],
        'retained_primitive_controls': 18, 'retained_terminal_counts': CAP+1,
        'all_clock_records_and_timeouts_retained': True,
        'work_controls_act_on_data_not_only_reference': True}


def exact_certificate():
    probabilities = distribution()
    timeout = sum(probabilities[:-1])
    independent = sum(F(math.comb(HORIZON, k)*11**(HORIZON-k), 12**HORIZON) for k in range(CAP))
    if timeout != independent or sum(probabilities) != 1:
        raise ArithmeticError('complete raw-clock recurrence and exact binomial tail disagree')
    return {'horizon': HORIZON, 'update_cap': CAP, 'bias': str(BIAS), 'step': str(STEP),
        'data_drive_parameter': str(PARAMETER), 'trigger': list(cooling.TRIGGER),
        'terminal_probabilities_exact': [str(x) for x in probabilities],
        'timeout_probability_exact': str(timeout),
        'update_rule': 'x <- x-(1/2)G_p x on each actual C trigger, capped at16',
        'control_primitive_menu': 'the same complete HS-orthonormal u(3)+u(3) frame for every terminal count',
        'target_twelve_menu_or_rank_used_in_update': False,
        'phase_augmentation': 'the native fivefold class sum is available as a work actuator, not just a reference bias',
        'all_raw_histories_including_timeout_execute_their_actual_filtered_control': True,
        'same_trigger_resets_reference_and_updates_controller': True,
        'completed_clock_implies_reference_reset_in_ideal_protocol': True,
        'whole_output': 'history_output_block sequentially contracts every retained cooling stage on the SAME reference leg, derives terminal_count from those same raw records, and applies work_output to the data; arbitrary shared quantum-reference columns remain',
        'long_history_amplitude_export': 'scaled=True retains a mantissa and integer binary exponent; unscaled export refuses vanished nonzero components rather than reporting a zero branch',
        'ideal_complete_history_bound': '2 epsilon |t| +2 Pr(timeout), with identical clock/bath outputs on both comparison sides',
        'comparison_domain': 'one unit-HS full-frame primitive at t=1/8, arbitrary work quantum reference; products use the sum of uniform per-stage drive bounds',
        'finite_exact_current_dimension_is_twelve': False,
        'finite_feedback_remains_full_rank': True,
        'bounded_actuation_required': 'uniform coefficient norm and drive-time budget; unlimited amplification recovers all eighteen finite currents',
        'exact_finite_SM_theorem_discharged': False,
        'new_architecture': ['full eighteen-coefficient classical disturbance-gradient controller',
            'native class-sum work phase', 'all work drives mediated by that controller',
            'fresh reference baths and counters', 'bounded count and coherent pulse actuator'],
        'unchanged_C_E_autonomously_selects_this_controller': False,
        'physical_gauge_group_and_measured_particle_parameters_derived': False}


def execute(reference_arrays, cooling_arrays, *, reference_report, cooling_report):
    # This public producer also owns the custody of the parent resource and
    # error contracts. Aggregate callers supply the shared factors; fresh
    # owner replays here validate them rather than substitute new samples.
    checked_reference, checked_reference_arrays = reference.execute()
    reference.compare_execution(reference_report, reference_arrays,
                                checked_reference, checked_reference_arrays)
    checked_cooling, checked_cooling_arrays = cooling.execute()
    cooling.compare_execution(cooling_report, cooling_arrays,
                              checked_cooling, checked_cooling_arrays)
    coefficients = _coefficients()
    words = [[_word(_generator(coefficients[k], mu)) for mu in range(18)] for k in range(CAP+1)]
    controls = np.asarray([[compiler.execute(word, 3) for word in row] for row in words])
    arrays = {'coefficients': coefficients, 'controls': controls}
    own = metrics(arrays, reference_arrays, cooling_arrays)
    preparation = [event for stage in transport._preparation_stages() for event in stage]
    inverse = compiler.inverse(preparation)
    native = reference_report['resources']['native_pulses_per_local_law']
    pulses = max(len(word) for row in words for word in row)+len(preparation)+len(inverse)+sum(native.values())
    export = reference._export_diagnostic([word for row in words for word in row]+[preparation, inverse], pulses)
    # Carry all raw-macro and intervening source norm growth on this path.
    _, macros = reference._native_costs(reference_report['resources']['native_factors'])
    raw = math.exp(pulses*math.log1p(2e-11))*math.expm1(sum(macros.values())*math.log1p(5.1e-11))
    successful = (1+SOURCE_ADMISSION)**6*(1+raw)-1
    work_error = successful+export['first_failure_amplitude_bound']
    cooling_error = cooling_report['error_budget']['one_cycle_channel_error']
    numerical = math.expm1(HORIZON*math.log1p(cooling_error)
                          +math.log1p(work_error*(2+work_error)))
    bounds = dynamics.feedback_bound(_source_data()[0], CAP, STEP)
    epsilon = float(bounds['coefficient_operator_error'])
    timeout = float(sum(distribution()[:-1]))
    ideal_error = 2*float(PARAMETER)*epsilon+2*timeout
    if not math.isfinite(numerical) or numerical > .01:
        raise ArithmeticError('complete finite history exceeds numerical admission')
    return {'exact': exact_certificate(), 'execution': own,
        'resources': {'clock_and_bath_attempts': HORIZON, 'maximum_controller_updates': CAP,
            'real_coefficient_registers': 18, 'stored_Gram_real_entries': 324,
            'dense_update_multiplications': 324, 'dense_update_additions': 306,
            'update_matrix_precomputed_from_full_source_Gram': True,
            'count_bits': CAP.bit_length(), 'timeout_flag_bits': 1,
            'raw_clock_record_bits': 3*HORIZON,
            'raw_data_repair_record_bits': 7,
            'complete_native_record_bits': 3*HORIZON+7,
            'first_export_failure_archive_bits': (HORIZON*cooling_report['resources']['maximum_native_pulses_per_cycle']+pulses).bit_length()+1,
            'retained_bath_and_counter_qubits': 15*HORIZON,
            'initial_unknown_reference_qubits': 6, 'initial_unknown_data_qubits': 3,
            'data_geometric_blank_qubits': 3, 'data_native_encoding_blank_qubits': 12,
            'data_coherent_high_syndrome_qubits': 12, 'data_blank_return_qubits': 3,
            'decoded_history_quantum_output_qubits': 15*HORIZON+24,
            'maximum_work_native_pulses': pulses,
            'maximum_history_native_pulses': HORIZON*cooling_report['resources']['maximum_native_pulses_per_cycle']+pulses,
            'classical_arithmetic': 'binary64 replay checked against exact matrix powers; bit-cost and physical precision are additional services',
            'physical_error_budget': '512 delta_cooling + delta_work + delta_controller + delta_storage/records, for supplied complete-instrument service bounds'},
        'error_budget': {'coefficient_source_admission': COEFFICIENT_ADMISSION,
            'component_source_admission': SOURCE_ADMISSION,
            'proper_code_exports': export, 'work_complete_operator_error': work_error,
            'source_clock_and_work_channel_error': numerical,
            'completed_controller_coefficient_error': epsilon,
            'ideal_timeout_probability': timeout,
            'ideal_complete_history_channel_error': ideal_error,
            'complete_source_history_vs_limiting_control_channel_error': numerical+ideal_error,
            'numerical_not_outward_rounded_or_physical': True,
            'reference_and_bath_outputs_matched_not_discarded_in_this_bound': True,
            'finite_protection_ground_error_not_needed_for_compiled_projected_reference': True}}, arrays


def compare_execution(report, arrays, fresh_report, fresh_arrays, *, reference_arrays,
                      cooling_arrays, fresh_reference_arrays, fresh_cooling_arrays):
    if type(report) is not dict or report.keys() != fresh_report.keys():
        raise ValueError('complete controller report differs')
    reference._owned_compare(metrics(arrays, reference_arrays, cooling_arrays), report['execution'])
    cooling._exact_compare(fresh_report['exact'], report['exact'])
    reference._fixed_compare(fresh_report['resources'], report['resources'])
    reference._fixed_compare(fresh_report['error_budget'], report['error_budget'])
    _schema(fresh_arrays)
    for key in arrays:
        if np.linalg.norm(arrays[key]-fresh_arrays[key]) > reference.REPLAY_ADMISSION:
            raise ValueError('retained controller source differs: '+key)


def source_bindings():
    paths = set(dynamics.source_bindings())|set(cooling.source_bindings())|{
        'model/common_source/native_response_controller.py',
        'model/common_source/test_native_response_controller.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}
