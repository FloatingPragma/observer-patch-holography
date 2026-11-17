"""Two native carriers, two interventions and complete C/E repair records.

All amplitudes come from the existing source words on the full nine-bit
carriers. The nineteen-bit collective word is evaluated by local tensor
contractions, not by a compiler call beyond its sixteen-bit workspace.
Numerical column compression retains an explicit discarded-norm allowance;
it does not replace a carrier by its ideal code or reset a syndrome.

The independent comparison is the fixed two-site correlation/exchange action.
For recorded dilations V,W, ||Ad(V)-Ad(W)||_diamond is bounded by
(||V||+||W||)||V-W||, including an arbitrary untouched reference. This is
the elementary forward Stinespring estimate; the nearest general continuity
theorem is Kretschmann/Schlingemann/Werner, quant-ph/0605009. The added result
is this actual finite source composition, not a new continuity theorem.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from tools import strict_json
from . import record_code as code, repair_field_attachment as field
from .mass_response import history as detector
from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit


ROOT = Path(__file__).resolve().parents[2]
TIMES = (F(3), F(5))
SEAMS = {'C': (2, 3), 'E': (2, 10)}
RANK = 4
COMPRESSION_ADMISSION = 1e-8
SOURCE_ADMISSION = 1e-7
REPLAY_ADMISSION = 2e-8
DIAGNOSTIC_REPLAY = 2e-10


def _freeze(value):
    raw = np.asarray(value, dtype=np.complex128)
    return np.frombuffer(raw.tobytes(), np.complex128).reshape(raw.shape)


def _seed_word():
    # X(pi/2)=iX; the scalar pulse retains the absolute preparation phase.
    return [compiler.pulse('X', math.pi/2, [0]), compiler.pulse('I', -math.pi/2, [0])]


def _dense(a, b):
    """carrier0 is the low nine bits; columns are untouched input references."""
    return np.einsum('ar,brc->bac', a, b, optimize=True).reshape(1 << 18, 4)


def _compress(value):
    """Column-space approximation of actual rows, with a measured residual.

    Fixed lexical column order avoids LAPACK singular-vector conventions.
    Tiny failure branches are compressed relative to their own amplitude.
    Their unnormalised coefficients and discarded norm both remain present.
    """
    matrix = value.reshape(512, 512, 4).transpose(1, 0, 2).reshape(512, 2048)
    scale = float(np.max(np.linalg.norm(matrix, axis=0)))
    basis = np.zeros((512, RANK), complex)
    pivots = []
    threshold = max(1e-30, scale*1e-8)
    for column in range(matrix.shape[1]):
        v = matrix[:, column].copy()
        for _ in range(2):
            if pivots:
                q = basis[:, :len(pivots)]
                v -= q@(q.conj().T@v)
        norm = float(np.linalg.norm(v))
        if norm > threshold:
            basis[:, len(pivots)] = v/norm
            pivots.append(column)
            if len(pivots) == RANK:
                break
    coefficients = basis.conj().T@matrix
    residual = float(np.linalg.norm(matrix-basis@coefficients, 'fro'))
    if not math.isfinite(residual) or residual > COMPRESSION_ADMISSION:
        raise ArithmeticError('actual native tensor compression exceeds fixed admission')
    other = coefficients.reshape(RANK, 512, 4).transpose(1, 0, 2)
    return _freeze(basis), _freeze(other), residual


def _gram(a, b):
    return np.einsum('rs,brc,bsd->cd', a.conj().T@a, b.conj(), b, optimize=True)


def _onsite(a, b, duration):
    word = field._walsh_diagonal_word(field._orbit_data()[2], duration)
    return (compiler.execute(word, 9, a),
            compiler.execute(word, 9, b.reshape(512, -1)).reshape(b.shape))


def _exchange_raw(a, b, duration, macros):
    """Execute every scratch row of the actual nineteen-bit factor word."""
    value = np.zeros((1 << 19, 4), complex)
    value[:1 << 18] = _dense(a, b)
    theta = float(duration*field.logical_parameters()['exchange_strength'])
    z = macros['Z'] if 'Z' in macros else compiler.execute([compiler.pulse('Z', theta, [0])], 1)
    scalar = macros['scalar'] if 'scalar' in macros else compiler.execute([compiler.pulse('I', -theta, [0])], 1)
    for stage in (0, 1):
        value = field._little_apply(macros['H'], (18,), 19, value)
        for bit in range(9):
            value = field._little_apply(macros['fredkin'], (bit, 9+bit, 18), 19, value)
        value = field._little_apply(macros['H'], (18,), 19, value)
        if stage == 0:
            value = field._little_apply(z, (18,), 19, value)
    value = field._little_apply(scalar, (18,), 19, value)
    return value


def _exchange(a, b, duration, macros):
    value = _exchange_raw(a, b, duration, macros)
    return [_compress(value[branch*(1 << 18):(branch+1)*(1 << 18)]) for branch in (0, 1)]


def _phase_controls(duration):
    theta = float(duration*field.logical_parameters()['exchange_strength'])
    return {'Z': compiler.execute([compiler.pulse('Z', theta, [0])], 1),
        'scalar': compiler.execute([compiler.pulse('I', -theta, [0])], 1)}


def _read(a):
    """Actual inverse encoder, pointer instrument and encoder, all512 rows."""
    decoded = code.decoding(a)
    charge = [0]*8
    for label in field.logical_parameters()['occupied_labels']:
        charge[label] = 1
    reading = detector.diagonal_read(decoded, [1-x for x in charge], targets=(0, 1, 2))
    return tuple(_freeze(code._encode(value)) for value in reading.outputs), reading.resources


def logical_target(*, first_read=True):
    p = field.logical_parameters()
    n = np.diag([0., 1., 1., 2.])
    swap = np.asarray(field.carrier_swap(0, 1, 2), complex)
    h = (2*float(p['onsite_baselines'][0])*np.eye(4)+float(p['positive_occupation_gaps'][0])*n).astype(complex)
    h += float(p['exchange_strength'])*(np.eye(4)-swap)
    unitaries = [expm(-1j*float(t)*h) for t in TIMES]
    projections = (np.diag([1., 0., 1., 0.]), np.diag([0., 1., 0., 1.]))
    if first_read:
        return {(u, v): projections[v]@unitaries[1]@projections[u]@unitaries[0]
                for u in (0, 1) for v in (0, 1)}
    return {v: projections[v]@unitaries[1]@unitaries[0] for v in (0, 1)}


def number_identifiability_control():
    """The chosen number history cannot determine the uniform onsite gap.

    This is a comparison of declared target actions, not additional source
    executions or a claim that the complete four-input dilations are equal.
    On fixed N=1 the uniform onsite term is scalar; exchange and both reads
    preserve that sector, so it contributes only one common phase.
    """
    import sympy as s
    number = s.diag(0, 1, 1, 2)
    one = s.diag(0, 1, 1, 0)
    swap = field.carrier_swap(0, 1, 2)
    projectors = (s.diag(1, 0, 1, 0), s.diag(0, 1, 0, 1))
    if (number*one != one or number*swap != swap*number
            or any(number*p != p*number for p in projectors)):
        raise ArithmeticError('exact fixed-charge number-history scalarization failed')
    j = float(field.logical_parameters()['exchange_strength'])
    laws = []
    for delta in (F(0), F(1, 256), F(1)):
        h = j*(np.eye(4)-np.asarray(swap, complex))+float(delta)*np.asarray(number, complex)
        u, v = (expm(-1j*float(t)*h) for t in TIMES)
        laws.append(np.array([[np.linalg.norm(np.asarray(q, complex)@v@np.asarray(p, complex)@u[:, 1])**2
                               for q in projectors] for p in projectors]))
    if max(np.max(np.abs(value-laws[0])) for value in laws) > 1e-12:
        raise ArithmeticError('uniform onsite number-history control failed')
    return {'fixed_charge': 1, 'uniform_onsite_restriction_exact': '(2*E0+Delta)*I_N1',
        'exchange_and_both_number_reads_preserve_this_sector_exactly': True,
        'comparison_uniform_onsite_coefficients': ['0', '1/256', '1'],
        'comparison_is_target_matrix_control_not_additional_native_source_execution': True,
        'number_law_depends_only_on_J_and_declared_intervention_parameters': True,
        'complete_four_input_dilations_still_contain_relative_charge_phases': True,
        'reported_number_history_identifies_native_mass_gap': False}


def _target_for(label):
    # stop, firstread, Ccoin,Craw,Ecoin,Eraw, finalread, firstancilla, secondancilla
    stop, u, cc, cr, ec, er, v, _, _ = map(int, label)
    if stop:
        return np.zeros((4, 4), complex)
    p = code.branch_weight('C', SEAMS['C'], cc, cr)*code.branch_weight('E', SEAMS['E'], ec, er)
    return math.sqrt(float(p))*logical_target()[(u, v)]


def _factor_target_bound(a, b, target):
    """Stable Frobenius upper estimate, avoiding subtraction of near-one norms."""
    if not np.any(target):
        return math.sqrt(max(0., float(np.trace(_gram(a, b)).real)))
    e = code.ideal_code_matrix()[:, (0, 1)]
    x = e.conj().T@a
    y = np.einsum('bi,brc->irc', e.conj(), b, optimize=True)
    projection = np.einsum('ar,brc->bac', x, y, optimize=True).reshape(4, 4)
    acode = e@x
    bperp = b-np.einsum('bi,irc->brc', e, y, optimize=True)
    return float(np.linalg.norm(projection-target, 'fro')
        +np.linalg.norm(a-acode, 2)*np.linalg.norm(b.reshape(512, -1), 'fro')
        +np.linalg.norm(acode, 2)*np.linalg.norm(bperp.reshape(512, -1), 'fro'))


def _schema(arrays):
    menu = {'A', 'B', 'terminal', 'labels', 'compression_residuals',
            'exchange_input_A', 'exchange_input_B', 'exchange_output_A', 'exchange_output_B',
            'fredkin_macro', 'ancilla_H', 'exchange_Z', 'exchange_scalar', 'prepared_seed'}
    if type(arrays) is not dict or set(arrays) != menu:
        raise ValueError('complete native tensor/record menu differs')
    a, b, terminal, labels, residuals = (arrays[k] for k in ('A', 'B', 'terminal', 'labels', 'compression_residuals'))
    for value in arrays.values():
        if np.ma.isMaskedArray(value) or type(value) is not np.ndarray or not np.isfinite(value).all():
            raise ValueError('original finite unmasked native arrays required')
    if (a.dtype != np.dtype(complex) or b.dtype != np.dtype(complex) or a.ndim != 3
            or a.shape[1:] != (512, RANK) or b.ndim != 4 or b.shape[1:] != (512, RANK, 4)
            or terminal.dtype != np.dtype('int64') or terminal.ndim != 2 or terminal.shape[1] != 2
            or labels.dtype != np.dtype('int64') or labels.shape != (len(terminal), 9)
            or residuals.dtype != np.dtype(float) or residuals.shape != (6,)
            or np.any(residuals < 0) or np.any(residuals > COMPRESSION_ADMISSION)
            or np.any(terminal < 0) or np.any(terminal[:, 0] >= len(a)) or np.any(terminal[:, 1] >= len(b))):
        raise ValueError('complete native tensor/record schema differs')
    if (not np.array_equal(labels, record_menu())
            or set(terminal[:, 0]) != set(range(len(a))) or set(terminal[:, 1]) != set(range(len(b)))):
        raise ValueError('complete native sequential record menu or tensor authority differs')
    for name, shape in {'exchange_input_A': (3, 512, RANK),
            'exchange_input_B': (3, 512, RANK, 4), 'exchange_output_A': (6, 512, RANK),
            'exchange_output_B': (6, 512, RANK, 4), 'fredkin_macro': (8, 8), 'ancilla_H': (2, 2),
            'exchange_Z': (3, 2, 2), 'exchange_scalar': (3, 2, 2)}.items():
        if arrays[name].dtype != np.dtype(complex) or arrays[name].shape != shape:
            raise ValueError('native retained stage/macro schema differs: '+name)
    seed = arrays['prepared_seed']
    if (seed.dtype != np.dtype(complex) or seed.shape != (4,)
            or np.linalg.norm(seed-np.array([0., 1., 0., 0.])) > field.LOCAL_MACRO_GATE):
        raise ValueError('source seed preparation differs, including its absolute phase')


def record_menu():
    rows = [(1, -1, -1, -1, -1, -1, -1, 1, -1)]
    for u in (0, 1):
        rows.append((1, u, -1, -1, -1, -1, -1, 0, 1))
        for cc in (0, 1):
            for cr in range(4):
                if cr not in circuit.record_values('C', cc):
                    rows.append((2, u, cc, cr, -1, -1, -1, 0, 0))
                    continue
                for ec in (0, 1):
                    for er in range(8):
                        if er not in circuit.record_values('E', ec):
                            rows.append((3, u, cc, cr, ec, er, -1, 0, 0))
                        else:
                            for v in (0, 1):
                                rows.append((0, u, cc, cr, ec, er, v, 0, 0))
    return np.asarray(rows, np.int64)


def metrics(arrays):
    _schema(arrays)
    a, b, indices, labels = (arrays[k] for k in ('A', 'B', 'terminal', 'labels'))
    complete = np.zeros((4, 4), complex)
    weights = np.zeros((len(indices), 4))
    errors = []
    joint = np.zeros((2, 2))
    failure = 0.
    padding = 0.
    native_rows = np.zeros(512, bool)
    native_rows[[circuit.encode(state) for state in field.native.states()]] = True
    seed = arrays['prepared_seed']
    for index, ((i, j), label) in enumerate(zip(indices, labels, strict=True)):
        aa, bb = a[i], b[j]
        gram = _gram(aa, bb)
        complete += gram
        weights[index] = np.maximum(0., np.diag(gram).real)
        seed_weight = max(0., float(np.vdot(seed, gram@seed).real))
        errors.append(_factor_target_bound(aa, bb, _target_for(label)))
        if label[0]:
            failure += seed_weight
        else:
            joint[label[1], label[6]] += seed_weight
        invalid_first = _gram(aa[~native_rows], bb)
        invalid_second = _gram(aa[native_rows], bb[~native_rows])
        padding += float(np.vdot(seed, (invalid_first+invalid_second)@seed).real)
    upper = math.sqrt(math.fsum(value*value for value in errors))
    defect = float(np.linalg.norm(complete-np.eye(4), 2))
    if upper > SOURCE_ADMISSION or defect > SOURCE_ADMISSION:
        raise ValueError('complete native history exceeds unchanged input/reference admission')
    ideal_joint = np.array([[float(np.vdot(value[:, 1], value[:, 1]).real)
                             for (a0, b0), value in logical_target().items() if a0 == u]
                            for u in (0, 1)])
    no_first = np.array([float(np.vdot(value[:, 1], value[:, 1]).real)
                        for value in logical_target(first_read=False).values()])
    return {'complete_input_columns': 4, 'terminal_records': len(indices),
        'complete_input_Gram_defect_diagnostic': defect,
        'all_record_dilation_Frobenius_error_upper_diagnostic': upper,
        'source_record_probability_mass': float(np.vdot(seed, complete@seed).real),
        'detector_joint_probabilities': joint.tolist(), 'independent_target_detector_joint': ideal_joint.tolist(),
        'retained_halt_probability_diagnostic': float(failure),
        'all512_row_padding_probability_diagnostic': max(0., padding),
        'first_read_backaction_final_marginal_L1_control': float(np.sum(np.abs(ideal_joint.sum(axis=0)-no_first))),
        'compression_residual_sum_diagnostic': float(np.sum(arrays['compression_residuals']))}


def execute():
    actual = code.code_matrix()[:, (0, 1)]
    a = np.zeros((512, RANK), complex); a[:, :2] = actual
    b = np.zeros((512, RANK, 4), complex)
    for column in range(4):
        b[:, column % 2, column] = actual[:, column//2]
    macros = {'fredkin': field.fredkin_source()['matrix'],
              'H': compiler.execute([compiler.pulse('H', math.pi/2, [0])], 1)}
    pool_a, pool_b, terminal, labels, residuals = [], [], [], [], []
    stage_input_a, stage_input_b, stage_output_a, stage_output_b = [], [], [], []
    stage_z, stage_scalar = [], []
    known_a, known_b = {}, {}

    def save(aa, bb, label):
        ids = []
        for value, pool, known in ((aa, pool_a, known_a), (bb, pool_b, known_b)):
            identity = hashlib.sha256(np.ascontiguousarray(value).tobytes()).digest()
            if identity not in known:
                known[identity] = len(pool)
                pool.append(value)
            ids.append(known[identity])
        terminal.append(ids); labels.append(label)

    a, b = _onsite(a, b, TIMES[0])
    stage_input_a.append(a); stage_input_b.append(b)
    controls = _phase_controls(TIMES[0]); stage_z.append(controls['Z']); stage_scalar.append(controls['scalar'])
    first = _exchange(a, b, TIMES[0], {**macros, **controls})
    stage_output_a.extend(item[0] for item in first); stage_output_b.extend(item[1] for item in first)
    residuals.extend(item[2] for item in first)
    save(*first[1][:2], (1, -1, -1, -1, -1, -1, -1, 1, -1))
    readings, pointer_resources = _read(first[0][0])
    for u, measured in enumerate(readings):
        aa, bb = _onsite(measured, first[0][1], TIMES[1])
        stage_input_a.append(aa); stage_input_b.append(bb)
        controls = _phase_controls(TIMES[1]); stage_z.append(controls['Z']); stage_scalar.append(controls['scalar'])
        second = _exchange(aa, bb, TIMES[1], {**macros, **controls})
        stage_output_a.extend(item[0] for item in second); stage_output_b.extend(item[1] for item in second)
        residuals.extend(item[2] for item in second)
        save(*second[1][:2], (1, u, -1, -1, -1, -1, -1, 0, 1))
        ca, cb = second[0][:2]
        cbranches = code.repair_branches('C', SEAMS['C'], ca)
        ebranches = code.repair_branches('E', SEAMS['E'], cb.reshape(512, -1))
        valid_c = [branch for branch in cbranches if not branch['abort']]
        combined = np.concatenate([branch['amplitude'] for branch in valid_c], axis=1)
        final_readings, _ = _read(combined)
        final_index = {id(branch): index for index, branch in enumerate(valid_c)}
        for cbranch in cbranches:
            cc, cr = cbranch['coin'], cbranch['record_code']
            if cbranch['abort']:
                save(cbranch['amplitude'], cb, (2, u, cc, cr, -1, -1, -1, 0, 0))
                continue
            offset = RANK*final_index[id(cbranch)]
            for ebranch in ebranches:
                ec, er = ebranch['coin'], ebranch['record_code']
                eb = ebranch['amplitude'].reshape(512, RANK, 4)
                if ebranch['abort']:
                    save(cbranch['amplitude'], eb, (3, u, cc, cr, ec, er, -1, 0, 0))
                else:
                    for v in (0, 1):
                        save(final_readings[v][:, offset:offset+RANK], eb, (0, u, cc, cr, ec, er, v, 0, 0))
    arrays = {'A': np.asarray(pool_a), 'B': np.asarray(pool_b),
        'terminal': np.asarray(terminal, np.int64), 'labels': np.asarray(labels, np.int64),
        'compression_residuals': np.asarray(residuals, float),
        'exchange_input_A': np.asarray(stage_input_a), 'exchange_input_B': np.asarray(stage_input_b),
        'exchange_output_A': np.asarray(stage_output_a), 'exchange_output_B': np.asarray(stage_output_b),
        'fredkin_macro': macros['fredkin'], 'ancilla_H': macros['H'],
        'exchange_Z': np.asarray(stage_z), 'exchange_scalar': np.asarray(stage_scalar),
        'prepared_seed': np.kron([1., 0.], compiler.execute(_seed_word(), 1, [1., 0.]))}
    report = {'schema': 'oph.native_field_history.v1',
        'classification': 'complete bounded native source composition; no measured particle inputs',
        'compatibility': 'independently specified tensor product of the matching codes at completed control boundaries',
        'carrier': {'native_bits_per_carrier': 9, 'carriers': 2, 'scratch_bits_per_exchange': 1,
            'logical_labels': [0, 1], 'arbitrary_untouched_reference_included': True,
            'reported_seed': 'actual source X and scalar-phase preparation of |10> from two logical blanks, followed by actual encoders',
            'actual512_row_factor_composition': True, 'ideal_encoded_arrays_used_as_source': False},
        'protocol': {'clock_parameter': 'supplied common model-coordinate parameter, not physical seconds',
            'intervals': list(map(str, TIMES)), 'native_seams': {k: list(v) for k, v in SEAMS.items()},
            'order': ['source preparation', 'onsite and exchange', 'first detector',
                'onsite and exchange', 'complete C repair/recovery carrier0',
                'complete E repair/recovery carrier1', 'second detector carrier0'],
            'ancilla_failure_and_raw_abort': 'retained unnormalised halt outputs; later unperformed records use -1',
            'no_syndrome_reset_or_code_projection': True,
            'compression': 'fixed-rank column-space approximation of actual native arrays with explicit discarded-norm allowance'},
        'admissions': {'all_record_source_input_reference_error': SOURCE_ADMISSION,
            'each_tensor_compression_Frobenius_allowance': COMPRESSION_ADMISSION,
            'fresh_complete_instrument_operator_replay': REPLAY_ADMISSION,
            'reported_seed_vector_preparation': field.LOCAL_MACRO_GATE},
        'composable_bound': {'theorem': '||(V tensor I)rho(V* tensor I)-(W tensor I)rho(W* tensor I)||1 <= (||V||+||W||)||V-W||',
            'compressed_outputs_reference_error_admission': SOURCE_ADMISSION,
            'six_compressions_additional_dilation_allowance': 6*COMPRESSION_ADMISSION,
            'fixed_complete_history_diamond_allowance': (2+SOURCE_ADMISSION+6*COMPRESSION_ADMISSION)*(SOURCE_ADMISSION+6*COMPRESSION_ADMISSION),
            'fixed_prepared_seed_output_trace_allowance':
                (2+SOURCE_ADMISSION+6*COMPRESSION_ADMISSION+(1+SOURCE_ADMISSION+6*COMPRESSION_ADMISSION)*field.LOCAL_MACRO_GATE)
                *(SOURCE_ADMISSION+6*COMPRESSION_ADMISSION+(1+SOURCE_ADMISSION+6*COMPRESSION_ADMISSION)*field.LOCAL_MACRO_GATE),
            'discarded_norm_adds_by_triangle_and_isometric_continuation': True,
            'numerical_admissions_are_outward_certificates_or_physical_service_errors': False,
            'additional_physical_error_contract': 'if every required proper service has diamond error eta, add at most required_service_events*eta; preparation/readout/memory errors add separately',
            'primary_source': 'https://arxiv.org/abs/quant-ph/0605009'},
        'execution': metrics(arrays), 'resources': resource_counts(pointer_resources),
        'number_history_identifiability_control': number_identifiability_control(),
        'scope': {'one_bounded_two_carrier_history': True,
            'reported_number_history_identifies_native_mass_gap': False,
            'C_E_repair_selects_the_correlation_exchange_controls_or_schedule': False,
            'complete_native_hardware_protection_load_export_history_executed': False,
            'proper_control_services_blank_preparation_and_readout_are_supplied': True,
            'physical_clock_spacetime_SM_or_observed_mass_law_derived': False}}
    return report, arrays


def resource_counts(pointer_resources):
    encoders = 2*code.encoding_resources()['native_pulses']+compiler.resources(_seed_word())['native_pulses']
    onsite = 2*sum(compiler.resources(field._walsh_diagonal_word(field._orbit_data()[2], t))['native_pulses'] for t in TIMES)
    exchange = 2*(18*field.fredkin_source()['resources']['native_pulses']+6)
    detector_cost = 2*(code.decoding_resources()['native_pulses']+pointer_resources['native_pulses']+code.encoding_resources()['native_pulses'])
    repair = {}
    for law, seam in SEAMS.items():
        recover = max(code.recovery_resources(law, seam, coin, raw)['native_pulses']
                      for coin in (0, 1) for raw in circuit.record_values(law, coin))
        repair[law] = code.repair_resources(law, seam)['native_pulses']+recover
    total = encoders+onsite+exchange+detector_cost+sum(repair.values())
    return {'basis': 'one successful physical source-word path; all-input host queries and alternative records are not separate trials',
        'preparation_native_pulses': encoders, 'onsite_native_pulses': onsite,
        'exchange_native_pulses': exchange, 'detector_native_pulses': detector_cost,
        'repair_and_record_selected_recovery_native_pulses_upper': repair,
        'one_history_native_pulses_upper': total,
        'required_proper_code_load_export_service_events_upper': 2*total,
        'proper_load_export_services_executed_in_this_producer': False,
        'supplied_logical_input_register_qubits': 6,
        'fixed_zero_padding_bits_in_logical_input_registers': 4,
        'reported_seed_logical_input_blanks': 2,
        'fresh_blank_code_qubits': 12, 'fresh_exchange_ancilla_qubits': 2,
        'fresh_detector_pointer_qubits': 2, 'fresh_native_coin_and_raw_record_qubits': 7,
        'maximum_retained_classical_record_bits': 11,
        'earlier_halt_paths_require_no_later_source_operations': True}


def _factor_distance(a, b, c, d):
    """Basis-independent upper bound between two represented four-column maps."""
    q, _ = np.linalg.qr(c, mode='reduced')
    projection = q.conj().T@a
    x = np.einsum('ir,brc->bic', projection, b, optimize=True)
    y = np.einsum('ir,brc->bic', q.conj().T@c, d, optimize=True)
    return float(np.linalg.norm((x-y).reshape(512, -1), 'fro')
        +np.linalg.norm(a-q@projection, 2)*np.linalg.norm(b.reshape(512, -1), 'fro'))


def verify_compressions(arrays):
    """Bind each discarded norm to its own retained input/output binary arrays."""
    _schema(arrays)
    macros = {'fredkin': arrays['fredkin_macro'], 'H': arrays['ancilla_H']}
    if (np.linalg.norm(macros['fredkin']-field._fredkin_target(), 2) > field.LOCAL_MACRO_GATE
            or np.linalg.norm(macros['H']-1j*np.array([[1., 1.], [1., -1.]])/math.sqrt(2), 2) > field.LOCAL_MACRO_GATE):
        raise ValueError('retained complete local source macro exceeds unchanged operator admission')
    for stage, duration in enumerate((TIMES[0], TIMES[1], TIMES[1])):
        theta = float(duration*field.logical_parameters()['exchange_strength'])
        controls = {'Z': arrays['exchange_Z'][stage], 'scalar': arrays['exchange_scalar'][stage]}
        if (np.linalg.norm(controls['Z']-np.diag([np.exp(1j*theta), np.exp(-1j*theta)]), 2) > field.LOCAL_MACRO_GATE
                or np.linalg.norm(controls['scalar']-np.exp(-1j*theta)*np.eye(2), 2) > field.LOCAL_MACRO_GATE):
            raise ValueError('retained complete phase control exceeds unchanged operator admission')
        incoming_norm = float(np.linalg.norm(_dense(arrays['exchange_input_A'][stage], arrays['exchange_input_B'][stage]), 'fro'))
        value = _exchange_raw(arrays['exchange_input_A'][stage], arrays['exchange_input_B'][stage], duration, {**macros, **controls})
        for outcome in (0, 1):
            index = 2*stage+outcome
            raw = value[outcome*(1 << 18):(outcome+1)*(1 << 18)]
            actual = _dense(arrays['exchange_output_A'][index], arrays['exchange_output_B'][index])
            residual = float(np.linalg.norm(raw-actual, 'fro'))
            # Cancellation into a tiny failed-return output still incurs
            # arithmetic at the unconditioned input scale, not its tiny norm.
            provenance = 1024*np.finfo(float).eps*max(incoming_norm, 1e-20)
            if residual > COMPRESSION_ADMISSION or abs(residual-arrays['compression_residuals'][index]) > provenance:
                raise ValueError('retained compression residual differs from its own binary source arrays')


def source_bindings():
    names = set(field.source_bindings())|{
        'model/common_source/native_field_history.py', 'model/common_source/test_native_field_history.py',
        'tools/strict_json.py', 'model/observer/source_discrimination/experiment.py'}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual):
        raise ValueError(path+' type differs')
    if type(expected) is dict:
        if expected.keys() != actual.keys():
            raise ValueError(path+' complete keys differ')
        for key in expected:
            _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual):
            raise ValueError(path+' complete list differs')
        for index, (left, right) in enumerate(zip(expected, actual, strict=True)):
            _compare(left, right, path+f'[{index}]')
    elif type(expected) is float:
        exact = path.startswith(('receipt.admissions.', 'receipt.composable_bound.'))
        if (not math.isfinite(actual) or (expected != actual if exact
                else not math.isclose(expected, actual, rel_tol=2e-9, abs_tol=2e-10))):
            raise ValueError(path+' numerical value differs')
    elif expected != actual:
        raise ValueError(path+' exact value differs')


def compare_execution(report, arrays, fresh, fresh_arrays):
    own = metrics(arrays)
    for key, value in own.items():
        if not key.endswith('_diagnostic'):
            continue
        retained = report['execution'].get(key)
        positive_sum = key in ('retained_halt_probability_diagnostic',
            'all512_row_padding_probability_diagnostic', 'compression_residual_sum_diagnostic')
        tolerance = 0. if positive_sum else 64*np.finfo(float).eps
        if (type(retained) is not float or not math.isfinite(retained)
                or not math.isclose(value, retained, rel_tol=2e-10 if positive_sum else 0., abs_tol=tolerance)):
            raise ValueError('retained diagnostic differs from its own binary arrays: '+key)
    _compare(own, report['execution'], 'retained.execution')
    # These are diagnostics of their respective binary arrays. Their fixed
    # admissions and complete targets are checked separately on both arrays.
    left, right = deepcopy(report), deepcopy(fresh)
    for key in own:
        if key.endswith('_diagnostic'):
            left['execution'][key] = right['execution'][key] = 0.
    _compare(right, left)
    error = []
    for (i, j), (k, ell) in zip(arrays['terminal'], fresh_arrays['terminal'], strict=True):
        error.append(_factor_distance(arrays['A'][i], arrays['B'][j], fresh_arrays['A'][k], fresh_arrays['B'][ell]))
    if math.sqrt(math.fsum(value*value for value in error)) > REPLAY_ADMISSION:
        raise ValueError('complete retained source dilation differs in operator/reference norm')
    for prefix in ('exchange_input_', 'exchange_output_'):
        for a, b, c, d in zip(arrays[prefix+'A'], arrays[prefix+'B'], fresh_arrays[prefix+'A'], fresh_arrays[prefix+'B'], strict=True):
            if _factor_distance(a, b, c, d) > REPLAY_ADMISSION:
                raise ValueError('retained native stage source operator differs')


def check(path):
    path = Path(path)
    retained = strict_json.load(path)
    bindings = source_bindings()
    if type(retained) is not dict or retained.get('source_bindings') != bindings:
        raise ValueError('complete native history closed source custody differs')
    artifact = path.with_suffix('.npz')
    if retained.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('complete native history binary custody differs')
    with np.load(artifact, allow_pickle=False) as saved:
        if len(saved.files) != len(set(saved.files)):
            raise ValueError('duplicate native history array names')
        arrays = {key: saved[key] for key in saved.files}
    verify_compressions(arrays)
    fresh, fresh_arrays = execute()
    scientific = {key: value for key, value in retained.items() if key not in ('source_bindings', 'artifact')}
    compare_execution(scientific, arrays, fresh, fresh_arrays)
    if source_bindings() != bindings:
        raise RuntimeError('native history source changed during replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('complete two-carrier native field/repair/detector history replay: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('fresh native history receipt and binary paths required')
    bindings = source_bindings()
    report, arrays = execute()
    verify_compressions(arrays)
    if source_bindings() != bindings:
        raise RuntimeError('native history source changed during execution')
    with artifact.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    report.update(source_bindings=bindings,
        artifact={'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()})
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print('complete two-carrier native field/repair/detector history produced')


if __name__ == '__main__':
    main()
