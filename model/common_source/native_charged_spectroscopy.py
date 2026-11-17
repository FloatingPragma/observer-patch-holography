"""Charged poles inferred from complete native-carrier interference records.

A dimer and a fermionic reference have six states at total occupation one or
two. They fit the existing six-dimensional recorded-repair code. Fixed local
occupation preparation and even transfer quadratures access both dimer
frequencies without preparing a spectral eigenvector. The source evolution
uses onsite and ordinary-exchange words; no pole occurs in its pulse angles.

The finite-step pole is a Floquet energy. A separately stated Strang bound
controls its difference from the continuous dimer pole. The clock, encoded
action, coherent services and reference arm remain declared apparatus inputs.
Neither occupation charge nor a spectator SM representation identifies a
physical particle species or measured mass.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_charged_spectrum as spectrum
from . import record_code as code
from .mass_response import history
from ..observer.source_scalar import compiler
from ..observer.source_discrimination import circuit


ROOT = Path(__file__).resolve().parents[2]
FOCK_STATES = (1, 2, 3, 4, 5, 6)
SPACING = F(8)
SAMPLE_NUMBERS = tuple(range(6))
STEPS = (F(1, 4), F(1, 16), F(1, 64))
SEAMS = {'C': (2, 3), 'E': (2, 10)}
SOURCE_ADMISSION = 2e-9
COMPRESSION_ADMISSION = 1e-9
REPLAY_ADMISSION = 1e-9
RANK = 6


def _freeze(value):
    raw = np.asarray(value, dtype=np.complex128)
    return np.frombuffer(raw.tobytes(), np.complex128).reshape(raw.shape)


def _local(matrix, state):
    """Actual low-three-bit source matrix, retaining all sixty-four blocks."""
    shape = state.shape
    return np.einsum('ab,sbc->sac', matrix, state.reshape(64, 8, -1),
                     optimize=True).reshape(shape)


def frame_target():
    matrix = np.zeros((8, 8))
    matrix[(*FOCK_STATES, 0, 7), np.arange(8)] = 1.
    return matrix


@lru_cache(maxsize=1)
def frame_word():
    return tuple(compiler.compile_real_orthogonal(frame_target()))


def preparation_target():
    # Fock3=|dimer full, reference empty>, Fock5=|site0, reference full>.
    value = np.zeros(8)
    value[FOCK_STATES.index(3)] = value[FOCK_STATES.index(5)] = 1/math.sqrt(2)
    return value


@lru_cache(maxsize=1)
def preparation_word():
    return tuple(compiler.compile_real_state(preparation_target()))


def cycle_word(mass, exchange, step):
    """Strang step from m(n0-n1) and J(I-S01), without a spectral gap."""
    if any(type(x) not in (int, F) or x <= 0 for x in (mass, exchange, step)):
        raise ValueError('positive exact rational source controls required')
    half = [compiler.pulse('Z', float(mass*step/4), [0]),
            compiler.pulse('Z', -float(mass*step/4), [1])]
    middle = [compiler.pulse('I', -float(exchange*step/2), [0])]
    for label in ('XXI', 'YYI', 'ZZI'):
        middle += compiler.pauli_rotation(label, float(exchange*step/2))
    return half+middle+half


def hamiltonian(mass, exchange):
    """Independent full-Fock action used only for validation."""
    result = np.zeros((8, 8), complex)
    for x in range(8):
        n0, n1 = x & 1, (x >> 1) & 1
        result[x, x] = float(mass)*(n0-n1)
        if n0 != n1:
            result[x, x] += float(exchange)
            result[x ^ 3, x] -= float(exchange)
    return result


@lru_cache(maxsize=9)
def ideal_cycle(mass, exchange, step):
    # The exact two-state exchange block is exp(-i J dt)(cos(J dt) I
    # + i sin(J dt) X). Round the analytic entries once, independently of
    # the platform's matrix-exponential implementation. Powers amplify a
    # one-ulp difference enough to change a near-zero norm diagnostic.
    jt, mt = map(s.Rational, (exchange*step, mass*step))
    phase = s.exp(-s.I*jt)
    diagonal = {sign: complex((phase*s.cos(jt)*s.exp(-s.I*mt*sign)).evalf(60))
                for sign in (-1, 1)}
    off = complex((s.I*phase*s.sin(jt)).evalf(60))
    result = np.eye(8, dtype=complex)
    for x in range(8):
        sign = (x & 1)-((x >> 1) & 1)
        if sign:
            result[x, x], result[x ^ 3, x] = diagonal[sign], off
    return _freeze(result)


def transfer_operator(axis):
    """Q=T+T* or -i(T-T*), T=c1* c_reference, hence [Q,N]=0."""
    if axis not in ('X', 'Y'):
        raise ValueError('an X or Y transfer quadrature required')
    t = np.zeros((8, 8), complex)
    for x in range(8):
        if x & 4 and not x & 2:
            t[x ^ 6, x] = 1.  # Adjacent modes1,2; spectator parity cancels.
    return t+t.conj().T if axis == 'X' else -1j*(t-t.conj().T)


def detector_basis_word(axis):
    if axis == 'X':
        requests = (('XY', math.pi/8), ('YX', -math.pi/8))
    elif axis == 'Y':
        requests = (('XX', -math.pi/8), ('YY', -math.pi/8))
    else:
        raise ValueError('an X or Y transfer quadrature required')
    return [event for label, angle in requests
            for event in compiler.pauli_rotation(label, angle)]


def ideal_detectors(axis):
    q = transfer_operator(axis)
    common = np.eye(8)/math.sqrt(2)+(.5-1/math.sqrt(2))*(q@q)
    return common+q/2, common-q/2


def _prefix():
    """Actual encoding, complete C/E instruments and unnormalised aborts."""
    raw = np.zeros((512, 7), complex)
    raw[:6, :6] = np.eye(6)
    raw[:8, 6] = compiler.execute(preparation_word(), 3, np.eye(8)[:, 0])
    initial = code._encode(raw)
    crows = code.repair_branches('C', SEAMS['C'], initial)
    valid = [row for row in crows if not row['abort']]
    combined = np.concatenate([row['amplitude'] for row in valid], axis=1)
    erows = code.repair_branches('E', SEAMS['E'], combined)
    success, halt, slabels, hlabels = [], [], [], []
    for row in crows:
        if row['abort']:
            halt.append(row['amplitude'])
            hlabels.append((1, row['coin'], row['record_code'], -1, -1))
    for i, crow in enumerate(valid):
        for erow in erows:
            value = erow['amplitude'][:, 7*i:7*(i+1)]
            label = (crow['coin'], crow['record_code'], erow['coin'], erow['record_code'])
            if erow['abort']:
                halt.append(value); hlabels.append((2, *label))
            else:
                success.append(value); slabels.append(label)
    joined = np.concatenate(success, axis=1)
    decoded = code.decoding(joined).reshape(512, len(success), 7).transpose(1, 0, 2)
    # A fixed, lexical Gram-Schmidt compression of the actual decoded arrays.
    matrix = decoded.transpose(1, 0, 2).reshape(512, -1)
    basis = np.zeros((512, RANK), complex)
    found = 0
    for column in range(matrix.shape[1]):
        v = matrix[:, column].copy()
        for _ in range(2):
            v -= basis[:, :found]@(basis[:, :found].conj().T@v)
        norm = float(np.linalg.norm(v))
        if norm > 1e-10:
            basis[:, found] = v/norm
            found += 1
            if found == RANK:
                break
    if found != RANK:
        raise ArithmeticError('complete six-input prefix rank was lost')
    coefficients = (basis.conj().T@matrix).reshape(RANK, len(success), 7).transpose(1, 0, 2)
    residual = float(np.linalg.norm(matrix-basis@(basis.conj().T@matrix), 'fro'))
    if residual > COMPRESSION_ADMISSION:
        raise ArithmeticError('native prefix compression exceeds fixed admission')
    return {'prefix_decoded': np.asarray(decoded), 'prefix_basis': basis,
        'prefix_coefficients': coefficients, 'prefix_halts': np.asarray(halt),
        'success_labels': np.asarray(slabels, np.int64), 'halt_labels': np.asarray(hlabels, np.int64),
        'prepared_input': raw[:, 6], 'compression_residual': np.array([residual])}


def _labels():
    success, halt = [], []
    for cc in (0, 1):
        for cr in range(4):
            if cr not in circuit.record_values('C', cc):
                halt.append((1, cc, cr, -1, -1))
            else:
                for ec in (0, 1):
                    for er in range(8):
                        if er in circuit.record_values('E', ec):
                            success.append((cc, cr, ec, er))
                        else:
                            halt.append((2, cc, cr, ec, er))
    # Actual executor stores first-step halts before second-step halts.
    halt.sort(key=lambda row: row[0])
    return np.asarray(success, np.int64), np.asarray(halt, np.int64)


def _source_operators():
    j, masses = spectrum.native_parameters()
    readings = [history.diagonal_read(np.eye(8, dtype=complex), [.5, 1., 0., .5],
                targets=(1, 2), basis_word=detector_basis_word(axis)) for axis in ('X', 'Y')]
    arrays = {'encoder': code._encode(np.eye(512, dtype=complex)),
        'frame': compiler.execute(frame_word(), 3),
        'inverse_frame': compiler.execute(compiler.inverse(frame_word()), 3),
        'cycles': np.asarray([[compiler.execute(cycle_word(m, j, dt), 3)
                              for dt in STEPS] for m in masses]),
        'detectors': np.asarray([read.outputs for read in readings])}
    return arrays, readings[0].resources


def infer_poles(contrasts, *, sampling_interval, contrast_error=0.):
    """Infer two poles from four records; later records are held out.

    The inputs contain neither Hamiltonian parameters nor target masses.
    The energy window is the open principal Nyquist interval. Coefficients
    solve s[k+2]=p*s[k+1]+q*s[k]. Root disks follow a direct Rouche bound.
    Numerical evaluation of this bound is an admission calculation, not an
    outward floating certificate. contrast_error bounds each complex sample.
    """
    if np.ma.isMaskedArray(contrasts):
        raise ValueError('original unmasked contrast records required')
    values = np.asarray(contrasts)
    if (values.shape != (6,) or values.dtype.kind not in 'fc'
            or not np.isfinite(values).all() or type(sampling_interval) not in (int, float)
            or not math.isfinite(sampling_interval) or sampling_interval <= 0
            or type(contrast_error) not in (int, float) or not math.isfinite(contrast_error)
            or contrast_error < 0):
        raise ValueError('six finite contrasts, positive sampling interval and nonnegative error required')
    matrix = np.array([[values[1], values[0]], [values[2], values[1]]])
    try:
        inverse = float(np.linalg.norm(np.linalg.inv(matrix), 2))
    except np.linalg.LinAlgError as exc:
        raise ValueError('record recurrence is singular') from exc
    margin = 1-2*inverse*contrast_error
    if margin <= 0:
        raise ValueError('record uncertainty destroys recurrence invertibility')
    p, q = np.linalg.solve(matrix, values[2:4])
    roots = np.roots([1., -p, -q])
    separation = float(abs(roots[0]-roots[1]))
    if separation < .25 or np.min(np.abs(roots)) < .5:
        raise ValueError('record poles are degenerate or outside the admitted annulus')
    coefficient_error = inverse*(math.sqrt(2)+2*float(np.linalg.norm([p, q])))*contrast_error/margin
    radius = 4*coefficient_error*(float(np.max(np.abs(roots)))+1)/separation
    if radius > min(separation/4, float(np.min(np.abs(roots)))/2):
        raise ValueError('record uncertainty does not isolate both spectral roots')
    if np.max(np.abs(np.abs(roots)-1)) > max(radius, 1e-10):
        raise ValueError('record poles are inconsistent with unitary source evolution')
    energies = -np.angle(roots)/sampling_interval
    order = np.argsort(energies)
    roots, energies = roots[order], energies[order]
    error = 2*radius/(float(np.min(np.abs(roots)))*sampling_interval)
    if np.max(np.abs(energies))+error >= math.pi/sampling_interval:
        raise ValueError('energy uncertainty crosses the Nyquist boundary')
    if energies[0]+error >= 0 or energies[1]-error <= 0:
        raise ValueError('record uncertainty does not establish one negative and one positive energy')
    weights = np.linalg.solve(np.array([[1., 1.], roots]), values[:2])
    predicted = np.array([weights@(roots**k) for k in SAMPLE_NUMBERS])
    held_out = float(np.max(np.abs(predicted[4:]-values[4:])))
    if held_out > max(1e-10, 64*contrast_error):
        raise ValueError('held-out records reject the two-pole recurrence')
    zero_gap_roots = np.array([1.+0j, roots[1]])
    design = np.array([zero_gap_roots**k for k in range(4)])
    zero_weights = np.linalg.lstsq(design, values[:4], rcond=None)[0]
    zero_predictions = np.array([zero_weights@(zero_gap_roots**k) for k in SAMPLE_NUMBERS])
    return {'energies': [float(e) for e in energies],
        'charged_gap': float(-energies[0]), 'root_moduli': [float(abs(z)) for z in roots],
        'record_uncertainty_excludes_zero_charged_gap': True,
        'residues': [[float(z.real), float(z.imag)] for z in weights],
        'training_sample_indices': [0, 1, 2, 3], 'held_out_sample_indices': [4, 5],
        'held_out_maximum_complex_residual': held_out,
        'zero_gap_control_refits_residues_on_training_only': True,
        'zero_gap_control_held_out_maximum_residual': float(np.max(np.abs(zero_predictions[4:]-values[4:]))),
        'one_sample_complex_error_admission': float(contrast_error),
        'record_only_energy_error_bound': float(error),
        'energies_interval': [[float(e-error), float(e+error)] for e in energies],
        'sampling_interval': float(sampling_interval),
        'energy_window': [-math.pi/sampling_interval, math.pi/sampling_interval],
        'no_source_J_mass_or_target_pole_supplied_to_inference': True}


def continuum_bound(mass, exchange, step):
    # ||[A,[A,B]]||=4m²J and ||[B,[A,B]]||=4mJ² in the nontrivial sector.
    local = step**3*(mass*mass*exchange/6+mass*exchange*exchange/3)
    if local >= 1 or step*(2*exchange+mass) >= 1:
        raise ValueError('Strang eigenphase comparison lies outside its unaliased domain')
    return {'step': str(step), 'single_cycle_operator_error_bound': str(local),
        'nested_commutator_operator_norms_exact': [str(4*mass*mass*exchange), str(4*mass*exchange*exchange)],
        'continuous_energy_error_bound': str(F(11, 7)*local/step),
        'bound': '2 asin(epsilon/2)/delta <= (11/7)*epsilon/delta',
        'Strang_bound_primary_source': 'https://arxiv.org/html/1912.08854v3#S5.SS1',
        'Strang_bound_locator': 'Childs et al., Eq.150; A=m(n0-n1), B=J(I-S01)',
        'independent_nyquist_operator_norm_bound': str(2*exchange+mass),
        'clock_parameter': 'supplied common model coordinate; no physical seconds'}


def _schema(arrays):
    shapes = {'encoder': (512, 512), 'frame': (8, 8), 'inverse_frame': (8, 8),
        'cycles': (3, 3, 8, 8), 'detectors': (2, 2, 8, 8),
        'prefix_basis': (512, 6), 'prepared_input': (512,), 'compression_residual': (1,)}
    slabels, hlabels = _labels()
    shapes.update(prefix_decoded=(len(slabels), 512, 7),
                  prefix_coefficients=(len(slabels), 6, 7), prefix_halts=(len(hlabels), 512, 7),
                  success_labels=slabels.shape, halt_labels=hlabels.shape)
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete charged spectroscopy array menu differs')
    for key, shape in shapes.items():
        a = arrays[key]
        kind = np.dtype('int64') if key.endswith('_labels') else np.dtype(float) if key == 'compression_residual' else np.dtype(complex)
        if (type(a) is not np.ndarray or a.shape != shape or a.dtype != kind
                or not np.isfinite(a).all()):
            raise ValueError('original finite charged spectroscopy array schema differs: '+key)
    if (not np.array_equal(arrays['success_labels'], slabels)
            or not np.array_equal(arrays['halt_labels'], hlabels)):
        raise ValueError('complete raw C/E record menu differs')
    residual = float(np.linalg.norm(arrays['prefix_decoded']-
        np.einsum('ar,nrc->nac', arrays['prefix_basis'], arrays['prefix_coefficients'])))
    if (residual > COMPRESSION_ADMISSION
            or abs(residual-float(arrays['compression_residual'][0])) > 128*np.finfo(float).eps):
        raise ValueError('prefix discarded norm differs from its own binary arrays')
    expected_seed = np.zeros(512, complex)
    expected_seed[:8] = preparation_target()
    if np.linalg.norm(arrays['prepared_input']-expected_seed) > SOURCE_ADMISSION:
        raise ValueError('fixed source preparation differs from its independent target')
    j, masses = spectrum.native_parameters()
    for i, m in enumerate(masses):
        for k, dt in enumerate(STEPS):
            if np.linalg.norm(arrays['cycles'][i, k]-ideal_cycle(m, j, dt), 2) > 1e-10:
                raise ValueError('complete local source cycle differs on unprepared inputs')
    for i, axis in enumerate(('X', 'Y')):
        for k, target in enumerate(ideal_detectors(axis)):
            if np.linalg.norm(arrays['detectors'][i, k]-target, 2) > 1e-10:
                raise ValueError('complete local source detector differs')
    if (np.linalg.norm(arrays['frame']-frame_target(), 2) > 1e-10
            or np.linalg.norm(arrays['inverse_frame']-frame_target().T, 2) > 1e-10):
        raise ValueError('source occupation-frame operator differs')


def _actual_output_basis(arrays, family, refinement, number, axis, outcome):
    exponent = int(SPACING*number/STEPS[refinement])
    dynamic = _compensated_power(arrays['cycles'][family, refinement], exponent)
    local = arrays['inverse_frame']@arrays['detectors'][axis, outcome]@dynamic@arrays['frame']
    return arrays['encoder']@_local(local, arrays['prefix_basis'])


def _ideal_map(mass, exchange, step, number, axis, outcome):
    dynamic = _compensated_power(ideal_cycle(mass, exchange, step), int(SPACING*number/step))
    frame = frame_target()
    return (frame.T@ideal_detectors(axis)[outcome]@dynamic@frame)[:6, :6]


def _rounded_complex_product(left, right):
    """Round four real products before addition, forbidding fused multiply-add.

    CPython complex multiplication may contract a multiply/subtract on ARM
    while x86 performs two rounded operations. Boxing each intermediate as a
    Python float fixes that distinction before long cycle powers amplify it.
    """
    ar, ai, br, bi = float(left.real), float(left.imag), float(right.real), float(right.imag)
    rr, ii = ar*br, ai*bi
    ri, ir = ar*bi, ai*br
    return complex(rr-ii, ri+ir)


def _compensated_power(matrix, exponent):
    """Fixed small-matrix products for backend-independent repeated cycles.

    The retained cycle itself is unchanged. In particular, this does not
    project it to a unitary or remove its accumulated source error.
    """
    if type(exponent) is not int or exponent < 0:
        raise ValueError('a nonnegative integer cycle count is required')
    value = np.asarray(matrix, dtype=complex)
    if value.ndim != 2 or value.shape[0] != value.shape[1] or not np.isfinite(value).all():
        raise ValueError('a finite square cycle matrix is required')
    def product(left, right):
        rows, columns = left.tolist(), right.T.tolist()
        result = np.empty_like(left)
        for i, row in enumerate(rows):
            for j, column in enumerate(columns):
                terms = [_rounded_complex_product(a, b) for a, b in zip(row, column, strict=True)]
                result[i, j] = complex(math.fsum(z.real for z in terms),
                                       math.fsum(z.imag for z in terms))
        return result
    result = np.eye(len(value), dtype=complex)
    while exponent:
        if exponent & 1:
            result = product(result, value)
        exponent //= 2
        if exponent:
            value = product(value, value)
    return result


def _compensated_gram(columns):
    """Small column Gram with backend-independent compensated row sums.

    Long sparse native records otherwise make the diagnostic depend on the
    BLAS reduction order at the scale of its near-zero norm defect.
    """
    columns = np.asarray(columns)
    vectors = columns.T.tolist()
    result = np.empty((columns.shape[1], columns.shape[1]), complex)
    for i in range(columns.shape[1]):
        for j in range(i, columns.shape[1]):
            terms = [_rounded_complex_product(a.conjugate(), b)
                     for a, b in zip(vectors[i], vectors[j], strict=True)]
            value = complex(math.fsum(z.real for z in terms), math.fsum(z.imag for z in terms))
            result[i, j], result[j, i] = value, value.conjugate()
    return result


def metrics(arrays):
    _schema(arrays)
    j, masses = spectrum.native_parameters()
    coefficients = arrays['prefix_coefficients']
    ideal_encoder = code.ideal_code_matrix()
    weights = [math.sqrt(float(code.branch_weight('C', SEAMS['C'], cc, cr)
                               *code.branch_weight('E', SEAMS['E'], ec, er)))
               for cc, cr, ec, er in arrays['success_labels'].tolist()]
    halt = arrays['prefix_halts']
    halt_gram = _compensated_gram(halt.reshape(-1, 7))
    coefficient_gram = _compensated_gram(coefficients.reshape(len(coefficients), -1)).reshape(6, 7, 6, 7)
    halt_norm = float(np.linalg.norm(halt[:, :, :6]))**2
    errors, seed_errors, defects, seed_defects, records, contrasts = [], [], [], [], [], []
    for family, mass in enumerate(masses):
        family_rows, family_contrasts = [], []
        for refinement, step in enumerate(STEPS):
            refinement_rows, refinement_contrasts = [], []
            for number in SAMPLE_NUMBERS:
                time_rows, quadratures = [], []
                for axis, name in enumerate(('X', 'Y')):
                    gram = halt_gram.copy()
                    error = halt_norm
                    seed_error = float(np.linalg.norm(halt[:, :, 6]))**2
                    probabilities = []
                    for outcome in (0, 1):
                        out = _actual_output_basis(arrays, family, refinement, number, axis, outcome)
                        target = ideal_encoder@_ideal_map(mass, j, step, number, name, outcome)
                        outputs = np.einsum('ar,nrc->nac', out, coefficients, optimize=True)
                        error += float(np.linalg.norm(outputs[:, :, :6]-np.asarray(weights)[:, None, None]*target))**2
                        seed_error += float(np.linalg.norm(outputs[:, :, 6]-np.asarray(weights)[:, None]*(target@preparation_target()[:6])))**2
                        # Contract the exact retained factors before summing
                        # native records, avoiding a 61,440-row BLAS reduction.
                        gram += np.einsum('rs,rcsd->cd', _compensated_gram(out),
                                          coefficient_gram, optimize=False)
                        probabilities.append(np.sum(np.abs(outputs[:, :, 6])**2, axis=1))
                    errors.append(math.sqrt(error))
                    seed_errors.append(math.sqrt(seed_error))
                    defects.append(float(np.linalg.norm(gram[:6, :6]-np.eye(6), 2)))
                    seed_defects.append(float(abs(gram[6, 6]-1)))
                    raw = np.asarray(probabilities).T
                    quadratures.append(float(np.sum(raw[:, 0]-raw[:, 1])))
                    time_rows.append(raw.tolist())
                refinement_rows.append(time_rows)
                refinement_contrasts.append(complex(*quadratures))
            family_rows.append(refinement_rows); family_contrasts.append(refinement_contrasts)
        records.append(family_rows); contrasts.append(family_contrasts)
    if max(errors+seed_errors) > SOURCE_ADMISSION or max(defects+seed_defects) > SOURCE_ADMISSION:
        raise ValueError('complete charged spectroscopy input/reference source admission failed')
    combined = SOURCE_ADMISSION+COMPRESSION_ADMISSION
    channel = (2+combined)*combined
    inference = [[infer_poles(row, sampling_interval=float(SPACING),
                             contrast_error=math.sqrt(2)*channel)
                  for row in family] for family in contrasts]
    return {'all_six_unknown_input_columns_and_untouched_references': True,
        'terminal_records_per_protocol': len(halt)+2*len(coefficients),
        'successful_prefix_records': len(coefficients), 'halt_records': len(halt),
        'maximum_complete_dilation_error_diagnostic': max(errors),
        'maximum_actual_seed_output_target_error_diagnostic': max(seed_errors),
        'maximum_complete_input_Gram_defect_diagnostic': max(defects),
        'maximum_prepared_source_probability_defect_diagnostic': max(seed_defects),
        'prepared_seed_retained_halt_probability_diagnostic': float(halt_gram[6, 6].real),
        'prefix_discarded_norm_diagnostic': float(arrays['compression_residual'][0]),
        'complex_contrasts': [[[[float(z.real), float(z.imag)] for z in row] for row in family] for family in contrasts],
        'complete_success_record_probabilities': records,
        'record_only_pole_inference': inference,
        'all_settings_retain_unnormalized_halts_without_conditional_rescaling': True,
        'probabilities_are_source_instrument_laws_not_experimental_sample_counts': True,
        'fixed_history_channel_error_allowance': channel,
        'numerical_admissions_are_outward_certificates_or_physical_hardware_errors': False}


def charge_certificate():
    number = s.diag(*[x.bit_count() for x in FOCK_STATES])
    j, masses = spectrum.native_parameters()
    for m in masses:
        h = s.Matrix(hamiltonian(m, j)).extract(FOCK_STATES, FOCK_STATES)
        if h*number != number*h:
            raise ArithmeticError('dimer action violates declared occupation charge')
    for axis in ('X', 'Y'):
        q = s.Matrix(transfer_operator(axis)).extract(FOCK_STATES, FOCK_STATES)
        if q*number != number*q:
            raise ArithmeticError('even transfer detector violates occupation charge')
    return {'occupation_number_eigenvalues': [x.bit_count() for x in FOCK_STATES],
        'action_preparation_and_detector_preserve_total_occupation': True,
        'prepared_total_occupation': 2,
        'occupation_charge_identified_with_SM_electric_charge': False,
        'encoded_occupation_is_distinct_from_the_native_total_port_load_two': True,
        'odd_CAR_operators_closed_on_this_six_level_carrier': False,
        'reference_is_same_encoded_fermion_mode_with_zero_declared_onsite_action': True,
        'charged_pole_reference_energy': 'empty and full dimer energies are zero for the declared ordinary-exchange action',
        'complete_three_mode_Fock_space': False,
        'retained_sector': 'N=1 and N=2; vacuum and completely filled three-mode sectors are outside this code'}


def sm_spectator_certificate():
    raw = strict_json.load(ROOT/'model/observer_particles/results/representation_identification_operators.json')
    subsets = [tuple(raw['basis_subsets'][i]) for i in raw['even_indices']]
    where = {x: i for i, x in enumerate(subsets)}
    matrices = []
    for item in raw['spin_lifted_one_body_generators']:
        one = s.SparseMatrix(*item['shape'], {(i, j): s.sympify(v) for i, j, v in item['entries']})
        lifted = s.zeros(15)
        for col, occupied in enumerate(subsets):
            for old in occupied:
                remaining = tuple(x for x in occupied if x != old)
                for new in range(5):
                    if new not in remaining and one[new, old] != 0:
                        sign = (-1)**(occupied.index(old)+sum(x < new for x in remaining))
                        lifted[where[tuple(sorted((*remaining, new)))], col] += sign*one[new, old]
        matrices.append(s.SparseMatrix(lifted))
    q = s.Matrix(transfer_operator('X')).extract(FOCK_STATES, FOCK_STATES)
    for g in matrices:
        # Exact component identity proves the rectangular source-Kraus case too.
        left = s.kronecker_product(q, s.eye(15))*s.kronecker_product(s.eye(6), g)
        right = s.kronecker_product(s.eye(6), g)*s.kronecker_product(q, s.eye(15))
        if left != right:
            raise ArithmeticError('supplied SM spectator intertwining failed')
    return {'source_representation': 'existing nonvacuum even exterior M15 identification',
        'actual_imported_source_generators': len(matrices), 'spectator_dimension': 15,
        'intertwining_identity': '(K_r tensor I15)(I6 tensor rho(g))=(I512 tensor rho(g))(K_r tensor I15)',
        'all_retained_Kraus_maps_allow_this_algebraic_extension': True,
        'ninety_state_SM_source_execution_performed': False,
        'spectator_extension_is_second_quantization_of_fifteen_internal_species': False,
        'SM_charge_Higgs_Yukawa_or_species_attachment_derived_by_this_extension': False}


def resource_counts(pointer):
    j, masses = spectrum.native_parameters()
    repair = {law: code.repair_resources(law, seam)['native_pulses']+
        max(code.recovery_resources(law, seam, coin, raw)['native_pulses']
            for coin in (0, 1) for raw in circuit.record_values(law, coin))
        for law, seam in SEAMS.items()}
    prepare = compiler.resources(preparation_word())['native_pulses']
    frame = 2*compiler.resources(frame_word())['native_pulses']
    encode = 2*code.encoding_resources()['native_pulses']+code.decoding_resources()['native_pulses']
    cycling = max(int(SPACING*max(SAMPLE_NUMBERS)/dt)*compiler.resources(cycle_word(m, j, dt))['native_pulses']
                  for m in masses for dt in STEPS)
    total = prepare+frame+encode+sum(repair.values())+pointer['native_pulses']+cycling
    return {'one_successful_protocol_native_pulses_upper': total,
        'preparation_native_pulses': prepare, 'occupation_frame_native_pulses': frame,
        'encoding_decoding_native_pulses': encode, 'C_E_repair_recovery_native_pulses_upper': repair,
        'evolution_native_pulses_upper': cycling, 'detector_native_pulses': pointer['native_pulses'],
        'required_proper_load_export_service_events_upper': 2*total,
        'physical_load_export_protection_and_memory_services_executed': False,
        'source_local_matrix_reuse_is_charged_for_every_pulse_invocation': True,
        'distinct_source_protocols': 3*3*6*2,
        'finite_shot_campaign_executed': False,
        'additional_physical_service_error_contract': 'add required_service_events*eta for per-service diamond error eta, plus state/readout/memory errors'}


def execute():
    arrays = _prefix()
    operators, pointer = _source_operators()
    arrays.update(operators)
    report = {'schema': 'oph.native_charged_spectroscopy.v1',
        'classification': 'native recorded carrier to finite interacting charged poles; no measured mass inputs',
        'protocol': {'frame_logical_to_Fock': list(FOCK_STATES),
            'native_seams': {k: list(v) for k, v in SEAMS.items()},
            'order': ['fixed source preparation', 'native matching-code encoding',
                'complete C repair and record-selected recovery', 'complete E repair and record-selected recovery',
                'full native decoder without syndrome selection', 'fixed occupation frame',
                'onsite/exchange Strang source words', 'fixed even transfer X or Y pointer read', 'inverse occupation frame', 'source reencoding'],
            'sampling_interval': str(SPACING), 'sample_numbers': list(SAMPLE_NUMBERS),
            'source_steps': list(map(str, STEPS)),
            'preparation_and_detector_do_not_depend_on_J_m_or_the_pole': True,
            'evolution_angles_use_only_native_J_m_and_declared_model_time': True,
            'same_unknown_source_pole_not_put_into_a_phase_gate': True,
            'independent_spectral_eigenvectors_not_used_to_prepare_the_source': True},
        'admissions': {'complete_source_dilation': SOURCE_ADMISSION,
            'prefix_compression': COMPRESSION_ADMISSION, 'fresh_array_replay': REPLAY_ADMISSION},
        'execution': metrics(arrays), 'charge_interface': charge_certificate(),
        'SM_spectator_interface': sm_spectator_certificate(),
        'continuum': [[continuum_bound(m, spectrum.native_parameters()[0], dt) for dt in STEPS]
                      for m in spectrum.native_parameters()[1]],
        'independent_charged_pole_certificates': [spectrum.dimer_certificate(m, spectrum.native_parameters()[0])
                                                 for m in spectrum.native_parameters()[1]],
        'resources': resource_counts(pointer),
        'scope': {'native_C_E_raw_records_and_512_rows_executed': True,
            'encoded_action_reference_and_controls_remain_declared': True,
            'record_law_identifies_finite_dimer_poles_in_the_supplied_clock': True,
            'measured_particle_spectrum_or_physical_SM_mass_attachment_derived': False,
            'interacting_three_dimensional_continuum_limit_derived': False}}
    return report, arrays


def source_bindings():
    paths = set(spectrum.source_bindings()) | {
        'model/common_source/native_charged_spectroscopy.py',
        'model/common_source/test_native_charged_spectroscopy.py',
        'model/common_source/mass_response/history.py', 'tools/strict_json.py',
        'model/observer_particles/representation_identification.py',
        'model/observer_particles/results/representation_identification_operators.json'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(left, right, path='receipt'):
    if type(left) is not type(right):
        raise ValueError(path+' original JSON type differs')
    if type(left) is dict:
        if left.keys() != right.keys():
            raise ValueError(path+' complete keys differ')
        for key in left:
            _compare(left[key], right[key], path+'.'+key)
    elif type(left) is list:
        if len(left) != len(right):
            raise ValueError(path+' complete length differs')
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            _compare(a, b, path+f'[{i}]')
    elif type(left) is float:
        exact = path.startswith('receipt.admissions.') or path.endswith(('one_sample_complex_error_admission', 'fixed_history_channel_error_allowance'))
        if not math.isfinite(right) or (left != right if exact else not math.isclose(left, right, abs_tol=2e-10, rel_tol=2e-8)):
            raise ValueError(path+' numerical value differs')
    elif left != right:
        raise ValueError(path+' exact value differs')


def compare_execution(report, arrays, fresh, outputs):
    if type(arrays) is not dict or type(outputs) is not dict or set(arrays) != set(outputs):
        raise ValueError('retained and fresh complete native array menus differ')
    metrics(outputs)
    own = metrics(arrays)
    for key, value in own.items():
        if not key.endswith('_diagnostic'):
            continue
        old = report['execution'].get(key)
        tolerance = 128*np.finfo(float).eps
        if type(old) is not float or not math.isclose(value, old, abs_tol=tolerance, rel_tol=2e-10):
            raise ValueError('retained diagnostic differs from its own binary arrays: '
                +key+f'; actual={value!r}; expected={old!r}; '
                +f'delta={abs(value-old) if type(old) is float else None!r}; atol={tolerance!r}')
    # Factor bases can differ by a unitary gauge. Compare their represented maps.
    for key in arrays:
        if key in ('prefix_basis', 'prefix_coefficients', 'compression_residual'):
            continue
        if arrays[key].dtype.kind in 'iu':
            if not np.array_equal(arrays[key], outputs[key]):
                raise ValueError('exact retained record menu differs')
        elif float(np.linalg.norm(arrays[key]-outputs[key])) > REPLAY_ADMISSION:
            raise ValueError('complete native source replay differs: '+key)
    original = np.einsum('ar,nrc->nac', arrays['prefix_basis'], arrays['prefix_coefficients'])
    other = np.einsum('ar,nrc->nac', outputs['prefix_basis'], outputs['prefix_coefficients'])
    if float(np.linalg.norm(original-other)) > REPLAY_ADMISSION:
        raise ValueError('represented complete prefix differs')
    normalized = deepcopy(report)
    for key in own:
        if key.endswith('_diagnostic'):
            normalized['execution'][key] = fresh['execution'][key]
    _compare(fresh, normalized)


def check(path):
    path = Path(path)
    report = strict_json.load(path)
    pins = source_bindings()
    if type(report) is not dict or report.get('source_bindings') != pins:
        raise ValueError('native charged spectroscopy source bindings differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('native charged spectroscopy artifact custody differs')
    with np.load(artifact, allow_pickle=False) as saved:
        if len(saved.files) != len(set(saved.files)):
            raise ValueError('duplicate binary artifact keys')
        arrays = {key: saved[key] for key in saved.files}
    fresh, outputs = execute()
    compare_execution({k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')}, arrays, fresh, outputs)
    if source_bindings() != pins:
        raise ValueError('native charged spectroscopy source changed during replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check)
        print('Native charged-pole spectroscopy and complete source records: PASS')
        return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists():
        parser.error('fresh receipt and artifact paths required')
    pins = source_bindings()
    report, arrays = execute()
    if source_bindings() != pins:
        raise ValueError('native charged spectroscopy source changed during build')
    report['source_bindings'] = pins
    with artifact.open('xb') as out:
        np.savez_compressed(out, **arrays)
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    with args.output.open('x') as out:
        json.dump(report, out, indent=2, sort_keys=True, allow_nan=False)
        out.write('\n')


if __name__ == '__main__':
    main()
