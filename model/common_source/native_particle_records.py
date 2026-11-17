"""Joint unbroken-charge and running-gap records on one native matter register.

The specified twelve-port current supplies the chosen spin/exterior matter
interface. One coherent register has four family addresses, a Nambu bit and
five exterior-mask bits. A neutral-Higgs mass action uses the explicitly
calibrated MSbar M_Z parameters. Three families, the fourth padding family,
all exterior complements and all detector outcomes are retained. Fixed probes
and readout recover charge and running gaps without a species-name table.

The Weyl spin block, calibration, clock and ideal continuous-current service
are supplied architecture. This is not a physical pole-mass prediction,
autonomous particle selection or a derivation of spacetime spin/statistics.
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
import sympy as s

from tools import strict_json
from . import native_matter_attachment as matter
from . import calibrated_matter as calibration
from .chiral_candidate import action
from ..observer.source_scalar import compiler

ROOT = Path(__file__).resolve().parents[2]
TIMES = (F(1, 512), F(1, 4), F(128), F(129))
MASS_WINDOW = (0., 256.)
SOURCE_ADMISSION = 2e-8
RECORD_ERROR = 1e-6
REPLAY_ADMISSION = 3e-9


def _freeze(value):
    value = np.asarray(value, complex)
    return np.frombuffer(value.tobytes(), complex).reshape(value.shape)


@lru_cache(maxsize=1)
def field_architecture():
    """Exact current/neutral-Higgs map; the fixed probe contains no masses."""
    rep = matter.representation()
    selected = rep['selected']
    old = action._representation()
    subsets = [tuple(i for i in range(5) if mask >> i & 1) for mask in selected]
    P = s.zeros(15)
    for i, subset in enumerate(old[0]): P[i, subsets.index(subset)] = 1
    for p, g in enumerate(rep['matter_currents']):
        if not matter.zero(old[2][rep['native']['native_to_fixture_port'][p]]*P-P*g):
            raise ArithmeticError('joint field register has no exact same-current map')
    # This is the same fixed internal singlet rephasing as the calibrated
    # weak-to-mass frame, before its separately supplied CKM family rotation.
    signs = s.ones(15, 1)
    labels = action.ordering()['all_left_components']
    for name, tensors in old[4].items():
        for (i, j), value in tensors['components'][1].todok().items():
            if (labels[i]['internal_type'] == calibration.LEFT_TYPE[name]
                    and labels[j]['internal_type'] == calibration.RIGHT_TYPE[name]):
                signs[j] = value
    R = s.diag(*signs)
    fixed = {name: matter.clean(P.T*R.T*data['components'][1]*R*P)
             for name, data in old[4].items()}
    if any(v not in (-1, 0, 1) for matrix in fixed.values() for v in matrix):
        raise ArithmeticError('the fixed internal Higgs frame is not a signed pairing')
    kernels = []
    for family in range(3):
        kernels.append(matter.clean(sum((matrix*calibration.mass_coordinates(name)[family]
            for name, matrix in fixed.items()), s.zeros(15))))
    # Exact bridge into the existing complete CKM-calibrated field action.
    W = calibration.mass_basis_change()
    P45 = s.kronecker_product(P, s.eye(3))
    full = P45.T*W.T*calibration.mass_kernel()*W*P45
    expected = s.zeros(45)
    for family, matrix in enumerate(kernels):
        for (i, j), value in matrix.todok().items(): expected[3*i+family, 3*j+family] = value
    if not calibration._zero(full-expected):
        raise ArithmeticError('joint register is not the calibrated action in its declared mass frame')
    unit_kernel = sum(fixed.values(), s.zeros(15))
    K64 = s.zeros(64)
    for (i, j), value in unit_kernel.todok().items():
        K64[selected[i], 32+selected[j]] = value
        K64[32+selected[j], selected[i]] = value
    pairs = [(i, j, int(K64[i, j])) for i in range(64) for j in range(i+1, 64) if K64[i, j] != 0]
    used = [i for pair in pairs for i in pair[:2]]
    if len(pairs) != 14 or len(set(used)) != 28:
        raise ArithmeticError('neutral Higgs support is not fourteen disjoint Nambu pairs')
    images = tuple(used+[i for i in range(64) if i not in used])
    phase_signs = tuple(x for _, _, sign in pairs for x in (1, sign))+(1,)*36
    permutation = np.eye(64)[:, list(images)]
    signed = permutation@np.diag(phase_signs)
    q32 = rep['operators']['T3']+rep['operators']['q']/6
    c32 = rep['operators']['color_C2']
    q15 = q32.extract(list(selected), list(selected))
    c15 = c32.extract(list(selected), list(selected))
    old_q = action.electric_charge(families=False)
    if not matter.zero(old_q*P-P*q15):
        raise ArithmeticError('native and existing field electric charges differ')
    for operator in (old_q, P*c15*P.T):
        on_families = s.kronecker_product(operator, s.eye(3))
        if not matter.zero(on_families*W-W*on_families):
            raise ArithmeticError('supplied family mass frame changes the read unbroken charge')
    if any(not calibration._zero(r.H*r-s.eye(3)) for r in calibration.mixing_rotations()):
        raise ArithmeticError('the supplied weak-to-mass family frame is not unitary')
    q64 = s.diag(q32, -q32.T)
    c64 = s.diag(c32, c32.T)
    masses = []
    targets = []
    for family, kernel in enumerate(kernels):
        H = s.zeros(64)
        for (i, j), value in kernel.todok().items():
            H[selected[i], 32+selected[j]] = value
            H[32+selected[j], selected[i]] = value
        if not matter.zero(q64*H-H*q64) or not matter.zero(c64*H-H*c64):
            raise ArithmeticError('mass action does not preserve the measured unbroken charges')
        values = tuple(s.simplify(H[a, b]*sign) for a, b, sign in pairs)
        if any(x.is_positive is not True or float(x) >= MASS_WINDOW[1] for x in values):
            raise ArithmeticError('positive calibrated action exceeds the declared observer window')
        masses.append(values); targets.append(H)
    return {'P15': P, 'P45': P45, 'weak_to_mass': W, 'fixed_higgs_components': fixed,
        'images': images, 'phase_signs': phase_signs, 'signed_pair_frame': signed,
        'pairs': tuple(pairs), 'masses': tuple(masses), 'targets': tuple(targets),
        'charge': tuple(q64.diagonal()), 'color': tuple(c64.diagonal()),
        'selected': selected}


def _flip_word():
    # Phase-exact six-bit controlled X, with every Walsh phase retained.
    bits, target = 6, 5
    angle = -math.pi/64
    word = [compiler.pulse('Y', math.pi/4, [target]), compiler.pulse('I', angle, [0])]
    for pivot in range(bits):
        count = 1 << pivot
        for i in range(count):
            word.append(compiler.pulse('Z', (-1 if pivot == target else 1)*angle, [pivot]))
            difference = (i ^ (i >> 1)) ^ (((i+1) % count) ^ (((i+1) % count) >> 1))
            if difference: word.extend(compiler.cnot(difference.bit_length()-1, pivot))
    word.append(compiler.pulse('Y', -math.pi/4, [target]))
    return word


def _edges(images):
    seen, pairs = set(), []
    for start in range(64):
        if start in seen: continue
        cycle, current = [], start
        while current not in seen:
            seen.add(current); cycle.append(current); current = images[current]
        pairs.extend((start, other) for other in cycle[1:])
    edges = []
    for a, b in pairs:
        path = [a]
        for bit in range(6):
            if (a ^ b) >> bit & 1: path.append(path[-1] ^ (1 << bit))
        adjacent = list(zip(path[:-1], path[1:], strict=True))
        for x, y in adjacent+adjacent[-2::-1]:
            bit = (x ^ y).bit_length()-1
            edges.append((bit, x & ~(1 << bit)))
    return tuple(edges)


def _permutation_source(images, macro):
    state = np.eye(64, dtype=complex)
    indices = np.arange(64)
    for target, pattern in _edges(images):
        word = [compiler.pulse('X', math.pi/2, [j]) for j in range(6) if pattern >> j & 1]
        state = matter.execute_word(word, 6, state)
        wires = [j for j in range(6) if j != target]+[target]
        chart = sum(((indices >> wire) & 1) << local for local, wire in enumerate(wires))
        state = macro[np.ix_(chart, chart)]@state
        state = matter.execute_word(compiler.inverse(word), 6, state)
    return state


def _frame_source():
    data = field_architecture()
    macro = matter.execute_word(_flip_word(), 6)
    images = data['images']
    forward = _permutation_source(images, macro)
    inverse = _permutation_source(tuple(images.index(i) for i in range(64)), macro)
    phases = [0 if sign == 1 else math.pi for sign in data['phase_signs']]
    signs = matter.execute_word(compiler.compile_diagonal(phases), 6)
    inverse_signs = matter.execute_word(compiler.compile_diagonal([-x for x in phases]), 6)
    return forward@signs, inverse_signs@inverse, macro


def _phases(time):
    phases = np.zeros(256)
    for family, values in enumerate(field_architecture()['masses']):
        for pair, mass in enumerate(values):
            angle = math.remainder(float(mass)*float(time), 2*math.pi)
            phases[64*family+2*pair] = -angle
            phases[64*family+2*pair+1] = angle
    return phases


def pointer_probabilities(states, local_kraus):
    """Exact tensor contraction of all eight pointer bits and quantum rows."""
    effects = np.asarray([(k.conj().T@k).T for k in local_kraus])
    density = np.einsum('ic,jc->ijc', states, states.conj()).reshape((2,)*16+(states.shape[1],))
    density = density.transpose(tuple(j for i in range(8) for j in (i, 8+i))+(16,))
    for _ in range(8):
        density = np.tensordot(effects, density, axes=((1, 2), (0, 1)))
        density = np.moveaxis(density, 0, -2)
    result = density.reshape(256, states.shape[1])
    if np.max(abs(result.imag)) > 1e-10 or np.min(result.real) < -1e-12:
        raise ArithmeticError('complete pointer contraction is not a probability law')
    return result.real


def record_dictionary():
    data = field_architecture(); inverse = np.argsort(data['images'])
    return tuple({'charge': str(data['charge'][r % 64]), 'color_C2': str(data['color'][r % 64]),
        'pair_sign': 1 if inverse[r % 64] % 2 == 0 else -1,
        'selected_exterior': r % 32 in data['selected'], 'padding_family': r // 64 == 3}
        for r in range(256))


def decode_records(probabilities, *, dictionary, times=TIMES,
                   mass_window=MASS_WINDOW, record_error=RECORD_ERROR):
    """No action, mass parameter, family label or species table enters here."""
    if type(probabilities) is not np.ndarray:
        raise ValueError('original unmasked complete record array required')
    values = probabilities
    if values.shape != (4, 2, 256) or values.dtype.kind not in 'if' or not np.isfinite(values).all():
        raise ValueError('four complete two-quadrature256-outcome laws are required')
    if (len(dictionary) != 256 or len(times) != 4
            or not all(float(t) > 0 for t in times)
            or not all(float(a) < float(b) for a, b in zip(times[:-1], times[1:], strict=True))
            or type(record_error) is not float or not math.isfinite(record_error) or not 0 < record_error < .01):
        raise ValueError('complete declared apparatus dictionary and error budget required')
    if np.max(abs(values.sum(axis=2)-1)) > record_error or np.min(values) < -record_error:
        raise ValueError('all outcomes must remain unnormalized and complete')
    scores = np.asarray([row['pair_sign'] for row in dictionary])
    z = values@scores
    phases = z[:, 0]+1j*z[:, 1]
    radius = math.sqrt(2)*record_error
    if np.max(abs(abs(phases)-1)) > radius:
        return {'resolved': False, 'reason': 'records do not resolve one coherent frequency'}
    angular = 2*math.asin(radius)
    low, high = map(float, mass_window)
    if not 0 <= low < high or not 2*high*float(times[0]) < math.pi:
        raise ValueError('a fixed positive initially unaliased clock window is required')
    for measured, time in zip(phases[:3], times[:3], strict=True):
        t = float(time)
        angle = math.atan2(measured.imag, measured.real)
        first = math.ceil((2*t*low-angle-angular)/(2*math.pi))
        last = math.floor((2*t*high-angle+angular)/(2*math.pi))
        if last-first > 4:
            return {'resolved': False, 'reason': 'record interval permits multiple clock windings'}
        candidates = []
        for winding in range(first, last+1):
            a = max(low, (angle+2*math.pi*winding-angular)/(2*t))
            b = min(high, (angle+2*math.pi*winding+angular)/(2*t))
            if a <= b: candidates.append((a, b))
        if len(candidates) != 1:
            return {'resolved': False, 'reason': 'record interval does not select one clock winding'}
        low, high = candidates[0]
    centre = (low+high)/2
    prediction = complex(math.cos(2*centre*float(times[3])), math.sin(2*centre*float(times[3])))
    heldout = float(abs(prediction-phases[3]))
    if heldout > radius+float(times[3])*(high-low):
        return {'resolved': False, 'reason': 'held-out source records reject the recovered gap'}
    groups = sorted({(row['charge'], row['color_C2']) for row in dictionary})
    weights = [float(np.mean(values[:, :, [i for i, row in enumerate(dictionary)
        if (row['charge'], row['color_C2']) == group]].sum(axis=2))) for group in groups]
    chosen = [g for g, weight in zip(groups, weights, strict=True) if weight >= 1-record_error]
    if len(chosen) != 1:
        return {'resolved': False, 'reason': 'records do not select a single unbroken-charge sector'}
    return {'resolved': True, 'electric_Q': chosen[0][0], 'color_C2': chosen[0][1],
        'running_gap_over_common_GeV_reference_interval': [low, high],
        'heldout_complex_record_error': heldout,
        'identified_particle_name_or_family_label': None}


def _evolution(arrays, number):
    frame = np.kron(np.eye(4), arrays['frame'])
    inverse = np.kron(np.eye(4), arrays['inverse_frame'])
    rotate = np.kron(np.eye(4), arrays['rotate'])
    inverse_rotate = np.kron(np.eye(4), arrays['inverse_rotate'])
    return frame@inverse_rotate@arrays['phases'][number]@rotate@inverse


def _ideal_evolution(time):
    data = field_architecture(); P = data['signed_pair_frame']
    result = np.eye(256, dtype=complex)
    for family, values in enumerate(data['masses']):
        value = np.eye(64, dtype=complex)
        for pair, mass in enumerate(values):
            phase = float(mass)*float(time)
            c, d = math.cos(phase), -1j*math.sin(phase)
            value[2*pair:2*pair+2, 2*pair:2*pair+2] = [[c, d], [d, c]]
        result[64*family:64*(family+1), 64*family:64*(family+1)] = P@value@P.T
    return result


def _source_resources():
    data = field_architecture()
    flip = matter.resources(_flip_word())['native_pulses']
    forward = len(_edges(data['images']))*flip+2*sum(pattern.bit_count() for _, pattern in _edges(data['images']))
    inverse_images = tuple(data['images'].index(i) for i in range(64))
    inverse = len(_edges(inverse_images))*flip+2*sum(pattern.bit_count() for _, pattern in _edges(inverse_images))
    sign = matter.resources(compiler.compile_diagonal([0 if x == 1 else math.pi for x in data['phase_signs']]))['native_pulses']
    forward += sign; inverse += sign
    phase_counts = [matter.resources(compiler.compile_diagonal(_phases(t)))['native_pulses'] for t in TIMES]
    # Initial encoding, evolution decode/reencode, detector decode/reencode,
    # two energy-frame rotations, one detector rotation, eight pointer CNOTs.
    total = 3*forward+2*inverse+max(phase_counts)+3+8*3
    return {'frame_native_pulses': forward, 'inverse_frame_native_pulses': inverse,
        'controlled_family_phase_native_pulses_per_time': phase_counts,
        'pointer_native_pulses': 24, 'one_protocol_native_pulses_upper': total,
        'required_proper_load_export_events_upper': 2*total,
        'logical_work_qubits': 8, 'sequential_pointer_qubits': 1,
        'source_services_executed': 'ideal continuous current and declared proper-code transfer',
        'full_finite64_native_pulse_backend_executed': False,
        'complete_native_C_E_matter_branch_history_executed': False,
        'required_matter_branch_lift': 'sqrt(w) times rho(g_lift) for the chosen lifted control path; the scalar is retained once and the endpoint U6 does not fix spin holonomy',
        'physical_state_clock_memory_and_readout_services_certified': False,
        'physical_error_obligation': 'add N*delta_current+2*N*eta_transfer plus clock/preparation/storage/pointer errors for N source pulses; the continuous-current and proper-transfer diamond errors are separate physical contracts'}


def _schema(arrays):
    shapes = {'frame': (64, 64), 'inverse_frame': (64, 64), 'macro': (64, 64),
        'rotate': (64, 64), 'inverse_rotate': (64, 64), 'detector': (2, 64, 64),
        'phases': (4, 256, 256), 'pointer': (2, 2, 2), 'probe_indices': (42,),
        'probabilities': (42, 4, 2, 256)}
    if type(arrays) is not dict or set(arrays) != set(shapes):
        raise ValueError('complete joint native particle array menu required')
    for key, shape in shapes.items():
        value = arrays[key]
        dtype = np.dtype('int64') if key == 'probe_indices' else np.dtype(float) if key == 'probabilities' else np.dtype(complex)
        if type(value) is not np.ndarray or value.dtype != dtype or value.shape != shape or not np.isfinite(value).all():
            raise ValueError('complete original array required: '+key)
    indices = np.array([64*f+2*p for f in range(3) for p in range(14)], np.int64)
    if not np.array_equal(indices, arrays['probe_indices']):
        raise ValueError('the fixed mass-independent probe menu differs')


def metrics(arrays):
    _schema(arrays)
    data = field_architecture(); target = data['signed_pair_frame']
    frame_error = float(np.linalg.norm(arrays['frame']-target))
    inverse_error = float(np.linalg.norm(arrays['inverse_frame']-target.T))
    rotation = np.array([[1, 1], [-1, 1]], complex)/math.sqrt(2)
    rotate_target = np.kron(np.eye(32), rotation)
    read_target = np.kron(np.eye(32), (np.eye(2)+1j*np.array([[0, 1], [1, 0]]))/math.sqrt(2))
    targets = [target@np.eye(64)@target.T, target@read_target@target.T]
    macro_target = np.eye(64, dtype=complex)
    macro_target[:, [0, 32]] = macro_target[:, [32, 0]]
    errors = [frame_error, inverse_error, float(np.linalg.norm(arrays['macro']-macro_target)),
        float(np.linalg.norm(arrays['rotate']-rotate_target)),
        float(np.linalg.norm(arrays['inverse_rotate']-rotate_target.T)),
        float(np.linalg.norm(arrays['detector']-np.asarray(targets)))]
    pointer_target = np.array([np.diag([1., 0.]), np.diag([0., 1.])], complex)
    pointer_error = float(np.linalg.norm(arrays['pointer']-pointer_target))
    errors.append(pointer_error)
    if max(errors) > SOURCE_ADMISSION:
        raise ValueError('complete fixed source apparatus admission exceeded')
    local_gram = sum((k.conj().T@k for k in arrays['pointer']), np.zeros((2, 2), complex))
    pointer_gram = local_gram
    for _ in range(7): pointer_gram = np.kron(pointer_gram, local_gram)
    prepared = np.kron(np.eye(4), arrays['frame'])
    defects, output_errors, probabilities = [], [], []
    for n, t in enumerate(TIMES):
        U = _evolution(arrays, n); target_U = _ideal_evolution(t)
        errors.append(float(np.linalg.norm(U-target_U)))
        if errors[-1] > SOURCE_ADMISSION:
            raise ValueError('complete unknown-input mass evolution admission exceeded')
        rows = []
        for q in range(2):
            output = np.kron(np.eye(4), arrays['detector'][q])@U@prepared
            ideal = np.kron(np.eye(4), targets[q])@target_U@np.kron(np.eye(4), target)
            output_errors.append(float(np.linalg.norm(output-ideal)))
            defects.append(float(np.linalg.norm(output.conj().T@pointer_gram@output-np.eye(256), 2)))
            rows.append(pointer_probabilities(output[:, arrays['probe_indices']], arrays['pointer']).T)
        probabilities.append(rows)
    calculated = np.asarray(probabilities).transpose(2, 0, 1, 3)
    if np.max(abs(calculated-arrays['probabilities'])) > 2e-12:
        raise ValueError('joint records differ from their own complete source factors')
    if max(errors+output_errors+defects) > SOURCE_ADMISSION:
        raise ValueError('complete unknown-input/reference source admission exceeded')
    # Four full source maps (preparation/evolution/detector/pointer), each
    # within its stated numerical admission; eight local pointer factors.
    dilation_allowance = (1+SOURCE_ADMISSION)**11-1
    channel_allowance = (2+dilation_allowance)*dilation_allowance
    if channel_allowance > RECORD_ERROR:
        raise ArithmeticError('record inverse budget does not cover source admission')
    dictionary = record_dictionary()
    decoded = [decode_records(row, dictionary=dictionary) for row in arrays['probabilities']]
    if not all(row['resolved'] for row in decoded):
        raise ValueError('complete joint native records did not resolve every declared probe')
    for index, row in zip(arrays['probe_indices'], decoded, strict=True):
        family, pair = int(index)//64, int(index)%64//2
        expected_mass = float(data['masses'][family][pair])
        low, high = row['running_gap_over_common_GeV_reference_interval']
        original = data['images'][2*pair]
        if not low <= expected_mass <= high or (row['electric_Q'], row['color_C2']) != (str(data['charge'][original]), str(data['color'][original])):
            raise ValueError('blind joint records disagree with the supplied same-register action')
    return {'complete_register_input_columns': 256, 'terminal_pointer_outcomes': 256,
        'fixed_probes': 42, 'common_times': list(map(str, TIMES)),
        'complete_source_operator_error_diagnostic': max(errors),
        'complete_protocol_operator_error_diagnostic': max(output_errors),
        'complete_unknown_input_Gram_defect_diagnostic': max(defects),
        'fixed_dilation_allowance': dilation_allowance, 'fixed_channel_allowance': channel_allowance,
        'record_error_admission': RECORD_ERROR, 'decoded_unordered_charge_gap_records': decoded,
        'record_laws_are_analytic_not_new_experimental_counts': True,
        'family_or_species_targets_supplied_to_inverse': False,
        'all_family_padding_exterior_and_nambu_records_retained': True,
        'any_output_postselected_or_probabilities_renormalized': False}


def _calibration_packet_name():
    return calibration.PACKET.relative_to(ROOT).as_posix()


def execute():
    frame, inverse, macro = _frame_source()
    rotate = matter.execute_word([compiler.pulse('Y', math.pi/4, [0])], 6)
    inverse_rotate = matter.execute_word([compiler.pulse('Y', -math.pi/4, [0])], 6)
    read = matter.execute_word([compiler.pulse('X', math.pi/4, [0])], 6)
    pointer_full = matter.execute_word(compiler.cnot(0, 1), 2)
    arrays = {'frame': frame, 'inverse_frame': inverse, 'macro': macro,
        'rotate': rotate, 'inverse_rotate': inverse_rotate,
        'detector': np.asarray([frame@inverse, frame@read@inverse]),
        'phases': np.asarray([matter.execute_word(compiler.compile_diagonal(_phases(t)), 8) for t in TIMES]),
        'pointer': np.asarray([pointer_full[:2, :2], pointer_full[2:, :2]]),
        'probe_indices': np.array([64*f+2*p for f in range(3) for p in range(14)], np.int64)}
    prepared = np.kron(np.eye(4), frame)[:, arrays['probe_indices']]
    laws = []
    for n in range(4):
        evolved = _evolution(arrays, n)@prepared
        laws.append([pointer_probabilities(np.kron(np.eye(4), detector)@evolved, arrays['pointer']).T
            for detector in arrays['detector']])
    arrays['probabilities'] = np.asarray(laws).transpose(2, 0, 1, 3)
    report = {'schema': 'oph.native_particle_records.v1',
        'classification': 'joint native same-register unbroken-charge and calibrated running-gap component',
        'calibration': {'parameter_packet': _calibration_packet_name(),
            'convention': 'MSbar running action at M_Z in the supplied neutral Higgs background',
            'masses_predicted_without_empirical_parameters': False,
            'experimental_or_interacting_pole_masses_recovered': False},
        'architecture': {'native_carrier_ports': 12, 'family_register_levels': 4,
            'calibrated_families': 3, 'retained_padding_family': 3,
            'exterior_mask_levels': 32, 'nambu_levels': 2,
            'nambu_current': 'diag(G,conjugate(G)) for anti-Hermitian matter G',
            'nambu_charge': 'diag(Q,-transpose(Q)) for Hermitian unbroken Q',
            'basis_order': '64*family+32*particle_hole+exterior_mask',
            'selected_spin_block': 'particle spin up and conjugate hole spin down; common supplied S_z=+1/2',
            'complementary_spin_block': 'particle spin down and conjugate hole spin up; mass off-diagonal sign reverses',
            'physical_spin_factor_supplied_not_weak_SU2': True,
            'spin_or_statistics_derived': False,
            'field_intertwiner': 'P maps native exterior order into existing action; I_spin tensor(P tensor I3), with supplied CKM weak-to-mass frame',
            'mass_frame_depends_on_calibrated_masses': False,
            'mass_parameters_in_evolution_word': True,
            'clock': 'one supplied inverse-GeV action parameter for all modes',
            'gap_convention': 'positive rest-frame Nambu energy m is half the observed2m Bohr frequency of the declared symmetric mass block',
            'pointer_dictionary': list(record_dictionary())},
        'execution': metrics(arrays), 'resources': _source_resources(),
        'scope': {'same_current_neutral_higgs_field_intertwiner_checked_exactly': True,
            'mass_frame_preserves_read_unbroken_currents_exactly': True,
            'coherent_family_register_and_complement_retained': True,
            'CKM_propagation_or_full_interacting_field_simulated_here': False,
            'autonomous_selection_or_measured_particle_spectrum_derived': False,
            'numerical_source_admissions_are_outward_or_hardware_certificates': False}}
    return report, arrays


def source_bindings():
    paths = set(matter.source_bindings())|set(calibration.source_bindings())|{
        'model/common_source/native_particle_records.py', 'model/common_source/test_native_particle_records.py',
        'tools/strict_json.py'}
    return {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in sorted(paths)}


def _compare(expected, actual, path='receipt'):
    if type(expected) is not type(actual): raise ValueError(path+' original type differs')
    if type(expected) is dict:
        if expected.keys() != actual.keys(): raise ValueError(path+' keys differ')
        for key in expected: _compare(expected[key], actual[key], path+'.'+key)
    elif type(expected) is list:
        if len(expected) != len(actual): raise ValueError(path+' length differs')
        for i, (a, b) in enumerate(zip(expected, actual, strict=True)): _compare(a, b, path+f'[{i}]')
    elif type(expected) is float:
        exact = path.endswith(('allowance', 'admission'))
        if not math.isfinite(actual) or (expected != actual if exact else not math.isclose(expected, actual, rel_tol=1e-8, abs_tol=1e-10)):
            raise ValueError(path+f' numerical value differs; actual={actual!r}; '
                +f'expected={expected!r}; delta={abs(actual-expected)!r}')
    elif expected != actual: raise ValueError(path+' exact value differs')


def compare_execution(report, arrays, fresh, outputs):
    """Bind stored diagnostics to stored arrays, then compare fresh execution.

    A fresh backend may execute a slightly different floating source matrix
    within the unchanged source and replay admissions. Its error norm belongs
    to that fresh matrix, rather than to the retained matrix's diagnostic.
    """
    own = metrics(arrays)
    for key, value in own.items():
        if key.endswith('_diagnostic'):
            retained = report['execution'].get(key)
            tolerance = 128*np.finfo(float).eps
            if type(retained) is not float or not math.isclose(value, retained, rel_tol=2e-10, abs_tol=tolerance):
                raise ValueError('diagnostic differs from its own binary arrays: '+key
                    +f'; actual={value!r}; expected={retained!r}; '
                    +f'delta={abs(value-retained) if type(retained) is float else None!r}; atol={tolerance!r}')
    _compare(own, report['execution'])
    _schema(outputs)
    if set(arrays) != set(outputs): raise ValueError('retained/fresh array menu differs')
    for key in arrays:
        if np.linalg.norm(arrays[key]-outputs[key]) > REPLAY_ADMISSION:
            raise ValueError('complete source replay differs: '+key)
    normalized = deepcopy(report)
    for key in own:
        if key.endswith('_diagnostic'):
            normalized['execution'][key] = fresh['execution'][key]
    _compare(fresh, normalized)


def check(path):
    path = Path(path); pins = source_bindings(); report = strict_json.load(path)
    if type(report) is not dict or report.get('source_bindings') != pins: raise ValueError('joint-record source pins differ')
    artifact = path.with_suffix('.npz')
    if report.get('artifact') != {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}:
        raise ValueError('joint-record binary custody differs')
    with np.load(artifact, allow_pickle=False) as saved:
        if len(saved.files) != len(set(saved.files)): raise ValueError('duplicate binary array names')
        arrays = {key: saved[key] for key in saved.files}
    fresh, outputs = execute()
    compare_execution({k: v for k, v in report.items() if k not in ('source_bindings', 'artifact')}, arrays, fresh, outputs)
    if source_bindings() != pins: raise ValueError('sources changed during joint-record replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--output', type=Path); mode.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        check(args.check); print('Joint native charge and calibrated running-gap records: PASS'); return
    artifact = args.output.with_suffix('.npz')
    if args.output.exists() or artifact.exists(): parser.error('fresh receipt and binary paths required')
    pins = source_bindings(); report, arrays = execute()
    if source_bindings() != pins: raise ValueError('sources changed during joint-record execution')
    with artifact.open('xb') as out: np.savez_compressed(out, **arrays)
    report['source_bindings'] = pins
    report['artifact'] = {'filename': artifact.name, 'sha256': hashlib.sha256(artifact.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
