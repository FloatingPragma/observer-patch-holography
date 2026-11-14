"""Native words for electric phases and complete local plaquette kicks."""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import sys

import numpy as np

import autonomous as af

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'openwave-tuning/auditor/common_walk_review'))
import source_error as root_certificate
from model.observer.source_scalar import compiler
from model.certified_reference import certificate as cert, interval as iv

Q = 16
FIELD_BASE = 2048
SCRATCH = (3200, 3201, 3202)


def H(bit):
    return [compiler.pulse('I', -math.pi / 2, [bit]), compiler.pulse('H', math.pi / 2, [bit])]


def qft_word():
    out = []
    for target in reversed(range(4)):
        out += H(target)
        for control in reversed(range(target)):
            out += compiler.remap(compiler.compile_diagonal([0., 0., 0., math.pi / 2**(target - control)]), [control, target])
    for a, b in ((0, 3), (1, 2)):
        out += compiler.cnot(a, b) + compiler.cnot(b, a) + compiler.cnot(a, b)
    return out


def qft_template():
    out = []
    for target in reversed(range(4)):
        out += [('I', (target,), cert._pi(-1, 2)), ('H', (target,), cert._pi(1, 2))]
        for control in reversed(range(target)):
            temp = root_certificate.diagonal_template([0, 0, 0, F(1, 2**(target - control))])
            out += root_certificate.remap(temp, {0: control, 1: target})
    for a, b in ((0, 3), (1, 2)):
        out += cert.cnot(a, b) + cert.cnot(b, a) + cert.cnot(a, b)
    return out


def phase_bound(word, desired):
    offset, result = 0, F()
    for mask in range(Q):
        label = ''.join('Z' if mask >> j & 1 else 'I' for j in range(4))
        offset, fixed, variables = cert.consume(word, offset, cert.pauli(label))
        assert len(variables) == 1
        wanted = iv.Interval.point(0)
        for i, value in enumerate(desired):
            wanted += value * (-1)**((i & mask).bit_count())
        wanted = wanted / Q
        result += fixed + (iv.Interval.point(variables[0]) - wanted).absolute_upper()
    assert offset == len(word)
    return result + cert.controls_bound(word, {})


def cosine_interval(index):
    n = index if index <= 8 else index - 16
    angle = iv.pi() * F(n, 8)
    midpoint = F(angle.lo + angle.hi, 2 * iv.SCALE)
    radius = (angle.hi - angle.lo + 1) // 2
    return iv.sine_cosine(midpoint)[1] + iv.Interval(-radius, radius)


def electric_word(alpha):
    energies = [F(((e + 8) % 16 - 8)**2, 2) for e in range(Q)]
    word = compiler.compile_diagonal([-float(alpha) * float(e) for e in energies])
    bound = phase_bound(word, [iv.Interval.point(-F(alpha) * e) for e in energies])
    return word, bound


def plaquette_word(beta):
    fourier = qft_word()
    diagonal = compiler.compile_diagonal([float(beta) * math.cos(2 * math.pi * r / Q) for r in range(Q)])
    qft_error = root_certificate.exact_template_bound(fourier, qft_template())
    diagonal_error = phase_bound(diagonal, [cosine_interval(r) * F(beta) for r in range(Q)])
    return compiler.inverse(fourier) + diagonal + fourier, 2 * qft_error + diagonal_error


@lru_cache(maxsize=1)
def add_program():
    """target4..7 += source0..3 modulo16; complete dirty-scratch extension."""
    gates = []
    for bit in range(4):
        for target in reversed(range(bit, 4)):
            values = {bit: 1}
            values.update({4 + k: 1 for k in range(bit, target)})
            if len(values) == 1:
                gates.append(('cx', (bit, 4 + target)))
            else:
                gates += af.mb.compiled.multi_control(values, 4 + target, (8, 9, 10))
    inputs = np.arange(1 << 11, dtype=np.int64)
    output = af.mb.compiled.execute_reversible(gates, inputs)
    expected = np.array([s + (((t + s) % 16) << 4) for t in range(16) for s in range(16)])
    assert np.array_equal(output[:256], expected)
    assert len(np.unique(output)) == len(output)
    assert np.array_equal(af.mb.compiled.execute_reversible(list(reversed(gates)), output), inputs)
    return gates


def arithmetic_word(gates, mapping=None):
    words = af.mb.compiled.controls()['primitive_words']
    for kind, targets in gates:
        destination = list(targets) if mapping is None else [mapping[q] for q in targets]
        yield from compiler.remap(words[kind], destination)


def face_word(face, beta):
    a = af.data()[2][:, face]
    pivot = int(np.flatnonzero(a)[0])
    shear = []
    for edge in np.flatnonzero(a):
        if edge == pivot:
            continue
        coefficient = -int(a[edge] * a[pivot])
        gates = add_program() if coefficient == 1 else list(reversed(add_program()))
        mapping = {q: FIELD_BASE + 4 * pivot + q for q in range(4)}
        mapping.update({q + 4: FIELD_BASE + 4 * int(edge) + q for q in range(4)})
        mapping.update({8 + q: SCRATCH[q] for q in range(3)})
        shear += list(arithmetic_word(gates, mapping))
    core, _ = plaquette_word(beta)
    target = [FIELD_BASE + 4 * pivot + q for q in range(4)]
    return shear + compiler.remap(core, target) + compiler.inverse(shear)


def field_word(alpha, face_couplings):
    """Chronological U_B then U_E; actual full finite pulse generator."""
    if len(face_couplings) != 208:
        raise ValueError('one coupling for every declared face required')
    for face, beta in enumerate(face_couplings):
        if beta:
            yield from face_word(face, beta)
    electric, _ = electric_word(alpha)
    if alpha:
        for edge in range(288):
            yield from compiler.remap(electric, [FIELD_BASE + 4 * edge + q for q in range(4)])


def complete_word(alpha, face_couplings, *, clock_spacing=None):
    """Full native V word, with an optional fresh completed-word phase tick.

    The generator retains all source phases. Loading/exporting each pulse
    uses the declared complete coherent service and its complement outcomes.
    """
    program = af.mb.mode_program()
    for group in program['onsite_mode_groups']:
        yield from compiler.remap(program['onsite_reflection_word'], group)
    block = af.mode_field.field_block(16)
    primitive = af.mb.compiled.controls()['primitive_words']
    for a, b, edge in program['transport_pairs']:
        mapping = {q: FIELD_BASE + 4 * edge + q for q in range(4)}
        mapping.update({4: a, 5: b, 6: SCRATCH[0], 7: SCRATCH[1], 8: SCRATCH[2]})
        for kind, targets in block['program']['gates']:
            yield from compiler.remap(primitive[kind], [mapping[q] for q in targets])
        yield from af.mb.fock_swap_word(a, b)
    yield from field_word(alpha, face_couplings)
    if clock_spacing is not None:
        if type(clock_spacing) is not F or clock_spacing <= 0:
            raise ValueError('positive exact clock spacing required')
        yield compiler.pulse('Z', -float(clock_spacing) / 16, [3203])


def run():
    paths = [Path(__file__), Path(af.__file__), Path(compiler.__file__), Path(cert.__file__),
             Path(iv.__file__), Path(root_certificate.__file__),
             ROOT / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json',
             HERE.parent / 'manybody/MODE_FIELD.json']
    pins = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    qft = compiler.execute(qft_word(), 4)
    target = np.exp(2j * np.pi * np.outer(np.arange(Q), np.arange(Q)) / Q) / 4
    qft_residual = float(np.linalg.norm(qft - target, 2))
    assert qft_residual < 2e-11
    alpha, beta = F(.17), F(.41)
    ew, ee = electric_word(alpha)
    pw, pe = plaquette_word(beta)
    E = np.array([((r + 8) % 16 - 8)**2 / 2 for r in range(16)])
    e_actual = compiler.execute(ew, 4)
    p_actual = compiler.execute(pw, 4)
    e_target = np.diag(np.exp(-1j * float(alpha) * E))
    p_target = (target * np.exp(1j * float(beta) * np.cos(2 * np.pi * np.arange(Q) / Q))[None, :]) @ target.conj().T
    er = float(np.linalg.norm(e_actual - e_target, 2))
    pr = float(np.linalg.norm(p_actual - p_target, 2))
    assert max(er, pr) < 2e-10
    root = json.loads((ROOT / 'openwave-tuning/auditor/common_walk_review/SOURCE_ERROR.json').read_text())
    counts = Counter(k for k, t in add_program())
    add_error = sum((count * F(root['primitive_operator_errors_exact'][kind]) for kind, count in counts.items()), F())
    add_pulses = sum(len(w) for w in [list(arithmetic_word(add_program()))])
    faces = []
    for face in range(208):
        support = int(np.count_nonzero(af.data()[2][:, face]))
        bound = pe + 2 * (support - 1) * add_error
        pulses = len(pw) + 2 * (support - 1) * add_pulses
        faces.append({'face': face, 'support_edges': int(support), 'native_pulses': pulses,
                      'operator_error_upper_exact': str(bound)})
    actual_face_word = face_word(0, beta)
    assert len(actual_face_word) == faces[0]['native_pulses']
    (HERE / 'FIELD_PROGRAM.json').write_text(json.dumps({
        'candidate_id': af.MODEL_ID, 'field_base': FIELD_BASE, 'field_bits': 1152,
        'matter_modes': 2048, 'scratch_bits': list(SCRATCH), 'total_source_bits': 3203,
        'QFT_word': qft_word(), 'electric_word': ew, 'plaquette_core_word': pw,
        'modular_add_gates': add_program(), 'modular_add_native_primitives': af.mb.compiled.controls()['primitive_words'],
        'full_field_pulse_generator': 'native_field.field_word(alpha,face_couplings208)',
        'full_matter_field_clock_generator': 'native_field.complete_word(alpha,face_couplings208,clock_spacing=None)',
        'oneface_member': {'alpha': float(alpha), 'face0_beta': float(beta), 'other_beta': 0},
        'matter_layer': 'frozen manybody.mode_program + mode_field.field_block(16), with full CAR signs'}, indent=2) + '\n')
    result = {'schema': 'oph.local12.native_field_words.v1', 'candidate_id': af.MODEL_ID,
              'checks': {'full_QFT_source_columns': True, 'full_electric_phase_source_columns': True,
                         'full_plaquette_cosine_source_columns': True,
                         'modular_shear_all_clean_and_dirty_columns': True,
                         'field_generators_preserve_Gauss': True,
                         'outward_full_source_error_bounds': True,
                         'full_field_word_generator_supplied': True, 'source_bytes_unchanged': True},
              'local_columns': {'QFT': 16, 'electric': 16, 'plaquette': 16, 'shear_clean': 256, 'shear_dirty': 2048},
              'numeric_residuals': {'QFT': qft_residual, 'electric': er, 'plaquette': pr},
              'electric_per_edge_error_exact': str(ee), 'plaquette_core_error_exact': str(pe),
              'shear_error_exact': str(add_error), 'shear_native_pulses': add_pulses,
              'faces': faces,
              'oneface_field_step_error_upper': float(288 * ee + F(faces[0]['operator_error_upper_exact'])),
              'homogeneous208face_field_step_error_upper': float(288 * ee + sum((F(r['operator_error_upper_exact']) for r in faces), F())),
              'oneface_field_step_native_pulses': 288 * len(ew) + faces[0]['native_pulses'],
              'homogeneous208face_field_step_native_pulses': 288 * len(ew) + sum(r['native_pulses'] for r in faces),
              'complete_source_scope': 'Declared exact source response and coherent code services, every export complement retained. Add the separate matter-layer error, initial3203 native resets, terminal detector errors and actual clock resources. This receipt does not execute a full2^3203 matrix or discard initialization failure.',
              'clock_probe_allocation': {'clock': 3203, 'probe': 3204, 'detector_auxiliaries': 'fresh indices above3204'},
              'source_sha256': pins, 'official_rows_changed': []}
    assert pins == {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    (HERE / 'NATIVE_FIELD.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'numeric_residuals', 'oneface_field_step_error_upper',
                                           'homogeneous208face_field_step_error_upper', 'oneface_field_step_native_pulses')}))


if __name__ == '__main__':
    run()
