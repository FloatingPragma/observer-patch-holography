"""Complete coherent C3 source and unconditioned weak phase feedback."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import math

import numpy as np
import sympy as s

from model.common_source import native_family_record_memory as memory
from model.common_source import native_repair_mass as repair
from model.common_source import native_matter_attachment as matter
from model.observer.source_scalar import compiler


P = s.ImmutableMatrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
WEIGHTS = (s.Rational(10, 72), s.Rational(5, 72), s.Rational(5, 72),
           s.Rational(2, 72), s.Rational(50, 72))
NATIVE_ADMISSION = 2e-8


def _clean(a):
    return s.ImmutableMatrix(a.applyfunc(s.simplify))


def _phase(z):
    z = s.sympify(z)
    if (z.has(s.Float) or z.free_symbols or s.re(z).is_Rational is not True
            or s.im(z).is_Rational is not True or s.expand(z*s.conjugate(z)) != 1):
        raise ValueError('an exact rational unit-circle phase is required')
    return z


def native_partition():
    rows = []; totals = [s.S.Zero]*5
    for cc, cr, ec, er, _, weight in repair.native_records():
        plus, minus, b = cc == 0 and cr == 1, cc == 1 and cr == 2, ec == 0 and er in (4, 5)
        group = (0 if b and not (plus or minus) else 1 if plus and not b
                 else 2 if minus and not b else 3 if b and (plus or minus) else 4)
        totals[group] += weight
        rows.append({'raw': [cc, cr, ec, er], 'group': group, 'probability': str(weight)})
    if tuple(totals) != WEIGHTS or len(rows) != 128:
        raise ArithmeticError('the full native exclusive partition differs')
    return rows


@lru_cache(maxsize=1)
def _native_rotation():
    data, (gc, ge) = memory.geometry()
    g = tuple(ge[gc[j]] for j in range(12))
    b = data['frame']; metric = (b.T*b)[0, 0]
    r = _clean(b.T*memory.objects.permutation_matrix(g)*b/metric)
    if _clean(r**3) != s.eye(3) or s.trace(r) != 0 or _clean(r.T*r) != s.eye(3):
        raise ArithmeticError('the actual C/E half-turn product is not regular C3')
    axis = (r-s.eye(3)).nullspace()
    a = axis[0]; a = _clean(a/s.sqrt((a.T*a)[0]))
    for j in range(3):
        v = s.eye(3)[:, j]-a*a[j]
        if (v.T*v)[0] != 0:
            v = _clean(v/s.sqrt((v.T*v)[0])); break
    v = _clean(a/s.sqrt(3)+s.sqrt(s.Rational(2, 3))*v)
    o = _clean(s.Matrix.hstack(v, r*v, r*r*v))
    if _clean(o.T*o) != s.eye(3) or _clean(r*o-o*P) != s.zeros(3):
        raise ArithmeticError('native and cyclic family frames do not intertwine')
    return r, o


def exact_source(z=1):
    z = _phase(z)
    a, b = s.sqrt(WEIGHTS[0]), s.sqrt(WEIGHTS[1])
    c = _clean(a*s.eye(3)+b*(z*P+s.conjugate(z)*P**2))
    y = _clean(c*c/3)
    if s.simplify(s.trace(c*c)-s.Rational(5, 6)) != 0 or s.trace(y) != s.Rational(5, 18):
        raise ArithmeticError('complete coherent source normalization differs')
    positive = all(s.simplify(c[:k, :k].det()).is_positive is True for k in (1, 2, 3))
    q = s.simplify(s.trace(c*c)/s.trace(c)**2)
    if q != s.Rational(2, 3):
        raise ArithmeticError('the signed balanced identity differs')
    return {'z': z, 'C_exclusive': c, 'action_coefficient': y,
            'positive_root_chamber_certified': positive,
            'signed_root_Koide_ratio': q,
            'physical_action_Koide_ratio_certified': positive}


@lru_cache(maxsize=8)
def analyzer_factors(z=1):
    """Exact full padding: eight flags and four family coordinates."""
    z = _phase(z)
    amp = s.Matrix([s.sqrt(w) for w in WEIGHTS]+[0, 0, 0])
    v = s.eye(8)[:, 0]-amp
    prep = _clean(s.eye(8)-2*v*v.T/(v.T*v)[0])
    w = -s.Rational(1, 2)+s.sqrt(3)*s.I/2
    fourier = s.eye(8)
    for i, j in product(range(3), repeat=2):
        fourier[i, j] = s.expand(w**(-i*j))/s.sqrt(3)
    fourier = _clean(fourier)
    r = s.eye(4); r[:3, :3] = P
    blocks = (s.eye(4), z*r, s.conjugate(z)*r*r)+tuple(s.eye(4) for _ in range(5))
    if _clean(prep*prep) != s.eye(8) or _clean(fourier.H*fourier) != s.eye(8):
        raise ArithmeticError('the complete flag service is not unitary')
    return prep, fourier, tuple(s.ImmutableMatrix(b) for b in blocks)


@lru_cache(maxsize=256)
def analyzer_block(out_flag, in_flag, z=1):
    if any(type(x) is not int or not 0 <= x < 8 for x in (out_flag, in_flag)):
        raise ValueError('all eight flag indices must be explicit')
    prep, fourier, blocks = analyzer_factors(z)
    return _clean(sum((fourier[out_flag, j]*prep[j, in_flag]*blocks[j]
                       for j in range(8)), s.zeros(4)))


def exact_analyzer_records(z=1):
    records = []
    for flag in range(8):
        k = analyzer_block(flag, 0, z)[:3, :3]
        prob = s.simplify(s.trace(k.H*k)/3)
        records.append({'flag': flag, 'probability': str(prob)})
    expected = [s.Rational(5, 54)]*3+[s.Rational(1, 36), s.Rational(25, 36)]+[0]*3
    if [s.sympify(x['probability']) for x in records] != expected:
        raise ArithmeticError('complete original tracial analyzer records differ')
    return records


def weak_kraus(z=1, tangent=s.Rational(1, 64)):
    """Exact kick angle h=2 atan(tangent); all eight reset flags retained."""
    tangent = s.sympify(tangent)
    if tangent.is_Rational is not True or tangent <= 0:
        raise ValueError('a positive rational half-angle tangent is required')
    co = (1-tangent*tangent)/(1+tangent*tangent)
    si = 2*tangent/(1+tangent*tangent)
    x = s.Matrix([[0, 1], [1, 0]])
    e = (co-1)*s.eye(2)-s.I*si*x
    first = analyzer_block(0, 0, z)[:3, :3]
    return tuple(_clean((s.eye(6) if flag == 0 else s.zeros(6))+
        s.kronecker_product(e, analyzer_block(0, flag, z)[:3, :3].H*first))
        for flag in range(8))


def weak_records(z=1):
    ks = weak_kraus(z)
    probabilities = [s.simplify((k.H*k)[0, 0]) for k in ks]
    if s.simplify(sum(probabilities)-1) != 0:
        raise ArithmeticError('weak-kick original records lost probability')
    return [{'flag': flag, 'probability': str(prob)} for flag, prob in enumerate(probabilities)]


def feedback_bound(steps, total_time, native_dilation_error=F(0)):
    if type(steps) is not int or steps <= 0 or type(total_time) not in (int, F) or total_time < 0:
        raise ValueError('positive integer steps and exact nonnegative total time required')
    t = F(total_time)
    if type(native_dilation_error) not in (int, F) or native_dilation_error < 0:
        raise ValueError('an exact certified nonnegative dilation error is required')
    return {'ideal_nonzero_archive_probability_upper': min(F(1), t*t/(4*steps)),
            'nonzero_archive_probability_upper': min(F(1), t*t/(4*steps)+steps*native_dilation_error),
            'complete_unhalved_trace_distance_upper': min(F(2), 9*t*t/(4*steps)+2*steps*native_dilation_error),
            'zero_record_postselection_used': False,
            'arbitrary_target_reference_inputs': True}


def _swap_word(a, b, bits):
    diff = a ^ b
    if diff == 0 or diff & (diff-1):
        raise ValueError('one-bit address transposition required')
    target = diff.bit_length()-1; word = []
    for mask in range(1 << bits):
        if mask & diff: continue
        sign = (-1)**((a & mask).bit_count())
        label = ['Z' if mask >> j & 1 else 'I' for j in range(bits)]
        angle = sign*math.pi/(1 << bits)
        word += compiler.pauli_rotation(label, angle)
        label[target] = 'X'
        word += compiler.pauli_rotation(label, -angle)
    return word


def analyzer_word(z=1):
    z = _phase(z); prep, _, _ = analyzer_factors(z)
    word = compiler.remap(compiler.compile_real_orthogonal(np.asarray(prep, float)), [2, 3, 4])
    phase = math.atan2(float(s.im(z)), float(s.re(z)))
    word += compiler.remap(compiler.compile_diagonal([0, phase, -phase, 0, 0, 0, 0, 0]), [2, 3, 4])
    for flag, pairs in ((1, ((0, 1), (0, 2))), (2, ((0, 2), (0, 1)))):
        for a, b in pairs: word += _swap_word(4*flag+a, 4*flag+b, 5)
    # F3 has real eigenvectors and eigenvalues 1,-1,-i.
    c = math.sqrt((1+1/math.sqrt(3))/2); d = math.sqrt((1-1/math.sqrt(3))/2)
    o = np.eye(8); o[:3, :3] = [[c, -d, 0], [d/math.sqrt(2), c/math.sqrt(2), 1/math.sqrt(2)],
                                              [d/math.sqrt(2), c/math.sqrt(2), -1/math.sqrt(2)]]
    basis = compiler.compile_real_orthogonal(o)
    word += compiler.remap(compiler.inverse(basis)+compiler.compile_diagonal(
        [0, math.pi, -math.pi/2, 0, 0, 0, 0, 0])+basis, [2, 3, 4])
    return word


def native_execution(z=1):
    """Literal same-current ideal native service; numerical diagnostics only."""
    word = analyzer_word(z)
    u = matter.execute_word(word, 5)
    inverse = matter.execute_word(compiler.inverse(word), 5)
    prep, f, blocks = analyzer_factors(z)
    controlled = s.diag(*blocks)
    ideal = np.asarray(s.kronecker_product(f, s.eye(4))*controlled*
                       s.kronecker_product(prep, s.eye(4)), complex)
    error = float(np.linalg.norm(u-ideal, 2))
    inverse_error = float(np.linalg.norm(inverse@u-np.eye(32), 2))
    h = 2*math.atan(1/64); kick = []
    for mask in range(8):
        label = ['I']*6; label[5] = 'X'
        for j in range(3):
            if mask >> j & 1: label[j+2] = 'Z'
        kick += compiler.pauli_rotation(label, -h/8)
    v = np.kron(np.eye(2), inverse)@matter.execute_word(kick, 6)@np.kron(np.eye(2), u)
    columns = [32*n+family for n in range(2) for family in range(3)]
    ks = [v[np.ix_([32*n+4*flag+family for n in range(2) for family in range(3)], columns)]
          for flag in range(8)]
    wanted = weak_kraus(z)
    weak_error = float(np.linalg.norm(np.vstack(ks)-np.vstack([np.asarray(k, complex) for k in wanted]), 2))
    q = np.zeros((32, 32)); q[:4, :4] = np.eye(4)
    projection = ideal.conj().T@q@ideal
    co, si = float(s.Rational(4095, 4097)), float(s.Rational(128, 4097))
    full_target = (np.eye(64)+(co-1)*np.kron(np.eye(2), projection)
                   -1j*si*np.kron([[0, 1], [1, 0]], projection))
    full_weak_error = float(np.linalg.norm(v-full_target, 2))
    gram_error = float(np.linalg.norm(v.conj().T@v-np.eye(64), 2))
    if max(error, inverse_error, weak_error, full_weak_error, gram_error) > NATIVE_ADMISSION:
        raise ArithmeticError('the complete native coherent word missed its ideal instrument')
    resources = matter.resources(word); kick_resources = matter.resources(kick)
    failures = sum(6-(4 if item['kind'] == 'J' else 2) for item in word+kick+compiler.inverse(word))
    return {'analyzer_native_pulses': len(word), 'weak_kick_native_pulses': len(kick),
            'one_round_native_pulses': 2*len(word)+len(kick),
            'complete_analyzer_columns': 32, 'complete_weak_round_columns': 64,
            'all_64_weak_columns_compared_with_exact_projection_formula': True,
            'retained_reset_flags': 8, 'indexed_first_service_failure_slots': failures,
            'valid_code_failure_weights_exact': '0',
            'primitive_source': 'native_matter_attachment.execute_word',
            'native_continuous_control_and_code_services_supplied': True,
            'floating_execution_passed_fixed_admission': True,
            'diagnostic_operator_admission': str(NATIVE_ADMISSION),
            'numerical_admission_is_a_rigorous_roundoff_certificate': False,
            'analyzer_single_code_pulses': resources['single_code_pulses'],
            'analyzer_pair_code_pulses': resources['pair_code_pulses'],
            'kick_single_code_pulses': kick_resources['single_code_pulses']}


def certificate():
    rows = native_partition(); r, o = _native_rotation()
    controls = []
    for z in (s.S.One, (s.Integer(99)+20*s.I)/101):
        d = exact_source(z)
        if not d['physical_action_Koide_ratio_certified']:
            raise ArithmeticError('the stated control is outside the positive chamber')
        controls.append({'phase_unit_complex_exact': str(z),
            'phase_selected_by_native_records': False,
            'C_exclusive_exact': [[str(x) for x in row] for row in d['C_exclusive'].tolist()],
            'action_coefficient_exact': [[str(x) for x in row] for row in d['action_coefficient'].tolist()],
            'positive_root_Koide_ratio_exact': '2/3',
            'action_trace_exact': '5/18',
            'analyzer_original_tracial_records': exact_analyzer_records(z),
            'weak_original_records_on_basis_input': weak_records(z)})
    return {'schema': 'oph.standard_model_source_map.koide_source.v1',
        'native_exclusive_partition': rows,
        'native_order_three_rotation_exact': [[str(x) for x in row] for row in r.tolist()],
        'native_to_cyclic_family_frame_exact': [[str(x) for x in row] for row in o.tolist()],
        'configured_phase_controls': controls,
        'literal_native_execution_at_zero_phase': native_execution(),
        'complete_feedback_bound': {k: str(v) if isinstance(v, F) else v
                                   for k, v in feedback_bound(1024, F(1)).items()},
        'same_M15_Higgs_family_handoff': 'Y_e=K_0^dagger K_0 on the same gauge-trivial C3 factor',
        'physical_CAR_action_handoff_is_a_configured_coefficient_map': True,
        'scope': {'different_configured_source_from_current_gated_Z_D': True,
            'complete_instrument_uses_postselection': False,
            'active_probability_divided_out_of_action': False,
            'conditional_preparation_and_action_coefficient_conflated': False,
            'balanced_mean_family_generator_proved': True,
            'full_physical_CAR_weak_kick_compiled': False,
            'physical_common_clock_selected': False,
            'weak_kick_angle_time_identified_with_global_native_attempt_time': False,
            'physical_pole_mass_or_running_to_pole_map_proved': False,
            'phase_or_mass_ratio_selected': False,
            'measured_mass_inputs_consumed': [], 'physical_particle_assignments': []}}
