"""A fixed two-rotation family memory driven by complete native C/E records.

The E event is selected from actual nonzero, moving source branches before
computing a spectrum. Native edge half-turns act on a complete twelve-port
memory; no state is removed if it is unreachable from the reset. Coherent
projector feedback is a chosen effective service, not a physical mass map.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np
from scipy.linalg import expm
import sympy as s

from tools import strict_json
from . import native_repair_mass as feedback
from . import native_code_transport as transport
from . import native_gauge_object as objects

ROOT = objects.ROOT
BENCHMARK = 'model/common_source/chiral_candidate/results.json'
STEPS = 64
TIME = 1
ADMISSION = 2e-10
SCHEMA = {'cycle_transition': ((12, 12), 'float64'),
          'cycle_stationary': ((12,), 'float64'),
          'cycle_Y': ((3, 3), 'float64'),
          'cycle_pulses': ((12, 6, 6), 'complex128'),
          'cycle_labels': ((1536, 6), 'int64'),
          'cycle_raw_weights': ((128,), 'float64'),
          'cycle_finite_channel': ((36, 36), 'complex128')}


def _clean(value): return objects._clean(value)
def _zero(value): return not _clean(value).todok()
def _entries(value): return [[str(x) for x in row] for row in _clean(value).tolist()]


def moving_E_labels():
    """Inspect primitive source outputs, not corrected identity or mass data."""
    code = transport.PAIR_CODES['E']; supported, moving = [], []
    for raw in range(8):
        if code.branch_weight('E', (2, 10), 0, raw) <= 0: continue
        supported.append(raw)
        for matching in code.MATCHINGS:
            for pair in matching:
                state = tuple(int(i in pair) for i in range(12))
                out, emitted = transport._raw('E', state, (2, 10), 0)
                if emitted == raw and out != state: moving.append(raw)
    if not supported: raise ArithmeticError('no supported native E event')
    return tuple(supported), tuple(sorted(set(moving)))


@lru_cache(maxsize=1)
def architecture():
    data = objects.geometry(); B = data['frame']; rotations = []
    for a, b in ((2, 3), (2, 10)):
        candidates = [g for g in data['rotations']
                      if objects.permutation_order(g) == 2 and g[a] == b and g[b] == a]
        if len(candidates) != 1: raise ArithmeticError('the native edge half-turn is not unique')
        rotations.append(candidates[0])
    gC, gE = rotations
    if tuple(gC[gE[i]] for i in range(12)) == tuple(gE[gC[i]] for i in range(12)):
        raise ArithmeticError('the selected native edge responses commute')
    supported, moving = moving_E_labels(); er = (moving or supported)[0]
    records = feedback.native_records(); labels = []; T = s.zeros(12)
    for cc, cr, ec, raw, advance, weight in records:
        for old in range(12):
            new = 2 if not advance else gC[old]
            if advance and ec == 0 and raw == er: new = gE[new]
            T[new, old] += weight
            labels.append((cc, cr, ec, raw, old, new))
    p = sum(row[-1] for row in records if row[-2])
    q = sum(row[-1] for row in records if row[2:4] == (0, er))
    coincidence = sum(row[-1] for row in records if row[-2] and row[2:4] == (0, er))
    if coincidence != p*q: raise ArithmeticError('source event independence changed')
    C, E = map(objects.permutation_matrix, rotations)
    S = (1-q)*C+q*E*C
    reset = s.zeros(12); reset[2, :] = s.ones(1, 12)
    if T != (1-p)*reset+p*S or any(x < 0 for x in T) or s.ones(1, 12)*T != s.ones(1, 12):
        raise ArithmeticError('complete memory transition disagrees with source records')
    stationary = (1-p)*(s.eye(12)-p*S).inv()*s.eye(12)[:, 2]
    if T*stationary != stationary or sum(stationary) != 1 or any(x < 0 for x in stationary):
        raise ArithmeticError('stationary reset memory failed')
    projectors = tuple(_clean(B[i, :].T*B[i, :]/(B[i, :]*B[i, :].T)[0]) for i in range(12))
    Y = _clean(sum((stationary[i]*projectors[i] for i in range(12)), s.zeros(3)))
    if s.trace(Y) != 1 or any(not _zero(P*P-P) for P in projectors):
        raise ArithmeticError('source family kernel lost normalization')
    return {'gC': gC, 'gE': gE, 'supported_E': supported, 'moving_E': moving, 'E_event': (0, er),
            'p': p, 'q': q, 'coincidence': coincidence, 'T': s.ImmutableMatrix(T),
            'S': s.ImmutableMatrix(S), 'stationary': s.ImmutableMatrix(stationary),
            'projectors': projectors, 'Y': Y, 'labels': tuple(labels),
            'rank': Y.rank(), 'zero_E_control_weight': sum(row[-1] for row in records if row[2:4] == (0, 1))}


def _sign(value):
    """Exact sign in Q(sqrt5), using a checked rational square-root bracket."""
    value = s.expand(s.radsimp(value)); b = value.coeff(s.sqrt(5)); a = s.expand(value-b*s.sqrt(5))
    if not a.is_Rational or not b.is_Rational: raise ArithmeticError('expected Q(sqrt5) coefficient')
    lo = s.Rational(2236067977499789696409173668731, 10**30); hi = lo+s.Rational(1, 10**30)
    if not lo**2 < 5 < hi**2: raise ArithmeticError('square-root enclosure changed')
    lower, upper = (a+b*lo, a+b*hi) if b >= 0 else (a+b*hi, a+b*lo)
    if lower > 0: return 1
    if upper < 0: return -1
    if a == b == 0: return 0
    raise ArithmeticError('sign not resolved by exact rational enclosure')


def spectral_certificate():
    a = architecture(); Y = a['Y']; x = s.Symbol('x')
    polynomial = s.Poly(s.expand(Y.charpoly(x).as_expr()), x)
    if a['rank'] != 3:
        return {'rank': a['rank'], 'three_positive_couplings': False,
                'characteristic_polynomial': str(polynomial.as_expr()), 'root_intervals_exact': [], 'ratios_after_one_exact': []}
    # Rational Sturm bisection makes retained endpoints independent of an
    # eigensolver, its ordering, and platform rounding near a decimal grid.
    sequence = s.Poly(polynomial.as_expr(), x, extension=s.sqrt(5)).sturm()
    coefficients = []
    for row in sequence:
        pairs = []
        for value in row.all_coeffs():
            value = s.expand(s.radsimp(value)); b = value.coeff(s.sqrt(5)); c = s.expand(value-b*s.sqrt(5))
            if not c.is_Rational or not b.is_Rational: raise ArithmeticError('Sturm coefficient left Q(sqrt5)')
            pairs.append((F(c), F(b)))
        coefficients.append(pairs)
    sqrtlo = F(2236067977499789696409173668731, 10**30); sqrthi = sqrtlo+F(1, 10**30)
    if not sqrtlo**2 < 5 < sqrthi**2: raise ArithmeticError('sqrt5 bracket failed')

    def variations(point):
        signs = []
        for row in coefficients:
            c = b = F(0)
            for dc, db in row: c, b = point*c+dc, point*b+db
            lower, upper = (c+b*sqrtlo, c+b*sqrthi) if b >= 0 else (c+b*sqrthi, c+b*sqrtlo)
            if lower > 0: signs.append(1)
            elif upper < 0: signs.append(-1)
            elif c or b: raise ArithmeticError('Sturm sign enclosure unresolved')
        return sum(a != b for a, b in zip(signs, signs[1:]))

    left, right = variations(F(0)), variations(F(1))
    if left-right != 3: raise ArithmeticError('source cubic must have three roots in (0,1)')
    pending = [(F(0), F(1), left, right)]; intervals = []
    while pending:
        lo, hi, vlo, vhi = pending.pop()
        if vlo == vhi: continue
        if vlo-vhi == 1 and hi-lo <= F(1, 10**12):
            intervals.append((s.Rational(lo), s.Rational(hi))); continue
        mid = (lo+hi)/2; vmid = variations(mid)
        pending.extend(((lo, mid, vlo, vmid), (mid, hi, vmid, vhi)))
    intervals.sort()
    for i, (lo, hi) in enumerate(intervals):
        if lo <= 0 or _sign(polynomial.eval(lo)) != (-1 if i % 2 == 0 else 1) or _sign(polynomial.eval(hi)) != (1 if i % 2 == 0 else -1):
            raise ArithmeticError('positive simple source-root endpoint certificate failed')
    ratios = [(lo/intervals[0][1], hi/intervals[0][0]) for lo, hi in intervals[1:]]
    return {'rank': 3, 'three_positive_couplings': True,
            'characteristic_polynomial': str(polynomial.as_expr()),
            'root_intervals_exact': [list(map(str, row)) for row in intervals],
            'ratios_after_one_exact': [list(map(str, row)) for row in ratios],
            'ratios_diagnostic': [list(map(float, row)) for row in ratios],
            'proof': 'deterministic rational Sturm bisection in Q(sqrt5); three positive sign-changing intervals exhaust the cubic'}


def _pulses():
    X = np.asarray([[0, 1], [1, 0]], complex); h = TIME/STEPS
    return np.asarray([np.eye(6)+(math.cos(h)-1)*np.kron(np.eye(2), np.asarray(P, complex))
                       -1j*math.sin(h)*np.kron(X, np.asarray(P, complex)) for P in architecture()['projectors']])


def finite_channel(T, pulses):
    blocks = [np.kron(U, U.conj()) for U in pulses]
    channel = np.block([[T[new, old]*blocks[new] for old in range(12)] for new in range(12)])
    initial = np.vstack([float(pi)*np.eye(36) for pi in architecture()['stationary']])
    return (np.linalg.matrix_power(channel, STEPS)@initial).reshape(12, 36, 36).sum(axis=0)


def validate_arrays(arrays):
    if set(arrays) != set(SCHEMA): raise ValueError('complete cycle-memory arrays required')
    for key, (shape, dtype) in SCHEMA.items():
        value = arrays[key]
        if type(value) is not np.ndarray or value.shape != shape or value.dtype != np.dtype(dtype) or not np.isfinite(value).all():
            raise ValueError('invalid cycle-memory array '+key)
    a = architecture()
    if not np.array_equal(arrays['cycle_labels'], np.asarray(a['labels'], np.int64)):
        raise ValueError('old memory, new memory or raw labels changed')
    for key, value in [('cycle_transition', a['T']), ('cycle_stationary', list(a['stationary'])),
                       ('cycle_Y', a['Y']), ('cycle_raw_weights', [row[-1] for row in feedback.native_records()])]:
        if np.max(abs(arrays[key]-np.asarray(value, float))) > ADMISSION: raise ValueError('source cycle coefficient differs: '+key)
    if np.max(abs(arrays['cycle_pulses']-_pulses())) > ADMISSION: raise ValueError('feedback pulse differs')
    if np.max(abs(arrays['cycle_finite_channel']-finite_channel(arrays['cycle_transition'], arrays['cycle_pulses']))) > ADMISSION:
        raise ValueError('complete correlated finite channel differs')


def metrics(arrays):
    validate_arrays(arrays)
    U = expm(-1j*TIME*np.kron([[0, 1], [1, 0]], arrays['cycle_Y']))
    difference = arrays['cycle_finite_channel']-np.kron(U, U.conj())
    choi = difference.reshape(6, 6, 6, 6).transpose(0, 2, 1, 3).reshape(36, 36)/6
    # ||Phi||diamond <= ||J_unnormalized||1 <= 6||J_unnormalized||F.
    diamond_upper = 36*float(np.linalg.norm(choi))
    p = float(architecture()['p']); bound = 4*(1+p)*TIME**2/(STEPS*(1-p))
    grams = np.zeros((12, 6, 6), complex)
    for index, row in enumerate(arrays['cycle_labels']):
        U = arrays['cycle_pulses'][row[5]]
        grams[row[4]] += arrays['cycle_raw_weights'][index//12]*(U.conj().T@U)
    defects = [float(np.linalg.norm(U.conj().T@U-np.eye(6), 2)) for U in arrays['cycle_pulses']]
    if max(defects+[float(np.max(abs(grams-np.eye(6))))]) > ADMISSION or diamond_upper > bound:
        raise ValueError('complete instrument or controlled averaging gate failed')
    return {'pulse_unitarity_defect_diagnostic': max(defects),
            'averaged_channel_diamond_upper_diagnostic': diamond_upper,
            'stationary_product_averaging_bound': bound,
            'steps': STEPS, 'dimensionless_model_time': TIME,
            'complete_stationary_product_finite_CQ_evolution': True,
            'iid_replacement_of_correlated_memory': False}


def execute():
    a = architecture(); spectrum = spectral_certificate(); pulses = _pulses()
    arrays = {'cycle_transition': np.asarray(a['T'], float), 'cycle_stationary': np.asarray(list(a['stationary']), float),
              'cycle_Y': np.asarray(a['Y'], float), 'cycle_pulses': pulses,
              'cycle_labels': np.asarray(a['labels'], np.int64),
              'cycle_raw_weights': np.asarray([row[-1] for row in feedback.native_records()], float)}
    arrays['cycle_finite_channel'] = finite_channel(arrays['cycle_transition'], pulses)
    execution = metrics(arrays)
    # This is the first access to empirical comparison data.
    comparison = strict_json.load(ROOT/BENCHMARK)['separate_empirical_comparison']
    windows = [tuple(F(x) for x in pair) for pair in comparison['conservative_ratio_windows']]
    if windows != [(F(206), F(207)), (F(3470), F(3485))]: raise ArithmeticError('retained comparison windows changed')
    ratios = [tuple(F(x) for x in pair) for pair in spectrum['ratios_after_one_exact']]
    passes = [max(lo, a) <= min(hi, b) for (lo, hi), (a, b) in zip(ratios, windows)] if ratios else [False, False]
    report = {'rule': {'C_event': [0, 1], 'E_event': list(a['E_event']), 'reset_port': 2,
              'E_selection': 'first positive coin0 branch moving a primitive matching state; otherwise first positive coin0 branch',
              'supported_E_raw_labels': list(a['supported_E']), 'moving_E_raw_labels': list(a['moving_E']),
              'C_probability_exact': str(a['p']), 'E_probability_exact': str(a['q']), 'coincidence_exact': str(a['coincidence']),
              'unsupported_E_0_1_negative_control_weight': str(a['zero_E_control_weight']),
              'C_edge_rotation': list(a['gC']), 'E_edge_rotation': list(a['gE']),
              'update': 'C failure resets2; C success applies gC, then gE only on selected E event',
              'transition_exact': _entries(a['T']), 'stationary_exact': list(map(str, a['stationary'])),
              'all12_memory_inputs_retained': True, 'raw_record_labels': 128, 'complete_record_memory_labels': 1536,
              'unreachable_stationary_ports_retained': [i for i, pi in enumerate(a['stationary']) if pi == 0]},
              'kernel': {'Y_exact': _entries(a['Y']), **spectrum}, 'execution': execution,
              'averaging': {'zero_sum_memory_contraction': 'T eta=p S eta; S is a convex sum of two permutations, ||T eta||1<=p||eta||1',
              'all_initial_classical_memory_reference_bound': '[4pt+4(1+p)t^2]/[N(1-p)]',
              'stationary_product_removes_linear_term': True},
              'separate_comparison': {'source': comparison['source'], 'ratio_window_intersections': passes,
              'joint_pass': all(passes), 'windows_are_experimental_standard_uncertainties': False},
              'scope': {'empirical_mass_inputs_to_source': False, 'source_native_pulse_history_executed': False,
              'fixed_rule_is_unchanged_C_E_alone': False, 'physical_fermion_mass_map_derived': False,
              'sector_dependent_Higgs_couplings_derived': False, 'measured_spectrum_claimed': False}}
    return report, arrays


def source_bindings():
    names = set(transport.source_bindings()) | set(objects.source_bindings()) | {
        'model/common_source/native_repair_mass.py', 'model/common_source/test_native_repair_mass.py',
        'model/common_source/native_mass_record_cycles.py', 'model/common_source/test_native_mass_record_cycles.py', BENCHMARK}
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}
