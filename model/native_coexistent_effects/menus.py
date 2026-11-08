"""Sequential weak filters on the unchanged 66-dimensional source carrier."""
from __future__ import annotations
from fractions import Fraction as F
import numpy as np

from model.native_effect_synthesis import intervals as iv, synthesis as s

DIMENSION = 66
MAX_STAGES = 8
WORD_PULSES = 58392
ZERO = (F(), F())
ONE = (F(1), F())


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def mul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]


def conjugate(a):
    return a[0], -a[1]


def scale(a, x):
    return a[0]*x, a[1]*x


def identity():
    return {(i, i): ONE for i in range(DIMENSION)}


def adjoint(matrix):
    return {(j, i): conjugate(x) for (i, j), x in matrix.items()}


def matrix_sum(*matrices):
    out = {}
    for matrix in matrices:
        for key, value in matrix.items():
            out[key] = add(out.get(key, ZERO), value)
    return {key: value for key, value in out.items() if value != ZERO}


def product(left, right):
    rows = {}
    for (j, k), value in right.items():
        rows.setdefault(j, []).append((k, value))
    out = {}
    for (i, j), value in left.items():
        for k, other in rows.get(j, ()):
            key = i, k
            out[key] = add(out.get(key, ZERO), mul(value, other))
    return {key: value for key, value in out.items() if value != ZERO}


def specification(stages):
    """Validate a nonempty bounded menu without silently dropping zero tests."""
    if type(stages) not in (list, tuple) or not 1 <= len(stages) <= MAX_STAGES:
        raise ValueError('one through eight native filter stages required')
    out = []
    for stage in stages:
        if type(stage) not in (list, tuple) or len(stage) != 2:
            raise ValueError('each stage contains a vector and survival amplitude')
        vector, survival, _ = s.target(stage[0], stage[1])
        out.append((vector, survival))
    return tuple(out)


def projector(vector):
    vector, _, norm = s.target(vector)
    return {(i, j): scale(mul(x, conjugate(y)), 1/norm)
            for i, x in enumerate(vector) if x != ZERO
            for j, y in enumerate(vector) if y != ZERO}


def ideal_menu(stages):
    """Exact Gaussian-rational effects, success maps and residual instrument."""
    stages = specification(stages)
    prefix = identity(); effects = []; prefixes = []
    for vector, survival in stages:
        p = projector(vector)
        local = matrix_sum(identity(), {key: scale(x, survival-1) for key, x in p.items()})
        effect = product(product(adjoint(prefix), p), prefix)
        effects.append({key: scale(x, 1-survival**2) for key, x in effect.items()
                        if scale(x, 1-survival**2) != ZERO})
        prefixes.append(prefix)
        prefix = product(local, prefix)
    effects.append(product(adjoint(prefix), prefix))
    if matrix_sum(*effects) != identity():
        raise ArithmeticError('the complete sequential menu lost its residual')
    return {'effects': tuple(effects), 'prefixes': tuple(prefixes), 'residual': prefix}


def coarse_effect(effects, labels):
    if type(labels) not in (tuple, list) or len(set(labels)) != len(labels):
        raise ValueError('distinct complete menu labels required')
    if any(type(i) is not int or not 0 <= i < len(effects) for i in labels):
        raise ValueError('a complete menu outcome label required')
    return matrix_sum(*(effects[i] for i in labels))


def preparation(vector, padding_weight):
    vector, _, norm = s.target(vector)
    if type(padding_weight) is not F or not 0 <= padding_weight <= 1:
        raise ValueError('a rational padding weight in [0,1] required')
    codes = s.tomo.codes(); padding = next(i for i in range(512) if i not in codes)
    weights = [F()]*512; phases = [0]*512
    for code, (x, y) in zip(codes, vector, strict=True):
        weights[code] = (1-padding_weight)*(x*x+y*y)/norm
        phases[code] = iv.phase_angle(x, y)
    weights[padding] = padding_weight
    return s._state_word(weights)+s._phase_word(phases)


def phase_word(vector, padding_weight, angle_numerator):
    """Rank-one phase in the same 512 addresses, with separate H66 export."""
    if type(angle_numerator) is not int or abs(angle_numerator) > 4*iv.AQ:
        raise ValueError('a bounded exact grid phase numerator required')
    word = preparation(vector, padding_weight)
    phases = [0]*512; phases[0] = angle_numerator
    result = s.compiler.inverse(word)+s._phase_word(phases)+word
    if len(result) != WORD_PULSES:
        raise ArithmeticError('the admitted rank-one phase resource count changed')
    s.interpreter.source_word(result, 9)
    return result


def filter_words(vector, survival):
    """Two admitted words implementing a weak filter and its phase correction."""
    vector, survival, _ = s.target(vector, survival)
    theta = iv.split_angle(survival**2, 1-survival**2)
    return (phase_word(vector, F(1, 2), 2*theta),
            phase_word(vector, F(), -theta))


def compile_menu(stages):
    return tuple(filter_words(vector, survival) for vector, survival in specification(stages))


def word_bounds(pulses, target_grid_units):
    if type(pulses) is not int or not 1 <= pulses <= WORD_PULSES:
        raise ValueError('a bounded native phase word required')
    if type(target_grid_units) is not int or target_grid_units not in (23, 24):
        raise ValueError('the proved source phase target gate required')
    x = pulses*s.PULSE_GATE
    e = x/(1-x)+F(target_grid_units, iv.AQ)
    failure = pulses*s.PULSE_GATE**2/(1-x)**2
    return {'operator_error': e, 'processor_failure': failure,
            'halfdiamond': (e*(2+e)+failure)/2,
            'gain': 1/(1-x)**2+failure}


def menu_bounds(stage_count):
    if type(stage_count) is not int or not 1 <= stage_count <= MAX_STAGES:
        raise ValueError('one through eight complete filter stages required')
    delta = F(); gain = F(1)
    for _ in range(stage_count):
        for units in (24, 23):
            bound = word_bounds(WORD_PULSES, units)
            delta += gain*bound['halfdiamond']
            gain *= bound['gain']
    return {'complete_halfdiamond_upper': delta, 'complete_gain_upper': gain,
            'native_words': 2*stage_count, 'native_pulses': 2*WORD_PULSES*stage_count}


def source_columns(value):
    if (not isinstance(value, np.ndarray) or value.ndim != 2 or value.shape[0] != DIMENSION
            or value.shape[1] < 1 or value.dtype.kind not in 'iufc'
            or not np.isfinite(value).all()):
        raise ValueError('finite source by untouched-reference amplitude columns required')
    return np.array(value, dtype=complex, copy=True)


def successful_word(word, columns):
    """Diagnostic exact-stored-product approximation, with no renormalization.

    First processor failures are queried separately by failure_branch; this
    routine does not assert that the returned success/padding blocks are TP.
    """
    columns = source_columns(columns)
    s.interpreter.source_word(word, 9)
    work = np.zeros((512, columns.shape[1]), complex)
    work[list(s.tomo.codes())] = columns
    for event in word:
        full = s.compiler.native_pulse(event['kind'], event['angle'])
        dimension = 1 << len(event['targets'])
        work = s.compiler._apply(work, full[:dimension, :dimension], event['targets'], 9)
    padding = sorted(set(range(512))-set(s.tomo.codes()))
    # Rows are amplitudes of a retained reference in a bipartite pure input.
    abort = work[padding].T@work[padding].conj()
    return work[list(s.tomo.codes())], abort


def diagnostic_output(words, columns):
    if type(words) not in (tuple, list) or not 1 <= len(words) <= MAX_STAGES:
        raise ValueError('a bounded compiled menu required')
    current = source_columns(columns); aborts = []
    for pair in words:
        if type(pair) not in (tuple, list) or len(pair) != 2:
            raise ValueError('both actual phase words must be retained')
        for word in pair:
            current, abort = successful_word(word, current)
            aborts.append(abort)
    return tuple(aborts), current


def failure_branch(words, stage, component, pulse, complement, columns):
    """One terminal first-processor-failure map, including its complete prefix."""
    if type(stage) is not int or not 0 <= stage < len(words):
        raise ValueError('a native stage label required')
    if type(component) is not int or component not in (0, 1):
        raise ValueError('a weak or correction word label required')
    current = source_columns(columns)
    for i, pair in enumerate(words):
        for j, word in enumerate(pair):
            if (i, j) == (stage, component):
                return s.processor_failure(word, pulse, complement, current)
            current, _ = successful_word(word, current)
    raise ArithmeticError('a first-failure address was lost')


def padding_branch(words, stage, component, address, columns):
    """One internal padding Kraus row behind its single observed word flag."""
    if type(stage) is not int or not 0 <= stage < len(words):
        raise ValueError('a native stage label required')
    if type(component) is not int or component not in (0, 1):
        raise ValueError('a weak or correction word label required')
    if type(address) is not int or not 0 <= address < 512 or address in s.tomo.codes():
        raise ValueError('an internal padding Kraus address required')
    current = source_columns(columns)
    for i, pair in enumerate(words):
        for j, word in enumerate(pair):
            if (i, j) == (stage, component):
                return s.work_row(word, address)[list(s.tomo.codes())]@current
            current, _ = successful_word(word, current)
    raise ArithmeticError('a padding Kraus address was lost')


def outcomes(words):
    rows = []
    for stage, pair in enumerate(words):
        for component, word in enumerate(pair):
            admitted = s.interpreter.source_word(word, 9)
            rows.append({'stage': stage, 'component': component,
                         'processor_first_failure_labels': sum(6-(1 << len(e['targets'])) for e in word),
                         'observed_padding_abort_flags': 1, 'internal_padding_kraus_rows': 446,
                         'resources': admitted['declared_resources']})
    return {'words': rows, 'final_quantum_residual_labels': 1,
            'native_tree_measurements': 0, 'terminal_failure_stops_menu': True,
            'zero_failure_rows_retained': True, 'source_success_is_coherent': True}


def dense(matrix):
    out = np.zeros((DIMENSION, DIMENSION), complex)
    for (i, j), (x, y) in matrix.items():
        out[i, j] = complex(float(x), float(y))
    return out
