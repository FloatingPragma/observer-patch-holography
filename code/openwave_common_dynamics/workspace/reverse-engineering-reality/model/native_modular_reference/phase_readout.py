"""A signed quadrature readout from two complete native current records.

Prepare the usual shared reference, with one unread CORE attempt as in
reference_tomography. A cap Ramsey stage uses either the real analyzer or
exp(-i*pi*X_ab/4), compiled through the existing source Y and Z controls.
After the cap current (c,j), read CORE=(0,1), obtaining coin d/current k.
The score 2*j*1[d=c,k=-c] isolates spectator occupation n0=0. Its two ideal
means are A*cos(s*Delta) and A*sin(s*Delta), where
A=(9/2)*success_mass*(a1-a0) and Delta=log(w_hat1/w_hat0).
The amplitude is unknown; no reference weights enter the score or controls.

Every zero, unmatched coin, preparation failure and export abort remains.
The second read acts on the first read's disturbed carrier. The complete
quantum instrument works with untouched references. The population helper
is explicitly the ideal-source instrument: its native read destroys the
only coherences created from a diagonal input. Raw source arrays are
executed and certified separately, without a dephasing or normalization fix.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math

import numpy as np

from . import cap_local_clock as cap, learned_clock as learned
from . import reference_tomography as tomography

common, full, local = cap.common, cap.full, cap.local
history, interface, compiler = cap.history, cap.interface, cap.compiler
QUADRATURES = ('cos', 'sin')
ROTATION_GATE, WORD_GATE = cap.ROTATION_GATE, cap.WORD_GATE
PREPARATION_FLAGS = ('herald_failure', 'warmup_timeout', 'completion_timeout')


def _setting(index, quadrature):
    cap._probe(index)
    if type(quadrature) is not str or quadrature not in QUADRATURES:
        raise ValueError('an original cos/sin quadrature is required')
    return index, quadrature


def _law(law):
    if type(law) is not str or law not in ('C', 'E'):
        raise ValueError('an original native C/E law is required')
    return law


def quadrature_pair_word(a, b):
    """Existing pair Y word conjugated by two existing pivot Z pulses.

    Chronological Z(+pi/4), R_Y, Z(-pi/4) gives exp(-i*pi*X/4).
    The pivot orientation accounts for the original full code addresses.
    The Z factors cancel on every address outside this one pair.
    """
    if (type(a) is not int or type(b) is not int or not 0 <= a < 512
            or not 0 <= b < 512 or a == b):
        raise ValueError('two distinct original source addresses required')
    pivot = ((a ^ b) & -(a ^ b)).bit_length()-1
    orientation = 1 if not (a >> pivot & 1) else -1
    angle = orientation*math.pi/4
    return ([compiler.pulse('Z', angle, [pivot])]
            +cap.pair_word(a, b, 1)
            +[compiler.pulse('Z', -angle, [pivot])])


def quadrature_word(index):
    cap._probe(index)
    return [pulse for _, _, a, b in cap.pairs(index) for pulse in quadrature_pair_word(a, b)]


def ideal_quadrature(index):
    """Full 512-address target, identity on all inactive and padding states."""
    cap._probe(index)
    matrix = np.eye(512, dtype=complex)
    for _, _, a, b in cap.pairs(index):
        matrix[np.ix_((a, b), (a, b))] = np.array([[1, -1j], [-1j, 1]])/math.sqrt(2)
    return matrix


@lru_cache(maxsize=2)
def _quadrature(index):
    matrix = np.eye(512, dtype=complex)
    for _, _, a, b in cap.pairs(index):
        matrix = compiler.execute(quadrature_pair_word(a, b), 9)@matrix
    return matrix.dtype.str, matrix.tobytes()


def quadrature_rotation(index):
    cap._probe(index)
    dtype, data = _quadrature(index)
    return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)


def certify_quadrature(index, array=None):
    """Exact rational enclosure of every unmodified original source entry."""
    cap._probe(index)
    array = quadrature_rotation(index) if array is None else array
    if type(array) is not np.ndarray:
        raise ValueError('an original base source array is required')
    array = full.exact._array(array, (512, 512))
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    middle, radius = (lo+hi)/2, (hi-lo)/2
    entries = {}
    for _, _, a, b in cap.pairs(index):
        entries.update({(a, a): (middle, F()), (b, b): (middle, F()),
                        (a, b): (F(), -middle), (b, a): (F(), -middle)})
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            re, im = full.exact._component(value)
            target_re, target_im = entries.get((i, j), (F(int(i == j)), F()))
            squared += (re-target_re)**2+(im-target_im)**2
    return full.exact._sqrt_interval(squared)[1]+full.exact._sqrt_interval(F(40))[1]*radius


def source_bounds():
    error, gain = F(), F(1)
    for epsilon in (ROTATION_GATE, WORD_GATE, ROTATION_GATE):
        error = (1+epsilon)**2*error+epsilon*(2+epsilon)/2
        gain *= (1+epsilon)**2
    return error, gain


@lru_cache(maxsize=16)
def _certificate(weights, s, index, quadrature):
    full.corner.report._gate(cap.certify_rotation(index, -1), ROTATION_GATE, 'phase preparation')
    analyzer_error = (cap.certify_rotation(index, 1) if quadrature == 'cos'
                      else certify_quadrature(index))
    full.corner.report._gate(analyzer_error, ROTATION_GATE, 'phase quadrature analyzer')
    unitary = learned.clock(weights, s)
    full.corner.report._gate(common.certify_diagonal(s, learned.ratios(weights), unitary),
                            WORD_GATE, 'record-inferred phase word')
    word = (cap.rotation_word(index, -1)
            +common.diagonal_word(s, learned.ratios(weights))
            +(cap.rotation_word(index, 1) if quadrature == 'cos' else quadrature_word(index)))
    error, gain = source_bounds()
    return tuple({'source_history_halfdiamond_upper': error,
                  'source_channel_norm_upper': gain, 'source_pulses': len(word),
                  'proper_code_service_events': 2*len(word),
                  'native_current_reads': 2, 'native_exports': 1,
                  'source_word_sha256': hashlib.sha256(full.corner.canonical(word).encode()).hexdigest(),
                  'intermediate_pulses_proved_spatially_local': False}.items())


def certificate(weights, s=F(1), *, index=0, quadrature='cos'):
    weights, s = learned._weights(weights), full.parameter(s)
    _setting(index, quadrature)
    # Fraction's deepcopy returns the same instance, despite its writable slots.
    # Detach every exact scalar before exposing the cached certificate.
    return {key: F(value.numerator, value.denominator) if type(value) is F else value
            for key, value in _certificate(weights, s, index, quadrature)}


def complete_trial_budget():
    """Uniform scalars; source-array admission is a separate required premise.

    A score in [-2,2] costs at most 4T for half-trace error T. Its square has
    norm4, so 8T is safe even for non-TP floating source arrays. For the
    ideal charge family E[X^2]=9*q*(a0+a1)<=9*q/20.
    """
    preparation, q, maximum = tomography._record_contract()
    source_error, gain = source_bounds()
    total = gain*preparation+q*source_error
    return {'complete_trial_halftrace_upper': total,
            'score_mean_bias_upper': 4*total,
            'score_second_moment_upper': 9*q/20+8*total,
            'ideal_score_second_moment_upper': 9*q/20,
            'ideal_preparation_success_probability': F(q.numerator, q.denominator),
            'source_history_halfdiamond_upper': source_error,
            'source_channel_norm_upper': gain,
            'native_attempts_per_trial_upper': maximum+3,
            'unread_core_attempts': 1,
            'native_current_reads': 2, 'native_exports': 1,
            'raw_source_certification_required': True}


def stage(weights, s=F(1), *, index=0, quadrature='cos', raw=True):
    weights, s = learned._weights(weights), full.parameter(s)
    _setting(index, quadrature)
    if type(raw) is not bool:
        raise ValueError('an original Boolean source switch is required')
    if raw:
        certificate(weights, s, index=index, quadrature=quadrature)
        analyzer = cap.rotation(index, 1) if quadrature == 'cos' else quadrature_rotation(index)
        return analyzer@learned.clock(weights, s)@cap.rotation(index, -1)
    ratios = learned.ratios(weights)
    phases = {ratio: common.phase_enclosure(s, ratio) for ratio in set(ratios)}
    unitary = np.diag([complex(float(phases[r][0]), float(phases[r][1])) for r in ratios])
    analyzer = cap.ideal_rotation(index, 1) if quadrature == 'cos' else ideal_quadrature(index)
    return analyzer@unitary@cap.ideal_rotation(index, -1)


def record_alphabet(law, index=0, *, prepared=False):
    _law(law)
    cap._probe(index)
    if type(prepared) is not bool:
        raise ValueError('an original Boolean preparation switch is required')
    records = [('abort', index)]
    for first in interface.fixed_maps(law, cap.SEAMS[index]):
        for second in interface.fixed_maps(law, local.CORE):
            records.append(('read', first.coin, first.label, first.current,
                            second.coin, second.label, second.current))
    if prepared:
        records += [('preparation', flag) for flag in PREPARATION_FLAGS]
    return tuple(records)


def score(record):
    if type(record) is not tuple or not record:
        raise ValueError('one original complete phase-read record is required')
    if record == ('numerical_abort',):
        return 0
    if record[0] == 'abort' and len(record) == 2 and type(record[1]) is int and record[1] in (0, 1):
        return 0
    if record[0] == 'preparation' and len(record) == 2 and type(record[1]) is str and record[1] in PREPARATION_FLAGS:
        return 0
    if (record[0] != 'read' or len(record) != 7 or any(type(x) is not int for x in record[1:])
            or record[1] not in (0, 1) or record[4] not in (0, 1)
            or any(record[i] not in (-1, 0, 1) for i in (2, 3, 5, 6))):
        raise ValueError('the two original native coin/label/current records are required')
    _, coin, _, current, later_coin, _, later_current = record
    return 2*current if later_coin == coin and later_current == -coin else 0


def execute(density, weights, law, *, index=0, quadrature='cos', s=F(1), raw=True):
    """Two actual native reads after a certified stage, retaining poststates."""
    density = history._density(density)
    _law(law)
    _setting(index, quadrature)
    first, abort = cap._read_stage(density, index, law,
        stage(weights, s, index=index, quadrature=quadrature, raw=raw))
    result = {('abort', index): abort}
    ref = len(density)//66
    for record, rho in first.items():
        for branch in interface.fixed_maps(law, local.CORE):
            matrix = np.kron(branch.matrix, np.eye(ref))
            result[('read', *record, branch.coin, branch.label, branch.current)] = matrix@rho@matrix.conj().T
    if set(result) != set(record_alphabet(law, index)) or any(not np.isfinite(r).all() for r in result.values()):
        raise ArithmeticError('complete quadrature history lost a record or finite output')
    return result


def prepared(weights, law, *, initial=0, index=0, quadrature='cos', s=F(1), raw=True):
    """Finite preparation plus unread CORE, then the two-read phase probe."""
    if type(initial) is not int or not 0 <= initial < 66:
        raise ValueError('an original initial native index is required')
    _law(law)
    maps = local.population_maps(law=law)
    result = {('preparation', flag): np.diag(maps[flag][:, initial]).astype(complex)
              for flag in PREPARATION_FLAGS}
    successful = (tomography.population_counts()[0]/2)@maps['success'][:, initial]
    result.update(execute(np.diag(successful).astype(complex), weights, law,
                          index=index, quadrature=quadrature, s=s, raw=raw))
    return result


def _combined_branches(law, index):
    result = []
    for first in interface.fixed_maps(law, cap.SEAMS[index]):
        for second in interface.fixed_maps(law, local.CORE):
            later = dict(second.mapping)
            mapping = tuple((x, later[y]) for x, y in first.mapping if y in later)
            record = ('read', first.coin, first.label, first.current,
                      second.coin, second.label, second.current)
            result.append((record, mapping))
    return tuple(result)


def population_instrument(weights, s=F(1), *, index=0, quadrature='cos', law='C'):
    """Ideal-source record/state kernel before preparation or unread CORE.

    Rows are consecutive 66-state blocks in returned record order; columns
    are initial native basis states. The export block is zero. This is an
    exact diagonal-closure model of ideal operations, evaluated numerically.
    It is NOT a dephased substitute for the continuing raw-source channel.
    """
    _law(law)
    _setting(index, quadrature)
    records = record_alphabet(law, index)
    codes = np.array(full.support()[0])
    ideal = stage(weights, s, index=index, quadrature=quadrature, raw=False)
    transition = np.abs(ideal[np.ix_(codes, codes)])**2
    arrays = [np.zeros((66, 66))]
    for record, mapping in _combined_branches(law, index):
        matrix = np.zeros((66, 66))
        for before, after in mapping:
            matrix[after] += transition[before]/4
        arrays.append(matrix)
    kernel = np.vstack(arrays)
    if (not np.isfinite(kernel).all() or np.any(kernel < 0)
            or np.max(np.abs(kernel.sum(axis=0)-1)) > 2e-13):
        raise ArithmeticError('ideal quadrature population instrument is incomplete')
    return records, kernel


def response_enclosure(weights, reference_weights, s=F(1), *, quadrature='cos', success_mass=F(1)):
    """Diagnostic target formula with explicitly supplied reference weights.

    The decoder must not call this validation helper with an analytical
    reference. Controls and the score do not use reference_weights.
    """
    weights, reference_weights = learned._weights(weights), learned._weights(reference_weights)
    _, quadrature = _setting(0, quadrature)
    s = full.parameter(s)
    if (type(success_mass) is not F or not 0 <= success_mass <= 1
            or sum(g*w for g, w in zip((45, 20, 1), reference_weights, strict=True)) != 1):
        raise ValueError('normalized exact charge weights and a probability are required')
    amplitude = F(9, 2)*success_mass*(reference_weights[1]-reference_weights[0])
    re, im, radius = common.phase_enclosure(s, weights[1]/weights[0])
    center = amplitude*(re if quadrature == 'cos' else -im)
    error = abs(amplitude)*radius
    return center-error, center+error


def _pauli_certificate():
    """Exact Gaussian-integer conjugations; normalized gates divide by two."""
    def matrix(rows):
        return tuple(tuple((int(complex(z).real), int(complex(z).imag)) for z in row) for row in rows)
    def multiply(a, b):
        return tuple(tuple((sum(a[i][k][0]*b[k][j][0]-a[i][k][1]*b[k][j][1] for k in range(2)),
                            sum(a[i][k][0]*b[k][j][1]+a[i][k][1]*b[k][j][0] for k in range(2)))
                           for j in range(2)) for i in range(2))
    def adjoint(a):
        return tuple(tuple((a[j][i][0], -a[j][i][1]) for j in range(2)) for i in range(2))
    def twice(a):
        return tuple(tuple((2*x, 2*y) for x, y in row) for row in a)
    x = matrix(((0, 1), (1, 0)))
    y = matrix(((0, -1j), (1j, 0)))
    z = matrix(((1, 0), (0, -1)))
    prepare = matrix(((1, -1), (1, 1)))
    cosine = matrix(((1, 1), (-1, 1)))
    sine = matrix(((1, -1j), (-1j, 1)))
    for gate, before, after in ((prepare, z, x), (cosine, x, z),
                                (cosine, y, y), (sine, y, z), (sine, x, x)):
        if multiply(multiply(gate, before), adjoint(gate)) != twice(after):
            raise ArithmeticError('exact quadrature Pauli response or sign changed')
    return True


def _audit_native_branches(law, seam):
    """Bind integer effects to the complete arrays and actual native records."""
    states = local.native.states(binary=True)
    branches = interface.fixed_maps(law, seam)
    expected = {}
    for coin in (0, 1):
        for i, state in enumerate(states):
            output, record = local.native.native(state, seam, coin)
            current = state[seam[0]]-output[seam[0]]
            label = current if law == 'C' else state[seam[0]]-state[seam[1]]
            if record != (coin, current):
                raise ArithmeticError('native source record changed')
            expected.setdefault((coin, label, current), []).append((i, states.index(output)))
    if len(branches) != len(expected):
        raise ArithmeticError('native read omitted or duplicated a complete branch')
    seen = set()
    for branch in branches:
        key = branch.coin, branch.label, branch.current
        if key in seen or key not in expected or branch.mapping != tuple(expected[key]):
            raise ArithmeticError('native read mapping or complete record changed')
        seen.add(key)
        target = np.zeros((66, 66), complex)
        for before, after in expected[key]:
            target[after, before] = 1/math.sqrt(2)
        if type(branch.matrix) is not np.ndarray or not np.array_equal(branch.matrix, target):
            raise ArithmeticError('native read array is not its declared fair-coin partial injection')
    return branches


def exact_effects():
    """Integer source proof of score effects, variance, locality and closure."""
    if local.CORE != (0, 1) or any(seam[0] != 1 for seam in cap.SEAMS):
        raise ArithmeticError('phase readout requires the declared shared-core wiring')
    _pauli_certificate()
    states = local.native.states(binary=True)
    rows = []
    for index, seam in enumerate(cap.SEAMS):
        generator = {address: 1 for i, j, _, _ in cap.pairs(index)
                     for address in ((i, j), (j, i))}
        if not cap._in_algebra(generator, seam):
            raise ArithmeticError('quadrature generator is not in its bare seam algebra')
        supports = {i: {i} for i in range(66)}
        for i, j, _, _ in cap.pairs(index):
            supports[i] = supports[j] = {i, j}
        for law in ('C', 'E'):
            effects, second = [F()]*66, [F()]*66
            audited = _audit_native_branches(law, seam)
            _audit_native_branches(law, local.CORE)
            for branch in audited:
                domain = {i for i, _ in branch.mapping}
                if any(len(domain & support) > 1 for support in supports.values()):
                    raise ArithmeticError('ideal phase read does not close on diagonal inputs')
            branches = _combined_branches(law, index)
            for record, mapping in branches:
                inputs, outputs = zip(*mapping, strict=True) if mapping else ((), ())
                if len(set(inputs)) != len(inputs) or len(set(outputs)) != len(outputs):
                    raise ArithmeticError('two-native-read branch is not a partial injection')
                for i in inputs:
                    effects[i] += F(score(record), 4)
                    second[i] += F(score(record)**2, 4)
            wanted = tuple(F((1-x[0])*(x[seam[0]]-x[seam[1]]), 2) for x in states)
            wanted_second = tuple(F((1-x[0])*abs(x[seam[0]]-x[seam[1]])) for x in states)
            if tuple(effects) != wanted or tuple(second) != wanted_second:
                raise ArithmeticError('actual native score lost its spectator sector')
            # Every recorded map is injective: all off-diagonal effects vanish.
            for left in range(66):
                for right in range(66):
                    for record, mapping in branches:
                        image = dict(mapping)
                        if left != right and left in image and right in image and image[left] == image[right]:
                            raise ArithmeticError('an off-diagonal score effect survived')
            active = [(i, j) for i, j, _, _ in cap.pairs(index) if states[i][0] == 0]
            if len(active) != 9:
                raise ArithmeticError('the common-gap score must include nine spectators')
            rows.append({'index': index, 'law': law, 'matrix_units': 4356,
                         'native_read_branches': len(branches), 'selected_spectator_pairs': len(active),
                         'diagonal_input_closes_after_native_read': True,
                         'actual_native_maps_and_complete_arrays_bound': True})
    return {'rows': rows, 'score_values': (-2, 0, 2),
            'native_read_input_score_effect': '(1-N0)*(N1-Nb)/2',
            'native_read_input_squared_score_effect': '(1-N0)*abs(N1-Nb)',
            'ideal_mean': '(9/2)*success_mass*(a1-a0)*cos_or_sin(s*log(w_hat1/w_hat0))',
            'ideal_second_moment': '9*success_mass*(a0+a1)',
            'exact_gaussian_integer_pauli_conjugations_checked': True,
            'native_current_and_coin_only': True, 'both_reads_disturb_the_same_carrier': True,
            'hidden_spectator_selection_or_projective_occupation_read': False,
            'completed_quadrature_in_bare_seam_algebra': True,
            'intermediate_source_pulse_locality_established': False}

