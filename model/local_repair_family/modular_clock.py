"""One product-reference clock through the complete biased repair instrument.

For every recorded partial injection K, sigma(output)/sigma(input) is a
constant r_K. Thus U_s K = r_K**(-is) K U_s, U_s=sigma**(-is). Scalars
multiply along arbitrary mixed C/E histories and disappear in each CP block,
also with an untouched reference. The same word therefore works before,
between, or after repairs. This is covariance, not preservation of probe
visibility. An overlapping E repair below destroys the chosen coherence.

The source phase and balanced probe/analyzer remain declared coherent
services. The exact full positive-seam scheduler response is separate from
the explicitly requested two-seam source execution; no selected history is
mistaken for the full scheduler. Preparation bounds concern its reduced
source/reference output, not an independent full preparation archive.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
import math

import numpy as np

from . import local_bias as law
from ..native_modular_reference import common_caps as common, cap_local_clock as cap
from ..native_modular_reference import local_history
from ..observer.source_discrimination import circuit

PROBE = (1, 7)  # Existing full bare-cap source mixer, all ten spectators.
REQUESTED_SEAMS = ((0, 5), (5, 9))
WORD_GATE = F(1, 10**9)
ROTATION_GATE = F(1, 10**8)
NATIVE_GATE = F(1, 10**12)


def parameters():
    return law.Parameters((25, 25)+(8,)*10, (1,)*30, F(1, 2), F(1))


def _parameters(value):
    if type(value) is not law.Parameters:
        raise ValueError('original bounded family Parameters required')
    fractions = []
    for x in (value.c_weight, value.repair_probability):
        if (type(x) is not F or type(x.numerator) is not int or type(x.denominator) is not int
                or x.denominator <= 0 or math.gcd(x.numerator, x.denominator) != 1):
            raise ValueError('canonical original rational family probabilities required')
        fractions.append(F(x.numerator, x.denominator))
    return law.Parameters(value.port_weights, value.seam_weights, *fractions)


def _pairs():
    states = law.basis(); index = {x: i for i, x in enumerate(states)}
    a, b = PROBE
    return tuple((index[tuple(int(p in (a, k)) for p in range(12))],
                  index[tuple(int(p in (b, k)) for p in range(12))])
                 for k in range(12) if k not in PROBE)


def exact_certificate(value=None):
    p = _parameters(parameters() if value is None else value)
    sigma = law.reference(p); checked = 0
    for name in ('C', 'E'):
        for edge in law.SEAMS:
            for coin, label, pairs in law.raw_maps(name, edge):
                ratios = {sigma[j]/sigma[i] for i, j in pairs}
                current = law.records.record_current(name, (*edge, coin, label))
                expected = F(p.port_weights[edge[1]], p.port_weights[edge[0]])**current
                if ratios and ratios != {expected}:
                    raise ArithmeticError('one complete record has incompatible modular frequencies')
                checked += len(pairs)
    states = law.basis(); index = {x: i for i, x in enumerate(states)}
    for a, b in combinations(range(12), 2):
        ratios = {sigma[index[tuple(int(p in (a, k)) for p in range(12))]] /
                  sigma[index[tuple(int(p in (b, k)) for p in range(12))]]
                  for k in range(12) if k not in (a, b)}
        if ratios != {F(p.port_weights[a], p.port_weights[b])}:
            raise ArithmeticError('a spectator changes the bare regional modular action')
    return {'parameters': p.description(), 'complete_native_columns_checked': checked,
            'all_66_bare_two_port_algebras_invariant': True,
            'every_number_preserving_regional_algebra_invariant': True,
            'compiled_global_clock_is_region_local': False,
            'generator_modulo_scalar': '-sum_i log(theta_i) n_i',
            'complete_history_identity': 'U_s K_h = product_b(r_b)^(-is) K_h U_s',
            'arbitrary_untouched_reference': True,
            'mixed_law_tags_and_record_only_stopping_allowed': True,
            'branch_probabilities_or_quantum_outputs_normalized': False,
            'adaptive_stationarity_claimed': False,
            'coherent_phase_actuator_selected_by_bare_repair': False}


def _ratios():
    result = [F(1)]*512
    p = parameters()
    for address, state in zip(common.full.support()[0], law.basis(), strict=True):
        result[address] = F(math.prod(w**n for w, n in zip(p.port_weights, state, strict=True)))
    return tuple(result)


def source_program():
    return {'probe': cap.rotation_word(0, 1), 'analyzer': cap.rotation_word(0, -1),
            'half_clock': common.diagonal_word(F(1, 2), _ratios()),
            'clock': common.diagonal_word(F(1), _ratios())}


def _program_snapshot():
    """Detach and validate the literal exported words before execution or cost."""
    program = source_program()
    names = ('probe', 'analyzer', 'half_clock', 'clock')
    if type(program) is not dict or set(program) != set(names):
        raise ValueError('the four original source program words required')
    result = []
    for name in names:
        word = program[name]
        if type(word) is not list:
            raise ValueError('an original chronological source pulse list required')
        detached = []
        for pulse in word:
            if type(pulse) is not dict or set(pulse) != {'kind', 'angle', 'targets'}:
                raise ValueError('the original source pulse schema required')
            kind, angle, targets = pulse['kind'], pulse['angle'], pulse['targets']
            if (type(kind) is not str or kind not in common.full.compiler.KINDS
                    or type(angle) is not float or not math.isfinite(angle)
                    or type(targets) is not list or len(targets) != (2 if kind == 'J' else 1)
                    or any(type(bit) is not int or not 0 <= bit < 9 for bit in targets)
                    or len(set(targets)) != len(targets)):
                raise ValueError('finite original source pulses on distinct nine-bit targets required')
            detached.append((kind, angle, tuple(targets)))
        result.append((name, tuple(detached)))
    return tuple(result)


def _program(snapshot):
    return {name: [{'kind': kind, 'angle': angle, 'targets': list(targets)}
                   for kind, angle, targets in word] for name, word in snapshot}


@lru_cache(maxsize=8)
def _executed_word(name, word):
    """Execute exactly these literal pulses on every original address column.

    Bounded chunks avoid long propagation of floating subnormals. The chunks
    are slices of the exported word, not independently regenerated factors;
    their unmodified arrays are multiplied in chronological order.
    """
    tape = _program(((name, word),))[name]
    array = np.eye(512, dtype=complex)
    for start in range(0, len(tape), 2048):
        array = common.full.compiler.execute(tape[start:start+2048], 9)@array
    if name in ('half_clock', 'clock'):
        parameter = F(1, 2) if name == 'half_clock' else F(1)
        if common.certify_diagonal(parameter, _ratios(), array) > WORD_GATE:
            raise ArithmeticError('literal full source phase array exceeds its rational certificate')
    elif cap.certify_rotation(0, 1 if name == 'probe' else -1, array) > ROTATION_GATE:
        raise ArithmeticError('literal full bare-region mixer exceeds its rational certificate')
    return name, array.dtype.str, array.tobytes()


@lru_cache(maxsize=1)
def _source_data(snapshot):
    return tuple(_executed_word(name, word) for name, word in snapshot)


def source_arrays():
    """Fresh copies of original executed 512-square arrays, including padding."""
    return {name: np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512).copy()
            for name, dtype, data in _source_data(_program_snapshot())}


def _ideal_arrays():
    return {'probe': cap.ideal_rotation(0, 1), 'analyzer': cap.ideal_rotation(0, -1),
            'half_clock': np.diag([complex(float(v))**(-.5j) for v in _ratios()]),
            'clock': np.diag([complex(float(v))**(-1j) for v in _ratios()])}


def native_branches(name, edge, value=None):
    """Every actual raw slot of the new biased law, including zero maps."""
    name = law.records.candidate(name); edge = law.records.seam(edge)
    p = _parameters(parameters() if value is None else value)
    return tuple((label, pairs, root, F(n, d))
                 for label, pairs, root, n, d in _native_rows(name, edge, p.port_weights))


@lru_cache(maxsize=256)
def _native_rows(name, edge, weights):
    p = law.Parameters(weights, (1,)*30, F(1, 2), F(1))
    maps = {(c, label): pairs for c, label, pairs in law.raw_maps(name, edge)}
    result = []
    for coin in (0, 1):
        probability = law.coin_probability(p, edge, coin)
        root = math.sqrt(float(probability))
        lo, hi = common.full.exact._sqrt_interval(probability)
        if max(abs(F.from_float(root)-lo), abs(F.from_float(root)-hi))*2 > NATIVE_GATE:
            raise ArithmeticError('biased native square root exceeds its frame allowance')
        decoded = circuit.record_values(name, coin)
        for raw in range(1 << circuit.RECORD_BITS[name]):
            label = decoded.get(raw)
            pairs = maps.get((coin, label), ())
            current = law.records.record_current(name, (*edge, coin, label)) if label is not None else None
            result.append(((coin, raw, label, current), pairs, root,
                           probability.numerator, probability.denominator))
    return tuple(result)


def _native(density, pairs, probability):
    r = len(density)//66
    result = np.zeros_like(density, dtype=complex)
    source = [r*i+a for i, _ in pairs for a in range(r)]
    target = [r*j+a for _, j in pairs for a in range(r)]
    result[np.ix_(target, target)] = float(probability)*density[np.ix_(source, source)]
    return result


def _validate_input(density):
    # This numerical input guard does not decide exact rank or remove tiny
    # positive populations. Component certificates give arbitrary-R bounds;
    # the executable density interface is explicitly R dimension one or two.
    if type(density) is np.ndarray and density.ndim == 2 and np.any(np.diag(density).real < 0):
        raise ValueError('negative original input population')
    return local_history._density(density).astype(complex, copy=True)


def _export(array, density):
    """The same complete address export, contracting the traced padding first."""
    r = len(density)//66
    codes = np.asarray(common.full.support()[0]); outside = np.setdiff1d(np.arange(512), codes)
    block, escaped = array[np.ix_(codes, codes)], array[np.ix_(outside, codes)]
    k = np.kron(block, np.eye(r))
    # Trace only the escaped address. Off-diagonal untouched-reference entries
    # survive this exact tensor contraction; no padding amplitude is dropped.
    abort = np.einsum('ji,irjs->rs', escaped.conj().T@escaped,
                      density.reshape(66, r, 66, r), optimize=True)
    return k@density@k.conj().T, abort


def execute(density, name, *, order='interleaved', raw=True, seams=REQUESTED_SEAMS):
    """Complete two-requested-seam outputs, with lazy exact-zero raw slots.

    Nonzero native maps are retained even when their particular input mass
    is zero. ``raw_output`` supplies every structurally zero raw history.
    Export aborts retain the untouched reference, not a normalized survivor.
    This declared schedule is not the full positive-seam scheduler.
    """
    density = _validate_input(density); name = law.records.candidate(name)
    if type(raw) is not bool or type(order) is not str or order not in ('interleaved', 'endpoint', 'zero'):
        raise ValueError('declared clock order and original Boolean source selector required')
    if type(seams) is not tuple or len(seams) != 2:
        raise ValueError('exactly two requested native seams required')
    seams = tuple(law.records.seam(e) for e in seams)
    arrays = source_arrays() if raw else _ideal_arrays()
    output = {}; r = len(density)//66
    rho, output['abort', 'probe'] = _export(arrays['probe'], density)
    frames = {(): rho}
    for step, edge in enumerate(seams, 1):
        following = {}
        for history, state in frames.items():
            for label, pairs, root, probability in native_branches(name, edge):
                if not pairs:
                    continue
                key = (*history, label)
                state2 = _native(state, pairs, root*root if raw else probability)
                if order == 'interleaved' or order == 'endpoint' and step == 2:
                    word = 'half_clock' if order == 'interleaved' else 'clock'
                    state2, abort = _export(arrays[word], state2)
                else:
                    abort = np.zeros((r, r), complex)
                output['abort', 'clock', step, key] = abort
                following[key] = state2
        frames = following
    for history, state in frames.items():
        state, output['abort', 'analyzer', history] = _export(arrays['analyzer'], state)
        for label, pairs, root, probability in native_branches(name, PROBE):
            if pairs:
                output['read', history, label] = _native(state, pairs, root*root if raw else probability)
    return output


def raw_output(outputs, name, slots):
    """Read a complete three-record raw word, refusing a missing live block."""
    name = law.records.candidate(name)
    if type(slots) is not tuple or len(slots) != 3:
        raise ValueError('three raw (coin, register) slots required')
    labels, zero = [], False
    for slot, edge in zip(slots, (*REQUESTED_SEAMS, PROBE), strict=True):
        if (type(slot) is not tuple or len(slot) != 2 or any(type(i) is not int for i in slot)
                or slot[0] not in (0, 1) or not 0 <= slot[1] < 1 << circuit.RECORD_BITS[name]):
            raise ValueError('original integer coin/raw slots required')
        found = next((row for row in native_branches(name, edge) if row[0][:2] == slot), None)
        zero |= not found[1]
        labels.append(found[0])
    if zero:
        r = len(outputs['abort', 'probe'])
        return np.zeros((66*r, 66*r), complex)
    return outputs['read', tuple(labels[:2]), labels[2]].copy()


def scheduler_response(value=None, attempts=2):
    """Exact full fixed-family mean current A+B*cos(s log(theta1/theta7)).

    The diagonal populations follow the actual heat-bath kernel. Off-diagonal
    pairs retain their difference vector under every common native record;
    their ten spectator coefficients follow the complete killed pair kernel.
    All law/seam/coin labels and holds are summed, with no conditioning.
    """
    p = _parameters(parameters() if value is None else value)
    if type(attempts) is not int or not 0 <= attempts <= 16:
        raise ValueError('bounded integer scheduler attempts0..16 required')
    sigma = law.reference(p); pairs = _pairs(); pair_index = {v: i for i, v in enumerate(pairs)}
    diagonal = list(sigma); coherence = []
    for i, j in pairs:
        diagonal[i] = diagonal[j] = (sigma[i]+sigma[j])/2
        coherence.append((sigma[j]-sigma[i])/2)
    transfer = [[F() for _ in pairs] for _ in pairs]
    for i in range(10):
        transfer[i][i] = 1-p.repair_probability
    for name, mixture in (('C', p.c_weight), ('E', 1-p.c_weight)):
        for edge, weight in zip(law.SEAMS, p.seam_weights, strict=True):
            choice = p.repair_probability*mixture*F(weight, sum(p.seam_weights))
            for coin, _, mapping in law.raw_maps(name, edge):
                mapping = dict(mapping); probability = choice*law.coin_probability(p, edge, coin)
                for k, (i, j) in enumerate(pairs):
                    if i in mapping and j in mapping:
                        target = (mapping[i], mapping[j])
                        if target not in pair_index:
                            raise ArithmeticError('a complete record changed the conserved pair difference')
                        transfer[k][pair_index[target]] += probability
    population = law.kernel(p)
    for _ in range(attempts):
        diagonal = [sum((diagonal[i]*population[i][j] for i in range(66)), F()) for j in range(66)]
        coherence = [sum((coherence[i]*transfer[i][j] for i in range(10)), F()) for j in range(10)]
    one_load = sum((diagonal[i]+diagonal[j] for i, j in pairs), F())
    plus = law.coin_probability(p, PROBE, 1); minus = law.coin_probability(p, PROBE, 0)
    constant, cosine = (minus-plus)*one_load/2, -sum(coherence, F())
    if p.port_weights[PROBE[0]] > p.port_weights[PROBE[1]] and cosine <= 0:
        raise ArithmeticError('full positive-seam scheduler lost every phase-sensitive path')
    return {'attempts': attempts, 'constant': str(constant), 'cosine': str(cosine),
            'frequency_ratio': str(F(p.port_weights[PROBE[0]], p.port_weights[PROBE[1]])),
            'full_seam_scheduler': True, 'C_weight': str(p.c_weight),
            'no_postselection': True, 'all_held_and_law_tagged_events_included': True}


def _unread(density, p):
    result = math.sqrt(float(1-p.repair_probability))**2*density
    for name, mixture in (('C', p.c_weight), ('E', 1-p.c_weight)):
        for edge, weight in zip(law.SEAMS, p.seam_weights, strict=True):
            choice = p.repair_probability*mixture*F(weight, sum(p.seam_weights))
            for _, pairs, root, _ in native_branches(name, edge, p):
                amplitude = math.sqrt(float(choice))*root
                result = result+_native(density, pairs, amplitude*amplitude)
    return result


def scheduler_frame_error(value):
    """Rational full flagged-frame bound on the actual scheduler coefficients.

    For each choice/coin exactly one raw record is attained per input column.
    Their environment labels are orthogonal. Consequently the squared frame
    error is the sum below, independent of the column and reference.
    """
    p = _parameters(value); squared = F()
    choices = [(1-p.repair_probability, math.sqrt(float(1-p.repair_probability)))]
    for name, mixture in (('C', p.c_weight), ('E', 1-p.c_weight)):
        for edge, weight in zip(law.SEAMS, p.seam_weights, strict=True):
            choice = p.repair_probability*mixture*F(weight, sum(p.seam_weights))
            for coin in (0, 1):
                probability = law.coin_probability(p, edge, coin)
                choices.append((choice*probability,
                                math.sqrt(float(choice))*math.sqrt(float(probability))))
    for probability, amplitude in choices:
        if not probability:
            if amplitude != 0:
                raise ArithmeticError('a zero scheduler choice gained amplitude')
            continue
        lo, hi = common.full.exact._sqrt_interval(probability)
        error = max(abs(F.from_float(amplitude)-lo), abs(F.from_float(amplitude)-hi))
        squared += error*error
    result = common.full.exact._sqrt_interval(squared)[1]
    if result > NATIVE_GATE:
        raise ArithmeticError('complete biased scheduler frame exceeds its declared gate')
    return F(result.numerator, result.denominator)


def detector_frame_error(value):
    p = _parameters(value); squared = F()
    for mixture in (p.c_weight, 1-p.c_weight):
        for coin in (0, 1):
            probability = law.coin_probability(p, PROBE, coin)
            exact = mixture*probability
            if not exact:
                continue
            actual = math.sqrt(float(mixture))*math.sqrt(float(probability))
            lo, hi = common.full.exact._sqrt_interval(exact)
            squared += max(abs(F.from_float(actual)-lo), abs(F.from_float(actual)-hi))**2
    result = common.full.exact._sqrt_interval(squared)[1]
    if result > NATIVE_GATE:
        raise ArithmeticError('complete biased terminal detector frame exceeds its gate')
    return F(result.numerator, result.denominator)


def scheduler_execute(value=None, *, raw=True, order='interleaved'):
    """Actual full-family marginal from raw source factors, both repair events.

    Only the displayed statistic discards the native tape. The component
    certificate bounds the complete taped instrument on arbitrary references.
    """
    p = _parameters(parameters() if value is None else value)
    if (p.port_weights != parameters().port_weights or type(raw) is not bool
            or type(order) is not str or order not in ('interleaved', 'endpoint', 'zero')):
        raise ValueError('canonical source reference, Boolean source flag, and clock order required')
    scheduler_frame_error(p)
    detector_frame_error(p)
    arrays = source_arrays() if raw else _ideal_arrays()
    rho = np.diag(list(map(float, law.reference(p)))).astype(complex)
    rho, abort = _export(arrays['probe'], rho); abort_mass = float(np.trace(abort).real)
    for step in (1, 2):
        rho = _unread(rho, p)
        if order == 'interleaved' or order == 'endpoint' and step == 2:
            rho, abort = _export(arrays['half_clock' if order == 'interleaved' else 'clock'], rho)
            abort_mass += float(np.trace(abort).real)
    rho, abort = _export(arrays['analyzer'], rho); abort_mass += float(np.trace(abort).real)
    expectation = mass = 0.
    for name, mixture in (('C', p.c_weight), ('E', 1-p.c_weight)):
        for label, pairs, root, _ in native_branches(name, PROBE, p):
            amplitude = math.sqrt(float(mixture))*root
            probability = float(np.trace(_native(rho, pairs, amplitude*amplitude)).real)
            mass += probability
            expectation += (label[3] or 0)*probability
    return {'mean_current': expectation, 'complete_mass': mass+abort_mass, 'export_abort_mass': abort_mass}


def _processor_failure_bound(pulses):
    pulse_gate = F(2, 10**11)
    if type(pulses) is not int or pulses < 0 or pulses*pulse_gate >= 1:
        raise ArithmeticError('processor pulse count exceeds the first-failure gain domain')
    return pulse_gate, pulses*pulse_gate**2*(1+NATIVE_GATE)**6/(1-pulses*pulse_gate)**2


def _processor_bounds(array, dimension):
    """Outward dyadic bounds for the original six-level pulse's two blocks.

    Frobenius norms bound the complement and A* A-I. The latter bounds
    max(0, ||A||-1) without rounding a singular value near one. Every stored
    real/imaginary component enters exact rational arithmetic first.
    """
    if type(array) is not np.ndarray or type(dimension) is not int or dimension not in (2, 4):
        raise ValueError('an original six-level pulse and proper code dimension required')
    array = common.full.exact._array(array, (6, 6))
    components = tuple(tuple(common.full.exact._component(value) for value in row[:dimension])
                       for row in array)
    complement_squared = sum((re*re+im*im for row in components[dimension:] for re, im in row), F())
    gram_squared = F()
    for i in range(dimension):
        for j in range(dimension):
            real, imag = -F(int(i == j)), F()
            for row in components[:dimension]:
                ar, ai = row[i]; br, bi = row[j]
                real += ar*br+ai*bi
                imag += ar*bi-ai*br
            gram_squared += real*real+imag*imag
    bounds = tuple(common.full.exact._sqrt_interval(value)[1]
                   for value in (complement_squared, gram_squared))
    return tuple(F(bound.numerator, bound.denominator) for bound in bounds)


def source_budget():
    snapshot = _program_snapshot()
    _source_data(snapshot)
    error, gain = F(), F(1)
    for eps in (ROTATION_GATE, NATIVE_GATE, WORD_GATE, NATIVE_GATE, WORD_GATE, ROTATION_GATE, NATIVE_GATE):
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    program = _program(snapshot)
    pulses = len(program['probe'])+max(2*len(program['half_clock']), len(program['clock']))+len(program['analyzer'])
    # The complete processor exporter has exactly zero complement on the
    # invariant source code. Keep the floating implementation's nonzero
    # complement separately, with first-failure labels, rather than silently
    # interpreting its compressed success product as a complete channel.
    pulse_gate, failure_mass = _processor_failure_bound(pulses)
    menu = {(e['kind'], e['angle'], len(e['targets'])) for word in program.values() for e in word}
    for kind, angle, bits in menu:
        actual = common.full.compiler.native_pulse(kind, angle); d = 1 << bits
        leakage, excess = _processor_bounds(actual, d)
        if max(leakage, excess) > pulse_gate:
            raise ArithmeticError('proper-code processor complement or gain exceeds its fixed gate')
    # Sum of P first-failure CP masses. Their source-code prefixes have norm
    # <=(1+g)^P<=1/(1-Pg); the three native frames add their independent gain.
    return {'complete_instrument_halfdiamond_upper': str(error+failure_mass/2),
            'successful_processor_path_halfdiamond_upper': str(error),
            'proper_code_first_failure_mass_upper': str(failure_mass),
            'proper_code_per_pulse_complement_and_gain_gate': str(pulse_gate),
            'actual_unique_processor_pulses_checked': len(menu),
            'channel_norm_upper': str(gain+failure_mass),
            'maximum_source_pulses': pulses, 'coherent_load_export_calls': 2*pulses,
            'native_continuation_attempts': 2, 'native_detector_attempts': 1,
            'quantum_address_qubits': 9, 'original_source_dimension': 512,
            'source_program_bytes': len(json.dumps(program, sort_keys=True, separators=(',', ':')).encode()),
            'source_position_counter_bits': max(map(len, program.values())).bit_length(),
            'source_position_counter_scope': 'live within-word interpreter only',
            'coherent_word_stage_bits': 2,
            'first_failed_load_or_export_ordinal_bits': (2*pulses).bit_length(),
            'processor_complement_environment_dimension_upper': 6,
            'processor_complement_index_is_accessible_record': False,
            'address_export_abort_position_and_status_bits': 3,
            'processor_failure_outputs': 'code-or-failure flag and first-failure ordinal retained; complement Kraus index traced as environment; exact zero on invariant proper codes, numerical aggregate bounded above',
            'requested_history_raw_bits_upper': 3*(1+3),
            'full_scheduler_history_raw_bits_upper': 2*(1+1+5+1+3)+(1+1+3),
            'source_export_abort_flags_retained': True,
            'physical_pulse_bound': 'operator error delta=P*epsilon/(1-P*epsilon), P*epsilon<1; CQ halfdiamond<=delta*(2+delta)/2; transfer/native errors additional',
            'hardware_and_physical_calibration_are_external_obligations': True,
            'supplied_service_assumptions': ['coherent proper-code load/export', 'source pulse actuator',
                                            'classical selector and biased placement coin',
                                            'record storage/export', 'terminal native read'],
            'numerical_certificates_are_not_physical_precision_claims': True}


def physical_allowance(pulse_error=F(), transfer_error=F(), native_attempt_error=F()):
    """Conditional exact bound from supplied full-isometry operator errors."""
    values = []
    for x in (pulse_error, transfer_error, native_attempt_error):
        if (type(x) is not F or type(x.numerator) is not int or type(x.denominator) is not int
                or x.denominator <= 0 or math.gcd(x.numerator, x.denominator) != 1 or x.numerator < 0):
            raise ValueError('canonical nonnegative original rational error required')
        values.append(F(x.numerator, x.denominator))
    snapshot = _program_snapshot()
    _source_data(snapshot)
    program = _program(snapshot)
    pulses = len(program['probe'])+max(2*len(program['half_clock']), len(program['clock']))+len(program['analyzer'])
    total = pulses*values[0]+2*pulses*values[1]+3*values[2]
    if total >= 1:
        raise ValueError('sum of supplied operator error allowances must be below one')
    delta = total/(1-total)
    return {'operator_upper': str(delta), 'halfdiamond_upper': str(delta*(2+delta)/2),
            'supplied_physical_errors_not_measured_here': True}


def source_sha256():
    paths = set(common.SOURCE_PATHS) | set(cap.interface.SOURCE_PATHS) | set(cap.local.SOURCE_PATHS)
    paths.update({'model/local_repair_family/local_bias.py',
                  'model/local_repair_family/modular_clock.py',
                  'model/local_repair_family/test_modular_clock.py',
                  'model/native_modular_reference/cap_local_clock.py',
                  'model/native_modular_reference/local_history.py',
                  'model/observer/source_discrimination/circuit.py',
                  'model/observer/record_recovery.py'})
    return law.source_bindings(sorted(paths))


def report():
    certificate = exact_certificate(); budget = source_budget()
    response = {}
    for label, mixture in (('C', F(1)), ('E', F()), ('mixed', F(1, 2))):
        p = law.Parameters(parameters().port_weights, (1,)*30, mixture, F(1))
        exact = scheduler_response(p); actual = scheduler_execute(p)
        expected = float(F(exact['constant']))+float(F(exact['cosine']))*math.cos(math.log(25/8))
        allowance = float(F(budget['complete_instrument_halfdiamond_upper']))
        if (not all(math.isfinite(actual[k]) for k in ('mean_current', 'complete_mass', 'export_abort_mass'))
                or abs(actual['mean_current']-expected) > 2*allowance or abs(actual['complete_mass']-1) > 2*allowance):
            raise ArithmeticError('source-executed full-family readout exceeds the complete bound')
        re, _, radius = common.phase_enclosure(F(1), F(25, 8))
        signal = tuple(F(exact['cosine'])*(v-1) for v in (re-radius, re+radius))
        scale = 10**9
        enclosure = [str(F(math.floor(signal[0]*scale), scale)),
                     str(F(math.ceil(signal[1]*scale), scale))]
        log_lo, log_hi = common.logarithm_interval(F(25, 8))
        x = F(113, 200)
        derivative = F(exact['cosine'])*F(113, 100)*(x-x**3/6)
        if not F(113, 100) < log_lo <= log_hi < F(115, 100) or derivative <= F(1, 32):
            raise ArithmeticError('full-scheduler mean clock lost its finite inverse bound')
        response[label] = {**exact, 'source_readout_within_complete_bound': True,
                           'mean_current_s1_minus_s0_interval': enclosure,
                           'mean_parameter_interval': ['1/2', '1'],
                           'inverse_mean_Lipschitz_upper': '32',
                           'individual_random_record_claimed_monotone': False}
    requested = {}
    density = np.diag(list(map(float, law.reference(parameters())))).astype(complex)
    allowance = float(F(budget['complete_instrument_halfdiamond_upper']))
    for name in ('C', 'E'):
        interleaved = execute(density, name)
        endpoint = execute(density, name, order='endpoint')
        target = execute(density, name, raw=False)
        if interleaved.keys() != endpoint.keys() or interleaved.keys() != target.keys():
            raise ArithmeticError('clock placement changed the complete output alphabet')
        source_error = sum(np.linalg.norm(interleaved[k]-target[k], 'nuc') for k in target)/2
        ordering_error = sum(np.linalg.norm(interleaved[k]-endpoint[k], 'nuc') for k in target)/2
        mass = sum(float(np.trace(x).real) for x in interleaved.values())
        if (not all(math.isfinite(v) for v in (source_error, ordering_error, mass))
                or source_error > allowance or ordering_error > 2*allowance or abs(mass-1) > 2*allowance):
            raise ArithmeticError('complete continued source clock violates its channel allowance')
        requested[name] = {'materialized_live_and_export_blocks': len(interleaved),
                           'three_attempt_full_raw_words': (2 << circuit.RECORD_BITS[name])**3,
                           'structural_zero_words_have_reference_bearing_zero_outputs': True,
                           'source_vs_ideal_halftrace_upper': str(F(budget['complete_instrument_halfdiamond_upper'])),
                           'interleaved_vs_endpoint_halftrace_upper': str(2*F(budget['complete_instrument_halfdiamond_upper'])),
                           'all_complete_output_gates_passed': True}
    return {'schema': 'oph.local_repair_family.modular_clock.v1', 'theorem': certificate,
            'requested_source_schedule': [list(e) for e in REQUESTED_SEAMS],
            'requested_schedule_is_full_scheduler': False,
            'full_scheduler_response': response, 'source_budget': budget,
            'requested_complete_source_execution': requested,
            'probe_region': list(PROBE), 'parameter': 'two half-steps, total s=1',
            'fresh_reference_reset_between_events': False,
            'physical_time_or_integrated_gate_claimed': False, 'source_sha256': source_sha256()}
