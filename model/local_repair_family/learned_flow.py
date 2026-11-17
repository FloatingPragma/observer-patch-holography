"""Compile the product-reference flow from an inferred port vector alone.

On the binary two-load carrier, a recorded native partial injection changes
the product weight by (w_b/w_a)**current, independently of the spectator.
Consequently its complete CP block intertwines the inferred modular flow,
even when these inferred weights differ from the apparatus's true weights.
This covariance does not assert stationarity for a misspecified reference.

The phase compiler receives no apparatus parameter or prepared state. The
separate executable experiment uses the existing biased native instrument,
two declared repair seams, and a terminal native read. All source-address
export failures and raw native slots remain in the output. Numerical source
certificates concern the supplied coherent service, not its physical fidelity
or the statistical accuracy of the inferred vector.
"""
from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
import json
import math

import numpy as np

from . import modular_clock as clock

law, common, cap = clock.law, clock.common, clock.cap
PROBE, REQUESTED_SEAMS = clock.PROBE, clock.REQUESTED_SEAMS


def _weights(values):
    """Snapshot original rational slots before comparison or cache lookup."""
    if type(values) is not tuple or len(values) != 12:
        raise ValueError('twelve original positive rational relative weights required')
    result = []
    for x in values:
        if type(x) is int:
            n, d = x, 1
        elif type(x) is F:
            n, d = x.numerator, x.denominator
        else:
            raise ValueError('original integer or canonical Fraction weights required')
        if (type(n) is not int or type(d) is not int or n <= 0 or d <= 0
                or math.gcd(n, d) != 1):
            raise ValueError('canonical strictly positive rational weights required')
        result.append(F(n, d))
    return tuple(x/result[0] for x in result)


def _key(weights):
    return tuple((x.numerator, x.denominator) for x in _weights(weights))


def _source_chart():
    """Bind cached charts to the original native state order and encoding.

    The inherited regional API exposes a cached dataclass. Its state tuple
    is therefore checked afresh, before any program or ideal-array cache is
    trusted. The fixed probe's targets are reconstructed independently.
    """
    states = law.native.states(binary=True)
    codes = tuple(clock.circuit.encode(state) for state in states)
    cached_codes = common.full.support()[0]
    cached_states = common.regional.source()[0].states
    if (law.basis() != states or type(cached_codes) is not tuple
            or any(type(x) is not int for x in cached_codes) or cached_codes != codes
            or type(cached_states) is not tuple
            or any(type(x) is not tuple or any(type(n) is not int for n in x) for x in cached_states)
            or cached_states != states):
        raise ArithmeticError('cached source chart differs from the committed native state order or encoding')
    index = {state: i for i, state in enumerate(states)}
    a, b = PROBE; pairs = []
    for i, state in enumerate(states):
        if state[a] == 1 and state[b] == 0:
            moved = list(state); moved[a], moved[b] = 0, 1
            j = index[tuple(moved)]
            pairs.append((i, j, codes[i], codes[j]))
    pairs = tuple(pairs)
    if cap.pairs(0) != pairs:
        raise ArithmeticError('inherited probe pairs differ from the committed native probe')
    return states, codes, pairs


def ratios(relative_weights):
    """Fresh full-address diagonal, identity on the 446 padding addresses.

    Dividing all port weights by w_0 removes only a scalar on the fixed-load
    source sector. Padding remains an explicitly specified identity extension.
    """
    weights = _weights(relative_weights)
    states, codes, _ = _source_chart()
    result = [F(1) for _ in range(512)]
    for address, state in zip(codes, states, strict=True):
        result[address] = math.prod(w**n for w, n in zip(weights, state, strict=True))
    return tuple(result)


def exact_certificate(relative_weights):
    weights = _weights(relative_weights)
    states = law.basis(); index = {x: i for i, x in enumerate(states)}
    products = tuple(math.prod(w**n for w, n in zip(weights, x, strict=True)) for x in states)
    checked = 0
    for name in ('C', 'E'):
        for a, b in law.SEAMS:
            for coin, label, pairs in law.raw_maps(name, (a, b)):
                current = law.records.record_current(name, (a, b, coin, label))
                factor = (weights[b]/weights[a])**current
                if any(products[j] != factor*products[i] for i, j in pairs):
                    raise ArithmeticError('the native branch does not intertwine the inferred flow')
                checked += len(pairs)
    for a in range(12):
        for b in range(a+1, 12):
            for k in range(12):
                if k in (a, b):
                    continue
                i = index[tuple(int(p in (a, k)) for p in range(12))]
                j = index[tuple(int(p in (b, k)) for p in range(12))]
                if products[i]/products[j] != weights[a]/weights[b]:
                    raise ArithmeticError('a spectator changed the inferred regional flow')
    return {'normalized_relative_weights': list(map(str, weights)),
            'complete_native_columns_checked': checked,
            'all_66_bare_two_port_algebras_invariant': True,
            'global_compiled_word_claimed_physically_region_local': False,
            'generator_modulo_scalar': '-sum_i log(inferred_weight_i) n_i',
            'complete_history_identity': 'U_s K_h = product_b(r_b)^(-is) K_h U_s',
            'mixed_law_tags_and_record_only_stopping_covered': True,
            'arbitrary_untouched_reference': True,
            'stationarity_of_misspecified_reference_claimed': False,
            'true_weights_are_compiler_inputs': False,
            'coherent_control_is_supplied_service': True}


@lru_cache(maxsize=4)
def _program_data(weight_key):
    diagonal = ratios(tuple(F(n, d) for n, d in weight_key))
    return _validate_program({'probe': cap.rotation_word(0, 1), 'analyzer': cap.rotation_word(0, -1),
                              'half_clock': common.diagonal_word(F(1, 2), diagonal),
                              'clock': common.diagonal_word(F(1), diagonal)})


def source_program(relative_weights):
    _source_chart()
    return clock._program(_program_data(_key(relative_weights)))


def _validate_program(program):
    names = ('probe', 'analyzer', 'half_clock', 'clock')
    if type(program) is not dict or set(program) != set(names):
        raise ValueError('the four original learned source program words required')
    result = []
    for name in names:
        if type(program[name]) is not list:
            raise ValueError('original chronological pulse lists required')
        word = []
        for pulse in program[name]:
            if type(pulse) is not dict or set(pulse) != {'kind', 'angle', 'targets'}:
                raise ValueError('original source pulse schema required')
            kind, angle, targets = pulse['kind'], pulse['angle'], pulse['targets']
            if (type(kind) is not str or kind not in common.full.compiler.KINDS
                    or type(angle) is not float or not math.isfinite(angle)
                    or type(targets) is not list or len(targets) != (2 if kind == 'J' else 1)
                    or any(type(bit) is not int or not 0 <= bit < 9 for bit in targets)
                    or len(set(targets)) != len(targets)):
                raise ValueError('finite original source pulses on distinct nine-bit targets required')
            word.append((kind, angle, tuple(targets)))
        result.append((name, tuple(word)))
    return tuple(result)


def _snapshot(relative_weights):
    # The public literal program is still reread before every cache lookup.
    _source_chart()
    return _validate_program(source_program(relative_weights))


@lru_cache(maxsize=512)
def _logarithm_payload(numerator, denominator):
    """Independent exact log authority; no inherited mutable Fraction cache."""
    ratio = F(numerator, denominator)
    exponent = numerator.bit_length()-denominator.bit_length()
    reduced = ratio/F(2)**exponent
    if reduced < 1:
        reduced *= 2; exponent -= 1
    def local(value):
        x = (value-1)/(value+1)
        lower = 2*sum((x**(2*j+1)/(2*j+1) for j in range(32)), F())
        tail = 2*x**65/(65*(1-x*x))
        return lower, lower+tail
    lo, hi = local(reduced)
    shift = sorted(exponent*x for x in local(F(2)))
    return tuple((x.numerator, x.denominator) for x in (lo+shift[0], hi+shift[1]))


def _sqrt_interval(value):
    """Fresh outward integer-isqrt endpoints, including subnormal scales."""
    n, d = value.numerator, value.denominator
    rn, rd = math.isqrt(n), math.isqrt(d)
    if rn*rn == n and rd*rd == d:
        return F(rn, rd), F(rn, rd)
    exponent = (n.bit_length()-d.bit_length())//2
    scale = 1 << max(0, 192-exponent)
    lower = math.isqrt(n*scale*scale//d)
    return F(lower, scale), F(lower+1, scale)


@lru_cache(maxsize=512)
def _phase_payload(sn, sd, numerator, denominator):
    """Exact enclosure with bounded dyadic centers, not giant log powers.

    exp(i*x) is one-Lipschitz on the real line. The degree40 Taylor
    remainder is at most |x|**41/41! by its integral formula. Each squaring
    changes radius e to e*(2+e), since the exact phase has modulus one.
    Every dyadic rounding is charged by its exact componentwise L1 error.
    """
    parameter = F(sn, sd)
    lo, hi = (F(n, d) for n, d in _logarithm_payload(numerator, denominator))
    theta, halvings = -parameter*(lo+hi)/2, 0
    while abs(theta) > 1:
        theta /= 2
        halvings += 1
    scale = 1 << 128
    center = F(math.floor(theta*scale), scale)
    radius = abs(parameter)*(hi-lo)/(2*2**halvings)+abs(theta-center)
    real, imag, pr, pi = F(1), F(), F(1), F()
    for k in range(1, 41):
        pr, pi = -center*pi/k, center*pr/k
        real += pr; imag += pi
    radius += abs(center)**41/math.factorial(41)
    for step in range(halvings+1):
        re, im = F(math.floor(real*scale), scale), F(math.floor(imag*scale), scale)
        radius += abs(real-re)+abs(imag-im)
        real, imag = re, im
        # An outward dyadic radius also prevents its denominator from growing
        # exponentially during repeated squaring. No uncertainty is removed.
        radius = F(math.ceil(radius*scale), scale)
        if step < halvings:
            real, imag = real*real-imag*imag, 2*real*imag
            radius *= 2+radius
    return tuple((x.numerator, x.denominator) for x in (real, imag, radius))


def _phase_enclosure(parameter, ratio):
    return tuple(F(n, d) for n, d in _phase_payload(parameter.numerator, parameter.denominator,
                                                  ratio.numerator, ratio.denominator))


def certify_phase(parameter, relative_weights, array):
    """Exact original-entry full512 Frobenius majorant, padding included."""
    if (type(parameter) is not F or type(parameter.numerator) is not int
            or type(parameter.denominator) is not int or parameter.denominator <= 0
            or math.gcd(parameter.numerator, parameter.denominator) != 1
            or parameter not in (F(1, 2), F(1))):
        raise ValueError('the original exact half-step or full-step parameter required')
    diagonal = ratios(relative_weights)
    if type(array) is not np.ndarray:
        raise ValueError('an original full512 source array required')
    array = common.full.exact._array(array, (512, 512))
    phases = {ratio: _phase_enclosure(parameter, ratio) for ratio in set(diagonal)}
    squared = F()
    for i, row in enumerate(array):
        re, im, _ = phases[diagonal[i]]
        for j, value in enumerate(row):
            ar, ai = common.full.exact._component(value)
            if i == j:
                ar, ai = ar-re, ai-im
            squared += ar*ar+ai*ai
    radius = max(x[2] for x in phases.values())
    return _sqrt_interval(squared)[1]+_sqrt_interval(F(512))[1]*radius


@lru_cache(maxsize=4)
def _checked_rotation(name, word):
    _, _, pairs = _source_chart()
    payload = clock._executed_word(name, word)
    _, dtype, data = payload
    array = np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)
    sign = 1 if name == 'probe' else -1
    lo, hi = _sqrt_interval(F(1, 2)); center, radius = (lo+hi)/2, (hi-lo)/2
    special = {}
    for _, _, a, b in pairs:
        special.update({(a, a): center, (b, b): center,
                        (a, b): sign*center, (b, a): -sign*center})
    squared = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            re, im = common.full.exact._component(value)
            squared += (re-special.get((i, j), F(int(i == j))))**2+im*im
    if _sqrt_interval(squared)[1]+_sqrt_interval(F(40))[1]*radius > clock.ROTATION_GATE:
        raise ArithmeticError('independent learned probe/analyzer certificate failed')
    return payload


def _processor_bounds(array, dimension):
    if type(array) is not np.ndarray:
        raise ValueError('an original full six-level pulse array required')
    array = common.full.exact._array(array, (6, 6))
    entries = tuple(tuple(common.full.exact._component(x) for x in row[:dimension]) for row in array)
    leakage = sum((re*re+im*im for row in entries[dimension:] for re, im in row), F())
    gram = F()
    for i in range(dimension):
        for j in range(dimension):
            re, im = -F(int(i == j)), F()
            for row in entries[:dimension]:
                ar, ai = row[i]; br, bi = row[j]
                re += ar*br+ai*bi; im += ar*bi-ai*br
            gram += re*re+im*im
    return _sqrt_interval(leakage)[1], _sqrt_interval(gram)[1]


@lru_cache(maxsize=8)
def _executed(name, word, weight_key):
    if name in ('probe', 'analyzer'):
        return _checked_rotation(name, word)
    tape = clock._program(((name, word),))[name]
    array = np.eye(512, dtype=complex)
    for start in range(0, len(tape), 2048):
        array = common.full.compiler.execute(tape[start:start+2048], 9)@array
    parameter = F(1, 2) if name == 'half_clock' else F(1)
    if certify_phase(parameter, tuple(F(n, d) for n, d in weight_key), array) > clock.WORD_GATE:
        raise ArithmeticError('literal learned full source phase exceeds its fixed certificate')
    return name, array.dtype.str, array.tobytes()


def source_arrays(relative_weights):
    """Detached full512 arrays from the actual exported literal words."""
    weights = _weights(relative_weights); key = _key(weights)
    snapshot = _snapshot(weights)
    return {name: np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512).copy()
            for name, dtype, data in (_executed(name, word, key) for name, word in snapshot)}


@lru_cache(maxsize=4)
def _ideal_data(weight_key):
    diagonal = ratios(tuple(F(n, d) for n, d in weight_key))
    arrays = {'probe': cap.ideal_rotation(0, 1), 'analyzer': cap.ideal_rotation(0, -1)}
    for name, parameter in (('half_clock', F(1, 2)), ('clock', F(1))):
        phases = {}
        for ratio in set(diagonal):
            re, im, _ = _phase_enclosure(parameter, ratio)
            phases[ratio] = complex(float(re), float(im))
        arrays[name] = np.diag([phases[x] for x in diagonal])
    return tuple((name, array.dtype.str, array.tobytes()) for name, array in arrays.items())


def _ideal_arrays(relative_weights):
    _source_chart()
    return {name: np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512).copy()
            for name, dtype, data in _ideal_data(_key(relative_weights))}


def _native_branches(name, edge):
    rows = clock.native_branches(name, edge)
    for _, _, root, probability in rows:
        lo, hi = _sqrt_interval(probability)
        if (not math.isfinite(root)
                or 2*max(abs(F.from_float(root)-lo), abs(F.from_float(root)-hi)) > clock.NATIVE_GATE):
            raise ArithmeticError('independent native frame certificate failed')
    return rows


def execute(density, relative_weights, name, *, order='interleaved', raw=True):
    """Two actual biased repairs and one terminal read, no conditioning.

    The declared native apparatus is used only here, never by the compiler.
    The numerical density interface permits untouched reference dimension one
    or two; component operator certificates cover arbitrary finite references.
    """
    weights = _weights(relative_weights)
    density = clock._validate_input(density); name = law.records.candidate(name)
    if type(raw) is not bool or type(order) is not str or order not in ('interleaved', 'endpoint', 'zero'):
        raise ValueError('original Boolean source selector and declared clock order required')
    arrays = source_arrays(weights) if raw else _ideal_arrays(weights)
    output = {}; r = len(density)//66
    rho, output['abort', 'probe'] = clock._export(arrays['probe'], density)
    frames = {(): rho}
    for step, edge in enumerate(REQUESTED_SEAMS, 1):
        following = {}
        for history, state in frames.items():
            for label, pairs, root, probability in _native_branches(name, edge):
                if not pairs:
                    continue
                key = (*history, label)
                state2 = clock._native(state, pairs, root*root if raw else probability)
                if order == 'interleaved' or order == 'endpoint' and step == 2:
                    word = 'half_clock' if order == 'interleaved' else 'clock'
                    state2, abort = clock._export(arrays[word], state2)
                else:
                    abort = np.zeros((r, r), complex)
                output['abort', 'clock', step, key] = abort
                following[key] = state2
        frames = following
    for history, state in frames.items():
        state, output['abort', 'analyzer', history] = clock._export(arrays['analyzer'], state)
        for label, pairs, root, probability in _native_branches(name, PROBE):
            if pairs:
                output['read', history, label] = clock._native(state, pairs, root*root if raw else probability)
    return output


def raw_output(outputs, name, slots):
    return clock.raw_output(outputs, name, slots)


def source_budget(relative_weights):
    weights = _weights(relative_weights); key = _key(weights)
    snapshot = _snapshot(weights)
    for name, word in snapshot:
        _executed(name, word, key)
    program = clock._program(snapshot)
    error, gain = F(), F(1)
    for eps in (clock.ROTATION_GATE, clock.NATIVE_GATE, clock.WORD_GATE,
                clock.NATIVE_GATE, clock.WORD_GATE, clock.ROTATION_GATE, clock.NATIVE_GATE):
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    pulses = len(program['probe'])+max(2*len(program['half_clock']), len(program['clock']))+len(program['analyzer'])
    pulse_gate, failure_mass = clock._processor_failure_bound(pulses)
    menu = {(p['kind'], p['angle'], len(p['targets'])) for word in program.values() for p in word}
    for kind, angle, bits in menu:
        leakage, excess = _processor_bounds(common.full.compiler.native_pulse(kind, angle), 1 << bits)
        if max(leakage, excess) > pulse_gate:
            raise ArithmeticError('learned processor complement or gain exceeds its fixed gate')
    return {'complete_instrument_halfdiamond_upper': str(error+failure_mass/2),
            'successful_processor_path_halfdiamond_upper': str(error),
            'proper_code_first_failure_mass_upper': str(failure_mass),
            'proper_code_failure_outputs_materialized': False,
            'proper_code_failure_mass_included_in_complete_bound': True,
            'channel_norm_upper': str(gain+failure_mass),
            'full512_phase_operator_gate': str(clock.WORD_GATE),
            'full512_rotation_operator_gate': str(clock.ROTATION_GATE),
            'native_flagged_frame_operator_gate': str(clock.NATIVE_GATE),
            'proper_code_per_pulse_complement_and_gain_gate': str(pulse_gate),
            'actual_unique_processor_pulses_checked': len(menu),
            'maximum_source_pulses': pulses, 'coherent_load_export_calls': 2*pulses,
            'native_continuation_attempts': 2, 'native_detector_attempts': 1,
            'quantum_address_qubits': 9, 'original_source_dimension': 512,
            'source_program_bytes': len(json.dumps(program, sort_keys=True, separators=(',', ':')).encode()),
            'source_position_counter_bits': max(map(len, program.values())).bit_length(),
            'coherent_word_stage_bits': 2,
            'first_failed_load_or_export_ordinal_bits': (2*pulses).bit_length(),
            'address_export_abort_position_and_status_bits': 3,
            'processor_complement_environment_dimension_upper': 6,
            'processor_complement_index_is_accessible_record': False,
            'requested_history_raw_bits_upper': 12,
            'source_export_abort_flags_retained': True,
            'hardware_and_physical_calibration_are_external_obligations': True,
            'numerical_certificates_are_not_statistical_or_physical_error_claims': True}


def physical_allowance(relative_weights, pulse_error=F(), transfer_error=F(), native_attempt_error=F()):
    """Conditional full-isometry errors; no hardware values are inferred."""
    values = []
    for x in (pulse_error, transfer_error, native_attempt_error):
        if (type(x) is not F or type(x.numerator) is not int or type(x.denominator) is not int
                or x.denominator <= 0 or math.gcd(x.numerator, x.denominator) != 1 or x.numerator < 0):
            raise ValueError('canonical nonnegative original rational error required')
        values.append(F(x.numerator, x.denominator))
    pulses = source_budget(relative_weights)['maximum_source_pulses']
    total = pulses*values[0]+2*pulses*values[1]+3*values[2]
    if total >= 1:
        raise ValueError('sum of supplied operator error allowances must be below one')
    delta = total/(1-total)
    return {'operator_upper': str(delta), 'halfdiamond_upper': str(delta*(2+delta)/2),
            'supplied_physical_errors_not_measured_here': True}


def source_sha256():
    paths = set(clock.source_sha256())
    paths.update({'model/local_repair_family/learned_flow.py',
                  'model/local_repair_family/test_learned_flow.py'})
    return law.source_bindings(sorted(paths))


def report(relative_weights):
    weights = _weights(relative_weights)
    theorem = exact_certificate(weights); budget = source_budget(weights)
    # A physical diagnostic, explicitly separate from the inferred program.
    density = np.diag(list(map(float, law.reference(clock.parameters())))).astype(complex)
    allowance = F(budget['complete_instrument_halfdiamond_upper']); results = {}
    for name in ('C', 'E'):
        actual = execute(density, weights, name)
        endpoint = execute(density, weights, name, order='endpoint')
        ideal = execute(density, weights, name, raw=False)
        if actual.keys() != endpoint.keys() or actual.keys() != ideal.keys():
            raise ArithmeticError('learned clock placement changed the complete output alphabet')
        if any(not np.all(np.isfinite(block)) for outputs in (actual, endpoint, ideal)
               for block in outputs.values()):
            raise ArithmeticError('complete learned source output contains nonfinite entries')
        errors = [sum(np.linalg.norm(actual[k]-target[k], 'nuc') for k in actual)/2
                  for target in (ideal, endpoint)]
        mass = sum(float(np.trace(x).real) for x in actual.values())
        if (not all(math.isfinite(x) for x in (*errors, mass))
                or errors[0] > float(allowance) or errors[1] > float(2*allowance)
                or abs(mass-1) > float(2*allowance)):
            raise ArithmeticError('complete learned source output exceeds its fixed allowance')
        results[name] = {'materialized_live_and_export_blocks': len(actual),
                         'three_attempt_full_raw_words': (2 << clock.circuit.RECORD_BITS[name])**3,
                         'structural_zero_raw_words_retained_lazily': True,
                         'source_vs_ideal_halftrace_upper': str(allowance),
                         'interleaved_vs_endpoint_halftrace_upper': str(2*allowance),
                         'complete_output_checks_passed': True}
    return {'schema': 'oph.local_repair_family.learned_flow.v1', 'theorem': theorem,
            'source_budget': budget, 'requested_complete_source_execution': results,
            'requested_source_schedule': [list(x) for x in REQUESTED_SEAMS],
            'requested_schedule_is_full_scheduler': False,
            'diagnostic_input': 'apparatus reference; not supplied to phase compilation',
            'parameter': 'two half-steps, total s=1',
            'inference_accuracy_is_separate': True, 'source_sha256': source_sha256()}
