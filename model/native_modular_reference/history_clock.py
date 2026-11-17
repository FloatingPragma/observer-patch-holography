"""Execute the transported modular clock across core-crossing native repairs.

The same finitely prepared reference supplies the input clock and every
record-conditioned output clock. A two-address probe rotation precedes two
record-adaptive native attempts; a conditioned analyzer and a complete native
current read follow them. The work is never reset during this history.
All original 512 source coordinates and export aborts are retained. Reference
normalization is only used to identify a modular channel, never to normalize
an actual observed branch. Coherent source apparatus remains declared.
"""
from __future__ import annotations

import argparse
import gzip
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import hashlib
import math
from pathlib import Path

import numpy as np

from tools import strict_json
from . import (history_transport as transport, full_flow as full,
               common_caps as common, cap_local_clock as cap,
               local_history as previous, local_preparation as preparation)
from ..native_interface import construction as interface
from ..observer import native, record_recovery as records
from ..observer import adaptive_record_recovery as adaptive
from ..observer.source_discrimination import circuit
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name('history_clock_results.json.gz')
SCHEMA = 'oph.native_modular_reference.history_clock.v2'
PARAMETER = F(4)
FIRST_EDGE = (0, 5)
READ_EDGE = (1, 5)
WORD_GATE = F(1, 10**9)
ROTATION_GATE = F(1, 10**8)
NATIVE_GATE = F(1, 10**12)
NUMERICAL_DIAGNOSTIC_GATE = 2e-10
NUMERICAL_RESIDUAL_CEILING = '1/5000000000'
RESIDUAL_INTERVAL_FIELDS = (
    'source_complete_halftrace_error_interval',
    'before_after_covariance_error_interval',
    'untouched_reference_covariance_error_interval',
    'untouched_reference_source_error_interval',
    'complete_preparation_history_error_interval',
)


def probe_pair():
    basis = native.states(binary=True)
    return tuple(basis.index(tuple(int(p in occupied) for p in range(12)))
                 for occupied in ((0, 1), (0, 2)))


def next_edge(law, record):
    return (0, 7) if records.record_current(law, record) > 0 else (1, 8)


def policy(law):
    """An explicitly bounded controller with no input other than its tape."""
    law = records.candidate(law)
    result = {(): FIRST_EDGE}
    basis = native.states(binary=True)
    for branch in transport.word_instrument(law, (FIRST_EDGE,)):
        edge = next_edge(law, branch.history[0])
        result[branch.history] = edge
        for output in branch.outputs:
            for coin in (0, 1):
                _, record = records.recorded_step(law, basis[output], edge, coin)
                result[(*branch.history, record)] = None
    return result


def branches(law):
    return transport.instrument(law, policy(law))


def ratios(branch=None, *, law=None):
    """Faithful permuted native reference; identity on non-native padding.

    A current-record port permutation extends each binary branch. Its full
    faithful reference restricts to the actual corner reference up to the
    common coin scalar. Only the core location changes, so histories sharing
    that location reuse the same word. No port permutation is executed as an
    extra service, and no graph automorphism is asserted.
    """
    sigma, codes = transport.reference(), full.support()[0]
    value = [F(1)]*512
    if branch is not None:
        sigma = transport.extended_reference(law, branch, sigma)
    elif law is not None:
        raise ValueError('a branch is required for a law-conditioned extension')
    for code, weight in zip(codes, sigma, strict=True):
        value[code] = 60*weight
    return tuple(value)


def _parameter(value):
    value = transport._fraction(value)
    if abs(value) > 16:
        raise ValueError('a rational parameter in [-16,16] required')
    return F(value.numerator, value.denominator)


@lru_cache(maxsize=64)
def _clock(s, values):
    array = full.compiler.execute(common.diagonal_word(s, values), 9)
    return array.dtype.str, array.tobytes()


def clock(s, values, raw):
    s = _parameter(s)
    if type(raw) is not bool:
        raise ValueError('explicit raw-source Boolean required')
    if type(values) is not tuple or len(values) != 512:
        raise ValueError('512 exact positive original ratios required')
    values = tuple(F(v.numerator, v.denominator) for v in
                   (transport._fraction(v, positive=True) for v in values))
    if raw:
        dtype, data = _clock(s, values)
        return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)
    return np.diag([complex(float(w))**(-1j*float(s)) for w in values])


def _pair(pair, sign):
    if (type(pair) is not tuple or len(pair) != 2 or pair[0] == pair[1]
            or any(type(i) is not int or not 0 <= i < 66 for i in pair)
            or type(sign) is not int or sign not in (-1, 1)):
        raise ValueError('two distinct original native indices and an integer sign required')
    return pair, sign


def rotation_word(pair, sign):
    _pair(pair, sign)
    codes = full.support()[0]
    return cap.pair_word(codes[pair[0]], codes[pair[1]], sign)


@lru_cache(maxsize=32)
def _rotation(pair, sign):
    array = full.compiler.execute(rotation_word(pair, sign), 9)
    return array.dtype.str, array.tobytes()


def rotation(pair, sign, raw):
    _pair(pair, sign)
    if type(raw) is not bool:
        raise ValueError('explicit raw-source Boolean required')
    if raw:
        dtype, data = _rotation(pair, sign)
        return np.frombuffer(data, dtype=np.dtype(dtype)).reshape(512, 512)
    target = np.eye(512, dtype=complex)
    a, b = (full.support()[0][i] for i in pair)
    target[np.ix_((a, b), (a, b))] = np.array([[1, sign], [-sign, 1]])/math.sqrt(2)
    return target


def certify_rotation(pair, sign, array):
    """Every original binary component, including all off-code directions."""
    _pair(pair, sign)
    if type(array) is not np.ndarray:
        raise ValueError('an original base array required')
    array = full.exact._array(array, (512, 512))
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    mid, rad = (lo+hi)/2, (hi-lo)/2
    a, b = (full.support()[0][i] for i in pair)
    expected = {(a, a): mid, (b, b): mid, (a, b): sign*mid, (b, a): -sign*mid}
    square = F()
    for i, row in enumerate(array):
        for j, value in enumerate(row):
            re, im = full.exact._component(value)
            wanted = expected.get((i, j), F(int(i == j)))
            square += (re-wanted)**2+im*im
    return full.exact._sqrt_interval(square)[1]+2*rad


def certify_native(law, edge):
    """Original native Kraus frame against independent exact endpoint maps."""
    maps = interface.fixed_maps(law, edge)
    basis, positions = native.states(binary=True), {x: i for i, x in enumerate(native.states(binary=True))}
    expected = {}
    for coin in (0, 1):
        for i, state in enumerate(basis):
            output, record = records.recorded_step(law, state, edge, coin)
            expected.setdefault((coin, record[-1]), set()).add((positions[output], i))
    if (len(maps) != len(expected) or set(expected) != {(b.coin, b.label) for b in maps}
            or any(type(v) is not int for b in maps for v in (b.coin, b.label, b.current))):
        raise ArithmeticError('native frame lost a record')
    lo, hi = full.exact._sqrt_interval(F(1, 2))
    mid, rad, square, entries = (lo+hi)/2, (hi-lo)/2, F(), 0
    for branch in maps:
        array = full.exact._array(branch.matrix, (66, 66))
        marked = expected[branch.coin, branch.label]
        mapping = tuple(branch.mapping)
        if (len(mapping) != len(marked)
                or any(type(i) is not int or type(j) is not int for i, j in mapping)
                or {(j, i) for i, j in mapping} != marked
                or branch.current != records.record_current(law, (*edge, branch.coin, branch.label))):
            raise ArithmeticError('native frame mapping or current differs from its exact matrix support')
        entries += len(marked)
        for i, row in enumerate(array):
            for j, value in enumerate(row):
                re, im = full.exact._component(value)
                square += (re-(mid if (i, j) in marked else 0))**2+im*im
    return full.exact._sqrt_interval(square)[1]+full.exact._sqrt_interval(F(entries))[1]*rad


def source_budget():
    values, pairs = {ratios()}, {probe_pair()}
    for law in ('C', 'E'):
        for branch in branches(law):
            values.add(ratios(branch, law=law))
            if all(i in branch.inputs for i in probe_pair()):
                pairs.add(tuple(branch.outputs[branch.inputs.index(i)] for i in probe_pair()))
        for edge in (FIRST_EDGE, (0, 7), (1, 8), READ_EDGE):
            full.corner.report._gate(certify_native(law, edge), NATIVE_GATE, 'native record frame')
    for value in values:
        full.corner.report._gate(common.certify_diagonal(PARAMETER, value, clock(PARAMETER, value, True)),
                                WORD_GATE, 'transported modular source word')
    for pair in pairs:
        for sign in (-1, 1):
            full.corner.report._gate(certify_rotation(pair, sign, rotation(pair, sign, True)),
                                    ROTATION_GATE, 'transported probe analyzer')
    error, gain = F(), F(1)
    # Either clock position has the same six-stage majorant. For the
    # adaptive instrument use the maximum branch word error, not their sum.
    for eps in (ROTATION_GATE, WORD_GATE, NATIVE_GATE, NATIVE_GATE, ROTATION_GATE, NATIVE_GATE):
        error = (1+eps)**2*error+eps*(2+eps)/2
        gain *= (1+eps)**2
    programs = [common.diagonal_word(PARAMETER, v) for v in sorted(values)]
    max_clock = max(map(len, programs))
    rotation_programs = [rotation_word(p, sign) for p in sorted(pairs) for sign in (-1, 1)]
    max_pulses = 2*max(map(len, rotation_programs))+max_clock
    programs += rotation_programs
    return {'complete_history_halfdiamond_upper': str(error), 'channel_norm_upper': str(gain),
            'distinct_clock_words_verified': len(values), 'distinct_analyzer_pairs_verified': len(pairs),
            'maximum_source_pulses_per_successful_history': max_pulses,
            'maximum_coherent_load_export_service_calls': 2*max_pulses,
            'native_repair_attempts_after_preparation': 2, 'native_detector_attempts': 1,
            'native_export_checks': 3, 'quantum_address_qubits': 9,
            'fine_history_and_detector_record_bits_upper': 3*(5+1+3),
            'reference_core_location_bits_upper': 8,
            'abort_location_and_status_bits': 3,
            'source_program_utf8_bytes': sum(len(full.corner.canonical(p).encode()) for p in programs),
            'controller_policy_utf8_bytes': sum(len(full.corner.canonical(
                [{'history': h, 'next_seam': edge} for h, edge in policy(law).items()]).encode())
                for law in ('C', 'E')),
            'source_position_counter_bits': max(map(len, programs)).bit_length(),
            'host_interpreter_and_physical_hardware_workspace_bounded': False,
            'physical_apparatus_error_bound': None}


def raw_menu(law):
    """Every coin/register slot, including exact zeros on binary66.

    The admitted raw register is the existing source instrument's register.
    A valid full78 label may still be impossible on the binary66 shell.
    """
    law = records.candidate(law)
    attained = {(coin, records.recorded_step(law, state, FIRST_EDGE, coin)[1][-1])
                for coin in (0, 1) for state in native.states(binary=True)}
    return tuple((coin, raw, circuit.record_values(law, coin).get(raw),
                  (coin, circuit.record_values(law, coin).get(raw)) in attained)
                 for coin in (0, 1) for raw in range(1 << circuit.RECORD_BITS[law]))


def raw_read_key(law, word):
    """Semantic terminal key, or None for a provably zero raw slot."""
    law = records.candidate(law)
    if (type(word) is not tuple or len(word) != 3 or any(
            type(slot) is not tuple or len(slot) != 2 or type(slot[0]) is not int
            or slot[0] not in (0, 1) or type(slot[1]) is not int
            or not 0 <= slot[1] < 1 << circuit.RECORD_BITS[law] for slot in word)):
        raise ValueError('three complete integer coin/raw-register slots required')
    menu = {(coin, raw): (label, attained) for coin, raw, label, attained in raw_menu(law)}
    decoded = [menu[slot] for slot in word]
    if not all(attained for _, attained in decoded):
        return None
    h1 = (*FIRST_EDGE, word[0][0], decoded[0][0])
    edge = next_edge(law, h1)
    h = (h1, (*edge, word[1][0], decoded[1][0]))
    coin, label = word[2][0], decoded[2][0]
    return ('read', h, coin, label, records.record_current(law, (*READ_EDGE, coin, label)))


def raw_output(outputs, law, word):
    """Lazy zero extension; a missing attained block is an error, never zero."""
    key = raw_read_key(law, word)
    ref = len(outputs['abort', 'probe'])
    return np.zeros((66*ref, 66*ref), complex) if key is None else outputs[key].copy()


def raw_alphabet_certificate(law, outputs):
    menu = raw_menu(law)
    slots = tuple((coin, raw) for coin, raw, _, _ in menu)
    decoded = [raw_read_key(law, word) for word in product(slots, repeat=3)]
    keys = [key for key in decoded if key is not None]
    if len(set(keys)) != len(keys) or set(keys) != {key for key in outputs if key[0] == 'read'}:
        raise ArithmeticError('complete raw alphabet lost or duplicated an attained record')
    return {'coin_register_slots_per_attempt': len(slots), 'three_attempt_raw_words': len(decoded),
            'zero_extension_raw_words': decoded.count(None),
            'attained_semantic_record_slots': len(keys),
            'raw_menu': [list(row) for row in menu],
            'zero_blocks_are_defined_on_work_and_untouched_reference': True,
            'export_abort_blocks_retained_separately': True}


def _export(array, density):
    ref, codes = len(density)//66, np.array(full.support()[0])
    outside = np.array([i for i in range(512) if i not in set(codes)])
    k = np.kron(array[np.ix_(codes, codes)], np.eye(ref))
    e = np.kron(array[np.ix_(outside, codes)], np.eye(ref))
    escaped = (e@density@e.conj().T).reshape(len(outside), ref, len(outside), ref)
    return k@density@k.conj().T, np.einsum('aras->rs', escaped)


def _native_step(density, branch, raw):
    ref = len(density)//66
    if raw:
        k = np.kron(branch.matrix, np.eye(ref))
        return k@density@k.conj().T
    value = np.zeros_like(density)
    pairs = branch.mapping
    ix = [i*ref+r for i, _ in pairs for r in range(ref)]
    ox = [j*ref+r for _, j in pairs for r in range(ref)]
    value[np.ix_(ox, ox)] = density[np.ix_(ix, ix)]/2
    return value


def execute(density, law, *, mode='transported', raw=True, s=PARAMETER):
    """Full quantum outputs; each export failure keeps the untouched reference."""
    density = previous._density(density)
    if np.any(np.diag(density).real < 0):
        raise ValueError('negative original input population')
    law, s = records.candidate(law), _parameter(s)
    if type(raw) is not bool or mode not in ('before', 'transported', 'fixed', 'zero'):
        raise ValueError('a declared clock placement and Boolean source switch required')
    if raw and s != PARAMETER:
        raise ValueError('raw execution uses the retained certified parameter')
    ref = len(density)//66
    zero = lambda: np.zeros((ref, ref), complex)
    work, abort = _export(rotation(probe_pair(), 1, raw), density)
    result = {('abort', 'probe'): abort, ('abort', 'input_clock'): zero()}
    if mode == 'before':
        work, result['abort', 'input_clock'] = _export(clock(s, ratios(), raw), work)
    admitted = {b.history: b for b in branches(law)}
    for first in interface.fixed_maps(law, FIRST_EDGE):
        h1 = (*FIRST_EDGE, first.coin, first.label)
        middle = _native_step(work, first, raw)
        edge = next_edge(law, h1)
        for second in interface.fixed_maps(law, edge):
            h = (h1, (*edge, second.coin, second.label))
            branch = admitted.get(h)
            final = _native_step(middle, second, raw)
            result['abort', 'output_clock', h] = zero()
            result['abort', 'analyzer', h] = zero()
            if branch is not None and mode in ('transported', 'fixed'):
                value = ratios(branch, law=law) if mode == 'transported' else ratios()
                final, result['abort', 'output_clock', h] = _export(clock(s, value, raw), final)
            if branch is not None and all(i in branch.inputs for i in probe_pair()):
                pair = tuple(branch.outputs[branch.inputs.index(i)] for i in probe_pair())
                final, result['abort', 'analyzer', h] = _export(rotation(pair, -1, raw), final)
            for read in interface.fixed_maps(law, READ_EDGE):
                result['read', h, read.coin, read.label, read.current] = _native_step(final, read, raw)
    if any(not np.isfinite(block).all() for block in result.values()):
        raise ArithmeticError('nonfinite complete history output')
    return result


def _probabilities(outputs):
    return {key: float(np.trace(block).real) for key, block in outputs.items()}


def _variation(left, right):
    if set(left) != set(right):
        raise ValueError('complete identical record alphabets required')
    return sum(abs(left[k]-right[k]) for k in left)/2


def prepared_target(law):
    """Complete ideal target with the ORIGINAL preparation failure branches."""
    raw_reference = previous.regional_preparation.reference(common.CORE)[0]
    failure = np.diag([float(F(1, 66)-F(5, 33)*w) for w in raw_reference]).astype(complex)
    value = {('preparation', 'herald_failure'): failure,
             ('preparation', 'warmup_timeout'): np.zeros((66, 66), complex),
             ('preparation', 'completion_timeout'): np.zeros((66, 66), complex)}
    value.update(execute(float(F(5, 33))*np.diag(list(map(float, transport.reference()))).astype(complex),
                         law, raw=False))
    return value


def _entangled_probe():
    # A channel-comparison diagnostic, not an uncharged source preparation.
    i, j = probe_pair()
    vector = np.zeros(132, complex)
    vector[2*i], vector[2*j+1] = 1/math.sqrt(2), 1j/math.sqrt(2)
    return np.outer(vector, vector.conj())


def _joined_budget(source):
    prep = preparation.budget()
    error, gain = F(source['complete_history_halfdiamond_upper']), F(source['channel_norm_upper'])
    complete = gain*F(prep['complete_flagged_reference_halftrace_upper'])+F(5, 33)*error
    return prep, complete


def _residual_interval(value, scientific_bound):
    """Admit a fresh numerical residual, retaining only its fixed ceiling.

    This interval summarizes the computed observation. It is not an interval
    arithmetic certificate or a replacement for the separate source bound.
    """
    if (type(value) not in (int, float) or not math.isfinite(value)
            or not 0 <= F(value) <= min(F(NUMERICAL_RESIDUAL_CEILING), scientific_bound)):
        raise ArithmeticError('fresh clock residual exceeds its numerical ceiling or scientific bound')
    return {'lower': '0', 'upper': NUMERICAL_RESIDUAL_CEILING}


def run():
    budget = source_budget()
    prep_budget, joined_bound = _joined_budget(budget)
    initial = np.diag(list(map(float, transport.reference()))).astype(complex)
    histories = []
    for law in ('C', 'E'):
        outputs = {mode: execute(initial, law, mode=mode, raw=False)
                   for mode in ('before', 'transported', 'fixed', 'zero')}
        raw = execute(initial, law)
        alphabet = raw_alphabet_certificate(law, raw)
        error = previous.distance(raw, outputs['transported'])
        covariance = previous.distance(outputs['before'], outputs['transported'])
        source_interval = _residual_interval(error, F(budget['complete_history_halfdiamond_upper']))
        covariance_interval = _residual_interval(covariance, F(NUMERICAL_RESIDUAL_CEILING))
        wrong = _variation(_probabilities(outputs['fixed']), _probabilities(outputs['transported']))
        signal = _variation(_probabilities(outputs['zero']), _probabilities(outputs['transported']))
        if not all(math.isfinite(v) and v > 1e-5 for v in (wrong, signal)):
            raise ArithmeticError('transported source clock lacks its positive or negative control')
        # Existing preparation executes its bounded population process, with
        # all four original terminal flags, from basis input0.
        prep = preparation.population_maps(law=law)
        prepared = {('preparation', flag): np.diag(m[:, 0]).astype(complex)
                    for flag, m in prep.items() if flag != 'success'}
        prepared.update(execute(np.diag(prep['success'][:, 0]).astype(complex), law))
        prepared_error = previous.distance(prepared, prepared_target(law))
        prepared_interval = _residual_interval(prepared_error, joined_bound)
        reference_probe = _entangled_probe()
        reference_before = execute(reference_probe, law, mode='before', raw=False)
        reference_after = execute(reference_probe, law, raw=False)
        reference_raw = execute(reference_probe, law)
        reference_covariance = previous.distance(reference_before, reference_after)
        reference_error = previous.distance(reference_raw, reference_after)
        reference_covariance_interval = _residual_interval(reference_covariance, F(NUMERICAL_RESIDUAL_CEILING))
        reference_source_interval = _residual_interval(reference_error, F(budget['complete_history_halfdiamond_upper']))
        prepared_mass = sum(_probabilities(prepared).values())
        if not math.isfinite(prepared_mass) or abs(prepared_mass-1) > NUMERICAL_DIAGNOSTIC_GATE:
            raise ArithmeticError('complete preparation/reference clock history exceeds its bound')
        histories.append({'law': law,
            'complete_raw_record_alphabet': alphabet,
            'exact_support_transport': transport.phase_ratio_certificate(law, branches(law), transport.reference()),
            'reference_extension_cores': sorted({tuple(sorted(transport.port_permutation(law, b)[i]
                                               for i in common.CORE)) for b in branches(law)}),
            'source_complete_halftrace_error_interval': source_interval,
            'before_after_covariance_error_interval': covariance_interval,
            'wrong_fixed_reference_record_variation_diagnostic': wrong,
            'nonzero_clock_record_variation_diagnostic': signal,
            'untouched_reference_covariance_error_interval': reference_covariance_interval,
            'untouched_reference_source_error_interval': reference_source_interval,
            'complete_preparation_history_error_interval': prepared_interval,
            'complete_prepared_mass_diagnostic': prepared_mass,
            'records': [{'record': list(k), 'probability': v} for k, v in _probabilities(prepared).items()]})
    paths = set(transport.source_bindings()) | set(adaptive.source_bindings()) | set(interface.SOURCE_PATHS)
    paths.update(('model/native_modular_reference/history_clock.py',
                  'model/native_modular_reference/test_history_clock.py',
                  'model/native_modular_reference/cap_local_clock.py',
                  'model/native_modular_reference/local_history.py'))
    return {'schema': SCHEMA,
            'source_sha256': bound_pins({p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}),
            'parameter': str(PARAMETER), 'source_budget': budget,
            'finite_preparation': transport.preparation_certificate(PARAMETER),
            'complete_preparation_history_halftrace_upper': full.compiled._outward(joined_bound),
            'preparation_resource_budget': {
                key: prep_budget[key] for key in ('total_native_attempts_upper', 'controller_register_bits_subtotal',
                                                'auxiliary_dimension', 'initial_work_sector')},
            'modular_transport': transport.exact_certificate(),
            'adaptive_observability': adaptive.report(), 'histories': histories,
            'scope': {'clock_reference_reprepared_between_repairs': False,
                      'all_post_preparation_native_records_retained': True,
                      'preparation_fine_tape_is_coarse_grained_by_existing_herald': True,
                      'source_error_bound_is_terminal_logical_word_and_native_export': True,
                      'native_repairs_executed_as_declared_kraus_interface': True,
                      'native_repair_source_pulse_words_reexecuted': False,
                      'faithful_extension_is_a_port_permutation_of_native_reference': True,
                      'extension_is_executed_routing_or_native_graph_symmetry': False,
                      'exact_proper_code_per_pulse_abort_is_zero': True,
                      'floating_per_pulse_apparatus_error_bound_claimed': False,
                      'fixed_bare_cap_algebra_claimed': False,
                      'passive_repair_identified_with_modular_unitary': False,
                      'laboratory_time_energy_or_particle_identification_derived': False}}


def compare_receipt(supplied, actual):
    full.corner.compare_receipt(supplied, actual)
    for supplied_row, actual_row in zip(supplied['histories'], actual['histories'], strict=True):
        variation = math.fsum(abs(a['probability']-b['probability'])
                              for a, b in zip(supplied_row['records'], actual_row['records'], strict=True))/2
        if variation > NUMERICAL_DIAGNOSTIC_GATE:
            raise ValueError('complete clock probability replay exceeds its total-variation allowance')
    for value in (supplied, actual):
        if value['schema'] != SCHEMA:
            raise ValueError('canonical clock residual-interval schema required')
        for bound in (value['source_budget']['complete_history_halfdiamond_upper'],
                      value['complete_preparation_history_halftrace_upper']):
            try:
                compatible = type(bound) is str and F(bound) >= F(NUMERICAL_RESIDUAL_CEILING)
            except (ValueError, ZeroDivisionError):
                compatible = False
            if not compatible:
                raise ValueError('scientific bound must contain the retained numerical residual interval')
        for row in value['histories']:
            for key in RESIDUAL_INTERVAL_FIELDS:
                if (type(row[key]) is not dict
                        or row[key] != {'lower': '0', 'upper': NUMERICAL_RESIDUAL_CEILING}):
                    raise ValueError('canonical numerical residual interval required')
            p = [r['probability'] for r in row['records']]
            if (any(type(v) not in (int, float) or not math.isfinite(v) or v < 0 for v in p)
                    or abs(sum(p)-1) > NUMERICAL_DIAGNOSTIC_GATE
                    or abs(sum(p)-row['complete_prepared_mass_diagnostic']) > NUMERICAL_DIAGNOSTIC_GATE
                    or min(row['wrong_fixed_reference_record_variation_diagnostic'],
                           row['nonzero_clock_record_variation_diagnostic']) <= 1e-5):
                raise ValueError('complete clock probabilities or distinguishing controls changed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', type=Path)
    group.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = strict_json.loads(full.corner.canonical(run()))
    if args.check:
        raw = args.check.read_bytes()
        if args.check.suffix == '.gz':
            raw = gzip.decompress(raw)
        compare_receipt(strict_json.loads(raw.decode('utf-8')), result)
        print('Complete transported native modular history and observability limits verified.')
    elif args.write:
        if args.write.suffix == '.gz':
            with args.write.open('xb') as stream:
                with gzip.GzipFile(filename='', mode='wb', fileobj=stream, mtime=0) as archive:
                    archive.write(full.corner.canonical(result).encode('utf-8'))
        else:
            with args.write.open('x', encoding='utf-8', newline='\n') as stream:
                stream.write(full.corner.canonical(result))
    else:
        print(full.corner.canonical(result))


if __name__ == '__main__':
    main()
