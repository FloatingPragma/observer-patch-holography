"""Prepare and read the native one-port phase contrast through source words.

Local words are executed on every local input, with their complete export
failure families. Tensor contractions retain all twelve occupation bits.
Classical preparation records and untouched ports are traced only when
forming the declared marginal readout. No branch is success-normalized.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import math
from pathlib import Path

import numpy as np

from ..native_fixed_algebra import occupation_service
from ..observer import native, quantum_mixing
from ..observer.source_discrimination import circuit
from ..observer.source_scalar import compiler
from ..reference_recovery import recovery
from .lineage import CLASSIFIER_PATH, classify_phase_loss


ROOT = Path(__file__).resolve().parents[2]
BITS = 12
SELECTOR_BITS = 11
SELECTOR_ATTEMPTS = 16
TOLERANCE = 2e-10
STATES = native.states(binary=True)
MASKS = tuple(sum(value << port for port, value in enumerate(state)) for state in STATES)
SEAMS = tuple(quantum_mixing.seams())


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _compact_word(word, active):
    positions = {port: index for index, port in enumerate(active)}
    _require(len(positions) == len(active), 'distinct local source coordinates required')
    return [dict(event, targets=[positions[target] for target in event['targets']])
            for event in word]


def _blank(bits, columns=1):
    result = np.zeros((1 << bits, columns), complex)
    result[:columns] = np.eye(columns)
    return result


def _source_run(word, bits, columns=1):
    return recovery.execute_complete(word, bits, _blank(bits, columns))


def source_runs():
    """Contract emitted M6 source pulses, including every first export failure."""
    beam = _compact_word(occupation_service._beam(), (0, 5))
    point = [compiler.pulse('X', math.pi / 2, [port]) for port in (0, 1)]
    pair = [compiler.pulse('X', math.pi / 2, [port]) for port in (0, 2)] + beam
    words = {
        'point': (point, 2, 1),
        'pair': (pair, 3, 1),
        'reference_phase': ([compiler.pulse('I', 0., [0]),
                             compiler.pulse('Z', 0., [0])], 1, 2),
        'twisted_phase': ([compiler.pulse('I', math.pi / 2, [0]),
                           compiler.pulse('Z', -math.pi / 2, [0])], 1, 2),
        'analyzer': (compiler.inverse(beam), 2, 4),
        'selector': ([compiler.pulse('Y', -math.pi / 4, [0])], 1, 1),
    }
    for law in ('C', 'E'):
        active = (0, 5, 12, 13, 14) + ((15,) if law == 'E' else ())
        word = _compact_word(occupation_service.repair_word(law, 1), active)
        words[law] = (word, len(active), 4)
    return {name: _source_run(*spec) for name, spec in words.items()}


def measurement_maps():
    """Keep all six labelled outputs of each source computational readout."""
    maps = compiler.transfer_kraus(2, coherent=False)
    _require(len(maps) == 6, 'six source measurement records required')
    total = np.zeros((2, 2), complex)
    for label, full in enumerate(maps):
        _require(full.shape == (3, 6) and np.isfinite(full).all(),
                 'finite source measurement map required')
        total += full[:, :2].conj().T @ full[:, :2]
        _require(np.count_nonzero(full) == 1 and abs(full[min(label, 2), label]) == 1,
                 'source computational readout changed')
    _require(np.linalg.norm(total - np.eye(2)) < TOLERANCE,
             'complete source measurement required')
    return maps


def selector(run, maps, attempts=SELECTOR_ATTEMPTS):
    """A bounded fair-bit selector retains source failures and exhausted retries."""
    probabilities = [float(np.linalg.norm(matrix[:, :2] @ run.success) ** 2)
                     for matrix in maps]
    bit_success = math.fsum(probabilities[:2])
    bit_failure = run.failure_probability + math.fsum(probabilities[2:])
    slots = 30 * 56
    draws = [math.prod(probabilities[(index >> bit) & 1] for bit in range(SELECTOR_BITS))
             for index in range(1 << SELECTOR_BITS)]
    retry = math.fsum(draws[slots:])
    series = math.fsum(retry ** step for step in range(attempts))
    accepted = [value * series for value in draws[:slots]]
    failed_draw = bit_failure * math.fsum(bit_success ** step for step in range(SELECTOR_BITS))
    failure = failed_draw * series
    timeout = retry ** attempts
    balance = math.fsum(accepted) + failure + timeout
    _require(abs(balance - 1) < TOLERANCE, 'selector complete probability budget failed')
    return accepted, {
        'bits_per_draw': SELECTOR_BITS, 'maximum_draws': attempts,
        'accepted_slots': slots, 'retry_slots': (1 << SELECTOR_BITS) - slots,
        'accepted_probability': math.fsum(accepted),
        'source_failure_probability': failure, 'timeout_probability': timeout,
        'normalization_error': abs(balance - 1),
        'bit_record_probabilities': probabilities,
        'slot_probabilities': accepted,
        'archive': 'literal bit records, retry indices, first source failure or timeout',
    }


def preparation_components(weights):
    """Describe the rank-56 symmetric basis for each actual source seam."""
    components = defaultdict(lambda: {'weight': 0., 'selector_slots': []})
    index = 0
    for a, b in SEAMS:
        choices = [('point', tuple(port for port, value in enumerate(state) if value))
                   for state in STATES if state[a] == state[b]]
        choices += [('pair', (a, b, spectator)) for spectator in range(BITS)
                    if spectator not in (a, b)]
        _require(len(choices) == 56, 'complete symmetric preparation basis required')
        for key in choices:
            components[key]['weight'] += weights[index]
            components[key]['selector_slots'].append(index)
            index += 1
    _require(index == len(weights), 'complete preparation selector required')
    return [{'kind': kind, 'ports': list(ports), **data}
            for (kind, ports), data in sorted(components.items())]


def prepare(components, runs):
    """Embed full source-executed local columns without a number-sector cut."""
    factor = np.zeros((1 << BITS, len(components)), complex)
    failure = 0.
    for column, component in enumerate(components):
        run = runs[component['kind']]
        ports = component['ports']
        for local, amplitude in enumerate(run.success[:, 0]):
            address = sum(((local >> bit) & 1) << port for bit, port in enumerate(ports))
            factor[address, column] = math.sqrt(component['weight']) * amplitude
        failure += component['weight'] * run.failure_probability
    return factor, failure


def _apply_local(factor, matrix, ports, bits=BITS):
    # A contracted word uses little-endian local bits. compiler._apply packs
    # its first listed target as the most significant local coordinate.
    return compiler._apply(factor, matrix, tuple(reversed(ports)), bits)


def regional_factor(factor, ports):
    """Trace untouched ports and preparation labels, preserving every QR row."""
    axes = [BITS - 1 - port for port in reversed(ports)]
    axes += [axis for axis in range(BITS) if axis not in axes] + [BITS]
    columns = factor.reshape([2] * BITS + [factor.shape[1]]).transpose(axes)
    columns = columns.reshape(1 << len(ports), -1)
    _, triangular = np.linalg.qr(columns.conj().T, mode='reduced')
    result = triangular.conj().T
    error = float(np.linalg.norm(result @ result.conj().T - columns @ columns.conj().T))
    _require(error < TOLERANCE, 'regional marginal contraction failed')
    return result, error


def _failure_probability(gram, factor):
    value = float(np.trace(factor.conj().T @ gram @ factor).real)
    _require(math.isfinite(value) and value >= -1e-25, 'invalid export failure probability')
    return value


def _repair_matrix(run, law, coin, raw):
    plus, minus = raw & 1, (raw >> 1) & 1
    parity = (raw >> 2) & 1 if law == 'E' else 0
    high = coin + 2 * plus + 4 * minus + 8 * parity
    return run.success[4 * high:4 * high + 4]


def readout(factor, ports, analyzer, maps):
    """Execute inverse beam and retain all 36 two-bit measurement records."""
    analyzed = _apply_local(factor, analyzer.success, ports, factor.shape[0].bit_length() - 1)
    local = _pair_restriction(analyzed, ports)
    failed = _failure_probability(analyzer.failure_gram, _pair_restriction(factor, ports))
    values = []
    for first, second in product(range(6), repeat=2):
        matrix = np.kron(maps[second][:, :2], maps[first][:, :2])
        probability = float(np.linalg.norm(matrix @ local) ** 2)
        values.append(probability)
    return values, failed


def _pair_restriction(factor, ports):
    bits = factor.shape[0].bit_length() - 1
    axes = [bits - 1 - port for port in reversed(ports)]
    axes += [axis for axis in range(bits) if axis not in axes] + [bits]
    columns = factor.reshape([2] * bits + [factor.shape[1]]).transpose(axes).reshape(4, -1)
    return columns


def _run_summary(run):
    word = [compiler.pulse(kind, angle, targets) for kind, angle, targets in run.word]
    return {
        'bits': run.bits, 'input_columns': run.initial.shape[1],
        'resources': compiler.resources(word), 'pulse_count': len(word),
        'all_input_export_failure_trace': run.failure_probability,
        'complete_gram_error': float(np.linalg.norm(run.success.conj().T @ run.success
                                                   + run.failure_gram
                                                   - run.initial.conj().T @ run.initial)),
        'first_failure_slots': sum(6 - (1 << len(targets)) for _, _, targets in run.word),
        'word': word,
    }


def _dependencies():
    names = [
        'model/native_phase_experiment/experiment.py',
        'model/native_phase_experiment/independent.py',
        'model/native_phase_experiment/verify.py',
        'model/native_phase_experiment/lineage.py',
        'model/native_phase_experiment/__main__.py',
        CLASSIFIER_PATH,
        'model/native_fixed_algebra/occupation_service.py',
        'model/reference_recovery/recovery.py',
        'model/reference_observation/observation.py',
        'model/observer/source_scalar/compiler.py',
        'model/observer/source_discrimination/circuit.py',
        'model/observer/quantum_mixing.py', 'model/observer/native.py',
        'model/common_source/slow_geometry.py',
        'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json',
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
    ]
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names}


def report(*, keep_arrays=False):
    runs = source_runs()
    maps = measurement_maps()
    weights, selection = selector(runs['selector'], maps)
    components = preparation_components(weights)
    prepared, preparation_failure = prepare(components, runs)
    # One-qubit phase service is contracted on the full twelve-bit factor.
    phase_local, _ = regional_factor(prepared, (0,))
    phase_gram = phase_local @ phase_local.conj().T
    phase_inputs = {}
    phase_failures = {}
    for setting in ('reference', 'twisted'):
        run = runs[setting + '_phase']
        phase_inputs[setting] = _apply_local(prepared, run.success, (0,))
        phase_failures[setting] = float(np.trace(run.failure_gram @ phase_gram).real)
    result = {
        'schema': 'oph.native_phase_experiment.v1',
        'target': {'occupation_bits': BITS, 'binary_total': 2, 'phase_port': 0,
                   'analyzer_ports': [0, 1], 'seams': [list(edge) for edge in SEAMS]},
        'pointer_records': [
            {'first': first, 'second': second, 'abort': first >= 2 or second >= 2,
             'score': int((first, second) == (1, 0)) - int((first, second) == (0, 1))}
            for first, second in product(range(6), repeat=2)],
        'selector': selection,
        'preparation': {'components': components, 'source_failure_probability': preparation_failure,
                        'success_probability': float(np.linalg.norm(prepared) ** 2),
                        'densities': {}},
        'source_execution': {name: _run_summary(run) for name, run in runs.items()},
        'source_sha256': _dependencies(),
        'protected_obstruction_reuse': classify_phase_loss(ROOT), 'trials': [],
        'resource_envelope': {
            'maximum_selector_source_pulses': SELECTOR_BITS * SELECTOR_ATTEMPTS,
            'maximum_selector_measurements': SELECTOR_BITS * SELECTOR_ATTEMPTS,
            'terminal_computational_measurements': 2,
            'maximum_trial_source_pulses': {
                law: SELECTOR_BITS * SELECTOR_ATTEMPTS + len(runs['pair'].word)
                + len(runs['twisted_phase'].word) + len(runs[law].word)
                + len(runs['analyzer'].word) for law in ('C', 'E')},
            'proper_code_load_export_events_per_source_pulse': 2,
        },
    }
    outside = sorted(set(range(1 << BITS)) - set(MASKS))
    for setting, factor in phase_inputs.items():
        body = factor[list(MASKS)]
        density = body @ body.conj().T
        result['preparation']['densities'][setting] = {
            'real': density.real.tolist(), 'imag': density.imag.tolist(),
            'outside_total_two_probability': float(np.linalg.norm(factor[outside]) ** 2),
        }
    for setting, original in phase_inputs.items():
        for law in ('before', 'C', 'E'):
            records = []
            repair_failure = analyzer_failure = 0.
            contraction_error = 0.
            seam_menu = SEAMS if law != 'before' else ((0, 1),)
            for seam in seam_menu:
                region = tuple(sorted(set((0, 1, *seam))))
                factor, error = regional_factor(original, region)
                contraction_error = max(contraction_error, error)
                read_ports = tuple(region.index(port) for port in (0, 1))
                repair_ports = tuple(region.index(port) for port in seam)
                if law == 'before':
                    branches = [(None, None, False, factor)]
                else:
                    local = _pair_restriction(factor, repair_ports)
                    repair_failure += _failure_probability(runs[law].failure_gram, local) / 30
                    branches = []
                    for coin, raw in product(range(2), range(1 << circuit.RECORD_BITS[law])):
                        output = _apply_local(factor, _repair_matrix(runs[law], law, coin, raw),
                                              repair_ports, len(region)) / math.sqrt(30)
                        branches.append((coin, raw, raw not in circuit.record_values(law, coin), output))
                for coin, raw, abort, output in branches:
                    outcomes, failed = readout(output, read_ports, runs['analyzer'], maps)
                    analyzer_failure += failed
                    records.append({'seam': list(seam), 'coin': coin, 'raw': raw,
                                    'raw_abort': abort, 'outcomes': outcomes})
            measured = math.fsum(value for record in records for value in record['outcomes'])
            score = math.fsum(value * pointer['score'] for record in records
                              for value, pointer in zip(record['outcomes'], result['pointer_records'], strict=True))
            failures = (selection['source_failure_probability'] + selection['timeout_probability']
                        + preparation_failure + phase_failures[setting]
                        + repair_failure + analyzer_failure)
            balance_error = abs(measured + failures - 1.)
            _require(balance_error < TOLERANCE, 'complete trial probability accounting failed')
            result['trials'].append({
                'setting': setting, 'law': law, 'records': records,
                'score': score, 'measured_probability': measured,
                'phase_source_failure_probability': phase_failures[setting],
                'repair_source_failure_probability': repair_failure,
                'analyzer_source_failure_probability': analyzer_failure,
                'complete_probability_error': balance_error,
                'regional_factorization_error': contraction_error,
            })
    if keep_arrays:
        arrays = {'prepared': prepared, **phase_inputs}
        for name, run in runs.items():
            for field in ('initial', 'success', 'failure_gram', 'failure_weights'):
                arrays[name + '_' + field] = getattr(run, field)
        return result, arrays
    return result
