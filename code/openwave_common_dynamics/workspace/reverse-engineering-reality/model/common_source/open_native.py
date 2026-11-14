"""Complete native instruments around the canonical mixed work history.

The executable input is a native/reference state independent of the declared
work preparation. Each output is represented by a native/reference amplitude
and an unnormalized mixed work density, never by an invented work Kraus map.
Environment alternatives have already been discarded by a complete CP sum.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import math

import numpy as np

from . import encoded_memory as memory
from . import native_routing as routing
from . import open_environment as environment
from . import open_protocol as mixed
from . import protocol, quantum_native as native


@dataclass(frozen=True)
class OpenNativeBranch:
    record: tuple
    abort: bool
    native_map: np.ndarray
    work_density: np.ndarray


@dataclass(frozen=True)
class OpenNativeHistory:
    law: str
    dimension: int
    steps: int
    backend: str
    noise_backend: str
    native_mode: str
    environment: environment.Environment
    work_histories: tuple[mixed.MixedHistory, ...]
    branches: tuple[OpenNativeBranch, ...]
    seal: str = ''


def _header(history):
    return {'law': history.law, 'dimension': history.dimension, 'steps': history.steps,
            'backend': history.backend, 'noise_backend': history.noise_backend,
            'native_mode': history.native_mode,
            'environment': [{'after': e.after, 'site': e.site, 'variance': str(e.variance),
                             'order': e.order} for e in history.environment.events],
            'work_seals': [h.seal for h in history.work_histories],
            'records': [(b.record, b.abort) for b in history.branches]}


def _seal(history):
    digest = hashlib.sha256(protocol._canonical(_header(history)))
    # Different visible records reuse the same immutable matrices. Their
    # contents are still bound at every occurrence, without repeatedly hashing
    # a large shared matrix on the host.
    cached = {}
    for branch in history.branches:
        for array in (branch.native_map, branch.work_density):
            key = id(array)
            if key not in cached:
                cached[key] = hashlib.sha256(np.ascontiguousarray(array, dtype=complex).tobytes()).digest()
            digest.update(protocol._canonical(list(array.shape)))
            digest.update(cached[key])
    return digest.hexdigest()


def _check(history, *, seal=True):
    if type(history) is not OpenNativeHistory:
        raise ValueError('an executed OpenNativeHistory is required')
    native._law(history.law)
    native._mode(history.native_mode)
    base = protocol.continuation_program(history.dimension, history.steps)
    env = environment.validated(base, history.environment)
    if (type(history.backend) is not str or history.backend not in protocol.BACKENDS
            or type(history.noise_backend) is not str or history.noise_backend not in mixed.NOISE_BACKENDS
            or type(history.work_histories) is not tuple or len(history.work_histories) != 3
            or type(history.branches) is not tuple):
        raise ValueError('complete typed mixed/native execution metadata is required')
    work = {}
    for current, h in zip((-1, 0, 1), history.work_histories, strict=True):
        mixed._check(h)
        expected = protocol.continuation_program(history.dimension, history.steps, native_current=current)
        if (h.program != expected or h.environment != env or h.cursor != len(expected.operations)
                or h.backend != history.backend or h.noise_backend != history.noise_backend
                or h.reference_dimension != 1 or h.preparation_mode != 'declared_vacuum'):
            raise ValueError('native attachment requires the complete independently prepared work program')
        work[current] = {b.record: b.density for b in h.branches}
    prepared = history.work_histories[1].initial_density
    if any(not np.array_equal(h.initial_density, prepared) for h in history.work_histories):
        raise ValueError('all currents must act on the same initial work preparation')
    menu, seen, arrays = native._menu(history.law), set(), set()
    for branch in history.branches:
        if type(branch) is not OpenNativeBranch:
            raise ValueError('typed native/mixed branches are required')
        key = branch.record
        if (type(key) is not tuple or len(key) != 5 or type(key[2]) is not tuple
                or any(type(x) is not int for x in (*key[:2], *key[2]))
                or any(x is not None and type(x) is not int for x in key[3:])
                or key not in menu or key in seen or type(branch.abort) is not bool
                or branch.abort != menu[key]):
            raise ValueError('every complete raw native outcome must have its original label')
        seen.add(key)
        if not isinstance(branch.native_map, np.ndarray) or not isinstance(branch.work_density, np.ndarray):
            raise ValueError('executed native maps and work densities must be numeric arrays')
        if id(branch.native_map) not in arrays:
            native._array(branch.native_map, (512, native.input_dimension(history.native_mode)), 'native map')
            arrays.add(id(branch.native_map))
        if id(branch.work_density) not in arrays:
            native._array(branch.work_density, (base.work_dimension,)*2, 'work density')
            arrays.add(id(branch.work_density))
        current = native._record_current(history.law, key[0], key[1])
        expected = prepared if not key[2] else work[current][key[2]]
        if not np.array_equal(branch.work_density, expected):
            raise ValueError('native outcome lost or replaced its retained mixed work branch')
    if seen != set(menu):
        raise ValueError('raw abort, zero and padding outcomes cannot be omitted')
    if seal and (not history.seal or history.seal != _seal(history)):
        raise ValueError('native/mixed history integrity digest changed')
    return history


def execute(env=None, law='C', dimension=2, steps=512, *, backend='source',
            noise_backend='source', native_mode='protected_qubit'):
    """Execute native repair, a complete noisy work history, then native repair.

    The first raw abort halts with the un-kicked prepared field retained.
    Otherwise the first decoded current selects the work kick; smooth record
    r2 selects the final seam. No work or native state is prepared again.
    """
    law, native_mode = native._law(law), native._mode(native_mode)
    base = protocol.continuation_program(dimension, steps)
    env = environment.validated(base, environment.Environment() if env is None else env)
    workers = tuple(mixed.execute(protocol.continuation_program(dimension, steps, native_current=j),
                                  env, backend=backend, noise_backend=noise_backend)
                    for j in (-1, 0, 1))
    by_current = dict(zip((-1, 0, 1), workers, strict=True))
    prepared = workers[1].initial_density
    branches = []
    for c1, raw1, current, abort, first, endings in native._native_tree(law, native_mode):
        if abort:
            branches.append(OpenNativeBranch((c1, raw1, (), None, None), True, first, prepared))
            continue
        for middle in by_current[current].branches:
            for c2, raw2, stopped, complete in endings[middle.record[2]]:
                branches.append(OpenNativeBranch((c1, raw1, middle.record, c2, raw2),
                                                stopped, complete, middle.density))
    result = OpenNativeHistory(law, dimension, steps, backend, noise_backend,
                              native_mode, env, workers, tuple(branches))
    _check(result, seal=False)
    return replace(result, seal=_seal(result))


def output(history, initial):
    """Stream (raw record, native/reference amplitude, mixed work density).

    Each tuple represents |amplitude><amplitude| tensor work_density.
    Both factors remain unnormalized. The work input is independently prepared;
    this API does not accept an initially correlated native/work state.
    """
    _check(history)
    initial = native._input(initial, native.input_dimension(history.native_mode))
    snapshots = {}

    def retained(array):
        key = id(array)
        if key not in snapshots:
            owner = array
            while isinstance(owner, np.ndarray):
                owner = owner.base
            # The normal source/runtime arrays have a bytes owner, which
            # cannot be made writable. A caller's equal-content replacement
            # may instead have mutable storage; detach that object once.
            # A fresh view also detaches mutable ndarray shape/stride metadata
            # while retaining the original immutable storage without a copy.
            snapshots[key] = array.view() if isinstance(owner, bytes) else native._immutable(array)
        return snapshots[key]

    # Frozen branch envelopes remain mutable through __dict__/object.__setattr__.
    # Bind their fields now, rather than after returning the lazy iterator.
    rows = tuple((b.record, retained(b.native_map), retained(b.work_density))
                 for b in history.branches)
    return ((record, native._immutable(matrix@initial), density)
            for record, matrix, density in rows)


def _weight(density):
    return float(np.trace(density).real)


def audit(history):
    """Check full native maps and all-input prepared-work trace preservation.

    These are floating diagnostics, not directed numerical or physical bounds.
    Independent endpoint operators check every complex native column. The work
    maps and environment channels have their separate mixed-runtime audits.
    """
    _check(history)
    size = native.input_dimension(history.native_mode)
    targets = (native._ideal_native_tree(history.law) if history.native_mode == 'full_carrier'
               else native._ideal_encoded_tree(history.law, history.native_mode))
    zero = np.zeros((512, size), complex)
    effect = np.zeros((size, size), complex)
    maximum, abort_mass, padding_mass = 0., 0., 0.
    padding = np.array([index for index in range(512) if index not in set(native.CODES)])
    residuals = {}
    for branch in history.branches:
        c1, raw1, bits, c2, raw2 = branch.record
        target = zero if not bits else targets.get((c1, raw1, bits[2], c2, raw2), zero)
        key = (id(branch.native_map), id(target))
        if key not in residuals:
            residuals[key] = float(np.linalg.norm(branch.native_map-target, 'fro'))
        maximum = max(maximum, residuals[key])
        weight = _weight(branch.work_density)
        gram = branch.native_map.conj().T@branch.native_map
        effect += weight*gram
        mass = weight*float(np.trace(gram).real)/size
        abort_mass += mass*branch.abort
        padding_mass += weight*float(np.linalg.norm(branch.native_map[padding], 'fro')**2)/size
    error = float(np.linalg.norm(effect-np.eye(size), 2))
    if any(not math.isfinite(x) or x > native.TOLERANCE for x in (maximum, error, abort_mass, padding_mass)):
        raise ArithmeticError('complete mixed/native output differs from its independent native contract')
    return {'raw_outcome_count_including_abort': len(history.branches),
            'first_abort_outcome_count': sum(not b.record[2] for b in history.branches),
            'final_abort_outcome_count': sum(b.abort and bool(b.record[2]) for b in history.branches),
            'maximum_native_map_frobenius_error_diagnostic': maximum,
            'prepared_work_native_input_completeness_error_diagnostic': error,
            'maximally_mixed_native_abort_mass_diagnostic': abort_mass,
            'maximally_mixed_native_padding_mass_diagnostic': padding_mass,
            'all_native_output_rows': 512, 'work_output_is_mixed_density': True,
            'arbitrary_initial_native_work_correlations_executed': False,
            'floating_or_physical_error_certified': False}


def _encoded(history):
    _check(history)
    if history.native_mode not in ('protected_qubit', 'passive_qubit'):
        raise ValueError('an encoded logical input mode is required')


def phase_distance(history):
    """Complete pre-read CQ trace distance of the independent logical +/- inputs."""
    _encoded(history)
    flip = np.array([[0., 1.], [1., 0.]])
    terms, plus, minus = [], [], []
    for b in history.branches:
        weight = _weight(b.work_density)
        _, triangular = np.linalg.qr(b.native_map, mode='reduced')
        delta = triangular@flip@triangular.conj().T
        terms.append(weight*float(np.abs(np.linalg.eigvalsh(delta)).sum())/2)
        plus.append(weight*float(np.linalg.norm(b.native_map[:, 0]+b.native_map[:, 1])**2)/2)
        minus.append(weight*float(np.linalg.norm(b.native_map[:, 0]-b.native_map[:, 1])**2)/2)
    distance = math.fsum(terms)
    if (not math.isfinite(distance) or not 0 <= distance <= 1+native.TOLERANCE
            or any(abs(math.fsum(values)-1) > native.TOLERANCE for values in (plus, minus))):
        raise ArithmeticError('phase distance requires complete normalized outputs')
    return distance


def logical_entanglement_fidelity(history):
    """Bell fidelity with the encoded target, after tracing work and all records."""
    _encoded(history)
    bell = np.eye(2, dtype=complex)/math.sqrt(2)
    target = native._logical_code()@bell
    terms, masses = [], []
    for _, amplitude, density in output(history, bell):
        weight = _weight(density)
        terms.append(weight*abs(np.vdot(target, amplitude))**2)
        masses.append(weight*float(np.vdot(amplitude, amplitude).real))
    fidelity, mass = math.fsum(terms), math.fsum(masses)
    if (not math.isfinite(fidelity) or not 0 <= fidelity <= 1+native.TOLERANCE
            or abs(mass-1) > native.TOLERANCE):
        raise ArithmeticError('logical/reference diagnostic omitted or renormalized an outcome')
    return fidelity


def memory_witness(active, passive):
    """Instantiate the existing memory theorem for the actual mixed work outputs."""
    _encoded(active)
    _encoded(passive)
    fields = ('law', 'dimension', 'steps', 'backend', 'noise_backend', 'environment')
    if (active.native_mode != 'protected_qubit' or passive.native_mode != 'passive_qubit'
            or any(getattr(active, field) != getattr(passive, field) for field in fields)):
        raise ValueError('matched active and passive noisy observers are required')
    probabilities = {h.program.native_current: math.fsum(b.probability for b in h.branches
                                                       if b.record[2] == 0)
                     for h in active.work_histories}
    qbar = (probabilities[-1]+2*probabilities[0]+probabilities[1])/4
    distances = {'active': phase_distance(active), 'passive': phase_distance(passive)}
    fidelity = logical_entanglement_fidelity(active)
    if (abs(distances['active']-1) > native.TOLERANCE
            or abs(distances['passive']-(1-qbar)) > native.TOLERANCE
            or abs(fidelity-1) > native.TOLERANCE):
        raise ArithmeticError('native memory differs from the recorded-channel theorem')
    return {'second_smooth_plus_probability_diagnostic': qbar,
            'pre_read_phase_distances_diagnostic': distances,
            'active_encoded_bell_fidelity_diagnostic': fidelity,
            'passive_distance_formula': '1 - qbar',
            'passive_any_decoder_error_lower_diagnostic': qbar/2,
            'existing_work_only_CPTP_memory_theorem_instantiated': True,
            'native_storage_or_record_fault_tolerance_certified': False,
            'terminal_destructive_read_executed': False}


def resources(history):
    """Maximum compatible realized path, before a destructive native read.

    The supplied unknown input is not prepared by diagnostic Bell/phase pulses.
    The encoded modes pay the actual encoder. All field/environment resources
    are inherited from the corresponding continuing mixed work branch.
    """
    _check(history)
    if history.backend != 'source' or history.noise_backend != 'source':
        raise ValueError('source field and environment execution are required for source resources')
    keys = mixed.RESOURCE_KEYS
    empty = dict.fromkeys(keys, 0)
    encoder = memory.encoding_resources() if history.native_mode != 'full_carrier' else {}
    first = routing.routed_repair_resources(history.law, (0, 1))
    common = {key: encoder.get(key, 0)+first.get(key, 0) for key in keys}
    work = {h.program.native_current: {b.record: dict(b.resources) for b in h.branches}
            for h in history.work_histories}
    prepared = dict(history.work_histories[1].preparation_resources)
    recovered = history.native_mode == 'protected_qubit'
    cache = {}
    final_costs = {bit: routing.routed_repair_resources(
        history.law, native.final_seam(history.native_mode, bit)) for bit in (0, 1)}

    def recovery(seam, coin, raw):
        if not recovered:
            return empty
        key = seam, coin, raw
        if key not in cache:
            cache[key] = memory.recovery_resources(history.law, seam, coin, raw)
        return cache[key]

    maxima, witness = dict.fromkeys(keys, 0), None
    for b in history.branches:
        c1, raw1, bits, c2, raw2 = b.record
        if not bits:
            costs = {key: common[key]+prepared[key] for key in keys}
        else:
            current = native._record_current(history.law, c1, raw1)
            seam = native.final_seam(history.native_mode, bits[2])
            last = final_costs[bits[2]]
            first_recovery = recovery((0, 1), c1, raw1)
            last_recovery = empty if b.abort else recovery(seam, c2, raw2)
            costs = {key: common[key]+work[current][bits][key]
                     +sum(stage.get(key, 0) for stage in (first_recovery, last, last_recovery))
                     for key in keys}
        if witness is None or costs['native_pulses'] > maxima['native_pulses']:
            witness = b.record
        maxima = {key: max(maxima[key], costs[key]) for key in keys}
    recorded_witness = [witness[0], witness[1], list(witness[2]), witness[3], witness[4]]
    return {'maximum_realized_path_resources': maxima, 'native_pulse_maximizing_record': recorded_witness,
            'maximum_resources_are_componentwise_over_compatible_paths': True,
            'native_input_qubits': 9,
            'supplied_unknown_logical_qubits': int(history.native_mode != 'full_carrier'),
            'native_encoder_count': int(history.native_mode != 'full_carrier'),
            'maximum_native_coin_initializations': 2, 'maximum_native_record_initializations': 2,
            'retained_classical_bits_including_two_status_bits':
                2*(1+native.circuit.RECORD_BITS[history.law])+4+2,
            'native_logical_storage_through_work_required': history.native_mode != 'full_carrier',
            'diagnostic_input_preparation_and_terminal_read_pulses': 0,
            'environment_draws_are_observer_records': False,
            'source_counts_include_every_repeated_invocation': True,
            'host_tree_storage_or_physical_service_errors_counted': False}
