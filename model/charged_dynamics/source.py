"""Complete source execution of the uniform four-site interacting action.

The known dimer ground and logical spin input are supplied preparations. Each
local word is contracted on its entire input algebra. The immutable archive
retains the actual six-level primitives and the input, so every first failed
export can be reconstructed with its separate complement Kraus outputs. No
failure or field sector is clipped, and no success state is renormalized.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import math
from types import MappingProxyType
import weakref

import numpy as np
from scipy.linalg import expm

from . import chain
from ..observer.source_scalar import compiler
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..full_fock_field.local_simulation import split_bound
from ..full_fock_field.field import FullFock
from ..interacting_field import source as local_source


SMALL_GATE = F(1, 10**10)
ENCODER_GATE = F(1, 10**7)
CONTROL_DECIMALS = 12
CONTROL_ZERO = 1e-12
COUNT_KEYS = ('native_pulses', 'single_code_pulses', 'pair_code_pulses',
              'zero_angle_requests', 'identity_phase_pulses',
              'proper_code_load_events', 'proper_code_export_events',
              'proper_code_service_events', 'integrated_processor_norm_upper')
SEEDS = {'hole': (1, 2), 'particle': (7, 11)}
_MACRO_AUTHORITY = {}
_ARCHIVE_AUTHORITY = {}


def _array(value, rows, *, columns=256):
    if np.ma.isMaskedArray(value):
        raise ValueError('masked quantum input is not admitted')
    value = np.asarray(value)
    if (value.ndim not in (1, 2) or value.shape[0] != rows
            or value.dtype.kind not in 'iufc' or not np.isfinite(value).all()
            or (value.ndim == 2 and not 1 <= value.shape[1] <= columns)):
        raise ValueError('finite bounded complete input/reference columns required')
    return np.array(value, dtype=complex, copy=True)


def _time_steps(duration, steps):
    if type(duration) not in (int, F):
        raise ValueError('exact rational duration in [0,4] required')
    if type(duration) is F and (type(duration.numerator) is not int
            or type(duration.denominator) is not int or duration.denominator <= 0):
        raise ValueError('canonical exact duration required')
    duration = F(duration)
    if not 0 <= duration <= 4 or type(steps) is not int or not 1 <= steps <= 256:
        raise ValueError('duration in [0,4] and integer split count in 1..256 required')
    return duration, steps


def _apply(state, matrix, targets, bits):
    vector = state.ndim == 1
    result = compiler._apply(state.reshape(1 << bits, -1), matrix, list(targets), bits)
    return result[:, 0] if vector else result


def _arrange(state, targets, bits):
    state = state.reshape(1 << bits, -1)
    axes = [bits-1-j for j in targets]
    axes += [j for j in range(bits) if j not in axes]
    axes += [bits]
    return state.reshape([2]*bits+[state.shape[1]]).transpose(axes).reshape(1 << len(targets), -1)


@dataclass(frozen=True)
class Event:
    kind: str
    angle: float
    targets: tuple
    primitive: np.ndarray


@dataclass(frozen=True)
class Macro:
    bits: int
    events: tuple
    success: np.ndarray
    failure_gram: np.ndarray
    resources: object
    gate: F
    numerical_error: float

    def failure(self, event_index):
        """All complement Kraus matrices of this local first-failure event."""
        _admit_macro(self)
        if type(event_index) is not int or not 0 <= event_index < len(self.events):
            raise ValueError('a retained primitive event index required')
        prefix = np.eye(1 << self.bits, dtype=complex)
        for event in self.events[:event_index]:
            d = 1 << len(event.targets)
            prefix = _apply(prefix, event.primitive[:d, :d], event.targets, self.bits)
        event = self.events[event_index]
        d = 1 << len(event.targets)
        failed = event.primitive[d:, :d]@_arrange(prefix, event.targets, self.bits)
        failed = failed.reshape(6-d, 1 << (self.bits-len(event.targets)), 1 << self.bits)
        return tuple(chain.frozen(row) for row in failed)


def _macro_fingerprint(macro):
    """Bind a compiled object's complete consumed data, including warm caches."""
    primitive_hashes = {}
    events = []
    for event in macro.events:
        if id(event.primitive) not in primitive_hashes:
            primitive_hashes[id(event.primitive)] = (event.primitive.shape, event.primitive.dtype.str,
                                                    hashlib.sha256(event.primitive.tobytes()).hexdigest())
        events.append((event.kind, event.angle, event.targets, primitive_hashes[id(event.primitive)]))
    return (macro.bits, tuple(events), macro.success.shape, macro.success.dtype.str,
            hashlib.sha256(macro.success.tobytes()).hexdigest(), macro.failure_gram.shape,
            macro.failure_gram.dtype.str, hashlib.sha256(macro.failure_gram.tobytes()).hexdigest(),
            tuple(macro.resources.items()),
            (type(macro.gate), type(macro.gate.numerator), macro.gate.numerator,
             type(macro.gate.denominator), macro.gate.denominator), macro.numerical_error)


def _admit_macro(macro):
    retained = _MACRO_AUTHORITY.get(id(macro))
    if (type(macro) is not Macro or retained is None or retained[0]() is not macro
            or _macro_fingerprint(macro) != retained[1]):
        raise ValueError('an unchanged, owned complete instrument from compile_macro is required')


def compile_macro(word, bits, target, gate=SMALL_GATE):
    """Execute a whole local instrument; admission is on every input column.

    The squared Frobenius norm of the complete Stinespring difference is
    ||success-target||_F^2 + sum_fail ||K_fail||_F^2. It bounds operator
    error for arbitrary references. The calculation is a floating diagnostic;
    its fixed admission and subsequent exact propagation are not interval or
    physical hardware certificates.
    """
    if (type(bits) is not int or not 1 <= bits <= 4 or type(gate) is not F
            or type(gate.numerator) is not int or type(gate.denominator) is not int
            or gate.denominator <= 0 or gate not in (SMALL_GATE, ENCODER_GATE)):
        raise ValueError('bounded local source macro required')
    size = 1 << bits
    target = _array(target, size)
    if target.shape != (size, size) or np.linalg.norm(target.conj().T@target-np.eye(size), 2) > 2e-11:
        raise ValueError('complete unitary macro target required')
    word = compiler._validate(word, bits)
    events, primitives = [], {}
    success = np.eye(size, dtype=complex)
    gram = np.zeros((size, size), complex)
    for request in word:
        key = request['kind'], request['angle']
        if key not in primitives:
            full = _array(compiler.native_pulse(*key), 6)
            if full.shape != (6, 6):
                raise ArithmeticError('source primitive must retain all six levels')
            defect = float(np.linalg.norm(full.conj().T@full-np.eye(6), 2))
            if not math.isfinite(defect) or defect > 2e-11:
                raise ArithmeticError('source primitive failed its full unitary admission')
            primitives[key] = chain.frozen(full)
        full = primitives[key]
        targets = tuple(request['targets']); d = 1 << len(targets)
        event = Event(*key, targets, full)
        arranged = _arrange(success, targets, bits)
        failure = (full[d:, :d]@arranged).reshape(-1, size)
        gram += failure.conj().T@failure
        success = _apply(success, full[:d, :d], targets, bits)
        events.append(event)
    error = math.sqrt(float(np.linalg.norm(success-target, 'fro'))**2+float(np.trace(gram).real))
    if not math.isfinite(error) or error > float(gate):
        raise ArithmeticError('complete local source instrument exceeds its fixed admission')
    macro = Macro(bits, tuple(events), chain.frozen(success), chain.frozen(gram),
                  MappingProxyType(compiler.resources(word)), F(gate.numerator, gate.denominator), error)
    identity = id(macro)
    _MACRO_AUTHORITY[identity] = (weakref.ref(macro, lambda _: _MACRO_AUTHORITY.pop(identity, None)),
                                  _macro_fingerprint(macro))
    return macro


@dataclass(frozen=True)
class Placement:
    macro: Macro
    targets: tuple


@dataclass(frozen=True)
class Archive:
    """Immutable, lossless factorization of every complete first-abort output.

    Local failure matrices are reconstructed from owned primitive snapshots;
    global prefixes are reconstructed from owned success maps and input. The
    complement labels remain separate Kraus values in each coarse abort slot.
    This representation makes no assertion that quantum outputs are recoverable
    from an observer's classical tape.
    """
    bits: int
    initial: np.ndarray
    schedule: tuple
    success: np.ndarray
    failure_probability: float
    resources: object

    def failure(self, macro_index, event_index):
        _admit_archive(self)
        if type(macro_index) is not int or not 0 <= macro_index < len(self.schedule):
            raise ValueError('a retained macro index required')
        state = self.initial.copy()
        for placement in self.schedule[:macro_index]:
            state = _apply(state, placement.macro.success, tuple(reversed(placement.targets)), self.bits)
        placement = self.schedule[macro_index]
        local = placement.macro.failure(event_index)
        arranged = _arrange(state, tuple(reversed(placement.targets)), self.bits)
        reference = 1 if state.ndim == 1 else state.shape[1]
        output = []
        event = placement.macro.events[event_index]
        local_remaining = [placement.targets[j] for j in reversed(range(placement.macro.bits))
                           if j not in event.targets]
        outside = [j for j in reversed(range(self.bits)) if j not in placement.targets]
        retained = local_remaining+outside
        permutation = [retained.index(j) for j in sorted(retained, reverse=True)]+[len(retained)]
        for operator in local:
            value = operator@arranged
            value = value.reshape([2]*len(retained)+[reference]).transpose(permutation).reshape(-1, reference)
            output.append(chain.frozen(value[:, 0] if state.ndim == 1 else value))
        return tuple(output)


def _archive_fingerprint(archive):
    """Bind every consumed archive field, not only its retained local macros."""
    def array(value):
        return (type(value), id(value), value.shape, value.dtype.str,
                value.flags.writeable, hashlib.sha256(value.tobytes()).hexdigest())

    return ((type(archive.bits), archive.bits), array(archive.initial),
            type(archive.schedule), tuple((type(placement), id(placement),
                id(placement.macro), type(placement.targets),
                tuple((type(target), target) for target in placement.targets))
                for placement in archive.schedule), array(archive.success),
            (type(archive.failure_probability), archive.failure_probability),
            type(archive.resources), tuple((type(key), key, type(value), value)
                                           for key, value in archive.resources.items()))


def _admit_archive(archive):
    retained = _ARCHIVE_AUTHORITY.get(id(archive))
    if type(archive) is not Archive or retained is None or retained[0]() is not archive:
        raise ValueError('an unchanged, owned complete archive from execute_schedule is required')
    try:
        fingerprint = _archive_fingerprint(archive)
    except (AttributeError, TypeError, ValueError):
        raise ValueError('an unchanged, owned complete archive from execute_schedule is required') from None
    if fingerprint != retained[1]:
        raise ValueError('an unchanged, owned complete archive from execute_schedule is required')
    seen = set()
    for placement in archive.schedule:
        if id(placement.macro) not in seen:
            _admit_macro(placement.macro)
            seen.add(id(placement.macro))


def execute_schedule(schedule, bits, initial):
    if type(bits) is not int or not 1 <= bits <= 16:
        raise ValueError('one through sixteen source work bits required')
    initial = _array(initial, 1 << bits)
    if initial.size > 2**20 or type(schedule) not in (list, tuple) or len(schedule) > 20000:
        raise ValueError('bounded complete source execution required')
    schedule = tuple(schedule)
    seen = set()
    for placement in schedule:
        if (type(placement) is not Placement or type(placement.macro) is not Macro
                or type(placement.targets) is not tuple
                or len(placement.targets) != placement.macro.bits
                or any(type(i) is not int or not 0 <= i < bits for i in placement.targets)
                or len(set(placement.targets)) != len(placement.targets)):
            raise ValueError('distinct in-register literal macro targets required')
        if id(placement.macro) not in seen:
            _admit_macro(placement.macro)
            seen.add(id(placement.macro))
    state = initial.copy(); failures = []
    totals = {key: 0 for key in COUNT_KEYS}
    error_sum = F()
    for placement in schedule:
        macro, targets = placement.macro, placement.targets
        # A macro's local bit zero is its first physical target. _apply packs
        # targets most-significant first, so reverse the local little-endian map.
        arranged = _arrange(state, tuple(reversed(targets)), bits)
        probability = float(np.vdot(arranged, macro.failure_gram@arranged).real)
        if not math.isfinite(probability) or probability < 0:
            raise ArithmeticError('complete abort mass is nonfinite or negative')
        failures.append(probability)
        state = _apply(state, macro.success, tuple(reversed(targets)), bits)
        for key in COUNT_KEYS:
            totals[key] += macro.resources[key]
        error_sum += macro.gate
    if error_sum >= 1:
        raise ValueError('source composition exceeds the declared small-error domain')
    # product(1+epsilon_j)-1 <= exp(sum epsilon_j)-1 <= s/(1-s).
    totals.update(macro_calls=len(schedule), retained_first_export_failure_slots=totals['native_pulses'],
                  initialized_work_bits=bits,
                  conditional_complete_Stinespring_error_upper=str(error_sum/(1-error_sum)),
                  numerical_admission_is_outward_certificate=False,
                  physical_service_error_included=False)
    archive = Archive(bits, chain.frozen(initial), schedule, chain.frozen(state), math.fsum(failures),
                      MappingProxyType(totals))
    identity = id(archive)
    _ARCHIVE_AUTHORITY[identity] = (weakref.ref(archive, lambda _: _ARCHIVE_AUTHORITY.pop(identity, None)),
                                   _archive_fingerprint(archive))
    return archive


def literal_word(archive):
    """Actual pulse tape, with every scalar-phase target literally remapped."""
    _admit_archive(archive)
    return tuple({'kind': event.kind, 'angle': event.angle,
                  'targets': [placement.targets[j] for j in event.targets]}
                 for placement in archive.schedule for event in placement.macro.events)


def composition_bound(archives):
    """Exact conditional bound accounting for numerical norm growth.

    If each admitted macro Stinespring differs by at most epsilon from its
    ideal isometry, composition is bounded by product(1+epsilon)-1. With
    s=sum(epsilon)<1, the rational upper bound s/(1-s) follows. Local tests
    are floating diagnostics, so this implication is not an outward audit of
    roundoff in the source matrices or of physical service errors.
    """
    if type(archives) not in (tuple, list) or not archives or any(type(a) is not Archive for a in archives):
        raise ValueError('a nonempty executed archive list required')
    for archive in archives:
        _admit_archive(archive)
    total = sum((placement.macro.gate for archive in archives for placement in archive.schedule), F())
    if total >= 1:
        raise ValueError('small-error composition requires sum epsilon below one')
    return total/(1-total)


def _encoder_data(axis):
    """Number-preserving completion of the actual main dimer columns."""
    chain._axis(axis)
    phases = np.array([1j**((mask >> 2).bit_count()) if axis != 1 else 1
                       for mask in range(16)], complex)
    fixed = {3: chain.ground_column(axis)}
    for kind in SEEDS:
        for seed, column in zip(SEEDS[kind], chain.local_columns(kind, axis).T, strict=True):
            fixed[seed] = column
    real = np.zeros((16, 16)); column_phases = np.ones(16, complex)
    for seed, value in fixed.items():
        candidate = phases.conj()*value
        # A first substantial coordinate fixes the phase chart. Choosing the
        # largest entry is discontinuous at the equal-amplitude spin columns.
        support = np.flatnonzero(np.abs(candidate) > CONTROL_ZERO)
        if not len(support):
            raise ArithmeticError('nonzero canonical dimer column required')
        pivot = int(support[0])
        phase = candidate[pivot]/abs(candidate[pivot])
        candidate = candidate/phase
        if np.linalg.norm(candidate.imag) > 2e-12:
            raise ArithmeticError('the declared dimer basis lost its real site-gauge chart')
        real[:, seed] = candidate.real; column_phases[seed] = phase
    numbers = np.array([mask.bit_count() for mask in range(16)])
    for number in range(5):
        seeds = [j for j in range(16) if numbers[j] == number]
        assigned = [j for j in seeds if j in fixed]
        remaining = iter(j for j in seeds if j not in fixed)
        for address in seeds:
            candidate = np.eye(16)[:, address]
            for j in assigned:
                candidate -= real[:, j]*np.dot(real[:, j], candidate)
            norm = float(np.linalg.norm(candidate))
            if norm > 1e-10:
                try: destination = next(remaining)
                except StopIteration: break
                real[:, destination] = candidate/norm
                assigned.append(destination)
    target = phases[:, None]*real*column_phases[None, :]
    residuals = [np.linalg.norm(target.conj().T@target-np.eye(16), 2)]
    residuals += [np.linalg.norm(target[:, seed]-value) for seed, value in fixed.items()]
    residuals += [np.linalg.norm((numbers[:, None]-numbers[None, :])*target, 2)]
    if any(not math.isfinite(float(r)) or r > 2e-11 for r in residuals):
        raise ArithmeticError('complete number-preserving dimer encoder failed')
    word = (compiler.compile_diagonal(np.angle(column_phases))
            +_canonical_real_word(real)+compiler.compile_diagonal(np.angle(phases)))
    # Quantize control parameters only. Every request remains, and all raw
    # quantum outputs are later retained and compared with the original target.
    word = [{**event, 'angle': round(event['angle'], CONTROL_DECIMALS)} for event in word]
    return chain.frozen(target), word


def _canonical_real_word(matrix):
    """Fixed-precision control synthesis prevents arbitrary empty-pivot turns.

    This is the existing Gray elimination schedule, with control-only decimal
    canonicalization before/after each elimination. Without this, atan2 of
    two roundoff remnants can emit a large rotation and change exact resource
    counts under a one-ULP input perturbation. The full original target is
    still checked by compile_macro; no quantum state is rounded or clipped.
    """
    bits, size = 4, 16
    gray = [j ^ (j >> 1) for j in range(size)]
    residual = np.round(matrix[np.ix_(gray, gray)], CONTROL_DECIMALS)
    rotations = []
    for column in range(size-1):
        for row in range(size-1, column, -1):
            x, y = residual[row-1, column], residual[row, column]
            if abs(x) < CONTROL_ZERO: x = 0.
            if abs(y) < CONTROL_ZERO: y = 0.
            theta = 0. if x == 0 and y == 0 else math.atan2(y, x)
            c, s = math.cos(theta), math.sin(theta)
            old = residual[[row-1, row], :].copy()
            residual[row-1], residual[row] = c*old[0]+s*old[1], -s*old[0]+c*old[1]
            residual = np.round(residual, CONTROL_DECIMALS)
            rotations.append((gray[row-1], gray[row], theta))
    signs = np.where(np.diag(residual) >= 0, 1., -1.)
    if np.linalg.norm(residual-np.diag(signs), 2) > 2e-10:
        raise ArithmeticError('canonical control elimination lost its complete basis')
    phases = np.zeros(size)
    for j, sign in enumerate(signs): phases[gray[j]] = 0. if sign == 1 else math.pi
    word = compiler.compile_diagonal(phases)
    for a, b, theta in reversed(rotations):
        word.extend(compiler._gray_rotation(a, b, -theta, bits))
    return word


@lru_cache(maxsize=6)
def _encoder(axis, inverse):
    target, word = _encoder_data(axis)
    if inverse: target, word = target.conj().T, compiler.inverse(word)
    return compile_macro(word, 4, target, ENCODER_GATE)


def encoder(axis=0, *, inverse=False):
    chain._axis(axis)
    if type(inverse) is not bool:
        raise ValueError('Boolean decoder selection required')
    return _encoder(axis, inverse)


@lru_cache(maxsize=512)
def _local(kind, coefficient=0j, duration=0.):
    if kind == 'contact':
        _, info = local_source.contact(duration, coefficient.real)
        target = np.diag(np.exp(-1j*duration*coefficient.real*np.array([0, 0, 0, 1.])))
        return compile_macro(info['word'], 2, target)
    if kind == 'cnot':
        target = np.array([[1,0,0,0], [0,0,0,1], [0,0,1,0], [0,1,0,0]], complex)
        return compile_macro(compiler.cnot(0, 1), 2, target)
    if kind == 'x':
        # exp(-i*pi/2) exp(i*pi*X/2)=X, including its scalar pulse.
        word = [compiler.pulse('I', -math.pi/2, [0]), compiler.pulse('X', math.pi/2, [0])]
        return compile_macro(word, 1, np.array([[0, 1], [1, 0]], complex))
    word = gaussian.local_word(kind, coefficient, duration)
    if kind == 'fswap':
        target = np.array([[1,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,-1]], complex)
    else:
        h = np.zeros((4, 4), complex); h[1, 2] = coefficient; h[2, 1] = coefficient.conjugate()
        target = expm(-1j*duration*h)
    return compile_macro(word, 2, target)


def evolution_schedule(axis=0, duration=F(4), steps=256):
    chain._axis(axis); duration, steps = _time_steps(duration, steps)
    geometry = chain.geometry(axis)
    half = float(duration)/(2*steps)
    terms = tuple(term for group in geometry.groups for term in group)
    swap = _local('fswap'); contact = _local('contact', 1+0j, half)
    half_schedule = []
    for p, q, z in terms:
        route = tuple(range(q-1, p, -1))
        half_schedule += [Placement(swap, (left, left+1)) for left in route]
        half_schedule.append(Placement(_local('hopping', z, half), (p, p+1)))
        half_schedule += [Placement(swap, (left, left+1)) for left in reversed(route)]
    # Reverse the order of terms, retaining each term's complete forward route.
    reverse = []
    for p, q, z in reversed(terms):
        route = tuple(range(q-1, p, -1))
        reverse += [Placement(swap, (left, left+1)) for left in route]
        reverse.append(Placement(_local('hopping', z, half), (p, p+1)))
        reverse += [Placement(swap, (left, left+1)) for left in reversed(route)]
    contacts = [Placement(contact, (2*site, 2*site+1)) for _ in range(2) for site in range(4)]
    return tuple((half_schedule+contacts+reverse)*steps)


def evolve(initial, axis=0, duration=F(4), steps=256):
    """All 256 field inputs and arbitrary bounded reference-column arrays."""
    return execute_schedule(evolution_schedule(axis, duration, steps), 8, initial)


def prepare(kind='hole', axis=0, spin=None):
    """Encode an unknown supplied spin/reference, with all other bits blank.

    Counted X/CNOT seed operations supply the declared initial occupation;
    they are apparatus preparation and not closed charge-preserving dynamics.
    The map is linear: unnormalized input columns are accepted and remain
    unnormalized. A probability experiment requires squared input norm one.
    """
    chain._kind(kind); chain._axis(axis)
    spin = np.array([1, 0], complex) if spin is None else _array(spin, 2)
    logical_bit = 1 if kind == 'hole' else 3
    initial = np.zeros((256,)+spin.shape[1:], complex)
    initial[0], initial[1 << logical_bit] = spin[0], spin[1]
    occupied = (0, 4, 5) if kind == 'hole' else (0, 1, 2, 4, 5)
    schedule = [Placement(_local('x'), (j,)) for j in occupied]
    schedule.append(Placement(_local('cnot'), (logical_bit, 0 if kind == 'hole' else 2)))
    schedule.extend((Placement(encoder(axis), (0, 1, 2, 3)), Placement(encoder(axis), (4, 5, 6, 7))))
    return execute_schedule(schedule, 8, initial)


def read(initial, axis=0, *, basis='occupation'):
    """All 256 eight-bit outcomes; up to sixteen untouched reference columns.

    The dimension limit is a host-allocation limit, not a mathematical limit
    on the complete instrument or its arbitrary-reference bound.
    """
    chain._axis(axis)
    if type(basis) is not str or basis not in ('occupation', 'candidate'):
        raise ValueError('occupation or candidate read basis required')
    initial = _array(initial, 256, columns=16)
    decoder = None
    if basis == 'candidate':
        decoder = execute_schedule((Placement(encoder(axis, inverse=True), (0, 1, 2, 3)),
                                    Placement(encoder(axis, inverse=True), (4, 5, 6, 7))), 8, initial)
        initial = decoder.success
    blank = np.zeros((65536,)+initial.shape[1:], complex); blank[:256] = initial
    pointer = execute_schedule(tuple(Placement(_local('cnot'), (j, 8+j)) for j in range(8)), 16, blank)
    outcomes = tuple(chain.frozen(pointer.success[256*j:256*(j+1)]) for j in range(256))
    return decoder, pointer, outcomes


def run(kind='hole', axis=0, duration=F(4), steps=256, *, spin=None, basis='occupation'):
    """A probability experiment on a normalized supplied spin/reference.

    The squared Frobenius input norm must be one within 1e-12; no input
    normalization is performed. The lower-level APIs remain linear maps.
    """
    chain._kind(kind); chain._axis(axis); duration, steps = _time_steps(duration, steps)
    if type(basis) is not str or basis not in ('occupation', 'candidate'):
        raise ValueError('occupation or candidate read basis required')
    if spin is not None:
        spin = _array(spin, 2, columns=16)
        with np.errstate(over='ignore', invalid='ignore'):
            squared_norm = float(np.vdot(spin, spin).real)
        if not math.isfinite(squared_norm) or abs(squared_norm-1) > 1e-12:
            raise ValueError('normalized supplied spin/reference required')
    preparation = prepare(kind, axis, spin)
    evolution = evolve(preparation.success, axis, duration, steps)
    decoder, pointer, outcomes = read(evolution.success, axis, basis=basis)
    archives = (preparation, evolution)+(() if decoder is None else (decoder,))+(pointer,)
    probabilities = tuple(float(np.vdot(value, value).real) for value in outcomes)
    charge_bins = tuple(math.fsum(p for mask, p in enumerate(probabilities)
                         if (mask & 15).bit_count() == a and (mask >> 4).bit_count() == b)
                        for a in range(5) for b in range(5))
    candidate = None
    if basis == 'candidate':
        a = math.fsum(probabilities[seed+16*3] for seed in SEEDS[kind])
        b = math.fsum(probabilities[3+16*seed] for seed in SEEDS[kind])
        others = math.fsum(p for mask, p in enumerate(probabilities)
            if mask not in tuple(seed+48 for seed in SEEDS[kind])+tuple(3+16*seed for seed in SEEDS[kind]))
        candidate = (a, b, others)
    return {'preparation': preparation, 'evolution': evolution, 'decoder': decoder,
            'pointer': pointer, 'outcomes': outcomes, 'probabilities': probabilities,
            'block_number_bins': charge_bins, 'candidate_location': candidate,
            'failure_probability': math.fsum(archive.failure_probability for archive in archives),
            'conditional_source_composition_bound': composition_bound(archives),
            'resources': {key: sum(archive.resources[key] for archive in archives) for key in COUNT_KEYS},
            'initialization_resources': {
                'field_work_bits': 8, 'supplied_unknown_internal_qubits': 1,
                'supplied_field_blank_bits': 7, 'fresh_pointer_blank_bits': 8,
                'maximum_live_work_bits': 16,
                'maximum_live_count_excludes_untouched_reference_and_six_level_processor': True,
                'physical_blank_preparation_cost_known': False},
            'ideal_split_operator_bound': split_bound(chain.geometry(axis).sides, duration, steps, refined=True),
            'basis': basis, 'initial_reference': 'supplied product charged dimer plus pure dimer ground',
            'coupled_ground_reference_prepared': False}


def certificate():
    """Bounded fresh source replay; mathematical bounds and diagnostics differ."""
    # The executable APIs reuse owned macros. A scientific receipt instead
    # recomputes the actual primitive inputs, including their data files.
    for cached in (_encoder, _local, compiler._source, compiler._native, local_source._contact):
        cached.cache_clear()
    source_paths = (compiler.RESPONSE, compiler.CONTROL_RECEIPT)
    source_pins = {path: hashlib.sha256((compiler.ROOT/path).read_bytes()).hexdigest()
                   for path in source_paths}
    if compiler.source_summary()['source_sha256'] != source_pins:
        raise ArithmeticError('source control data changed during fresh admission')
    rows = []
    for kind in ('hole', 'particle'):
        geometry = chain.geometry(0)
        target = chain.bare_columns(kind, 0)[:, 0]
        ideal = FullFock(8).exact_evolve(geometry, target, 4)
        for basis in ('occupation', 'candidate'):
            result = run(kind, basis=basis)
            preparation_error = float(np.linalg.norm(result['preparation'].success-target))
            evolution_error = float(np.linalg.norm(result['evolution'].success-ideal))
            total_mass = math.fsum(result['probabilities'])+result['failure_probability']
            for value, bound in ((preparation_error, 1e-7), (abs(total_mass-1), 1e-7)):
                if not math.isfinite(value) or value > bound:
                    raise ArithmeticError('full source history failed its mass/preparation admission')
            split = result['ideal_split_operator_bound']
            source = result['conditional_source_composition_bound']
            if not math.isfinite(evolution_error) or evolution_error > float(split+source):
                raise ArithmeticError('full source history exceeds its conditional evolution budget')
            expected = ideal if basis == 'occupation' else (
                np.kron(_encoder_data(0)[0].conj().T, _encoder_data(0)[0].conj().T)@ideal)
            probability_error = .5*math.fsum(abs(p-abs(z)**2)
                for p, z in zip(result['probabilities'], expected, strict=True))
            if not math.isfinite(probability_error) or probability_error > float(split+source):
                raise ArithmeticError('complete source outcome law failed its conditional bound')
            # The lower bound is a statement about ideal source words. Floating
            # execution is additionally admitted above, not certified outward.
            lower = None
            if basis == 'candidate':
                from .certification import short_time_transfer
                amplitude_lower = F(short_time_transfer()['opposite_bare_code_amplitude_lower'])
                if amplitude_lower <= split+source:
                    raise ArithmeticError('positive transfer margin was lost')
                lower = (amplitude_lower-split-source)**2
            rows.append({'kind': kind, 'axis': 0, 'read_basis': basis,
                'duration': '4', 'split_steps': 256,
                'raw_occupation_outcomes': 256, 'block_number_bins': list(result['block_number_bins']),
                'candidate_location_probabilities': (None if result['candidate_location'] is None
                                                     else list(result['candidate_location'])),
                'complete_output_probability': total_mass,
                'retained_export_failure_mass_freshly_included': True,
                'resources': result['resources'],
                'initialization_resources': result['initialization_resources'],
                'ideal_split_operator_and_half_diamond_upper': str(split),
                'conditional_numerical_source_Stinespring_upper': str(source),
                'conditional_positive_candidate_transfer_probability_lower': None if lower is None else str(lower),
                'fresh_numerical_output_checks_passed': True,
                'numerical_mass_and_preparation_gate': '1e-7',
                'numerical_propagation_is_outward_certificate': False,
                'all_number_sectors_and_export_failures_retained': True})
    if {path: hashlib.sha256((compiler.ROOT/path).read_bytes()).hexdigest()
            for path in source_paths} != source_pins:
        raise ArithmeticError('source control data changed during fresh execution')
    return {'schema': 'oph.charged_dynamics.source.v1', 'rows': rows,
            'fresh_native_control_data_sha256': source_pins,
            'carrier': 'complete eight-mode Fock space, dimension 256',
            'initial_reference': 'supplied pure ground in each decoupled dimer; one charged input',
            'initial_reference_is_coupled_ground': False,
            'preparation_accepts_unknown_internal_qubit_and_untouched_reference': True,
            'probability_experiment_input_squared_norm_gate': '1e-12',
            'maximum_executed_reference_columns_in_complete_readout': 16,
            'preparation_seed_occupation_supplied_by_counted_X_CNOT_apparatus': True,
            'read_refinements_are_separate_experiments': True,
            'failure_archive': 'owned immutable local primitives, input, literal remaps, and prefix replay',
            'local_complete_instrument_gate': str(SMALL_GATE),
            'encoder_complete_instrument_gate': str(ENCODER_GATE),
            'encoder_control_canonicalization': {
                'decimal_places': CONTROL_DECIMALS, 'empty_pivot_threshold': '1e-12',
                'all_pulse_and_service_requests_retained': True,
                'original_full_target_admission_unchanged': True,
                'quantum_outputs_clipped_or_normalized': False},
            'physical_hardware_error_or_clock_calibration_included': False,
            'native_passive_protection_or_topological_defect_claimed': False}
