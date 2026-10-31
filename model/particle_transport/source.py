"""Complete source preparation and supplied transport of the existing dimer modes.

Fixed-number occupation inputs are supplied charge resources. Completed encoders
and swaps conserve CAR number; individual proper-code processor pulses are not
asserted to be physical CAR operations. No passive native protection is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
import itertools
import math
from types import MappingProxyType

import numpy as np

from . import algebra
from ..observer.source_scalar import compiler
from ..common_source.chiral_candidate import gaussian_source as gaussian
from ..field_thermal.thermal import execute_with_exports, failure_probability


ADMISSION = 1e-7
SEEDS = {"ground": (3,), "hole": (1, 2), "particle": (7, 11)}
# Internal arithmetic convention, not a caller-selected support or spectral
# threshold. It is below the existing 2e-11 full-unitary and 3e-11 elimination
# checks. The unchanged full source-matrix/column admission still applies
# against the requested target after this arithmetic convention is applied.
_ELIMINATION_ZERO = 2.**-44


def _frozen(value):
    value = np.asarray(value)
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def _axis(axis):
    if type(axis) is not int or axis not in (0, 1, 2):
        raise ValueError("axis must be the integer 0, 1 or 2")
    return axis


def _kind(kind):
    if type(kind) is not str or kind not in SEEDS:
        raise ValueError("ground, hole or particle preparation required")
    return kind


def _pauli_image(label, address):
    target, phase = address, 1+0j
    for bit, letter in enumerate(label):
        occupied = (address >> bit) & 1
        if letter in "XY":
            target ^= 1 << bit
        if letter == "Y":
            phase *= -1j if occupied else 1j
        elif letter == "Z" and occupied:
            phase = -phase
    return target, phase


def number_rotation_word(first, second, angle, bits=4):
    """exp(i angle Y_first,second), identity on every other occupation.

    Equal-number endpoints give even Pauli strings. Their complete commuting
    expansion implements the selective Givens rotation, without a Gray path
    through intermediate occupation sectors. Number conservation applies to
    the completed rotation; parity also holds for each expanded Pauli factor.
    """
    if (type(bits) is not int or not 1 <= bits <= 4
            or type(first) is not int or type(second) is not int
            or not 0 <= first < second < 1 << bits
            or first.bit_count() != second.bit_count()
            or type(angle) not in (int, float) or not math.isfinite(angle)):
        raise ValueError("bounded equal-number endpoints and finite real angle required")
    difference = first ^ second
    choices = ["XY" if difference >> bit & 1 else "IZ" for bit in range(bits)]
    word = []
    for labels in itertools.product(*choices):
        label = "".join(labels)
        target, phase = _pauli_image(label, first)
        if target != second:
            raise ArithmeticError("selective Pauli expansion changed its support")
        coefficient = phase.imag/(1 << (bits-1))
        if coefficient:
            word.extend(compiler.pauli_rotation(label, angle*coefficient))
    return word


def _targets(kind, axis):
    _kind(kind); _axis(axis)
    if kind == "ground":
        return np.asarray(algebra.ground_column(axis), complex).reshape(16, 1)
    return np.asarray(algebra.candidate_columns(kind, axis), complex)


def _canonical_residual(value):
    """Remove cancellation-sized real/imaginary components and signed zeros."""
    result = np.array(value, dtype=complex, copy=True)
    result.real[np.abs(result.real) <= _ELIMINATION_ZERO] = 0.
    result.imag[np.abs(result.imag) <= _ELIMINATION_ZERO] = 0.
    return result


def _canonical_angle(angle):
    """Canonicalize roundoff near integer multiples of pi/2."""
    quarter = round(angle/(math.pi/2))
    exact = quarter*(math.pi/2)
    return exact if abs(angle-exact) <= _ELIMINATION_ZERO else angle


def _phase(value):
    # A negative real entry has phase +pi, never a platform-dependent choice
    # of +/-pi from an imaginary cancellation residue or a negative zero.
    real = 0. if abs(value.real) <= _ELIMINATION_ZERO else float(value.real)
    imag = 0. if abs(value.imag) <= _ELIMINATION_ZERO else float(value.imag)
    return _canonical_angle(math.atan2(imag, real)) if real or imag else 0.


def _encoder(kind, axis):
    target = _targets(kind, axis)
    seeds = SEEDS[kind]
    number = seeds[0].bit_count()
    addresses = [j for j in range(16) if j.bit_count() == number]
    if (target.shape != (16, len(seeds)) or not np.isfinite(target).all()
            or np.linalg.norm(target.conj().T@target-np.eye(len(seeds))) > 1e-11
            or np.count_nonzero(target[[j for j in range(16) if j not in addresses]]) != 0):
        raise ArithmeticError("candidate columns must be complete orthonormal fixed-number columns")
    # Complete only the already specified number sector. The other sectors
    # remain literal identities, rather than being removed from the unitary.
    columns = [target[:, j].copy() for j in range(len(seeds))]
    for address in addresses:
        value = np.eye(16, dtype=complex)[:, address]
        for _ in range(2):
            for previous in columns:
                value -= previous*np.vdot(previous, value)
        value = _canonical_residual(value)
        norm = float(np.linalg.norm(value))
        if norm > 1e-11:
            columns.append(value/norm)
    if len(columns) != len(addresses):
        raise ArithmeticError("fixed-number basis completion failed")
    unitary = np.eye(16, dtype=complex)
    ordering = list(seeds)+[j for j in addresses if j not in seeds]
    unitary[:, ordering] = np.column_stack(columns)
    if np.linalg.norm(unitary.conj().T@unitary-np.eye(16), 2) > 2e-11:
        raise ArithmeticError("full encoder is not unitary")
    residual = _canonical_residual(unitary)
    eliminations = []
    for column_position, pivot in enumerate(addresses[:-1]):
        for row in reversed(addresses[column_position+1:]):
            residual = _canonical_residual(residual)
            x, y = residual[pivot, pivot], residual[row, pivot]
            phase_x, phase_y = _phase(x), _phase(y)
            residual[pivot] *= np.exp(-1j*phase_x)
            residual[row] *= np.exp(-1j*phase_y)
            angle = _canonical_angle(math.atan2(abs(y), abs(x))) if x != 0 or y != 0 else 0.
            cosine, sine = math.cos(angle), math.sin(angle)
            old = residual[[pivot, row]].copy()
            residual[pivot] = cosine*old[0]+sine*old[1]
            residual[row] = -sine*old[0]+cosine*old[1]
            # These two entries have exact algebraic values. Re-evaluating
            # their phase from floating cancellation would choose an arbitrary
            # future phase on an empty pivot and change zero-angle requests.
            residual[pivot, pivot] = math.hypot(abs(x), abs(y))
            residual[row, pivot] = 0.
            eliminations.append((pivot, row, angle, phase_x, phase_y))
    phases = np.array([_phase(value) for value in np.diag(residual)])
    if np.linalg.norm(residual-np.diag(np.exp(1j*phases)), 2) > 3e-11:
        raise ArithmeticError("number-preserving elimination failed")
    word = compiler.compile_diagonal(phases)
    for first, second, angle, phase_x, phase_y in reversed(eliminations):
        word += number_rotation_word(first, second, -angle)
        phases = np.zeros(16)
        phases[first], phases[second] = phase_x, phase_y
        word += compiler.compile_diagonal(phases)
    # Walsh summation can leave tiny nonzero coefficients even after the
    # elimination phases have been fixed. Canonicalize the actual emitted
    # angles, retaining every pulse/service slot and counting these literal
    # requests. The angle-L1 change is at most len(word)*this threshold, below
    # one hundredth of ADMISSION on this fixed domain. This is an angle budget,
    # not a bound on full native generators; full matrix admission is separate.
    if len(word)*_ELIMINATION_ZERO >= ADMISSION/100:
        raise ArithmeticError("encoder angle canonicalization exceeds its fixed angle budget")
    for event in word:
        if abs(event["angle"]) <= _ELIMINATION_ZERO:
            event["angle"] = 0.
    initial = np.eye(16, dtype=complex)[:, seeds]
    return initial, target, unitary, word


def preparation_word(kind="ground", axis=0):
    """A fresh literal word; its input is SEEDS[kind], not a vacuum substitute."""
    return _encoder(_kind(kind), _axis(axis))[3]


@dataclass(frozen=True)
class Preparation:
    state: np.ndarray
    target: np.ndarray
    initial: np.ndarray
    failures: tuple
    resources: object
    word: tuple
    error: float


def _distance(actual, target, failures):
    value = math.sqrt(float(np.linalg.norm(actual-target))**2+failure_probability(failures))
    if not math.isfinite(value) or value > ADMISSION:
        raise ArithmeticError("complete source columns exceed numerical Stinespring admission")
    return value


def _execute(word, bits, initial):
    """Execute and validate every labeled zero or nonzero source outcome."""
    actual, failures = execute_with_exports(word, bits, initial)
    shape = np.asarray(initial).shape
    if (np.ma.isMaskedArray(actual) or np.asarray(actual).shape != shape
            or np.asarray(actual).dtype.kind not in "iufc" or not np.isfinite(actual).all()):
        raise ArithmeticError("finite full success/reference columns required")
    if type(failures) not in (tuple, list) or len(failures) != len(word):
        raise ArithmeticError("complete first-export failure event grid required")
    for event, branches in zip(word, failures, strict=True):
        dimension = 1 << len(event["targets"])
        if type(branches) not in (tuple, list) or len(branches) != 6-dimension:
            raise ArithmeticError("complete labeled export complement required")
        expected = ((1 << bits)//dimension,)+shape[1:]
        for branch in branches:
            if (np.ma.isMaskedArray(branch) or np.asarray(branch).shape != expected
                    or np.asarray(branch).dtype.kind not in "iufc" or not np.isfinite(branch).all()):
                raise ArithmeticError("finite full failure/reference columns required")
    return actual, failures


def _counts(resources):
    """Retain stable exact resource counts, after validating float diagnostics."""
    for key in ("integrated_processor_norm_upper", "maximum_numerical_valid_code_leakage"):
        if not math.isfinite(resources[key]) or resources[key] < 0:
            raise ArithmeticError("finite nonnegative source resource diagnostics required")
    result = {key: value for key, value in resources.items() if type(value) in (int, bool)}
    # A conservative integer budget avoids platform-dependent residual bytes.
    result["integrated_processor_norm_upper_integer"] = math.ceil(
        resources["integrated_processor_norm_upper"])
    return result


def prepare(kind="ground", axis=0):
    """Execute the encoder on a fixed-number input or unknown-spin isometry.

    For hole/particle, the two columns represent an arbitrary supplied dual-rail
    spin input, including entanglement with an untouched reference. Failure
    outputs are separate Kraus columns; no successful branch is normalized.
    """
    initial, target, _, word = _encoder(_kind(kind), _axis(axis))
    if kind == "ground":
        initial, target = initial[:, 0], target[:, 0]
    actual, failures = _execute(word, 4, initial)
    error = _distance(actual, target, failures)
    return Preparation(_frozen(actual), _frozen(target), _frozen(initial), failures,
                       MappingProxyType(compiler.resources(word)),
                       tuple({**e, "targets": tuple(e["targets"])} for e in word), error)


def controlled_pauli_rotation(label, angle, pointer):
    """Phase-exact exp(i angle P) on the neutral pointer's one arm."""
    if (type(label) is not str or not label or any(x not in "IXYZ" for x in label)
            or type(pointer) is not int or not len(label) <= pointer <= 8
            or type(angle) not in (int, float) or not math.isfinite(angle)):
        raise ValueError("bounded Pauli string, separate neutral pointer and finite angle required")
    first = label+"I"*(pointer+1-len(label))
    second = first[:pointer]+"Z"
    return compiler.pauli_rotation(first, angle/2)+compiler.pauli_rotation(second, -angle/2)


def controlled_fswap_word(left, right, pointer, bits=8):
    if (type(bits) is not int or not 2 <= bits <= 8
            or type(left) is not int or type(right) is not int
            or not 0 <= left < right < bits or right != left+1
            or type(pointer) is not int or not bits <= pointer <= 8):
        raise ValueError("adjacent field modes and separate bounded pointer required")
    word = []
    for pauli, angle in (("XX", -math.pi/4), ("YY", -math.pi/4),
                         ("II", math.pi/4), ("ZI", -math.pi/4),
                         ("II", math.pi/4), ("IZ", -math.pi/4)):
        label = ["I"]*bits
        label[left], label[right] = pauli
        word += controlled_pauli_rotation("".join(label), angle, pointer)
    return word


def _ordinary_swap_word(left, right, pointer, bits):
    word = []
    for pauli, angle in (("II", math.pi/4), ("XX", -math.pi/4),
                         ("YY", -math.pi/4), ("ZZ", -math.pi/4)):
        label = ["I"]*bits
        label[left], label[right] = pauli
        word += controlled_pauli_rotation("".join(label), angle, pointer)
    return word


def block_swap_word(*, controlled=False, ordinary=False):
    if type(controlled) is not bool or type(ordinary) is not bool or ordinary and not controlled:
        raise ValueError("Boolean controlled/ordinary requests; ordinary control requires a pointer")
    order = list(range(8)); route = []
    for position, wanted in enumerate((4, 5, 6, 7, 0, 1, 2, 3)):
        location = order.index(wanted)
        for left in reversed(range(position, location)):
            route.append((left, left+1))
            order[left], order[left+1] = order[left+1], order[left]
    if len(route) != 16:
        raise ArithmeticError("four-mode block route lost a crossing")
    word = []
    for left, right in route:
        if controlled:
            word += (_ordinary_swap_word(left, right, 8, 8) if ordinary else
                     controlled_fswap_word(left, right, 8, 8))
        else:
            word += compiler.remap(gaussian.local_word("fswap"), [left, right])
    return word


def _gram(state, failures):
    matrix = np.asarray(state).reshape(len(state), -1)
    gram = matrix.conj().T@matrix
    for event in failures:
        for branch in event:
            value = np.asarray(branch).reshape(len(branch), -1)
            gram += value.conj().T@value
    if np.linalg.norm(gram-np.eye(gram.shape[0]), 2) > ADMISSION:
        raise ArithmeticError("complete source preparation/transport lost input mass")
    return gram


def _failure_slots(failures):
    return sum(len(event) for event in failures)


def _hadamard(pointer):
    return compiler.remap(compiler.compile_real_orthogonal(
        np.array([[1., 1.], [1., -1.]])/math.sqrt(2)), [pointer])


def transport(kind="hole", axis=0):
    """Transport the same unknown charged doublet between two supplied dimers."""
    _kind(kind)
    if kind not in ("hole", "particle"):
        raise ValueError("a charge-one hole or particle doublet required")
    local_initial, local_target, _, local_word = _encoder(kind, _axis(axis))
    ground_initial, ground_target, _, ground_word = _encoder("ground", axis)
    initial = np.kron(ground_initial[:, 0, None], local_initial)
    prepared_target = np.kron(ground_target[:, 0, None], local_target)
    target = np.kron(local_target, ground_target[:, 0, None])
    word = local_word+compiler.remap(ground_word, [4, 5, 6, 7])+block_swap_word()
    actual, failures = _execute(word, 8, initial)
    error = _distance(actual, target, failures)
    _gram(actual, failures)
    if np.linalg.norm(algebra.block_swap()@prepared_target-target) > 2e-11:
        raise ArithmeticError("graded transport phase changed against an even ground exterior")
    generator = np.asarray(algebra.pair_generator(axis))
    energy = np.trace(prepared_target.conj().T@generator@prepared_target)/2
    if max(np.linalg.norm(generator@prepared_target-energy*prepared_target),
           np.linalg.norm(generator@target-energy*target)) > 2e-10:
        raise ArithmeticError("transport did not preserve the same candidate energy")
    return dict(state=_frozen(actual), target=_frozen(target), initial=_frozen(initial),
                prepared_target=_frozen(prepared_target), failures=failures,
                error=error, resources=compiler.resources(word),
                first_failure_kraus_slots=_failure_slots(failures),
                unknown_input_dimension=2, arbitrary_untouched_reference=True,
                initialized_field_number=SEEDS[kind][0].bit_count()+2)


def presence_readout_word(kind="hole", axis=0):
    """Nondestructive lower-doublet read on four field bits and one pointer.

    The actual inverse encoder resolves the charged eigenspace into its two
    occupation seeds. A diagonal seed-controlled pointer flip reads only their
    union; re-encoding restores the unknown internal state on either outcome.
    """
    _kind(kind); _axis(axis)
    if kind not in ("hole", "particle"):
        raise ValueError("a charged lower doublet required for presence readout")
    word = preparation_word(kind, axis)
    phases = np.zeros(32)
    phases[[16+seed for seed in SEEDS[kind]]] = math.pi
    hadamard = _hadamard(4)
    return compiler.inverse(word)+hadamard+compiler.compile_diagonal(phases)+hadamard+word


def _presence_certificate(kind, axis=0):
    word = presence_readout_word(kind, axis)
    actual, failures = _execute(word, 5, np.eye(32, dtype=complex))
    columns = _targets(kind, axis)
    projector = columns@columns.conj().T
    other = np.eye(16)-projector
    target = np.block([[other, projector], [projector, other]])
    _distance(actual, target, failures); _gram(actual, failures)
    return dict(full_Fock_and_pointer_columns=32, resolved_doublet_dimension=2,
                resources=_counts(compiler.resources(word)),
                first_failure_kraus_slots=_failure_slots(failures),
                all_input_sectors_and_pointer_arms_checked=True,
                unknown_internal_state_and_reference_preserved=True)


def transport_readout(kind="hole", axis=0):
    """One source history: prepare, move, and read the same charged doublet."""
    _kind(kind); _axis(axis)
    if kind not in ("hole", "particle"):
        raise ValueError("a charged lower doublet required for transport and readout")
    local_initial, local_target, _, local_word = _encoder(kind, axis)
    ground_initial, ground_target, _, ground_word = _encoder("ground", axis)
    initial = np.zeros((512, 2), complex)
    initial[:256] = np.kron(ground_initial[:, 0, None], local_initial)
    target = np.kron(local_target, ground_target[:, 0, None])
    word = (local_word+compiler.remap(ground_word, [4, 5, 6, 7])+block_swap_word()
            +compiler.remap(presence_readout_word(kind, axis), [4, 5, 6, 7, 8]))
    actual, failures = _execute(word, 9, initial)
    expected = np.zeros_like(initial); expected[256:] = target
    error = _distance(actual, expected, failures); _gram(actual, failures)
    outcomes = tuple(dict(pointer=pointer, state=_frozen(actual[pointer*256:(pointer+1)*256]))
                     for pointer in (0, 1))
    return dict(outcomes=outcomes, state=_frozen(actual), target=_frozen(expected),
                initial=_frozen(initial), failures=failures, error=error,
                resources=compiler.resources(word),
                first_failure_kraus_slots=_failure_slots(failures),
                unknown_input_dimension=2, arbitrary_untouched_reference=True,
                initialized_field_number=SEEDS[kind][0].bit_count()+2,
                neutral_pointer_bits_initialized=1, neutral_pointer_reads=1,
                supplied_occupation_register_bits=8)


def exchange(kind="hole", axis=0, *, ordinary=False, dephase=False):
    """Coherent block-permutation read on two identically prepared candidates.

    This is a supplied mode permutation with a neutral pointer, not autonomous
    spatial braiding. Which-path readouts and every source export abort remain.
    """
    _kind(kind); _axis(axis)
    if type(ordinary) is not bool or type(dephase) is not bool:
        raise ValueError("Boolean exchange controls required")
    local_initial, local_target, _, prep = _encoder(kind, axis)
    initial = np.kron(local_initial[:, 0], local_initial[:, 0])
    target = np.kron(local_target[:, 0], local_target[:, 0])
    prep_word = prep+compiler.remap(prep, [4, 5, 6, 7])
    prepared, prep_failures = _execute(prep_word, 8, initial)
    _distance(prepared, target, prep_failures)
    pointer_input = np.zeros(512, complex); pointer_input[:256] = prepared
    hadamard = _hadamard(8)
    split, split_failures = _execute(hadamard, 9, pointer_input)
    # The pointer is neutral: its read does not change field CAR number.
    branches = []
    if dephase:
        for path in (0, 1):
            state = np.zeros(512, complex); state[path*256:(path+1)*256] = split[path*256:(path+1)*256]
            branches.append((path, state))
    else:
        branches.append((None, split))
    tail = block_swap_word(controlled=True, ordinary=ordinary)+hadamard
    outcomes, continued_failures = [], []
    sign = 1 if ordinary or kind == "ground" else -1
    squared_error = 0.
    for path, state in branches:
        actual, failures = _execute(tail, 9, state)
        continued_failures.extend(failures)
        if path is None:
            ideal = np.concatenate(((1+sign)*target/2, (1-sign)*target/2))
        elif path == 0:
            ideal = np.concatenate((target/2, target/2))
        else:
            ideal = np.concatenate((sign*target/2, -sign*target/2))
        squared_error += float(np.linalg.norm(actual-ideal))**2
        for pointer, output in enumerate(actual.reshape(2, 256)):
            outcomes.append(dict(which_path=path, pointer=pointer, state=_frozen(output),
                                 probability=float(np.vdot(output, output).real)))
    failures = tuple(prep_failures)+tuple(split_failures)+tuple(continued_failures)
    abort = failure_probability(failures)
    error = math.sqrt(squared_error+abort)
    total = sum(row["probability"] for row in outcomes)+abort
    if not math.isfinite(error) or error > ADMISSION or abs(total-1) > ADMISSION:
        raise ArithmeticError("complete exchange interferometer failed source admission")
    # The tail is the same conditional service after either which-path record.
    # Counts are per trial, not multiplied by the number of alternatives.
    resources = compiler.resources(prep_word+hadamard+tail)
    return dict(outcomes=tuple(outcomes), failures=failures, error=error,
                complete_probability=total, resources=resources,
                first_failure_kraus_slots=_failure_slots(failures),
                neutral_pointer_reads=2 if dephase else 1,
                initialized_field_number=2*SEEDS[kind][0].bit_count(),
                ideal_pointer_probabilities=[.5, .5] if dephase else
                    ([1., 0.] if sign == 1 else [0., 1.]))


def certificate():
    """Fresh source admission followed by stable counts and exact scope claims."""
    compiler._native.cache_clear()
    compiler._source.cache_clear()
    swap_certificate = _swap_certificate()
    preparations = {}
    for kind in SEEDS:
        initial, target, unitary, word = _encoder(kind, 0)
        full, failures = _execute(word, 4, np.eye(16, dtype=complex))
        _distance(full, unitary, failures); _gram(full, failures)
        if np.linalg.norm(full@initial-target) > ADMISSION:
            raise ArithmeticError("full encoder did not prepare its charged/ground columns")
        preparations[kind] = dict(fixed_number=SEEDS[kind][0].bit_count(),
            occupation_input_addresses=list(SEEDS[kind]), full_Fock_columns=16,
            encoded_columns=len(SEEDS[kind]), resources=_counts(compiler.resources(word)),
            supplied_occupation_register_bits=4, neutral_pointer_bits_initialized=0,
            first_failure_kraus_slots=_failure_slots(failures))
    transports = {}
    for kind in ("hole", "particle"):
        result = transport(kind)
        transports[kind] = {key: result[key] for key in
            ("resources", "first_failure_kraus_slots", "unknown_input_dimension",
             "arbitrary_untouched_reference", "initialized_field_number")}
        transports[kind]["resources"] = _counts(result["resources"])
        transports[kind]["supplied_occupation_register_bits"] = 8
        transports[kind]["neutral_pointer_bits_initialized"] = 0
    presence_reads, joined_reads = {}, {}
    for kind in ("hole", "particle"):
        presence_reads[kind] = _presence_certificate(kind)
        result = transport_readout(kind)
        joined_reads[kind] = {key: result[key] for key in
            ("first_failure_kraus_slots", "unknown_input_dimension",
             "arbitrary_untouched_reference", "initialized_field_number",
             "neutral_pointer_bits_initialized", "neutral_pointer_reads",
             "supplied_occupation_register_bits")}
        joined_reads[kind].update(resources=_counts(result["resources"]),
            complete_pointer_slots=2, ideal_pointer_probabilities=[0, 1])
    reads = []
    for kind, ordinary, dephase in (("hole", False, False), ("particle", False, False),
                                    ("ground", False, False), ("hole", True, False),
                                    ("hole", False, True)):
        result = exchange(kind, ordinary=ordinary, dephase=dephase)
        reads.append(dict(kind=kind, ordinary_swap_control=ordinary,
                          which_path_control=dephase,
                          ideal_pointer_probabilities=result["ideal_pointer_probabilities"],
                          complete_pointer_slots=len(result["outcomes"]),
                          first_failure_kraus_slots=result["first_failure_kraus_slots"],
                          neutral_pointer_reads=result["neutral_pointer_reads"],
                          supplied_occupation_register_bits=8,
                          neutral_pointer_bits_initialized=1,
                          initialized_field_number=result["initialized_field_number"],
                          resources=_counts(result["resources"])))
    return dict(source_execution_axis=0,
                preparations=preparations, transports=transports, exchange_reads=reads,
                presence_reads=presence_reads, transport_and_readout=joined_reads,
                swap_source=swap_certificate,
                source_numerical_Stinespring_Frobenius_admission="1/10000000",
                full_preparation_columns_and_all_abort_Kraus_columns_checked=True,
                transport_unknown_spin_reference_columns_checked=True,
                number_preserving_completed_field_operations=True,
                microscopic_processor_pulses_identified_as_CAR_operations=False,
                fixed_number_occupation_inputs_are_supplied_charge_resources=True,
                source_sender_reset_on_abort_and_run_terminated=True,
                supplied_disconnected_identical_dimer_reference=True,
                passive_native_protection_claimed=False,
                autonomous_spatial_exchange_claimed=False,
                physical_electron_or_positron_identified=False)


report = certificate


def _swap_certificate():
    """Bind both local complement channels and the emitted complete route.

    The local checks retain every failed export column. The independent global
    check executes all256/512 successful columns, including all number sectors,
    and catches a wrong route that local gate checks alone would not detect.
    """
    local = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                      [0, 1, 0, 0], [0, 0, 0, -1]], complex)
    controlled = np.eye(8, dtype=complex); controlled[4:, 4:] = local
    for bits, target, word in (
            (2, local, gaussian.local_word("fswap")),
            (3, controlled, controlled_fswap_word(0, 1, 2, bits=2))):
        actual, failures = _execute(word, bits, np.eye(1 << bits, dtype=complex))
        _distance(actual, target, failures); _gram(actual, failures)
    full_target = np.asarray(algebra.block_swap())
    full_controlled = np.eye(512, dtype=complex); full_controlled[256:, 256:] = full_target
    counted = {}
    for bits, target, controlled in ((8, full_target, False), (9, full_controlled, True)):
        word = block_swap_word(controlled=controlled)
        counted[controlled] = compiler.resources(word)["native_pulses"]
        actual = compiler.execute(word, bits)
        if not np.isfinite(actual).all() or np.linalg.norm(actual-target) > ADMISSION:
            raise ArithmeticError("emitted full-Fock route disagrees with graded block permutation")
    return dict(adjacent_FSWAPs=16, uncontrolled_native_pulses=counted[False],
                controlled_native_pulses=counted[True],
                local_complete_input_columns=4, controlled_local_complete_input_columns=8,
                global_full_Fock_columns=256, controlled_global_columns=512,
                all_local_source_export_complements_retained=True,
                full_emitted_route_checked_against_independent_permutation=True,
                global_complete_route_Stinespring_composition_bound="16/10000000",
                global_success_Frobenius_admission="1/10000000",
                individual_local_complete_Stinespring_admission="1/10000000",
                arbitrary_untouched_reference_extension=True)
