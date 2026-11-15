"""Small complete instruments compiled through the existing native source.

This module is an adapter for an experiment, not another history executor.
All words use ``observer.source_scalar.compiler`` and act on the caller's
retained, possibly unnormalized work/reference columns. A read appends
blank helpers, executes its word and returns BOTH unnormalized pointer
outputs. Scratch-failure outputs, spectators and reference columns remain
present rather than being selected away.

The compiled proper-code controls, pointer blanks, analogue precision and
quantum readout remain declared services. Pulse counts are not a physical
clock calibration or a certificate of hardware accuracy.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from model.observer.source_scalar import compiler as c


MAX_PREPARATION_QUBITS = 8
MAX_MEASURED_QUBITS = 7  # The diagonal word also contains its pointer bit.
TOLERANCE = 2e-11


@dataclass(frozen=True)
class Readout:
    """Complete source outputs, with logical outcomes in order zero, one."""

    outputs: tuple[np.ndarray, np.ndarray]
    word: tuple[dict, ...]
    resources: dict
    failure_outputs: tuple[np.ndarray, ...] = ()
    failure_records: tuple[tuple[int, int, int], ...] = ()
    scratch_return_error: float = 0.


def _bits(size, maximum):
    if type(size) is not int or size < 2 or size & (size - 1) or size > 1 << maximum:
        raise ValueError(f"power-of-two dimension from 2 through {1 << maximum} required")
    return size.bit_length() - 1


def _original_values(value):
    if np.ma.isMaskedArray(value) or isinstance(value, (bool, np.bool_)):
        raise ValueError("masked or Boolean quantum data are inadmissible")
    if isinstance(value, (tuple, list)):
        for item in value:
            _original_values(item)


def _state(value):
    _original_values(value)
    raw = np.asarray(value)
    if raw.ndim not in (1, 2) or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all():
        raise ValueError("finite numeric work/reference columns required")
    bits = _bits(int(raw.shape[0]), 16)
    out = np.asarray(raw, dtype=complex).reshape(raw.shape[0], -1)
    if out.shape[1] < 1:
        raise ValueError("at least one reference column required")
    with np.errstate(over="ignore", invalid="ignore"):
        mass = float(np.vdot(out, out).real)
    if not math.isfinite(mass):
        raise ValueError("the retained probability mass must be finite")
    # Branches are intentionally NOT required to have unit norm.
    return out, bits


def _vector(value):
    _original_values(value)
    raw = np.asarray(value)
    if raw.ndim != 1 or raw.dtype.kind not in "iufc" or not np.isfinite(raw).all():
        raise ValueError("a finite numeric preparation vector required")
    bits = _bits(int(raw.size), MAX_PREPARATION_QUBITS)
    vector = np.asarray(raw, dtype=complex)
    with np.errstate(over="ignore", invalid="ignore"):
        mass = float(np.vdot(vector, vector).real)
    if not math.isfinite(mass) or abs(mass - 1) > 2e-12:
        raise ValueError("the supplied preparation vector must already be normalized")
    return vector, bits


def _targets(bits, work_bits, targets):
    try:
        targets = tuple(range(bits)) if targets is None else tuple(targets)
    except TypeError as exc:
        raise ValueError("measured targets must be a finite bit sequence") from exc
    if (len(targets) != bits
            or any(type(bit) is not int or not 0 <= bit < work_bits for bit in targets)
            or len(set(targets)) != bits):
        raise ValueError("distinct measured bits within the work register required")
    return targets


def _hadamard(bit):
    # The native H request is exp(i angle H), so remove its scalar i.
    return [c.pulse("I", -math.pi / 2, [bit]), c.pulse("H", math.pi / 2, [bit])]


def preparation_word(vector):
    """Prepare a known complex vector from zero on at most eight work bits.

    Absolute amplitudes use the existing real preparation compiler; the
    existing diagonal compiler supplies their phases. This is an emitted
    unitary, not a replacement of the caller's unknown state.
    """
    vector, _ = _vector(vector)
    return c.compile_real_state(np.abs(vector)) + c.compile_diagonal(np.angle(vector))


def execute_word(tape, qubits, state_matrix):
    """Apply native pulses to retained columns and return output plus costs.

    A vector input is returned as a one-column matrix. Reference columns
    are untouched. No norm repair or source-target polar correction occurs.
    """
    state, bits = _state(state_matrix)
    if type(qubits) is not int or qubits != bits:
        raise ValueError("qubit count differs from the retained work register")
    output = c.execute(tape, qubits, state=state)
    resources = {**c.resources(tape), "work_qubits": qubits,
                 "reference_columns": state.shape[1],
                 "physical_service_error_certified": False}
    return output, resources


def diagonal_read(state_matrix, effect0, *, targets=None, basis_word=()):
    """Execute the positive binary instrument for a diagonal effect.

    ``effect0`` contains probabilities in local little-endian bit order.
    ``targets`` defaults to those lowest work bits; other work bits are
    spectators. ``basis_word`` is a local word B executed before the read
    and undone afterward. The ideal maps are B* sqrt(E) B and
    B* sqrt(I-E) B, where * denotes adjoint. Both source outputs retain
    their full work/reference shape and their original probability weight.
    """
    state, work_bits = _state(state_matrix)
    if work_bits + 1 > 16:
        raise ValueError("pointer read exceeds the sixteen-qubit source workspace")
    _original_values(effect0)
    raw = np.asarray(effect0)
    if (raw.ndim != 1 or raw.dtype.kind not in "iuf" or not np.isfinite(raw).all()
            or np.any(raw < 0) or np.any(raw > 1)):
        raise ValueError("a real diagonal effect with every entry in [0,1] required")
    measured_bits = _bits(int(raw.size), MAX_MEASURED_QUBITS)
    targets = _targets(measured_bits, work_bits, targets)
    # Revalidation through remap prevents basis requests from addressing the
    # pointer or an undeclared spectator.
    basis = c.remap(basis_word, list(range(measured_bits)))
    pointer = measured_bits
    phi = np.arccos(np.sqrt(raw.astype(float)))
    phases = np.concatenate((phi, -phi))
    # H exp(i phi Z) H gives cos(phi)|0> + i sin(phi)|1>.
    # The final diagonal(1,-i) removes the known outcome-one phase so the
    # two maps are the positive square roots, including endpoints 0 and 1.
    local = (basis + _hadamard(pointer) + c.compile_diagonal(phases)
             + _hadamard(pointer)
             + [c.pulse("I", -math.pi / 4, [pointer]),
                c.pulse("Z", math.pi / 4, [pointer])]
             + c.inverse(basis))
    word = c.remap(local, [*targets, work_bits])
    prepared = np.zeros((2 * len(state), state.shape[1]), complex)
    prepared[:len(state)] = state
    output = c.execute(word, work_bits + 1, state=prepared)
    branches = (output[:len(state)].copy(), output[len(state):].copy())
    # The input may be a branch of weight less than one, or identically zero.
    incoming = float(np.vdot(state, state).real)
    outgoing = sum(float(np.vdot(branch, branch).real) for branch in branches)
    if abs(incoming - outgoing) > TOLERANCE * max(1., incoming):
        raise ArithmeticError("complete pointer outputs changed the incoming probability mass")
    resources = {**c.resources(word), "fresh_pointer_blanks": 1,
                 "committed_classical_bits": 1, "work_preparations": 0,
                 "work_qubits": work_bits, "reference_columns": state.shape[1],
                 "physical_service_error_certified": False}
    return Readout(branches, tuple(word), resources)


def projector_read(state_matrix, vector, *, targets=None):
    """Read P=|vector><vector| and I-P, retaining clock/reference spectators.

    The supplied known-state preparation V is inverted before the zero
    flag and restored afterward. Outcome zero is P, outcome one is I-P.
    There is no preparation of the work register during this read.

    Eight measured bits use two four-bit zero flags and a Toffoli pointer;
    the scratch flags are uncomputed. Their six nonzero-scratch outputs
    are retained separately as terminal failures, with records ordered
    (pointer, scratch_A, scratch_B). They vanish for the ideal circuit.
    """
    vector, bits = _vector(vector)
    if bits == 8:
        return _eight_bit_projector(state_matrix, vector, targets)
    effect = np.zeros(len(vector))
    effect[0] = 1
    return diagonal_read(state_matrix, effect, targets=targets,
                         basis_word=c.inverse(preparation_word(vector)))


def _zero_flag_word(bits):
    """Flip a blank flag iff every selected data bit is zero."""
    phases = np.zeros(1 << (bits + 1))
    phases[1 << bits] = math.pi
    return _hadamard(bits) + c.compile_diagonal(phases) + _hadamard(bits)


def _eight_bit_projector(state_matrix, vector, targets):
    state, work_bits = _state(state_matrix)
    if work_bits + 3 > 16:
        raise ValueError("eight-bit projector exceeds the sixteen-qubit source workspace")
    targets = _targets(8, work_bits, targets)
    pointer, scratch_a, scratch_b = work_bits, work_bits + 1, work_bits + 2
    prepare = c.remap(preparation_word(vector), list(targets))
    flag = _zero_flag_word(4)
    flag_a = c.remap(flag, [*targets[:4], scratch_a])
    flag_b = c.remap(flag, [*targets[4:], scratch_b])
    # CCZ on the two flags and pointer, conjugated by H on the pointer.
    phases = np.zeros(8)
    phases[7] = math.pi
    toffoli = (_hadamard(pointer)
               + c.remap(c.compile_diagonal(phases), [scratch_a, scratch_b, pointer])
               + _hadamard(pointer))
    # The raw zero flag is one on the selected vector. Flip its label so
    # callers retain the same outcome-zero=P convention as the small read.
    label_flip = [c.pulse("I", -math.pi / 2, [pointer]),
                  c.pulse("X", math.pi / 2, [pointer])]
    word = (c.inverse(prepare) + flag_a + flag_b + toffoli
            + c.inverse(flag_b) + c.inverse(flag_a) + label_flip + prepare)
    prepared = np.zeros((8 * len(state), state.shape[1]), complex)
    prepared[:len(state)] = state
    output = c.execute(word, work_bits + 3, state=prepared)
    blocks = tuple(output[index * len(state):(index + 1) * len(state)].copy()
                   for index in range(8))
    incoming = float(np.vdot(state, state).real)
    outgoing = sum(float(np.vdot(branch, branch).real) for branch in blocks)
    if abs(incoming - outgoing) > TOLERANCE * max(1., incoming):
        raise ArithmeticError("scratch-inclusive instrument changed incoming probability mass")
    failures = blocks[2:]
    records = tuple((index & 1, (index >> 1) & 1, (index >> 2) & 1)
                    for index in range(2, 8))
    error = math.sqrt(sum(float(np.vdot(branch, branch).real) for branch in failures))
    resources = {**c.resources(word), "fresh_pointer_blanks": 3,
                 "scratch_blanks": 2, "committed_classical_bits": 3,
                 "work_preparations": 0, "work_qubits": work_bits,
                 "reference_columns": state.shape[1],
                 "scratch_failure_outputs_retained": True,
                 "ideal_scratch_failure_probability": 0,
                 "physical_service_error_certified": False}
    return Readout(blocks[:2], tuple(word), resources, failures, records, error)


def pauli_read(state_matrix, bit, axis="Z"):
    """QND read of one X/Y/Z, with zero denoting its +1 eigenspace."""
    state, work_bits = _state(state_matrix)
    if type(bit) is not int or not 0 <= bit < work_bits:
        raise ValueError("a chosen work qubit required")
    if type(axis) is not str or axis not in ("X", "Y", "Z"):
        raise ValueError("clock axis must be X, Y or Z")
    if axis == "X":
        basis = _hadamard(0)
    elif axis == "Y":
        # B=H S*, so B* Z B=Y. The scalar-corrected phase is S*.
        basis = ([c.pulse("I", -math.pi / 4, [0]),
                  c.pulse("Z", math.pi / 4, [0])] + _hadamard(0))
    else:
        basis = []
    return diagonal_read(state, np.array([1., 0.]), targets=(bit,), basis_word=basis)

