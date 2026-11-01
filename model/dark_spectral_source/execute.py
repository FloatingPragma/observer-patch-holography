"""Bounded13-bit source execution with every first-export failure retained."""
from types import MappingProxyType
import math

import numpy as np

from ..observer.source_scalar import compiler
from ..reference_observation.observation import frozen
from ..reference_recovery.recovery import CompleteRun, _failure_vectors


ADMISSION = 3e-9
MAX_PULSES = 32_000


def complete(word, bits, initial):
    """Same M6/code export as the earlier9-bit evaluator; new bounded domain.

    Every literal pulse acts on its full six-level source before exporting.
    The inherited CompleteRun stores immutable source matrices and prefix
    tape, so all rectangular first-failure vectors remain reconstructible.
    """
    if type(bits) is not int or bits not in (12, 13):
        raise ValueError('twelve field bits with at most one explicit pointer required')
    if np.ma.isMaskedArray(initial):
        raise ValueError('masked input forbidden')
    a = np.asarray(initial)
    if a.ndim == 1:
        a = a[:, None]
    if (a.ndim != 2 or a.shape[0] != 1 << bits or not 1 <= a.shape[1] <= 4
            or a.dtype.kind not in 'fiuc' or not np.isfinite(a).all()):
        raise ValueError('complete bounded field/reference columns required')
    a = np.asarray(a, complex)
    raw = compiler._validate(word, bits)
    if len(raw) > MAX_PULSES:
        raise ValueError('native tape exceeds the declared bounded execution domain')
    tape = tuple((e['kind'], e['angle'], tuple(e['targets'])) for e in raw)
    matrices = {}
    for kind, angle, _ in tape:
        if (kind, angle) not in matrices:
            value = compiler.native_pulse(kind, angle)
            if value.shape != (6, 6) or np.linalg.norm(value.conj().T@value-np.eye(6)) > 2e-11:
                raise ArithmeticError('actual full source pulse failed unitarity admission')
            matrices[kind, angle] = frozen(value)
    state = a.copy()
    weights = np.zeros((len(tape), 4))
    gram = np.zeros((a.shape[1], a.shape[1]), complex)
    for event, (kind, angle, targets) in enumerate(tape):
        full = matrices[kind, angle]
        for slot, branch in enumerate(_failure_vectors(state, full, targets, bits)):
            gram += branch.conj().T@branch
            weights[event, slot] = np.vdot(branch, branch).real
        size = 1 << len(targets)
        state = compiler._apply(state, full[:size, :size], targets, bits)
    error = np.linalg.norm(state.conj().T@state+gram-a.conj().T@a)
    if (not np.isfinite(state).all() or not math.isfinite(error)
            or error > ADMISSION*max(1., np.linalg.norm(a)**2)
            or abs(weights.sum()-np.trace(gram).real) > ADMISSION):
        raise ArithmeticError('complete source output/reference accounting failed')
    return CompleteRun(tape, bits, frozen(a), MappingProxyType(matrices), frozen(state),
                       frozen(weights), frozen(gram))


def pointer_outputs(state):
    """Actual incoherent M6 export; both pointer values and four aborts."""
    if state.shape[0] != 8192:
        raise ValueError('complete thirteen-bit state required')
    kraus = compiler.transfer_kraus(2, coherent=False)
    if len(kraus) != 6 or not np.array_equal(sum(k.conj().T@k for k in kraus), np.eye(6)):
        raise ArithmeticError('pointer export omitted an outcome')
    arranged = state.reshape(2, -1)
    outputs = tuple(k[:, :2]@arranged for k in kraus)
    return outputs
