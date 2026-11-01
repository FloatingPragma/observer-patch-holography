"""Exact finite-block composition of the complete native source instrument.

Every block contains literal source pulses on at most three work bits.
Its success and EVERY first-export failure matrix are composed locally.
Only that composition is contracted before applying the block to all
131072 field/pointer rows. No number, eta, or spectral projection occurs.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np

from ..observer.source_scalar import compiler
from ..reference_recovery.recovery import _failure_vectors
from ..reference_observation.observation import frozen


ADMISSION = 2e-8


@dataclass(frozen=True)
class LocalInstrument:
    bits: int
    tape: tuple
    success: np.ndarray
    failures: tuple
    effects: np.ndarray
    complete_error: float


def key(word):
    return tuple((event['kind'], event['angle'], tuple(event['targets'])) for event in word)


@lru_cache(maxsize=128)
def compile_block(tape, bits):
    if type(bits) is not int or not 1 <= bits <= 3: raise ValueError('one through three local work bits required')
    word = [dict(kind=k, angle=a, targets=list(t)) for k, a, t in tape]
    compiler._validate(word, bits)
    if len(tape) > 256: raise ValueError('bounded local native word required')
    d = 1 << bits; state = np.eye(d, dtype=complex); failures = []; effects = []
    for kind, angle, targets in tape:
        full = compiler.native_pulse(kind, angle)
        if full.shape != (6, 6) or np.linalg.norm(full.conj().T@full-np.eye(6)) > 2e-11:
            raise ArithmeticError('actual six-level source lost unitarity admission')
        branches = tuple(frozen(x) for x in _failure_vectors(state, full, targets, bits))
        failures.append(branches)
        effect = np.zeros((4, d, d), complex)
        for slot, branch in enumerate(branches): effect[slot] = branch.conj().T@branch
        effects.append(effect)
        size = 1 << len(targets)
        state = compiler._apply(state, full[:size, :size], targets, bits)
    effects = np.array(effects).reshape(len(tape), 4, d, d)
    error = float(np.linalg.norm(state.conj().T@state+effects.sum(axis=(0, 1))-np.eye(d)))
    if error > 2e-11: raise ArithmeticError('complete local first-failure composition failed')
    return LocalInstrument(bits, tape, frozen(state), tuple(failures), frozen(effects), error)


def arranged(state, targets, bits):
    axes = [bits-1-j for j in reversed(targets)]
    axes += [j for j in range(bits) if j not in axes]+[bits]
    value = state.reshape([2]*bits+[state.shape[1]]).transpose(axes)
    return value.reshape(1 << len(targets), -1, state.shape[1]), axes, value.shape


def apply_success(state, instrument, targets, bits):
    value, axes, shape = arranged(state, targets, bits)
    result = (instrument.success@value.reshape(len(instrument.success), -1)).reshape(shape)
    return result.transpose(np.argsort(axes)).reshape(state.shape)


@dataclass(frozen=True)
class Run:
    bits: int
    blocks: tuple
    initial: np.ndarray
    success: np.ndarray
    failure_weights: np.ndarray
    failure_gram: np.ndarray
    block_offsets: tuple
    normalization_error: float
    rounded_negative_weight_count: int
    rounded_negative_weight_total: float
    failure_sum_error: float
    failure_sum_admission: float

    @property
    def failure_probability(self): return float(np.trace(self.failure_gram).real)

    def failure(self, block, event, slot):
        if type(block) is not int or not 0 <= block < len(self.blocks): raise ValueError('valid block index required')
        instrument, targets = self.blocks[block]
        if (type(event) is not int or not 0 <= event < len(instrument.tape)
                or type(slot) is not int or not 0 <= slot < len(instrument.failures[event])):
            raise ValueError('a genuine retained source outcome, not padding, required')
        state = np.array(self.initial)
        for previous, old_targets in self.blocks[:block]:
            state = apply_success(state, previous, old_targets, self.bits)
        value, _, _ = arranged(state, targets, self.bits)
        failed = instrument.failures[event][slot]@value.reshape(len(instrument.success), -1)
        return failed.reshape(-1, state.shape[1])


def complete(blocks, bits, initial):
    if type(bits) is not int or not 1 <= bits <= 17: raise ValueError('one through seventeen complete work bits required')
    if np.ma.isMaskedArray(initial): raise ValueError('masked inputs forbidden')
    state = np.asarray(initial)
    if state.ndim == 1: state = state[:, None]
    if (state.ndim != 2 or state.shape[0] != 1 << bits or not 1 <= state.shape[1] <= 4
            or state.dtype.kind not in 'fiuc' or not np.isfinite(state).all()):
        raise ValueError('complete bounded field/pointer/reference columns required')
    state = np.asarray(state, complex); initial = frozen(state); state = state.copy()
    blocks = tuple(blocks); weights = []; gram = np.zeros((state.shape[1], state.shape[1]), complex); offsets = []
    offset = 0; rounded_count = 0; rounded_total = 0.
    for instrument, targets in blocks:
        if (type(targets) is not tuple or len(targets) != instrument.bits or len(set(targets)) != len(targets)
                or any(type(j) is not int or not 0 <= j < bits for j in targets)):
            raise ValueError('distinct valid physical targets for every local work bit required')
        value, axes, shape = arranged(state, targets, bits)
        flat = value.reshape(len(instrument.success), -1)
        density = flat@flat.conj().T
        current = np.einsum('pqab,ba->pq', instrument.effects, density).real
        if np.min(current, initial=0) < -1e-40: raise ArithmeticError('negative native failure weight exceeds arithmetic admission')
        # Tiny cancellation below1e-40 is a floating evaluation artifact;
        # every slot remains in the table and the exact local map is retained.
        rounded_count += int(np.count_nonzero(current < 0))
        rounded_total += float(-current[current < 0].sum())
        current = np.maximum(current, 0)
        total_effect = instrument.effects.sum(axis=(0, 1))
        gram += np.einsum('ari,ab,brj->ij', value.conj(), total_effect, value, optimize=True)
        weights.append(current); offsets.append(offset); offset += len(instrument.tape)
        updated = (instrument.success@flat).reshape(shape)
        state = updated.transpose(np.argsort(axes)).reshape(state.shape)
    table = np.concatenate(weights, axis=0) if weights else np.zeros((0, 4))
    error = float(np.linalg.norm(state.conj().T@state+gram-initial.conj().T@initial))
    failure_sum_error = abs(float(table.sum()-np.trace(gram).real))
    failure_sum_admission = float(1e-20+8*np.finfo(float).eps*max(1, offset)*
                                  max(abs(float(table.sum())), abs(float(np.trace(gram).real))))
    if error > ADMISSION*max(1., np.linalg.norm(initial)**2) or failure_sum_error > failure_sum_admission:
        raise ArithmeticError('complete contracted source output/reference accounting failed')
    return Run(bits, blocks, initial, frozen(state), frozen(table), frozen(gram), tuple(offsets), error,
               rounded_count, rounded_total, failure_sum_error, failure_sum_admission)
