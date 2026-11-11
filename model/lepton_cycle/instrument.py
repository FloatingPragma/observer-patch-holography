"""Complete source instruments; every nonreturn source state is retained."""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
import math

import numpy as np

from model.lepton_register_round import carrier, rounds

LABELS = ('yes_yes', 'x_no', 'y_no', 'both_no')
P = F(2, 9)
TRACE_EFFECT = F(5, 18)


@lru_cache(maxsize=1)
def source_data():
    """Read the actual native double-load event and its exact return weights."""
    walk = carrier.carrier_walk()
    records, exact = rounds.records_by_address(walk)
    weights = rounds.complete_record(walk)
    if not exact or len(records) != 60 or weights['outside'] != 0:
        raise ArithmeticError('native event or returned source support changed')
    if weights['overlap'] != weights['norm'] or weights['overlap'] != P:
        raise ArithmeticError('exclusive native read effect changed')
    if tuple(weights[k] for k in LABELS) != (F(4, 81), F(11, 324), F(125, 1296), F(55, 1296)):
        raise ArithmeticError('complete native return weights changed')
    double = rounds.register_branches(walk)
    neighbours = carrier.neighbours()[rounds.PORT]
    other = min(port for port in range(12) if port != rounds.PORT and port not in neighbours)
    pair = rounds.native_pair_branches(walk, (rounds.PORT, other))
    if (double['trace'] != TRACE_EFFECT or pair['trace'] != F(2, 9)
            or double['a2'] / pair['a2'] != F(5, 4)
            or double['b2'] / pair['b2'] != F(5, 4)
            or double['b2_minus'] / pair['b2_minus'] != F(5, 4)):
        raise ArithmeticError('native analyzer normalization comparison changed')
    return records, weights


def reflection(state):
    """Two native uniform-address reflections on axes zero and one."""
    value = np.asarray(state, dtype=np.complex128)
    value = 2 * value.mean(axis=0, keepdims=True) - value
    return 2 * value.mean(axis=1, keepdims=True) - value


def return_parts(state):
    """The four orthogonal outcomes of the two address-consensus tests."""
    value = np.asarray(state, dtype=np.complex128)
    if value.ndim < 2 or value.shape[:2] != (60, 60):
        raise ValueError('two complete sixty-address source registers required')
    x_yes = np.broadcast_to(value.mean(axis=0, keepdims=True), value.shape)
    y_yes = np.broadcast_to(value.mean(axis=1, keepdims=True), value.shape)
    both = np.broadcast_to(value.mean(axis=(0, 1), keepdims=True), value.shape)
    return {'yes_yes': both.copy(), 'x_no': y_yes - both,
            'y_no': x_yes - both, 'both_no': value - x_yes - y_yes + both}


def preparation():
    """Declared fresh blank double at the chosen port, uniform addresses."""
    return np.full((60, 60), 1 / 60, dtype=np.complex128)


def address_preparation_unitary():
    """An explicit declared Householder gate prepares a uniform address."""
    uniform = np.full(60, 1 / math.sqrt(60))
    first = np.zeros(60); first[0] = 1
    vector = first - uniform
    vector /= np.linalg.norm(vector)
    return np.eye(60) - 2 * np.outer(vector, vector)


def _step_size(h):
    if isinstance(h, bool) or not np.isscalar(h) or not np.isreal(h) or not np.isfinite(h):
        raise ValueError('finite real kick angle required')
    return float(h)


def weak_instrument(h):
    """Source amplitudes for each current eigenvalue (-1,0,1), all returns."""
    h = _step_size(h)
    records, _ = source_data()
    f = np.array([r == 2 for r in records])
    g = np.array([r in (2, 3) for r in records])
    exclusive = f[:, None] ^ g[None, :]
    source = preparation()
    x = reflection(exclusive * source)
    phases = np.exp(-1j * h * np.array([-1, 0, 1]))
    output = source[:, :, None] + x[:, :, None] * (phases - 1)
    return return_parts(output)


def analyzer_data(theta):
    """Pull the actual three equal analyzer flags back through the native tick."""
    theta = _step_size(theta)
    records, _ = source_data()
    plus = np.array([r == 2 for r in records])
    minus = np.array([r == 3 for r in records])
    blank = ~(plus | minus)
    masks = (blank[:, None] & (~blank)[None, :],
             plus[:, None] & blank[None, :], minus[:, None] & blank[None, :])
    phi = sum(mask / math.sqrt(int(mask.sum())) for mask in masks) / math.sqrt(3)
    phases = theta + 2 * math.pi * np.arange(3) / 3
    tick = np.ones((60, 60, 3), dtype=np.complex128)
    tick[plus, :, :] = np.exp(1j * phases)[None, None, :]
    tick[minus, :, :] = np.exp(-1j * phases)[None, None, :]
    source = preparation()
    amplitudes = np.sum(phi[:, :, None] * tick * source[:, :, None], axis=(0, 1))
    pulled = reflection(tick.conj() * phi[:, :, None])
    x = pulled * amplitudes
    effects = abs(amplitudes) ** 2
    expected = math.sqrt(5 / 108) + 2 * math.sqrt(5 / 216) * np.cos(phases)
    if np.max(abs(amplitudes - expected)) > 2e-13:
        raise ArithmeticError('native analyzer compression differs from the register')
    return x, effects, amplitudes


def mass_instrument(theta, h):
    """Actual native analyzer kick, inverse tick and four return projections."""
    h = _step_size(h)
    x, _, _ = analyzer_data(theta)
    output = preparation()[:, :, None] + np.expm1(-1j * h) * x
    return return_parts(output)


def gram_parts(parts):
    """Kraus coefficients: G[r]_(k,l)=sum_a v[r]_(a,k) conjugate(v[r]_(a,l))."""
    if set(parts) != set(LABELS):
        raise ValueError('every return outcome is required')
    return {name: block.reshape(-1, 3).T @ block.reshape(-1, 3).conj()
            for name, block in parts.items()}


def apply_with_reference(rho, grams):
    """Keep the complete classical record for a register and arbitrary reference."""
    rho = np.asarray(rho, dtype=np.complex128)
    if rho.ndim != 2 or rho.shape[0] != rho.shape[1] or rho.shape[0] % 3:
        raise ValueError('a square three-mode register/reference matrix required')
    if set(grams) != set(LABELS):
        raise ValueError('every return outcome is required')
    reference = rho.shape[0] // 3
    shaped = rho.reshape(3, reference, 3, reference)
    return {name: (matrix[:, None, :, None] * shaped).reshape(rho.shape)
            for name, matrix in grams.items()}


@dataclass(frozen=True)
class ResetBranch:
    """Factorized SWAP output: fresh active source and complete old-source archive."""
    record: str
    fresh_source: np.ndarray
    archived_source_register: np.ndarray


def reset_branch(record, branch):
    """Swap a fresh prepared pair in; the old pair and its correlations survive."""
    if record not in LABELS or np.asarray(branch).shape[:2] != (60, 60):
        raise ValueError('complete returned source branch required')
    return ResetBranch(record, preparation(), np.array(branch, copy=True))


def unconditional_weak_multiplier(h):
    h = _step_size(h)
    current = np.array([-1, 0, 1])
    difference = current[:, None] - current[None, :]
    return 1 - float(P) + float(P) * np.exp(-1j * h * difference)


def unconditional_mass_multiplier(theta, h):
    h = _step_size(h)
    _, effects, _ = analyzer_data(theta)
    u = np.expm1(-1j * h)
    result = 1 + u * effects[:, None] + u.conjugate() * effects[None, :]
    np.fill_diagonal(result, 1)
    return result
