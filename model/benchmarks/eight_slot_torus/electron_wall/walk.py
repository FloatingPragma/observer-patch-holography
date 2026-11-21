"""Local charged-wall walk. Array axes after the port are x, y, z, orbital, wall, spin."""
from __future__ import annotations

import numpy as np

PAULI = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
ROOT = np.sqrt([1 / 24] * 12 + [4 / 12, 1 / 12, 1 / 12])
WEIGHT = ROOT[:, None, None, None, None, None, None]


def prepare(psi):
    return WEIGHT * psi


def project(state):
    return (WEIGHT * state).sum(axis=0)


def links(shape):
    return np.ones((3, *shape), dtype=complex)


def route(state, connection, charge, *, spectator=False):
    """Each port sends to its inverse. Reversing an oriented edge takes the adjoint."""
    if state.ndim != 7 or state.shape[0] != 15 or state.shape[-3::2] != (2, 2):
        raise ValueError("expected 15 ports on three public axes, two orbitals and two spins")
    if connection.shape != (3, *state.shape[1:4]) or charge not in (-1, 1):
        raise ValueError("invalid connection or charge")
    if not np.isfinite(connection).all() or np.max(np.abs(np.abs(connection) - 1)) > 1e-12:
        raise ValueError("link phases must have unit modulus")
    out = np.empty(state.shape, dtype=complex)
    for axis in range(3):
        phase = connection[axis] if charge == 1 else connection[axis].conj()
        phase = phase[..., None, None, None]
        for orientation, sign in enumerate((1, -1)):
            port = 2 * axis + orientation
            factor = phase if sign == 1 else np.roll(phase.conj(), 1, axis=axis)
            kinetic = state[port] if spectator else state[port] @ PAULI[axis].T
            kinetic = sign * 1j * kinetic
            kinetic[..., 1, :, :] *= -1
            out[port ^ 1] = np.roll(factor * kinetic, sign, axis=axis)
            wilson = -state[6 + port][..., ::-1, :, :]
            out[6 + (port ^ 1)] = np.roll(factor * wilson, sign, axis=axis)
    out[12] = state[12][..., ::-1, :, :]
    for port, boundary in ((13, 1), (14, -1)):
        source, dest = state[port], out[port]
        dest[..., 0, 1:, :] = -source[..., 1, :-1, :]
        dest[..., 1, :-1, :] = -source[..., 0, 1:, :]
        dest[..., 0, 0, :] = -boundary * source[..., 0, 0, :]
        dest[..., 1, -1, :] = -boundary * source[..., 1, -1, :]
    return out


def step(state, connection, charge, *, spectator=False):
    return route(2 * prepare(project(state)) - state, connection, charge, spectator=spectator)


def dimer(length, a, periodic=False):
    shift = np.eye(length, k=-1)
    if periodic:
        shift[0, -1] = 1
    return a * np.eye(length) - 2 * shift


def symbol(k, length, *, periodic=False, wilson=True):
    """Independent Fourier construction, ordered orbital, private cell, spin."""
    k = np.asarray(k)
    a = 1 + (np.sum(1 - np.cos(k)) if wilson else 0)
    m = dimer(length, a, periodic)
    zero = np.zeros_like(m)
    wall = np.block([[zero, m], [m.T, zero]])
    kinetic = np.einsum("a,aij->ij", np.sin(k), PAULI)
    return (np.kron(wall, np.eye(2)) + np.kron(np.kron(PAULI[2], np.eye(length)), kinetic)) / 12


def wall_mode(length, p, spin, *, spectator=False):
    """Lowest positive wall eigenvector at momentum (p,0,0), spin prepared at rest."""
    u, values, vh = np.linalg.svd(dimer(length, 2 - np.cos(p)))
    left, right = u[:, -1], vh[-1]
    # The singular pair has one arbitrary common sign. Fix it in the two arms.
    phase = np.sign(left[0]) or 1
    plus = phase * np.stack((left, right)) / np.sqrt(2)
    minus = plus.copy()
    minus[1] *= -1
    theta = np.arctan2(np.sin(p), values[-1])
    basis = np.eye(2)[spin]
    rotated = basis if spectator else PAULI[0] @ basis
    vector = (np.cos(theta / 2) * plus[..., None] * basis
              + np.sin(theta / 2) * minus[..., None] * rotated)
    return vector, np.hypot(np.sin(p), values[-1])


def wave(shape, p, vector):
    factors = np.exp(1j * p * np.arange(shape[0]))[:, None, None]
    factors = np.broadcast_to(factors, shape) / np.sqrt(np.prod(shape))
    return factors[..., None, None, None] * vector


def pointer(amplitude):
    """Four phase ports and the complete unmatched part of the signal arm."""
    ports = np.abs(amplitude + np.array([1, 1j, -1, -1j])) ** 2 / 8
    if abs(amplitude) > 1 + 1e-12:
        raise ValueError("noncontractive overlap")
    return np.r_[ports, max(0., 1 - abs(amplitude) ** 2) / 2]


def read_grating(shape, length, charge, strength, spin, *, spectator=False):
    p = 2 * np.pi / shape[0]
    incident, _ = wall_mode(length, -p, spin, spectator=spectator)
    result = []
    connection = links(shape)
    grating = strength * np.cos(2 * p * np.arange(shape[0]))
    connection[1] *= np.exp(1j * grating)[:, None, None]
    signal = project(step(prepare(wave(shape, -p, incident)), connection, charge, spectator=spectator))
    for output_spin in (spin, 1 - spin):
        detector, _ = wall_mode(length, p, output_spin, spectator=spectator)
        result.append(np.vdot(wave(shape, p, detector), signal))
    return np.array(result)


def forecast(period, length, charge, strength, spin, *, spectator=False):
    """Both field parities. The off-diagonal spin port has zero amplitude."""
    p = 2 * np.pi / period
    coords = np.arange(period)
    fourier = np.exp(-2j * p * coords)
    field = strength * np.cos(2 * p * coords)
    # The y hopping changes by -tau_x (cos(A)-1) and -q tau_z sigma_y sin(A).
    incident, _ = wall_mode(length, -p, spin, spectator=spectator)
    detector, _ = wall_mode(length, p, spin, spectator=spectator)
    even = -np.mean((np.cos(field) - 1) * fourier) * np.vdot(detector, incident[::-1])
    _, r = wall_mode(length, p, spin)
    odd = 0j if spectator else -1j * charge * (1 - 2 * spin) * np.sin(p) / r
    return (even + odd * np.mean(np.sin(field) * fourier)) / 12


def rest_properties(length):
    left, values, right = np.linalg.svd(dimer(length, 1.))
    mass = values[-1]
    eta = left[:, -1] @ right[-1]
    denom = np.sqrt(144 - mass ** 2)
    return dict(singular_mass=float(mass), rest_rate=float(np.arcsin(mass / 12)),
                curvature=float((1 / mass + eta) / (2 * denom)),
                magnetic_coefficient=float(1 / (2 * mass * denom)),
                g=float(2 / (1 + mass * eta)),
                wall_endpoint_weight=float((left[0, -1] ** 2 + right[-1, -1] ** 2) / 2))
