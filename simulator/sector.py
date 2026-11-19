"""Exact evolution inside one private Fourier sector, and the scalar field reading.

A state of the form psi(x, p) = exp(2 pi i m p / P) chi(x) keeps its private
order m under the law, because every update is invariant under translation
of the private register. Routing a slot one position along the private axis
then multiplies chi by exp(-i s theta) with theta = 2 pi m / P and s the
slot's direction. ``sector_step`` executes exactly this, so a carrier with a
private register of any size P costs the memory of the three spatial axes.
``test_simulator.py`` compares it with the literal four-axis update.
"""
from __future__ import annotations

import numpy as np

from .law import EQUAL, PARTNER, SLOT_AXIS, SLOT_SIGN, SLOTS, Law, coin


def sector_step(chi: np.ndarray, law: Law, private_angle: float) -> np.ndarray:
    """One update of chi with shape (N1, N2, N3, 8) in the private sector of angle theta."""
    if chi.ndim != 4 or chi.shape[-1] != SLOTS:
        raise ValueError("a sector state has shape (N1, N2, N3, 8)")
    coined = coin(chi, law.equal_phase)
    out = np.empty_like(coined)
    for slot in range(SLOTS):
        axis, sign = int(SLOT_AXIS[slot]), int(SLOT_SIGN[slot])
        source = coined[..., slot]
        if axis < 3:
            moved = np.roll(source, sign, axis=axis)
        else:
            moved = source * np.exp(-1j * sign * private_angle)
        out[..., PARTNER[slot]] = moved
    return out


def field(chi: np.ndarray) -> np.ndarray:
    """The equal-slot projection <e|chi(x)>, the scalar field the observers' interferometers read."""
    return chi @ EQUAL.conj()


def prepare_field(initial: np.ndarray) -> np.ndarray:
    """Release the real field ``initial`` at rest: chi = f(x) e, so that <e|chi> = f."""
    return initial[..., None] * EQUAL


def plane_modes(size: int, modes) -> np.ndarray:
    """Sum of amplitude * cos(2 pi K.x) on a chart of side one sampled at N points per axis."""
    x = np.arange(size) / size
    grid = np.stack(np.meshgrid(x, x, x, indexing="ij"), axis=-1)
    total = np.zeros((size,) * 3)
    for amplitude, wave in modes:
        total += amplitude * np.cos(2 * np.pi * grid @ np.asarray(wave, float))
    return total


def continuum_field(modes, private_order: int, chart_times, points) -> np.ndarray:
    """Klein-Gordon solution with c = 1/2 and rest frequency |Q|/2, released at rest.

    phi(t, x) = sum_j a_j cos(omega_j t) cos(2 pi K_j.x), omega_j = sqrt(|2 pi K_j|^2 + Q^2) / 2,
    Q = 2 pi * private_order. ``points`` are chart coordinates with shape (..., 3).
    """
    q = 2 * np.pi * private_order
    t = np.asarray(chart_times, float)[:, None]
    points = np.asarray(points, float)
    total = 0.0
    for amplitude, wave in modes:
        wave = np.asarray(wave, float)
        omega = 0.5 * np.hypot(2 * np.pi * np.linalg.norm(wave), q)
        spatial = np.cos(2 * np.pi * points @ wave).reshape(1, -1)
        total = total + amplitude * np.cos(omega * t) * spatial
    return np.asarray(total).reshape((len(t),) + points.shape[:-1])
