"""Finite electric register on the three directed x-links of a periodic ring."""
from itertools import product

import numpy as np

from .walk import links, prepare, step


def flux_register(charge, modulus=5, services=3):
    """Fourier sectors of three link shifts, transformed back before reading flux."""
    if modulus < 3 or services < 1 or charge not in (-1, 1):
        raise ValueError("test requires distinguishable charges and a positive service count")
    shape = (3, 1, 1)
    seed = np.zeros((*shape, 2, 2, 2), complex)
    seed[0, 0, 0, 0, 0, 0] = 1
    seed = prepare(seed)
    amplitudes = np.empty((modulus,) * 3 + (3, 15, 2, 2, 2), complex)
    for sectors in product(range(modulus), repeat=3):
        conn = links(shape)
        conn[0, :, 0, 0] = np.exp(2j * np.pi * np.array(sectors) / modulus)
        state = seed
        for _ in range(services):
            state = step(state, conn, charge)
        amplitudes[sectors] = state[:, :, 0, 0].swapaxes(0, 1)
    amplitudes = np.fft.fftn(amplitudes, axes=(0, 1, 2)) / modulus ** 3
    probabilities = (abs(amplitudes) ** 2).sum(axis=(4, 5, 6, 7))
    bad = 0.
    for field in product(range(modulus), repeat=3):
        divergence = np.array(field) - np.roll(field, 1)
        for position in range(3):
            density = np.zeros(3, int)
            density[position] += charge
            density[0] -= charge
            if np.any((divergence + density) % modulus):
                bad += probabilities[(*field, position)]
    displaced = probabilities[..., 1:].sum()
    if bad > 1e-20 or abs(probabilities.sum() - 1) > 1e-12 or displaced < .001:
        raise ArithmeticError("Gauss test failed or excitation did not move")
    return probabilities, float(bad)
