"""Addressed native pulses that prepare the one-charge interferometer states."""

import numpy as np
from scipy.linalg import expm

from .measure import source


def contact(terms, dimension, a, b):
    diag, hop = np.array(terms) / 8
    h = np.zeros((dimension, dimension))
    h[a, a] = h[b, b] = diag
    h[a, b] = h[b, a] = hop
    return h


def phase_gates(angle):
    terms = source.contacts()
    blank = contact(terms["blank"], 3, 0, 1)
    marked = contact(terms["even"], 3, 1, 2)
    # A change of basis makes both native full transfers unsigned permutations.
    gauge = np.diag([1, 1, -1])
    marked = gauge @ marked @ gauge
    x01 = expm(-1j * np.pi * blank)
    x12 = expm(-1j * np.pi / 8 * marked)
    cycle = x12 @ x01

    def pair_phase(theta):
        t = (-theta / 8) % (2 * np.pi)
        # After the permutation both commuting generators act on the same pair.
        return (
            cycle @ expm(-1j * t / 8 * marked) @ cycle.conj().T @ expm(-1j * t * blank)
        )

    pos, neg = pair_phase(angle / 2), pair_phase(-angle / 2)
    # The three factors have supports 01, 02 and 12 respectively.
    return pos @ (x12 @ pos @ x12.conj().T) @ (cycle.conj().T @ neg @ cycle)


def endpoint_pulse(t):
    # The first physical blank moves both winding endpoints at once.
    term = source.contacts()["blank"]
    left = contact(term, 4, 0, 1)
    right = contact(term, 4, 2, 3)
    free = expm(-1j * t / 2 * (left + right))
    z = np.diag([1, 1, -1, 1])
    correction = np.diag([1, 1, *([np.exp(1j * 8.5 * t)] * 2)])
    return correction @ free @ z @ free @ z


def compile_state(target, terms):
    """Native partial exchanges set all moduli, then the address phases set the arguments."""
    out = np.zeros_like(target, dtype=complex)
    out[0] = 1
    for i in range(len(target) - 1):
        d, h = np.array(terms["blank" if i % 2 == 0 else "even"]) / 8
        angle = np.arctan2(np.linalg.norm(target[i + 1 :]), abs(target[i]))
        t = angle / abs(h)
        u, v = out[i : i + 2]
        c, s = np.cos(h * t), -1j * np.sin(h * t)
        out[i : i + 2] = np.exp(-1j * d * t) * np.array([c * u + s * v, s * u + c * v])
    angle = np.angle(target) - np.angle(out)
    return out * np.exp(1j * angle)
