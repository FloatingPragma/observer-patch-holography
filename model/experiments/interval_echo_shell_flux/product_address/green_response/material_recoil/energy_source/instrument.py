"""Four native exchange links, a local mark and an independent free-flight spectrometer."""

from functools import cache

import numpy as np
from scipy.linalg import expm

from ..common_receiver import receiver

# Binary order; the physical cycle passes through 00, 01, 11, 10.
PORTS = (61, 65, 69, 64)
EDGES = ((0, 1), (1, 3), (3, 2), (2, 0))
BITS = np.array([0, 1, 1, 2])


@cache
def instruments(load):
    board = receiver.original.rails.layout()[0]
    k = receiver.original.rails.candidate.KAPPA
    diag, isolated = receiver.original.terms(load, k)
    if receiver.swap_error(k) > 1e-12:
        raise ArithmeticError("the phase-parking swaps introduce a material phase")

    def energy(pos):
        return (6 - board.degree[pos]) * k * load**2 + sum(
            diag[load, 0] for _ in board.neighbours[pos]
        )

    energies = np.array([energy(p) for p in PORTS])
    flight = np.diag(energies)
    phases = []
    for port in PORTS:
        park = next(v for v in board.neighbours[port] if board.degree[v] == 6)
        # Swap the chosen component out, hold the whole native diagonal, swap back.
        phases.append(np.exp(-1j * np.pi * (energy(park) - energy(port))))
    projectors = []
    mix = expm(-0.5j * np.pi * (isolated - isolated[0, 0] * np.eye(2)))
    for a, b in EDGES:
        if PORTS[b] not in board.neighbours[PORTS[a]]:
            raise ArithmeticError("missing ring seam")
        flight[a, b] = flight[b, a] = isolated[0, 1]
        row = np.zeros(4, complex)
        row[[a, b]] = mix[1] * [1, phases[b]]
        projectors.append(np.outer(row.conj(), row))

    # Drop the common phase after summing all native diagonal terms.
    flight -= energies[0] * np.eye(4)
    prephase = np.diag(np.array(phases) ** BITS)
    sorter = expm(-0.5j * np.pi * flight) @ prephase
    return flight, np.asarray(projectors), sorter


def preparations(load):
    _, _, sorter = instruments(load)
    sharp = np.eye(4, dtype=complex)
    # Sharp inputs to inverse sorting carry a load-dependent overall phase only.
    stocks = {
        name: sorter.conj().T[:, j]
        for name, j in (("ground", 0), ("one", 1), ("top", 3))
    }
    for name, a, b, phase in (
        ("first dimer", 0, 1, -1),
        ("opposite dimer", 3, 2, -1),
        ("quadrature dimer", 0, 1, 1j),
    ):
        stocks[name] = (sharp[a] + phase * sharp[b]) / np.sqrt(2)
    return stocks


def beat_probability(load, flips, time):
    flight, _, _ = instruments(load)
    j = np.arange(4)
    # In binary positions this preparation is uniform on the parity-even
    # ports for the selected bits: a coherent dimer, made by a half exchange
    # and phase parking (with a routed endpoint for the diagonal pair).
    parity = j % 2 if flips == 1 else BITS % 2
    initial = (parity == 0).astype(float) / np.sqrt(2)
    final = expm(-1j * time * flight) @ initial
    return float(np.sum(abs(final[parity == 0]) ** 2))
