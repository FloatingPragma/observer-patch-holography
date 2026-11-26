"""Ternary reference echo and free vacancy-edge release on the registered guides."""

from collections import deque
from functools import cache
from itertools import product

import numpy as np
from scipy.linalg import expm

from .. import receiver as original

TRIT_TIME = 4 * np.pi / 3
RELEASE = np.pi / np.hypot(1, 0.5)
WORDS = tuple(product((0, 1), repeat=3))


@cache
def apparatus(equal_degree=False):
    board, guides, flag = original.rails.layout()
    occupied = set(flag).union(*guides)
    for left in range(board.n):
        for right in sorted(map(int, board.neighbours[left])):
            if board.degree[left] != 6 or board.degree[right] != (
                6 if equal_degree else 5
            ):
                continue
            around = {left, right, *board.neighbours[left], *board.neighbours[right]}
            if occupied & around:
                continue
            for pad in sorted(
                set(board.neighbours[left]) - {right} - set(board.neighbours[right])
            ):
                pad = int(pad)
                paths = {pad: (pad,)}
                queue = deque([pad])
                while queue:
                    pos = queue.popleft()
                    for dest in map(int, board.neighbours[pos]):
                        if dest not in paths and dest not in (left, right):
                            paths[dest] = (*paths[pos], dest)
                            queue.append(dest)
                outside = sorted(paths.keys() - around - occupied)
                if not outside or not set().union(*guides) <= paths.keys():
                    continue
                reference = outside[0]
                routes = [paths[p] for p in (*sum(guides, ()), reference)]
                for route in routes:
                    steps = list(zip(route, route[1:]))
                    labels = list(range(board.n))
                    for a, b in steps + steps[-2::-1]:
                        labels[a], labels[b] = labels[b], labels[a]
                    target = list(range(board.n))
                    a, b = route[0], route[-1]
                    target[a], target[b] = target[b], target[a]
                    if labels != target:
                        raise ArithmeticError("a pad route moved an unrequested label")
                return left, right, pad, reference, tuple(routes)
    raise ArithmeticError("no independent blank release edge")


@cache
def native_terms(probe, stiffness, equal_degree=False):
    """Sum the native onsite and incident terms before taking the two-rail restriction."""
    board, _, _ = original.rails.layout()
    left, right, pad, _, _ = apparatus(equal_degree)
    pair = {left, right}
    seams = [(a, b) for a, b in board.seams if a in pair or b in pair]
    diagonal, kinetic = original.terms(probe, stiffness)
    energies = []
    for pad_load in range(3):
        readings = []
        for pos in (left, right):
            loads = {pos: probe, pad: pad_load}
            energy = stiffness * (6 - board.degree[pos]) * probe**2
            energy += sum(diagonal[loads.get(a, 0), loads.get(b, 0)] for a, b in seams)
            readings.append(energy)
        energies.append(readings)
    flight = np.diag(energies[0]) + np.array([[0, kinetic[0, 1]], [kinetic[1, 0], 0]])
    return np.asarray(energies), flight


@cache
def echo(probe, stiffness, sense, reference, trit, equal_degree):
    """In the original rail order the second hold takes the two reversed energies."""
    energies, _ = native_terms(probe, stiffness, equal_degree)
    sums = energies[sense] + energies[reference, ::-1]
    return np.exp(-1j * (trit * TRIT_TIME) * (sums - sums[0]))


@cache
def standing_wave(probe, stiffness, sign):
    """Prepare opposite real coherences on the blank unequal-degree edge."""
    _, seam = original.terms(probe, stiffness)
    psi = expm(-0.5j * np.pi * (seam - seam[0, 0] * np.eye(2)))[:, 0]
    blank = native_terms(probe, stiffness)[1].diagonal()
    phase_time = (3 if sign == 1 else 1) * np.pi
    return psi * np.exp(-1j * phase_time * (blank - blank[0]))


@cache
def right_probability(
    probe, stiffness, flags, source, blind, sign, dwell, parked, equal_degree
):
    _, flight = native_terms(probe, stiffness, equal_degree)
    if parked:
        flight = np.diag(flight.diagonal())
    unitary = expm(-1j * dwell * (flight - np.eye(2) * flight.trace() / 2))
    # The averaged single-coordinate density channels compose because the trits
    # and reference draws are independent. Keep both populations for every trial.
    psi = standing_wave(probe, stiffness, sign)
    rho = np.outer(psi, psi.conj())
    for matched, tag in zip(flags, source, strict=True):
        sense = tag if matched else 0
        standards = (1, 2) if blind else (tag or 1,)
        states = [
            echo(probe, stiffness, sense, ref, k, equal_degree)
            for ref in standards
            for k in range(3)
        ]
        rho *= sum(np.outer(v, v.conj()) for v in states) / len(states)
    prob = (unitary @ rho @ unitary.conj().T)[1, 1].real
    return float(prob)


def swap_error(stiffness):
    law = original.rails.candidate.native_law
    exchange, contrast = law.local_term_parts()
    generator = exchange + stiffness * contrast
    src = law.SRC
    whole = expm(-1j * np.pi * generator)
    errors = []
    for n, (a, b) in enumerate(src):
        desired = src.index((b, a))
        errors.append(np.max(abs(whole[:, n] - np.eye(len(whole))[:, desired])))
    return float(max(errors))
