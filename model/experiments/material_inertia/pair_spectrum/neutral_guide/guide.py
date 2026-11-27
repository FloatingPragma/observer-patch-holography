"""Blank-seam physical rings and complete relative-coordinate propagation."""

from fractions import Fraction

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model.carrier import federation
from ..native_source_probe import circuit, measure as source, physical
from model.experiments.interval_echo_consensus.signed_effects import fixed_pointer as pointer


def first_path(neighbours, start, stop, edges, allowed):
    def visit(path):
        if len(path) == edges + 1:
            return path if path[-1] == stop else None
        for other in neighbours[path[-1]]:
            if other in allowed and other not in path and (other != stop or len(path) == edges):
                result = visit(path + [other])
                if result:
                    return result
        return None
    return visit([start])


def physical_rings():
    board = federation(0)
    bridge = {(a // 12, b // 12): (a, b) for a, b in board.external}
    bridge.update({(j, i): (b, a) for (i, j), (a, b) in list(bridge.items())})
    dual = [sorted(j for i, j in bridge if i == owner) for owner in range(board.cells)]
    candidates = []
    for first in range(board.cells):
        for last in dual[first]:
            path = first_path(dual, first, last, 4, range(board.cells))
            if path:
                candidates.append(path)
    owners = min(candidates)
    apparatus = source.routes(board, owners)
    paths = dict(zip(apparatus['target_ports'], apparatus['endpoint_transposition_paths']))
    rings = {}
    for steps in (3, 5, 7):
        ring = []
        for i, owner in enumerate(owners):
            entry = bridge[owners[i - 1], owner][1]
            exit_port = bridge[owner, owners[(i + 1) % len(owners)]][0]
            path = first_path(board.neighbours, entry, exit_port, steps, range(12 * owner, 12 * owner + 12))
            if not path:
                raise ValueError(f'Missing physical path of length {steps} in owner {owner}')
            ring.extend(path)
        if len(set(ring)) != len(ring) or any(b not in board.neighbours[a] for a, b in zip(ring, ring[1:] + ring[:1])):
            raise ValueError('physical ring is not simple and closed')
        swaps = [physical.validate_route(board, paths[p], apparatus['source_park'], apparatus['source_cross']) for p in ring]
        rings[len(ring)] = dict(ports=ring, ramp_routing_contacts_per_pass=2 * sum(swaps),
                                target_to_station_paths=[paths[p] for p in ring])
    station = {k: apparatus[k] for k in ('interaction_port', 'source_cross', 'source_park')}
    return owners, rings, station


def relative_hamiltonian(size, character, kappa):
    """Assemble the Fourier block by acting on every unlike native contact of a representative word."""
    k = 2 * np.pi * character / size
    h = np.zeros((size - 1, size - 1), complex)
    local = law.seam_term(kappa)
    for r in range(1, size):
        loads = np.ones(size, int)
        loads[0], loads[r] = 0, 2
        for a in sorted({0, size - 1, r, r - 1}):
            b = (a + 1) % size
            n, m = loads[[a, b]]
            src = 3 * n + m
            h[r - 1, r - 1] += local[src, src]
            if n != m:
                exchanged = {a: b, b: a}
                x, y = exchanged.get(0, 0), exchanged.get(r, r)
                h[(y - x) % size - 1, r - 1] += local[3 * m + n, src] * np.exp(-1j * k * x)
    np.testing.assert_allclose(h, h.conj().T, atol=3e-14, rtol=0)
    return h


def source_vector(size, character, orientation):
    v = np.zeros(size - 1, complex)
    v[0], v[-1] = 1, orientation * np.exp(-2j * np.pi * character / size)
    return v / np.sqrt(2)


class Flight:
    def __init__(self, size, character, kappa, orientation):
        self.initial = source_vector(size, character, orientation)
        self.h = relative_hamiltonian(size, character, kappa)
        self.energy, self.basis = eigh(self.h)
        self.coefficients = self.basis.conj().T @ self.initial

    def at(self, time):
        return self.basis @ (np.exp(-1j * time * self.energy) * self.coefficients)


def signed_pulse_audit(size):
    """Check the unit ramp's action on all three loads, including leakage into the other 24 words."""
    results = {}
    for kind, denominator, second in [('signed', 65, -1), ('unsigned', 31, 1)]:
        dwell = Fraction(-4, size * denominator)
        # Cross initially zero. A full material swap reverses the two stock rails.
        swap = [pointer.pulse('contact', 1, 2, 1)]
        word = pointer.diagonal_word(0, 2, dwell % 4) + swap
        word += pointer.diagonal_word(0, 2, (second * dwell) % 4) + swap
        indices = [idx[0] for idx in circuit.CODE]
        actual = circuit.apply(word, np.eye(27)[:, indices])
        expected = np.zeros_like(actual)
        q = np.arange(-1, 2)
        phases = np.exp(2j * np.pi / size * (q if kind == 'signed' else q*q))
        expected[indices, range(3)] = phases * actual[indices[1], 1]
        err = float(np.max(abs(actual - expected)))
        if err > 1e-11:
            raise ArithmeticError('native impulse fails ternary phase or returns wrong rails')
        results[kind] = dict(pulses=len(word), max_column_entry_error=err)
    return results
