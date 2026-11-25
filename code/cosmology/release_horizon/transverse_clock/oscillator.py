"""Two native exchange contacts per edge implement the Fourier incidence lift."""

from functools import cache
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from angular_transfer.rays import CELL, PORTS, Transport, align, apply, frame
from angular_transfer.polarimeter import Polarimeter

COUNT = 30
COUPLING = 1 / np.sqrt(2 * COUNT)


class Instrument:
    def __init__(self):
        self.pol = Polarimeter()
        self.net = Transport()
        _, words, directions = frame()
        self.vertices = directions[:12]
        self.edges = sorted(tuple(sorted(e)) for e in CELL.edges)
        self.deltas = np.array(
            [self.vertices[j] - self.vertices[i] for i, j in self.edges]
        )
        assert len(self.edges) == COUNT and all(len(w) == 1 for w in words[:12])
        # The central rail at port 3 couples to a borrowed edge rail at port 2.
        # The vertex of the second copy is on a different physical owner.
        self.hubs = (0, 32)
        self.routes = []
        for hub in self.hubs:
            pairs = [(hub, hub + i + 1) for i in range(COUNT)]
            self.net.compile(pairs, PORTS[0])
            sizes = []
            for home, edge in pairs:
                a, b = 12 * home + PORTS[0], 12 * edge + PORTS[0]
                ops = self.net.ops[a, b]
                labels = {p: p for gate in ops for p in gate}
                apply(labels, ops)
                assert labels[a] == b and labels[b] == a
                assert all(labels[p] == p for p in labels if p not in (a, b))
                sizes.append(len(ops))
            self.routes.append(sizes)
        for owner in (
            *self.hubs,
            *(h + i + 1 for h in self.hubs for i in range(COUNT)),
        ):
            ports = PORTS if owner in self.hubs else PORTS[:1]
            for p in ports:
                self.net.read_ops[p] = align([p, 11], [p, min(CELL[p])])
                a = 12 * owner + p
                self.net.word[a] = 2
                assert self.net.read(owner, p) == 1
                self.net.word[a] = 0
                assert self.net.read(owner, p) == 0

    def joint_polarimeter(self):
        # Bring the second vertex to the scratch rail beside the first vertex.
        # Read diagonal and circular polarizations coherently on this pair.
        hub, other = self.hubs
        self.net.compile([(hub, other)], PORTS[0])
        route = [(12 * other + 2, 12 * other + 3)] + self.net.ops[
            12 * other + 2, 12 * hub + 2
        ]
        assert 3 in CELL[2]
        labels = list(range(self.net.board.n))
        apply(labels, route)
        assert [labels[12 * hub + 3], labels[12 * hub + 2]] == [
            12 * hub + 3,
            12 * other + 3,
        ]
        apply(labels, route[::-1])
        assert labels == list(range(self.net.board.n))
        return len(route)

    def clock(self, steps):
        owner = self.net.clock_owner
        assert owner not in range(64)
        addr = 12 * owner + 6
        word = self.net.word
        before = self.net.read(owner, 6)
        marks = [before]
        for _ in range(steps):
            a, b = self.net.flip
            word[addr : addr + 2] = [word[addr + a], word[addr + b]]
            marks.append(self.net.read(owner, 6))
        assert all(a != b for a, b in zip(marks, marks[1:]))
        return marks

    @cache
    def bare_contact(self, angle):
        # Native Givens rotation on [scratch, vertex] transformed to
        # exp(-i angle X); public coordinates below are [vertex, scratch].
        pol = self.pol
        shift = pol.phase(1, -np.pi / 2)
        g = shift @ pol.givens(0, 1, angle) @ shift.conj().T
        return g[:2, :2][::-1, ::-1]

    def sweep(self, phases, dt, broken=False):
        pairs = []
        g = self.bare_contact(COUPLING * dt / 2)
        for i, phi in enumerate(phases, start=1):
            for angle in (0.0, phi + (0 if broken else np.pi)):
                # A parked native phase on the borrowed rail precedes its return.
                d = self.pol.phase(0, angle)[0, 0]
                pair = g.copy()
                pair[1, 0] *= d
                pair[0, 1] *= d.conjugate()
                pairs.append((i, pair))
        out = np.eye(COUNT + 1, dtype=complex)
        for i, pair in pairs + pairs[::-1]:
            out[[0, i]] = pair @ out[[0, i]]
        return out


def generator(phases, broken=False):
    b = COUPLING * (1 + (1 if broken else -1) * np.exp(1j * phases))
    h = np.zeros((COUNT + 1, COUNT + 1), complex)
    h[1:, 0], h[0, 1:] = b, b.conj()
    return h


def detector(amplitude, pol):
    # Source and parked reference each have half the original weight.
    # Interference reads both ports and the outside event in either basis.
    pair = np.array([1, amplitude], complex) / np.sqrt(2)
    outcomes = []
    for phase in (0.0, -np.pi / 2):
        u = pol.givens(0, 1, np.pi / 4).conj().T @ pol.phase(1, phase)
        hit = abs(u[:2, :2] @ pair) ** 2
        outcomes.append([*hit, max(0.0, 1 - sum(hit))])
    return np.array(outcomes)


def rounded(x):
    return (np.round(np.asarray(x), 9) + 0.0).tolist()


def oscillator(rng):
    inst = Instrument()
    v = inst.vertices[0]
    oblique = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
    marker = inst.clock(3)
    joint_route = inst.joint_polarimeter()
    records = []
    shot = 2**20
    for label, k, broken in [
        ("zero", 0 * v, False),
        ("port", 1.7 * v, False),
        ("oblique", 2.4 * oblique, False),
        ("wrong_sign", 0 * v, True),
    ]:
        phases = inst.deltas @ k
        h = generator(phases, broken)
        lam = (1 - np.cos(phases).mean()) if not broken else (1 + np.cos(phases).mean())
        np.testing.assert_allclose((h @ h)[0, 0], lam, atol=2e-14)
        # The sign of velocity reverses on independently prepared imaginary
        # edge launches. Pick a non-dark edge from its geometrical phase alone.
        edge = int(np.argmax(abs(h[0, 1:]))) + 1
        # Check the actual analyzer for all three known reference phases.
        for a in (1.0, -1.0, 1j):
            z = detector(a, inst.pol)
            np.testing.assert_allclose(z[:, 0] - z[:, 1], [a.real, a.imag], atol=3e-13)
        refines = []
        for n in (8, 32, 128):
            dt = 1 / n
            step = inst.sweep(phases, dt, broken)
            per_tick = np.linalg.matrix_power(step, n)
            exact = expm(-1j * h)
            one_step = np.linalg.norm(step - expm(-1j * dt * h), 2)
            # Telescoping of exact unitaries gives this sharper numerical
            # estimate. The uniform analytic bound is independent of spectra.
            err = np.linalg.norm(per_tick - exact, 2)
            assert err <= n * one_step + 2e-11
            refines.append([n, *rounded([err, n * one_step])])
        rows = []
        for t in (1, 2, 3):
            unitary = np.linalg.matrix_power(per_tick, t)
            for preparation, col, amp in [
                ("vertex", 0, 1),
                ("edge_plus", edge, 1j),
                ("edge_minus", edge, -1j),
            ]:
                a = amp * unitary[0, col]
                expected = (
                    np.cos(t * np.sqrt(lam))
                    if col == 0
                    else -1j * t * np.sinc(t * np.sqrt(lam) / np.pi) * h[0, col] * amp
                )
                assert abs(a - expected) <= t * err + 2e-11
                probs = detector(a, inst.pol)
                # Independent repeats on the two disjoint rail banks.
                shots = [
                    [
                        rng.multinomial(shot, p / p.sum()).tolist()
                        for p in np.round(probs, 11)
                    ]
                    for _ in inst.hubs
                ]
                rows.append(
                    dict(
                        t=t,
                        input=preparation,
                        response=rounded([a.real, a.imag]),
                        counts=shots,
                    )
                )
        # Both transverse amplitudes see identical schedules on disjoint rails.
        # Route and recombine the surviving vertices before a sharp read.
        axis = k / np.linalg.norm(k) if np.any(k) else v
        basis = inst.pol.frame(axis)[:, :2]
        np.testing.assert_allclose(basis.conj().T @ axis, 0, atol=2e-13)
        vertex = per_tick[0, 0]
        joined = []
        for analyzer in inst.pol.axes[1:]:
            emit = analyzer[:2, 0]
            for sign in (1, -1):
                hits = abs(analyzer[:2, :2].conj().T @ (vertex * emit * [1, sign])) ** 2
                np.testing.assert_allclose(
                    hits,
                    abs(vertex) ** 2 * np.array([sign == 1, sign == -1]),
                    atol=2e-12,
                )
                joined.append(rounded([*hits, max(0.0, 1 - sum(hits))]))
        records.append(
            dict(
                mode=label,
                coherent_polarization=joined,
                wavevector=rounded(k),
                stiffness=round(float(lam), 11),
                edge=edge - 1,
                refinements=refines,
                reads=rows,
            )
        )
    return dict(
        edges=inst.edges,
        native_vertices=rounded(inst.vertices),
        rail_hubs=inst.hubs,
        route_lengths=inst.routes,
        polarization_route=joint_route,
        marker_owner=inst.net.clock_owner,
        marker=marker,
        shots=shot,
        modes=records,
    )
