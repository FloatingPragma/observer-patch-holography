"""Two-rail energy transport and independent clock/range instruments on one carrier."""

from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from radial_transport.carrier_routes import AXES, POINTS, Transport, coords, grid_edges
from relay import apply, dimer, observations

PROBES = (12, 6, 35, 33)


class Bench(Transport):
    def __init__(self):
        super().__init__()
        self.word[12 * self.clock_owner + 6] = 0
        self.compile(grid_edges(), 3)
        self.slot = max(map(len, self.ops.values()))
        transit = {p // 12 for ops in self.ops.values() for edge in ops for p in edge}
        self.clock_owner = min(set(range(self.board.cells)) - transit)
        self.word[12 * self.clock_owner + 6] = 2
        site, decode, _ = dimer.gates("native")
        self.energy = site[:2, :2]
        self.decoder = decode[:2, :2]
        self.eigenstates = self.decoder.conj().T
        np.testing.assert_allclose(
            self.energy @ self.eigenstates,
            self.eigenstates * [32, 33],
            rtol=0,
            atol=2e-13,
        )
        law = observations.native_law
        idx = [law.SRC.index(p) for p in ((2, 0), (0, 2))]
        parked = observations.KAPPA * law.local_term_parts()[1]
        np.testing.assert_allclose(parked[np.ix_(idx, idx)], 32 * np.eye(2), atol=1e-13)
        self.parked = parked[
            np.ix_(idx, idx)
        ]  # Both energies carry the same idle phase.

    def idle_ticks(self, steps):
        clock = 12 * self.clock_owner + 6
        last = self.read(self.clock_owner, 6)
        for _ in range(steps):
            self.word[clock], self.word[clock + 1] = (
                self.word[clock + 1],
                self.word[clock],
            )
            obs = self.read(self.clock_owner, 6)
            self.service_clock += int(obs != last)
            last = obs

    def volume_returns(self, ticks):
        # Single-step time changes at fixed material endpoints. Read a persistent
        # reflector and the returning packet, then restore the calibration bank.
        returns = []
        for endpoint in AXES:
            self.word[12 * endpoint + 4] = self.word[2] = 2
            start = self.service_clock
            for a, b in ((0, endpoint), (endpoint, 0)):
                self.hop(a, b, 2, metered=True)
                self.idle_ticks(ticks - self.slot)
                assert self.read(b, 2) and not self.read(a, 2)
                if b:
                    assert self.read(b, 4)
            returns.append(self.service_clock - start)
            self.word[12 * endpoint + 4] = self.word[2] = 0
        return returns

    def energy_exchange(self, a, b, cut=False):
        # The exact native pi columns carry a common phase, including spectator
        # occupations. Follow distinct endpoint labels under their full permutations.
        rails = [12 * a + 2, 12 * a + 3, 12 * b + 2, 12 * b + 3]
        ops = self.ops[rails[0], rails[2]] + (
            [] if cut else self.ops[rails[1], rails[3]]
        )
        labels = {p: p for edge in ops for p in edge} | {p: p for p in rails}
        apply(labels, ops)
        assert all(p == label for p, label in labels.items() if p not in rails)
        take = [rails.index(labels[p]) for p in rails]
        # The full 16 sharp words include empty and doubly occupied destinations.
        # Embed the two dimers' four coherent basis columns in that space.
        sector = [5, 9, 6, 10]  # (L0,R0),(L0,R1),(L1,R0),(L1,R1)
        lift = np.eye(16)[:, sector]
        perm = [sum(((x >> take[k]) & 1) << k for k in range(4)) for x in range(16)]
        propagated = np.zeros_like(lift)
        propagated[perm] = lift
        amplitudes = propagated @ np.kron(self.eigenstates, self.eigenstates)
        decoded = np.kron(self.decoder, self.decoder) @ lift.T @ amplitudes
        return decoded, 1 - (abs(decoded) ** 2).sum(axis=0)

    def marker_stroke(self):
        changes = np.zeros(len(POINTS), dtype=int)
        for v in POINTS:
            before = self.read(v, 8)
            a = 12 * v + 8
            self.word[a], self.word[a + 1] = self.word[a + 1], self.word[a]
            changes[v] = int(before != self.read(v, 8))
        return changes

    def frame(self, axis, bank, live=True):
        stride = AXES[axis]
        pairs = [(a, a + stride) for a in POINTS if coords(a)[axis] in (0, 2)]
        assert len(pairs) == 50 and len({v for e in pairs for v in e}) == 100
        begin = self.service_clock
        old = bank.copy()
        changes = np.zeros(len(POINTS), dtype=int)
        completed = 0
        for a, b in pairs:
            decoded, leak = self.energy_exchange(a, b)
            np.testing.assert_allclose(
                decoded, np.eye(4)[:, [0, 2, 1, 3]], rtol=0, atol=4e-13
            )
            np.testing.assert_allclose(leak, 0, rtol=0, atol=5e-13)
            for port in (2, 3):
                self.hop(a, b, port, cut=not live, metered=True)
                completed += 1
                if completed % 25 == 0:
                    changes += self.marker_stroke()
            if live:
                # The decoded quantum channel swaps complete registers even with
                # a reference, so each column may challenge a different input.
                bank[[a, b]] = old[[b, a]]
        return self.service_clock - begin, changes

    def spectroscopy(self, active, total):
        # Remove the common parked phase analytically before exponentiating.
        gap = self.energy - self.parked
        u = expm(-1j * active * gap)
        ref = expm(-0.5j * np.pi * gap)
        position = np.eye(2)[:, 0]
        bright = [
            float(abs((u @ prep @ position)[1]) ** 2) for prep in (np.eye(2), ref)
        ]
        phase = float(np.arctan2(2 * bright[1] - 1, 1 - 2 * bright[0]))
        measured = abs(self.decoder @ u @ self.eigenstates) ** 2
        np.testing.assert_allclose(measured, np.eye(2), atol=3e-13, rtol=0)
        start = self.service_clock
        self.idle_ticks(total)
        elapsed = self.service_clock - start
        assert elapsed == total
        return (
            np.round(bright, 12).tolist(),
            phase / (np.pi * elapsed),
            measured.tolist(),
        )
