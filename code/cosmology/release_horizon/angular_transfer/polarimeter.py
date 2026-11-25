"""Three-amplitude native seam analyzer with a nondestructive rejected rail."""

import numpy as np
from scipy.linalg import expm

if __package__:
    from .rays import PORTS
else:
    from rays import PORTS
from relay import CELL, align, observations


class Polarimeter:
    def __init__(self):
        law = observations.native_law
        exchange, contrast = law.local_term_parts()
        self.blocks = {}
        for a, b in ((0, 1), (0, 2)):
            assert PORTS[b] in CELL[PORTS[a]]
            loads = 2 * np.eye(3, dtype=int)
            idx = [law.SRC.index((v[a], v[b])) for v in loads]
            h = exchange + observations.KAPPA * contrast
            np.testing.assert_allclose(np.delete(h[:, idx], idx, axis=0), 0, atol=2e-13)
            self.blocks[a, b] = h[np.ix_(idx, idx)]
        self.rate = float(
            observations.KAPPA * contrast[law.SRC.index((2, 0)), law.SRC.index((2, 0))]
        )
        np.testing.assert_allclose(self.rate, 32, atol=2e-13, rtol=0)
        assert all(set(CELL[p]) - set(PORTS) for p in PORTS)
        parked = observations.KAPPA * contrast
        for x, want in ((0, 0), (2, self.rate)):
            col = law.SRC.index((x, 0))
            np.testing.assert_allclose(
                parked[:, col], np.eye(law.DIM)[:, col] * want, rtol=0, atol=2e-13
            )
        diagonal = self.givens(0, 1, np.pi / 4)
        self.axes = [np.eye(3), diagonal, self.phase(1, np.pi / 2) @ diagonal]
        self.launch = np.column_stack(
            (np.eye(3)[:, :2], self.axes[1][:, 0], self.axes[2][:, 0])
        )

    def phase(self, rail, angle):
        # Park the addressed rail against a local blank. Occupancy two costs
        # 32; every other component has zero load on this seam.
        t = (-angle) % (2 * np.pi) / self.rate
        gate = np.eye(3, dtype=complex)
        gate[rail, rail] = np.exp(-1j * self.rate * t)
        return gate

    def givens(self, a, b, theta):
        t = (2 * theta) % (4 * np.pi)
        h = self.blocks[a, b]
        scalar = float(h[a, a].real) * t
        correction = self.phase(a, scalar) @ self.phase(b, scalar)
        gate = (
            self.phase(b, -np.pi / 2)
            @ correction
            @ expm(-1j * t * h)
            @ self.phase(b, np.pi / 2)
        )
        rotation = np.eye(3)
        rotation[np.ix_([a, b], [a, b])] = [
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]
        np.testing.assert_allclose(gate, rotation, atol=2e-13, rtol=0)
        return gate

    def frame(self, n):
        basis = self.givens(0, 1, np.arctan2(n[1], n[0])) @ self.givens(
            0, 2, -np.arccos(np.clip(n[2], -1, 1))
        )
        np.testing.assert_allclose(basis[:, 2], n, atol=2e-13, rtol=0)
        return basis

    def reads(self, frames, net, paths):
        # The verified full port permutations preserve the ordered three rails
        # with common phase. Bind the Born effects to real arrival pointers.
        for rail in PORTS:
            net.read_ops[rail] = align([rail, 11], [rail, min(CELL[rail])])
        for owner in {0, *(b for _, b in paths)}:
            for rail in PORTS:
                net.word[12 * owner + rail] = 2
                assert [net.read(owner, p) for p in PORTS] == [
                    int(p == rail) for p in PORTS
                ]
                net.word[12 * owner + rail] = 0
        records = np.empty((len(frames), len(frames), 4, 3, 3))
        for s, emit in enumerate(frames):
            state = emit @ self.launch
            for t, analyze in enumerate(frames):
                out = analyze.conj().T @ state
                # One QND normal-rail query. Both transverse components have
                # the identical dark local input and therefore retain coherence.
                rejected = abs(out[2]) ** 2
                kept = out.copy()
                kept[2] = 0
                for k, basis in enumerate(self.axes):
                    branch = abs(basis.conj().T @ kept) ** 2
                    records[s, t, :, k] = np.column_stack(
                        (branch[0], branch[1], rejected)
                    )
        np.testing.assert_allclose(records.sum(axis=-1), 1, rtol=0, atol=4e-13)
        return records
