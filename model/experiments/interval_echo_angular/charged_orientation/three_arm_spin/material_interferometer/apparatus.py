"""Charge-path splitter, terminal windings and a local phase from blank contacts."""

from functools import cache

import numpy as np
from scipy.linalg import expm

from ...live_trimer.measure import LAW, blocks
from ..apparatus import fixture

# Durations were synthesized solely on a local three-port blank chain to
# return a vacancy to its first port with phase exp(2 pi i/3). All entries
# are reduced modulo the exact 2 pi native contact period before use.
PHASE_PULSES = (
    1.3990479646043,
    1.287065241844444,
    3.147103136350763,
    1.3487577241810094,
    0.02368613343709569,
    5.856326422907152,
    4.243939393275632,
    0.9472865671404831,
)
PERIOD = 2 * np.pi


def inverse(gates):
    return [(a, b, (-dt) % PERIOD) for a, b, dt in reversed(gates)]


@cache
def seam(marked, reverse, duration):
    term = blocks()[(1, 1) if marked else "blank"]
    n = 2 if marked else 1
    words = []
    for n_left, n_right in ((0, 1), (1, 0)):
        pair = (n_right, n_left) if reverse else (n_left, n_right)
        words.extend(
            9 * (bit ^ int(reverse and marked)) + LAW.SRC.index(pair)
            for bit in range(n)
        )
    rest = np.setdiff1d(np.arange(len(term)), words)
    assert not term[np.ix_(rest, words)].any()
    filled = np.arange(LAW.SRC.index((1, 1)), len(term), 9)
    h0 = term[filled[0], filled[0]]
    assert np.array_equal(
        term[:, filled], (h0 * np.eye(len(term), dtype=int))[:, filled]
    )
    h = (term[np.ix_(words, words)] - h0 * np.eye(len(words))) / 8
    assert np.max(abs(np.linalg.eigvalsh(h) - np.rint(np.linalg.eigvalsh(h)))) < 1e-12
    return expm(-1j * duration * h)


class Terminal:
    def __init__(self, side):
        f = fixture()
        self.hub = f[side]
        self.arms = [leg[:3] if side == "hub" else leg[:-4:-1] for leg in f["legs"]]
        self.ports = [self.hub] + sum(self.arms, [])
        self.blank = set(f["blank"]) | {(b, a) for a, b in f["blank"]}
        self.marks = {(a, g): bit for bit, (a, c, g) in enumerate(self.arms)}
        assert set(self.marks) <= set(f["marked"])

    def gate(self, state, a, b, dt):
        i, j = self.ports.index(a), self.ports.index(b)
        if (a, b) in self.blank:
            state[[i, j]] = seam(False, a > b, dt) @ state[[i, j]]
            return
        # Listed marked direction a->g fixes the stored coordinate convention.
        bit = self.marks[a, b]
        masks = np.arange(8)
        off = masks[(masks & (1 << bit)) == 0]
        on = off | (1 << bit)
        addresses = np.array([i, i, j, j])[:, None], np.array([off, on, off, on])
        state[addresses] = seam(True, a > b, dt) @ state[addresses]

    def run(self, gates, mask=0):
        state = np.zeros((10, 8), complex)
        state[0, mask] = 1
        for a, b, dt in gates:
            self.gate(state, a, b, dt)
        return state

    def recipe(self, sense):
        split = [
            (self.hub, a, dt)
            for (a, c, g), dt in zip(
                self.arms,
                (2 * np.arcsin(1 / np.sqrt(3)), np.pi / 2, np.pi),
                strict=True,
            )
        ]
        gates = split.copy()
        for index, (a, c, g) in enumerate(self.arms):
            blank_pair = ((a, c), (c, g))
            phase = [(*blank_pair[j % 2], t) for j, t in enumerate(PHASE_PULSES)]
            # A single charge at this corner accumulates the cube-root phase.
            corner = np.outer(np.eye(10)[self.ports.index(a)], np.eye(8)[0])
            calibrated = corner.astype(complex)
            for u, v, dwell in phase:
                self.gate(calibrated, u, v, dwell)
            np.testing.assert_allclose(
                calibrated, corner * np.exp(2j * np.pi / 3), atol=6e-14
            )
            if index:
                gates += (
                    phase if sense * (1 if index == 1 else -1) == 1 else inverse(phase)
                )
            gates += [(a, g, np.pi / 8), (g, c, np.pi), (c, a, np.pi)]
        return gates + inverse(split)

    def zero_flags(self):
        """Nondestructive all-unwound check by three directed crossing heralds."""
        for bit, (a, c, g) in enumerate(self.arms):
            go = [(self.hub, a, np.pi), (a, g, np.pi / 8)]
            for mask in range(8):
                y = self.run(go, mask)
                hit = self.ports.index(g)
                passed = (mask & (1 << bit)) == 0
                np.testing.assert_allclose(
                    np.sum(abs(y[hit]) ** 2), int(passed), atol=1e-12
                )
                if passed:
                    y[np.arange(10) != hit] = 0
                    for u, v, t in inverse(go):
                        self.gate(y, u, v, t)
                    np.testing.assert_allclose(y[0, mask], 1, atol=1e-12)
        return True

    def ports_and_maps(self):
        assert self.zero_flags()
        maps = []
        for sense in (1, -1):
            before = self.run(self.recipe(sense))
            # The inverse pulse sequence returns the full state to initial
            # stock. Its local success bra is the adjoint source column.
            back = before.copy()
            for a, b, t in inverse(self.recipe(sense)):
                self.gate(back, a, b, t)
            np.testing.assert_allclose(
                back,
                np.outer(np.eye(10)[0], np.eye(8)[0]),
                atol=5e-13,
            )
            maps.append(before[0].copy())
        return np.array(maps)
