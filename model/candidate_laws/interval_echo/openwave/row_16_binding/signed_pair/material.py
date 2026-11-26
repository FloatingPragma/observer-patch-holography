"""Native signed word sectors, their contact spectra and threefold preparation."""

from itertools import combinations

import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation

from ...row_14_profile import measure as profile
from model.executable_model import carrier


def vertices():
    v, edges = carrier.icosahedron_from_coordinates()
    return v[carrier.isomorphism(carrier.cell().seams, edges)] / np.linalg.norm(v[0])


def axes(direction):
    z = direction / np.linalg.norm(direction)
    x = np.eye(3)[np.argmin(abs(z))]
    x = x - z * (z @ x)
    x /= np.linalg.norm(x)
    return np.array([z, x, np.cross(z, x)])


def charge_words(plus, minus):
    words = []
    for p in combinations(range(12), plus):
        remaining = sorted(set(range(12)) - set(p))
        for m in combinations(remaining, minus):
            w = np.zeros(12, dtype=np.int8)
            w[list(p)] = 1
            w[list(m)] = -1
            words.append(w)
    return np.asarray(words)


def quantities(w):
    mask, _ = profile.native_geometry()
    a, b = profile.EDGES.T
    d = w[:, a] - w[:, b]
    potential = (d.astype(int) ** 2).sum(axis=1) + (((w[:, mask] + 1) ** 2) - 1).sum(
        axis=1
    )
    changes = (d != 0).sum(axis=1)
    return potential, changes


def bound(plus, minus):
    w = charge_words(plus, minus)
    potential, changes = quantities(w)
    return int(8 * potential.min()), float((8 * potential + changes / 2).min())


class Assembly:
    def __init__(self, plus, minus, contrast=8):
        self.w, self.swaps = profile.sector(plus, minus)
        potential, changes = quantities(self.w)
        self.v = contrast * potential
        self.diag = self.v + changes / 2
        self.i = int(np.argmin(8 * potential + changes / 2))
        idx = np.arange(len(self.w))
        rows, cols = [], []
        for perm in self.swaps:
            take = idx != perm
            rows.extend(idx[take])
            cols.extend(perm[take])
        self.h = sparse.diags(self.diag).tocsr() + sparse.csr_matrix(
            (np.full(len(rows), -0.5), (rows, cols)), shape=(len(idx), len(idx))
        )
        self.bracket = (int(self.v.min()), int((self.v + changes).max()))
        self.contrast = contrast
        v = vertices()
        mask, _ = profile.native_geometry()
        self.axes = axes(v[~mask].sum(axis=0))
        turn = Rotation.from_rotvec((2 * np.pi / 3) * self.axes[0]).apply(v)
        perm = np.argmax(turn @ v.T, axis=1)
        assert np.allclose(turn, v[perm]) and np.array_equal(mask, mask[perm])
        original = self.w[self.i].copy()
        index = {row.tobytes(): i for i, row in enumerate(self.w)}
        self.orbit = []
        for _ in range(3):
            self.orbit.append(index[original.tobytes()])
            original = original[np.argsort(perm)]
        assert np.array_equal(original, self.w[self.i])

    def sources(self):
        seeds = np.zeros((len(self.w), 3), complex)
        seeds[self.i, 0] = 1
        for k in (0, 1):
            for j, idx in enumerate(self.orbit):
                seeds[idx, k + 1] += np.exp(-2j * np.pi * j * k / 3) / 3
        return seeds

    def contact_bank(self, state, rng, shots):
        """Select one of the 39 native terms; return every possible energy outcome."""
        energies = {}

        def add(energy, weight):
            energies[float(energy)] = energies.get(float(energy), 0.0) + weight

        for (a, b), sw in zip(profile.EDGES, self.swaps):
            mag = abs(self.w[:, a] - self.w[:, b])
            for difference in range(3):
                ids = mag == difference
                for parity in (0, 1):
                    amp = (state[ids] + (-1) ** parity * state[sw[ids]]) / 2
                    add(self.contrast * difference**2 + parity, np.vdot(amp, amp).real)
        p = abs(state) ** 2
        mask, _ = profile.native_geometry()
        for site in np.flatnonzero(mask):
            for q in (-1, 0, 1):
                add(self.contrast * ((q + 1) ** 2 - 1), p[self.w[:, site] == q].sum())
        levels = sorted(energies)
        weights = np.array([energies[y] for y in levels]) / 39
        assert abs(weights.sum() - 1) < 1e-10
        assert abs(39 * (weights @ levels) - np.vdot(state, self.h @ state).real) < 1e-9
        counts = rng.multinomial(shots, weights / weights.sum())
        return {
            "levels": levels,
            "counts": counts.tolist(),
            "energy": float(39 * counts @ np.asarray(levels) / shots),
        }


def decompositions(w):
    """Census connected defect subsets before looking at amplitudes."""
    powers = 1 << np.arange(12)
    positive = (w > 0) @ powers
    support = (w != 0) @ powers
    links = carrier.cell().neighbours
    connected = np.zeros(2**12, bool)
    for mask in range(1, len(connected)):
        todo = [int((mask & -mask).bit_length() - 1)]
        reached = 0
        while todo:
            v = todo.pop()
            if reached & (1 << v):
                continue
            reached |= 1 << v
            todo.extend(
                u for u in links[v] if mask & (1 << u) and not reached & (1 << u)
            )
        connected[mask] = reached == mask
    multiplicities = np.zeros(len(w), int)
    for selected in combinations(range(12), 3):
        bits = int(sum(1 << p for p in selected))
        if connected[bits]:
            multiplicities += (
                ((support & bits) == bits)
                & ((w[:, selected] > 0).sum(axis=1) == 2)
                & connected[support ^ bits]
            )
    a, b = profile.EDGES.T
    prod = w[:, a] * w[:, b]
    return (
        np.array(
            [
                connected[support],
                multiplicities > 0,
                connected[positive],
                connected[support ^ positive],
            ]
        ),
        (prod == 1).sum(axis=1),
        (prod == -1).sum(axis=1),
    )


def competing_orbit():
    """Energy and Gram matrices on simultaneous rotated sharp fragment pairs."""
    sites = 1 << np.arange(12)
    grams, energies, masks = [], [], []
    for stock in ((4, 0), (0, 3)):
        assembly = Assembly(*stock)
        orbit = np.array(assembly.orbit)
        grams.append(np.equal.outer(orbit, orbit))
        energies.append(assembly.h[orbit][:, orbit].toarray())
        masks.append((abs(assembly.w[orbit]) @ sites).tolist())
    gram = grams[0] * grams[1]
    compression = energies[0] * grams[1] + grams[0] * energies[1]
    doubled = np.rint(2 * compression).astype(int)
    assert np.array_equal(2 * compression, doubled)
    return {
        "masks": masks,
        "gram": gram.astype(int).tolist(),
        "doubled_total_energy": doubled.tolist(),
    }
