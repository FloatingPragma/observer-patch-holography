from itertools import combinations

import numpy as np
from scipy.linalg import expm
from scipy.special import jv

from ..measure import ring_modes
from .measure import bessel_abs_upper, run


def test_native_return_before_ring_transit():
    data = run()
    rows = {r["time"]: r for r in data["measurements"]}
    assert rows[0]["phase_flip_change"] > 0.63
    assert rows[96]["maximum_change_over_one_station_channels"] < 0.0073
    for r in rows.values():
        np.testing.assert_allclose(abs(r["phase_flip_change"]), r["maximum_change_over_one_station_channels"], atol=1e-13)
    for r in data["certified_windows"]:
        t = r["start"]
        assert float(bessel_abs_upper(0, t)) >= abs(jv(0, t)) - 1e-15
        assert 2 * np.hypot(jv(0, t), jv(1, t)) <= r["universal_bond_error_upper"]
        assert r["propagator_winding_error_upper"] < 6e-69


def test_phase_flip_in_the_actual_fixed_load_spin_basis():
    n, m = 6, 3
    words = [sum(1 << j for j in word) for word in combinations(range(n), m)]
    at = {w: i for i, w in enumerate(words)}
    edges = []
    for j in range(n):
        h = np.zeros((len(words), len(words)))
        mask = (1 << j) | (1 << ((j + 1) % n))
        for col, w in enumerate(words):
            if (w & mask).bit_count() == 1:
                h[at[w ^ mask], col] = -0.5
        edges.append(h)
    total = sum(edges)
    values, vectors = np.linalg.eigh(total)
    assert values[1] > values[0] + 0.1
    vacuum = vectors[:, 0]
    flipped = vacuum * [(-1)**(w & 1) for w in words]
    k, basis, full = ring_modes(n)
    for t in (0, 2.37, 5):
        direct = expm(-1j * t * total) @ flipped
        change = np.vdot(direct, edges[0] @ direct).real - np.vdot(vacuum, edges[0] @ vacuum).real
        phase = np.exp(1j * t * np.cos(k))
        a = basis[:2] @ (phase * basis[0].conj())
        d = basis[:2] @ (phase * (full - 0.5) * basis[0].conj())
        z = -0.5 * (d[0].conjugate() * a[1] + d[1].conjugate() * a[0])
        np.testing.assert_allclose(change, -4 * z.real, atol=1e-13)
