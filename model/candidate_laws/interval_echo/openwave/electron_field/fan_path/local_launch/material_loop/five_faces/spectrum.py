"""Rotation-resolved stationary lines, bond counters and flux calibration."""

import numpy as np
import sympy as sp
from scipy.linalg import expm

from .wheel import FACES, RIM, rotation, wheel_matrices


def sectors(basis, h):
    perm = rotation(basis)
    np.testing.assert_array_equal(h[np.ix_(perm, perm)], h)
    seed = np.array([i for i, (x, _) in enumerate(basis) if x == RIM[0]])
    modes = {}
    for k in (0, 1, -1, 2, -2):
        f = np.zeros((30, 6), complex)
        idc = seed
        for j in range(5):
            f[idc, np.arange(6)] = np.exp(-2j * np.pi * k * j / 5) / np.sqrt(5)
            idc = perm[idc]
        np.testing.assert_array_equal(idc, seed)
        hk = f.conj().T @ h @ f
        e, v = np.linalg.eigh(hk)
        full = f @ v
        assert np.linalg.norm(h @ full - full * e, ord=2) < 1e-12
        modes[k] = e, full
    if not np.any(h.imag):
        for k in (1, 2):
            np.testing.assert_allclose(modes[k][0], modes[-k][0], atol=5e-14, rtol=0)
    return modes


def certified_lines(h):
    scaled = np.rint(4 * h.real).astype(int)
    np.testing.assert_array_equal(scaled / 4, h)
    # The spectrum and its multiplicities use the full microscopic rational matrix.
    poly = sp.Poly(sp.Matrix(scaled).charpoly().as_expr())
    return [
        (float(a / 4) - 1e-13, float(b / 4) + 1e-13, int(m))
        for (a, b), m in poly.intervals(eps=sp.Rational(1, 10**12))
    ]


def lines(basis, g):
    h, j = wheel_matrices(basis, g)
    modes = sectors(basis, h)
    rows = []
    for lo, hi, mult in certified_lines(h):
        tags = [
            (k, a)
            for k, (en, _) in modes.items()
            for a, e in enumerate(en)
            if lo <= e <= hi
        ]
        assert len(tags) == mult
        v = np.column_stack([modes[k][1][:, a] for k, a in tags])
        p = v @ v.conj().T
        assert abs(p.imag).max() < 3e-13
        rows.append(
            (sum(modes[k][0][a] for k, a in tags) / mult, p.real, tags, [lo, hi])
        )
    assert min(b[0] - a[0] for a, b in zip(rows, rows[1:])) > 0.009
    return h, j, modes, rows


def readers(basis, j):
    # Each enabled two-word block of this side is separately distinguishable
    # from a disabled contact. The two signs within a block remain unresolved.
    out = []
    for s, ay in enumerate(j):
        inside = np.abs(ay).sum(axis=1) > 0
        u, v = FACES[s][1:]
        src = np.diag(
            [
                int(active and x == u)
                for active, (x, _) in zip(inside, basis, strict=True)
            ]
        )
        dest = np.diag(
            [
                int(active and x == v)
                for active, (x, _) in zip(inside, basis, strict=True)
            ]
        )
        turn = expm(1j * np.pi / 4 * np.abs(ay) / 4)
        a, b = (turn.conj().T @ p @ turn for p in (src, dest))
        np.testing.assert_allclose(4 * (b - a), ay, atol=6e-15, rtol=0)
        out.append((a, b))
    return out


def spectral_filter(h, rows, at):
    e, p, _, _ = rows[at]
    # All Kraus factors commute. Each undesired line is rejected by its own wait.
    en, v = np.linalg.eigh(h)
    weights = np.ones(len(en), complex)
    for j, (a, *_) in enumerate(rows):
        if j != at:
            weights *= 0.5 * (1 + np.exp(-1j * np.pi * (en - e) / abs(a - e)))
    np.testing.assert_allclose((v * weights) @ v.conj().T, p, atol=5e-12, rtol=0)
    return p


def flux_test(basis, g, modes, j):
    at_zero = modes[1][1][:, 0]
    primary = float(np.vdot(at_zero, j[0] @ at_zero).real)
    for side in j:
        np.testing.assert_allclose(
            np.vdot(at_zero, side @ at_zero), primary, atol=2e-14
        )
    rows = []
    for phi in (1 / 32, 1 / 64):
        plus = sectors(basis, wheel_matrices(basis, g, phi)[0])
        minus = sectors(basis, wheel_matrices(basis, g, -phi)[0])
        for k in plus:
            np.testing.assert_allclose(plus[k][0], minus[-k][0], atol=6e-14, rtol=0)
        eplus, eminus = plus[1][0][0], plus[-1][0][0]
        derivative = (eplus - eminus) / (2 * phi)
        if g:
            assert abs(derivative / (-primary) - 1) < 0.01
        else:
            assert abs(derivative) < 2e-12 and abs(primary) < 2e-14
        rows.append(
            {
                "positive_phase": phi,
                "energies_plus_one_minus_one": [eplus, eminus],
                "signed_split_over_twice_flux": derivative,
            }
        )
    return {"single_character_directed_current": primary, "phase_probes": rows}


def closed_form(basis, h, g, flux):
    """Independent six-site Jacobi reduction for every angular character."""
    modes = sectors(basis, h)
    for k, (en, _) in modes.items():
        b = abs(g + 4 * np.exp(1j * (2 * np.pi * k + flux) / 5))
        diag = np.array([1, 2, 2, 2, 2, 1]) * (4 + g) + 8
        jac = np.diag(diag) + np.diag(np.full(5, -b), 1) + np.diag(np.full(5, -b), -1)
        np.testing.assert_allclose(np.linalg.eigvalsh(jac), en, atol=6e-14, rtol=0)
