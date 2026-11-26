"""Physical alternating seam contacts, coherently prepared clock and reflected sharp signal."""

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.row_30.experiment import select_route
from model.executable_model import federation

T = 24
HOLDS = (np.pi, np.pi - 0.4)


def gates(exchange, contrast):
    law = candidate.native_law
    a, b = law.local_term_parts()
    h = exchange * a + contrast * b
    basis = np.eye(99)[:, [law.SRC.index((2, 0)), law.SRC.index((0, 2))]]
    np.testing.assert_allclose(h[:, law.SRC.index((0, 0))], 0.0, atol=4e-13)
    reduced = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ reduced, atol=4e-13)
    return [expm(-1j * t * reduced) for t in HOLDS]


def step(psi, layers):
    """Last two indices are physical cell and A/B rail (arbitrary batch axes)."""
    for layer, u in enumerate(layers):
        if layer == 0:
            psi = psi @ u.T
        else:
            left = psi[..., 1]
            right = np.roll(psi[..., 0], -1, axis=-1)
            psi = np.stack(
                (
                    np.roll(u[1, 0] * left + u[1, 1] * right, 1, axis=-1),
                    u[0, 0] * left + u[0, 1] * right,
                ),
                axis=-1,
            )
    return psi


def spin(layers, k, sign, rail, common):
    u, v = layers
    cell = (
        np.array(
            [[v[1, 1], v[1, 0] * np.exp(-1j * k)], [v[0, 1] * np.exp(1j * k), v[0, 0]]]
        )
        @ u
    )
    eig, mat = np.linalg.eig(cell / common)
    index = np.flatnonzero(sign * np.angle(eig) < 0)
    if len(index) != 1:
        raise ArithmeticError("Energy branches unresolved")
    vec = mat[:, index[0]]
    return vec * np.exp(-1j * np.angle(vec[rail]))


def arms(n, kid, layers, exchange, contrast, plane=False):
    common = -np.exp(-1j * (4 * contrast + exchange / 2) * sum(HOLDS))
    rail = int(kid < 0)
    out = []
    offsets = [0] if plane else np.arange(-2, 3)
    lobes = np.exp(-(np.asarray(offsets) ** 2) / (4 * 0.65**2))
    for sign in (1, -1):
        psi = np.zeros((n, 2), complex)
        for m, weight in zip(offsets, lobes, strict=True):
            k = sign * 2 * np.pi * (kid + m) / n
            psi += (
                np.exp(1j * k * np.arange(n))[:, None]
                * spin(layers, k, sign, rail, common)
                * weight
            )
        out.append(psi / np.linalg.norm(psi))
    z = np.array(out)
    np.testing.assert_allclose(np.vdot(*z), 0, atol=5e-12)
    return z


def exact(n, kid, exchange):
    a, b = np.asarray(HOLDS) * exchange / 2
    k = 2 * np.pi * kid / n
    m = np.arccos(np.sin(a) * np.sin(b) - np.cos(a) * np.cos(b))
    e = np.arccos(np.sin(a) * np.sin(b) * np.cos(k) - np.cos(a) * np.cos(b))
    v = np.sin(a) * np.sin(b) * np.sin(k) / np.sin(e)
    return dict(
        k=k,
        rest=m,
        band_energy=e,
        group_velocity=v,
        ideal_comoving_fraction=(e - k * v) / m,
        ideal_lorentz=np.sqrt(1 - v * v),
    )


def draw(p, shots, rng):
    a = np.ravel(p).real
    if np.min(a) < -1e-11 or abs(a.sum() - 1) > 3e-9:
        raise ValueError("Material detector outcomes incomplete")
    a = a.clip(0)
    return rng.multinomial(shots, a / a.sum())


def address(counts):
    n = len(counts)
    z = counts @ np.exp(2j * np.pi * np.arange(n) / n)
    return int(np.rint(n * np.angle(z) / (2 * np.pi)))


def sharp(n, cell, rail):
    a = np.zeros((n, 2), complex)
    a[cell % n, rail] = 1.0
    return a


def go(a, n, layers):
    for _ in range(n):
        a = step(a, layers)
    return a


def radar(n, x0, x1, layers, shots, rng):
    dist = abs(x1 - x0)
    outgoing = int(x1 < x0)
    if 2 * dist >= n:
        raise ValueError("Ambiguous simultaneous return round the opposite side")
    start = sharp(n, x0, outgoing)
    arrival = go(start, dist, layers)
    p = float(abs(arrival[x1 % n, outgoing]) ** 2)
    initial = sharp(n, x1, 1 - outgoing)
    q = float(abs(go(initial, dist, layers)[x0 % n, 1 - outgoing]) ** 2)
    parked = float(
        abs(go(sharp(n, x1, outgoing), dist, layers)[x0 % n, 1 - outgoing]) ** 2
    )
    early = (
        float(np.linalg.norm(go(start, dist - 1, layers)[x1 % n]) ** 2) if dist else 0.0
    )

    def joint(return_chance):
        return draw(
            [1 - p, p * (1 - return_chance), p * return_chance], shots, rng
        ).tolist()

    return dict(
        emission=T - dist,
        turn=T,
        receipt=T + dist,
        sites=[x0, x1],
        forward_then_return=joint(q),
        parked_turn=joint(parked),
        early_clicks=int(rng.binomial(shots, min(1.0, early))),
        empty_clicks=0,
    )


def reference_hand(n, appointments):
    # A disjoint copy of the route, with two complete addressed native swaps per appointment.
    full = gates(1.0, 8.0)[0]
    state = sharp(n, 0, 0)
    seen = {}
    for t in range(max(appointments) + 1):
        if t in appointments:
            k = int(np.argmax(abs(state[:, 0])))
            np.testing.assert_allclose(abs(state[k, 0]), 1, atol=1e-11)
            seen[t] = k
        state = step(state, [full, full])
    return [seen[t] for t in appointments]


def mount(level):
    board = federation(level)
    route, _ = select_route(board)
    assert all(
        b in board.neighbours[a]
        for a, b in zip(route, route[1:] + route[:1], strict=True)
    )
    return list(map(int, route))
