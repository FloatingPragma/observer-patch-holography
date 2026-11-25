"""Free signed two-constituent evolution with the original writable gateway heads."""

from functools import cache
import numpy as np
from scipy.linalg import eigh
from .....live_trimer.measure import LAW, blocks
from .geometry import layout

TIMES = (0, 32, 64, 128, 256, 512, 1024, 2048, 4096)


@cache
def edges(written=True):
    ring, marks, leaves, _, _ = layout()
    return [(*e, None) for e in zip(ring, ring[1:] + ring[:1])] + [
        (a, b, j) for j, (a, b) in enumerate(marks) if written
    ]


def transitions(state, written):
    vacancy, excess, mask = state
    move = {}
    diag = 0.0
    for a, b, j in edges(written):
        loads = tuple(
            0 if v == vacancy else 2 if v == excess else 1 for v in sorted((a, b))
        )
        if loads == (1, 1):
            continue
        mat = blocks()["blank" if j is None else (1, 1)]
        coin = 0 if j is None else ((mask >> j) & 1) ^ int(a > b)
        col = 9 * coin + LAW.SRC.index(loads)
        host = 9 * coin + LAW.SRC.index((1, 1))
        diag += (mat[col, col] - mat[host, host]) / 8
        for row in np.flatnonzero(mat[:, col]):
            if row == col:
                continue
            assert LAW.SRC[row % 9] == loads[::-1]
            swap = {a: b, b: a}
            newmask = mask if j is None else mask ^ (1 << j)
            assert j is None or row // 9 == 1 - coin
            adj = (swap.get(vacancy, vacancy), swap.get(excess, excess), newmask)
            move[adj] = move.get(adj, 0.0) + mat[row, col] / 8
    move[state] = diag
    return move


@cache
def chart(start, written=True):
    states = [start]
    seen = {start}
    for x in states:
        for y in transitions(x, written):
            if y not in seen:
                seen.add(y)
                states.append(y)
    index = {x: j for j, x in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for j, x in enumerate(states):
        for y, v in transitions(x, written).items():
            h[index[y], j] = v
    np.testing.assert_array_equal(h, h.T)
    return tuple(states), h


def signed_pair(sign, init):
    return init if sign == -1 else init[::-1]


def frame(sign, mask, init=(160, 167), written=True):
    ring, marks, leaves, _, _ = layout()
    n = len(ring)
    vacancy, excess = signed_pair(sign, init)
    states, h = chart((vacancy, excess, mask), written)
    near = np.array(
        [
            u in ring
            and v in ring
            and (ring.index(u) - ring.index(v)) % n in (1, n - 1)
            for u, v, _ in states
        ]
    )
    # Each written side is a bridge; a sharp starting side record fixes a
    # conserved difference between its leaf load and stored direction bit.
    assert len({(u, v) for u, v, _ in states}) == len(states)
    vols = np.diag(h)
    r = np.sum(abs(h), axis=1) - abs(vols)
    low = float(max((vols + r)[~near]))
    high = float(min((vols - r)[near]))
    assert low <= 36 and high >= 40
    separate = h[np.ix_(~near, ~near)]
    contact = h[np.ix_(near, near)]
    off = h[np.ix_(~near, near)]
    min_contact = float(
        np.min(np.diag(contact) - np.sum(abs(contact), axis=1) + abs(np.diag(contact)))
    )
    max_sep = float(
        np.max(
            np.diag(separate) + np.sum(abs(separate), axis=1) - abs(np.diag(separate))
        )
    )
    cross_norm = np.sqrt(
        np.max(np.sum(abs(off), axis=0)) * np.max(np.sum(abs(off), axis=1))
    )
    residual = np.linalg.norm(h[0] - np.eye(len(h))[0] * vols[0])
    certified = (
        max(
            0.0,
            1
            - (residual / (vols[0] - low) + cross_norm / (min_contact - max_sep)) ** 2,
        )
        if near[0]
        else 0.0
    )
    ev, U = eigh(h, driver="evd")
    upper = ev > 38
    assert sum(upper) == sum(near) == 2 * n
    upper_weight = float(np.sum(U[0, upper] ** 2))
    outside = float(np.linalg.norm(U[~near][:, upper], 2))
    lower = max(
        0.0, 1 - (np.sqrt(1 - upper_weight) + outside * np.sqrt(upper_weight)) ** 2
    )
    prop = (np.exp(-1j * np.outer(TIMES, ev)) * U[0]) @ U.T
    p = abs(prop) ** 2
    np.testing.assert_allclose(p.sum(axis=1), 1, atol=1e-12)
    state = np.array(states)
    p_near = p[:, near].sum(axis=1)
    sensors = {"opposite": (207, 206), "second_owner": (106, 96)}
    field = {}
    for name, (a, b) in sensors.items():
        probs = []
        for t, w in zip(TIMES, p):
            # Independently measuring the two physical ternary loads gives nine
            # possible outcomes. On these sectors (0,0),(2,2) cannot occur.
            va, ve = state[:, :2].T
            aa = np.where(va == a, 0, np.where(ve == a, 2, 1))
            bb = np.where(va == b, 0, np.where(ve == b, 2, 1))
            probs.append(np.bincount(3 * aa + bb, weights=w, minlength=9))
        field[name] = np.asarray(probs)
    return (
        dict(
            component_size=len(states),
            diagonal_exclusion=[low, high],
            upper_dimension=int(sum(upper)),
            spectral_gap=[float(max(ev[~upper])), float(min(ev[upper]))],
            high_energy_weight=upper_weight,
            separation_projection_norm=outside,
            adjacent_all_time_lower=lower,
            rational_bounds=[
                low,
                min_contact,
                max_sep,
                float(cross_norm**2),
                float(residual**2),
                float(vols[0]),
            ],
            analytic_all_time_floor=certified,
        ),
        p_near,
        field,
    )
