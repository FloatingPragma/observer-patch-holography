"""Range and independently metered cubical shells on a three-address native register."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import product
import json
from math import lcm
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from sympy import ZZ
from sympy.polys.matrices import DomainMatrix

from model.candidate_laws.interval_echo import candidate
from model.executable_model import CELL_SEAMS, federation
from model.experiments.interval_echo_consensus.refinement_cooling.dimer_band.measure import (
    pointer,
)

HERE = Path(__file__).resolve().parent
CENTRES = (5, 6)
READY = (-1, -1, -1)
DEPTH, SHOTS = 5, 65536


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


@cache
def layout():
    board = federation(0)
    links = {frozenset(edge) for edge in CELL_SEAMS}

    def extend(path):
        if len(path) == 12:
            return path
        for v in range(12):
            if v not in path and (not path or frozenset((v, path[-1])) in links):
                found = extend((*path, v))
                if found:
                    return found
        return ()

    local = extend(())
    axes = tuple(tuple(12 * o + v for v in local) for o in range(3))
    flag = (36, 37)
    for chain in (*axes, flag):
        require(
            all(v in board.neighbours[u] for u, v in zip(chain, chain[1:])),
            "missing rail seam",
        )
    return board, axes, flag


def layer(pos):
    if pos == READY:
        return 0
    return 1 + max(min(abs(x - c) for c in CENTRES) for x in pos)


def neighbors(pos):
    for axis in range(3):
        for delta in (-1, 1):
            dest = list(pos)
            dest[axis] += delta
            if 0 <= dest[axis] < 12:
                yield tuple(dest)


def edges():
    for corner in product(CENTRES, repeat=3):
        yield READY, corner
    for p in product(range(12), repeat=3):
        for q in neighbors(p):
            if p < q:
                yield p, q


@cache
def read_contact(parked, side):
    """Locate the sharp double after one actual blank two-port pulse."""
    hop, contrast = candidate.native_law.local_term_parts()
    h = (0 if parked else hop) + candidate.KAPPA * contrast
    labels = (
        candidate.native_law.SRC.index((2, 0)),
        candidate.native_law.SRC.index((0, 2)),
    )
    col = expm(-1j * np.pi * h)[:, labels[side]]
    probs = abs(col) ** 2
    chosen = int(np.argmax(probs))
    np.testing.assert_allclose(probs, np.eye(len(h))[chosen], atol=2e-13, rtol=0)
    return labels.index(chosen)


def pulse(pair, before, parked=False):
    side = pair.index(before)
    return pair[read_contact(parked, side)]


def release_corner(corner):
    _, axes, flag = layout()
    require(pulse(flag, flag[0]) == flag[1], "source flag did not fire")
    end = []
    for axis, target in zip(axes, corner, strict=True):
        p = axis[5]
        if target == 6:
            p = pulse((axis[5], axis[6]), p)
        end.append(axis.index(p))
    return tuple(end)


def calibrate():
    hop, contrast = candidate.native_law.local_term_parts()
    src = candidate.native_law.SRC
    targets = [src.index((b, a)) for a, b in src]
    for name, h in (
        ("native", hop + candidate.KAPPA * contrast),
        ("exchange_off", candidate.KAPPA * contrast),
    ):
        columns = expm(-1j * np.pi * h)[:, : len(src)]
        expected = targets if name == "native" else list(range(len(src)))
        np.testing.assert_allclose(
            columns, np.eye(len(h))[:, expected], atol=2e-13, rtol=0
        )
    pointer()
    board, axes, _ = layout()
    # Each upper bit may be read against an empty neighbouring port by routing
    # a fresh pointer. Direct sharp occupation is the supplied address read.
    return [
        dict(
            ports=list(axis),
            directed_arrivals=[
                [32 * int(pulse((u, v), u) == v), 32 * int(pulse((u, v), v) == u)]
                for u, v in zip(axis, axis[1:])
            ],
            parked_arrivals=[
                [
                    32 * int(pulse((u, v), u, True) == v),
                    32 * int(pulse((u, v), v, True) == u),
                ]
                for u, v in zip(axis, axis[1:])
            ],
        )
        for axis in axes
    ]


def timing():
    _, axes, _ = layout()
    surveys = []
    positions = list(product(range(12), repeat=3))
    for corner in product(CENTRES, repeat=3):
        rows = []
        for target in positions:
            lamps = [axes[i][p] for i, p in enumerate(release_corner(corner))]
            expected = [axes[i][p] for i, p in enumerate(target)]
            ticks = 1
            while lamps != expected:
                for i, path in enumerate(axes):
                    old = path.index(lamps[i])
                    if lamps[i] != expected[i]:
                        # One surveyed full contact per foreign-disjoint coordinate.
                        direction = 1 if target[i] > old else -1
                        lamps[i] = pulse((lamps[i], path[old + direction]), lamps[i])
                ticks += 1
            rows.append(ticks)
        surveys.append(rows)
    measured = dict(zip(positions, map(min, zip(*surveys, strict=True)), strict=True))
    measured[READY] = 0
    return measured, surveys


def measured_cuts(ranges, cut):
    _, axes, _ = layout()
    census = Counter()
    arrivals = Counter()
    for a, b in edges():
        if a == READY:
            # A source tag selects a corner and the independent fired flag.
            # Keep its eight transactions distinct even at central (5,5,5).
            dst = release_corner(b)
        else:
            axis = next(i for i in range(3) if a[i] != b[i])
            crossed = (axes[axis][a[axis]], axes[axis][b[axis]])
            observed = pulse(crossed, crossed[0], parked=(cut == (axis, a[axis])))
            reverse = pulse(crossed, crossed[1], parked=(cut == (axis, a[axis])))
            require(
                (observed == crossed[1]) == (reverse == crossed[0]),
                "asymmetric measured bond",
            )
            word = list(a)
            word[axis] = axes[axis].index(observed)
            dst = tuple(word)
        arrivals[(a == READY, dst == b)] += 1
        small, large = sorted([ranges[a], ranges[dst]])
        for n in range(small, large):
            census[n] += 1
    return dict(
        outward=[census[n] for n in range(DEPTH + 1)],
        source_arrivals=arrivals[True, True],
        active_arrivals=arrivals[False, True],
        missed_active=arrivals[False, False],
    )


def orbit(p):
    return tuple(sorted(min(abs(x - c) for c in CENTRES) for x in p))


def representative(t):
    return tuple(6 + i for i in t)


def reward(p, q):
    """Face-centre crossings minus face-corner crossings at shell two."""
    lp, lq = layer(p), layer(q)
    if {lp, lq} != {2, 3}:
        return 0
    inside, outside = (p, q) if lp < lq else (q, p)
    tangents = tuple(
        sorted(
            min(abs(a - c) for c in CENTRES)
            for a, b in zip(inside, outside, strict=True)
            if a == b
        )
    )
    value = int(tangents == (0, 0)) - int(tangents == (1, 1))
    return value if lp < lq else -value


def solve_integer(A, rhs):
    dens = 1
    for value in rhs:
        dens = lcm(dens, value.denominator)
    b = DomainMatrix([[ZZ(int(v * dens))] for v in rhs], (len(rhs), 1), ZZ)
    num, den = A.solve_den(b)
    return [Fraction(int(v), int(den) * dens) for v in num.to_Matrix()]


def steady():
    positions = list(product(range(12), repeat=3))
    types = sorted({orbit(p) for p in positions if layer(p) <= DEPTH})
    index = {t: i for i, t in enumerate(types)}
    A = [[0] * len(types) for _ in types]
    for i, t in enumerate(types):
        A[i][i] = 6
        for q in neighbors(representative(t)):
            if layer(q) <= DEPTH:
                A[i][index[orbit(q)]] -= 1
    matrix = DomainMatrix(
        [[ZZ(x) for x in row] for row in A], (len(types), len(types)), ZZ
    )
    phi = solve_integer(
        matrix, [Fraction(1, 8) if t == (0, 0, 0) else Fraction(0) for t in types]
    )

    def at(values, p):
        return values[index[orbit(p)]] if layer(p) <= DEPTH else Fraction(0)

    # Verify against every unsymmetrized vertex and bond on the actual graph.
    outward = Counter()
    class_currents = {}
    net = Fraction(0)
    for p in positions:
        if layer(p) > DEPTH:
            continue
        balance = Fraction(0)
        for q in neighbors(p):
            flux = at(phi, p) - at(phi, q)
            balance += flux
            if layer(q) > layer(p):
                outward[layer(p)] += flux
                if reward(p, q):
                    tag = str(reward(p, q))
                    if tag in class_currents:
                        require(class_currents[tag] == flux, "a cube orbit split")
                    class_currents[tag] = flux
                    net += reward(p, q) * flux
        if layer(p) == 1:
            balance -= Fraction(1, 8)  # incoming source-flag transaction
        require(balance == 0, "noncentral sourced vertex")
    require(all(outward[n] == 1 for n in range(1, DEPTH + 1)), "flux discontinuity")
    # First and raw second moment of the signed per-launch contrast. Idle
    # terms cancel from the Dirichlet equations (coin tags six and seven).
    means = solve_integer(
        matrix,
        [
            Fraction(
                sum(reward(representative(t), q) for q in neighbors(representative(t)))
            )
            for t in types
        ],
    )
    second = solve_integer(
        matrix,
        [
            sum(
                Fraction(reward(representative(t), q) ** 2)
                + 2 * reward(representative(t), q) * at(means, q)
                for q in neighbors(representative(t))
            )
            for t in types
        ],
    )
    require(means[0] == net, "independent expectation calculation disagrees")
    return dict(
        grounded_orbits=len(types),
        full_interior=sum(layer(p) <= DEPTH for p in positions),
        source_vertex_outflow="1",
        outward=[str(outward[n]) for n in range(1, DEPTH + 1)],
        per_edge={k: str(j) for k, j in class_currents.items()},
        mean=str(means[0]),
        variance=str(second[0] - means[0] ** 2),
    )


def sample():
    rng = np.random.default_rng(771903)
    released = rng.integers(0, 8, SHOTS)
    pos = 5 + ((released[:, None] >> np.arange(3)) & 1)
    active = np.arange(SHOTS)
    contrast = np.zeros(SHOTS, dtype=int)
    completed = Counter()
    # Matrix-free signed reward lookup on 1,728 vertices.
    dims = (12, 12, 12)
    positions = list(product(range(12), repeat=3))
    events = np.zeros((12**3, 8), dtype=int)
    destinations = np.repeat(np.arange(12**3)[:, None], 8, axis=1)
    for i, p in enumerate(positions):
        if layer(p) == 6:
            continue
        for k in range(6):
            q = list(p)
            q[k // 2] += 1 if k % 2 else -1
            events[i, k] = reward(p, tuple(q))
            destinations[i, k] = np.ravel_multi_index(q, dims)
    weights = np.array([layer(p) for p in positions])
    words = np.ravel_multi_index(pos.T, dims)
    ticks = 0
    while active.size:
        ticks += 1
        tags = rng.integers(8, size=len(active))
        old = words[active]
        contrast[active] += events[old, tags]
        words[active] = destinations[old, tags]
        done = weights[words[active]] == 6
        if done.any():
            completed[ticks] = int(done.sum())
        active = active[~done]
    return dict(
        released=np.bincount(released, minlength=8).tolist(),
        completions=sorted(completed.items()),
        contrast=sorted(Counter(contrast.tolist()).items()),
    )


def report():
    banks = calibrate()
    ranges, survey = timing()
    require(all(ranges[p] == layer(p) for p in ranges), "counted range disagreement")
    return dict(
        axes=banks,
        flag=layout()[2],
        root=READY,
        native_range_reads=survey,
        surveyed=[measured_cuts(ranges, c) for c in (None, (0, 5), (0, 6))],
        steady=steady(),
        walk=sample(),
        shots=SHOTS,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = json.dumps(report(), separators=(",", ":")) + "\n"
    out = HERE / "records.json"
    if args.write:
        out.write_text(data)
    elif data != out.read_text():
        raise SystemExit("Product-address survey differs")
