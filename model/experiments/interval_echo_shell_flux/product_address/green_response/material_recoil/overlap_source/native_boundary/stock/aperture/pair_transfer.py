"""Original-contact work and a held transit after matching two ordinary probe energies."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json

import numpy as np
from scipy.optimize import linprog

from . import pair_work as lab

HERE = lab.HERE
pair = lab.energy.pair


def original_contact(board, home, marks, h, w, cells):
    graph, _ = lab.mount(board, home, marks)
    base = np.zeros((len(w), board.n), int)
    base[:, home] = w[:, :4]
    e = lab.whole(board, marks, w, base, h)

    def place(x, y):
        flat = base.copy()
        flat[:, home[:2]] = 0
        flat[:, [x, y]] = w[:, :2]
        return lab.whole(board, marks, w, flat, h)

    moves = [list(graph[v]) for v in home[:2]]
    for x in moves[0]:
        for y in moves[1]:
            assert len(set([*home[:2], x, y])) == 4
            np.testing.assert_array_equal(
                place(x, y) - place(x, home[1]) - place(home[0], y) + e, 0
            )
    x, y = w[:, :2].T
    assert (
        marks[0] == (10, 20)
        and graph.has_edge(8, 6)
        and graph.has_edge(6, 10)
        and graph.has_edge(20, 13)
    )
    measured = place(10, 20) - place(6, 20) - place(10, 13) + place(6, 13)
    diag = h.diagonal().reshape(11, 3, 3)
    k = w[:, 4]
    np.testing.assert_array_equal(
        measured, diag[k, x, y] - diag[k, x, 0] - diag[k, 0, y] + diag[k, 0, 0]
    )
    np.testing.assert_array_equal(measured, -65 * (cells >= 2))
    return dict(
        one_step_pairs=int(np.prod(list(map(len, moves)))),
        mixed_work_one_step=0,
        left_stroke=[8, 6, 10],
        right_withdrawal=[20, 13],
        marked_mixed_work=str(F(int(measured.sum()), 4 * len(w))),
        changes_at_fixed_head_marginal=0,
    )


def scale(a, k):
    ends = [v * k for v in a]
    return min(ends), max(ends)


def divide(a, b):
    assert b[0] * b[1] > 0
    vals = [x / y for x in a for y in b]
    return min(vals), max(vals)


def stationary(w, cells, codes, h, ids, groups, v):
    f = (w[:, 1] - 1) ** 2
    g = (w[:, 0] - 1) ** 2
    rows = np.vstack((pair.constraints(w, cells, codes, h, "energy"), w[:, :4].T, f, g))
    A = pair.group_matrix(rows, ids, groups)
    b = np.array(A.sum(1), object)
    b[-2:] = [F(5, 28) * len(w), v * len(w)]
    sol = linprog(np.zeros(groups), A_eq=A, b_eq=np.array(b, float), bounds=(0, None))
    assert sol.success
    z = np.array(pair.rational(sol.x), object)
    assert all(z >= 0) and np.all(A @ z == b)
    return z


def run():
    _, w, cells, _, codes, h, ids, groups = pair.setup()
    board, home, marks, *_ = lab.energy.returned(w)
    result = original_contact(board, home, marks, h, w, cells)
    _, path = lab.mount(board, home, marks)
    degree = int(board.degree[path[1]])
    # For any eligible sharp distribution with both means one the pair table
    # is completely determined by f and v. Eliminate the ordinary gain.
    f0, f1, v0 = F(3, 28), F(5, 28), F(11, 28)
    tables = {}
    gain = {}
    coeff = []
    for t in (1, 2):
        vals = [[lab.reader(t, h, degree, i, j) for j in range(3)] for i in range(3)]
        lev = [[v[0] for v in row] for row in vals]
        mid = [(i + j) / 2 for i, j in zip(lev[0][1], lev[2][1], strict=True)]
        gain[t] = lab.sub(lev[1][0], lev[0][0])
        coeff.append(
            [
                divide(lab.sub(mid, lev[1][1]), gain[t]),
                divide(lab.sub(lev[0][2], mid), gain[t]),
            ]
        )
        tables[t] = vals
    dc = [lab.sub(a, b) for a, b in zip(*coeff, strict=True)]
    dv = divide(scale(dc[1], f0 - f1), dc[0])
    root = [v0 + x for x in dv]
    target = F(17, 50)
    endpoints = [F(3399, 10000), F(3401, 10000)]
    assert endpoints[0] < root[0] <= root[1] < endpoints[1]
    weights = {
        v: stationary(w, cells, codes, h, ids, groups, v) for v in [*endpoints, target]
    }

    # By convexity the exact algebraic matching value is attainable. Check
    # that the third observable disagrees on the entire rational bracket.
    def table(f, v):
        joint = np.full((3, 3), F(), object)
        joint[1, 1] = 1 - v
        joint[0, 1] = joint[2, 1] = (v - f) / 2
        joint[0, 2] = joint[2, 0] = f / 2
        return joint

    p = weights[target][ids] / len(w)
    actual = np.array(
        [[sum(p[(w[:, 0] == i) & (w[:, 1] == j)]) for j in range(3)] for i in range(3)]
    )
    assert np.array_equal(table(f1, target), actual)
    rng = np.random.default_rng(692173)
    arms = {}
    for j, t in enumerate((1, 2)):
        cf = coeff[j]

        def mass(v):
            a, b = scale(cf[0], v - v0), scale(cf[1], f1 - f0)
            return (a[0] + b[0], a[1] + b[1])

        massband = [
            float(min(mass(v)[0] for v in root)),
            float(max(mass(v)[1] for v in root)),
        ]
        data = np.array([row[1:] for cells_ in tables[t] for row in cells_])
        gaps, amplitudes, arrive = data.T
        ref = float(np.asarray(table(f0, v0), float).ravel() @ arrive)
        ordinary = arrive[3] - arrive[0]
        assert ordinary < 0
        expected = [
            float(np.asarray(table(f1, v), float).ravel() @ arrive) for v in endpoints
        ]
        bracket = [min(expected) - ref, max(expected) - ref]
        trial = rng.choice(9, lab.SHOTS, p=np.array(actual, float).ravel())
        time = rng.uniform(0, lab.DURATION, lab.SHOTS)
        hits = int(
            sum(
                rng.random(lab.SHOTS)
                < amplitudes[trial] * np.sin(gaps[trial] * time / 2) ** 2
            )
        )
        baseline = json.loads((HERE / "pair_work.json").read_bytes())["mobile"][str(t)][
            "original"
        ]["arrivals"]
        bandwidth = 2 * np.sqrt(np.log(8 / 0.01) / (2 * lab.SHOTS))
        diff = (hits - baseline) / lab.SHOTS
        held = [diff - bandwidth, diff + bandwidth]
        assert massband[0] > 0 and held[0] > 0 and bracket[0] > 0
        arms[str(t)] = dict(
            matching_mass=massband,
            rational_arm_mass=list(map(float, mass(target))),
            held_expected_transit_change=bracket,
            sparse_ordinary_prediction=[massband[1] * ordinary, massband[0] * ordinary],
            rational_arm_arrivals=hits,
            held_measured_change_99=held,
        )
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "pair_work.json", HERE / "closed_pairing.json")
        },
        original=result,
        matching_partner_variance=list(map(float, root)),
        middle_variance=str(f1),
        stationary_group_weights={
            str(v): list(map(str, z)) for v, z in weights.items()
        },
        ordinary_transfer=arms,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    result = run()
    dst = HERE / "pair_transfer.json"
    if args.write:
        dst.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        lab.check(json.loads(dst.read_bytes()), result)
    print(
        result["original"],
        result["matching_partner_variance"],
        result["ordinary_transfer"],
    )
