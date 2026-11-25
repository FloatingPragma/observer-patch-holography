"""Directional, one-way first-hit clock responses on the formed and equilibrium collars."""

from itertools import combinations, permutations, product
from pathlib import Path
import json
import argparse

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

from model.experiments.interval_echo_einstein_station.collar_pairing import (
    measure as pairing,
)
from model.experiments.interval_echo_einstein_station.collar_pairing.ranging import (
    measure as ranging,
)
from ..full_geometry.measure import leg_kernels
from . import matched


HERE = Path(__file__).resolve().parent
PAIRS = list(permutations(range(4), 2))
SEND, RECEIVE = np.eye(4)[np.array(PAIRS).T]


def directions(x):
    total = x.sum(axis=0)
    _, right = pairing.formation.marginals(total)
    product = np.array([np.kron(pairing.formation.marginals(b)[0], right) for b in x])
    ordinary = np.array([total / 2 - b for b in x])
    return dict(state=x, unpair=product - x, ordinary=ordinary)


def source():
    raw = json.loads((HERE.parent / "full_algebra/results.json").read_text())
    x = np.array(raw["x_real"]) + 1j * np.array(raw["x_imag"])
    return directions(x)


def arrivals(u, step, cap):
    ops = []
    for i, j in PAIRS:
        first = ranging.first_hits(u[: cap + 1], i, j)
        effect = cap * step * np.eye(4, dtype=complex)
        for n in range(1, cap):
            effect += (n - cap) * step * (first[n].conj().T @ first[n])
        ops.append(effect)
    return ops


def rectangles(p, m):
    r = []
    for a, b in combinations(range(4), 2):
        j, k = sorted({0, 1, 2, 3} - {a, b})
        ids = [PAIRS.index(t) for t in ((a, j), (a, k), (b, k), (b, j))]

        def differences(vals):
            one, two, three, four = vals[ids]
            return [one - two, three - four]

        v, w = differences(p), differences(m)
        if np.sign(v[0]) == np.sign(v[1]) == -np.sign(w[0]) == -np.sign(w[1]):
            r.append(dict(senders=[a, b], receivers=[j, k], excess=v, vacancy=w))
    return r


def gain_feasibility(plus, minus):
    # Orthogonally eliminate all row/column offsets. Normalize the eight gains
    # to sum to one and maximize their smallest entry.
    annihilator = null_space(np.column_stack([SEND, RECEIVE]).T).T
    a = annihilator @ np.column_stack([SEND * plus[:, None], -SEND * minus[:, None]])
    eq = np.block([[a, np.zeros((len(a), 1))], [np.ones((1, 8)), np.zeros((1, 1))]])
    result = linprog(
        np.r_[np.zeros(8), -1],
        A_eq=eq,
        b_eq=np.r_[np.zeros(len(a)), 1],
        A_ub=np.column_stack([-np.eye(8), np.ones(8)]),
        b_ub=np.zeros(8),
        bounds=(None, None),
        method="highs",
    )
    if not result.success:
        return dict(feasible_positive=False)
    return dict(
        feasible_positive=bool(result.x[-1] > 0),
        best_minimum_gain=float(result.x[-1]),
        normalized_gains=result.x[:-1].tolist(),
    )


def run():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    x = source()
    x["dephased_unpair"] = np.stack([x["unpair"].sum(axis=0) / 2] * 2)
    rho, pi, _ = pairing.prepared_states(board, route, 1)
    x["formed_pairing"] = np.stack(
        [rho / rho.trace() - pi / pi.trace(), np.zeros((4, 4))]
    )
    waves = {
        (sg, odd): leg_kernels(ranging.generator(board, route, sg, 8, odd), route)
        for sg, odd in product((-1, 1), (False, True))
    }
    result = []
    for horizon in (2, 4, 6):
        ops = {
            sg: np.array(
                [arrivals(waves[sg, odd], 0.5, 2 * horizon) for odd in (False, True)]
            )
            for sg in (1, -1)
        }
        preparations = dict(x)
        if horizon == 6:
            preparations.update(matched.center(x))
        readings = {
            key: np.array(
                [np.einsum("bnij,bji->n", ops[sg], st).real for sg in (1, -1)]
            )
            for key, st in preparations.items()
        }
        for key, (excess, vacancy) in readings.items():
            result.append(
                dict(
                    horizon=horizon,
                    source=key,
                    charge_order=[1, -1],
                    readings=[excess.tolist(), vacancy.tolist()],
                    opposite_rectangles=rectangles(excess, vacancy),
                    gain_test=gain_feasibility(excess, vacancy),
                )
            )
    return dict(ports=route, edge_indices=PAIRS, tick=0.5, panels=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
