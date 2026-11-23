"""Complete, uniformly launched probe reads after unconstrained collar formation."""

import argparse
import hashlib
import json
from collections import deque
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..mobile_collar_source import measure as formation
from .measure import HERE

TIMES = (1.0, 2.0, 4.0)
TRIALS = 2**42
SEED = 390178
COHORTS = 864
TAIL = 0.01 / (2 * COHORTS)
READS = ("first", "second", "x_plus", "x_minus", "y_plus", "y_minus")
VECS = np.array([[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]])
VECS = VECS / np.linalg.norm(VECS, axis=1)[:, None]


def monitor_edges(board, route):
    home = route[1] // 12
    neighbours = {i: set() for i in range(board.cells)}
    for a, b in board.external:
        c, d = a // 12, b // 12
        neighbours[c].add(d)
        neighbours[d].add(c)
    depths = {home: 0}
    queue = deque([home])
    while queue:
        start = queue.popleft()
        for end in neighbours[start]:
            if end not in depths:
                depths[end] = depths[start] + 1
                queue.append(end)
    distant = min(
        (a, b)
        for a, b in board.external
        if {depths[a // 12], depths[b // 12]} == {2, 3}
    )
    return [tuple(route[1:3]), distant]


def prepared_states(board, route, s):
    h = formation.hamiltonian(board, route, s, 8)
    left, right = formation.starting_coins(route, s)
    initial = np.zeros(h.shape[0], complex)
    initial[4 * route[0] + 2 * left + right] = 1
    psi = expm_multiply(-4j * h, initial, traceA=0)
    rho = formation.reduce(psi)
    return (rho, np.kron(*formation.marginals(rho)), np.eye(4) / 4)


def backward_effects(h, edges):
    n = h.shape[0] // 4
    psi = np.zeros((n, 4, len(edges), 2, 4), complex)
    for j, (a, b) in enumerate(edges):
        psi[a, :, j, 0, :] = np.eye(4)
        psi[b, :, j, 1, :] = np.eye(4)
    waves = expm_multiply(
        1j * h, psi.reshape(4 * n, -1), start=0, stop=4, num=5, traceA=0
    )
    out = []
    for t in TIMES:
        vec = waves[int(t)].reshape(4 * n, len(edges), 2, 4)
        v = np.einsum("njac,ka->njkc", vec, VECS).reshape(
            n, 4, len(edges), len(VECS), 4
        )
        effects = np.einsum("xnjkc,xmjkc->jknm", v, v.conj(), optimize=True) / n
        np.testing.assert_allclose(
            np.trace(effects, axis1=-1, axis2=-2) / 4, 1 / n, rtol=0, atol=2e-12
        )
        out.append(effects)
    return np.array(out)


def intervals(passed):
    """Simultaneous Clopper–Pearson bounds, including exact endpoint events."""
    passed = np.asarray(passed)
    if np.any((passed < 0) | (passed > TRIALS)):
        raise ValueError("Binary counts outside the number of trials")
    lower = np.where(passed == 0, 0.0, beta.ppf(TAIL, passed, TRIALS - passed + 1))
    upper = np.where(passed == TRIALS, 1.0, beta.isf(TAIL, passed + 1, TRIALS - passed))
    return np.stack([lower, upper], axis=-1)


def run():
    rng = np.random.default_rng(SEED)
    rows = []
    for level in (0, 1):
        board = formation.federation(level)
        route = formation.region(board)
        edges = monitor_edges(board, route)
        states = {s: np.array(prepared_states(board, route, s)) for s in (-1, 1)}
        for q in (-1, 1):
            effects = backward_effects(formation.hamiltonian(board, route, q, 8), edges)
            for s, rhos in states.items():
                probs = np.einsum("tjkab,mba->mtjk", effects, rhos).real
                hits = rng.binomial(TRIALS, probs)
                bounds = intervals(hits)
                nominal = 1 / board.n
                delta = np.maximum(
                    0, np.maximum(bounds[..., 0] - nominal, nominal - bounds[..., 1])
                )
                signed = hits / TRIALS - nominal
                diff_lower = bounds[0, ..., 0] - bounds[1, ..., 1]
                diff_upper = bounds[0, ..., 1] - bounds[1, ..., 0]
                diff = np.maximum(0, np.maximum(diff_lower, -diff_upper))
                rows.append(
                    dict(
                        cells=board.cells,
                        ports=board.n,
                        edges=edges,
                        source=s,
                        probe=q,
                        plus_counts=hits.tolist(),
                        intervals=bounds.tolist(),
                        min_departure_from_uniform_by_arm=delta.max(
                            axis=(1, 2, 3)
                        ).tolist(),
                        significant_uniform_departures_by_arm=(delta > 0)
                        .sum(axis=(1, 2, 3))
                        .tolist(),
                        joint_product_max_lower_by_station=diff.max(
                            axis=(0, 2)
                        ).tolist(),
                        joint_product_significant_by_station=(diff > 0)
                        .sum(axis=(0, 2))
                        .tolist(),
                        max_abs_signed_joint_departure_by_station=np.abs(signed[0])
                        .max(axis=(0, 2))
                        .tolist(),
                    )
                )
    if len(rows) * 3 * len(TIMES) * 2 * len(READS) != COHORTS:
        raise ValueError("wrong cohort count")
    return dict(
        protocol_sha256=hashlib.sha256(
            (HERE / "FLIGHT_PROTOCOL.md").read_bytes()
        ).hexdigest(),
        seed=SEED,
        shots=TRIALS,
        tail_probability=TAIL,
        times=TIMES,
        axes=READS,
        arms=["generated", "product", "uniform"],
        rows=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    with args.output.open("x") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
    for r in result["rows"]:
        print(
            r["cells"],
            r["source"],
            r["probe"],
            r["significant_uniform_departures_by_arm"],
            r["joint_product_max_lower_by_station"],
            r["joint_product_significant_by_station"],
        )
