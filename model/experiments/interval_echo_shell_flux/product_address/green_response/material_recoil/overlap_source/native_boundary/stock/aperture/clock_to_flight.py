"""Predict cold flight from the local insertion clock calibrated on ordinary matter."""

import argparse
import json
from hashlib import sha256

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from . import mediated_release as lab

HERE = lab.HERE


def prepare(p):
    n = len(p)
    return (1 - np.eye(n)) * np.sqrt(p / (n - 1))


def kinetic_force(lap, potential, p, shells, gain):
    h = lap / 2 + sparse.diags(gain * potential, format="csr")
    wave = prepare(p)
    velocity = h @ wave
    accel = 2 * np.sum(velocity**2 - wave * (h @ velocity), axis=1)
    return np.bincount(shells, weights=accel).cumsum()


def forecast(lap, potential, p, shells, gain):
    h = lap / 2 + sparse.diags(gain * potential, format="csr")
    wave = prepare(p).astype(complex)
    laws = []
    last = 0
    for time in lab.RELEASE_TIMES:
        if time != last:
            wave = expm_multiply((-1j * (time - last)) * h, wave)
        last = time
        laws.append(np.bincount(shells, weights=(abs(wave) ** 2).sum(axis=1)).cumsum())
    return np.asarray(laws) - laws[0]


def potential(board, mask, p, g):
    m = np.bincount(np.array(mask, dtype=int).ravel(), minlength=board.n)
    # The observer's clock includes a kinetic diagonal d_live/2. Removing it
    # from the total frequency change puts half a unit per mark back.
    out = []
    for x, near in enumerate(board.neighbours):
        live = [v for v in near if tuple(sorted((x, v))) not in mask]
        out.append(-(g - 0.5) * m[x] - g * p[live].sum() / (1 - p[x]))
    return np.array(out)


def exterior_residual(data, board, shells, coefficient):
    """Bound true minus scalar change using just counted ordinary addresses.

    For a boundary edge inside u, outside v, the potential part of the
    scalar acceleration is C*(1-p_u-p_v)*(s_u/(1-p_u)-s_v/(1-p_v)),
    with s_x the probability of a neighbor and C=coefficient/(N-1).
    Separate its linear and quadratic terms before propagating count errors.
    """
    n = board.n
    lap = board.laplacian().astype(np.int64)
    wave = 1 - np.eye(n, dtype=np.int64)
    b = lap @ wave
    kinetic = b * b - wave * (lap @ b)
    by_shell = np.zeros((shells.max() + 1, n), dtype=np.int64)
    np.add.at(by_shell, shells, kinetic)
    kernel = (np.asarray(data["ordinary_integer_kernel"]) - by_shell.cumsum(axis=0)) / (
        2 * (n - 1)
    )
    C = coefficient / (n - 1)
    cases = []
    for r in range(1, shells.max() + 1):
        edges = [
            (u, v) if shells[u] <= r else (v, u)
            for u, v in board.seams
            if (shells[u] <= r) != (shells[v] <= r)
        ]
        linear = kernel[r].copy()
        for u, v in edges:
            linear[list(board.neighbours[u])] -= C
            linear[list(board.neighbours[v])] += C

        def enclosure(counts):
            p = np.asarray(counts) / sum(counts)
            sums = board.adjacency() @ p
            err = np.sqrt(lab.RISK / (2 * sum(counts)))
            lo = np.maximum(0, p - err)
            hi = np.minimum(1, p + err)
            s_lo = np.maximum(0, sums - err)
            s_hi = np.minimum(1, sums + err)
            assert max(hi) < 1

            def term(v, u):
                return np.array(
                    [lo[v] * s_lo[u] / (1 - lo[u]), hi[v] * s_hi[u] / (1 - hi[u])]
                )

            bounds = linear @ p + np.array([-1, 1]) * np.ptp(linear) * err
            for u, v in edges:
                bounds += C * lab.delta(term(v, u), term(u, v))
            return bounds

        change = lab.delta(
            enclosure(data["source_counts"]["written"]),
            enclosure(data["source_counts"]["blank"]),
        )
        cases.append((change + np.array([-1, 1]) * 1e-10).tolist())
    return cases


def run():
    paths = [
        HERE / f
        for f in (
            "mediated_release.json",
            "ordinary_release.json",
            "CLOCK_TO_FLIGHT.md",
        )
    ]
    data, standard = [json.loads(path.read_bytes()) for path in paths[:2]]
    board = lab.cold.federation(0)
    marks, _ = lab.cold.read_map(board, data["spectral_tape"])
    edge = marks[0]
    assert list(edge) == standard["seam"]
    shells = lab.cold.distances(board, edge)
    g = np.mean(data["decoded"]["ordinary_neighbor_gain"])
    weights = {k: np.asarray(c) / sum(c) for k, c in data["source_counts"].items()}
    p = np.zeros(board.n)
    p[list(edge)] = 0.5
    _, lap, _ = lab.cold.motion(board, [], 2, lab.K)
    v = potential(board, [], p, g)
    native = np.cumsum(
        np.array(standard["arms"]["ordinary"]["acceleration"])
        - standard["arms"]["blank"]["acceleration"]
    )
    # Density commutes with a scalar potential, hence the initial
    # acceleration is affine in its strength for every initial wave.
    model = kinetic_force(lap, v, p, shells, 1) - kinetic_force(lap, v, p, shells, 0)
    gain = native[1] / model[1]
    assert gain > 0
    out = {}
    for key, lam in (("unit", 1.0), ("ordinary_calibrated", gain)):
        a = []
        laws = []
        for medium, mask in (("blank", []), ("written", marks)):
            _, lap, _ = lab.cold.motion(board, mask, 2, lab.K)
            p = weights[medium]
            field = potential(board, mask, p, g)
            a.append(kinetic_force(lap, field, p, shells, lam))
            laws.append(forecast(lap, field, p, shells, lam))
        prediction = laws[1] - laws[0]
        observed = np.asarray(data["decoded"]["inside_changes"])
        out[key] = dict(
            gain=lam,
            acceleration=(a[1] - a[0]).tolist(),
            finite_probability_change=prediction.tolist(),
            bins_outside_position_count_band=np.argwhere(
                abs(observed - prediction) > data["decoded"]["two_arm_change_error"]
            ).tolist(),
        )
    return dict(
        sha256={p.name: sha256(p.read_bytes()).hexdigest() for p in paths},
        ordinary_first_exterior=native[1],
        unit_first_exterior=model[1],
        comparisons=out,
        exterior_native_minus_scalar_interval=exterior_residual(
            data, board, shells, gain * g
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record = run()
    path = HERE / "clock_to_flight.json"
    if args.write:
        path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    else:
        lab.cold.apparatus.check(json.loads(path.read_bytes()), record)
    for k, v in record["comparisons"].items():
        print(k, v["gain"], v["acceleration"][1], v["bins_outside_position_count_band"])
