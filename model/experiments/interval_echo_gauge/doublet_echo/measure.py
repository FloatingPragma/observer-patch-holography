"""Two-aperture reciprocal echo of the registered pair doublet."""

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import brentq
from math import log, log1p, sqrt

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import (
    PairSector,
    pair_index,
)
from model.candidate_laws.interval_echo import currents
from model.executable_model import cell, federation
from model.experiments.interval_echo_gauge.doublet_flight.measure import endpoint


HERE = Path(__file__).resolve().parent
TIMES = (8, 32, 64)
ANGLES = (0, np.pi / 2, np.pi, 3 * np.pi / 2, np.pi / 3, 4 * np.pi / 5)
SHOTS = 2**60
J = np.array([[0, -1], [1, 0]])
X = np.array([[0, 1], [1, 0]])
INPUTS = np.array([[1, 0], [0, 1], [1, 1], [1, -1], [1, 1j], [1, -1j]], complex)
INPUTS /= np.linalg.norm(INPUTS, axis=1)[:, None]


def apparatus():
    board = federation(0)
    p, c = board.external[0]
    near = {
        n for n in board.neighbours[c] if n // 12 == c // 12 and board.degree[n] == 5
    }
    edges = [e for e in board.seams if set(e) <= near]
    left, right = next(
        (a, b) for a, b in combinations(edges, 2) if set(a).isdisjoint(b)
    )
    sector = PairSector(board)
    ports = [sector.index(c, a) for a in (*left, *right)]
    local = PairSector(cell())
    home = pair_index(board.n, local.aa, local.bb)
    frame = endpoint(p)["frame"][:, 1:]
    return board, sector, (p, c, left, right), ports, home, frame


def transfer(sector, home, frame, ports, kappa):
    initial = np.zeros((sector.P, 2))
    initial[home] = frame
    h = sector.hamiltonian(kappa)
    states = expm_multiply(-1j * h, initial, start=0, stop=64, num=9)
    return [states[t // 8, ports] for t in TIMES]


def source_words(board, axis, frame):
    near = {n for n in board.neighbours[axis] if n // 12 == 0 and board.degree[n] == 5}
    a, b = next(edge for edge in board.seams if set(edge) <= near)
    small = PairSector(cell())
    rows = [small.index(axis, site) for site in (a, b)]
    # The coherent difference port applies (tag - I)/2 on all input amplitudes.
    alpha = (endpoint(axis)["phase"] - 1) / 2
    return (a, b), frame[rows].T, alpha


def turnaround(kappa, theta, signed=True):
    local = law.local_term_parts()
    h = local[0] + kappa * local[1]
    single = [law.SRC.index(pair) for pair in ((1, 0), (0, 1))]
    swap = expm(-1j * np.pi * h)[:, single]
    expected = np.eye(len(h))[:, single[::-1]]
    np.testing.assert_allclose(swap, expected, atol=2e-14)
    sign = (
        np.diag(np.exp(-1j * kappa * (np.pi / kappa) * np.array([0, 1])))
        if signed
        else np.eye(2)
    )
    gate = swap[single] @ sign
    tag = np.exp(-1j * kappa * ((-theta) % (2 * np.pi) / kappa))
    return gate, tag * gate


def experiment_matrices(t, kappa):
    pair = t.reshape(2, 2, 2)
    operations = {}
    for n, theta in enumerate(ANGLES):
        gates = turnaround(kappa, theta)
        parts = [m.T @ g @ m for m, g in zip(pair, gates, strict=True)]
        operations[f"phase{n}"] = [(sum(parts), [0, 1])]
    straight, _ = turnaround(kappa, 0, signed=False)
    operations["ordinary_swap"] = [(sum(m.T @ straight @ m for m in pair), [0, 1])]
    for i, m in enumerate(pair):
        operations[f"only{i}"] = [(m.T @ J @ m, [i])]
    operations["which"] = operations["only0"] + operations["only1"]
    return operations


def probabilities(t, branches, transform, alpha, source):
    prepared = alpha * transform @ source
    target = J @ source
    good = wrong = filtered = midpoint = 0.0
    for matrix, paths in branches:
        returned = alpha * matrix @ prepared
        decoded = transform.T @ returned
        good += abs(np.vdot(target, decoded)) ** 2
        wrong += (
            abs(
                np.vdot(
                    np.array([-target[1].conjugate(), target[0].conjugate()]), decoded
                )
            )
            ** 2
        )
        filtered += np.linalg.norm(returned) ** 2
        midpoint += sum(
            np.linalg.norm(t[2 * i : 2 * i + 2] @ prepared) ** 2 for i in paths
        )
    initial = np.linalg.norm(prepared) ** 2
    values = [
        good,
        wrong,
        max(0.0, filtered - good - wrong),
        max(0.0, midpoint - filtered),
        max(0.0, initial - midpoint),
        max(0.0, 1 - initial),
    ]
    np.testing.assert_allclose(sum(values), 1, atol=2e-14)
    return values


def stable_binomial(rng, n, p):
    """Direct small-mean inversion avoids subtracting a tiny p from one."""
    if p == 0:
        return np.zeros_like(n)
    if p >= 1e-10 or np.max(n) * p > 30:
        return rng.binomial(n, p)
    u = rng.random(len(n))
    term = np.exp(n * np.log1p(-p))
    total = term.copy()
    out = np.zeros_like(n)
    j = 0
    while np.any(pending := (u > total)):
        j += 1
        out[pending] += 1
        term[pending] *= ((n[pending] - j + 1) * p / j) / (1 - p)
        total[pending] += term[pending]
        if j > 160:
            raise ArithmeticError("Rare-bin inversion exceeded its numerical range")
    return out


def acquire(prob, batches, seed):
    """Draw in integer batches; preserve Python-integer totals beyond 64 bits."""
    rng = np.random.default_rng(seed)
    # Put the large complement last so rounding of unit-minus-tiny cannot erase a rare click.
    order = np.argsort(prob)
    p = np.array(prob)[order]
    p[-1] = 1 - sum(p[:-1])
    strata = np.zeros((batches, len(p)), dtype=np.int64)
    remaining = np.full(batches, SHOTS, dtype=np.int64)
    for j in range(len(p) - 1):
        draws = stable_binomial(rng, remaining, p[j] / p[j:].sum())
        strata[:, j] = draws
        remaining -= draws
    strata[:, -1] = remaining
    first = [int(v) for v in strata[0, np.argsort(order)]]
    # Summing row chunks in float would lose hundreds of counts in common bins.
    low = (strata & 0xFFFFFFFF).sum(axis=0)
    high = (strata >> 32).sum(axis=0)
    sums = [int(low[j]) + (int(high[j]) << 32) for j in np.argsort(order)]
    assert sum(sums) == SHOTS * batches
    return sums, first


def interval(k, n):
    # E[exp(tX)] <= exp(n p (exp(t)-1)) for every binomial X.
    # Inverting its Chernoff tails in the mean avoids beta quantiles at huge n.
    budget = log(2 * 2048 / 0.01)
    if k == 0:
        return np.array([0.0, min(1.0, budget / n)])

    def boundary(mu):
        return (mu - k) - k * log1p((mu - k) / k) - budget

    # The low bracket is strictly positive even for a single observed count.
    left = max(k * np.exp(-budget / k - 2.0), np.finfo(float).tiny)
    right = k + 2 * budget + 2 * sqrt(k * budget)
    while boundary(right) < 0:
        right *= 2
    low = brentq(boundary, left, k)
    high = brentq(boundary, k, right)
    return np.array([max(0.0, low * (1 - 1e-12) / n), min(1.0, high * (1 + 1e-12) / n)])


def held_test(panels, shots):
    pooled = {
        name: sum(row["counts"][0] for row in rows) for name, rows in panels.items()
    }
    n = 6 * shots
    intervals = [interval(pooled[f"phase{i}"], n) for i in range(len(ANGLES))]
    predictions = []
    for i in (4, 5):
        theta = ANGLES[i]
        coefficients = [
            (1 + np.cos(theta)) / 2,
            np.sin(theta) / 2,
            (1 - np.cos(theta)) / 2,
            -np.sin(theta) / 2,
        ]
        products = [
            np.sort(c * iv) for c, iv in zip(coefficients, intervals[:4], strict=True)
        ]
        prediction = np.maximum(0.0, np.sum(products, axis=0))
        predictions.append(
            dict(
                prediction=prediction.tolist(),
                held=intervals[i].tolist(),
                compatible=bool(
                    max(prediction[0], intervals[i][0])
                    <= min(prediction[1], intervals[i][1])
                ),
            )
        )
    return dict(
        pooled_target_counts=pooled,
        calibration_intervals=[v.tolist() for v in intervals[:4]],
        held=predictions,
        dark_vs_bright_resolved=bool(
            intervals[0][0] > intervals[2][1] or intervals[2][0] > intervals[0][1]
        ),
        quadratures_resolved=bool(
            intervals[1][0] > intervals[3][1] or intervals[3][0] > intervals[1][1]
        ),
    )


def summarize(matrices, kappa, home_transfer, alpha):
    results = []
    for t, matrix in zip(TIMES, matrices, strict=True):
        operations = experiment_matrices(matrix, kappa)
        determinants = np.linalg.det(matrix.reshape(2, 2, 2))
        calibrations = [("plane", np.eye(2), 1.0, 65 if t >= 32 else 1)]
        if t == 32:
            calibrations.append(("sharp", home_transfer, alpha, 2**18))
        for kind, transform, factor, batches in calibrations:
            panels = {}
            for j, (name, branches) in enumerate(operations.items()):
                rows = []
                for s, source in enumerate(INPUTS):
                    p = probabilities(matrix, branches, transform, factor, source)
                    seed = [74115, int(kappa), t, int(kind == "sharp"), j, s]
                    if name == "which":
                        path_bins = [
                            probabilities(matrix, [branch], transform, factor, source)[
                                :4
                            ]
                            for branch in branches
                        ]
                        raw, single = acquire(
                            path_bins[0] + path_bins[1] + p[4:], batches, seed
                        )
                        row = {
                            "labelled_path_counts": [raw[:4], raw[4:8]],
                            "first_labelled_path_counts": [single[:4], single[4:8]],
                        }

                        def discard_path_label(counts):
                            return [
                                a + b
                                for a, b in zip(counts[:4], counts[4:8], strict=True)
                            ] + counts[8:]

                        row["counts"] = discard_path_label(raw)
                        row["first_batch"] = discard_path_label(single)
                    else:
                        counts, first = acquire(p, batches, seed)
                        row = dict(counts=counts, first_batch=first)
                    rows.append(row)
                panels[name] = rows
            results.append(
                dict(
                    exposure=t,
                    instrument=kind,
                    batches_per_setting=batches,
                    mode_path_probabilities=(abs(determinants) ** 2).tolist(),
                    native_path_phase=float(
                        np.angle(determinants[1] / determinants[0])
                    ),
                    coherent_visibility=float(
                        2 * abs(np.prod(determinants)) / sum(abs(determinants) ** 2)
                    ),
                    complete_counts=panels,
                    forecast=held_test(panels, SHOTS * batches),
                )
            )
    return results


def run():
    board, sector, addresses, ports, home, frame = apparatus()
    edge, transform, alpha = source_words(board, addresses[0], frame)
    rows = []
    for kappa in (currents.KAPPA, currents.KAPPA / 2):
        matrices = transfer(sector, home, frame, ports, kappa)
        rows.append(
            dict(kappa=kappa, runs=summarize(matrices, kappa, transform, alpha))
        )
    return dict(
        protocol_sha256=hashlib.sha256((HERE / "PROTOCOL.md").read_bytes()).hexdigest(),
        board_cells=board.cells,
        addresses=addresses,
        launch_edge=edge,
        return_efficiency_from_sharp=float(
            abs(alpha) ** 4 * abs(np.linalg.det(transform)) ** 2
        ),
        batch_size=SHOTS,
        input_order=["0", "1", "+", "-", "+i", "-i"],
        bins=[
            "J_target",
            "orthogonal_logical",
            "other_sharp_words_after_plane_pass",
            "failed_home_or_plane_read",
            "failed_midpoint",
            "failed_source_filter",
        ],
        arms=rows,
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    with ap.parse_args().output.open("x") as stream:
        json.dump(run(), stream, indent=1, allow_nan=False)
        stream.write("\n")
