"""Read the initial acceleration of the same counter under ordinary and regional sources."""

import argparse
from fractions import Fraction
from hashlib import sha256
from math import factorial
import json

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.executable_model import federation
from model.experiments.interval_echo_einstein_station import volume
from . import word_routes as routes

HERE = routes.HERE
SEED = 28115739
STEP = Fraction(1, 1024)
WEIGHTS = (
    Fraction(-205, 72),
    Fraction(16, 5),
    Fraction(-2, 5),
    Fraction(16, 315),
    Fraction(-1, 280),
)


def select(board, inside, common):
    a, b = min(
        (a, b)
        for edge in board.external
        for a, b in (edge, edge[::-1])
        if inside[a] and not inside[b]
    )
    s = 12 * (a // 12) + int(np.flatnonzero(common[a % 12])[0])
    assert (min(s, b), max(s, b)) not in board.seams
    return s, a, b


def contact_terms():
    law = volume.candidate.native_law
    e, c = law.local_term_parts()
    h = e + 8 * c
    # All written memories on this control are sharp zero modes, including
    # when two moving loads arrive on opposite ends of a written edge.
    mark = np.kron(volume.label(0, 1)[:, None], np.eye(law.N_SRC))
    np.testing.assert_allclose(h @ mark, 0, atol=1e-13)
    blank = np.eye(law.DIM)[:, : law.N_SRC]
    np.testing.assert_allclose(h @ blank, blank @ law.seam_term(8), atol=1e-13)
    block = law.seam_term(8)

    def energy(x, y):
        i = law.SRC.index((x, y))
        return block[i, i]

    return {q: float(energy(1, q) - energy(1, 0) - energy(0, q)) for q in (1, 2)}


def two_loads(one, adjacency, contact, q):
    # Ordered tensor coordinates restricted to distinct occupied ports.
    # For identical particles the symmetric subspace is invariant. Two
    # equal neighbors have no exchange cost; this is included in contact.
    n = one.shape[0]
    generator = sparse.kronsum(one, one, format="csr")
    row, col = adjacency.nonzero()
    pairs = row * n + col
    diag = np.zeros(n * n)
    diag[pairs] = contact
    generator = generator + sparse.diags(diag)
    if q != 1:
        generator += sparse.csr_matrix(
            (np.full(len(pairs), -0.5), (pairs, col * n + row)), shape=generator.shape
        )
    x, y = np.indices((n, n)).reshape(2, -1)
    keep = np.flatnonzero(x != y)
    return generator[keep][:, keep].tocsr(), x[keep], y[keep]


def prepare(x, y, s, a, b, q):
    def branch(v):
        return (
            (((x == s) & (y == v)) | ((x == v) & (y == s))).astype(float) / np.sqrt(2)
            if q == 1
            else ((x == s) & (y == v)).astype(float)
        )

    left, right = branch(a), branch(b)
    return np.column_stack((left + right, left - right)) / np.sqrt(2)


def read(generator, starts, score, rng):
    # Removing a scalar leaves every occupation expectation invariant.
    norm = float(np.asarray(abs(generator).sum(axis=1)).max())
    h = float(STEP)
    remainder = (
        max(score)
        * (2 * norm) ** 10
        * h**8
        * sum(float(abs(w)) * j**10 for j, w in enumerate(WEIGHTS))
        / factorial(10)
    )
    gen = -1j * generator
    path = expm_multiply(
        gen, starts, start=0, stop=4 * h, num=5, traceA=gen.diagonal().sum()
    )
    direct = 2 * ((generator @ starts) ** 2 * score[:, None]).sum(axis=0) - 2 * (
        starts * score[:, None] * (generator @ (generator @ starts))
    ).sum(axis=0)
    readings = []
    # Centering by a constant reduces the interval width of the count read.
    centered = np.arange(int(max(score)) + 1) - 1
    for k in range(2):
        tape = []
        mid = width = 0.0
        for j, w in enumerate(WEIGHTS):
            prob = np.bincount(
                score, weights=abs(path[j, :, k]) ** 2, minlength=len(centered)
            )
            assert abs(prob.sum() - 1) < 1e-12
            if j:
                counts = routes.lab.counts(prob, rng)
                tape.append(counts)
                freqs = np.array(counts) / routes.lab.SHOTS
                err = routes.reg.src.probability_error(freqs, routes.lab.SHOTS) @ abs(
                    centered
                )
                value = freqs @ centered
            else:
                err = 0.0
                # At release each known preparation gives precisely half of
                # its counts at each of two adjacent occupation values.
                value = round(2 * (prob @ centered)) / 2
            mid += float(w) * value / h**2
            width += abs(float(w)) * err / h**2
        width += remainder + 2e-7
        readings.append(
            {
                "counts": tape,
                "held_acceleration": [mid - width, mid + width],
                "exact_acceleration": float(direct[k]),
            }
        )
    return {
        "centered_norm_bound": norm,
        "truncation_bound": float(remainder),
        "phases": readings,
    }


def energy_standards(one, board, site):
    # Compare the same prepared packet under native evolution and a bypass.
    # For the kinetic read the reference branch evolves in the lowest mode.
    laplacian = 2 * one
    eigen, modes = np.linalg.eigh(laplacian.toarray())
    assert abs(eigen[0]) < 1e-12 and 0 < eigen[1] < 1
    init = np.stack((modes[:, 0], np.eye(1, board.n, site).ravel())).T
    rng = np.random.default_rng(SEED + 1)
    h = 1 / 32768

    def acquire(moment):
        counts = routes.lab.counts([(1 + moment) / 2, (1 - moment) / 2], rng)
        f = counts[0] / routes.lab.SHOTS
        err = 2 * routes.reg.src.probability_error(np.asarray(f), routes.lab.SHOTS)
        return counts, routes.reg.band(2 * f - 1 - err, 2 * f - 1 + err)

    results = {}
    for q in (1, 2):
        degree = laplacian.diagonal()
        standard = 8 * q * q * (degree + 6 - board.degree)
        np.testing.assert_array_equal(standard, np.full(board.n, 40 * q * q))
        native = one + sparse.diags(standard)
        gen = -1j * native
        evolved = expm_multiply(h * gen, init, traceA=h * gen.diagonal().sum())
        norm = float(abs(native).sum(axis=1).max())
        error = norm**3 * h * h / 6
        reads = {}
        for j, name in enumerate(("rest", "localized")):
            tape, word = acquire(float(-(init[:, j] @ evolved[:, j]).imag))
            reads[name] = {
                "counts": tape,
                "energy": routes.bounds(word / h + routes.iv.mpf([-error, error])),
            }
        # A recombined, equal superposition of the lowest and first
        # excited mode measures their relative phase, canceling the rest gap.
        pair = modes[:, :2]
        evolved = expm_multiply(gen / 2, pair, traceA=gen.diagonal().sum() / 2)
        phases = np.sum(pair * evolved, axis=0)
        tape, word = acquire(float(-(phases[0].conjugate() * phases[1]).imag))
        # The calibrated one-mode branch lies within the principal sine interval.
        gap = 2 * routes.iv.atan2(word, routes.iv.sqrt(1 - word**2))
        reads["first_gap"] = {
            "counts": tape,
            "energy": routes.bounds(gap),
            "gap_per_laplacian_eigenvalue": routes.bounds(gap / eigen[1]),
        }
        results[str(q)] = reads
    return {
        "source_seed": SEED + 1,
        "bypass_time": h,
        "mode_time": 0.5,
        "clock_laplacian_first_eigenvalue": float(eigen[1]),
        "species": results,
    }


def run():
    routes.iv.dps = 35
    board = federation(0)
    protocol = json.loads((HERE / "word_routes.json").read_bytes())
    inside = np.repeat(np.isin(np.arange(board.cells), protocol["source_ball"]), 12)
    code, common = routes.codebook()
    cal, _ = routes.calibrate(code)
    adjacency, deg = routes.graph(board, cal)(-1, -1, 0, 0)
    one = (sparse.diags(deg, dtype=float) - adjacency) * 0.5
    s, a, b = select(board, inside, common)
    assert adjacency[s, a] and adjacency[a, b] and not adjacency[s, b]
    assert not set(adjacency[a].indices) & set(adjacency[b].indices)
    weights = contact_terms()
    rng = np.random.default_rng(SEED)
    rows = {}
    single = np.eye(board.n)[:, [a, b]]
    plus_minus = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    for q in (1, 2):
        blank = read(one, single @ plus_minus, inside.astype(int), rng)
        h, x, y = two_loads(one, adjacency, weights[q], q)
        occupied = read(
            h,
            prepare(x, y, s, a, b, q),
            inside[x].astype(int) + inside[y].astype(int),
            rng,
        )

        def contrast(data):
            p, m = data["phases"]
            return routes.iv.mpf(p["held_acceleration"]) - routes.iv.mpf(
                m["held_acceleration"]
            )

        measured = (contrast(occupied) - contrast(blank)) / 2
        prediction = -weights[q] / 2
        exact = sum(
            sign * arm["exact_acceleration"] / 2
            for data, ds in ((occupied, 1), (blank, -1))
            for arm, sign in zip(data["phases"], (ds, -ds), strict=True)
        )
        assert (
            abs(exact - prediction) < 1e-11
            and prediction in measured
            and measured.a > 0
        )
        rows[str(q)] = {
            "empty": blank,
            "ordinary": occupied,
            "contact_energy_shift": weights[q],
            "expected": prediction,
            "measured": routes.bounds(measured),
        }
    return {
        "seed": SEED,
        "starts_per_exposure": routes.lab.SHOTS,
        "step": str(STEP),
        "weights": list(map(str, WEIGHTS)),
        "ordinary_site": s,
        "launch_sites": [a, b],
        "rest_gap": 40,
        "phase_energy_standards": energy_standards(one, board, s),
        "localized_energy_above_vacuum": float(40 + deg[s] / 2),
        "source_cells": protocol["source_ball"],
        "modular_initial_acceleration_contrast": 0,
        "unit_and_double": rows,
        "inputs": {
            k: sha256((HERE / k).read_bytes()).hexdigest()
            for k in ("INITIAL_ACCELERATION.md", "word_routes.json")
        },
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    routes.iv.dps = 35
    data = run()
    path = HERE / "initial_acceleration.json"
    if args.write:
        path.write_text(json.dumps(data, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), data)
    for k, z in data["unit_and_double"].items():
        print(k, z["expected"], z["measured"])
