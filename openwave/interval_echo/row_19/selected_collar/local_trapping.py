"""Measure lasting ordinary occupation of cell modes after all receivers participate."""

from fractions import Fraction as F
from itertools import combinations, product
from math import factorial, sqrt

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.special import jv
import sympy as sp

from . import stationary_polarization as steady

exp = steady.exp
local = exp.local
SEED = 28122917
T = F(1, 2**19)


def projector(internal, gates):
    h = np.zeros((12, 12), dtype=object)
    for a, b in internal:
        h[a, a] += 1
        h[b, b] += 1
        h[a, b] -= 1
        h[b, a] -= 1
    block = np.eye(12, dtype=object)[:, list(gates)]
    seq = [block]
    for _ in range(11):
        seq.append(h @ seq[-1])
    kernel = sp.Matrix(np.column_stack(seq).T.tolist()).to_DM().nullspace().to_Matrix()
    proj = kernel.T * (kernel * kernel.T).inv() * kernel
    assert (
        proj * proj == proj
        and proj * sp.Matrix(h.tolist()) == sp.Matrix(h.tolist()) * proj
    )
    assert all(proj[j, j] == 0 for j in gates)
    assert proj * sp.ones(12, 1) == sp.zeros(12, 1)
    return np.array([[F(x) for x in row] for row in proj.tolist()], dtype=object)


def sign_mean(G, S):
    if not S:
        return F(1)
    a = sp.Matrix([[2 * G[i, j] - int(i == j) for j in S] for i in S])
    z = a.inv() * sp.ones(len(S), 1) * 2
    return 1 - sum(F(G[0, i]) * F(v) for i, v in zip(S, z, strict=True))


def arrangements(board, steps, codes, internal):
    projectors = {}
    out = []
    for c, z in enumerate(steps):
        gates = tuple(int(g) for g in np.flatnonzero(z != c))
        if gates not in projectors:
            projectors[gates] = projector(internal, gates)
        P = projectors[gates]
        partners = [next(a + b - g for a, b in codes[c] if g in (a, b)) for g in gates]
        region = (c, *(int(z[g]) for g in gates))
        scores = []
        for bits in product((0, 1), repeat=4):
            active = [j for j, b in zip(partners, bits[1:], strict=True) if not b]
            scores.append(sum(P[x, y] for x in active for y in active) * bits[0])
        assert min(scores) >= 0 and max(scores) > 0
        out.append((region, partners, P, scores))
    return out


def transfer_matrix(board, internal):
    edges = [
        (i, j)
        for i, j in board.seams
        if i // 12 != j // 12 or (i % 12, j % 12) in internal
    ]
    a, b = np.array(edges).T
    deg = np.bincount(np.r_[a, b], minlength=board.n)
    return (
        coo_matrix(
            (
                np.r_[deg, -np.ones(2 * len(a))],
                (np.r_[np.arange(board.n), a, b], np.r_[np.arange(board.n), b, a]),
            ),
            shape=(board.n, board.n),
        ).tocsr()
        / 2
    )


def local_rotation(desired, internal):
    # A real Givens rotation uses a calibrated blank partial swap and port
    # phases. Nonadjacent modes are brought together along a blank route
    # with full swaps and returned along the inverse permutation.
    rem = desired.copy()
    n = len(desired)
    total = np.eye(n)
    graph = nx.Graph(internal)
    swaps = rotations = 0
    for col in range(n):
        for row in range(n - 1, col, -1):
            a, b = rem[row - 1, col], rem[row, col]
            if abs(b) < 1e-14:
                continue
            norm = np.hypot(a, b)
            c = a / norm
            sine = b / norm
            # The native off-diagonal amplitude is i sin(angle); phase
            # settings on the inputs and outputs convert to real rotation.
            calibrated = np.array([[c, 1j * sine], [1j * sine, c]])
            rot = np.real(np.diag([1, 1j]) @ calibrated @ np.diag([1, -1j]))
            ids = [row - 1, row]
            rem[ids] = rot @ rem[ids]
            total[ids] = rot @ total[ids]
            swaps += 2 * (nx.shortest_path_length(graph, *ids) - 1)
            rotations += 1
    reader = total.T @ np.diag(np.sign(np.diag(rem)))
    np.testing.assert_allclose(reader, desired, rtol=0, atol=5e-13)
    return reader, {"full_swaps": swaps, "partial_swaps": rotations}


def port_reader(P, internal):
    values, vecs = np.linalg.eigh(np.array(P, float))
    reader, circuit = local_rotation(vecs.T, internal)
    selected = reader[values > 0.5]
    assert np.max(abs(selected.T @ selected - np.array(P, float))) < 3e-13
    return selected, {
        **circuit,
        "selected_modes": int(sum(values > 0.5)),
    }


def isolated_path_bounds(board, items, codes, internal, hops):
    n = board.n
    active = [
        (a, b)
        for a, b in board.seams
        if a // 12 != b // 12 or (a % 12, b % 12) in internal
    ]
    output = []
    order = 52
    error = F(8) ** (order + 1) / factorial(order + 1)
    den = 2**order * factorial(order)
    for cell, (_, _, P, _) in enumerate(items):
        ignored = {(12 * cell + i, 12 * cell + j) for i, j in codes[cell]}
        edges = np.array([e for e in active if e not in ignored])
        a, b = edges.T
        deg = np.bincount(np.r_[a, b], minlength=n)
        diag = (deg - 16 * (np.arange(n) // 12 == cell)).astype(object)
        power = np.ones(n, dtype=object)
        re = np.zeros_like(power)
        im = np.zeros_like(power)
        for k in range(order + 1):
            coeff = (-1) ** (k // 2) * (den // (2**k * factorial(k)))
            if k % 2:
                im -= coeff * power
            else:
                re += coeff * power
            z = diag * power
            np.add.at(z, a, -power[b])
            np.add.at(z, b, -power[a])
            power = z
        x = re[12 * cell : 12 * cell + 12]
        y = im[12 * cell : 12 * cell + 12]
        mass = (x @ P @ x + y @ P @ y) / (n * den**2)
        path = F(1, 32) * (F(31, 32) ** 2 * F(9, 144)) ** hops[cell]
        lower = path * (mass - 2 * error - error**2) / 144
        assert lower > 0
        output.append(
            {
                "hop": hops[cell],
                "path_probability": str(path),
                "accepted_compact": float(mass),
                "unconditional_excess_lower": float(lower) * (1 - 1e-13),
            }
        )
    return output


def propagation(board, steps, items, codes, internal, rng):
    # Each trajectory is a native courier start. Conditional projection
    # probabilities allow the much smaller correlated excess to be estimated.
    h = transfer_matrix(board, internal)
    mates = np.empty(board.n, int)
    for c, matching in enumerate(codes):
        for i, j in matching:
            mates[12 * c + i] = 12 * c + j
            mates[12 * c + j] = 12 * c + i
    projections = np.array([np.asarray(z[2], float) for z in items])
    meters = [port_reader(z[2], internal) for z in items]
    read = np.array([m[0] for m in meters])
    event_counts = np.zeros(board.cells + 1, dtype=np.int64)
    shots = 2**20
    total = 8192
    batch = 128
    first = np.zeros(board.cells)
    second = np.zeros(board.cells)
    score_sum = np.zeros(board.cells)
    for _ in range(total // batch):
        flags = np.zeros((board.cells, batch), bool)
        stations = np.zeros(batch, int)
        live = np.arange(batch)
        while len(live):
            flags[stations[live], live] ^= True
            live = live[rng.integers(32, size=len(live)) != 0]
            stations[live] = steps[stations[live], rng.integers(12, size=len(live))]
        for c, (region, other, P, scores) in enumerate(items):
            ix = np.array([8, 4, 2, 1]) @ flags[list(region)].astype(int)
            score_sum[c] += np.asarray(scores, float)[ix].sum()
        mark = np.repeat(flags, 12, axis=0)

        def apply(x):
            return h @ x - mark * (8 * x + (x - x[mates]) / 2)

        z = np.full((board.n, batch), 1 / sqrt(board.n))

        # Chebyshev polynomials on the enclosing spectral interval [-8,5].
        def scaled(x):
            return (apply(x) + 1.5 * x) / 6.5

        last = z
        cur = scaled(z)
        wave = jv(0, 6.5) * z + (-2j) * jv(1, 6.5) * cur
        for k in range(2, 43):
            last, cur = cur, 2 * scaled(cur) - last
            wave += 2 * (-1j) ** k * jv(k, 6.5) * cur
        shaped = wave.reshape(board.cells, 12, batch)
        compact = np.einsum("cij,cjb->cib", projections, shaped)
        chance = np.sum(abs(compact) ** 2, axis=1)
        assert np.max(abs(np.sum(abs(wave) ** 2, axis=0) - 1)) < 1e-12
        observed = np.sum(abs(np.einsum("cij,cjb->cib", read, shaped)) ** 2, axis=1)
        np.testing.assert_allclose(chance, observed, atol=5e-14)
        for prob in observed.T:
            event_counts += rng.multinomial(shots, np.r_[prob, 1 - prob.sum()])
        # Check the first three coefficients directly using native action.
        if _ == 0:
            third = apply(apply(apply(z))).reshape(board.cells, 12, batch)
            predicted = np.zeros_like(third)
            for c, (region, partners, _, _) in enumerate(items):
                for arm, neighbor in zip(partners, region[1:], strict=True):
                    predicted[c, arm] = (
                        -2 * (flags[c] & ~flags[neighbor]) / sqrt(board.n)
                    )
            np.testing.assert_allclose(
                np.einsum("cij,cjb->cib", projections, third - predicted), 0, atol=1e-13
            )
        first += chance.sum(axis=1)
        second += (chance**2).sum(axis=1)
    mean = first / total
    variance = (second - total * mean**2) / (total - 1)
    # Independent trajectory estimates of bounded [0,1] conditional chances.
    # Individual nonnegative layout contributions also give a rigorous finite
    # lower bound for the complete ensemble using one explicit path below.
    return {
        "paths": total,
        "shots_per_path": shots,
        "meter_counts": event_counts.tolist(),
        "local_circuits": [m[1] for m in meters],
        "duration": 1,
        "conditional_means": mean.tolist(),
        "conditional_variance": variance.tolist(),
        "short_coefficient_path_means": (
            score_sum / total / (9 * board.n * 144)
        ).tolist(),
        "truncation_probability_bound": 1e-12,
    }


def run():
    rng = np.random.default_rng(SEED)
    internal = sorted(
        set(local.survey.routes.cell().seams)
        - set(local.survey.routes.codebook()[0][0])
    )
    reports = []
    for lev in (0, 1, 2):
        b = local.survey.bulk.federation(lev)
        steps, hist, _ = exp.flight.routes.layout(b)
        codes = exp.flight.routes.receiver_codes(b, steps)
        items = arrangements(b, steps, codes, internal)
        graph = nx.Graph((i // 12, j // 12) for i, j in b.external)
        hops = nx.single_source_shortest_path_length(graph, 0)
        if lev == 0:
            exact = (sp.eye(b.cells) - sp.Rational(31, 384) * sp.Matrix(hist)).inv()
            G = np.array([[F(x) for x in row] for row in exact.tolist()], dtype=object)
            eps = F(0)
        else:
            approx, e = exp.flight.routes.green(hist, 5)
            G = np.vectorize(F, otypes=[object])(approx)
            eps = F(float(np.nextafter(e, np.inf)))
        cache = {}
        coeff = []
        for region, partners, P, scores in items:
            expectation = F(0)
            error = F(0)
            for k in range(5):
                for ii in combinations(range(4), k):
                    vertices = tuple(sorted(region[i] for i in ii))
                    if vertices not in cache:
                        cache[vertices] = sign_mean(G, vertices)
                    fourier = (
                        sum(
                            v * (-1) ** sum(bits[i] for i in ii)
                            for bits, v in zip(
                                product((0, 1), repeat=4), scores, strict=True
                            )
                        )
                        / 16
                    )
                    expectation += fourier * cache[vertices]
                    error += (
                        abs(fourier)
                        * (2 * k + 8 * 32 * k * k)
                        * eps
                        / (1 - 2 * k * eps)
                    )
            den = 9 * b.n * 144
            coeff.append((expectation / den, error / den))
        shell = [
            sum(x for v, (x, e) in enumerate(coeff) if hops[v] == i)
            for i in range(max(hops.values()) + 1)
        ]
        # Error in the projected vector divided by t^3 after retaining the
        # cubic term is <= 8^4 t/24 in the space of path-weighted vectors.
        lower = [
            max(0, sqrt(max(0, float(144 * (x - e)))) - 8**4 * float(T) / 24) ** 2 / 144
            for x, e in coeff
        ]
        smallest = min(x - e for x, e in coeff)
        assert smallest > 0
        duration = F(1)
        while (8**4 * duration / 24) ** 2 > 36 * smallest:
            duration /= 2
        row = {
            "cells": b.cells,
            "all_cell_duration": str(duration),
            "all_cell_excess_probability_lower": str(smallest * duration**6 / 4),
            "compact_dimensions": [int(np.trace(z[2])) for z in items],
            "per_cell_t6": [float(x) for x, e in coeff],
            "per_cell_error": [motion_error(x, e) for x, e in coeff],
            "shell_t6": [float(x) for x in shell],
            "positive_at_fixed_duration": sum(v > 0 for v in lower),
            "finite_positive_t6_lower": [
                max(0, v * (1 - 1e-12) - 1e-25) for v in lower
            ],
        }
        if lev == 0:
            row["unit_time_path_witnesses"] = isolated_path_bounds(
                b, items, codes, internal, hops
            )
            row["exact_shell_t6"] = list(map(str, shell))
            row["exact_cell_t6"] = [str(x) for x, e in coeff]
            row["accepted_trajectory_sample"] = propagation(
                b, steps, items, codes, internal, rng
            )
        reports.append(row)
    return {
        "seed": SEED,
        "duration": str(T),
        "surveys": reports,
        "inputs": exp.save_inputs(
            "LOCAL_TRAPPING.md", "routed_contacts.json", "material_polarization.json"
        ),
    }


def motion_error(x, e):
    return steady.motion.shells.output_error(x, e)


if __name__ == "__main__":
    z = run()
    exp.emit("local_trapping", z)
    for r in z["surveys"]:
        print(
            r["cells"],
            min(r["per_cell_t6"]),
            max(r["per_cell_error"]),
            r["positive_at_fixed_duration"],
            r["shell_t6"][:4],
        )
