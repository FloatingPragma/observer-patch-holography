"""Release a charged three-defect dressing in the registered blank material."""

import argparse
import json
from collections import Counter
from itertools import combinations, product
from math import factorial
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.linalg import expm
from scipy.sparse.csgraph import shortest_path
from scipy.sparse.linalg import expm_multiply
from scipy.optimize import brentq

from model.candidate_laws.interval_echo.openwave.row_04_compact import experiment
from model.candidate_laws.interval_echo.openwave.row_04_compact import (
    material_reference,
)
from model.executable_model import cell, federation

K = material_reference.K
LAW = material_reference.LAW
SHOTS = 2**20
HERE = Path(__file__).resolve().parent


def transitions(board, state):
    """Destinations of native swaps. Two indistinguishable zeros are stored in order."""
    a, b, x = state
    occupied = {a: 0, b: 0, x: 2}
    dest = []
    for i in state:
        for j in board.neighbours[i]:
            if j in occupied and (i > j or occupied[i] == occupied[j]):
                continue

            def move(p):
                return j if p == i else i if p == j else p

            u, v = sorted((move(a), move(b)))
            dest.append((u, v, move(x)))
    return dest


def potential(board, state):
    a, b, x = state
    holes = int(b in board.neighbours[a])
    contacts = int(a in board.neighbours[x]) + int(b in board.neighbours[x])
    # Host-subtracted sum of squared seam contrasts and degree references.
    return (
        18
        + 2 * (board.degree[a] + board.degree[b] - board.degree[x] - 6)
        + 2 * (contacts - holes)
    )


def restriction(board, states, strength):
    ids = {s: i for i, s in enumerate(states)}
    row, col, diag = [], [], []
    for i, state in enumerate(states):
        dest = transitions(board, state)
        diag.append(len(dest) / 2 + strength * potential(board, state))
        for s in dest:
            if s in ids:
                row.append(i)
                col.append(ids[s])
    matrix = sparse.csr_matrix(
        (np.full(len(row), -0.5), (row, col)), shape=(len(ids), len(ids))
    )
    matrix += sparse.diags(diag, format="csr")
    assert (matrix != matrix.T).nnz == 0
    return matrix


def row_basis(matrix, prime):
    a = matrix.copy() % prime
    pivot = 0
    for col in range(a.shape[1]):
        ids = np.flatnonzero(a[pivot:, col])
        if not len(ids):
            continue
        at = pivot + ids[0]
        a[[pivot, at]] = a[[at, pivot]]
        a[pivot] = (a[pivot] * pow(int(a[pivot, col]), -1, prime)) % prime
        tail = a[pivot + 1 :, col:]
        tail -= a[pivot + 1 :, col, None] * a[pivot, None, col:]
        tail %= prime
        pivot += 1
        if pivot == len(a):
            break
    return a[:pivot]


def compact_rank(board, gateways):
    # On configurations with all three defects off gateways, the interior matrix
    # coincides with that of a degree-five cell, including the degree reference.
    states = [
        (a, b, x)
        for a, b in combinations(range(12), 2)
        for x in range(12)
        if x not in (a, b)
    ]
    h = (2 * restriction(board, states, K)).astype(np.int64)
    interior = np.array([not set(gateways).intersection(s) for s in states])
    a = h[interior][:, interior]
    b = h[~interior][:, interior]
    current = b.toarray()
    current = current[np.any(current != 0, axis=1)]
    prime = 1009
    basis = row_basis(current, prime)
    ranks = [len(basis)]
    while len(basis) < a.shape[0]:
        next_rows = (a.T @ basis.T).T
        basis = row_basis(np.vstack((basis, next_rows)), prime)
        ranks.append(len(basis))
        if ranks[-1] == ranks[-2]:
            break
    return {"interior_states": int(a.shape[0]), "modulus": prime, "ranks": ranks}


def contact_check(strength):
    # In the compressed separated dynamics each of the two species can hop
    # into the other's site. Only their two diagonal half-rates survive compression.
    ids = [LAW.SRC.index(w) for w in ((0, 2), (2, 0))]
    seam = material_reference.anchor.free_seam(strength)
    h = seam[np.ix_(ids, ids)]
    isolated = seam[LAW.SRC.index((0, 1)), LAW.SRC.index((0, 1))]
    isolated += seam[LAW.SRC.index((1, 2)), LAW.SRC.index((1, 2))]
    block = h - isolated * np.eye(2)
    expected = (2 * strength - 0.5) * np.eye(2) - 0.5 * np.eye(2)[::-1]
    np.testing.assert_allclose(block, expected, atol=2e-13)
    return np.linalg.eigvalsh(block).tolist()


def layout(board):
    face = next(
        tuple(f)
        for f in board.faces()
        if max(f) < 12 and all(board.degree[p] == 5 for p in f)
    )
    graph = nx.Graph(board.seams)
    distance = shortest_path(
        board.adjacency(), directed=False, indices=np.arange(12), unweighted=True
    ).min(axis=0)
    targets = {
        d: min(
            np.flatnonzero(
                (board.degree == d) & (distance == max(distance[board.degree == d]))
            )
        )
        for d in (5, 6)
    }
    a, b, c = face
    graph.remove_nodes_from((a, b))
    path = nx.shortest_path(graph, c, int(targets[5]))
    return face, targets, distance, path


def ball(board, starts, radius):
    states = set(starts)
    shell = set(starts)
    for _ in range(radius):
        shell = {u for s in shell for u in transitions(board, s)} - states
        states |= shell
    return sorted(states)


def census():
    charts = Counter()
    bounds = []
    for level in range(3):
        board = federation(level)
        charts.update(tuple(c) for c in experiment.gateway_table(board))
        face, targets, distance, _ = layout(board)
        state = (*face[:2], int(targets[6]))
        value = potential(board, state) * K + len(transitions(board, state)) / 2
        bounds.append(
            {
                "cells": board.cells,
                "separated_trial": state,
                "mean_native_energy": float(value),
                "storage_distance_from_home": int(distance[state[2]]),
            }
        )
    return {
        "layouts": [
            {"gateways": c, "occurrences": num, **compact_rank(cell(), c)}
            for c, num in sorted(charts.items())
        ],
        "separated_comparisons": bounds,
    }


def measure(board, strength, rng):
    face, targets, _, path = layout(board)
    cluster = tuple(face)
    separate = (*face[:2], int(targets[5]))
    radius = 3
    states = ball(board, (cluster, separate), radius)
    h = restriction(board, states, strength)
    ports = [states.index(cluster), states.index(separate)]
    # A half contact splits the excess between the two endpoints. A chain of
    # full pi contacts sends its second component to storage and is undone
    # after free evolution. The two vacancies remain in place in either arm.
    split, combine = material_reference.optics(1)
    np.testing.assert_allclose(abs(split) ** 2, [0.5, 0.5], atol=2e-13)
    # For -8 the same binary half contact and full swaps have identical phases.
    ref_seam = material_reference.anchor.free_seam(strength)
    ids = [LAW.SRC.index(w) for w in ((2, 1), (1, 2))]
    np.testing.assert_allclose(
        expm(-1j * np.pi / 2 * ref_seam)[np.ix_(ids, ids)][:, 0], split, atol=2e-13
    )
    full_swap = expm(-1j * np.pi * ref_seam)
    permutation = [LAW.SRC.index((v, u)) for u, v in LAW.SRC]
    np.testing.assert_allclose(full_swap, np.eye(LAW.N_SRC)[permutation], atol=5e-14)
    rack = np.ones(board.n, dtype=int)
    rack[list(face[:2])] = 0
    rack[path[1]] = 2
    for p, q in zip(path[1:-1], path[2:], strict=True):
        assert q in board.neighbours[p]
        rack[p], rack[q] = rack[q], rack[p]
    assert (
        np.array_equal(np.flatnonzero(rack == 0), face[:2]) and rack[separate[2]] == 2
    )
    # At a degree-five contact the 2 and 1 arms differ in native reference
    # cost by 3*kappa. Every analyzer phase has a nonnegative physical wait.
    cost = 3 * strength * (6 - board.degree[face[2]])
    for reference_phase in (1, -1, 1j, -1j):
        duration = (
            (-np.angle(reference_phase) * np.sign(cost)) % (2 * np.pi) / abs(cost)
        )
        np.testing.assert_allclose(
            np.exp(-1j * duration * cost), reference_phase, atol=1e-13
        )
    inputs = np.zeros(len(states), complex)
    inputs[ports] = split
    deficit = float(h[ports[0], ports[0]] - h[ports[1], ports[1]])
    reduced = h - sparse.eye(len(states), format="csr") * h[ports[0], ports[0]]
    rows = []
    # The terminal success includes both vacancy-address checks. On the one-
    # excess-per-arm subspace their terminal signatures distinguish both words.
    for wait in (0, 1 / 128, 1 / 64, 1 / 32):
        sample = expm_multiply(-1j * wait * reduced, inputs)[ports]
        quadratures = []
        contrasts = []
        for phase in (1, -1, 1j, -1j):
            # A physical degree-five reference gives a relative quadrature;
            # fixed degree charge phase and the reverse path are part of the analyzer.
            out = combine @ (sample * np.array([phase, 1]))
            p = np.r_[abs(out) ** 2, max(0, 1 - np.linalg.norm(out) ** 2)].round(12)
            counts = rng.multinomial(SHOTS, p / p.sum())
            contrasts.append((counts[0] - counts[1]) / SHOTS)
            quadratures.append(counts.tolist())
        fringe = complex(contrasts[0] - contrasts[1], contrasts[3] - contrasts[2]) / 2
        rows.append(
            {
                "time": wait,
                "counts_by_phase": quadratures,
                "complex_fringe": [fringe.real, fringe.imag],
                "return_amplitude_error_bound": 2
                * (9 * wait) ** (2 * radius + 2)
                / float(factorial(2 * radius + 2)),
            }
        )
    phases = np.unwrap(np.angle([complex(*r["complex_fringe"]) for r in rows]))
    inferred = (
        -(phases[1:] - phases[0]) / np.array([r["time"] for r in rows[1:]])
    ).tolist()
    return {
        "sharp_triangle": cluster,
        "separated": separate,
        "route": path,
        "state_count": len(states),
        "full_board_diagonal_difference": deficit,
        "acquisitions": rows,
        "slope_estimates": inferred,
    }


def top_addresses(board):
    return {
        (a, b, x)
        for x in range(board.n)
        if board.degree[x] == 5
        for a, b in combinations(
            sorted(i for i in board.neighbours[x] if board.degree[i] == 6), 2
        )
        if b not in board.neighbours[a]
    }


def high_band(board):
    # Include even the formal degree/adjacency combinations not realized on a board.
    classes = []
    for a, b, x in product((5, 6), repeat=3):
        for h, u, v in product((0, 1), repeat=3):
            energy = K * (18 + 2 * (a + b - x - 6) + 2 * (u + v - h))
            unequal = a + b + x - 2 * h - u - v
            classes.append((energy, unequal))
    top = max(e for e, _ in classes)
    low_cap = max(e + n for e, n in classes if e < top)
    high = top_addresses(board)
    outside_counts = Counter()
    split = []
    moves = Counter()
    for s in high:
        assert potential(board, s) * K == top and len({p // 12 for p in s}) == 1
        dest = transitions(board, s)
        inside = set(dest) & high
        outside_counts.update(set(dest) - high)
        split.append((len(dest), len(inside)))
        for q in inside:
            moves["within_home" if q[2] // 12 == s[2] // 12 else "crosses_cell"] += 1
    centers = [top + n / 2 for n, _ in split]
    if max(centers) != min(centers):
        raise ValueError("Highest configurations do not share their sharp energy")
    rest = centers[0]
    low = min(top + (n - i) / 2 for n, i in split)
    up = max(top + (n + i) / 2 for n, i in split)
    cross = np.sqrt(max(n - i for n, i in split) * max(outside_counts.values())) / 2

    # Rayleigh quotient of the two-by-two matrix of block upper bounds.
    def bound(p):
        return up * p + low_cap * (1 - p) + 2 * cross * np.sqrt(p * (1 - p))

    floor = brentq(lambda p: bound(p) - rest, 0.0, 1.0)
    ceiling = np.linalg.eigvalsh([[up, cross], [cross, low_cap]])[-1]
    return {
        "board_cells": board.cells,
        "top_addresses": len(high),
        "top_potential": float(top),
        "sharp_energy": float(rest),
        "low_block_upper_bound": float(low_cap),
        "high_block_enclosure": [float(low), float(up)],
        "cross_norm_upper": float(cross),
        "whole_sector_upper_bound": float(ceiling),
        "intact_probability_all_exposures": float(floor),
        "highest_band_dimension": len(high),
        "upper_band_target_lower": float(np.sqrt(2) * low),
        "maximizing_address_hops": dict(moves),
    }


def top_flight(board, strength, rng):
    high = top_addresses(board)
    launch = min(high)
    states = ball(board, [launch], 8)
    inner = set(ball(board, [launch], 7))
    boundary = np.array([s not in inner for s in states])
    h = restriction(board, states, strength) - 18 * strength * sparse.eye(
        len(states), format="csr"
    )
    initial = np.zeros(len(states))
    initial[states.index(launch)] = 1
    step = 1 / 16
    timed = expm_multiply(
        -1j * h, initial, start=0.0, stop=1.0, num=33, traceA=-1j * h.diagonal().sum()
    )
    bins = []
    for s in states:
        a, b, x = s
        edges = (
            int(a in board.neighbours[b])
            + int(x in board.neighbours[a])
            + int(x in board.neighbours[b])
        )
        if s == launch:
            idx = 0
        elif s in high:
            idx = 1 if x // 12 == launch[2] // 12 else 2
        elif edges >= 2:
            idx = 3
        else:
            idx = 4
        bins.append(idx)
    # The full and restricted kinetic off-diagonal blocks have norm <= 9.
    # The radius-eight exterior couples only through the last shell. Integrate
    # its actual amplitude using local Hermitian Taylor bounds at 16 midpoints.
    radius = float(np.asarray(abs(h).sum(axis=1)).max())
    error = 0.0
    reads = []
    for n in range(1, 17):
        term = timed[2 * n - 1]
        integral = np.linalg.norm(term[boundary])
        for j in range(1, 19):
            term = (step / (2 * j)) * (h @ term)
            integral += np.linalg.norm(term[boundary]) / (j + 1)
        tail = (step * radius / 2) ** 19 / factorial(19)
        error += 9 * step * (integral + tail)
        if n in (4, 16):
            probability = np.bincount(bins, weights=abs(timed[2 * n]) ** 2, minlength=5)
            np.testing.assert_allclose(sum(probability), 1, atol=2e-12)
            counts = rng.multinomial(
                SHOTS, probability.round(12) / probability.round(12).sum()
            )
            reads.append(
                {
                    "time": float(n * step),
                    "counts": counts.tolist(),
                    "probabilities": probability.tolist(),
                    "full_board_bin_error_upper": float(2 * error),
                }
            )
    return {
        "launch": launch,
        "computed_addresses": len(states),
        "sector_addresses": int(board.n * (board.n - 1) * (board.n - 2) // 2),
        "detector_bins": [
            "original top address",
            "other top address in home",
            "top address in other cell",
            "other connected trimer",
            "disconnected",
        ],
        "exposures": reads,
    }


def run():
    rng = np.random.default_rng(178051)
    board = federation(0)
    return {
        "upper_band": [high_band(federation(i)) for i in range(3)],
        "high_release": top_flight(board, K, np.random.default_rng(178052)),
        "exchange_only_release": top_flight(board, 0, np.random.default_rng(178053)),
        "seed": 178051,
        "release_seeds": [178052, 178053],
        "shots_per_setting": SHOTS,
        "census": census(),
        "strengths": {
            str(k): {
                "unlike_contact_eigenvalues": contact_check(k),
                "transfer": measure(board, k, rng),
            }
            for k in (K, -K)
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run()
    if args.output:
        with args.output.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    else:
        from ..electron_field.replay import compare

        compare(
            json.loads(json.dumps(result)),
            json.loads((HERE / "results.json").read_text()),
            "dressed_charge",
        )
        print("Reproduced the three-defect contact, confinement and transfer tests.")
