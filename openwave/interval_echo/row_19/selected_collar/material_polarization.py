"""Let existing material fall in the overlap field, erase the record, and probe that material."""

from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial, log, sqrt

import networkx as nx
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

from . import contact_census as exp
from . import shell_receivers as shells
from . import spectator_census as spectators

local = exp.local
HERE = exp.HERE
SEED = 28122913
EXPOSURES = [F(1, 8), F(1, 4), F(1, 2), F(1), F(2), F(4)]
FLIGHTS = (1 / 32, 1 / 8)
N = 2**50


def all_boards():
    stock = local.inputs("routed_contacts")["profiles"]
    base = set(local.survey.routes.codebook()[0][0])
    surveys = []
    for level in (0, 1, 2):
        board = local.survey.bulk.federation(level)
        graph = nx.Graph((u // 12, v // 12) for u, v in board.external)
        hop = nx.single_source_shortest_path_length(graph, 0)
        opened = nx.Graph(
            (a, b)
            for a, b in board.seams
            if a // 12 != b // 12 or (a % 12, b % 12) not in base
        )
        for row in (r for r in stock if r["stations"] == board.cells):
            diag = list(map(F, row["diagonal_green"]))
            g = list(map(F, row["green_from_source"]))
            u = [a / (2 * d - 1) for a, d in zip(g, diag, strict=True)]
            err = (
                3
                * F(float(np.nextafter(row["green_rounding_bound"], np.inf)))
                / (2 * min(diag) - 1)
            )
            scale = F(8, 144 * board.n)
            density = [F(0)] * board.n
            for a, b in board.external:
                density[a] = scale * (u[a // 12] - u[b // 12])
                density[b] = -density[a]
            grouped = [sum(density[12 * i : 12 * i + 12]) for i in range(board.cells)]
            assert sum(grouped) == 0
            radial = []
            for k in range(max(hop.values())):
                cut = [
                    (a, b) if hop[a // 12] == k else (b, a)
                    for a, b in board.external
                    if {hop[a // 12], hop[b // 12]} == {k, k + 1}
                ]
                mass = sum(density[a] for a, b in cut)
                exterior = sum(
                    grouped[v] for v in range(board.cells) if hop[v] == k + 1
                )
                assert mass == sum(z for v, z in enumerate(grouped) if hop[v] <= k)
                read = sum(
                    sum(density[j] for j in opened[a] if j != b)
                    - sum(density[j] for j in opened[b] if j != a)
                    for a, b in cut
                ) / len(cut)
                bound = 2 * err * scale
                radial.append(
                    {
                        "hop": k,
                        "edges": len(cut),
                        "enclosed": float(mass),
                        "enclosed_error": shells.output_error(mass, bound * len(cut)),
                        "next_shell_mass": float(exterior),
                        "next_shell_error": shells.output_error(
                            exterior, bound * 3 * sum(h == k + 1 for h in hop.values())
                        ),
                        "secondary_available_gradient": float(read),
                        "secondary_error": shells.output_error(read, bound * 10),
                    }
                )
            surveys.append(
                {
                    "cells": board.cells,
                    "stop_bits": row["stop_bits"],
                    "shells": radial,
                    "density_by_cell": list(map(float, grouped)),
                    "per_cell_error": float(
                        np.nextafter(float(6 * err * scale), np.inf)
                    ),
                }
            )
    return surveys


def exposure_distribution(h, t):
    # Exact integer Taylor numerators on the unnormalised all-ones launch.
    # The integral remainder uses unitarity and the Hermitian infinity norm.
    norm = max(sum(abs(x) for x in row) for row in h)
    size = F(norm) * t / 2
    n = 0
    rem = size
    while rem > F(1, 2**66):
        n += 1
        rem *= size / (n + 1)
    den = (2 * t.denominator) ** n * factorial(n)
    power = np.ones(len(h), dtype=object)
    re = np.zeros_like(power)
    im = np.zeros_like(power)
    for k in range(n + 1):
        term = (
            (-1) ** (k // 2)
            * t.numerator**k
            * (den // ((2 * t.denominator) ** k * factorial(k)))
        )
        if k % 2:
            im -= term * power
        else:
            re += term * power
        if k < n:
            power = h @ power
    probabilities = [
        F(x * x + y * y, len(h) * den**2) for x, y in zip(re, im, strict=True)
    ]
    return np.array([float(v) for v in probabilities]), n, float(2 * rem + rem**2)


def second_instruments(ports, ab, degrees, edges):
    n = len(ports)
    # Outcomes are two or three ordinary returns and a separate blocked launch.
    tables = np.zeros((len(FLIGHTS), n, 8, 4))
    tables[:, list(ab), :, 3] = 1
    returns = {}
    for q in (1, 2):
        basis = list(combinations(range(n), 2) if q == 1 else permutations(range(n), 2))
        index = {v: i for i, v in enumerate(basis)}
        occupied = np.array(
            [sum(ports[x] // 12 == ports[ab[0]] // 12 for x in v) for v in basis]
        )
        row, col, val = [], [], []
        for i, (x, y) in enumerate(basis):
            state = [0] * n
            state[x], state[y] = 1, q
            for result, element in spectators.column(
                tuple(state), edges, degrees
            ).items():
                if q == 1:
                    key = tuple(j for j, z in enumerate(result) if z)
                else:
                    key = (result.index(1), result.index(q))
                row.append(index[key])
                col.append(i)
                val.append(element)
        h = coo_matrix((val, (row, col)), shape=(len(index), len(index))).tocsr()
        assert (h - h.T).nnz == 0
        allowed = [j for j in range(n) if j not in ab]
        ids = np.array(
            [
                [index[tuple(sorted((j, g))) if q == 1 else (j, g)] for g in ab]
                for j in allowed
            ]
        )
        launch = np.eye(len(index))[:, ids.ravel()]
        for it, t in enumerate(FLIGHTS):
            waves = expm_multiply(-0.5j * t * h, launch) * np.exp(40j * t * (1 + q * q))
            for i, j in enumerate(allowed):
                a, b = waves[:, 2 * i : 2 * i + 2].T
                for k, psi in enumerate(((a + b) / sqrt(2), (a - b) / sqrt(2))):
                    tables[it, j, 2 * (q - 1) + k, :3] = np.bincount(
                        occupied, weights=abs(psi) ** 2, minlength=3
                    )
                returns[it, j, q] = [a[ids[i, 0]], b[ids[i, 1]]]
    for it, j in product(range(len(FLIGHTS)), allowed):
        for side, (a, b) in enumerate(
            zip(returns[it, j, 1], returns[it, j, 2], strict=True)
        ):
            vis = a.conjugate() * b
            retained = (abs(a) ** 2 + abs(b) ** 2) / 2
            for k, signal in enumerate((vis.real, -vis.imag)):
                tables[it, j, 4 + 2 * side + k, :3] = [
                    (retained + signal) / 2,
                    (retained - signal) / 2,
                    1 - retained,
                ]
    np.testing.assert_allclose(tables.sum(axis=-1), 1, atol=3e-13)
    return tables


def outcomes(prob, score):
    return float(np.dot(prob, score))


def observe(pos, panels, inside, score, ab, rng, starts=exp.flight.N):
    sample = rng.multinomial(N, pos / pos.sum())

    def band(weights):
        mid = outcomes(sample / N, weights)
        half = np.ptp(weights) * sqrt(log(2 * 10000 / 0.01) / (2 * N)) + 1e-12
        return [mid - half, mid + half]

    reads = []
    for tbl, dt in zip(panels, FLIGHTS, strict=True):
        mixed = np.einsum("j,jmk->mk", pos, tbl)
        counts = [
            local.acquisition(
                vec, [0, 1, 2, 0] if k < 4 else [1, -1, 0, 0], rng, starts
            )
            for k, vec in enumerate(mixed)
        ]

        def diff(left, right, div):
            return [
                v / div
                for v in local.difference(counts[left]["mean"], counts[right]["mean"])
            ]

        reads.append(
            {
                "counts": [x["counts"] for x in counts],
                "scores": [
                    diff(0, 1, dt**2),
                    diff(2, 3, dt**2),
                    diff(5, 7, dt),
                ],
            }
        )
    return {
        "position_counts": sample.tolist(),
        "near": band(inside),
        "available": band(score),
        "blocked": band([int(j in ab) for j in range(len(pos))]),
        "held_flights": reads,
    }


def guarded():
    setup, *_ = exp.prepare()
    profile, _, ports, ab, _, degree, edges, second = setup
    rng = np.random.default_rng(SEED)
    n = len(ports)
    tree = nx.Graph()
    tree.add_nodes_from(range(n))
    tree.add_edges_from(e for e in edges if e not in second[0] and e not in second[1])
    launch = shells.prepare_shell(tree, [(i, i + 1) for i in range(0, n, 2)], 1)
    panels = second_instruments(ports, ab, degree, edges)
    prob = {}
    orders = []
    max_error = 0
    for i, j in product(range(2), repeat=2):
        allowed = [
            e for e in edges if not (i and e in second[0] or j and e in second[1])
        ]
        h = np.diag(16 * (6 - degree) - 80).astype(object)
        for a, b in allowed:
            h[a, a] += 17
            h[b, b] += 17
            h[a, b] -= 1
            h[b, a] -= 1
        for t in EXPOSURES:
            z, norder, error = exposure_distribution(h, t)
            prob[i, j, t] = z
            orders.append(norder)
            max_error = max(error, max_error)
    law = profile["joint_parity_and_stop"]
    accepted = {
        (i, j): float(sum(map(F, law[i][j]))) for i, j in product(range(2), repeat=2)
    }
    score = []
    for s in range(n):
        neighbors = [int(s not in ab and tuple(sorted((s, v))) in edges) for v in ab]
        score.append(neighbors[0] - neighbors[1])
    inside = np.array([float(p // 12 == ports[ab[0]] // 12) for p in ports])
    data = []
    for t in EXPOSURES:
        end = sum(w * prob[i, j, t] for (i, j), w in accepted.items())
        base = prob[0, 0, t]
        arms = {}
        for name, p in (("quiet", 0), ("fresh", 25 / 36), ("shared", 101 / 144)):
            pos = (1 - p) * base + p * end
            arms[name] = observe(pos, panels, inside, score, ab, rng)
        z = arms["shared"]
        ref = arms["fresh"]
        delta = (end - base) / 144
        source = outcomes(delta, score)
        panel = [
            [
                local.difference(x, y)
                for x, y in zip(left["scores"], right["scores"], strict=True)
            ]
            for left, right in zip(z["held_flights"], ref["held_flights"], strict=True)
        ]
        data.append(
            {
                "t": str(t),
                "arms": arms,
                "near_change": local.difference(z["near"], ref["near"]),
                "available_change": local.difference(z["available"], ref["available"]),
                "blocked_change": local.difference(z["blocked"], ref["blocked"]),
                "predicted_near": outcomes(delta, inside),
                "predicted_available": source,
                "predicted_secondary_release": [
                    source * c for c in (8.5, 16.25, -15.5)
                ],
                "finite_secondary": panel,
            }
        )
    return {
        "launch": launch,
        "source_starts": N,
        "secondary_starts": exp.flight.N,
        "secondary_times": list(FLIGHTS),
        "series_max_order": max(orders),
        "series_probability_error": max_error,
        "available_score_by_site": score,
        "ports": ports,
        "data": data,
    }


def run():
    return {
        "seed": SEED,
        "full_boards": all_boards(),
        "guarded": guarded(),
        "inputs": exp.save_inputs(
            "MATERIAL_POLARIZATION.md",
            "routed_contact_flight.json",
            "routed_contacts.json",
        ),
    }


if __name__ == "__main__":
    result = run()
    exp.emit("material_polarization", result)
    for x in result["guarded"]["data"]:
        print(x["t"], x["near_change"], x["available_change"], x["finite_secondary"][1])
