"""Compare a clock field carried by recorded walks with ordinary-source transport."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log, sqrt

import networkx as nx
import numpy as np
import sympy as sp

from . import shared_contact as local

HERE = local.HERE
SEED = 28122891
N = 2**17
DEN = 2**40


def layout(board):
    """Derive each transition from uniformly launched ports and the external swaps."""
    steps = np.broadcast_to(np.arange(board.cells)[:, None], (board.cells, 12)).copy()
    ports = {}
    for a, b in board.external:
        for u, v in ((a, b), (b, a)):
            steps[u // 12, u % 12] = v // 12
            ports[u // 12, v // 12] = u
    hist = np.array([np.bincount(row, minlength=board.cells) for row in steps])
    assert np.all(np.diag(hist) == 9) and (hist == 1).sum() == 3 * board.cells
    np.testing.assert_array_equal(hist, hist.T)
    return steps, hist, ports


def receiver_codes(board, steps):
    edges = sorted(local.survey.routes.cell().seams)
    base = set(local.survey.routes.codebook()[0][0])
    candidates = list(
        local.survey.routes.perfect_matchings(
            [e for e in edges if e not in base], set(range(12))
        )
    )
    receivers = []
    for v, row in enumerate(steps):
        gates = set(np.flatnonzero(row != v))
        receivers.append(
            next(
                m
                for m in candidates
                if all(a not in gates or b not in gates for a, b in m)
            )
        )
    return receivers


def green(hist, k):
    q = 2**k
    K = 12 * q * np.eye(len(hist), dtype=np.int64) - (q - 1) * hist
    values = np.rint(np.linalg.solve(K / (12 * q), np.eye(len(K))) * DEN).astype(
        np.int64
    )
    # Rowwise infinity error for each column, from the exact residual
    # and a stochastic Neumann inverse with row-sum norm q.
    residual = 12 * q * DEN * np.eye(len(K), dtype=np.int64) - K @ values
    error = F(int(abs(residual).max()), 12 * DEN)
    return values / DEN, float(error)


def visit_evidence(steps, k, edge, rng):
    n = len(steps)
    odd = np.zeros(n, dtype=int)
    both = np.zeros(n, dtype=int)
    end = odd.copy()
    stopodd = odd.copy()
    two = np.zeros((2, 2, 3), dtype=int)
    starts = 4096
    for _ in range(N // starts):
        sites = np.zeros(starts, dtype=int)
        rows = np.arange(starts)
        parity = np.zeros((starts, n), bool)
        while len(rows):
            parity[rows, sites[rows]] ^= True
            rows = rows[rng.integers(2**k, size=len(rows)) != 0]
            sites[rows] = steps[sites[rows], rng.integers(12, size=len(rows))]
        odd += parity.sum(axis=0)
        both += (parity[::2] & parity[1::2]).sum(axis=0)
        end += np.bincount(sites, minlength=n)
        stopodd += np.bincount(sites[parity[np.arange(starts), sites]], minlength=n)
        stop = np.where(sites == edge[0], 0, np.where(sites == edge[1], 1, 2))
        bins = (2 * parity[:, edge[0]].astype(int) + parity[:, edge[1]]) * 3 + stop
        two += np.bincount(bins, minlength=12).reshape(2, 2, 3)
    return {
        "odd_visits": odd.tolist(),
        "independent_both_odd": both.tolist(),
        "stops": end.tolist(),
        "odd_stops": stopodd.tolist(),
        "two_receivers": two.tolist(),
    }


def endpoint_law(hist, k, edge):
    n = len(hist)
    alpha = sp.Rational(2**k - 1, 2**k)
    signed = {}
    for i, j in np.ndindex(2, 2):
        signs = sp.eye(n)
        for bit, v in zip((i, j), edge, strict=True):
            signs[v, v] = (-1) ** bit
        termination = (
            (1 - alpha)
            * (sp.eye(n) - alpha * signs * sp.Matrix(hist) / 12).inv()
            * signs[:, 0]
        )
        signed[i, j] = [
            termination[edge[0]],
            termination[edge[1]],
            sum(termination) - sum(termination[v] for v in edge),
        ]
    law = [
        [
            [
                str(
                    sum(
                        (-1) ** (i * x + j * y) * signed[i, j][v]
                        for i, j in np.ndindex(2, 2)
                    )
                    / 4
                )
                for v in range(3)
            ]
            for y in range(2)
        ]
        for x in range(2)
    ]
    probs = [F(z) for x in law for y in x for z in y]
    assert min(probs) >= 0 and sum(probs) == 1
    return law


def run():
    rng = np.random.default_rng(SEED)
    data = []
    for level in (0, 1, 2):
        board = local.survey.bulk.federation(level)
        steps, hist, ports = layout(board)
        codes = receiver_codes(board, steps)
        graph, coords = local.survey.bulk.geometry._atlas(
            [local.survey.bulk.federation(i) for i in range(level + 1)]
        )
        graph = graph[level]
        points = coords[level]
        hops = np.array(
            [nx.shortest_path_length(graph, 0, v) for v in range(board.cells)]
        )
        np.testing.assert_array_equal(
            nx.to_numpy_array(graph, nodelist=range(board.cells)),
            hist - 9 * np.eye(board.cells),
        )
        edge = next((u, v) for u, v in sorted(ports) if hops[u] == 2 and hops[v] == 3)
        crossing = [ports[edge], ports[edge[::-1]]]
        park = []
        for g in crossing:
            m = codes[g // 12]
            park.append(
                next(12 * (g // 12) + a + b - g % 12 for a, b in m if g % 12 in (a, b))
            )
        chord = np.linalg.norm(points - points[0], axis=1)
        for k in (3, 5):
            G, eps = green(hist, k)
            D = np.diag(G)
            visits = G[0] / (2 * D - 1)
            alpha = 1 - 1 / 2**k
            stop = (1 - alpha) * G[0]
            stop_odd = (1 - alpha) * D * visits
            P = hist / 12
            expected = (np.eye(len(P)) - alpha * P) @ visits
            row = {
                "stations": len(P),
                "stop_bits": k,
                "green_rounding_bound": eps,
                "diagonal_green": D.tolist(),
                "green_from_source": G[0].tolist(),
                "odd_visit_probabilities": visits.tolist(),
                "odd_stop_probabilities": stop_odd.tolist(),
                "independent_two_reports": (visits**2).tolist(),
                "raw_exterior_equation_max_abs": float(abs(expected[1:]).max()),
                "return_corrected_equation_max_abs": float(
                    abs(
                        ((np.eye(len(P)) - alpha * P) @ ((2 * D - 1) * visits))[1:]
                    ).max()
                ),
                "held_edge": list(edge),
                "probe_ports": crossing,
                "ordinary_partners": park,
            }
            samples = visit_evidence(steps, k, edge, rng)
            error = sqrt(log(2 * 20000 / 0.01) / (2 * N))
            for name, forecast, size in (
                ("odd_visits", visits, N),
                ("stops", stop, N),
                ("odd_stops", stop_odd, N),
                ("independent_both_odd", visits**2, N // 2),
            ):
                assert abs(
                    np.array(samples[name]) / size - forecast
                ).max() < error * sqrt(N / size)
            row["sampled_paths"] = samples
            shells = []
            for h in range(1, int(hops.max()) + 1):
                group = hops == h
                shells.append(
                    {
                        "hops": h,
                        "sites": int(group.sum()),
                        "mean_odd": float(visits[group].mean()),
                        "odd_counts": int(np.array(samples["odd_visits"])[group].sum()),
                        "mean_inverse_chord": float((1 / chord[group]).mean()),
                    }
                )
            x, y = np.array(
                [[z["mean_inverse_chord"], z["mean_odd"]] for z in shells[:3]]
            ).T
            w = (x[1] - x[2]) / (x[0] - x[2])
            row["shells"] = shells
            row["middle_shell_inverse_radius_residual"] = float(
                y[1] - w * y[0] - (1 - w) * y[2]
            )
            gradient = float(visits[edge[0]] - visits[edge[1]])
            # A serial echo at a blank cell lowers the 2-versus-1 rest beat
            # by 24 native units, one fifth of the 120-unit blank beat.
            row["serial_near_minus_far_frequency_per_unit_acceptance"] = -24 * gradient
            row["serial_fraction_of_blank_beat_per_unit_acceptance"] = -gradient / 5
            if level == 0:
                exact = (
                    sp.eye(len(P)) - sp.Rational(2**k - 1, 2**k) * sp.Matrix(hist) / 12
                ).inv()
                assert len(set(exact.diagonal())) == 1
                row["exact_return_green"] = str(exact[0, 0])
                row["exact_odd_on_held_edge"] = [
                    str(exact[0, v] / (2 * exact[v, v] - 1)) for v in edge
                ]
                row["joint_parity_and_stop"] = endpoint_law(hist, k, edge)
            data.append(row)
    return {
        "seed": SEED,
        "starts": N,
        "profiles": data,
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in ("ROUTED_CONTACTS.md", "shared_contact.json")
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    z = run()
    path = HERE / "routed_contacts.json"
    if args.write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    for v in z["profiles"]:
        print(
            v["stations"],
            v["stop_bits"],
            v["serial_near_minus_far_frequency_per_unit_acceptance"],
            v["middle_shell_inverse_radius_residual"],
            v["raw_exterior_equation_max_abs"],
        )
