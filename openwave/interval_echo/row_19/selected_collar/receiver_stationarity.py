"""Measure clocked tetrahedra at the exterior receivers and compare equal volumes."""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import asin, atan2, prod, sqrt

from mpmath import iv
import networkx as nx
import numpy as np

from . import routed_contacts as routes

local = routes.local
HERE = local.HERE
SEED = 28122893
T = 1 / 128
INTERROGATIONS = 256


def point(x):
    x = F(x)
    return iv.mpf(x.numerator) / x.denominator


def bounds(x):
    return [
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    ]


def band(ab):
    return iv.mpf(ab)


def entropy(p):
    return -p * iv.log(p) - (1 - p) * iv.log(1 - p)


def clock(graph, port, original_degrees, rng):
    degree = np.array([graph.degree(v) for v in range(12)])
    L = np.diag(degree) - nx.to_numpy_array(graph, nodelist=range(12))
    amplitudes = []
    for q in (1, 2):
        H = L / 2 + np.diag(8 * q * q * (6 - original_degrees + degree))
        energies, states = np.linalg.eigh(H)
        amplitude = (
            np.exp(-1j * T * energies / INTERROGATIONS) @ states[port] ** 2
        ) ** INTERROGATIONS
        amplitudes.append(amplitude)
    a, b = amplitudes
    visibility = a.conjugate() * b
    retained = (abs(a) ** 2 + abs(b) ** 2) / 2
    reads = [
        local.acquisition(
            [(retained + x) / 2, (retained - x) / 2, 1 - retained], [1, -1, 0], rng
        )
        for x in (visibility.real, -visibility.imag)
    ]
    centers = [sum(r["mean"]) / 2 for r in reads]
    radius = max((r["mean"][1] - r["mean"][0]) / 2 for r in reads)
    # The centered single-step return differs from unity by at most
    # t² Var(H)/2. Telescoping both branch products gives this disk radius.
    zeno = degree[port] * T * T / (4 * INTERROGATIONS)
    angle = atan2(centers[1], centers[0])
    error = asin((sqrt(2) * radius + zeno + 1e-12) / np.linalg.norm(centers))
    return {
        "port": int(port),
        "reads": reads,
        "visibility_zeno_error": float(zeno),
        "relative_frequency": [
            (angle - error) / (120 * T),
            (angle + error) / (120 * T),
        ],
    }


def gram_data(t, D):
    a, *other = t
    g = [
        [int(D[a, x] ** 2 + D[a, y] ** 2 - D[x, y] ** 2) for x in other] for y in other
    ]
    x, y, z = g[0]
    _, u, v = g[1]
    w = g[2][2]
    det = x * (u * w - v * v) - y * (y * w - z * v) + z * (y * v - z * u)
    if min(x, x * u - y * y, det) <= 0:
        return None
    herons = []
    for i, j, k in combinations(t, 3):
        c, d, e = int(D[i, j]), int(D[j, k]), int(D[k, i])
        herons.append(prod([c + d + e, -c + d + e, c - d + e, c + d - e]))
    return [det, herons]


def geometry(raw, rho):
    det, herons = raw
    return (
        iv.sqrt(point(F(det, 8))) / 6 * rho**3,
        sum(iv.sqrt(point(h)) / 4 for h in herons) * rho**2,
    )


def compensate(left, right, pf, ps):
    def mix(values, p, i):
        return (1 - p) * values[0][i] + p * values[1][i]

    target = (mix(left, pf, 0) + mix(right, pf, 0)) / 2
    first = mix(left, ps, 0)
    denom = mix(right, ps, 0) - first
    if 0 in denom:
        return "unresolved_or_equal_denominator", None
    alpha = (target - first) / denom
    if alpha.a <= 0 or alpha.b >= 1:
        return "outside_or_touches_public_probability_boundary", None
    # Differentiate at the fresh preparation. The public fair coin is at
    # its entropy maximum, so its tangent entropy contribution vanishes.
    fresh_gap = mix(right, pf, 0) - mix(left, pf, 0)
    if 0 in fresh_gap:
        return "unresolved_or_equal_fresh_volumes", None
    dp = ps - pf
    alpha_dot = (
        -dp * (left[1][0] - left[0][0] + right[1][0] - right[0][0]) / (2 * fresh_gap)
    )
    area_dot = dp * (
        left[1][1] - left[0][1] + right[1][1] - right[0][1]
    ) / 2 + alpha_dot * (mix(right, pf, 1) - mix(left, pf, 1))
    entropy_dot = dp * iv.log((1 - pf) / pf)
    area = (
        (1 - alpha) * mix(left, ps, 1)
        + alpha * mix(right, ps, 1)
        - (mix(left, pf, 1) + mix(right, pf, 1)) / 2
    )
    ds = entropy(ps) - entropy(pf) + entropy(alpha) - iv.log(2)
    row = {
        "second_tetrahedron_probability": bounds(alpha),
        "equal_volume": bounds(target),
        "delta_area": bounds(area),
        "delta_entropy_nats": bounds(ds),
        "area_tangent": bounds(area_dot),
        "entropy_tangent_nats": bounds(entropy_dot),
    }

    def sign(x):
        return "positive" if x.a > 0 else ("negative" if x.b < 0 else "unresolved")

    return sign(area) + "_area_" + sign(ds) + "_entropy", row


def survey(cell, board, steps, extra, odd, rng):
    cell_edges = set(local.survey.routes.cell().seams)
    base = set(local.survey.routes.codebook()[0][0])
    port = int(np.flatnonzero(steps[cell] == cell)[0])
    readings, matrices = [], []
    originals = np.array(board.degree)[12 * cell : 12 * (cell + 1)]
    for marks in (base, base | set(extra)):
        g = nx.Graph()
        g.add_nodes_from(range(12))
        g.add_edges_from(sorted(cell_edges - marks))
        # Each recorded route launches on a reset word. Every blank gate is
        # the calibrated native pi swap, each unavailable gate is a sharp mark.
        D = np.array(
            [[nx.shortest_path_length(g, i, j) for j in range(12)] for i in range(12)]
        )
        matrices.append(D)
        readings.append(clock(g, port, originals, rng))
    ideal = [point(1), point(F(4, 5))]
    acquired = [band(r["relative_frequency"]) for r in readings]
    stock = {}
    for t in combinations(range(12), 4):
        raw = [gram_data(t, d) for d in matrices]
        if all(raw):
            stock[t] = raw
    pf, ps = F(25, 36) * odd, F(101, 144) * odd
    reads = [local.acquisition([float(p), float(1 - p)], [1, 0], rng) for p in (pf, ps)]
    stats = Counter()
    held = {}
    for x, y in combinations(stock, 2):
        if len(set(x) & set(y)) != 3:
            continue

        def shapes(rate):
            return [
                [geometry(g, f) for g, f in zip(stock[t], rate, strict=True)]
                for t in (x, y)
            ]

        key, data = compensate(*shapes(ideal), point(pf), point(ps))
        stats[key] += 1
        if key == "positive_area_positive_entropy" and not held:
            measured, interval = compensate(
                *shapes(acquired), *(band(r["mean"]) for r in reads)
            )
            held = {
                "stations": [list(x), list(y)],
                "ideal": data,
                "measured_class": measured,
                "measured": interval,
                "sharp_gram_and_heron": [stock[x], stock[y]],
            }
    return {
        "cell": cell,
        "earliest_swap_counts": [d.tolist() for d in matrices],
        "sharp_clocks": readings,
        "occupancy_reads": reads,
        "eligible_tetrahedra": [list(t) for t in stock],
        "examined": 495,
        "face_pair_counts": dict(stats),
        "first_positive_pair": held,
    }


def run():
    iv.dps = 50
    rng = np.random.default_rng(SEED)
    b = local.survey.bulk.federation(0)
    steps, _, _ = routes.layout(b)
    second = routes.receiver_codes(b, steps)
    profile = next(
        r
        for r in local.inputs("routed_contacts")["profiles"]
        if (r["stations"], r["stop_bits"]) == (20, 5)
    )
    return {
        "seed": SEED,
        "starts_per_acquisition": local.N,
        "duration": T,
        "interrogations": INTERROGATIONS,
        "surveys": [
            survey(c, b, steps, second[c], F(p), rng)
            for c, p in zip(
                profile["held_edge"], profile["exact_odd_on_held_edge"], strict=True
            )
        ],
        "inputs": {
            f: sha256((HERE / f).read_bytes()).hexdigest()
            for f in (
                "RECEIVER_STATIONARITY.md",
                "routed_contacts.json",
                "routed_contact_flight.json",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "receiver_stationarity.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    for row in z["surveys"]:
        print(row["cell"], row["face_pair_counts"], row["first_positive_pair"])
