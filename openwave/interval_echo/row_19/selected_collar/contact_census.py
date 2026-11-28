"""Separate written contacts from ordinary material around an open crossing."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log, sqrt

import numpy as np
from scipy.linalg import expm

from . import routed_contact_flight as flight

local = flight.local
HERE = local.HERE
N = 2**48
SEED = 28122911
TIMES = (*local.TIMES, 1 / 8, 1 / 4)
COEFFICIENTS = [
    [F(17, 4), F(17, 2), F(65, 4)],
    [F(65, 4), F(65, 4), F(65, 2)],
    [F(-24), F(-31, 2), F(-65, 2)],
    [F(1), F(1), F(1)],
]


def save_inputs(*files):
    return {f: sha256((HERE / f).read_bytes()).hexdigest() for f in files}


def native_pulse():
    law = local.survey.routes.volume.candidate.native_law
    exchange, cost = law.local_term_parts()
    unitary = expm(-1j * np.pi * (exchange + 8 * cost))
    blank = np.eye(11)[:, 0]
    mark = (law.echo() @ np.kron(blank, np.eye(9)[:, 0])).reshape(11, 9)[:, 0]
    truth = []
    for n, mem in enumerate((blank, mark)):
        table = []
        for q in range(3):
            rows = []
            for m in range(3):
                index = law.SRC.index((q, m))
                key = (q, m) if n else (m, q)
                target = np.kron(mem, np.eye(9)[:, law.SRC.index(key)])
                out = unitary @ np.kron(mem, np.eye(9)[:, index])
                overlap = np.vdot(target, out)
                np.testing.assert_allclose(out, overlap * target, atol=2e-13)
                assert abs(abs(overlap) - 1) < 2e-13
                rows.append(key[0])
            table.append(rows)
        truth.append(table)
    return truth


def prepare():
    profile, cells, ports, ab, sites, degrees, edges, second = flight.instrument()
    relay = local.inputs("shared_contact")["relay"]["serial"]
    pf, ps = [F(relay[v]["probability"]) for v in ("fresh", "shared")]
    u, v = map(F, profile["exact_odd_on_held_edge"])
    d = (ps - pf) * (u - v)
    open_sum = 2 - pf * (u + v)
    outward = 31 * d / (3 * open_sum)
    inward = 17 * d / (3 * open_sum)
    parity = profile["joint_parity_and_stop"]
    arms = {}
    for name, p, stock in (
        ("record_fresh", pf, [(F(1), 0, 0)]),
        ("record_shared", ps, [(F(1), 0, 0)]),
        (
            "ordinary_before",
            pf,
            [(1 - outward - inward, 0, 0), (outward, 1, 0), (inward, 2, 1)],
        ),
        (
            "ordinary_after",
            pf,
            [(1 - outward - inward, 0, 0), (outward, 1, 1), (inward, 2, 0)],
        ),
    ):
        weights = defaultdict(F)
        for i, j in np.ndindex(2, 2):
            path = sum(map(F, parity[i][j]))
            for passed, prob in enumerate((1 - p, p)):
                for fraction, m, k in stock:
                    weights[(i * passed, j * passed, m, k)] += path * prob * fraction
        assert sum(weights.values()) == 1 and min(weights.values()) >= 0
        arms[name] = weights
    return (
        (profile, cells, ports, ab, sites, degrees, edges, second),
        arms,
        d,
        outward,
        inward,
    )


def local_census(setup, key, truth):
    _, _, _, ab, sites, degrees, edges, second = setup
    first, last, m, k = key
    active = [
        e
        for e in edges
        if not any(z and e in w for z, w in zip((first, last), second, strict=True))
    ]
    freq = np.zeros((2, 2, 3), int)
    for i, g in enumerate(ab):
        near = [a + b - g for a, b in active if g in (a, b)]
        for q in (1, 2):
            freq[i, q - 1, truth[1][q][0]] += int(degrees[g]) - len(near)
            for j in near:
                freq[i, q - 1, truth[0][q][m if j == sites[k] else 0]] += 1
    return active, freq


def read_census(probs, rng):
    raw = np.array([[rng.multinomial(N, vector) for vector in site] for site in probs])
    # One bound for every load histogram, both probes, sides and preparations.
    margin = sqrt(log(2 * 1000 / 0.01) / (2 * N)) + 1e-14
    grads = 6 * (raw[0] / N - raw[1] / N)
    lo, hi = grads - 12 * margin, grads + 12 * margin
    coords = np.array(
        [
            [lo[0, 1] - hi[1, 1], hi[0, 1] - lo[1, 1]],
            [lo[1, 1], hi[1, 1]],
            [lo[0, 2], hi[0, 2]],
        ]
    )
    vacancy = [-hi[0, 0], -lo[0, 0]]
    return {
        "counts": raw.tolist(),
        "source_intervals": coords.tolist(),
        "vacancy_interval": vacancy,
    }


def linear(matrix, intervals):
    rows = []
    for vector in matrix:
        z = [
            sorted(float(c) * v for v in b)
            for c, b in zip(vector, intervals, strict=True)
        ]
        rows.append([sum(v[0] for v in z) - 1e-12, sum(v[1] for v in z) + 1e-12])
    return rows


def run():
    rng = np.random.default_rng(SEED)
    setup, arms, d, outward, inward = prepare()
    _, _, ports, ab, sites, degrees, _, _ = setup
    truth = native_pulse()
    cache = {}
    readings = {}
    for name, weights in arms.items():
        freqs = np.zeros((2, 2, 3))
        bank = []
        slopes = [F(0)] * 3
        for key, weight in weights.items():
            if not weight:
                continue
            allowed, counts = local_census(setup, key, truth)
            freqs += float(weight) * counts / 6
            _, _, mass, k = key
            if key not in cache:
                cache[key] = [
                    local.spectrum(
                        ports, allowed, degrees, (sites[k], *ab), mass, q, TIMES
                    )
                    for q in (1, 2)
                ]
            one, two = cache[key]
            bank.append((weight, [one, two]))
            diag = two["clock_diagonal"] - one["clock_diagonal"]
            for i, v in enumerate(
                [one["curvature"], two["curvature"], F(int(diag[0] - diag[1]), 2)]
            ):
                slopes[i] += weight * v
        readings[name] = {
            "census": read_census(freqs, rng),
            "exact_release": list(map(str, slopes)),
            "held": local.timed_reads(bank, rng, flight.N, TIMES),
        }
    contrasts = {}
    for name, (new, old) in {
        "record": ("record_shared", "record_fresh"),
        "redistribution": ("ordinary_after", "ordinary_before"),
    }.items():
        left, right = readings[new], readings[old]
        source = [
            local.difference(x, y)
            for x, y in zip(
                left["census"]["source_intervals"],
                right["census"]["source_intervals"],
                strict=True,
            )
        ]
        release = [
            F(x) - F(y)
            for x, y in zip(left["exact_release"], right["exact_release"], strict=True)
        ]
        held = [
            [
                local.difference(x, y)
                for x, y in zip(a["statistics"], b["statistics"], strict=True)
            ]
            for a, b in zip(left["held"], right["held"], strict=True)
        ]
        contrasts[name] = {
            "source_intervals": source,
            "predicted_release": linear(COEFFICIENTS[:3], source),
            "measured_vacancy": local.difference(
                left["census"]["vacancy_interval"], right["census"]["vacancy_interval"]
            ),
            "exact_release": list(map(str, release)),
            "finite": held,
        }
    return {
        "seed": SEED,
        "times": list(TIMES),
        "census_starts": N,
        "flight_starts": flight.N,
        "truth": truth,
        "ordinary_transfer_near_to_far_one": str(outward),
        "ordinary_transfer_far_to_near_two": str(inward),
        "extra_mark_gradient": str(d),
        "arms": readings,
        "contrasts": contrasts,
        "inputs": save_inputs(
            "CONTACT_CENSUS.md",
            "shared_contact.json",
            "routed_contacts.json",
            "routed_contact_flight.json",
        ),
    }


def emit(module, result):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    path = HERE / (module + ".json")
    if parser.parse_args().write:
        path.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), result)


if __name__ == "__main__":
    result = run()
    emit("contact_census", result)
    for name, row in result["contrasts"].items():
        print(
            name, row["exact_release"], row["source_intervals"], row["measured_vacancy"]
        )
