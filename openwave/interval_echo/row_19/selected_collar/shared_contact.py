"""Read the clock and free-probe force of a receiving overlap memory."""

import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import log, pi, sqrt

import numpy as np
from scipy.linalg import eigh, expm

from . import paired_ranging as survey

HERE = survey.HERE
N = 2**50
SEED = 28122885
TIMES = (1 / 64, 1 / 32)
# Simultaneous Hoeffding bound on every reported scalar outcome frequency.
ERROR = sqrt(log(2 * 5000 / 0.01) / (2 * N)) + 1e-12


def inputs(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def acquisition(probs, scores, rng, starts=N):
    p = np.asarray(probs, float)
    assert p.min() >= -1e-12 and abs(p.sum() - 1) < 1e-12
    p = p.clip(0)
    p /= p.sum()
    batch = min(starts, N)
    assert starts % batch == 0
    counts = rng.multinomial(batch, p, size=starts // batch).sum(axis=0)
    # Hoeffding directly on the per-start score, including zero-valued misses.
    scores = np.asarray(scores, float)
    avg = (counts / starts) @ scores
    half = (sqrt(log(2 * 5000 / 0.01) / (2 * starts)) + 1e-12) * np.ptp(scores)
    return {"counts": counts.tolist(), "mean": [float(avg - half), float(avg + half)]}


def difference(new, old):
    return [new[0] - old[1], new[1] - old[0]]


def divide(a, b):
    assert b[0] * b[1] > 0
    values = [x / y for x in a for y in b]
    return [min(values), max(values)]


def serial_reader(board, symbols, stocks, rng):
    import networkx as nx
    import scipy.linalg as la

    code, common = survey.routes.codebook()
    old = inputs("initial_acceleration")
    forbidden = {v // 12 for v in old["launch_sites"]}
    i, j = [c for c in range(board.cells) if c not in forbidden][:2]
    start, second_start = 12 * i + 3, 12 * j + 3
    first, end = 12 * i + 6, 12 * j + 6
    graph = nx.Graph()
    for a, b in sorted(board.seams):
        if a // 12 == b // 12:
            edge = (a % 12, b % 12)
            passes = common[edge] if a // 12 in (i, j) else edge not in code[0]
            if not passes:
                continue
        graph.add_edge(a, b)
    graph.remove_nodes_from([start, end])
    path = min(nx.all_shortest_paths(graph, first, second_start))
    pulses = [
        (start, first),
        *zip(path[:-1], path[1:], strict=True),
        (second_start, end),
    ]
    truth = np.zeros((9, 9), int)
    law = survey.routes.volume.candidate.native_law
    ex, cost = law.local_term_parts()
    U = la.expm(-1j * np.pi * (ex + 8 * cost))
    labels = [
        survey.routes.volume.label(0, 1),
        survey.routes.volume.label(2, 1),
        survey.routes.volume.label(2, -1),
    ]
    blank = np.eye(11)[0]
    localswap = {}
    for letter, mem in enumerate([blank, *labels]):
        # Project on the same memory after the pulse and check there is no leakage.
        columns = np.kron(
            mem[:, None],
            np.eye(9)[
                :, [law.SRC.index((1, 0)), law.SRC.index((0, 1)), law.SRC.index((0, 0))]
            ],
        )
        block = columns.conj().T @ U @ columns
        np.testing.assert_allclose(U @ columns, columns @ block, atol=3e-13)
        localswap[letter] = block[:2, :2] / block[2, 2]
    for w, z in np.ndindex(9, 9):
        tags = {}
        for c, word in ((i, w), (j, z)):
            letters = divmod(word, 3)
            for n, (a, b) in enumerate(code[word + 1]):
                tags[12 * c + a, 12 * c + b] = 1 + (letters[n] if n < 2 else 0)
        packet = np.zeros(board.n, complex)
        packet[start] = 1
        for a, b in pulses:
            block = localswap[tags.get(tuple(sorted((a, b))), 0)]
            packet[[a, b]] = block @ packet[[a, b]]
        chance = abs(packet[end]) ** 2
        truth[w, z] = round(float(chance))
        np.testing.assert_allclose(chance, symbols[w, 1] * symbols[z, 1], atol=3e-12)
    rows = {}
    for name, stock in zip(("opposite", "fresh", "shared"), stocks[:3], strict=True):
        p = F(int(np.sum(stock * truth)), 2592)
        rows[name] = {
            "probability": str(p),
            "read": acquisition([1 - float(p), float(p)], [0, 1], rng),
        }
    return rows, {
        "cells": [i, j],
        "pulses": [list(v) for v in pulses],
        "sharp_arrivals": truth.tolist(),
    }


def source_and_echo(rng):
    board = survey.bulk.federation(0)
    stocks, _ = survey.standards(board)
    local = inputs("pair_clock_information")["local_reader"]
    detector = local["calibrations"][local["primary"]]
    symbols = np.array(detector["sharp_records_in_calibration_order"])
    one = symbols[:, 1]
    odd = one[:, None] != one
    law = survey.routes.volume.candidate.native_law
    echo = law.echo()
    ex, cost = law.local_term_parts()
    h = ex + 8 * cost
    blank = np.eye(law.DIM)[:, 0]
    record = echo @ blank
    np.testing.assert_allclose(echo @ record, blank, atol=2e-13)
    mem = record.reshape(11, 9)[:, 0]
    mark = np.kron(mem[:, None], np.eye(9))
    np.testing.assert_allclose(h @ mark, 0, atol=2e-13)
    J = np.eye(law.DIM)[:, :9]
    np.testing.assert_allclose(h @ J, J @ law.seam_term(8), atol=2e-13)
    # The two returned alternatives (1,2) and (2,1) admit a physical
    # balanced splitter using the same native exchange.
    pair = [law.SRC.index(z) for z in ((1, 2), (2, 1))]
    splitter = expm(-0.5j * pi * h)[np.ix_(pair, pair)]
    np.testing.assert_allclose(abs(splitter) ** 2, 0.5, atol=1e-13)
    relay = {}
    for seed in (0, 1):
        parity = []
        for n in range(3):
            out = np.linalg.matrix_power(echo, n + seed) @ blank
            parity.append(round(float(np.sum(abs(out[9:]) ** 2))))
        seen = np.array(parity)[one[:, None] + one]
        assert np.array_equal(seen, np.logical_xor(odd, bool(seed)))
        rows = {}
        for tag, weights in zip(
            ("opposite", "fresh", "shared"), stocks[:3], strict=True
        ):
            count = int(weights[seen.astype(bool)].sum())
            p = F(count, 2592)
            rows[tag] = {
                "probability": str(p),
                "read": acquisition([1 - float(p), float(p)], [0, 1], rng),
            }
        relay[str(seed)] = rows
    relay["serial"], path = serial_reader(board, symbols, stocks, rng)
    return relay, path


def patch():
    board = survey.bulk.federation(0)
    old = inputs("initial_acceleration")
    s = old["ordinary_site"]
    a, b = old["launch_sites"]
    owners = [a // 12, b // 12]
    ports = [12 * i + j for i in owners for j in range(12)]
    loc = {v: i for i, v in enumerate(ports)}
    code = survey.routes.codebook()[0][0]
    blank = []
    for u, v in board.seams:
        if u not in loc or v not in loc:
            continue
        if u // 12 == v // 12 and (u % 12, v % 12) in code:
            continue
        if u // 12 != v // 12 and {u, v} != {a, b}:
            continue
        blank.append((loc[u], loc[v]))
    s, a, b = map(loc.get, (s, a, b))
    assert tuple(sorted((s, a))) in blank and tuple(sorted((a, b))) in blank
    assert tuple(sorted((s, b))) not in blank
    return ports, blank, np.array(board.degree)[ports], (s, a, b)


def spectrum(ports, active, degrees, sab, m, q, times=TIMES):
    """Assemble the sector by applying every native contact to its two loads."""
    s, a, b = sab
    n = len(ports)
    basis = (
        (list(combinations(range(n), 2)) if m == q else list(permutations(range(n), 2)))
        if m
        else [(v,) for v in range(n)]
    )

    def ordered(t):
        return tuple(sorted(t)) if m == q else tuple(t)

    index = {k: i for i, k in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=np.int32)
    inside = np.array([sum(ports[v] // 12 == ports[a] // 12 for v in k) for k in basis])
    for j, k in enumerate(basis):
        loads = dict(zip(k, (m, q) if m else (q,), strict=True))
        H[j, j] = 16 * sum((6 - degrees[v]) * z * z for v, z in loads.items())
        for u, v in active:
            d = loads.get(u, 0) - loads.get(v, 0)
            if d:
                H[j, j] += 1 + 16 * d * d
                other = ordered([v if x == u else u if x == v else x for x in k])
                H[index[other], j] -= 1
    np.testing.assert_array_equal(H, H.T)
    # All energies are twice their physical value. Removing the uniform
    # unguarded rest beat demodulates both clock branches before readout.
    offset = 80 * (m * m + q * q)
    H[np.diag_indices_from(H)] -= offset
    ids = [index[ordered([s, v] if m else [v])] for v in (a, b)]
    # The first two columns are prepared by the physical half-splitter.
    paths = np.eye(len(H))[:, ids]
    left, right = paths.T
    patterns = np.column_stack(
        [(left + right) / sqrt(2), (left - right) / sqrt(2), left, right]
    )
    energies, U = eigh(H * 0.5)
    coeff = U.T @ patterns
    panels = []
    for t in times:
        wave = U @ ((np.exp(-1j * energies * t))[:, None] * coeff)
        np.testing.assert_allclose(np.sum(abs(wave) ** 2, axis=0), 1, atol=2e-13)
        hist = [
            np.bincount(inside, weights=abs(wave[:, j]) ** 2, minlength=3)
            for j in (0, 1)
        ]
        panels.append(
            (np.array(hist), np.array([wave[v, j + 2] for j, v in enumerate(ids)]))
        )
    h2s = H @ paths
    # Twice the two-site coherence isolates the ab entry of -[H,[H,X]].
    cross = 2 * (h2s[:, 0] * inside) @ h2s[:, 1] - (inside[ids[0]] + inside[ids[1]]) * (
        h2s[:, 0] @ h2s[:, 1]
    )
    v = paths[:, 0] + paths[:, 1]
    Hv = H @ v
    plus_num = int(round((inside * Hv) @ Hv - (inside * v) @ (H @ Hv)))
    return {
        "plus_curvature": F(plus_num, 4),
        "curvature": F(int(cross), 4),
        "short_sines": np.sin(energies / 8192) @ (coeff[:, 2:] ** 2),
        "short_sine_remainder": (abs(energies) ** 3)
        @ (coeff[:, 2:] ** 2)
        / (6 * 8192**2),
        "clock_diagonal": H[ids, ids].copy(),
        "panels": panels,
    }


def timed_reads(bank, rng, starts=N, times=TIMES):
    reads = []
    for ti, t in enumerate(times):
        force = []
        curvature = []
        for q in (1, 2):
            hist = sum(float(f) * pair[q - 1]["panels"][ti][0] for f, pair in bank)
            data = [acquisition(row, [0, 1, 2], rng, starts) for row in hist]
            curvature.append(
                [v / t**2 for v in difference(data[0]["mean"], data[1]["mean"])]
            )
            force.append(data)
        clock = []
        freq = []
        for site in (0, 1):
            probs = np.zeros((2, 3))
            for f, pair in bank:
                z1, z2 = [pair[q - 1]["panels"][ti][1][site] for q in (1, 2)]
                vis = z1.conjugate() * z2
                retained = (abs(z1) ** 2 + abs(z2) ** 2) / 2
                for j, val in enumerate((vis.real, -vis.imag)):
                    probs[j] += float(f) * np.array(
                        [
                            (retained + val) / 2,
                            (retained - val) / 2,
                            1 - retained,
                        ]
                    )
            data = [acquisition(row, [1, -1, 0], rng, starts) for row in probs]
            freq.append([v / t for v in data[1]["mean"]])
            clock.append(data)
        reads.append(
            {
                "position_counts": force,
                "ramsey_counts": clock,
                "statistics": [*curvature, difference(*freq)],
            }
        )
    return reads


def calibrations(by_arm):
    def contrast(tag):
        return [
            F(a) - F(b)
            for a, b in zip(
                by_arm[tag]["derivatives"],
                by_arm["fresh"]["derivatives"],
                strict=True,
            )
        ]

    inc = contrast("shared")
    doses = {}
    for tag in ("ordinary_one", "ordinary_two"):
        ordinary = contrast(tag)
        observed = []
        for a, b, c in zip(
            by_arm["shared"]["reads"],
            by_arm[tag]["reads"],
            by_arm["fresh"]["reads"],
            strict=True,
        ):
            observed.append(
                [
                    divide(difference(va, vc), difference(vb, vc))
                    for va, vb, vc in zip(
                        a["statistics"],
                        b["statistics"],
                        c["statistics"],
                        strict=True,
                    )
                ]
            )
        doses[tag] = {
            "reference_change": list(map(str, ordinary)),
            "required_fractions": [
                str(a / b) for a, b in zip(inc, ordinary, strict=True)
            ],
            "finite_dose_intervals": observed,
        }
    return inc, doses


def rest_lapse(stock, arms, rng):
    both = {}
    for q in (1, 2):
        reads = {}
        for tag in ("fresh", "shared"):
            p = F(arms[tag]["probability"])
            freq = []
            for v in (0, 1):
                weights = (1 - p, p)
                amplitude = sum(
                    float(w) * stock[k, 0, q]["short_sines"][v]
                    for k, w in enumerate(weights)
                )
                err = sum(
                    float(w) * stock[k, 0, q]["short_sine_remainder"][v]
                    for k, w in enumerate(weights)
                )
                # A reference-arm quadrature on all starts reads -Im of
                # the return amplitude. Escapes have zero mean score.
                data = acquisition(
                    [(1 + amplitude) / 2, (1 - amplitude) / 2], [1, -1], rng
                )
                val = [
                    8192 * z + sgn * err
                    for z, sgn in zip(data["mean"], (-1, 1), strict=True)
                ]
                freq.append((data, val))
            # In this cage both gateways have the same base coordination.
            # A written focal contact removes precisely half a kinetic unit.
            field = difference(freq[0][1], freq[1][1])
            occupancy = arms[tag]["read"]["mean"]
            reads[tag] = {
                "phase_counts": [x[0] for x in freq],
                "fractional_rest_difference": [
                    (x + y / 2) / (40 * q * q)
                    for x, y in zip(field, occupancy, strict=True)
                ],
            }
        both[str(q)] = {
            "arms": reads,
            "source_change": difference(
                reads["shared"]["fractional_rest_difference"],
                reads["fresh"]["fractional_rest_difference"],
            ),
        }
    return both


def material(rng, relay):
    ports, edges, degrees, sab = patch()
    stock = {
        (occ, m, q): spectrum(
            ports,
            [e for e in edges if not occ or e != tuple(sorted(sab[:2]))],
            degrees,
            sab,
            m,
            q,
        )
        for occ in (0, 1)
        for m in (0, 1, 2)
        for q in (1, 2)
    }
    exact = {}
    for occ in (0, 1):
        for m in (0, 1, 2):
            first, second = [stock[occ, m, q] for q in (1, 2)]
            clock = second["clock_diagonal"] - first["clock_diagonal"]
            exact[occ, m] = [
                first["curvature"],
                second["curvature"],
                F(int(clock[0] - clock[1]), 2),
            ]
    assert exact[1, 0][0] - exact[0, 0][0] == F(17, 4)
    assert exact[1, 0][1] - exact[0, 0][1] == F(65, 4)
    assert exact[1, 0][2] - exact[0, 0][2] == -24
    out = {}
    for seed, arms in relay.items():
        by_arm = {}
        for tag, writer, m in (
            ("fresh", "fresh", 0),
            ("shared", "shared", 0),
            ("ordinary_one", "fresh", 1),
            ("ordinary_two", "fresh", 2),
        ):
            p = F(arms[writer]["probability"])
            exact_response = [
                (1 - p) * a + p * b
                for a, b in zip(exact[0, m], exact[1, m], strict=True)
            ]
            bank = [
                (f, [stock[k, m, q] for q in (1, 2)]) for k, f in enumerate((1 - p, p))
            ]
            reads = timed_reads(bank, rng)
            by_arm[tag] = {
                "derivatives": list(map(str, exact_response)),
                "reads": reads,
            }

        inc, doses = calibrations(by_arm)
        out[seed] = {
            "rest_lapse": rest_lapse(stock, arms, rng),
            "arms": by_arm,
            "overlap_change": list(map(str, inc)),
            "ordinary_doses": doses,
        }
    return {
        "ports": ports,
        "global_source_and_launch_ports": [ports[v] for v in sab],
        "blank_edges": [list(e) for e in edges],
        "full_board_degree": degrees.tolist(),
        "sharp_derivatives": {
            f"{k}_{m}": list(map(str, row)) for (k, m), row in exact.items()
        },
        "comparisons": out,
    }


def repeat(relay):
    out = []
    for n in (1, 2, 4, 8, 16, 32):
        rows = []
        for m in (0, 1):
            r = {}
            for tag in ("fresh", "shared"):
                p = F(relay["0"][tag]["probability"])
                r[tag] = (1 - (1 - 2 * m) * (1 - 2 * p) ** n) / 2
            rows.append(str(r["shared"] - r["fresh"]))
        out.append({"reports_pairs": n, "occupancy_changes": rows})
    return out


def run():
    rng = np.random.default_rng(SEED)
    relay, path = source_and_echo(rng)
    return {
        "seed": SEED,
        "starts_per_setting": N,
        "times": list(TIMES),
        "one_frequency_error": ERROR,
        "relay": relay,
        "serial_route": path,
        "independent_repetition": repeat(relay),
        "material": material(rng, relay),
        "inputs": {
            k: sha256((HERE / k).read_bytes()).hexdigest()
            for k in (
                "SHARED_CONTACT.md",
                "initial_acceleration.json",
                "pair_clock_information.json",
            )
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    z = run()
    path = HERE / "shared_contact.json"
    if parser.parse_args().write:
        path.write_text(json.dumps(z, indent=1, allow_nan=False) + "\n")
    else:
        from cosmology.rar_deep_regime.lensing_profile import check

        check(json.loads(path.read_bytes()), z)
    for k, v in z["material"]["comparisons"].items():
        print(k, v["overlap_change"], v["ordinary_doses"]["ordinary_one"])
