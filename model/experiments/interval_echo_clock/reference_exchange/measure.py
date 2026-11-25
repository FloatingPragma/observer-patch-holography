"""Two local reference readouts after complete exchange-word instruments."""

import argparse
import json
from collections import deque
from itertools import pairwise, product
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.stats import beta

from model.candidate_laws.interval_echo import clock
from model.executable_model import federation

from ..modular_reference import reference_reader

HERE = Path(__file__).resolve().parent
SHOTS = 2**19
TOL = 0.012
BASIS = list(product((0, 2), repeat=4))
SECTOR = [i for i, z in enumerate(BASIS) if sum(z[:2]) == sum(z[2:]) == 2]
PI = np.pi


def path(graph, a, b, forbidden):
    prev = {a: None}
    queue = deque([a])
    while b not in prev:
        if not queue:
            raise ValueError("no enrolled route")
        v = queue.popleft()
        for u in sorted(graph[v]):
            if u not in prev and not forbidden(v, u):
                prev[u] = v
                queue.append(u)
    q = [b]
    while q[-1] != a:
        q.append(prev[q[-1]])
    return list(reversed(q))


def contacts(car, cells):
    """Restrict the contact words that bring two remote rails together."""
    links = {frozenset(e) for e in pairwise(cells)}

    def forbidden(a, b):
        c, d = a // 12, b // 12
        return (
            c not in cells
            or d not in cells
            or (c != d and frozenset((c, d)) not in links)
        )

    routes = []
    for offset in (6, 7):
        a, b = (12 * c + offset for c in (cells[0], cells[-1]))
        route = path(car.neighbours, a, b, forbidden)
        # Bring a to b's neighbour, apply the selected exposure, and undo
        # the approach. This conjugates a single contact to the desired pair.
        ports = list(range(car.n))
        for i, j in pairwise(route[:-1]):
            ports[i], ports[j] = ports[j], ports[i]
        assert (ports[route[-2]], ports[route[-1]]) == (a, b)
        routes.append(route)
    return routes


def local_exposure(duration):
    law = clock.law
    h0, h1 = law.local_term_parts()
    h = h0 + clock.KAPPA * h1
    ids = [law.SRC.index(z) for z in product((0, 2), repeat=2)]
    basis = np.eye(law.DIM)[:, ids]
    reduced = basis.T @ h @ basis
    np.testing.assert_allclose(h @ basis, basis @ reduced, atol=1e-13)
    return expm(-1j * duration * reduced)


def gate(duration):
    u = local_exposure(duration)
    out = np.eye(len(BASIS), dtype=complex)
    for a, b in ((0, 2), (1, 3)):
        single = np.zeros_like(out)
        for j, z in enumerate(BASIS):
            for k, (aa, bb) in enumerate(product((0, 2), repeat=2)):
                s = list(z)
                s[a], s[b] = aa, bb
                single[BASIS.index(tuple(s)), j] = u[k, 2 * (z[a] // 2) + z[b] // 2]
        out = single @ out
    return out


def certify_route(words, duration, target):
    """Apply the actual port-level transpositions to all occupancy columns."""
    ends = [q[0] for q in words] + [q[-1] for q in words]
    ports = [frozenset(v for v, b in zip(ends, z) if b) for z in BASIS]
    index = {sites: i for i, sites in enumerate(ports)}
    routed = np.zeros_like(target)
    contact = local_exposure(duration)
    swap = local_exposure(PI)
    np.testing.assert_allclose(swap, np.eye(4)[:, [0, 2, 1, 3]], atol=1e-12, rtol=0)

    def permute(sites, a, b):
        if (a in sites) == (b in sites):
            return sites
        return sites ^ {a, b}

    for column, sites in enumerate(ports):
        state = {sites: 1.0}
        for word in words:
            hops = list(pairwise(word[:-1]))
            for a, b in hops:
                state = {permute(k, a, b): v for k, v in state.items()}
            a, b = word[-2:]
            evolved = {}
            for k, v in state.items():
                initial = 2 * (a in k) + (b in k)
                for final, amp in enumerate(contact[:, initial]):
                    if abs(amp) < 1e-14:
                        continue
                    new = (
                        (k - {a, b})
                        | ({a} if final & 2 else set())
                        | ({b} if final & 1 else set())
                    )
                    evolved[new] = evolved.get(new, 0) + amp * v
            state = evolved
            for a, b in reversed(hops):
                state = {permute(k, a, b): v for k, v in state.items()}
        for k, v in state.items():
            routed[index[k], column] = v
    np.testing.assert_allclose(routed, target, atol=1e-12, rtol=0)


def sources(car):
    launches, settings, pulse = reference_reader.instruments(car)
    refs = []
    for ready, state in launches:
        phases = [np.linalg.matrix_power(pulse, n) @ state for n in range(32)]
        refs.append((ready, sum(np.outer(s[:2], s[:2].conj()) for s in phases) / 32))
    # reference_reader orders one-double states [20,02]; BASIS reverses it.
    effect, before = settings[2]
    effect = (before.conj().T @ effect @ before)[:2, :2][::-1, ::-1]
    return [(r, x[::-1, ::-1]) for r, x in refs], effect


def terminals(rho, effect):
    one = rho[np.ix_(SECTOR, SECTOR)]
    effects = [np.eye(2) - effect, effect]
    hits = [np.trace(np.kron(a, b) @ one).real for a, b in product(effects, repeat=2)]
    bunched = [
        rho[BASIS.index(z), BASIS.index(z)].real for z in ((2, 2, 0, 0), (0, 0, 2, 2))
    ]
    law = np.array(hits + bunched)
    np.testing.assert_allclose(law.sum(), 1, atol=2e-11)
    if min(law) < -1e-13:
        raise ValueError("negative terminal probability")
    return law


def bank(rng, refs, effect, full, half, coins, second=1, coherent=False):
    ready, right = refs[second]
    source = np.zeros_like(full)
    source[np.ix_(SECTOR, SECTOR)] = np.kron(refs[0][1], right)
    mass = [1 - ready]
    words = list(product((0, 1), repeat=len(coins)))
    for word in words:
        weight = np.prod([p if bit else 1 - p for bit, p in zip(word, coins)])
        u = half if coherent else full if sum(word) % 2 else np.eye(16)
        mass.extend(weight * ready * terminals(u @ source @ u.conj().T, effect))
    # The finite categorical law is rounded to 13 decimal places before its
    # fixed-seed multinomial draw, including all complementary outcomes.
    mass = np.round(mass, 13)
    np.testing.assert_allclose(mass.sum(), 1, atol=2e-11, rtol=0)
    counts = rng.multinomial(SHOTS, mass / mass.sum())
    return {
        "failed": int(counts[0]),
        "by_exchange_word": [
            {"word": "".join(map(str, w)), "outcomes": ns.tolist()}
            for w, ns in zip(words, counts[1:].reshape(-1, 6), strict=True)
        ],
    }


def decode(b):
    counts = np.array([row["outcomes"] for row in b["by_exchange_word"]]).sum(axis=0)
    assert counts.sum() + b["failed"] == SHOTS
    p = counts[:4] / counts[:4].sum()
    return np.array([p[2] + p[3], p[1] + p[3]]), float(p[3]), counts[-2:]


def probe_intervals(banks, hits):
    """Simultaneous 99% coverage for two calibration and two probe counts."""

    def slot(hits, total):
        tail = 0.01 / 8
        return np.array(
            [
                beta.ppf(tail, hits, total - hits + 1) if hits else 0.0,
                beta.isf(tail, hits + 1, total - hits) if hits < total else 1.0,
            ]
        )

    def calibrated_gap(bank, side):
        counts = np.sum([r["outcomes"][:4] for r in bank["by_exchange_word"]], axis=0)
        return 2 * np.arctanh(4 * (slot(counts[side] + counts[3], sum(counts)) - 3 / 8))

    home = calibrated_gap(banks["direct_fair"], 2)
    other = calibrated_gap(banks["held_settled"], 1)
    lo, hi = PI / 2 * other / home[::-1]
    if not 0 < lo < hi < PI:
        raise ValueError("calibration does not resolve an unwrapped first half-turn")
    first = (0.5 + 0.5 * np.cos([hi, lo])).tolist()
    endpoints = 0.5 + 0.5 * np.cos(np.array([lo, hi]) + PI / 2)
    shifted = [
        0.0 if lo <= PI / 2 <= hi else float(min(endpoints)),
        float(max(endpoints)),
    ]
    measured = [slot(k, SHOTS) for k in hits]
    for (a, b), (c, d) in zip([first, shifted], measured, strict=True):
        if max(a, c) > min(b, d):
            raise ValueError(
                "calibration and independent probe count intervals are disjoint"
            )
    return [first, shifted]


def run():
    car = federation(0)
    graph = [set() for _ in range(car.cells)]
    for a, b in car.external:
        graph[a // 12].add(b // 12)
        graph[b // 12].add(a // 12)
    a, b = 0, min(graph[0])
    detour = path(graph, a, b, lambda x, y: {x, y} == {a, b})
    assert len(detour) > 2
    routes = {
        name: {"cells": cs, "port_words": contacts(car, cs)}
        for name, cs in [("direct", [a, b]), ("held", detour)]
    }
    # The station on each end uses the same native reader on its own rails.
    for cell in (a, b):
        for j, k in ((6, 7), (6, 8), (7, 8), (8, 9), (6, 10)):
            assert 12 * cell + k in car.neighbours[12 * cell + j]
    refs, effect = sources(car)
    full, half = gate(PI), gate(PI / 2)
    for route in routes.values():
        for duration, target in [(PI, full), (PI / 2, half)]:
            certify_route(route["port_words"], duration, target)
    _, rates, _ = clock._primitives()
    np.testing.assert_allclose(list(rates.values()), 0.5, atol=1e-12)
    p = float(np.round(rates[0], 13))
    rng = np.random.default_rng(754307)
    banks = {}
    for name, coins, second, coherent in [
        ("direct_zero", [], 1, False),
        ("direct_quarter", [p * p], 1, False),
        ("direct_fair", [p], 1, False),
        ("held_four", [p * p] * 4, 1, False),
        ("held_settled", [p * p] * 4 + [p], 1, False),
        ("held_two_sharp", [p], 0, False),
        ("held_coherent", [], 1, True),
    ]:
        banks[name] = bank(rng, refs, effect, full, half, coins, second, coherent)
    means = {k: decode(v)[0] for k, v in banks.items()}
    first = means["direct_zero"]
    contraction = np.diff(means["direct_quarter"])[0] / np.diff(first)[0]
    forecast = np.mean(first) + (first - np.mean(first)) * contraction**4
    for prediction, observed in [
        (forecast, means["held_four"]),
        (means["direct_fair"], means["held_settled"]),
        (np.repeat(first[0], 2), means["held_two_sharp"]),
        (np.repeat(first[0], 2), means["held_coherent"]),
    ]:
        if max(abs(observed - prediction)) > TOL:
            raise ValueError("held reference read differs from the calibration")
    means_set = means["held_settled"]
    if (
        abs(np.diff(means_set)[0]) > 0.008
        or max(abs(means_set - means["held_two_sharp"])) < 0.05
    ):
        raise ValueError("settled contrast or input dependence unresolved")
    for parity in (0, 1):
        words = banks["held_settled"]["by_exchange_word"]
        filtered = np.sum(
            [row["outcomes"] for row in words if row["word"].count("1") % 2 == parity],
            axis=0,
        )
        p4 = filtered[:4] / filtered[:4].sum()
        resolved = np.array([p4[2] + p4[3], p4[1] + p4[3]])
        if max(abs(resolved - (first[::-1] if parity else first))) > TOL:
            raise ValueError(
                "retaining swap parity failed to resolve the initial references"
            )
    _, joint, _ = decode(banks["held_settled"])
    if joint > np.prod(means_set) - 0.003:
        raise ValueError("joint click does not distinguish private independent mixing")
    pile = decode(banks["held_coherent"])[2]
    if not all(abs(pile / (SHOTS - banks["held_coherent"]["failed"]) - 0.125) < TOL):
        raise ValueError("coherent contact did not resolve both bunching alternatives")
    # Directly acquire terminal signs on fresh sharp hands at two settings.
    # Each exposure here is one accepted tick, independent of both reference banks.
    phi = PI / 16
    probe = clock._dimer_pulse(phi)
    state = np.linalg.matrix_power(probe, 8) @ np.array([1.0, 0.0])
    hits = [
        int(
            rng.binomial(
                SHOTS, np.round(abs((clock._dimer_pulse(phase) @ state)[0]) ** 2, 13)
            )
        )
        for phase in (0.0, PI / 2)
    ]

    def g(prob):
        return 2 * np.arctanh(4 * (prob - 3 / 8))

    s = 8 * phi / g(means["direct_fair"][0])
    pred = [(1 + np.cos(s * g(means_set[1]) + shift)) / 2 for shift in (0.0, PI / 2)]
    # Propagate both reference count errors; each p near 0.46 has SE 0.0008.
    if max(abs(np.array(hits) / SHOTS - pred)) > 0.07:
        raise ValueError("held modular phase misses the independent sharp probe")
    probe_intervals(banks, hits)
    return {
        "board": car.cells,
        "shots": SHOTS,
        "routes": routes,
        "banks": banks,
        "terminal_order": [
            "00",
            "01",
            "10",
            "11",
            "both doubles left",
            "both doubles right",
        ],
        "eight_tick_probe": hits,
        "held_four_forecast": np.round(forecast, 9).tolist(),
        "held_probe_forecast": np.round(pred, 9).tolist(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, help="write a fresh receipt")
    args = parser.parse_args()
    result = run()
    if args.record:
        with args.record.open("x") as f:
            json.dump(result, f, indent=2, allow_nan=False)
            f.write("\n")
    elif result != json.loads((HERE / "records.json").read_text()):
        raise SystemExit("reference-exchange receipt differs")
    for name, bank_ in result["banks"].items():
        print(name, bank_["failed"], *decode(bank_))
    print(
        "detour forecasts and probe counts",
        result["held_four_forecast"],
        result["held_probe_forecast"],
        result["eight_tick_probe"],
    )
