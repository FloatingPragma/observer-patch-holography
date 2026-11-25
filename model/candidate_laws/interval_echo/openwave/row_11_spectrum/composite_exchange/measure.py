"""Exchange and free flight of two charged three-port cores on an eight-port ring."""

import argparse
from itertools import combinations, product
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm

from model.executable_model import federation
from model.experiments.interval_echo_gauge.flag_measure.measure import meter
from ..conjugate_flight.measure import LAW, edge, seam_data

N = 8
SHOTS = 2**20
WAITS = (0, 1, 16, 256, 4096, 65536, 1048576)
REGIMES = ("native", "center_fermion", "extra_half_flux", "fermion_and_flux")


def source(copies):
    return (0, 2, 0, 1) + ((0, 2, 0, 1) if copies == 2 else (1, 1, 1, 1))


def cycle(board):
    def extend(path):
        if len(path) == N:
            return path if path[0] in board.neighbours[path[-1]] else None
        for port in sorted(board.neighbours[path[-1]]):
            if port < 12 and port not in path:
                result = extend(path + [int(port)])
                if result:
                    return result
        return None

    return extend([0])


def preparation(board, ring, copies):
    """Native π swaps present equal physical loads to every silent write."""
    word = list(source(copies))
    initial = word.copy()
    cuts = sorted(
        set(e for e in board.seams if set(e).intersection(ring))
        - {edge(ring[j - 1], ring[j]) for j in range(N)}
    )
    # Outside loads at the incident ports remain unit throughout the writes.
    excluded = set(ring).union(*(set(e) for e in cuts))
    stock = [p for p in reversed(range(board.n)) if p not in excluded][:copies]
    assert sum(v - 1 for v in word) + len(stock) == 0
    writes = []
    for kind, ports in [("private", [a]) for a in ring if board.degree[a] == 5] + [
        ("physical", list(e)) for e in cuts
    ]:
        targets = [ring.index(a) for a in ports if a in ring]
        load = 0 if kind == "private" else 1
        swaps = []
        # Stable adjacent permutation on the open ring path realizes any choice,
        # including a chord whose endpoints must both see unit load.
        desired = {}
        donors = iter(i for i, value in enumerate(initial) if value == load)
        for t in targets:
            desired[t] = next(donors)
        unused = iter(i for i in range(N) if i not in desired.values())
        order = [desired[i] if i in desired else next(unused) for i in range(N)]
        labels = list(range(N))
        for i, label in enumerate(order):
            for j in range(labels.index(label), i, -1):
                swaps.append(j - 1)
                labels[j - 1], labels[j] = labels[j], labels[j - 1]
                word[j - 1], word[j] = word[j], word[j - 1]
        assert all(word[t] == load for t in targets)
        if kind == "physical":
            assert all(port not in stock for port in ports)
        writes.append(dict(kind=kind, ports=ports, ring_edge_indices=swaps))
        for j in reversed(swaps):
            word[j], word[j + 1] = word[j + 1], word[j]
        assert word == initial
    return dict(ring=ring, charge_stock=stock, writes=writes)


def basis(copies):
    words = []
    for centers in combinations(range(N), copies):
        for vacancies in combinations(sorted(set(range(N)) - set(centers)), 2 * copies):
            word = [1] * N
            for i in centers:
                word[i] = 2
            for i in vacancies:
                word[i] = 0
            words.append(tuple(word))
    return words


def crossing(a, b, left, right, copies, regime):
    """Exterior sign of an excess crossing the cut past the other excesses."""
    if {a, b} != {0, N - 1} or left == right or 2 not in (left, right):
        return 1
    parity = (copies - 1) * (regime in ("center_fermion", "fermion_and_flux"))
    parity += regime in ("extra_half_flux", "fermion_and_flux")
    return (-1) ** parity


def generator(words, copies, regime, local):
    lookup = {word: i for i, word in enumerate(words)}
    h = np.zeros((len(words), len(words)))
    for col, word in enumerate(words):
        for a in range(N):
            b = (a + 1) % N
            left, right = word[a], word[b]
            p = LAW.SRC.index((left, right))
            h[col, col] += local[p, p]
            if left != right:
                other = list(word)
                other[a], other[b] = right, left
                h[lookup[tuple(other)], col] += local[
                    LAW.SRC.index((right, left)), p
                ] * crossing(a, b, left, right, copies, regime)
    return h


def exchange(copies, regime):
    # Each step advances both intact cores by one port, starting at the leading
    # edge. Every constituent hops into a unit sea site, preserving its order.
    word = list(source(copies))
    tags = [f"{i//4}:{i%4}" if x != 1 else "" for i, x in enumerate(word)]
    phase = 1
    circuit = []
    for step in range(N // copies):
        hops = []
        for c in range(copies):
            for offset in (2, 1, 0):
                a = (4 * c + step + offset) % N
                b = (a + 1) % N
                assert word[a] != 1 and word[b] == 1
                phase *= crossing(a, b, word[a], word[b], copies, regime)
                hops.append(a)
                word[a], word[b] = word[b], word[a]
                tags[a], tags[b] = tags[b], tags[a]
        assert all(
            word[(4 * c + step + 1 + i) % N] == x
            for c in range(copies)
            for i, x in enumerate((0, 2, 0))
        )
        circuit.append(hops)
    assert word == list(source(copies))
    final_tags = [[tags[4 * c + j] for j in range(3)] for c in range(copies)]
    return dict(macrostep_edges=circuit, final_tags=final_tags, phase=int(phase))


def odd_core_control(rng):
    """Physical exchanges compare an odd internal wave on each whole core."""
    left = (2, 0, 0)
    right = (0, 0, 2)
    reflection = [0, 1, 0]
    schedule = sum(exchange(2, "native")["macrostep_edges"], [])

    def apply_swaps(word, edges):
        w = list(word)
        for j in edges:
            k = (j + 1) % N
            w[j], w[k] = w[k], w[j]
        return tuple(w)

    # In each of two addressed blocks a controlled reflection from a sharp
    # left word prepares (left - right)/√2 at its dark output, with p=1/2.
    coeff = {}
    for a, b in product((left, right), repeat=2):
        coeff[a + (1,) + b + (1,)] = (1 if a == b else -1) / 2
    out = {}
    for label, circuit in (
        ("exchange", schedule),
        ("exchange_and_first_reflection", schedule + reflection),
    ):
        overlap = sum(
            v * coeff.get(apply_swaps(w, circuit), 0) for w, v in coeff.items()
        )
        selected = int(rng.binomial(SHOTS, 0.25))
        bright = int(rng.binomial(selected, (1 + overlap) / 2))
        out[label] = dict(
            overlap=overlap,
            bright=bright,
            dark=selected - bright,
            preparation_refused=SHOTS - selected,
        )
    return out


def record(rng, copies, regime, h, words):
    start = words.index(source(copies))
    vals, vec = eigh(h)
    shifted = [source(copies)[-j:] + source(copies)[:-j] for j in range(N)]
    intact = [words.index(w) for w in dict.fromkeys(shifted)]
    out = []
    for t in WAITS:
        psi = vec @ (np.exp(-1j * t * (vals - vals[0])) * vec[start, :])
        probs = abs(psi) ** 2
        counts = rng.multinomial(SHOTS, probs / sum(probs))
        out.append(
            dict(
                time=t,
                port_word_counts=counts.tolist(),
                intact=int(sum(counts[intact])),
                source=int(counts[start]),
                intact_away=int(sum(counts[intact]) - counts[start]),
            )
        )
    return out


def run():
    board = federation(0)
    ring = cycle(board)
    assert len(ring) == N
    local, _ = seam_data(8)
    # Check all blank inputs, including coherences, against actual native time π.
    swap = np.eye(9)[:, [LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    np.testing.assert_allclose(expm(-1j * pi * local), swap, atol=2e-13)
    assert tuple(meter()) == (0, 1, 2)
    rng = np.random.default_rng(20261125)
    runs = {}
    for copies in (1, 2):
        words = basis(copies)
        records = {}
        for regime in REGIMES:
            circuit = exchange(copies, regime)
            s = circuit["phase"]
            circuit["phase_counts"] = {}
            for label, brightness in (
                ("real", (1 + s) / 2),
                ("reversed_reference", (1 - s) / 2),
                ("quadrature", 0.5),
            ):
                bright = int(rng.binomial(SHOTS, brightness))
                circuit["phase_counts"][label] = [bright, SHOTS - bright]
            h = generator(words, copies, regime, local)
            records[regime] = dict(
                circuit=circuit, free=record(rng, copies, regime, h, words)
            )
        runs[str(copies)] = dict(
            preparation=preparation(board, ring, copies),
            port_words=["".join(map(str, w)) for w in words],
            regimes=records,
        )
    return dict(
        shots=SHOTS,
        waits=WAITS,
        ring=ring,
        runs=runs,
        odd_internal_control=odd_core_control(rng),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x") as fp:
        json.dump(run(), fp, separators=(",", ":"))
        fp.write("\n")
