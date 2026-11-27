"""Retain two-station parity measurements for the complete binary exclusion sectors."""

import argparse
from collections import Counter, deque
import gzip
import hashlib
from itertools import combinations
import json
from math import comb, log, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from ...seam_word_meter.finite_program import machine
from .. import apparatus

HERE = Path(__file__).resolve().parent
SHOTS = 65536
SEED = 1241721
TIMES = (0, 0.25, 2, 8)


def admissible(c, k, r):
    return 0 <= k <= 12 * c and 0 <= r <= min(c, k, 12 * c - k) and (k - r) % 2 == 0


def word_count(c, k):
    return sum(comb(c, r) for r in range(c + 1) if admissible(c, k, r))


def coefficient_columns():
    law = apparatus.candidate.native_law
    ex, contrast = law.local_term_parts()
    h = ex + apparatus.candidate.KAPPA * contrast
    diagonal = apparatus.candidate.KAPPA + 0.5
    for a in (0, 1):
        for b in (0, 1):
            j = law.SRC.index((a, b))
            column = np.zeros(law.DIM)
            if a != b:
                column[j] = diagonal
                column[law.SRC.index((b, a))] = -0.5
            np.testing.assert_allclose(h[:, j], column, atol=3e-13, rtol=0)
    return diagonal


def sector(config, k):
    ports = [12 * owner + j for owner in config["owners"] for j in range(12)]
    positions = {p: i for i, p in enumerate(ports)}
    graph = [[] for _ in ports]
    for a, b in apparatus.federation(0).seams:
        if a in positions and b in positions:
            i, j = positions[a], positions[b]
            graph[i].append(j)
            graph[j].append(i)
    states = list(combinations(range(len(ports)), k))
    ids = {s: i for i, s in enumerate(states)}
    target = []
    pointer = [0]
    for state in states:
        active = set(state)
        for u in state:
            for v in graph[u]:
                if v not in active:
                    target.append(ids[tuple(sorted((active - {u}) | {v}))])
        pointer.append(len(target))
    adjacency = csr_matrix(
        (np.ones(len(target)), target, pointer), shape=(len(states), len(states))
    )
    if (adjacency != adjacency.T).nnz:
        raise ValueError("Unsymmetric exclusion graph")
    h = diags(np.diff(pointer) * coefficient_columns()) - 0.5 * adjacency
    start = ids[tuple(range(k))]
    return ports, states, start, adjacency, h.tocsr()


def shortest_certificate(adj, start):
    """Exact natural-number counts of shortest paths in each leading coefficient."""
    d = [-1] * adj.shape[0]
    ways = [0] * len(d)
    d[start], ways[start] = 0, 1
    queue = deque([start])
    while queue:
        u = queue.popleft()
        for v in adj.indices[adj.indptr[u] : adj.indptr[u + 1]]:
            if d[v] == -1:
                d[v] = d[u] + 1
                queue.append(v)
            if d[v] == d[u] + 1:
                ways[v] += ways[u]
    if not all(ways):
        raise ValueError("A configuration has no nonzero leading coefficient")
    return dict(
        distances=d, shortest_walks=ways, distance_hist=sorted(Counter(d).items())
    )


def read_basis(config, ports, states, heads):
    # The compact binary rail configuration and the two calibrated columns encode
    # both complete ternary read tapes, including a sharp vacant exterior sentinel.
    keys = []
    for station in heads:
        if [line["port"] for line in station["paths"]] != ports + [
            config["boundary"][1]
        ]:
            raise ValueError("Declared head route does not address the enrolled word")
        labels = [apparatus.READ_CODES.index(v) for v in station["sharp"]]
        if labels != [0, 1, 2]:
            raise ValueError("Physical meter does not separate ternary digits")
        bits = []
        for state in states:
            word = 0
            for i in state:
                word ^= labels[1] << (i // 12)
            bits.append(word)
        keys.append(bits)
    if keys[0] != keys[1]:
        raise ValueError("Physical meter words disagree")
    return np.asarray(keys[0], dtype=np.int64)


def material_readback(config, ports, state, heads):
    """Independently operate each native route on one selected sharp outcome."""
    stock = [0] * apparatus.federation(config["level"]).n
    for i in state:
        stock[ports[i]] = 1
    for head in heads:
        zero, two = head["references"]
        stock[zero], stock[two] = 0, 2
    before = stock.copy()
    effects, _ = machine.digit_instruments()
    raw = []
    for head in heads:
        digits = []
        for leg in head["paths"]:
            swaps = apparatus.exchange_word(leg["path"])
            for a, b in swaps:
                stock[a], stock[b] = stock[b], stock[a]
            v = stock[head["head"]]
            digits.append(
                sum(
                    effects[v, stock[r]] << (1 - j)
                    for j, r in enumerate(head["references"])
                )
            )
            for a, b in reversed(swaps):
                stock[a], stock[b] = stock[b], stock[a]
        raw.append(digits)
    word, outside = apparatus.decode(raw, len(config["owners"]))
    if stock != before or outside != 0:
        raise ArithmeticError("Held remote measurement changed the sample or sentinel")
    return dict(state=state, word=word, raw=raw)


def displayed(p):
    return "<1e-12" if p < 1e-12 else format(float(p), ".5g")


def panel(probs, keys, allowed, rng):
    np.testing.assert_allclose(probs.sum(), 1, atol=3e-11, rtol=0)
    probs = probs / probs.sum()
    # Sample only supported states; a rounding remainder on an impossible last cell
    # cannot leak through the multinomial routine's last category.
    positive = np.flatnonzero(probs)
    counts = np.zeros(len(probs), dtype=np.int64)
    counts[positive] = rng.multinomial(SHOTS, probs[positive] / probs[positive].sum())
    word_probs = np.bincount(keys, weights=probs, minlength=64)
    word_counts = np.zeros(64, dtype=np.int64)
    np.add.at(word_counts, keys, counts)
    if word_counts.sum() != SHOTS or any(
        word_counts[i] for i in range(64) if i not in allowed
    ):
        raise ValueError("Forbidden parity or missing readout")
    miss = [i for i in allowed if word_counts[i] == 0]
    seen = [i for i in allowed if word_counts[i]]
    residual = max(abs(word_counts / SHOTS - word_probs))
    # 2 inventories × 2 Hamiltonians × 4 times × 64 bins.
    radius = sqrt(log(2 * 1024 / 0.01) / (2 * SHOTS))
    if residual > radius:
        raise ValueError("Frozen simultaneous concentration check failed")
    return dict(
        counts=counts.tolist(),
        words=word_counts.tolist(),
        probability_per_1e10=np.rint(word_probs * 10**10).astype(np.int64).tolist(),
        observed=len(seen),
        missing_mass=displayed(word_probs[miss].sum()),
        smallest_seen=displayed(min(word_probs[i] for i in seen)),
        smallest_unseen=displayed(min((word_probs[i] for i in miss), default=0))
        if miss
        else None,
        largest_unseen=displayed(max((word_probs[i] for i in miss), default=0))
        if miss
        else None,
        largest_residual=round(float(residual), 10),
    )


def report():
    config = apparatus.layout(0)
    heads = apparatus.calibrate(config)
    enrolled = []
    rng = np.random.Generator(np.random.PCG64(SEED))
    for k in (1, 2):
        ports, states, start, adj, h = sector(config, k)
        keys = read_basis(config, ports, states, heads)
        allowed = sorted(set(map(int, keys)))
        if allowed != [m for m in range(64) if admissible(6, k, m.bit_count())]:
            raise ValueError("Measured register image disagrees with finite inventory")
        initial = np.zeros(len(states), complex)
        initial[start] = 1
        generator = -1j * h
        trajectory = expm_multiply(
            generator,
            initial,
            start=0,
            stop=8,
            num=33,
            traceA=generator.diagonal().sum(),
        )
        # Keep the identical true native diagonal on the shuttered-exchange arm.
        parked = -1j * h.diagonal()
        panels = {}
        for mode in ("native", "diagonal"):
            rows = []
            for t in TIMES:
                state = (
                    trajectory[int(4 * t)]
                    if mode == "native"
                    else np.exp(t * parked) * initial
                )
                rows.append(panel(abs(state) ** 2, keys, allowed, rng))
            panels[mode] = rows
        enrolled.append(
            dict(
                inventory=k,
                ports=ports,
                basis=states,
                start=start,
                words=keys.tolist(),
                allowed=allowed,
                undirected_exchanges=adj.nnz // 2,
                leading=shortest_certificate(adj, start),
                material_readbacks=[
                    material_readback(
                        config, ports, states[int(np.flatnonzero(keys == w)[-1])], heads
                    )
                    for w in allowed
                ],
                panels=panels,
            )
        )
    files = [
        HERE / "PROTOCOL.md",
        Path(__file__),
        Path(apparatus.__file__),
        Path(machine.__file__),
        Path(apparatus.candidate.native_law.__file__),
    ]
    return dict(
        provenance="Simulated finite-shot complete joint sharp rail readings following the unswitched simultaneous native Hamiltonian.",
        config=config,
        stations=heads,
        shots=SHOTS,
        seed=SEED,
        times=TIMES,
        radius=sqrt(log(2 * 1024 / 0.01) / (2 * SHOTS)),
        samples_order="Each joint projection yields the same full binary load configuration at both heads.",
        sources={
            str(p.resolve().relative_to(HERE.parents[4])): hashlib.sha256(
                p.read_bytes()
            ).hexdigest()
            for p in files
        },
        enrolled=enrolled,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--output", type=Path)
    args = cli.parse_args()
    data = (json.dumps(report(), separators=(",", ":")) + "\n").encode()
    if args.output:
        with args.output.open("xb") as f:
            f.write(gzip.compress(data, mtime=0))
    elif gzip.decompress((HERE / "records.json.gz").read_bytes()) != data:
        raise ValueError("Fixed-inventory replay differs")
