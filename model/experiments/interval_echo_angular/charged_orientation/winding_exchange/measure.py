"""Exchange and one-body winding with the final physical directions matched."""

import argparse
from collections import Counter
from functools import cache
import json
from math import pi
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm

from ..live_trimer.measure import LAW, blocks
from ..two_winding_loops.measure import stock
from ..winding_memory import measure as winding

N = 15
B = N // 2
SPUR = N
RING = [(i, (i + 1) % N) for i in range(N)]
EDGES = RING + [(B, SPUR)]
SHOTS = 2**20
TIMES = tuple(2.0**j for j in range(-2, 11))


@cache
def apparatus():
    board, g, ring, _, _, _ = winding.flight.layout()
    # Stock deposits use a separate bank of two vacancy-excess pairs. Seal
    # the compensating excesses on different remote dimers.
    branch = next(b for b in sorted(board.neighbours[ring[B]]) if b not in ring)
    ports = ring + [branch]
    route = {tuple(sorted((ports[a], ports[b]))) for a, b in EDGES}
    bypass = g.copy()
    bypass.remove_edges_from(route)
    assert nx.is_connected(bypass)
    ref = next(p for p in sorted(g) if p not in ports and board.degree[p] == 5)
    for a, b in RING:
        test = bypass.copy()
        test.remove_node(ports[a])
        assert nx.has_path(test, ports[b], ref)
    distances = nx.multi_source_dijkstra_path_length(g, ports)
    pairs = []
    for a, b in board.seams:
        if min(distances[a], distances[b]) > 2 and all(
            z not in (a, b)
            and a not in board.neighbours[z]
            and b not in board.neighbours[z]
            for pair in pairs
            for z in pair
        ):
            pairs.append((a, b))
            if len(pairs) == 2:
                break
    assert len(pairs) == 2
    # Store each launch vacancy off its own shutter after writes, verify every
    # other incident seam can close while filled, then feed the two launches
    # back through the separate shutters. No ring edge aids stock writing.
    blocked = set(ports) | {p for pair in pairs for p in pair}
    launch = [ports[0], ports[B]]
    staging = []
    for p in launch:
        q = next(
            q
            for q in sorted(board.neighbours[p])
            if q not in blocked and all(q not in board.neighbours[z] for z in staging)
        )
        staging.append(q)
        blocked.add(q)
    cuts = {
        tuple(sorted((a, b)))
        for a, b in g.edges()
        if set((a, b)).intersection(set(ports) | {z for pair in pairs for z in pair})
        and tuple(sorted((a, b))) not in route | {tuple(sorted(p)) for p in pairs}
    }
    shutters = {tuple(sorted(ab)) for ab in zip(launch, staging)}
    sealed = g.copy()
    sealed.remove_edges_from(cuts - shutters)
    safe = [p for p in sorted(g) if p not in blocked and not any(p in e for e in cuts)]
    outside = sealed.subgraph(set(g) - set(ports))
    storage = []
    for start in staging:
        storage.append(
            next(p for p in safe if p not in storage and nx.has_path(outside, p, start))
        )

    def prepare_loads(count, coherent_tag=False):
        word = np.ones(board.n, int)
        for p, (a, b) in zip(storage[:count], pairs):
            word[p], word[a] = 0, 2
        initial = word.copy()

        def transpose(graph, left, right):
            path = nx.shortest_path(graph, left, right)
            for a, b in winding.flight.closed_walk(path):
                word[a], word[b] = word[b], word[a]

        for p in ports:
            if board.degree[p] == 5:
                transpose(bypass, storage[0], p)
                assert word[p] == 0
                transpose(bypass, p, storage[0])
        for i, j in RING:
            a, b = ports[i], ports[j]
            start = max(a, b) if coherent_tag and i == 1 else a
            transpose(bypass, storage[0], start)
            assert word[start] == 0 and word[a] + word[b] == 1
            # The ordinary even/odd source corrects either herald. On the
            # sharp signed-contrast input retain the requested oriented load
            # herald, leaving a vacancy at a, in both contrast parity sectors.
            word[a], word[b] = 0, 1
            transpose(bypass, a, storage[0])
        assert np.array_equal(word, initial)
        # Every exterior closure uses an actual equal-load write, protecting
        # each distant excess by its reserved native blank dimer.
        for a, b in sorted(cuts - shutters):
            parked = []
            for src, mate in pairs[:count]:
                if src in (a, b):
                    assert mate not in (a, b)
                    word[src], word[mate] = word[mate], word[src]
                    parked.append((src, mate))
            assert word[a] == word[b] == 1
            for src, mate in parked:
                word[src], word[mate] = word[mate], word[src]
        for j, p in enumerate(launch):
            if j < count:
                transpose(outside, storage[j], staging[j])
                q = staging[j]
                word[p], word[q] = word[q], word[p]
                park = ports[((0, B)[j] + 1) % N]
                word[p], word[park] = word[park], word[p]
                assert word[q] == word[p] == 1
                word[p], word[park] = word[park], word[p]
            else:
                assert word[p] == word[staging[j]] == 1
        assert sorted(np.flatnonzero(word == 0)) == sorted(launch[:count])
        assert sorted(np.flatnonzero(word == 2)) == sorted(a for a, b in pairs[:count])

    prepare_loads(1)
    prepare_loads(2)
    prepare_loads(2, coherent_tag=True)
    hard = g.copy()
    hard.remove_edges_from(cuts)
    assert set(nx.node_connected_component(hard, ports[0])) == set(ports)
    assert all(
        set(nx.node_connected_component(hard, a)) == set(pair)
        for pair in pairs
        for a in pair
    )
    stock(False)
    stock(True)
    return dict(
        ports=list(map(int, ports)),
        owners=sorted({p // 12 for p in ports}),
        compensator_dimers=[[int(x) for x in pair] for pair in pairs],
        phase_ref=int(ref),
        private_zero_writes=sum(int(board.degree[p] == 5) for p in ports),
        zero_physical_writes=len(cuts),
        staging=list(map(int, staging)),
        storage=list(map(int, storage)),
    )


@cache
def local(parity):
    raw = blocks()["blank" if parity == 0 else (1, parity)] / 8
    f = LAW.SRC.index((1, 1))
    host = raw[f::9, f::9]
    np.testing.assert_array_equal(host, np.eye(len(host)) * host[0, 0])
    return raw - host[0, 0] * np.eye(len(raw))


def sense(port, other):
    ports = apparatus()["ports"]
    return int(ports[port] > ports[other])


def prepared(n, start=None):
    initial = sum(sense(a, b) << j for j, (a, b) in enumerate(RING))
    holes = (0,) if n == 1 else (0, B)
    return (holes if start is None else tuple(start), initial)


def column(state, odd, parking, twist=False):
    holes, heads = state
    diag = 0.0
    moves = []
    frozen = set(holes)
    for j, (a, b) in enumerate(EDGES if parking else RING):
        hl = int(a not in frozen)
        hr = int(b not in frozen)
        order = sense(a, b)
        pair = (hr, hl) if order else (hl, hr)
        coin = (heads >> j) & 1 if j < N else 0
        native = local(((-1 if j in odd else 1) if j < N else 0))
        i = 9 * coin + LAW.SRC.index(pair)
        diag += native[i, i]
        for o in np.flatnonzero(native[:, i]):
            if o == i:
                continue
            new = LAW.SRC[o % 9]
            assert new == pair[::-1] and hl != hr
            src, dst = (a, b) if hl == 0 else (b, a)
            updated = tuple(sorted((frozen - {src}) | {dst}))
            nexthead = heads ^ ((coin ^ (int(o) // 9)) << j)
            parity = (-1) ** sum(min(a, b) < z < max(a, b) for z in frozen)
            link = -1 if twist and j == N - 1 else 1
            moves.append(((updated, nexthead), native[o, i] * link, parity))
    return diag, moves


def generators(n, odd=frozenset(), parking=True, twist=False):
    states = [prepared(n)]
    lookup = {states[0]: 0}
    entries = []
    diagonal = []
    for k, s in enumerate(states):
        d, m = column(s, odd, parking, twist)
        diagonal.append(d)
        for neighbor, v, swap in m:
            if neighbor not in lookup:
                lookup[neighbor] = len(states)
                states.append(neighbor)
            entries.append((lookup[neighbor], k, v, swap))
    pos = np.diag(diagonal)
    neg = pos.copy()
    for i, k, v, swap in entries:
        pos[i, k] += v
        neg[i, k] += v * swap
    np.testing.assert_array_equal(pos, pos.T)
    np.testing.assert_array_equal(neg, neg.T)
    return states, pos, neg


def programs():
    # Exchange moves the two charges alternately to one another's original
    # addresses. Reference winds A the same distance in total, while B parks.
    arms = [[], [(B, SPUR)]]
    for j in range(B + 1):
        if j < B:
            arms[0].append((j, j + 1))
        arms[0].append(((B + j) % N, (B + j + 1) % N))
    arms[1] += RING + [(SPUR, B)]
    return arms


def forward(state, a, b, time, odd, fermion=False, flux=False):
    holes, heads = state
    assert a in holes and b not in holes
    j = next(i for i, e in enumerate(EDGES) if {a, b} == set(e))
    if j < N:
        parity = -1 if j in odd else 1
        coin = (heads >> j) & 1
        out = 1 - coin
    else:
        coin = out = 0
    i = 9 * coin + LAW.SRC.index((1, 0) if sense(a, b) else (0, 1))
    o = 9 * out + LAW.SRC.index((0, 1) if sense(a, b) else (1, 0))
    u = pulse(j < N, parity if j < N else 0, time)
    amp = u[o, i]
    if fermion and sum(min(a, b) < z < max(a, b) for z in holes) % 2:
        amp = -amp
    if flux and j == N - 1:
        amp = -amp
    return (tuple(sorted((set(holes) - {a}) | {b})), heads ^ ((coin ^ out) << j)), amp


@cache
def pulse(live, parity, time):
    return expm(-1j * time * local(parity if live else 0))


def schedule(program, odd, fermion, flux, flipped=False):
    holes, heads = prepared(2)
    state = (holes, heads ^ ((1 << N) - 1) if flipped else heads)
    tags = {a: k for k, a in enumerate(state[0])}
    amp = 1 + 0j
    for a, b in program:
        state, z = forward(
            state, a, b, pi if SPUR in (a, b) else pi / 8, odd, fermion, flux
        )
        tags[b] = tags.pop(a)
        amp *= z
    return state, [tags[a] for a in state[0]], amp


def pulse_bank(rng):
    arms = programs()
    records = {}
    # One sharp (1,0) input to the native echo prepares a sharp signed
    # contrast, i.e. equal even/odd parity, on either direction herald.
    # Keep the required orientation on arm 1 with a load read at the write.
    p = LAW.SRC.index((1, 0))
    out = LAW.echo()[:, p]
    c = sense(*RING[1])
    addr = LAW.SRC.index((1, 0) if c else (0, 1))
    retained = out.reshape(11, 9)[1:, addr]
    expected = np.zeros(10)
    expected[LAW.LABELS.index((c, 1))] = 2**-0.5
    np.testing.assert_allclose(retained, expected, atol=1e-12)
    prepared_probability = float(retained @ retained)
    # Include a coherent superposition of even and odd contrast on one head.
    # The comparator sign is applied to the sorted-order transport of both
    # identical vacancies; its one-charge hopping and all contact energies agree.
    for model, fermion, flux in [
        ("native", False, False),
        ("fermion", True, False),
        ("half_flux", False, True),
        ("fermion_and_flux", True, True),
    ]:
        v = []
        for odd in (frozenset(), frozenset([1])):
            x = [schedule(prog, odd, fermion, flux) for prog in arms]
            assert x[0][0] == x[1][0] and x[0][1] == [1, 0] and x[1][1] == [0, 1]
            v.append([row[-1] for row in x])
        a, b = np.array(v).T / np.sqrt(2)
        logs = {}
        for name, phase in [("real", 1), ("imag", 1j), ("reversed", -1)]:
            bright = np.linalg.norm(a + phase * b) ** 2 / 4
            dark = np.linalg.norm(a - phase * b) ** 2 / 4
            probs = np.maximum(
                [
                    prepared_probability * bright,
                    prepared_probability * dark,
                    prepared_probability * (1 - bright - dark),
                    1 - prepared_probability,
                ],
                0,
            )
            counts = rng.multinomial(SHOTS, probs / sum(probs))
            logs[name] = counts.tolist()
        records[model] = dict(
            coherent_sector_transfer=[
                [[float(z.real), float(z.imag)] for z in pair] for pair in v
            ],
            banks_bright_dark_transfer_refused_source_refused=logs,
        )
    return records


def winding_qubit(rng):
    arms = programs()
    models = {}
    for model, fermion in (("native", False), ("fermion", True)):
        physical = []
        for arm in arms:
            reverse = [(b, a) for a, b in reversed(arm)]
            op = np.zeros((2, 2), complex)
            for bit in range(2):
                for prog in (arm, reverse):
                    state, _, amp = schedule(
                        prog, frozenset(), fermion, False, bool(bit)
                    )
                    final = int(state[1] != prepared(2)[1])
                    assert final == 1 - bit
                    # Prepare direction |+>, route, keep all forward-current
                    # heralds, recombine and keep its plus output.
                    op[final, bit] += amp / 2
            physical.append(op)
        ident = np.eye(2)
        x = np.array([[0, 1], [1, 0]])
        for op in physical:
            np.testing.assert_allclose(op.conj().T @ op, ident / 4, atol=3e-14)
            # Correct the known winding flip and its phase, then decode a
            # Bell pair between this code and an untouched local reference.
            np.testing.assert_allclose(abs(np.trace(x @ op)), 1.0, atol=3e-14)
        banks = {}
        for name, z in (("real", 1), ("reversed", -1), ("quadrature", 1j)):
            # Use half the identity on the winding qubit, as obtained by
            # tracing the untouched half of a prepared maximally entangled pair.
            plus, minus = ((physical[0] + k * z * physical[1]) / 2 for k in (1, -1))
            bp, dp = (np.linalg.norm(k) ** 2 / 2 for k in (plus, minus))
            probs = np.maximum([bp, dp, 1 - bp - dp], 0)
            banks[name] = rng.multinomial(SHOTS, probs / sum(probs)).tolist()
        odd = np.array([1, -1]) / np.sqrt(2)
        even = np.array([1, 1]) / np.sqrt(2)
        ratios = [(odd @ u @ odd) / (even @ u @ even) for u in physical]
        np.testing.assert_allclose(ratios, [-1, -1], atol=3e-14)
        models[model] = dict(
            transfers=[np.round(u.real, 13).tolist() for u in physical],
            intrinsic_odd_vs_even=[round(float(z.real), 12) for z in ratios],
            banks_bright_dark_refused=banks,
        )
    return models


def head_gauge(states):
    return np.array([(-1) ** ((head >> (N - 1)) & 1) for _, head in states])


def free_bank(n, odd, parking, rng):
    states, h, g = generators(n, odd, parking)
    # Trace every head in port-only returns. Store every port-pair count.
    # A half-flux at the stated ring cut changes the sign of exactly
    # those matrix elements which reverse its explicitly recorded head.
    d = head_gauge(states)
    twisted = generators(n, odd, parking, True)
    assert twisted[0] == states
    for actual, expected in zip(twisted[1:], (h, g)):
        np.testing.assert_array_equal(actual, d[:, None] * expected * d[None, :])
    if not parking and n == 2:
        np.testing.assert_array_equal(g, d[:, None] * h * d[None, :])
    fibers = {}
    for pair, heads in states:
        fibers.setdefault(pair, []).append(heads)
    for hs in fibers.values():
        assert len(hs) == 1 or (len(hs) == 2 and hs[0] ^ hs[1] == (1 << N) - 1)
    bins = sorted(fibers)
    address = {w: i for i, w in enumerate(bins)}
    group = np.array([address[w] for w, head in states])
    results = {}
    for name, m in [("native", h), ("fermion", g)]:
        ev, u = eigh(m)
        records = []
        for t in TIMES:
            p = abs(u @ (np.exp(-1j * (ev - ev[0]) * t) * u[0])) ** 2
            weights = np.bincount(group, weights=p, minlength=len(bins))
            obs = rng.multinomial(SHOTS, weights / sum(weights))
            records.append(
                dict(
                    time=t,
                    counts=obs.tolist(),
                    returns=float(weights[address[states[0][0]]]),
                )
            )
        results[name] = records
    return dict(
        dimension=len(h),
        fixed_address_fibers=dict(Counter(map(len, fibers.values()))),
        ordered_bins=bins,
        records=results,
    )


def run():
    rng = np.random.default_rng(20261126)
    return dict(
        shots=SHOTS,
        material=apparatus(),
        arms=programs(),
        pulse=pulse_bank(rng),
        winding_qubit=winding_qubit(rng),
        free={
            name: free_bank(n, odd, parking, rng)
            for name, n, odd, parking in (
                ("one", 1, frozenset(), True),
                ("two", 2, frozenset(), True),
                ("odd", 2, frozenset([1]), True),
                ("closed_spur", 2, frozenset(), False),
            )
        },
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    name = p.parse_args().output
    data = run()
    with name.open("x") as f:
        json.dump(data, f, separators=(",", ":"))
        f.write("\n")
