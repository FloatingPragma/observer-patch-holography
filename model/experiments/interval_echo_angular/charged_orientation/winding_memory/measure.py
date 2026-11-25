"""Free charge propagation and coherent winding of native unit-contrast heads."""

import argparse
import json
from fractions import Fraction
from itertools import pairwise
from math import factorial
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh, expm
from scipy.stats import beta

from model.candidate_laws.interval_echo.openwave.row_11_spectrum.conjugate_flight import (
    measure as flight,
)
from ..live_trimer.measure import LAW, blocks


HERE = Path(__file__).resolve().parent
SHOTS = 32768
TIMES = (0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)


def source():
    """Write both contrast signs in phase. Correct either load herald by native flight."""
    ex, con = LAW.local_term_parts()
    native = ex + flight.candidate.KAPPA * con
    input_index = [LAW.SRC.index(w) for w in ((0, 1), (1, 0))]
    block = native[np.ix_(input_index, input_index)]
    split = expm(-0.5j * np.pi * block)[:, 0]
    # A private zero reference gives the vacancy a +kappa phase per unit time
    # relative to the filled load. The phase below acts on the second rail.
    vacancy, host = (LAW.SRC.index((v, 0)) for v in (0, 1))
    rate = flight.candidate.KAPPA * (con[host, host] - con[vacancy, vacancy])
    dwell = 3 * np.pi / (2 * rate)
    split[1] *= np.exp(1j * rate * dwell)
    np.testing.assert_allclose(split[0], split[1], atol=1e-12)
    state = LAW.echo()[:, input_index] @ split
    frame = np.zeros((99, 2))
    for c, address in enumerate(input_index):
        for d in (-1, 1):
            frame[9 + 9 * LAW.LABELS.index((c, d)) + address, c] = 2**-0.5
    amplitudes = frame.T @ state
    np.testing.assert_allclose(frame @ amplitudes, state, atol=2e-12)
    np.testing.assert_allclose(abs(amplitudes) ** 2, [0.5, 0.5], atol=1e-12)
    pulse = expm(-1j * np.pi / flight.candidate.KAPPA * native)
    np.testing.assert_allclose(abs(frame[:, 0] @ pulse @ frame[:, 1]), 1, atol=1e-12)
    return dict(herald_probabilities=[0.5, 0.5], both_heralds_correctable=True)


def local_identity():
    """Check both material currents on every nonzero entry of both full m=1 blocks."""
    raw = blocks()
    triples = [(a, b, c) for c in (0, 1) for a, b in LAW.SRC]
    changing = {}
    for parity in (1, -1):
        h = raw[1, parity]
        n = 0
        for i, j in zip(*np.nonzero(h), strict=True):
            a, b, c = triples[i]
            x, y, z = triples[j]
            assert a - x == c - z == y - b
            n += int(i != j)
        changing[str(parity)] = n
    assert not np.any(raw[0, 1])  # zero-mark shutters for the cut and private contacts
    return dict(native_columns=99, marked_unit_columns=36, directed_moves=changing)


def apparatus():
    board, graph, ring, cuts, _, donor = flight.layout()
    edges = list(pairwise(ring + ring[:1]))
    bypass = graph.copy()
    bypass.remove_edges_from(edges)
    # Borrow a blank degree-five private phase reference outside the ring, keeping
    # the opposite endpoint fixed in the two-rail superposition.
    reference = next(
        p
        for p in sorted(graph)
        if p not in ring and board.degree[p] == 5 and p != donor
    )
    paths = []
    for a, b in edges:
        allowed = bypass.copy()
        allowed.remove_nodes_from((a, donor))
        paths.append(nx.shortest_path(allowed, b, reference))
    staging = next(p for p in board.neighbours[ring[0]] if p not in ring and p != donor)
    stock = bypass.copy()
    stock.remove_node(donor)
    # Stage a single vacancy through blank bypasses while collecting the heads.
    # A completed write, including either herald correction, leaves it at a.
    assert all(nx.has_path(stock, staging, p) for p in ring)
    material = np.ones(board.n, dtype=int)
    material[donor], material[staging] = 2, 0

    def relay(a, b):
        for v, w in pairwise(nx.shortest_path(stock, a, b)):
            material[v], material[w] = material[w], material[v]

    for p in ring:
        if board.degree[p] == 5:
            relay(staging, p)
            assert material[p] == 0  # Write (0,0) on the private zero reference.
            relay(p, staging)
    for a, b in edges:
        relay(staging, a)
        assert (material[a], material[b]) == (0, 1)
        relay(a, staging)
    # Store the vacancy away from ring cuts. Close every shutter except the last,
    # deliver it along the remaining blank path, and flip one native link to park
    # it off the shutter before the final equal-load write.
    shuttle = flight.edge(ring[0], staging)
    storage = next(
        v
        for v in stock
        if v not in ring
        and v not in {donor, staging}
        and all(v not in cut for cut in cuts)
        and nx.has_path(stock.subgraph(set(stock) - set(ring)), staging, v)
    )
    relay(staging, storage)
    assert all((material[a], material[b]) == (1, 1) for a, b in cuts)
    sealed = stock.copy()
    sealed.remove_edges_from(cut for cut in cuts if cut != shuttle)
    path = nx.shortest_path(sealed, storage, ring[0])
    for a, b in pairwise(path):
        material[a], material[b] = material[b], material[a]
    assert (material[ring[0]], material[ring[1]]) == (0, 1)
    material[ring[0]], material[ring[1]] = 1, 0
    assert all(material[p] == 1 for p in shuttle)
    # Reverse the same full native transfer after the shutter write.
    material[ring[0]], material[ring[1]] = 0, 1
    assert material[ring[0]] == 0 and material[donor] == 2
    return ring, dict(
        ring=ring,
        owners=list(dict.fromkeys(p // 12 for p in ring)),
        external_contacts=sum(a // 12 != b // 12 for a, b in edges),
        exterior_zero_writes=len(cuts),
        private_zero_writes=sum(board.degree[p] == 5 for p in ring),
        private_reference=reference,
        stock_staging=staging,
        stock_storage=storage,
        final_shutter=list(shuttle),
        compensating_excess=donor,
        maximum_phase_route=max(len(path) - 1 for path in paths),
    )


def sector(n):
    """Lift every reduced column through the native local matrix on every actual ring link."""
    words = [(x % n, (1 << x) - 1) for x in range(n + 1)]
    index = {w: i for i, w in enumerate(words)}
    local = blocks()[1, 1]
    idle = local[LAW.SRC.index((1, 1)), LAW.SRC.index((1, 1))]
    terms = []
    for a in range(n):
        b = (a + 1) % n
        term = np.zeros((n + 1, n + 1), dtype=int)
        for col, (hole, heads) in enumerate(words):
            coin = (heads >> a) & 1
            load = (int(a != hole), int(b != hole))
            start = 9 * coin + LAW.SRC.index(load)
            for end in np.flatnonzero(local[:, start]):
                pair = LAW.SRC[end % 9]
                new_hole = a if pair[0] == 0 else b if pair[1] == 0 else hole
                new_heads = heads ^ ((coin ^ (end // 9)) << a)
                row = index[
                    new_hole, new_heads
                ]  # any omitted outgoing word is an error
                term[row, col] += local[end, start]
        term -= idle * np.eye(n + 1, dtype=int)
        assert np.array_equal(term, term.T)
        terms.append(term / 8)
    return words, terms


def certify_remote_rank(terms, target):
    """Exact rational enclosure of the foreign two-port minor at t=3/2."""
    scaled = -1.5 * sum(terms)
    np.testing.assert_array_equal(scaled, np.rint(scaled))
    matrix = scaled.astype(int).astype(object)
    norm_bound = int(np.max(np.sum(abs(matrix), axis=1)))
    order = 96
    remainder = Fraction(norm_bound ** (order + 1), factorial(order + 1))
    assert remainder < Fraction(1, 10**16)
    power = np.eye(len(matrix), 2, dtype=int).astype(object)
    real = power.copy()
    imag = 0 * power
    for k in range(1, order + 1):
        power = matrix @ power
        real = k * real + (1, 0, -1, 0)[k % 4] * power
        imag = k * imag + (0, 1, 0, -1)[k % 4] * power
    a, b, c, d = real[target : target + 2].flat
    e, f, g, h = imag[target : target + 2].flat
    real_det = a * d - e * h - b * c + f * g
    imag_det = a * h + e * d - b * g - f * c
    det_squared = Fraction(real_det**2 + imag_det**2, factorial(order) ** 4)
    # Each exact unitary entry is bounded by one. The 2x2 determinant has
    # two products, each displaced by at most 2*remainder + remainder^2.
    det_error = 4 * remainder + 2 * remainder**2
    assert det_squared > (Fraction(1, 10) + det_error) ** 2
    return dict(
        time="3/2",
        order=order,
        operator_norm_upper=norm_bound,
        entry_remainder_upper="1e-16",
        determinant_lower="1/10",
    )


def counts(probability, rng):
    successes = int(rng.binomial(SHOTS, np.clip(probability, 0, 1)))
    # Exact marginal 99% Clopper-Pearson intervals. Counts are simulated launches.
    lo = beta.ppf(0.005, successes, SHOTS - successes + 1) if successes else 0
    hi = beta.ppf(0.995, successes + 1, SHOTS - successes) if successes < SHOTS else 1
    return dict(
        probability=round(float(probability), 12),
        successes=successes,
        launches=SHOTS,
        interval99=[float(lo), float(hi)],
    )


def run():
    preparation = source()
    identity = local_identity()
    ring, route = apparatus()
    n = len(ring)
    words, terms = sector(n)
    assert words[0][0] == words[-1][0]
    winding_support = (words[0][1] ^ words[-1][1]).bit_count()
    assert winding_support == n
    sharp = np.eye(n + 1, dtype=complex)[0]
    k = flight.candidate.KAPPA
    swaps = [expm(-1j * np.pi / k * h) for h in terms]
    assert flight.meter() == [0, 1, 2]
    flight.seam_data(k)  # All blank relays are native full swaps.
    # Two interferometer inputs have identical initial local load and head statistics.
    # Opposite first-splitter quadratures arise from two durations of the same native contact.
    inputs, seeds = [], []
    for duration in (0.5, 1.5):
        state = expm(-1j * np.pi * duration / k * terms[0]) @ sharp
        seeds.append(state.copy())
        for gate in swaps[1:]:
            state = gate @ state
        expected = np.zeros(n + 1)
        expected[[0, n]] = 0.5
        np.testing.assert_allclose(abs(state) ** 2, expected, atol=2e-13)
        inputs.append(state)
    assert abs(inputs[0].conj() @ inputs[1]) < 2e-13

    # Addressed deterministic winding and sense reversal are apparatus controls.
    wound = sharp
    for gate in swaps:
        wound = gate @ wound
    np.testing.assert_allclose(abs(wound), np.eye(n + 1)[n], atol=2e-13)
    unwound = wound
    for gate in reversed(swaps):
        unwound = gate @ unwound
    np.testing.assert_allclose(abs(unwound), abs(sharp), atol=2e-13)
    # At the launch port a single forward contact transfers the unflipped head
    # population to the neighbouring port. The fully wound state stays put.
    source_reader = [float(abs((swaps[0] @ v)[1]) ** 2) for v in (sharp, wound)]
    np.testing.assert_allclose(source_reader, [1, 0], atol=2e-13)
    # Keep a trajectory only if the vacancy starts at the launch and is read at
    # the next prescribed port after each swap. All filters are material reads.
    holes = np.array([hole for hole, _ in words])
    oriented = np.diag((holes == 0).astype(complex))
    for step, gate in enumerate(swaps, 1):
        oriented = (holes == step % n)[:, None] * (gate @ oriented)
    expected = np.zeros_like(oriented)
    expected[-1, 0] = 1
    np.testing.assert_allclose(abs(oriented), expected, atol=2e-13)
    np.testing.assert_allclose(oriented @ oriented, 0, atol=2e-13)
    # Coherent two-input readback reverses the routing, then reverses the first
    # half splitter (its period is 2pi/k after the common idle phase is removed).
    reader = expm(-1.5j * np.pi / k * terms[0])
    for gate in swaps[1:]:
        reader = reader @ gate
    ramsey = [float(abs((reader @ state)[0]) ** 2) for state in inputs]
    np.testing.assert_allclose(ramsey, [1, 0], atol=2e-13)

    eig, vec = eigh(sum(terms))
    rng = np.random.default_rng(551482)
    # Choose a distant material port by graph order, independently of all probabilities.
    target = n // 2
    assert ring[target] // 12 == ring[target + 1] // 12
    assert words[target][1] ^ words[target + 1][1] == 1 << target
    local_read = expm(-0.5j * np.pi / k * terms[target])
    report = []
    for t in TIMES:
        u = (vec * np.exp(-1j * t * eig)) @ vec.T
        probs = [abs(u @ v) ** 2 for v in (sharp, *inputs)]
        mixture = (abs(u[:, 0]) ** 2 + abs(u[:, -1]) ** 2) / 2
        np.testing.assert_allclose((probs[1] + probs[2]) / 2, mixture, atol=3e-13)
        # Charge and prefix-head predictions use the same born weights. Source
        # charge has two indistinguishable positions with different winding head.
        sharp_source = (u @ sharp).copy()
        sharp_source[1:-1] = 0  # Nondemolition source-vacancy herald.
        after_probe = swaps[0] @ sharp_source
        source_and_wound = abs(after_probe[0]) ** 2 + abs(after_probe[-1]) ** 2
        np.testing.assert_allclose(source_and_wound, probs[0][-1], atol=1e-13)
        # Compare a strictly local launch with a relative phase, on two runs.
        # The foreign owner reads signed arrival on its own pair and may mix it
        # before reading. No source or intervening head enters either detector.
        arrivals = []
        for seed in seeds:
            local = u @ seed
            mixed = local_read @ local
            accepted = sum(abs(local[target : target + 2]) ** 2)
            np.testing.assert_allclose(
                accepted, sum(abs(mixed[target : target + 2]) ** 2)
            )
            arrivals.append(
                dict(
                    arrival=counts(accepted, rng),
                    first_port=counts(abs(local[target]) ** 2, rng),
                    first_port_after_mix=counts(abs(mixed[target]) ** 2, rng),
                )
            )
        transfer = u[target : target + 2, :2]
        charge = [p[0] + p[-1] for p in probs]
        report.append(
            dict(
                time=t,
                sharp_at_target=counts(probs[0][target], rng),
                sharp_at_source=counts(charge[0], rng),
                sharp_source_and_wound=counts(source_and_wound, rng),
                local_launches=arrivals,
                remote_pair_determinant=round(float(abs(np.linalg.det(transfer))), 12),
                target_two_inputs=[counts(p[target], rng) for p in probs[1:]],
                target_dephased=counts(mixture[target], rng),
            )
        )
    return dict(
        preparation=preparation,
        seam_current=identity,
        geometry=route,
        remote_rank_certificate=certify_remote_rank(terms, target),
        cyclic_sector_dimension=len(words),
        head_sites_needed_for_stationary_winding_coherence=winding_support,
        target_port=ring[target],
        target_owner=ring[target] // 12,
        source_reader=[round(v, 12) for v in source_reader],
        repeated_directed_winding_probability=round(
            float(np.linalg.norm(oriented @ oriented) ** 2), 12
        ),
        addressed_ramsey=[counts(v, rng) for v in ramsey],
        free_records=report,
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, sort_keys=True, default=int) + "\n"
    if args.write:
        (HERE / "records.json").write_text(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
