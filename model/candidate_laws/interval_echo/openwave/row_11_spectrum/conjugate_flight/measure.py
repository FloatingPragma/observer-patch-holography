"""Conjugate material charges on a freely released route across owners."""

import argparse
from collections import Counter
from functools import cache
from itertools import combinations
import json
from math import atan2, pi, sqrt
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta
from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from model.experiments.interval_echo_gauge.flag_measure.measure import meter

HERE = Path(__file__).resolve().parent
LAW = candidate.native_law
SHOTS = 2**23
TIMES = np.linspace(0, 3, 13)


def edge(a, b):
    return (min(a, b), max(a, b))


def carrier_certificate():
    """The exact Gram determinants test injectivity from gateway-only profiles."""
    adjacency = cell().adjacency().toarray().astype(int)
    determinants = []
    for tips in combinations(range(12), 3):
        b = adjacency[np.ix_([s for s in range(12) if s not in tips], list(tips))]
        a, d, e, f, g, h, i, j, k = (b.T @ b).ravel().tolist()
        determinants.append(
            a * (g * k - h * j) - d * (f * k - h * i) + e * (f * j - g * i)
        )
    if min(determinants) <= 0:
        raise AssertionError("Unseen pure gateway channel")
    return dict(
        three_gateway_charts=len(determinants),
        gram_determinant_histogram=dict(sorted(Counter(determinants).items())),
    )


def charge_reflection():
    """Check the full echo, exchange and contrast symmetry via the integer write table."""
    load = [LAW.SRC.index((2 - a, 2 - b)) for a, b in LAW.SRC]
    memory = [LAW.LABELS.index((1 - c, -d)) for c, d in LAW.LABELS]
    perm = [9 * m + p for m in memory for p in load]
    write = LAW._write_pattern()
    np.testing.assert_array_equal(write[np.ix_(perm, load)], write)
    for operator in LAW._label_ops():
        np.testing.assert_array_equal(operator[np.ix_(memory, memory)], operator)
    total = [*load, *(9 + v for v in perm)]
    for term in (*LAW.local_term_parts(), LAW.echo()):
        np.testing.assert_allclose(term[np.ix_(total, total)], term, atol=2e-14)
    return dict(
        load_permutation=load, memory_permutation=memory, native_dimension=len(total)
    )


@cache
def seam_data(k):
    exchange, contrast = LAW.local_term_parts()
    h = exchange + k * contrast
    echo = LAW.echo()
    blank = np.eye(LAW.DIM)[:, :9]
    np.testing.assert_allclose(h @ blank, blank @ h[:9, :9], atol=1e-13)
    swap = np.eye(9)[:, [LAW.SRC.index((b, a)) for a, b in LAW.SRC]]
    np.testing.assert_allclose(expm(-1j * pi * h[:9, :9]), swap, atol=2e-13)
    for value in (0, 1):
        marked = echo[:, LAW.SRC.index((value, value))].reshape(11, 9)
        head = marked[:, LAW.SRC.index((value, value))]
        np.testing.assert_allclose(
            marked, np.outer(head, np.eye(9)[LAW.SRC.index((value, value))])
        )
        code = np.kron(head[:, None], np.eye(9))
        for term in (exchange, contrast):
            np.testing.assert_allclose(term @ code, 0, atol=2e-14)
        assert abs(head @ head - 1) < 1e-13
    private = contrast[:9, :9].diagonal()
    np.testing.assert_allclose(private, [(a - b) ** 2 for a, b in LAW.SRC], atol=1e-13)
    return h[:9, :9], private


def canonical_cycle(cycle):
    first = cycle.index(min(cycle))
    word = cycle[first:] + cycle[:first]
    return min(tuple(word), (word[0], *word[:0:-1]))


def layout():
    board = federation(0)
    graph = nx.Graph(board.seams)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    cycles = []
    for a, b in sorted(owners.edges()):
        opened = owners.copy()
        opened.remove_edge(a, b)
        cycles.extend(canonical_cycle(c) for c in nx.all_shortest_paths(opened, a, b))
    tour = min(cycles, key=lambda c: (len(c), c))
    gateways = {
        (a // 12, b // 12): a for x, y in board.external for a, b in [(x, y), (y, x)]
    }
    ring = []
    for n, home in enumerate(tour):
        left = gateways[home, tour[n - 1]]
        right = gateways[home, tour[(n + 1) % len(tour)]]
        ring.extend(
            min(
                nx.all_shortest_paths(
                    graph.subgraph(range(12 * home, 12 * home + 12)), left, right
                )
            )
        )
    assert len(set(ring)) == len(ring)
    active = {edge(a, b) for a, b in zip(ring, ring[1:] + ring[:1])}
    cuts = {e for e in board.seams if set(e).intersection(ring)} - active
    launch_gate, ref = next(
        (a, b)
        for a in ring
        for b in board.neighbours[a]
        if b not in ring and board.degree[b] == 5
    )
    offset = ring.index(launch_gate)
    ring = ring[offset:] + ring[:offset]
    distances = nx.multi_source_dijkstra_path_length(graph, ring)
    donor = max(
        (p for p in range(board.n) if p != ref), key=lambda p: (distances[p], -p)
    )
    assert distances[donor] > 2
    return board, graph, ring, sorted(cuts), ref, donor


def closed_walk(path):
    bonds = list(zip(path, path[1:]))
    return bonds + bonds[-2::-1] if len(bonds) > 1 else bonds


def prepare_seams(board, graph, ring, cuts, ref, donor, sign, dark_private):
    """Use a physical vacancy for each private write and two filled loads for each cut."""
    word = np.ones(board.n, dtype=int)
    start = ring[1]
    word[start], word[donor] = 0, 2

    def walk(bonds):
        for a, b in bonds:
            assert graph.has_edge(a, b)
            word[a], word[b] = word[b], word[a]

    sealed = []
    if dark_private:
        for p in ring:
            if board.degree[p] == 5:
                path = nx.shortest_path(graph, start, p)
                relay = closed_walk(path)
                walk(relay)
                assert (
                    word[p] == 0
                )  # The reference also stays at zero during its write.
                sealed.append(p)
                walk(reversed(relay))
    if sign == 1:
        walk(closed_walk(nx.shortest_path(graph, start, donor)))
    shutter = edge(ring[0], ref)
    open_graph = graph.copy()
    for cut in cuts:
        if cut == shutter:
            continue
        shift = []
        if start in cut:
            destination = next(p for p in (ring[0], ring[2]) if p not in cut)
            shift.append((start, destination))
        walk(shift)
        assert all(word[p] == 1 for p in cut)
        open_graph.remove_edge(*cut)
        walk(reversed(shift))
    assert word[start] == 1 + sign and word[donor] == 1 - sign
    routes = [nx.shortest_path(open_graph, p, ref) for p in ring]
    for p, path in zip(ring, routes, strict=True):
        labels = list(range(board.n))
        for a, b in closed_walk(path):
            labels[a], labels[b] = labels[b], labels[a]
        expected = list(range(board.n))
        expected[p], expected[ref] = ref, p
        assert labels == expected
    return dict(dark_private=sealed, phase_routes=routes, shutter=list(shutter))


def compile_state(ring, target, k, sign):
    """Physical exchanges on consecutive ports, with phases on the borrowed private reference."""
    n = len(ring)
    h, private = seam_data(k)
    pair = [LAW.SRC.index((1 + sign, 1)), LAW.SRC.index((1, 1 + sign))]
    block = h[np.ix_(pair, pair)]
    clock_rate = -k * (
        private[LAW.SRC.index((1 + sign, 0))] - private[LAW.SRC.index((1, 0))]
    )
    state = np.zeros(n, complex)
    state[1] = 1

    def phase(p, angle):
        dwell = (angle * np.sign(clock_rate) % (2 * pi)) / abs(clock_rate)
        state[p] *= np.exp(1j * clock_rate * dwell)

    # The last programming shutter adjoins port 0. The two-mode input has a node there.
    for a in range(1, n - 1):
        b = a + 1
        angle = atan2(np.linalg.norm(target[b:]), abs(target[a]))
        t = 2 * angle
        phase(b, pi / 2)
        state[[a, b]] = expm(-1j * t * block) @ state[[a, b]]
        phase(a, (k + 0.5) * t)
        phase(b, (k + 0.5) * t - pi / 2)
    for j, z in enumerate(target):
        phase(j, np.angle(z))
    np.testing.assert_allclose(state, target, atol=8e-13)
    assert state[0] == 0  # Both shutter endpoints have load one on every branch.
    return state


def binomial_interval(observed, total, risk):
    lower = beta.ppf(risk / 2, observed, total - observed + 1) if observed else 0
    upper = (
        beta.ppf(1 - risk / 2, observed + 1, total - observed)
        if observed < total
        else 1
    )
    return np.array([lower, upper])


def energy_intervals(final_read, n, parity, contrast, records):
    """Endpoint phase, material timer, and independent bond counts, at familywise 99%."""
    # Eight two-mode banks: two drive momenta, both charges and two contrast laws.
    # Bound the timer, parity, and two bounded Fourier quadratures in each bank.
    risk = 0.01 / 32
    clock_p = binomial_interval(records["clock_crossed"], SHOTS, risk)
    times = 2 * np.arcsin(np.sqrt(clock_p))
    fourier_error = sqrt(4 * np.log(2 / risk) / SHOTS)
    phase = -np.angle(-final_read)
    phase_error = np.arcsin(fourier_error / abs(final_read))
    gaps = [(phase - phase_error) / times[1], (phase + phase_error) / times[0]]
    bond = binomial_interval(int(parity[1]), SHOTS, risk)
    means = contrast + n * bond
    return dict(
        gap=gaps,
        lower_line=[means[0] - gaps[1] / 2, means[1] - gaps[0] / 2],
        upper_line=[means[0] + gaps[0] / 2, means[1] + gaps[1] / 2],
    )


def ring_hamiltonian(board, ring, preparation, k, sign):
    local, _ = seam_data(k)
    pair = [LAW.SRC.index((1 + sign, 1)), LAW.SRC.index((1, 1 + sign))]
    n = len(ring)
    h = local[np.ix_(pair, pair)]
    gen = np.zeros((n, n))
    private = np.zeros(n)
    for a, p in enumerate(ring):
        b = (a + 1) % n
        gen[np.ix_([a, b], [a, b])] += h
        if board.degree[p] == 5 and p not in preparation["dark_private"]:
            private[a] = k * ((1 + sign) ** 2 - 1)
    return gen + np.diag(private), private, h


def sharp_flight(board, ring, prep, k, sign, rng):
    gen, _, _ = ring_hamiltonian(board, ring, prep, k, sign)
    reads = []
    for t in [0.0, 1.0, 2.0, 3.0]:
        state = expm(-1j * t * gen)[:, 1]
        counts = rng.multinomial(SHOTS, abs(state) ** 2 / np.vdot(state, state).real)
        owners = sorted(set(p // 12 for p in ring))
        by_owner = [
            int(sum(c for p, c in zip(ring, counts) if p // 12 == home))
            for home in owners
        ]
        reads.append(
            dict(
                time=t,
                position=counts.tolist(),
                owners=owners,
                by_owner=by_owner,
                outside_launch_owner=SHOTS - by_owner[owners.index(ring[1] // 12)],
            )
        )
    return reads


def acquire(board, ring, preparation, num, k, sign, rng):
    n = len(ring)
    gen, degree_term, local = ring_hamiltonian(board, ring, preparation, k, sign)
    coord = np.exp(2j * pi * num * np.arange(n) / n)
    source = compile_state(ring, (1 - coord) / sqrt(2 * n), k, sign)
    vals, vecs = np.linalg.eigh(gen)
    time_records = []
    readings = []
    clock = []
    for t in TIMES:
        state = vecs @ (np.exp(-1j * t * vals) * (vecs.T @ source))
        counts = rng.multinomial(SHOTS, abs(state) ** 2 / np.vdot(state, state).real)
        nc = int(rng.binomial(SHOTS, abs(expm(-1j * t * local)[1, 0]) ** 2))
        time_records.append(dict(position=counts.tolist(), clock_crossed=nc))
        readings.append(np.vdot(coord, counts / SHOTS))
        clock.append(2 * np.arcsin(sqrt(nc / SHOTS)))
    slope = float(-np.polyfit(clock, np.unwrap(np.angle(-np.array(readings))), 1)[0])
    # One uniformly drawn physical edge per trial, one-occupied flag and native swap parity.
    antisym = float(np.sum(abs(source - np.roll(source, -1)) ** 2) / (2 * n))
    sym = 2 / n - antisym
    parity = rng.multinomial(SHOTS, [sym, antisym, 1 - 2 / n])
    # Every sharp terminal word has exactly two unit-contrast bonds.
    tags = meter()
    words = np.ones((n, n), dtype=int)
    words[np.diag_indices(n)] = 1 + sign
    decoded = np.asarray(tags)[words]
    loads = (decoded - 1).sum(axis=1)
    costs = np.sum((decoded - np.roll(decoded, -1, axis=1)) ** 2, axis=1)
    np.testing.assert_array_equal(loads, sign * np.ones(n, dtype=int))
    np.testing.assert_array_equal(costs, 2 * np.ones(n, dtype=int))
    first_counts = np.array(time_records[0]["position"])
    contrast = float(k * (costs @ first_counts) / SHOTS)
    # A complete position cohort supplies the private-reference contribution if present.
    mean = (
        contrast
        + n * parity[1] / SHOTS
        + degree_term @ (np.array(time_records[0]["position"]) / SHOTS)
    )
    owner = ring[0] // 12
    remote = np.array([p // 12 != owner for p in ring])
    return dict(
        mode=num,
        measured_route_charge=int(loads @ first_counts // SHOTS),
        load_contrast_energy=contrast,
        simultaneous_intervals=energy_intervals(
            readings[-1], n, parity, contrast, time_records[-1]
        )
        if preparation["dark_private"]
        else None,
        clock_and_position=time_records,
        bond_plus_minus_other=parity.tolist(),
        inferred_frequency=slope,
        mean_energy=float(mean),
        inferred_lower_and_upper=[float(mean - slope / 2), float(mean + slope / 2)]
        if preparation["dark_private"]
        else None,
        harmonic_modulus_range=[
            float(min(map(abs, readings))),
            float(max(map(abs, readings))),
        ],
        exterior_owner_change=float(
            (np.array(time_records[-1]["position"]) - time_records[0]["position"])
            @ remote
            / SHOTS
        ),
    )


def measure():
    board, graph, ring, cuts, ref, donor = layout()
    assert meter() == [0, 1, 2]
    cases = []
    for dark, k in [(True, 8), (False, 8), (True, 4)]:
        for sign in (-1, 1):
            prep = prepare_seams(board, graph, ring, cuts, ref, donor, sign, dark)
            rows = [
                acquire(
                    board,
                    ring,
                    prep,
                    m,
                    k,
                    sign,
                    np.random.default_rng(
                        26002 + 1000 * k + 100 * dark + 10 * sign + m
                    ),
                )
                for m in (1, 2)
            ]
            cases.append(
                dict(
                    private_marks=dark,
                    kappa=k,
                    sign=sign,
                    banks=rows,
                    sharp_launch=sharp_flight(
                        board,
                        ring,
                        prep,
                        k,
                        sign,
                        np.random.default_rng(31000 + 1000 * k + 100 * dark + sign),
                    ),
                )
            )
    return dict(
        ring=ring,
        owners=sorted({p // 12 for p in ring}),
        cuts=cuts,
        donor=donor,
        preparation=prepare_seams(board, graph, ring, cuts, ref, donor, -1, True),
        board_certificate=carrier_certificate(),
        native_reflection=charge_reflection(),
        shots_per_bank=SHOTS,
        cases=cases,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = json.dumps(measure(), indent=2) + "\n"
    if parser.parse_args().write:
        (HERE / "results.json").write_text(result)
    else:
        print(result, end="")
