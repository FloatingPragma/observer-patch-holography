"""Measure the same arriving winding channel with separately prepared probe charges."""

import argparse
import json
from functools import cache
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply
from scipy.stats import beta

from ..apparatus import LAW, blocks, fixture, native_graph, winding
from ..measure import BASE, ONE, bank, build, pt
from ..reciprocal import channel

HERE = Path(__file__).resolve().parent
FLIGHT = 32
DWELL = np.pi / 32
SHOTS = 2**39
DIGITS = list(product(range(3), repeat=3))


def physical_probe():
    """Find a separate waiting dimer and an empty route to a blank detector spoke."""
    layout = fixture()
    full = nx.Graph(winding.flight.federation(0).seams)
    quiet = nx.Graph(layout["blank"])
    hub = layout["tip"]
    quiet.remove_node(hub)
    occupied = set(layout["ports"] + layout["donor"])
    options = []
    for near in sorted(leg[-1] for leg in layout["legs"]):
        if not full.has_edge(hub, near):
            raise ValueError("The probe needs an actual blank spoke.")
        for inside in quiet:
            for x in set(full[inside]) - occupied:
                for y in set(full[x]) - occupied:
                    path = nx.shortest_path(quiet, near, inside)
                    options.append((len(path), near, tuple(path), int(x), int(y)))
    _, near, path, x, y = min(options)
    # The storage compensation occupies a separate, shielded dimer. Both
    # stock bodies can be parked on their partner during equal-load writes.
    protected = occupied | {x, y}
    distance = nx.multi_source_dijkstra_path_length(full, protected)
    stock = min(
        tuple(sorted((a, b)))
        for a, b in full.edges
        if min(distance[a], distance[b]) > 1
    )
    added = {x, y, *stock}
    dimer_edges = {frozenset((x, y)), frozenset(stock)}
    cuts = {
        frozenset((a, b))
        for a, b in full.edges
        if (a in added or b in added) and frozenset((a, b)) not in dimer_edges
    }
    # With sharp parked loads, ensure a partner accessible by an active
    # dimer bond is available on every equal-load exterior write.
    for edge in cuts:
        for site, spare in ((y, x), tuple(stock)):
            if site in edge:
                assert spare not in edge and full.has_edge(site, spare)
    free = full.copy()
    active = {frozenset(e) for e in layout["blank"] + layout["marked"]}
    guide = set(layout["ports"])
    free.remove_edges_from(
        [
            e
            for e in full.edges
            if frozenset(e) in cuts or ((guide & set(e)) and frozenset(e) not in active)
        ]
    )
    assert set(nx.node_connected_component(free, hub)) == guide
    assert set(nx.node_connected_component(free, x)) == {x, y}
    assert set(nx.node_connected_component(free, stock[0])) == set(stock)
    # Everything outside the guide is isolated during flight except for the
    # two internal dimer bonds. Zero writes have zero native operator.
    np.testing.assert_array_equal(blocks()[0, 1], np.zeros((18, 18)))
    # Stored endpoint y has a charged deviation. Return its spectator dimer
    # to its initial state before touching the bridge.
    blank = blocks()["blank"] / 8
    period = expm(-2j * np.pi * blank)
    np.testing.assert_allclose(period, np.eye(9), atol=8e-14)
    write = LAW.echo()
    filled = LAW.SRC.index((1, 1))
    initial = np.eye(99)[:, filled]
    stored = write @ initial
    assert np.count_nonzero(abs(stored) > 1e-13) == 2
    np.testing.assert_allclose(write @ stored, initial, atol=1e-13)
    # Restore precisely the two endpoint loads after all intermediate swaps.
    route = tuple(reversed((*path, x)))
    word = list(zip(route, route[1:]))
    labels = {p: p for p in route}
    for a, b in (*word, *reversed(word[:-1])):
        labels[a], labels[b] = labels[b], labels[a]
    assert labels[near] == x and labels[x] == near
    assert all(labels[p] == p for p in route[1:-1])
    return dict(
        hub=hub,
        near=near,
        route=route,
        probe_dimer=[x, y],
        stock=stock,
        probe_cut_seams=len(cuts),
        erased_probe_private_references=sum(int(full.degree[p] == 5) for p in added),
    )


def conjugation():
    """Check complementing every load, direction and contrast on the raw law."""
    perm = [LAW.SRC.index((2 - a, 2 - b)) for a, b in LAW.SRC]
    for c, d in LAW.LABELS:
        base = 9 + 9 * LAW.LABELS.index((1 - c, -d))
        perm.extend(base + LAW.SRC.index((2 - a, 2 - b)) for a, b in LAW.SRC)
    np.testing.assert_array_equal(np.array(perm)[perm], np.arange(99))
    for op in (LAW.echo(), *LAW.local_term_parts()):
        np.testing.assert_array_equal(op[np.ix_(perm, perm)], op)
    return len(perm)


def station_gate(near, operator):
    """Lift a physical two-port operation to (carrier, near probe, far probe)."""
    untouched = 2 if near == (0, 1) else 0
    target = np.zeros((27, 27), complex)
    for k, digits in enumerate(DIGITS):
        col = LAW.SRC.index(tuple(digits[j] for j in near))
        for i, pair in enumerate(LAW.SRC):
            out = list(digits)
            out[near[0]], out[near[1]] = pair
            assert out[untouched] == digits[untouched]
            target[DIGITS.index(tuple(out)), k] = operator[i, col]
    return target


@cache
def pulses(kappa):
    blank = blocks(kappa)["blank"] / 8
    splitter = station_gate((1, 2), expm(-0.5j * np.pi * blank))
    contact = station_gate((0, 1), expm(-1j * DWELL * blank))
    return splitter, contact


def local_assay(q, s, kappa, angle, contact=True):
    """Record the two probe rails and failure to retain the initial carrier load."""
    split, gate = pulses(kappa)
    state = split[:, DIGITS.index((1 + q, 1, 1 + s))]
    if contact:
        state = gate @ state
    else:
        # Both probe rails stay in the disconnected waiting dimer. The
        # carrier still sees the blank spoke with a neutral neighbor.
        h = blocks(kappa)["blank"] / 8
        i = LAW.SRC.index((1 + q, 1))
        state = expm(-1j * DWELL * h)[i, i] * state
    state = state * np.array([a == 1 + q for a, _, _ in DIGITS])
    state *= np.array([np.exp(1j * angle) if b == 1 + s else 1 for _, b, _ in DIGITS])
    out = split.conj().T @ state
    plus = abs(out[DIGITS.index((1 + q, 1, 1 + s))]) ** 2
    minus = abs(out[DIGITS.index((1 + q, 1 + s, 1))]) ** 2
    np.testing.assert_allclose(plus + minus, np.vdot(state, state).real, atol=1e-14)
    return np.array([plus, minus, max(0.0, 1 - plus - minus)])


def bin_interval(c, n):
    # Bonferroni over all 36 scan rows and both quadratures/two output rails.
    tail = 0.01 / 512
    return np.array(
        [
            beta.ppf(tail, c, n - c + 1) if c else 0,
            beta.isf(tail, c + 1, n - c) if c < n else 1,
        ]
    )


def phase_read(q, s, k, code, arrival, launch, contact, rng):
    banks = []
    quadratures = []
    for angle in (0, -np.pi / 2):
        response = local_assay(q, s, k, angle, contact)
        # Use every original trial. Apertures independent of probe rail since
        # the assay acts on loads only, and the retained hub charge is sharp.
        probs = [
            code * response[0],
            code * response[1],
            (launch * arrival - code) * sum(response[:2]),
            launch * arrival * response[2],
            launch * (1 - arrival),
            1 - launch,
        ]
        counts = rng.multinomial(SHOTS, np.maximum(probs, 0))
        accepted = int(counts[0] + counts[1])
        if accepted == 0:
            raise ValueError("No probe outcomes for the phase estimate")
        quadratures.append(2 * bin_interval(counts[0], accepted) - 1)
        banks.append(counts.tolist())
    x, y = quadratures
    # These experimental confidence rectangles all fit in a single angular
    # semicircle, including the isolated control whose phase is zero.
    corners = [float(np.arctan2(b, a)) for a in x for b in y]
    if max(corners) - min(corners) >= np.pi:
        raise ValueError("The measured phase wraps or is unresolved")
    return dict(
        counts=banks,
        phase99=[min(corners), max(corners)],
        accepted_prob=code * sum(local_assay(q, s, k, 0, contact)[:2]),
    )


def filtered_probe():
    """The independently frozen finite energy clock from the original source."""
    clock = json.loads((HERE.parent / "results.json").read_text())["clock"]
    e, amplitudes, weights, *_ = build()
    frequency = (e - clock["center"]) * clock["tick"] / (2 * np.pi)
    n = clock["number_delays"]
    response = (
        np.exp(-1j * np.pi * (n - 1) * frequency)
        * np.sinc(n * frequency)
        / np.sinc(frequency)
    )
    K = amplitudes @ response
    admitted = float(np.linalg.norm(weights * response) ** 2)
    rng = np.random.default_rng(269125)
    result = []
    for q, s in product((-1, 1), repeat=2):
        bright = local_assay(q, s, 8, -np.pi / 2)[0]
        # Accepted local head and probe instruments have a scalar action on
        # one another. Trace the unobserved starting-owner head in bank().
        observed = bank(K * np.sqrt(bright), admitted, rng)
        result.append(dict(charge=q, probe=s, **observed))
    return result


def run():
    geometry = physical_probe()
    rawdim = conjugation()
    rows = []
    flight = []
    rng = np.random.default_rng(326711)
    for k in (8, 4):
        minus = native_graph(-1, k)
        for q in (-1, 1):
            h = native_graph(q, k)
            assert (minus != h).nnz == 0
            source = np.zeros((h.shape[0], 2))
            source[ONE] = BASE
            final = expm_multiply(-1j * FLIGHT * h, source)
            K, _, hub, _ = channel(final, 0)
            yield_, lowest, _ = pt(K)
            arrival = np.linalg.norm(hub) ** 2 / 2
            flight.append(
                dict(
                    kappa=k,
                    charge=q,
                    arrival=arrival,
                    code_yield=yield_ / 9,
                    pt_min=lowest / 9,
                )
            )
            for s in (-1, 1):
                for enabled in (True, False):
                    row = phase_read(q, s, k, yield_ / 9, arrival, 1 / 3, enabled, rng)
                    rows.append(
                        dict(kappa=k, charge=q, probe=s, contact=enabled, **row)
                    )
        # Neutral reference, including its actual probe and source-retention
        # losses. This panel calibrates only the response of an empty hub.
        for s in (-1, 1):
            row = phase_read(0, s, k, 1, 1, 1, True, rng)
            rows.append(dict(kappa=k, charge=0, probe=s, contact=True, **row))
    shifts = []
    for row in rows:
        if not row["charge"] or not row["contact"]:
            continue
        ref = next(
            r
            for r in rows
            if r["kappa"] == row["kappa"]
            and r["charge"] == 0
            and r["probe"] == row["probe"]
        )
        lo, hi = np.subtract(row["phase99"], ref["phase99"][::-1])
        shifts.append(
            dict(
                kappa=row["kappa"],
                charge=row["charge"],
                probe=row["probe"],
                shift99=[lo, hi],
            )
        )
    return dict(
        raw_local_dimension=rawdim,
        geometry=geometry,
        flight_time=FLIGHT,
        probe_time=DWELL,
        shots_per_quadrature=SHOTS,
        columns=[
            "probe_plus_code",
            "probe_minus_code",
            "other_hub_records",
            "lost_hub_charge",
            "no_hub_arrival",
            "source_refusal",
        ],
        flight=flight,
        rows=rows,
        shifts=shifts,
        filtered_probe=filtered_probe(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    if parser.parse_args().write:
        (HERE / "results.json").write_text(json.dumps(run(), indent=2) + "\n")
    else:
        print(json.dumps(run()["shifts"], indent=2))
