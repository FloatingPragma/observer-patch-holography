"""Matched signed charge and head histories for the alternating 144-port material."""

import argparse
from functools import cache
import json
from math import lcm
from pathlib import Path

import networkx as nx
import numpy as np
from scipy.linalg import eigh_tridiagonal

from model.executable_model import federation
from model.experiments.interval_echo_angular.charged_orientation.winding_memory.measure import (
    source,
)
from ...conjugate_flight.measure import edge
from ...marked_bands import material, stock
from ..measure import local_sector
from .certify import DENOM, enclose

N = 144
RATIOS = np.array([206.7682827, 1776.93 / 0.51099895069])
SPACINGS = tuple(s for s in range(1, 49) if not 48 % s and 3 * lcm(2, 3 * s) <= N)
SHOTS = 2**20
WAITS = (0, 1, 16, 256, 4096, 65536)


@cache
def apparatus():
    ring, arms = material.geometry()
    b = federation(2)
    graph = nx.Graph(b.seams)
    active = [(ring[j], ring[j + 1]) for j in range(0, N, 2)]
    bypass = graph.copy()
    bypass.remove_edges_from(active + arms)
    assert nx.is_biconnected(bypass)
    # A remote stock pair per possible deposit, plus the source compensator.
    ds = nx.multi_source_dijkstra_path_length(graph, ring)
    pairs, reserved = [], set()
    for a, bp in sorted(graph.edges()):
        if min(ds[a], ds[bp]) > 2 and not reserved.intersection((a, bp)):
            pairs.append([a, bp])
            reserved.update((a, bp))
            if len(pairs) == 49:
                break
    assert len(pairs) == 49
    assert set(map(int, b.degree[ring[1::3]])) == {5}
    outside_ref = next(
        p
        for p in range(b.n)
        if p not in ring and p not in reserved and b.degree[p] == 5
    )
    # These native source routines check the full local columns of both deposits.
    source()
    for endpoint in (0, 1):
        stock.write_column(endpoint)
    route = {edge(a, bp) for a, bp in zip(ring, ring[1:] + ring[:1])}
    return ring, arms, graph, route, pairs, outside_ref


def preparation(s):
    ring, arms, g, route, stock_pairs, ref = apparatus()
    chosen = arms[::s] if s else []
    marks = {edge(*a) for a in chosen}
    store, rest = stock_pairs[-1]
    active = {edge(ring[j], ring[j + 1]) for j in range(0, N, 2)}
    # One writes the active and side heads using untouched bypasses, temporarily
    # transposing stock into the requested ports and undoing the surviving loads.
    # The common biconnected bypass supplies a local phase rail avoiding the
    # first endpoint. Closing unused edges is done only after these deposits.
    cuts = sorted(
        edge(a, b)
        for a, b in g.edges()
        if {a, b}.intersection(set(ring) | {store, rest})
        and edge(a, b) not in route | marks | {edge(store, rest)}
    )
    charge = ring[0]
    graph = g.copy()
    graph.remove_edges_from(marks)
    for a, b in cuts:
        for here, park in [(charge, ring[-1]), (store, rest)]:
            if here in (a, b):
                assert park not in (a, b) and graph.has_edge(here, park)
        graph.remove_edge(a, b)
    assert set(nx.node_connected_component(graph, charge)) == set(ring)
    assert set(nx.node_connected_component(graph, store)) == {store, rest}
    return dict(
        marked_sides=[list(pair) for pair in chosen],
        active_edges=len(active),
        private_zero_writes=48,
        exterior_zero_writes=len(cuts),
        phase_rail=ref,
        spent_vacancy_excess_pairs=stock_pairs[: len(chosen)],
        stock_dimer=[store, rest],
    )


@cache
def contacts(sign):
    active = local_sector(sign, (1, 1))
    passive = local_sector(sign, (2, 1))
    plus = np.array([[1, 0], [0, 1], [1, 0], [0, 1]]) / np.sqrt(2)
    np.testing.assert_allclose(passive @ plus, 32 * plus)
    blank = local_sector(sign, "blank")
    return blank, active


def finite(s, length, sign):
    """Expand the two incident native seam columns, including the inaccessible wall."""
    ring = apparatus()[0]
    blank, active = contacts(sign)
    # Physical ordered native addresses orient each coin. Both signs use the
    # same operator sector and opposite directed source head preparations.
    initial = 0
    for j in range(0, N, 2):
        side = int(ring[j] > ring[j + 1])
        before = next(h for h in (0, 1) if active[2 * h + side, 2 * h + side])
        initial |= before << (j // 2)
    if length < N:
        initial ^= 1 << (length // 2)
    labels = []
    heads = initial
    for x in range(-1, length + 1):
        labels.append((x % N, heads))
        if x >= 0 and not x % 2:
            heads ^= 1 << (x // 2)
    index = {w: i for i, w in enumerate(labels)}
    assert len(index) == length + 2
    matrix = np.zeros((len(labels), len(labels)))
    for col, (pos, heads) in enumerate(labels):
        for a in ((pos - 1) % N, pos):
            b = (a + 1) % N
            side = int(ring[pos] > ring[b if pos == a else a])
            if a % 2:
                matrix[col, col] += blank[side, side]
                other = b if a == pos else a
                matrix[index[other, heads], col] += blank[1 - side, side]
            else:
                h = (heads >> (a // 2)) & 1
                i = 2 * h + side
                for o in np.flatnonzero(active[:, i]):
                    other = pos if o % 2 == side else (b if pos == a else a)
                    new_heads = heads ^ ((h ^ (int(o) // 2)) << (a // 2))
                    matrix[index[other, new_heads], col] += active[o, i]
        if s and (pos - 1) % (3 * s) == 0:
            matrix[col, col] += 32
    np.testing.assert_array_equal(matrix, matrix.T)
    d, e = matrix.diagonal(), matrix.diagonal(1)
    np.testing.assert_array_equal(matrix, np.diag(d) + np.diag(e, 1) + np.diag(e, -1))
    return d, e


def bands(s):
    """For a scalar periodic Jacobi chain band edges lie at transfer trace ±2."""
    p = lcm(2, 3 * s) if s else 2
    out = []
    for twist in (1, -1):
        h = np.zeros((p, p))
        for a in range(p):
            b = (a + 1) % p
            t = 0.5 if a % 2 else 4.0
            h[a, a] += t
            h[b, b] += t
            h[a, b] -= t * (twist if b == 0 else 1)
            h[b, a] -= t * (twist if b == 0 else 1)
            if s and (a - 1) % (3 * s) == 0:
                h[a, a] += 32
        out.append(np.linalg.eigvalsh(h))
    edges = []
    if s:
        for twist, eig in zip((1, -1), out):
            diag = [9 + 64 * int((j - 1) % (3 * s) == 0) for j in range(p)]
            links = [-8 if j % 2 == 0 else -1 for j in range(p)]
            edges.append(
                [
                    enclose(diag, links[:-1], twist * links[-1], i, val)
                    for i, val in enumerate(eig)
                ]
            )
    return np.minimum(*out), np.maximum(*out), edges


def comparison(v):
    """Compare two independent heavier levels anchored to the fixed lowest energy."""
    ratios = v / v[0]
    # Joint minimax relative error, requiring ordered distinct heavy bands.
    upper = np.triu_indices(len(v) - 1, 1)
    j, k = upper[0] + 1, upper[1] + 1
    a = np.abs(ratios[j] / RATIOS[0] - 1)
    b = np.abs(ratios[k] / RATIOS[1] - 1)
    error = np.maximum(a, b)
    # Total relative error breaks minimax ties among two candidates with
    # the same limiting heavy line.
    i = np.lexsort((k, j, a + b, error))[0]
    return dict(
        band_indices=[0, int(j[i]), int(k[i])],
        ratios=ratios[[0, j[i], k[i]]].tolist(),
        joint_error=float(error[i]),
        first_three_ratios=ratios[:3].tolist(),
    )


def flight(s, vec, energies, lo, hi, rng):
    p = lcm(2, 3 * s)
    # Sharp source in the middle of the second barrier-free interval, reached
    # by directed native swaps from the initial launch. Detector counts the
    # entire next material period by addresses, tracing over all head outcomes.
    x0 = 1 + p // 2
    x = np.arange(-1, N + 1)
    mask = ((x % N) >= p + 1) & ((x % N) < 2 * p + 1)
    rows = []
    dt = 1 / 8
    separation = lo[1] - hi[0]
    clock_size = 2 ** int(np.ceil(np.log2(32 / (dt * separation))))
    center = (lo[0] + hi[0]) / 2
    angles = (energies - center) * dt / (2 * np.pi)
    filt = (
        np.exp(-1j * np.pi * (clock_size - 1) * angles)
        * np.sinc(clock_size * angles)
        / np.sinc(angles)
    )
    amplitudes = filt * vec[x0 + 1]
    admitted = float(np.vdot(amplitudes, amplitudes).real)
    first_band = energies < (lo[1] + hi[0]) / 2
    purity = float(sum(abs(amplitudes[first_band]) ** 2) / admitted)
    filtered = []
    # Compare the same prepared state immediately and after free flight. The
    # phase reference has a finite binary gate register; every refusal counts.
    for t in (*WAITS, 1048576, 16777216):
        ps = vec @ (amplitudes * np.exp(-1j * t * (energies - energies[0])))
        joint = float(sum(abs(ps[mask]) ** 2))
        counts = rng.multinomial(SHOTS, [joint, admitted - joint, 1 - admitted])
        filtered.append(
            dict(
                time=t,
                arrived=int(counts[0]),
                elsewhere=int(counts[1]),
                refused=int(counts[2]),
                joint_probability=joint,
            )
        )
    for t in WAITS:
        state = vec @ (np.exp(-1j * t * (energies - energies[0])) * vec[x0 + 1])
        prob = sum(abs(state[mask]) ** 2)
        count = int(rng.binomial(SHOTS, np.clip(prob, 0, 1)))
        rows.append(
            dict(time=t, arrivals=count, other=SHOTS - count, probability=float(prob))
        )
    assert not mask[x0 + 1]
    return dict(
        launch_ring_index=x0,
        clock=dict(
            step=dt,
            register_size=clock_size,
            center=center,
            admitted=admitted,
            first_band_weight=purity,
            banks=filtered,
        ),
        receiver_ring_indices=np.flatnonzero(
            (np.arange(N) >= p + 1) & (np.arange(N) < 2 * p + 1)
        ).tolist(),
        banks=rows,
    )


def run():
    rng = np.random.default_rng(20261126)
    report = {}
    for s in (0, *SPACINGS):
        p = lcm(2, 3 * s) if s else 2
        lo, hi, certified = bands(s)
        short = p * (N // (2 * p))
        lengths = sorted(set([short, 2 * short, N]))
        finite_rows = {}
        for length in lengths:
            d, e = finite(s, length, -1)
            other = finite(s, length, 1)
            np.testing.assert_array_equal(d, other[0])
            np.testing.assert_array_equal(e, other[1])
            vals, vec = eigh_tridiagonal(d, e)
            finite_rows[str(length)] = dict(
                states=len(vals),
                minimum=float(vals[0]),
                minus_bulk=float(vals[0] - lo[0]),
            )
            finite_rows[str(length)]["ground_interval_numerator"] = enclose(
                np.rint(2 * d).astype(int).tolist(),
                np.rint(2 * e).astype(int).tolist(),
                0,
                0,
                vals[0],
            )
            assert vals[0] >= lo[0] - 1e-11
        data = dict(
            period=p,
            periods_on_board=N // p,
            band_minima=lo.tolist(),
            band_tops=hi.tolist(),
            finite=finite_rows,
            preparation=preparation(s),
        )
        if s:
            data["periodic_edge_intervals_numerator"] = certified
            data["joint_family"] = comparison(lo)
            data["finite_family"] = comparison(vals)
            # Evaluate every distinct triple of bulk band bottoms, including
            # the case in which the lightest chosen member is excited.
            candidates = [(comparison(lo[i:]), i) for i in range(len(lo) - 2)]
            best, offset = min(candidates, key=lambda pair: pair[0]["joint_error"])
            data["any_bulk_triple"] = dict(starting_band=offset, **best)
            data["sharp_flight"] = flight(s, vec, vals, lo, hi, rng)
        report[str(s)] = data
    return dict(
        target_ratios=RATIOS.tolist(),
        certificate_denominator=DENOM,
        ports=apparatus()[0],
        spacing_units="eligible middle ports",
        shots=SHOTS,
        settings=report,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    path = p.parse_args().output
    result = run()
    with path.open("x") as stream:
        json.dump(result, stream, separators=(",", ":"))
        stream.write("\n")
