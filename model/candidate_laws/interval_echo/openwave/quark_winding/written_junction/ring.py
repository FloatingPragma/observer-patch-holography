"""Stationary and released wall sources on fully active closed native head strings."""

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.interval_echo import currents
from model.executable_model import carrier
from ..belt_walls.composite import gap_arc


def physical_cycles():
    board = carrier.federation(0)
    g = nx.Graph(carrier.cell().seams)
    owners = nx.Graph((a // 12, b // 12) for a, b in board.external)
    gateway = {
        (u // 12, v // 12): u % 12
        for a, b in board.external
        for u, v in ((a, b), (b, a))
    }
    paths = {}
    for tour in sorted(nx.cycle_basis(owners), key=lambda t: (len(t), t)):
        if len(tour) > 9:
            continue
        choices = []
        for j, k in enumerate(tour):
            a = gateway[k, tour[j - 1]]
            b = gateway[k, tour[(j + 1) % len(tour)]]
            options = {}
            for path in nx.all_simple_paths(g, a, b, cutoff=6):
                options.setdefault(len(path), [12 * k + i for i in path])
            choices.append(options)
        # Dynamic programming in total path length, one explicit witness per size.
        totals = {0: []}
        for opt in choices:
            updated = {}
            for n, p in totals.items():
                for m, part in sorted(opt.items()):
                    if n + m <= 18:
                        updated.setdefault(n + m, p + part)
            totals = updated
        for length in (12, 15, 18):
            if length in totals:
                paths.setdefault(length, totals[length])
        if len(paths) == 3:
            break
    if 12 not in paths:
        # Hamiltonian cycle of the isolated native cell at owner 0.
        path = next(p for p in nx.all_simple_paths(g, 0, 1) if len(p) == 12)
        paths[12] = path
    assert set(paths) == {12, 15, 18}
    for path in paths.values():
        assert len(set(path)) == len(path) and all(
            b in board.neighbours[a] for a, b in zip(path, path[1:] + path[:1])
        )
    return board, paths


def basis(gauss):
    length = len(gauss)
    words = []

    def extend(bits, last, first, i):
        if i == length:
            if last == first:
                words.append(bits)
            return
        for b in (0, 1):
            if gauss[i] + b - last in (0, 1):
                extend(bits | (b << i), b, first, i + 1)

    for start in (0, 1):
        extend(0, start, start, 0)
    w = np.array(sorted(words))
    heads = ((w[:, None] >> np.arange(length)) & 1).astype(np.int8)
    loads = gauss + heads - np.roll(heads, 1, axis=1)
    return w, heads, loads


def generator(w, loads, ref, extension):
    size, length = loads.shape
    index = {int(v): i for i, v in enumerate(w)}
    row, col = [], []
    for i, word in enumerate(w):
        for j in range(length):
            dest = index.get(int(word) ^ (1 << j))
            if dest is not None:
                row.append(i)
                col.append(dest)
    adj = sparse.csr_matrix((np.ones(len(row)), (row, col)), shape=(size, size))
    deg = np.asarray(adj.sum(axis=1)).ravel()
    # subtract common 8L for native, and constant host reference for both.
    diag = (4 if extension else -4) * deg + 24 * (loads @ ref)
    return sparse.diags(diag) - 4 * adj


def maxima(h, arcs):
    e, v = np.linalg.eigh(h.toarray())
    groups = np.split(np.arange(len(e)), np.where(np.diff(e) > 1e-8)[0] + 1)
    low = np.nanmin(arcs)
    out = {}
    for tag, mask in [
        ("any", np.isfinite(arcs)),
        ("tight", arcs == low),
        ("wider", arcs <= low + 3),
    ]:
        restricted = v[mask]
        weights = np.sum(abs(restricted) ** 2, axis=0)
        best = np.argmax(weights)
        maximum = weights[best]
        ground = None
        for ix in groups:
            if len(ix) == 1:
                weight = weights[ix[0]]
            else:
                b = restricted[:, ix]
                weight = np.linalg.norm(b, 2) ** 2
            if ix[0] == 0:
                ground = weight
            if weight > maximum:
                maximum = weight
                best = ix[0]
        out[tag] = dict(maximum=maximum, energy=e[best], ground_max=ground)
    return dict(minimum_arc=low, ground_energy=e[0], extrema=out)


def run():
    board, cycles = physical_cycles()
    result = []
    for length, path in sorted(cycles.items()):
        j = np.arange(length)
        n0 = (j % 3 == 0).astype(int)
        e0 = n0.copy()
        g0 = n0 - e0 + np.roll(e0, 1)
        # Exact geometric holding reference, including every port touching an external contact.
        ref = 6 - board.degree[path]
        for addition in (None, 1, 2):
            g = g0.copy()
            if addition is not None:
                g[addition] += 1
            w, heads, loads = basis(g)
            ones = 1 if addition is None else 2
            arcs = np.asarray(
                [
                    gap_arc(tuple(np.flatnonzero(n)), length, ones, 3 - ones)
                    for n in loads
                ],
                float,
            )
            for name in ("native", "positive_laplacian"):
                h = generator(w, loads, ref, name == "positive_laplacian")
                records = {}
                if addition is None:
                    idx = {int(word): i for i, word in enumerate(w)}
                    bulk = sum(1 << int(i) for i in j[e0 == 1])
                    for shift in (1, -2):
                        panels = []
                        for k in range(1, length // 3):
                            bits = bulk
                            for i in range(0, 3 * k, 3):
                                for step in range(abs(shift)):
                                    edge = (
                                        i + step if shift > 0 else i - step - 1
                                    ) % length
                                    bits ^= 1 << edge
                            p, q = idx[bits], idx[bulk]
                            diff = (loads[p] - loads[q]) @ ref
                            panels.append(
                                dict(
                                    periods=k,
                                    energy_change=h[p, p] - h[q, q],
                                    reference_change=24 * float(diff),
                                    contact_change=h[p, p] - h[q, q] - 24 * float(diff),
                                )
                            )
                            if k == length // 6:
                                psi = np.zeros((len(w), 2), complex)
                                psi[[p, q], [0, 1]] = 1
                                # Circular electrode centred on the left defect, disjoint from the right core at peak.
                                x = np.minimum(j, length - j)
                                f = (
                                    1
                                    + np.cos(
                                        np.pi * np.minimum(x / max(2.0, 3 * k / 2), 1)
                                    )
                                ) / 2
                                counts = loads @ f
                                potential_from_heads = heads @ (f - np.roll(f, -1))
                                for t in (0.0, 0.125, 1.0, 4.0):
                                    prob = abs(expm_multiply(-1j * t * h, psi)) ** 2
                                    diff = prob[:, 0] - prob[:, 1]
                                    value = counts @ diff
                                    measured = potential_from_heads @ diff
                                    assert abs(value - measured) < 1e-9
                                    panels[-1].setdefault("release", []).append(
                                        dict(time=t, material=value, heads=measured)
                                    )
                        records[str(shift)] = panels
                result.append(
                    dict(
                        length=length,
                        path=path,
                        extra_unit_at=addition,
                        law=name,
                        dimension=len(w),
                        wall_words=int(np.sum(np.isfinite(arcs))),
                        stationary=maxima(h, arcs),
                        separation=records,
                    )
                )
    return currents._portable(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    encoded = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as f:
        f.write(encoded)
