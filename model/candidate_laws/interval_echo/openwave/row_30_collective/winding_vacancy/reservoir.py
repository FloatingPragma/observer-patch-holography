"""Native all-filling reservoir and charge transfer spectra on finite winding stock."""

import argparse
from dataclasses import dataclass
from hashlib import sha256
import io
import json

import networkx as nx
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.special import logsumexp

from . import geometry, measure

HERE = measure.HERE
PROTOCOL = "ae830d040298b5b078225d12e7931c738cb69e6a334f457df1db204887b3e6ce"
SHOTS = 2**22


def native_contacts():
    src = measure.source.measure
    tables = src.blocks()
    local = {}
    for name, key in (("blank", "blank"), ("even", (1, 1)), ("idle", (0, 1))):
        block = tables[key]
        heads = range(len(block) // 9)
        sites = [
            9 * c + src.LAW.SRC.index((1 - a, 1 - b))
            for c in heads
            for a in range(2)
            for b in range(2)
        ]
        complement = [
            9 * (c ^ int(name == "even")) + src.LAW.SRC.index((a, b))
            for c in heads
            for a in range(2)
            for b in range(2)
        ]
        assert not block[
            np.ix_(np.setdiff1d(np.arange(len(block)), sites), sites)
        ].any()
        assert np.array_equal(
            block[np.ix_(sites, sites)], block[np.ix_(complement, complement)]
        )
        baseline = block[sites[0], sites[0]]
        local[name] = block[np.ix_(sites, sites)] - baseline * np.eye(
            len(sites), dtype=int
        )
        assert np.all(local[name][:, ::4] == 0)
    assert not local["idle"].any()
    return local


def integer_sector(n, number, tables):
    masks = np.arange(1 << n)
    masks = masks[np.array([int(m).bit_count() for m in masks]) == number]
    words = (masks + (np.arange(1 << (n // 2)) << n)[:, None]).ravel()
    return words, compile_words(n, words, tables)


def compile_words(n, words, tables):
    columns = np.arange(len(words))
    rr, cc, values = [], [], []
    for left in range(n):
        right = (left + 1) % n
        a = (words >> left) & 1
        b = (words >> right) & 1
        marked = left % 2
        j = left // 2
        coin = ((words >> (n + j)) & 1) if marked else 0
        sub = tables["even" if marked else "blank"]
        start = 4 * coin + 2 * a + b
        for end, old in zip(*np.nonzero(sub), strict=True):
            at = columns[start == old]
            if not len(at):
                continue
            new = (
                words[at]
                ^ ((((old // 2) ^ (end // 2)) & 1) << left)
                ^ (((old ^ end) & 1) << right)
            )
            if marked:
                new ^= ((old // 4) ^ (end // 4)) << (n + j)
            rows = np.searchsorted(words, new)
            np.testing.assert_array_equal(words[rows], new)
            rr.append(rows)
            cc.append(at)
            values.append(np.full(len(at), sub[end, old], dtype=np.int32))
    if not values:
        return coo_matrix((len(words), len(words)), dtype=int).tocsr()
    h = coo_matrix(
        (np.concatenate(values), (np.concatenate(rr), np.concatenate(cc))),
        shape=(len(words), len(words)),
    ).tocsr()
    h.eliminate_zeros()
    assert not (h - h.T).nnz
    return h


@dataclass
class Sector:
    words: np.ndarray
    eigenvalues: np.ndarray
    blocks: list
    components: int
    largest: int


def diagonalize(n, number, tables, retain):
    words, h = integer_sector(n, number, tables)
    count, tags = connected_components(h, directed=False)
    vals = []
    blocks = []
    largest = 0
    for t in range(count):
        ix = np.flatnonzero(tags == t)
        largest = max(largest, len(ix))
        local = h[np.ix_(ix, ix)].toarray() / 8
        if retain:
            e, v = eigh(local)
            blocks.append((ix, e, v))
        else:
            e = eigh(local, eigvals_only=True)
        vals.append(e)
    return Sector(words, np.sort(np.concatenate(vals)), blocks, count, largest)


def ensembles(n, sectors, raw, rng):
    rows = []
    fillings = np.arange(n + 1)
    tol = np.sqrt(np.log(2 * 100 / 0.01) / (2 * SHOTS))
    for beta in [1 / 32, 1 / 16, 1 / 8]:
        lz = np.array([logsumexp(-beta * s.eigenvalues) for s in sectors])
        appointments = [(0.0, 0.0), (4.5, 0.0)]
        if beta == 1 / 16:
            appointments += [(8.0, 0.0), (0.0, 1.75), (1.75, 1.75)]
        for mu, shift in appointments:
            lweights = lz + beta * (mu - shift) * fillings
            p = np.exp(lweights - logsumexp(lweights))
            y = rng.multinomial(SHOTS, p)
            raw[f"n{n}_b{beta}_mu{mu}_shift{shift}"] = y
            N = fillings @ (y / SHOTS)
            exact = fillings @ p / n
            number = n // 2
            work = -np.log(y[[number - 1, number + 1]] / y[number]) / beta
            expected = -(
                lz[[number - 1, number + 1]] - lz[number]
            ) / beta + shift * np.array([-1, 1])
            rows.append(
                dict(
                    beta=beta,
                    mu=mu,
                    uniform_onsite=shift,
                    measured_filling=N / n,
                    exact_filling=exact,
                    simultaneous_filling_radius=tol,
                    half_within_bound=abs(N / n - 0.5) <= tol,
                    exact_canonical_work_remove_add=expected.tolist(),
                    measured_canonical_work_remove_add=(
                        work + mu * np.array([-1, 1])
                    ).tolist(),
                )
            )
    return rows


def ground_input(sector):
    low = sector.eigenvalues[0]
    basis = []
    for ix, e, v in sector.blocks:
        for j in np.flatnonzero(e < low + 1e-9):
            x = np.zeros(len(sector.words))
            x[ix] = v[:, j]
            basis.append(x)
    return np.stack(basis, axis=1)


def charge_kick(n, k, parity, source, dest, ground, sign):
    L = n // 2
    y = np.zeros((len(dest.words), ground.shape[1]), dtype=complex)
    for c in range(L):
        v = 2 * c + parity
        take = ((source.words >> v) & 1) == int(sign < 0)
        output = source.words[take] ^ (1 << v)
        rows = np.searchsorted(dest.words, output)
        np.testing.assert_array_equal(dest.words[rows], output)
        y[rows] += np.exp(1j * k * c) / L * ground[take]
    return y


def response(n, sectors, raw):
    N = n // 2
    src = sectors[N]
    G = ground_input(src)
    rows = []
    for sign in [-1, 1]:
        dest = sectors[N + sign]
        for j in range(N):
            k = 2 * np.pi * j / N
            # k measures the cell carrier, and upper rest is at k=pi.
            q = abs(np.angle(-np.exp(1j * k)))
            expected = np.sqrt(3.5**2 + 8 * np.sin(q / 2) ** 2)
            for parity in (0, 1):
                y = charge_kick(n, k, parity, src, dest, G, sign)
                eall, wall = [], []
                for ix, e, v in dest.blocks:
                    eall.extend(e - src.eigenvalues[0])
                    wall.extend(np.sum(abs(v.T @ y[ix]) ** 2, axis=1) / G.shape[1])
                e = np.array(eall)
                w = np.array(wall)
                order = np.argsort(e)
                e, w = e[order], w[order]
                assert abs(w.sum() - np.sum(abs(y) ** 2) / G.shape[1]) < 1e-11
                raw[f"add_remove_{sign}_{j}_{parity}"] = np.column_stack((e, w))
                lines, inv = np.unique(np.round(e, 9), return_inverse=True)
                weights = np.bincount(inv, weights=w, minlength=len(lines))
                visible = weights > 1e-10
                peak = int(np.argmax(weights))
                rows.append(
                    dict(
                        sign=sign,
                        k=k,
                        q=q,
                        parity=parity,
                        kg_frequency=expected,
                        accepted=w.sum(),
                        visible_lines=int(sum(visible)),
                        first_visible=lines[visible][0],
                        dominant_line=lines[peak],
                        dominant_fraction=weights[peak] / w.sum(),
                        work_centroid=e @ w / w.sum(),
                        fraction_within_KG=weights[abs(lines - expected) < 0.02].sum()
                        / w.sum(),
                    )
                )
    return rows


def boards():
    board = geometry.federation(2)
    graph = nx.Graph(board.seams)
    _, reps = geometry.cycles(graph)
    survey = json.loads((HERE / "geometry.json").read_text())
    belt = next(r["route"] for r in survey["belts"] if r["qualified"])
    routes = [reps[4], reps[6], belt]
    assert [len(r) for r in routes] == [4, 6, 10]
    return [(r, geometry.qualify(board, graph, r)) for r in routes]


def run():
    assert sha256((HERE / "RESERVOIR.md").read_bytes()).hexdigest() == PROTOCOL
    tables = native_contacts()
    results = []
    raw = {}
    rng = np.random.default_rng(316144)
    for route, qual in boards():
        assert qual["qualified"]
        n = len(route)
        sectors = [
            diagonalize(n, N, tables, n == 10 and abs(N - n // 2) <= 1)
            for N in range(n + 1)
        ]
        for N in range(n + 1):
            np.testing.assert_allclose(
                sectors[N].eigenvalues, sectors[n - N].eigenvalues, atol=3e-12
            )
        minima = np.array([x.eigenvalues[0] for x in sectors])
        hist = ensembles(n, sectors, raw, rng)
        spect = response(n, sectors, raw) if n == 10 else []
        results.append(
            dict(
                ports=route,
                stock=qual,
                n=n,
                words=[len(c.words) for c in sectors],
                component_counts=[c.components for c in sectors],
                largest_component=[c.largest for c in sectors],
                sector_minima=minima.tolist(),
                ground_cost_remove_or_add=(
                    minima[[n // 2 - 1, n // 2 + 1]] - minima[n // 2]
                ).tolist(),
                half_ground_degeneracy=int(
                    sum(sectors[n // 2].eigenvalues < minima[n // 2] + 1e-9)
                ),
                equilibrium=hist,
                local_transfer_spectra=spect,
            )
        )
    zipbuf = io.BytesIO()
    np.savez_compressed(zipbuf, **raw)
    blob = zipbuf.getvalue()
    out = dict(
        protocol_sha256=PROTOCOL,
        counts_and_spectra_sha256=sha256(blob).hexdigest(),
        trials=SHOTS,
        rings=results,
        integer_blank=tables["blank"].tolist(),
        integer_marked=tables["even"].tolist(),
    )
    return json.dumps(
        out, indent=2, sort_keys=True, default=lambda x: x.item()
    ) + "\n", blob


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    opt = parser.parse_args()
    out, blob = run()
    for name, value in [
        ("reservoir.json", out.encode()),
        ("reservoir_counts.npz", blob),
    ]:
        if opt.write:
            (HERE / name).write_bytes(value)
        elif (HERE / name).read_bytes() != value:
            raise AssertionError(f"Reservoir replay mismatch: {name}")
    print("Verified all-filling winding reservoir and addition/removal spectra.")
