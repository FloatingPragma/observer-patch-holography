"""All mobile Gauss backgrounds on a fully open twelve-port written ring."""

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import currents
from ..belt_walls.pinned_string import core_read
from . import isolated_ring, ring


def bits(mask, length):
    return (int(mask) >> np.arange(length)) & 1


def mask(sites):
    return sum(1 << int(s) for s in sites)


def layout(length):
    for level in range(3):
        board, paths = isolated_ring.geometry(level)
        if length in paths:
            return board, paths[length]
    raise ValueError("No induced physical ring of the required length")


@lru_cache
def arcs(length, n):
    out = np.full(2**length, np.nan)
    for s in combinations(range(length), n):
        out[mask(s)] = core_read(s, length, n - length // 3 + 1)
    return out


def block(g, ref):
    length = len(g)
    w, heads, loads = ring.basis(g)
    n = g.sum()
    chart = arcs(length, n)
    a = chart[loads @ (1 << np.arange(length))]
    low, high = np.nanmin(chart), np.nanmax(chart)
    tight = a == low
    far = a >= (low + high) / 2
    # Cache the sparse flip matrix. A flip toggles the head and transports
    # one load in the direction permitted by its current orientation.
    left = []
    right = []
    for j in range(length):
        dest = w ^ (1 << j)
        ix = np.searchsorted(w, dest)
        ok = (ix < len(w)) & (w[np.minimum(ix, len(w) - 1)] == dest) & (w < dest)
        left.extend(np.flatnonzero(ok))
        right.extend(ix[ok])
    edges = (
        np.concatenate((left, right)).astype(int),
        np.concatenate((right, left)).astype(int),
    )
    adj = sparse.csr_matrix((np.ones(len(edges[0])), edges), shape=(len(w), len(w)))
    degree = np.asarray(adj.sum(axis=1)).ravel()
    return w, heads, loads, a, tight, far, adj, degree, 24 * (loads @ ref)


def generator(data, extension):
    adj, degree, ref = data[-3:]
    return sparse.diags(ref + (4 if extension else -4) * degree) - 4 * adj


def peak(e, v, tight):
    scores = np.sum(v[tight] ** 2, axis=0)
    clusters = np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1e-8) + 1)
    repeated = {}
    for group in clusters:
        if len(group) == 1:
            continue
        q = v[tight][:, group]
        s, rot = np.linalg.eigh(q.T @ q)
        scores[group] = -1
        scores[group[0]] = s[-1]
        # Fix any basis ambiguity in equally optimal states by projecting a
        # sharp word on the maximizing subspace, in enumeration order.
        sub = v[:, group] @ rot[:, s >= s[-1] - 1e-10]
        row = int(np.flatnonzero(np.sum(sub**2, axis=1) > 1e-12)[0])
        z = sub @ sub[row]
        repeated[int(group[0])] = z / np.linalg.norm(z)
    i = int(np.flatnonzero(scores >= scores.max() - 1e-10)[0])
    return float(e[i]), repeated.get(i, v[:, i])


def terminal(data, ext):
    w, heads, loads, a, tight, far, *_ = data
    h = generator(data, ext)
    e, v = eigh(h.toarray(), driver="evd", check_finite=False)
    energy, state = peak(e, v, tight)
    prob = state**2
    assert np.linalg.norm(h @ state - energy * state) < 2e-7
    far_values = np.linalg.eigvalsh(h[far][:, far].toarray()) if far.any() else []
    nonwall = h.diagonal()[~np.isfinite(a)]
    feature = dict(
        energy=energy,
        tight=float(prob[tight].sum()),
        any_wall=float(prob[np.isfinite(a)].sum()),
        stretched=float(prob[far].sum()),
        stretched_interval=[float(far_values[0]), float(far_values[-1])]
        if len(far_values)
        else None,
        nonwall_trial_interval=[float(nonwall.min()), float(nonwall.max())]
        if len(nonwall)
        else None,
    )
    return feature, state


def candidates(n, ref):
    length = len(ref)
    rows = []
    stats = Counter()
    for sites in combinations(range(length), n):
        key = mask(sites)
        g = bits(key, length)
        data = block(g, ref)
        tight, far = data[4:6]
        stats["backgrounds"] += 1
        stats["joint_words"] += len(data[0])
        stats["with_tight"] += int(tight.any())
        stats["with_stretched"] += int(far.any())
        stats["with_both"] += int(tight.any() and far.any())
        row = [key, len(data[0]), int(tight.sum()), int(far.sum())]
        for ext in (False, True):
            if tight.any():
                result, _ = terminal(data, ext)
                row.extend([round(result["tight"], 11), round(result["energy"], 11)])
            else:
                row.extend([0, None])
        rows.append(row)
    return dict(statistics=dict(stats), rows=sorted(rows))


def run():
    geometry = {}
    for length in (12, 15, 18):
        b, path = layout(length)
        geometry[str(length)] = isolated_ring.certificate(b, path)
    ref = np.array(geometry["12"]["reference"])
    scans = [candidates(n, ref) for n in (4, 5)]
    paired = {}
    for column, name in [(4, "native"), (6, "extension")]:
        plus = {r[0]: r for r in scans[1]["rows"] if r[2] and r[3]}
        possible = []
        for r in scans[0]["rows"]:
            if not (r[2] and r[3]):
                continue
            for site in range(12):
                p = r[0] | (1 << site)
                if p != r[0] and p in plus:
                    possible.append((min(r[column], plus[p][column]), r[0], p, site))
        # Scores within numerical tolerance share the same rank.
        best = max(s[0] for s in possible)
        _, minus, positive, site = min(
            (s for s in possible if s[0] >= best - 1e-9), key=lambda s: s[1:]
        )
        points = []
        for length in (12, 15):
            for n, g in [(length // 3, minus), (length // 3 + 1, positive)]:
                key = g if length == 12 else g | (1 << 12)
                data = block(
                    bits(key, length), np.array(geometry[str(length)]["reference"])
                )
                if data[4].any():
                    result, _ = terminal(data, name == "extension")
                else:
                    result = dict(tight=0)
                points.append(
                    dict(
                        length=length,
                        n=n,
                        gauss=key,
                        dimension=len(data[0]),
                        tight_rows=int(data[4].sum()),
                        stretched_rows=int(data[5].sum()),
                        **result,
                    )
                )
        paired[name] = dict(
            eligible_pairs=len(possible), additional_site=site, selected=points
        )
    return currents._portable(
        dict(
            geometry=geometry,
            census=scans,
            columns=[
                "Gauss mask",
                "dimension",
                "tight rows",
                "stretched rows",
                "native maximum",
                "native energy",
                "extension maximum",
                "extension energy",
            ],
            paired=paired,
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    r = run()
    with a.output.open("x") as f:
        f.write(json.dumps(r, indent=2, allow_nan=False) + "\n")
