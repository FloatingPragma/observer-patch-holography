"""Unforced transfer between the two winding records on one conserved background."""

import argparse
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import currents
from . import backgrounds, isolated_ring, ring


def patch(length, site):
    b, path = backgrounds.layout(length)
    region = set(path)
    choices = set(b.neighbours[path[site]]) - region
    z = min(choices, key=lambda v: (len(region.intersection(b.neighbours[v])), v))
    spokes = [i for i, p in enumerate(path) if p in b.neighbours[z]]
    vertices = path + (z,)
    edges = [(i, (i + 1) % length) for i in range(length)] + [
        (i, length) for i in spokes
    ]
    assert {tuple(sorted((vertices[a], vertices[b]))) for a, b in edges} == {
        e for e in b.seams if set(e) <= set(vertices)
    }
    idx = {tuple(sorted(e)): i for i, e in enumerate(edges)}
    faces = []
    for f in combinations(range(len(vertices)), 3):
        keys = list(combinations(f, 2))
        if all(e in idx for e in keys):
            faces.append(tuple(idx[e] for e in keys))
    cert = isolated_ring.shield(b, vertices)
    assert cert["triangles_by_ring_ports"][3] == len(faces)
    return z, spokes, edges, faces, cert


def configurations(g, spokes):
    length = len(g)
    words, materials, entry, tight = [], [], [], []
    for k, u in enumerate([-1] + spokes):
        background = g.copy()
        if u >= 0:
            background[u] += 1
        w, _, n = ring.basis(background)
        weights = n @ (1 << np.arange(length))
        a = backgrounds.arcs(length, int(background.sum()))
        words.extend(w if u < 0 else w | (1 << (length + k - 1)))
        materials.extend(weights | (int(u < 0) << length))
        entry.extend([u] * len(w))
        tight.extend(a[weights] == np.nanmin(a))
    return tuple(np.asarray(x) for x in (words, materials, entry, tight))


def generator(words, material, edges, faces, ref, extension):
    index = {int(w): i for i, w in enumerate(words)}
    rows, cols, weights = [], [], []

    def bond(i, newword, newload, weight):
        j = index[newword]
        assert material[j] == newload
        if i < j:
            rows.append(i)
            cols.append(j)
            weights.append(weight)

    for i, (word, load) in enumerate(zip(words, material)):
        w, n = int(word), int(load)
        tips = []
        for e, (u, v) in enumerate(edges):
            src, dst = (u, v) if w & (1 << e) else (v, u)
            tips.append(src)
            if n & (1 << src) and not n & (1 << dst):
                bond(i, w ^ (1 << e), n ^ (1 << src) ^ (1 << dst), 4)
        if extension:
            for f in faces:
                if len({tips[e] for e in f}) == 3:
                    bond(i, w ^ sum(1 << e for e in f), n, 0.5)
    adj = sparse.coo_matrix(
        (weights, (rows, cols)), shape=(len(words), len(words))
    ).tocsr()
    adj += adj.T
    loads = (material[:, None] >> np.arange(len(ref))) & 1
    diagonal = (
        24 * (loads @ ref)
        + (1 if extension else -1) * np.asarray(adj.sum(axis=1)).ravel()
    )
    return sparse.diags(diagonal, dtype=float) - adj


def outcomes(p, entry, tight, site):
    return {
        name: [float(p[mask & tight].sum()), float(p[mask & ~tight].sum())]
        for name, mask in (
            ("reservoir_full", entry == -1),
            ("target_spoke", entry == site),
            ("other_spokes", (entry >= 0) & (entry != site)),
        )
    }


def panel(row, site, extension):
    length = row["length"]
    z, spokes, edges, faces, cert = patch(length, site)
    g = backgrounds.bits(row["gauss"], length)
    words, loads, entry, tight = configurations(g, spokes)
    ref = np.asarray(cert["reference"])
    h = generator(words, loads, edges, faces, ref, extension)
    count = np.count_nonzero(entry == -1)
    data = backgrounds.block(g, ref[:-1])
    assert np.array_equal(data[0], words[:count])
    np.testing.assert_array_equal([int(x).bit_count() for x in loads], g.sum() + 1)
    result = dict(
        length=length,
        neutral_gauss=row["gauss"],
        reservoir_port=z,
        spokes=spokes,
        triangle_count=len(faces),
        dimension=len(words),
        boundary_heads=len(cert["boundary_heads"]),
        triangles_by_patch_ports=cert["triangles_by_ring_ports"],
    )
    if not data[4].any():
        return result | {"source_available": False}
    _, state = backgrounds.terminal(data, extension)
    v = np.pad(state, (0, len(words) - count))
    e, basis = eigh(h.toarray(), driver="evd", check_finite=False)
    coeff = basis.T @ v
    groups = np.split(np.arange(len(e)), np.flatnonzero(np.diff(e) > 1e-8) + 1)
    mean = np.zeros(len(e))
    best_weight = -1
    for q in groups:
        u = basis[:, q] @ coeff[q]
        wt = float(u @ u)
        mean += u**2
        if wt > best_weight + 1e-12:
            best_weight, best_energy, best_p = wt, float(e[q[0]]), u**2 / wt
    times = [0, 1 / 16, 1 / 4, 1, 4, 16, 32]
    table = []
    for t in times:
        p = abs(basis @ (coeff * np.exp(-1j * t * e))) ** 2
        table.append(dict(time=t, read=outcomes(p, entry, tight, site)))
    # Remove every crossing matrix element but preserve its diagonal energy.
    # The reservoir read then commutes with the control Hamiltonian.
    rows, cols = h.nonzero()
    cut = entry[rows] != entry[cols]
    assert np.all((entry[rows[~cut]] < 0) == (entry[cols[~cut]] < 0))
    coupling = h[count:, :count] @ state
    onset = float(coupling @ coupling)
    prediction = 16 * float((state**2) @ (1 - data[2][:, spokes]).sum(axis=1))
    np.testing.assert_allclose(onset, prediction, atol=1e-9)
    return result | dict(
        source_available=True,
        initial_neutral_tight=row["tight"],
        transfer_t2_coefficient=onset,
        blocked_control_transfer=0,
        blocked_control_removed_pairs=int(np.sum(cut) // 2),
        terminal=table,
        infinite_time=outcomes(mean, entry, tight, site),
        largest_overlap_energy=best_energy,
        largest_overlap=best_weight,
        largest_overlap_read=outcomes(best_p, entry, tight, site),
    )


def run(scans):
    return currents._portable(
        dict(
            columns=["tight", "other"],
            laws={
                name: [
                    panel(r, pair["additional_site"], name == "extension")
                    for r in pair["selected"]
                    if r["n"] == r["length"] // 3
                ]
                for name, pair in scans["paired"].items()
            },
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    s = json.loads(Path(__file__).with_name("backgrounds.json").read_text())
    d = run(s)
    with args.output.open("x") as f:
        f.write(json.dumps(d, indent=2, allow_nan=False) + "\n")
