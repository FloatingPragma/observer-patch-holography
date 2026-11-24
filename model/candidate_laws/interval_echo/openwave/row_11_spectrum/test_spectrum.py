"""Independent microscopic, variational, and meter checks for the charged minima."""

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation

HERE = Path(__file__).parent


def load(name):
    return json.loads((HERE / name).read_text())


def test_direct_native_two_vacancy_one_double():
    law = candidate.native_law
    trans, cost = law.local_term_parts()
    block = (trans + 8 * cost)[: law.N_SRC, : law.N_SRC]
    words = []
    for a in range(12):
        for b in range(a):
            for x in set(range(12)) - {a, b}:
                word = [1] * 12
                word[a] = word[b] = 0
                word[x] = 2
                words.append(tuple(word))
    ix = {w: i for i, w in enumerate(words)}
    h = np.zeros((len(ix), len(ix)))
    for i, w in enumerate(words):
        h[i, i] += 8 * sum(n * n - 1 for n in w)
        for a, b in cell().seams:
            x, y = w[a], w[b]
            j = law.SRC.index((x, y))
            h[i, i] += block[j, j]
            if x != y:
                dest = list(w)
                dest[a], dest[b] = y, x
                h[ix[tuple(dest)], i] += block[law.SRC.index((y, x)), j]
    np.testing.assert_allclose(h, h.T)
    v = eigh(h, subset_by_index=(0, 0), eigvals_only=True)[0]
    rec = load("results.json")["cell"]["8"]["compositions"][1]
    assert (
        abs(v - rec["min_energy"]) < 1e-10
        and abs(v - rec["formation_read"]["energy"]) < 1e-6
    )


def test_all_dressings_and_held_density():
    s = load("results.json")
    for n in ("belt", "cell"):
        k = s[n]["8"]
        assert sum(x["size"] for x in k["compositions"]) == (
            8350 if n == "belt" else 69576
        )
        assert max(x["offset_free"] for x in k["formation_triples"]) < 5
    p = load("transport.json")
    for c in p["cases"]:
        delta = s["belt"][str(c["k"])]["compositions"][c["p"]]["recoils"][0]["read"][
            "energy"
        ]
        for r in c["panels"]:
            count = np.array(r["port_012_and_failure"])
            assert np.all(count.sum(axis=1) == p["shots"])
            signed = (count[:, 2] - count[:, 0]) / count[:, :3].sum(axis=1)
            z = np.fft.fft(signed)[1]
            if r["time"] == 0:
                first = z
            assert (
                abs(z - first * np.exp(-1j * r["time"] * delta))
                <= r["simultaneous_error_bound"]
            )


def test_variational_word_envelope_and_owner_bound():
    record = load("bounds.json")
    v = record["arbitrary_three_gateways"]
    assert v["minimum_defects_for_ratio"] == {"207": 139, "3477": 2319}
    assert max(r["mass_ratio_ceiling"] for r in record["boards"]) < 760
    # Independent level-order traversal follows the physical federation adjacency.
    for level, r in enumerate(record["boards"]):
        b = federation(level)
        fringe = [0]
        layers = []
        seen = {0}
        while fringe:
            layers.extend(sorted(fringe))
            nxt = {
                q // 12
                for owner in fringe
                for i in range(12 * owner, 12 * owner + 12)
                for q in b.neighbours[i]
            } - seen
            seen.update(nxt)
            fringe = nxt
        order = [12 * i + j for i in layers for j in range(12)]
        edges = np.array(b.seams)
        for p in (0, b.n // 4, b.n // 2 - 1):
            n = np.ones(b.n, dtype=int)
            n[order[: p + 1]] = 0
            if p:
                n[order[-p:]] = 2
            a, c = n[edges[:, 0]], n[edges[:, 1]]
            twice_energy = 16 * (
                sum((a - c) ** 2) + np.dot(6 - b.degree, n * n - 1)
            ) + sum(a != c)
            assert twice_energy == r["twice_sharp_energies"][p]
        assert r["lowest_floor"] == 32 and r["minimum_defect_count"]["3477"] is None


def test_reported_inertia_and_cell_energy_span():
    d = load("results.json")
    belt = d["belt"]["8"]["compositions"]
    for m, expected in enumerate(((168.56, 3209.78), (164.91, 3301.44))):
        gap = [s["recoils"][m]["read"]["energy"] for s in belt]
        np.testing.assert_allclose(
            [gap[1] / gap[j] for j in (3, 4)], expected, atol=0.02
        )
    cell = d["cell"]["8"]["compositions"]
    assert (
        max(x["spectral_ceiling"] for x in cell)
        / min(x["spectral_floor"] for x in cell)
        < 22
    )


def test_owner_inequality_by_cuts_for_every_gateway_choice():
    # Count unlike-colour edges by two binary cuts. A vacancy/double edge costs four.
    s = np.indices((3,) * 12, dtype=np.int8).reshape(12, -1).T - 1
    v, w = s < 0, s > 0
    z, g = v.sum(axis=1), w.sum(axis=1)
    cuts = np.zeros(len(s), dtype=int)
    for i, j in cell().seams:
        cuts += v[:, i] != v[:, j]
        cuts += w[:, i] != w[:, j]
        cuts += 2 * ((v[:, i] & w[:, j]) | (v[:, j] & w[:, i]))
    correction = np.sort(np.where(v, -1, np.where(w, 3, 0)), axis=1)[:, :9].sum(axis=1)
    c = cuts + correction
    assert np.min(4 * c - 6 * (g - z) - 3 * (g + z)) == 0
    # Independent set-partition recursion permits more budget than any physical board.
    mins = {
        (i, j): int(min(c[(z == i) & (g == j)])) for i in range(5) for j in range(4)
    }
    best = {(0, 0): 0}
    for d in range(1, 8):
        for i in range(min(d, 4) + 1):
            j = d - i
            if j > 3:
                continue
            best[i, j] = min(
                u + best[i - a, j - b]
                for (a, b), u in mins.items()
                if a <= i and b <= j and a + b
            )
    assert [best[j + 1, j] for j in range(4)] == [4, 11, 14, 15]
