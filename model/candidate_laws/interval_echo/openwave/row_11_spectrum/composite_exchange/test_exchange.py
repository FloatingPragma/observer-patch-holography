"""Checks using explicitly labeled material centers and a direct contact propagator."""

from collections import Counter
import json
from math import sqrt
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix

from . import measure as ex

DATA = Path(__file__).with_name("results.json")


def tagged(words, copies, sign):
    """Full named-center cover then symmetric/antisymmetric compression."""
    lookup = {}
    tagged_words = []
    entries = []
    for j, word in enumerate(words):
        w = tuple(3 if x == 2 and i != word.index(2) else x for i, x in enumerate(word))
        reps = (w, tuple({2: 3, 3: 2}.get(x, x) for x in w)) if copies == 2 else (w,)
        for k, rep in enumerate(reps):
            lookup[rep] = len(tagged_words)
            tagged_words.append(rep)
            entries.append((len(tagged_words) - 1, j, sign**k / sqrt(len(reps))))
    rows, cols, vals = [], [], []
    for col, w in enumerate(tagged_words):
        diag = 0
        for j in range(ex.N):
            k = (j + 1) % ex.N
            x, y = min(w[j], 2), min(w[k], 2)
            diag += 8 * (x - y) ** 2
            if x == y:
                continue
            diag += 0.5
            other = list(w)
            other[j], other[k] = other[k], other[j]
            rows.append(lookup[tuple(other)])
            cols.append(col)
            vals.append(-0.5)
        rows.append(col)
        cols.append(col)
        vals.append(diag)
    h = coo_matrix((vals, (rows, cols)), shape=(len(tagged_words),) * 2).tocsr()
    u, v, a = zip(*entries)
    q = coo_matrix((a, (u, v)), shape=(len(tagged_words), len(words))).tocsr()
    block = (q.T @ h @ q).toarray()
    np.testing.assert_allclose(h @ q.toarray(), q @ block, atol=2e-12)
    return block


def test_winding_matches_every_species_and_edge():
    receipts = json.loads(DATA.read_text())
    pair = []
    for copies in (1, 2):
        steps = ex.exchange(copies, "native")
        word = list(ex.source(copies))
        crossings = Counter()
        for macro in steps["macrostep_edges"]:
            for j in macro:
                k = (j + 1) % ex.N
                assert word[k] == 1
                crossings[word[j], j] += 1
                word[j], word[k] = word[k], word[j]
        pair.append(crossings)
        expected = [[f"{(p-1)%copies}:{i}" for i in range(3)] for p in range(copies)]
        assert steps["final_tags"] == expected
    assert (
        pair[0]
        == pair[1]
        == Counter({(s, j): 1 if s == 2 else 2 for j in range(ex.N) for s in (0, 2)})
    )
    for regime in ex.REGIMES:
        phase = [ex.exchange(c, regime)["phase"] for c in (1, 2)]
        assert phase[0] * phase[1] == (-1 if "fermion" in regime else 1)
    odd = receipts["odd_internal_control"]
    assert odd["exchange"]["dark"] == 0
    assert odd["exchange_and_first_reflection"]["bright"] == 0
    for row in odd.values():
        assert row["dark"] + row["bright"] + row["preparation_refused"] == ex.SHOTS
        assert abs(row["preparation_refused"] / ex.SHOTS - 0.75) < 0.003


def test_tagged_center_cover_and_independent_free_probabilities():
    data = json.loads(DATA.read_text())
    local, _ = ex.seam_data(8)
    for copies in (1, 2):
        words = ex.basis(copies)
        run = data["runs"][str(copies)]
        assert ["".join(map(str, w)) for w in words] == run["port_words"]
        assert len(words) == (168 if copies == 1 else 420)
        for reg, sign in (("native", 1), ("center_fermion", -1)):
            h = tagged(words, copies, sign)
            np.testing.assert_allclose(
                h, ex.generator(words, copies, reg, local), atol=2e-12
            )
            # Independently propagate a complete sharp source without diagonalization.
            # Recursively double a short Taylor/Padé interval to 65,536.
            u = expm(-1j * (h - np.eye(len(h)) * h[0, 0]) / 16)
            for _ in range(20):
                u = u @ u
            src = ex.source(copies)
            p = abs(u[:, words.index(src)]) ** 2
            assert abs(sum(p) - 1) < 2e-8
            bins = [words.index(w) for w in {src[-j:] + src[:-j] for j in range(ex.N)}]
            bins.remove(words.index(src))
            obs = next(r for r in run["regimes"][reg]["free"] if r["time"] == 65536)
            assert abs(obs["intact_away"] / ex.SHOTS - sum(p[bins])) < 0.002
        for regime in run["regimes"].values():
            for row in regime["free"]:
                assert sum(row["port_word_counts"]) == ex.SHOTS
                assert row["intact"] == row["source"] + row["intact_away"]
    for reg, sign in (("native", 1), ("center_fermion", -1)):
        x = data["runs"]["2"]["regimes"][reg]
        r = next(v for v in x["free"] if v["time"] == 65536)
        assert r["intact"] > 0.98 * ex.SHOTS
        assert (r["intact_away"] > ex.SHOTS / 2) == (sign == 1)
