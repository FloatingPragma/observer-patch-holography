"""Independent capacity, leading-amplitude, held-route and inventory checks."""

from collections import Counter
import gzip
import hashlib
import json
from math import comb

from . import measure as m


def load():
    return json.loads(gzip.decompress((m.HERE / "records.json.gz").read_bytes()))


def test_inventory_support_formula_against_independent_local_population_enumeration():
    possible = {(0, 0)}
    for c in range(1, 7):
        possible = {
            (stock + number, bits | ((number % 2) << (c - 1)))
            for stock, bits in possible
            for number in range(13)
        }
        for k in range(12 * c + 1):
            observed = {word for stock, word in possible if stock == k}
            assert observed == {
                q for q in range(1 << c) if m.admissible(c, k, q.bit_count())
            }
            assert len(observed) == m.word_count(c, k) == m.word_count(c, 12 * c - k)
    for c in (20, 80, 320):
        assert m.word_count(c, 1) == c and m.word_count(c, 2) == 1 + comb(c, 2)


def test_every_sharp_target_has_a_nonzero_signed_leading_coefficient():
    data = load()
    for sec in data["enrolled"]:
        k = sec["inventory"]
        ports, states, start, adj, h = m.sector(data["config"], k)
        assert len(ports) == 72 and len(states) == comb(72, k)
        assert [list(x) for x in states] == sec["basis"]
        assert adj.nnz == 2 * sec["undirected_exchanges"]
        d = sec["leading"]["distances"]
        walks = sec["leading"]["shortest_walks"]
        assert (d[start], walks[start]) == (0, 1)
        assert Counter(d)[0] == 1
        for i in range(len(states)):
            adjacent = adj.indices[adj.indptr[i] : adj.indptr[i + 1]]
            assert all(abs(d[j] - d[i]) <= 1 for j in adjacent)
            if d[i]:
                assert (
                    walks[i] == sum(walks[j] for j in adjacent if d[j] == d[i] - 1) > 0
                )
            assert h[i, i] == 8.5 * len(adjacent)
            assert all(h[i, j] == -0.5 for j in adjacent)


def test_all_raw_simulated_outcomes_and_two_independent_meter_tapes():
    r = load()
    for path, sha in r["sources"].items():
        assert (
            hashlib.sha256((m.HERE.parents[4] / path).read_bytes()).hexdigest() == sha
        )
    assert sorted(r["config"]["auditors"]) == [18, 19]
    assert set(r["config"]["auditors"]).isdisjoint(r["config"]["owners"])
    for sec in r["enrolled"]:
        k = sec["inventory"]
        recon = []
        for state in sec["basis"]:
            q = 0
            for j in state:
                q ^= 1 << (j // 12)
            recon.append(q)
        assert recon == sec["words"]
        assert (
            set(recon)
            == set(sec["allowed"])
            == {x["word"] for x in sec["material_readbacks"]}
        )
        for held in sec["material_readbacks"]:
            raw = held["raw"]
            assert raw[0] == raw[1] and len(raw[0]) == 73
            assert all(port in (0, 1) for port in raw[0]) and raw[0][-1] == 1
            occupied = [j for j, val in enumerate(raw[1][:-1]) if val == 0]
            assert occupied == held["state"]
            assert (
                sum((sum(i // 12 == c for i in occupied) % 2) << c for c in range(6))
                == held["word"]
            )
        for mode, panels in sec["panels"].items():
            for row in panels:
                histogram = Counter()
                for word, n in zip(recon, row["counts"], strict=True):
                    histogram[word] += n
                assert row["words"] == [histogram[i] for i in range(64)]
                assert sum(row["words"]) == r["shots"]
                assert row["observed"] == sum(n > 0 for n in row["words"])
                assert row["largest_residual"] < r["radius"]
                if mode == "diagonal":
                    assert row["words"][recon[sec["start"]]] == r["shots"]
        assert [p["observed"] for p in sec["panels"]["native"]] == (
            [1, 3, 4, 5] if k == 1 else [1, 2, 4, 5]
        )
