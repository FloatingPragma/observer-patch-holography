"""Independent coordinate/ancestry census and count-only material readbacks."""

from itertools import permutations
import gzip
import json

import numpy as np
import pytest

from model.candidate_laws.interval_echo.openwave.row_29_maxwell.experiment import (
    geometry,
)
from .chart import graph_from_bank
from .five_roots import RADIUS, SHOTS, past_of_hub
from .measure import HERE, READ_CODES, decode


def receipt(file):
    if file.endswith(".gz"):
        with gzip.open(HERE / file, "rt") as src:
            return json.load(src)
    return json.loads((HERE / file).read_text())


def test_full_census_independently():
    """Run all starts; overlapping seams precede one another, then take closure."""
    flights = receipt("flights.json")
    edges = sorted(graph_from_bank(flights["sampled"][0]).edges)
    assert edges == sorted(graph_from_bank(flights["sampled"][1]).edges)
    with pytest.raises(ValueError, match="Contact standard"):
        graph_from_bank(flights["parked"])
    edges = np.asarray(edges)
    points = geometry()[edges].mean(axis=1)
    adjacent = np.any(edges[:, None, :, None] == edges[None, :, None, :], axis=(2, 3))
    dist = ((points[:, None] - points[None, :]) ** 2).sum(axis=2)
    expected = receipt("census.json")["lengths"]
    for n in range(1, 5):
        counts = np.zeros(3, dtype=int)
        for first in range(30):
            w = np.array(
                [
                    (first, *s)
                    for s in permutations([k for k in range(30) if k != first], n - 1)
                ]
            )
            m = len(w)
            past = np.zeros((m, n, n), dtype=bool)
            height = np.zeros((m, n), dtype=int)
            for i in range(n):
                for j in range(i):
                    contact = adjacent[w[:, j], w[:, i]]
                    past[:, i, j] |= contact
                    past[:, i, :] |= contact[:, None] & past[:, j, :]
                    height[:, i] = np.maximum(
                        height[:, i], np.where(contact, 1 + height[:, j], 0)
                    )
            low, high = np.zeros(m), np.full(m, np.inf)
            for i in range(n):
                for j in range(i):
                    dt = height[:, i] - height[:, j]
                    d2 = dist[w[:, i], w[:, j]]
                    assert np.all(d2 > 1e-8)
                    bound = d2 / np.maximum(1, dt**2)
                    low = np.maximum(low, np.where(past[:, i, j], bound, 0))
                    high = np.minimum(
                        high, np.where((dt != 0) & ~past[:, i, j], bound, np.inf)
                    )
            counts += [
                sum(low < high - 1e-10),
                sum(abs(low - high) <= 1e-10),
                sum(low > high + 1e-10),
            ]
        assert counts.tolist() == [
            expected[n - 1][k] for k in ("feasible", "tie", "reversed")
        ]


def test_auditors_and_native_controls():
    rows = receipt("records.json.gz")
    witnesses = receipt("census.json")["witnesses"]
    for name, panel in rows.items():
        expected = list(map(tuple, witnesses[name]["edges"]))
        for scans in panel["ordinary"]:
            word, ps, loads = decode(scans)
            assert word == expected
            assert all(sum(xs) == sum(loads[0]) for xs in loads)
            written = {}
            for k, ((a, b), parents) in enumerate(zip(word, ps, strict=True)):
                assert set(parents) == {written[v] for v in (a, b) if v in written}
                written[a] = written[b] = k
                total = loads[k][a] + loads[k][b]
                assert sorted([loads[k + 1][a], loads[k + 1][b]]) == [
                    total // 2,
                    (total + 1) // 2,
                ]
        _, _, loads = decode(panel["parked"])
        assert all(xs == loads[0] for xs in loads)
        assert any(decode(s)[2][-1] != decode(s)[2][0] for s in panel["ordinary"])
        with pytest.raises(ValueError, match="Incomplete endpoint"):
            decode(panel["shutter"])


def test_matched_local_comparisons():
    rows = receipt("matched.json.gz")

    def endpoint_tags(scans):
        word, _, _ = decode(scans)
        return [
            (scan[0]["stamp"][6 * v : 6 * v + 6], scan[0]["load"][v])
            for i in (0, 3)
            for scan in scans[i : i + 2]
            for v in word[i]
        ]

    common = endpoint_tags(rows["linked"]["empty"][0])
    for name, bank in rows.items():
        for scans in bank["empty"]:
            assert endpoint_tags(scans) == common
        for scans in bank["loaded"]:
            word, parents, load = decode(scans)
            assert parents[3] == (2,) and ((0 in parents[2]) == (name != "separate"))
            assert sum(load[-1][v] for v in word[-1]) == int(name == "linked")


def test_five_root_physical_contrasts():
    d = receipt("five_roots.json.gz")
    assert d["muted_trace"] == d["rows"][-1]["trace"]
    prefix = d["rows"][0]["trace"][:6]
    for mask, row in enumerate(d["rows"]):
        selected = [i for i in range(5) if mask & (1 << i)]
        word, roots, loads = past_of_hub(row["trace"])
        assert roots == selected and row["trace"][:6] == prefix
        assert all(xs == [1] * 12 for xs in loads)
        assert len(set(word)) == len(word) and len(word) == 6 + len(selected)
        assert word[-1] == tuple(d["hub"])
        for bank in row["banks"]:
            i = bank["root"]
            totals = np.array(
                [sum(READ_CODES.index(v) for v in pair) for pair in bank["bins"]]
            )
            weights = np.array(bank["branches"])
            assert weights.sum() == 2 ** len(word) and set(totals) <= {0, 1}
            expected = 0
            if i in selected:
                (side,) = set(d["spokes"][i]) & set(d["hub"])
                later = sum(side in d["spokes"][j] for j in selected if j > i)
                expected = 2.0 ** (-later)
            assert (totals @ weights) / weights.sum() == expected
            assert sum(bank["sampled"]) == SHOTS
            empirical = (totals @ bank["sampled"]) / SHOTS
            assert abs(empirical - expected) <= RADIUS
            assert (empirical > 2 * RADIUS) == (i in selected)
    assert d["muted"]["bins"] == [[READ_CODES[0]] * 2]
