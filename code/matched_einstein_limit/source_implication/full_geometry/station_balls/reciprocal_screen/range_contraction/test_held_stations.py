import json
from itertools import combinations, permutations

import numpy as np
from mpmath import iv

from . import held_stations as h


def test_count_conservation_euclidean_distance_geometry_and_signed_changes():
    record = json.loads(h.RECORD.read_text())
    iv.dps = 35
    raw = record["counts"]
    checked = h.decode(raw)
    for key in checked:
        assert checked[key] == record[key]
    assert record["station_ports"] == [0, 1, 38, 39, 5, 6]
    assert record["directed_routes"] == list(map(list, permutations(range(6), 2)))
    assert record["quartets"] == list(map(list, combinations(range(6), 4)))
    old = json.loads((h.survey.HERE / "results.json").read_text())["acquisitions"]
    for name, batches in raw.items():
        # Independent binary64 reconstruction from complete histograms. The
        # Cayley-Menger determinant gives volumes without choosing an origin.
        logs = np.zeros((6, 6))
        for (i, j), counts in zip(record["directed_routes"], batches, strict=True):
            assert sum(counts) == record["starts_per_route"]
            logs[i, j] = np.log(np.dot(counts, np.arange(2, 13)) / (4 * h.STARTS))
        skew = logs - logs.T
        anchor = skew[:4, :4].mean(axis=1)
        gains = np.r_[anchor, (skew[4:, :4] + anchor).mean(axis=1)]
        sq = np.exp(logs + logs.T - gains[:, None] - gains[None, :])
        np.fill_diagonal(sq, 0)
        center = np.eye(6) - np.ones((6, 6)) / 6
        eigen = np.linalg.eigvalsh(-center @ sq @ center / 2)
        assert abs(eigen[0]) < 1e-12 and eigen[1] > 3
        for q, (lo, hi) in zip(
            record["quartets"],
            record["preparations"][name]["volume_bands"],
            strict=True,
        ):
            cm = np.ones((5, 5))
            cm[0, 0] = 0
            cm[1:, 1:] = sq[np.ix_(q, q)]
            value = np.sqrt(np.linalg.det(cm) / 288)
            assert lo < value < hi
            if q == [0, 1, 2, 3]:
                a, b = old[name]["clocks"]["volume"]
                assert a < value < b
    ordinary = [r["sign"] for r in record["volume_changes"]["ordinary_write"]]
    correlation = [r["sign"] for r in record["volume_changes"]["independent_right"]]
    assert (ordinary.count("shrinks"), ordinary.count("expands")) == (12, 3)
    assert (correlation.count("shrinks"), correlation.count("expands")) == (7, 8)
