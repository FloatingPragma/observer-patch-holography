import json
from itertools import combinations

import numpy as np
import sympy as sp
from mpmath import iv

from . import measure as m


def test_all_matching_responses_and_independent_calibration_cancellation():
    rec = json.loads(m.RECORD.read_text())
    banks = rec["directed_mean_bands"]
    iv.dps = 35
    assert m.compare(banks) == rec["comparisons"]
    edges = list(combinations(range(6), 2))
    rng = np.random.default_rng(413)
    dist = {}
    for q, stocks in banks.items():
        edge_factor = np.exp(rng.normal(size=len(edges)))
        for stock, data in stocks.items():
            directed = {
                tuple(route): sum(ab) / (2 * m.SCALE)
                for route, ab in zip(rec["directed_routes"], data, strict=True)
            }
            gain = np.exp(rng.normal(size=6))
            dist[q, stock] = {
                (i, j): gain[i] * gain[j] * f * np.sqrt(directed[i, j] * directed[j, i])
                for (i, j), f in zip(edges, edge_factor, strict=True)
            }
    constraints = []
    for row in rec["comparisons"]:
        i, j, k, last = row["quartet"]
        a, b = row["matching_pair"]
        pairs = [[(i, j), (k, last)], [(i, k), (j, last)], [(i, last), (j, k)]]
        coef = [int(e in pairs[a]) - int(e in pairs[b]) for e in edges]
        constraints.append(coef)
        quotients = {}
        for label, ds in dist.items():
            top = np.prod([ds[e] for e in pairs[a]])
            bottom = np.prod([ds[e] for e in pairs[b]])
            quotients[label] = top / bottom
        changes = {
            q: quotients[q, "formed"] / quotients[q, "independent"] for q in banks
        }
        gap = changes["excess"] / changes["vacancy"]
        lo, hi = np.array(row["excess_over_vacancy"]) / m.SCALE
        assert lo <= gap <= hi
    equations = sp.Matrix(constraints)
    endpoint = sp.Matrix([[int(v in e) for v in range(6)] for e in edges])
    assert (
        equations.shape == (45, 15) and equations.rank() == 9 and endpoint.rank() == 6
    )
    assert equations * endpoint == sp.zeros(45, 6)


def test_every_excluded_shape_and_opposite_signs():
    rows = json.loads(m.RECORD.read_text())["comparisons"]
    assert len(rows) == 45
    assert all(
        not (x["log_excess_over_vacancy"][0] <= 0 <= x["log_excess_over_vacancy"][1])
        for x in rows
    )
    pairs = [x["log_cross_response"] for x in rows]
    assert (
        sum(
            (x["vacancy"][1] < 0 < x["excess"][0])
            or (x["excess"][1] < 0 < x["vacancy"][0])
            for x in pairs
        )
        == 18
    )
    witness = next(
        x for x in rows if x["quartet"] == [1, 2, 3, 5] and x["matching_pair"] == [0, 2]
    )
    assert witness["log_excess_over_vacancy"][0] > 4 * 252968400
