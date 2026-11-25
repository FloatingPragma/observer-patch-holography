"""Direct four-label propagation and independent checks of the paired readouts."""

import json

import numpy as np
from scipy.sparse import csr_matrix, kron, block_diag, eye
from scipy.sparse.linalg import expm_multiply
from scipy.special import eval_chebyt

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.executable_model import federation
from . import common_time, measure


def direct_matrix(k):
    board = federation(0)
    n = board.n
    h = sum(q * part for q, part in zip([1, k], law.local_term_parts()))
    # Site-major, then sharp (coin, signed contrast) basis. The source is
    # never propagated in a clock-energy eigenbasis in this check.
    tags = [law.LABELS.index((c, s)) for c in (0, 1) for s in (2, -2)]
    blocks = []
    ends = measure.ENDPOINTS
    for pos in range(n):
        loads = tuple(1 + int(pos == v) for v in ends)
        columns = [law.N_SRC * (1 + t) + law.SRC.index(loads) for t in tags]
        blk = h[np.ix_(columns, columns)]
        np.testing.assert_allclose(np.delete(h[:, columns], columns, axis=0), 0, atol=1e-13)
        blocks.append(csr_matrix(blk))
    position = np.diag(3 * k * (6 - board.degree)).astype(float)
    twin = [law.SRC.index(w) for w in ((2, 1), (1, 2))]
    bare = h[np.ix_(twin, twin)]
    for edge in board.seams:
        if tuple(edge) != ends:
            position[np.ix_(edge, edge)] += bare
    return block_diag(blocks, format="csr") + kron(csr_matrix(position), eye(4, format="csr"), format="csr")


def test_full_coin_and_excess_propagation_matches_counts():
    n = federation(0).n
    mask = np.zeros(n, dtype=bool)
    mask[list(measure.ENDPOINTS)] = True
    with np.load(measure.HERE / "counts.npz") as records:
        for k, pos, distance in ((8, 1, 1), (8, 2, 2), (0, 0, 0)):
            matrix = direct_matrix(k)
            for parity, name in ((1, "even"), (-1, "odd")):
                vectors = np.zeros((n, 4, 2))
                for coin in (0, 1):
                    vectors[pos, 2 * coin:2 * coin + 2, coin] = np.array([1, parity]) / np.sqrt(2)
                for appointment in (0, 8, 32, 79):
                    waves = expm_multiply(-1j * measure.TIMES[appointment] * matrix,
                                          vectors.reshape(n * 4, 2)).reshape(n, 2, 2, 2)
                    prob = np.sum(abs(waves) ** 2, axis=2)
                    same = (prob[:, 0, 0] + prob[:, 1, 1]) / 2
                    flip = (prob[:, 1, 0] + prob[:, 0, 1]) / 2
                    bins = [same[~mask].sum(), same[mask].sum(), flip[~mask].sum(), flip[mask].sum()]
                    np.testing.assert_allclose(sum(bins), 1, atol=2e-12)
                    observed = records[f"k{k}_{name}_free{distance}"][appointment] / measure.SHOTS
                    np.testing.assert_allclose(observed, bins, rtol=0, atol=measure.EPS)


def test_replay_and_calibrations():
    summary = json.loads((measure.HERE / "results.json").read_text())
    with np.load(measure.HERE / "counts.npz") as records:
        for key, row in summary["panels"].items():
            np.testing.assert_array_equal(records[key].sum(1), measure.SHOTS)
            if key.startswith("k0_even_") or key.endswith("_park0"):
                assert not np.any(records[key][:, 2:])
            if "_free" in key:
                for branch in ("sym", "anti"):
                    np.testing.assert_array_equal(records[key + "_" + branch].sum(1), measure.SHOTS)
        for entry in summary["source_only_energy_meter"].values():
            degree = entry["degree"]
            np.testing.assert_allclose(entry["formation_energy"], .5 * degree + 8 * (18 - 2 * degree),
                                       rtol=0, atol=degree * measure.EPS)
        assert summary["panels"]["k8_even_free1"]["backaction"]["absolute_difference"] > .04


def test_chebyshev_images_and_exclusions():
    # Every tested p from the exact binomial interval must map inside the
    # reported bounds, including very small p and interior turning points.
    from scipy.stats import binomtest
    for order in (32, 33):
        for hits in (0, 25, 510, 1000):
            a, b = common_time.chebyshev_image([hits], 1000, .001, order)[0]
            ci = binomtest(hits, 1000).proportion_ci(.998)
            values = eval_chebyt(order, 1 - 2 * np.linspace(ci.low, ci.high, 2001))
            assert a - 1e-8 <= min(values) and max(values) <= b + 1e-8
    result = common_time.audit()
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import compare
    compare(result, json.loads((measure.HERE / "common_time.json").read_text()), "paired clocks")
    for name, row in result["settings"].items():
        if name in {"free0", "free1", "free2"}:
            assert row["disjoint_appointments"] > 90
        else:
            assert row["disjoint_appointments"] == 0


def test_exact_common_time_series_including_far_source():
    from fractions import Fraction
    from math import factorial
    from . import series

    # Composition on a blank ideal clock equals cos(32*33*t) coefficient by
    # coefficient. On the moving-source records all five differences appear.
    coeff = [Fraction((-1)**j * 32**(2*j), factorial(2*j)) for j in range(17)]
    exact = [Fraction((-1)**j * 1056**(2*j), factorial(2*j)) for j in range(17)]
    assert series.chebyshev_coefficients(33, coeff) == exact
    for site, order, value in (
        (0, 6, -503360), (1, 6, -513040), (2, 8, Fraction(-8712847, 105)),
        (13, 10, Fraction(-25289, 144)), (156, 32, Fraction(-2058299, 185940969061374140350464000000)),
    ):
        assert series.first_difference(site, order) == (order, value)
