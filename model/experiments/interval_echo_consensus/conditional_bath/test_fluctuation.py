from fractions import Fraction as Fr
from itertools import product
import json
from math import log

import numpy as np

from .fluctuation import EPS, SHOTS, encode, forecast, make_repair, visible_repair
from .gibbs import exact_and_native
from .measure import HERE, WORDS


def test_visible_lumping_on_every_native_input_and_every_address():
    for j in range(4):
        k, den, physical = exact_and_native(make_repair(j))
        for i, w in enumerate(WORDS):
            x = w & 15
            target = dict(visible_repair(x, j))
            for y in range(16):
                cols = [h for h, v in enumerate(WORDS) if v & 15 == y]
                assert Fr(sum(k[i, h] for h in cols), den) == target.get(y, 0)
                assert abs(physical[i, cols].sum() - float(target.get(y, 0))) < 2e-13


def test_all_record_odds_their_integral_and_entropy_loss():
    r, p, q, f, b, odds = forecast()
    for j, x, y in product(range(4), range(16), range(16)):
        assert r[x] * dict(visible_repair(x, j)).get(y, 0) == r[y] * dict(
            visible_repair(y, j)
        ).get(x, 0)
    assert len(f) == 64 and set(f) == set(b)
    assert sum(f.values()) == sum(b.values()) == sum(q) == sum(p) == sum(r) == 1
    assert min(q) > 0
    for path in product(range(16), repeat=3):
        x, m, y = path
        assert f.get(path, 0) * q[y] * r[x] == b.get(path, 0) * p[x] * r[y]
    assert sum(v / odds[k] for k, v in f.items()) == 1

    def divergence(a):
        return sum(float(s) * log(s / ref) for s, ref in zip(a, r))

    mean = sum(float(v) * log(odds[k]) for k, v in f.items())
    assert abs(mean - (divergence(p) - divergence(q))) < 1e-13
    assert min(odds.values()) < 1 < max(odds.values()) and mean > 0


def test_observed_paths_and_calibrations_with_registered_controls():
    receipt = json.loads((HERE / "fluctuation.json").read_text())
    hist = {key: dict(rows) for key, rows in receipt["histograms"].items()}
    assert all(sum(row.values()) == SHOTS for row in hist.values())
    r, p, q, f, b, odds = forecast()
    for lab, target in (("R", r), ("P", p), ("Q", q)):
        assert (
            max(abs(hist[lab].get(x, 0) / SHOTS - v) for x, v in enumerate(target))
            < EPS
        )
    forward, reverse = hist["forward"], hist["reverse"]
    assert set(forward) == {encode(*k) for k in f}
    assert set(reverse) == {encode(*k[::-1]) for k in f}
    for k in f:
        v = float(odds[k])
        assert abs(forward[encode(*k)] / SHOTS - float(f[k])) < EPS
        assert abs(reverse[encode(*k[::-1])] / SHOTS - float(b[k])) < EPS
        assert (
            abs(forward[encode(*k)] - v * reverse[encode(*k[::-1])])
            < EPS * (1 + v) * SHOTS
        )
    for v, g in receipt["statistics"]["by_odds"].items():
        v = float(Fr(v))
        assert abs(g["forward_count"] - v * g["reverse_count"]) < EPS * (1 + v) * SHOTS
    # The inventory reference is stationary on each visible conditional fibre.
    for x, m, y in f:
        pred_eq = r[x] * dict(visible_repair(x, 0))[m] * dict(visible_repair(m, 2))[y]
        assert (
            abs(hist["equilibrium"].get(encode(x, m, y), 0) / SHOTS - float(pred_eq))
            < EPS
        )
    assert all(c // 256 == c // 16 % 16 == c % 16 for c in hist["parked"])
    dif = [v for k, v in f.items() if len(set(k)) > 1]
    assert sum(dif) > EPS
    stats = receipt["statistics"]
    assert stats["reverse_wrong_order_outside_support"] / SHOTS > EPS
    assert stats["flat_reference_residual"] > stats["flat_reference_radius"]
    assert stats["independent_odds_intervals_overlap"] == 64
    assert stats["negative_from_independent_calibration"] > 0
    assert (
        abs(stats["negative_count"] / SHOTS - float(Fr(stats["negative_forecast"])))
        < EPS
    )
    multipliers = np.array(list(map(float, odds.values())))
    assert abs(stats["exp_minus_sigma"] - 1) < EPS * np.ptp(1 / multipliers)
    assert abs(stats["mean_sigma"] - stats["kl_descent"]) < EPS * np.ptp(
        np.log(multipliers)
    )
