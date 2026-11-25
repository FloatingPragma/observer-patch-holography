from fractions import Fraction
from itertools import product
import json
from math import comb, exp, log, sqrt

import numpy as np
import sympy as sp

from .measure import HERE, SHOTS, count_law, expected_fraction, log_mgf, log_rate


def test_exact_lumping_for_distinct_arrangements():
    p0, p1 = Fraction(1, 3), Fraction(2, 3)
    for source in product((0, 1), repeat=4):
        count = [Fraction()] * 5
        for target in product((0, 1), repeat=4):
            weight = Fraction(1)
            for x, y in zip(source, target, strict=True):
                q = p1 if x else p0
                weight *= q if y else 1 - q
            count[sum(target)] += weight
        assert count == count_law(4, sum(source), p0, p1)
        assert sum(count) == 1


def test_declared_curvature_gauge_and_mode():
    x, y, z, a, s, c = sp.symbols("x y z a s c", real=True)
    bilinear = c + s * (x + y - 2 * x * y)
    lift = bilinear + a * y * (y - 1) / 2
    for endpoint in (0, 1):
        two_steps = lift.subs({x: endpoint, y: z}) + lift.subs({x: z, y: endpoint})
        assert sp.solve(sp.diff(two_steps, z).subs(z, endpoint), a) == [4 * s]
        assert (
            sp.expand(two_steps.subs(a, 4 * s) - 2 * c - 2 * s * (z - endpoint) ** 2)
            == 0
        )
    gauge = s * (y * (y - 1) - x * (x - 1))
    assert sp.expand(lift.subs(a, 4 * s) - (c + s * (y - x) ** 2) - gauge) == 0


def test_measured_dual_and_count_refinement():
    p0, p1 = 1 / 3, 2 / 3
    for x in (0, 0.25, 0.5, 1):
        for p in (-7, -1, 0, 3, 10):
            y = expected_fraction(x, p, p0, p1)
            np.testing.assert_allclose(
                log_rate(x, y, p0, p1) + log_mgf(x, p, p0, p1), p * y, atol=3e-14
            )
        assert abs(log_rate(x, expected_fraction(x, 0, p0, p1), p0, p1)) < 1e-14
        for y in (0, 0.25, 0.5, 0.75, 1):
            # Supermultiplicativity at doubled sizes and the exponential tilt upper bound.
            last = 1 + max(-log(p0), -log(1 - p1))
            for n in (4, 8, 16, 32):
                exact = count_law(n, int(n * x), Fraction(1, 3), Fraction(2, 3))[
                    int(n * y)
                ]
                finite = -log(exact) / n
                assert log_rate(x, y, p0, p1) <= finite + 1e-13 < last + 2e-13
                last = finite
    assert comb(64, 16) ** 4 > comb(64, 32) ** 3
    assert 32 > 27


def test_every_endpoint_and_history():
    result = json.loads((HERE / "results.json").read_text())
    assert result["shots_per_bank"] == SHOTS
    panels = [g for menu in result["menus"].values() for g in menu["fraction_banks"]]
    # Simultaneous 99% DKW across all count banks and the four independent prefix banks.
    threshold = sqrt(log(2 * (len(panels) + 4) / 0.01) / (2 * SHOTS))
    for bank in panels:
        assert len(bank["counts"]) == bank["copies"] + 1
        assert sum(bank["counts"]) == SHOTS
        assert bank["cdf_distance"] <= threshold
    assert result["menus"]["symmetric"]["initial_one"] == "1/2"
    assert result["menus"]["symmetric"]["kernel"] == [["2/3", "1/3"], ["1/3", "2/3"]]
    assert result["imported_chain"]["initial_one"] == "54356/61511"
    assert result["imported_chain"]["kernel"] == [
        ["1324/1431", "107/1431"],
        ["5/508", "503/508"],
    ]
    for menu in [*result["menus"].values(), result["imported_chain"]]:
        chain = menu["histories"]
        assert len(chain["counts_lexicographic"]) == 32
        assert sum(chain["counts_lexicographic"]) == SHOTS
        reference = np.cumsum([float(Fraction(p)) for p in chain["expected"]])
        assert (
            np.max(abs(np.cumsum(chain["counts_lexicographic"]) / SHOTS - reference))
            < threshold
        )
    held = result["held_test"]
    lo, hi = held["simultaneous_99"]
    assert 0 < lo < held["exact_finite"] < hi
    assert abs(exp(4 * held["limiting"]) - 32 / 27) < 5e-12
    common = result["shared_draw_control"]["counts"]
    assert sum(common) == SHOTS and sum(common[1:-1]) == 0


def test_fraction_action_symmetries():
    # Reflection is an exact relabeling symmetry of every finite symmetric count bank.
    n = 8
    p0, p1 = Fraction(1, 3), Fraction(2, 3)
    for m in range(n + 1):
        assert count_law(n, m, p0, p1) == count_law(n, n - m, p0, p1)[::-1]
    # On the edge fixing zero, even a nonidentity increasing transformation that
    # fixes the three mandatory points 0, p0 and 1 is detected by the cost.
    q = float(p0)
    y = q / 2
    changed = y * y / q
    assert (
        log_rate(0, changed, q, float(p1))
        > log_rate(0, y, q, float(p1))
        > log_rate(0, q, q, float(p1))
    )
    assert (
        log_rate(0, 1, q, float(p1))
        > log_rate(0, (1 + q) / 2, q, float(p1))
        > log_rate(0, q, q, float(p1))
    )
    # The asymmetric historical schedule distinguishes a rail relabeling.
    assert (
        count_law(1, 0, Fraction(107, 1431), Fraction(503, 508))
        != count_law(1, 1, Fraction(107, 1431), Fraction(503, 508))[::-1]
    )
