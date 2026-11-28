"""Check the coupled material and geometric inferences from deposited measurements."""

from fractions import Fraction as F
from itertools import combinations, permutations
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / (name + ".json")).read_bytes())


def test_interval_inverse_from_meter_standards():
    matched = read("six_meter_pointer")
    source = read("charge_selected_flight")["windows"]
    ordinary = read("ordinary_service_distribution")["windows"]
    times, labels = matched["times"], matched["labels"]
    k = sp.Rational(*matched["common_ratio"])

    def matrix(side):
        return sp.Matrix(
            [
                [
                    sp.Rational(z[side])
                    for t in times
                    for z in source[t]["meter_quadratures"][n]
                ]
                for n in labels
            ]
        )

    def target(side):
        return sp.Matrix(
            [
                sum(sp.Rational(ordinary[t]["held_ordinary"][n][side]) for t in times)
                for n in labels
            ]
        )

    low, high = matrix(0), matrix(1)
    left, right = target(0), target(1)
    inverse = ((low + high) / 2).inv()
    residual = inverse.applyfunc(abs) * ((high - low) / 2)
    resolvent = (sp.eye(6) - residual).inv()
    assert all(v >= 0 for v in resolvent)
    norm_weights = resolvent * sp.ones(6, 1)
    assert all(v >= 1 for v in norm_weights)
    norm = max(
        v / w for v, w in zip(residual * norm_weights, norm_weights, strict=True)
    )
    assert norm <= sp.Rational(matched["preconditioned_contraction_upper"]) < 1
    center = inverse * ((left + right) / 2)
    radius = resolvent * (
        inverse.applyfunc(abs) * ((right - left) / 2) + residual * center.applyfunc(abs)
    )
    for j, t in enumerate(times):
        box = matched["exact_phase_enclosures"][t]
        for axis, (lo, hi) in enumerate(box):
            i = 2 * j + axis
            assert (
                sp.Rational(lo)
                <= k * (center[i] - radius[i])
                <= k * (center[i] + radius[i])
                <= sp.Rational(hi)
            )
            assert lo <= matched["central_settings"][t][axis] <= hi
        assert sum(max(F(a) ** 2, F(b) ** 2) for a, b in box) < F(9, 16)
    for meter in matched["material_at_center"].values():
        assert all(sum(tape) == 2**72 for tape in meter["counts"].values())
        lo, hi = meter["observed"]
        a, b = meter["forecast_at_center"]
        assert lo < 0 < hi and max(lo, a) <= min(hi, b)


def test_shape_and_fixed_volume_areas_from_raw_clock_counts():
    clocks = read("six_meter_clocks")
    shape = read("pointer_shape")
    rates = read("conical_ball_radii")
    starts = clocks["starts"]
    routes = list(permutations(range(6), 2))
    arms = {}
    for arm, table in clocks["clock_counts"].items():
        assert len(table) == len(routes)
        assert all(
            sum(tape) == starts and min(tape) >= 0 and len(tape) == 11 for tape in table
        )
        arms[arm] = [
            F(sum(n * tick for tick, n in enumerate(tape, 2)), 4 * starts)
            for tape in table
        ]
    beta = F(*rates["ordinary_fraction"])
    base, change = [], []
    for r, d, h, o in zip(
        *(arms[a] for a in ("reference", "dense", "held", "ordinary")), strict=True
    ):
        base.append((255 * r + d + beta * (o - d)) / 256)
        change.append(((h - d) - beta * (o - d)) / 256)
    for row in shape["ratio_matrix"]:
        for sender in range(6):
            assert (
                sum(F(v) for (i, _), v in zip(routes, row, strict=True) if i == sender)
                == 0
            )

    def inside(interval, value):
        assert interval[0] <= value <= interval[1]

    with mp.workdps(46):

        def geometry(h):
            means = [a + h * b for a, b in zip(base, change, strict=True)]
            u = {
                route: mp.log(mp.mpf(v.numerator) / v.denominator)
                for route, v in zip(routes, means, strict=True)
            }

            def skew(i, j):
                return u[i, j] - u[j, i] if i != j else 0

            gains = [sum(skew(i, j) for j in range(4)) / 4 for i in range(4)]
            gains.extend(
                sum(skew(i, j) + gains[j] for j in range(4)) / 4 for i in (4, 5)
            )
            edge = {
                (i, j): mp.exp((u[i, j] + u[j, i] - gains[i] - gains[j]) / 2)
                for i, j in combinations(range(6), 2)
            }

            def length(i, j):
                return edge[tuple(sorted([i, j]))] if i != j else 0

            angles = {}
            for tet in combinations(range(6), 4):
                for i in tet:
                    rest = [j for j in tet if j != i]
                    g = mp.matrix(
                        [
                            [
                                (
                                    length(i, j) ** 2
                                    + length(i, k) ** 2
                                    - length(j, k) ** 2
                                )
                                / 2
                                for k in rest
                            ]
                            for j in rest
                        ]
                    )
                    den = mp.sqrt(g[0, 0] * g[1, 1] * g[2, 2]) + sum(
                        g[j, (j + 1) % 3] * mp.sqrt(g[(j + 2) % 3, (j + 2) % 3])
                        for j in range(3)
                    )
                    angles[tet, i] = 2 * mp.atan2(mp.sqrt(mp.det(g)), den)
            balls = {}
            for row in rates["balls"]:
                tri, i = tuple(row["triangle"]), row["vertex"]
                balls[tri, i] = sum(
                    v
                    for (tet, vertex), v in angles.items()
                    if vertex == i and len(set(tri) & set(tet)) == 2
                )
            return edge, balls

        before, a = geometry(0)
        after, b = geometry(1)
        for e in combinations(range(6), 2):
            if e != (0, 1):
                value = mp.log((after[e] / after[0, 1]) / (before[e] / before[0, 1]))
                inside(clocks["log_edge_ratio_differences"][str(e)], value)
        h = F(1, 1024)
        _, plus = geometry(h)
        _, minus = geometry(-h)
        for rate, finite in zip(
            rates["balls"], clocks["proper_volume_area_differences"], strict=True
        ):
            key = (tuple(rate["triangle"]), rate["vertex"])
            assert key == (tuple(finite["triangle"]), finite["vertex"])
            area = a[key] / 4 * (mp.root(b[key] / a[key], 3) - 1)
            inside(finite["difference"], area)
            tangent = (plus[key] - minus[key]) * 512 / 12
            for j, interval in enumerate(rate["fixed_volume_area_rates"]):
                inside(interval, tangent / (4**j))
            inside(rate["affine_source_contrast"], 3 * tangent / 8)
            inside(rate["half_minus_full_area_over_radius_fourth"], 48 * tangent)
            assert 0 < rate["minimum_maximum_residual_lower"] <= abs(tangent) / 16
