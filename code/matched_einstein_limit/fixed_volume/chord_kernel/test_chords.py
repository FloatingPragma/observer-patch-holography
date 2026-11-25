import hashlib
import json
import math

import numpy as np
import pytest

from . import measure


@pytest.fixture(scope="module")
def outcome():
    return measure.run()


def test_first_exit_flow_and_flat_roots():
    for r in (0.0, 0.3, 0.8):
        for mu in (-1.0, -0.7, 0.0, 0.5, 1.0):
            tick = 1 / 4096
            first = measure.clock_branch(1, r, mu, 0, tick)
            root = -r * mu + math.sqrt(1 - r * r + (r * mu) ** 2)
            assert -0.0000000001 <= first * tick - root <= tick + 1e-13
    for k in (-1.0, 1.0):
        ell, r, tick = 0.5, 0.25, 1e-8
        f, inv = (math.cosh, math.acosh) if k < 0 else (math.cos, math.acos)
        independent = inv(f(ell) / f(r))
        first = measure.clock_branch(ell, r, 0, k, tick)
        assert 0 <= first * tick - independent <= tick


def test_native_receivers_and_aperture(outcome):
    a = outcome["native_chart"]
    period = a["hand_periods_per_mesh_edge"]
    assert abs(a["native_hand_single_period"] - math.pi) < 1e-12
    for row in a["straight"]:
        site, axis = row["site"], row["axis"]
        b = math.ceil(
            math.sqrt(3.5**2 - sum(x * x for i, x in enumerate(site) if i != axis))
        )
        assert row["legs"] == [b - site[axis], b + site[axis]]
        for leg, ar in zip(row["legs"], row["trials"]):
            assert ar["hand_reads"] == leg * period
            assert all(x[1] for x in ar["trace"])
    assert all(x["hand_reads"] == 5 * period for x in a["turn_after_first"])
    assert all(x["hand_reads"] == 2 * period for x in a["anisotropic_aperture"])
    assert not a["contrast_only_data"]["trace"][-1][1]


def test_refinement_chart_controls_and_signed_collar(outcome):
    flat = outcome["flat_mesh"]
    assert all(
        0 < b["mass_relative_error"] < a["mass_relative_error"]
        for a, b in zip(flat, flat[1:])
    )
    for row in flat:
        assert row["maximum_weight_error_per_radius"] < 1 / row["points_per_radius"]
    np.testing.assert_allclose(
        flat[-1]["integral_per_radius_fourth"][1:], 4 * np.pi / 105, atol=0.006
    )
    for prop, end_bound in (
        ("maximum_product_relative_error", 0.00074),
        ("maximum_sectional_error", 0.001),
        ("maximum_contraction_area_disagreement", 0.00205),
    ):
        err = [max(c[prop] for c in p["cases"]) for p in outcome["comparison_charts"]]
        assert (
            all(0 < b < 0.28 * a for a, b in zip(err, err[1:])) and err[-1] < end_bound
        )
    c = outcome["formed_collar"]
    assert (
        len(c["completion"]) == 12
        and 0.0305 < c["maximum_convex_transfer_fraction"] < 0.0306
    )
    archive = json.loads((measure.HERE / "results.json").read_text())
    archived_chart = archive["native_chart"]
    digest = hashlib.sha256(json.dumps(archived_chart, sort_keys=True, separators=(",", ":")).encode())
    assert digest.hexdigest() == "b4e4a3c5c8f3d122059625e668c267560b81dc85e4ca3726cba6a873e97e2015"
    for name, value in archived_chart.items():
        if name == "native_hand_single_period":
            # Two independent numerical clock calibrations of the exact pi period.
            for period in (value, outcome["native_chart"][name]):
                assert abs(period - math.pi) <= 1e-13
        else:
            assert outcome["native_chart"][name] == value
    np.testing.assert_allclose(
        c["maximum_convex_transfer_fraction"],
        archive["formed_collar"]["maximum_convex_transfer_fraction"],
        rtol=1e-8,
    )


def test_anisotropic_product_planes_and_separate_area(outcome):
    panels = outcome["product_charts"]
    assert len(panels) == 8
    for panel in panels:
        k = panel["curved_factor"]
        ell = panel["radius"]
        for frame in panel["frames"]:
            ratios = (
                np.array([1, 0, 0])
                if frame["frame"] == "aligned"
                else np.array([4, 36, 9]) / 49
            )
            np.testing.assert_allclose(
                [p["expected_sectional"] for p in frame["planes"]],
                k * ratios,
                atol=1e-13,
            )
            assert (
                abs(frame["clock_contraction"] - panel["equal_volume_area_contraction"])
                < 0.13 * ell**2
            )
            for ray in frame["planes"]:
                assert (
                    abs(ray["expected_sectional"] - ray["inferred_sectional"])
                    < 0.17 * ell**2
                )
