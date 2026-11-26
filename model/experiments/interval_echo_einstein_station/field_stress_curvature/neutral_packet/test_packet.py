"""Check the separately counted neutral field and charged ensemble forecasts."""

from fractions import Fraction as F
import json

import numpy as np

from . import measure as m, transfer as t


def test_remote_interference_and_vacuum():
    faces, likelihoods, overlap = m.laws()
    omitted = [not {10, 11}.issubset(set(f)) for f in faces]
    assert np.max(abs(overlap[0, omitted])) < 1e-14
    assert np.max(abs(likelihoods[:, :, 4] - likelihoods[:, :, 5])) < 1e-14
    assert (
        np.max(abs(overlap) - np.sqrt(likelihoods[:, :, 6] * likelihoods[:, :, 7]))
        < 1e-13
    )
    summary = m.acquire(False)
    assert summary["total_banks"] == 640
    assert [r["resolved_vertex_disjoint_faces"] for r in summary["by_delay"]] == [
        0,
        0,
        0,
        2,
    ]
    assert summary["by_delay"][-1]["simultaneous_interval"][0] > 0.00249


def test_affine_transfer_and_mislabelled_charge_control():
    report = t.run()
    for held, truth in zip(report["held"], t.WEIGHTS, strict=True):
        lo, hi = held["field_inferred_weight"]
        assert lo < truth < hi and hi - lo < 0.011
        assert held["ordinary_assignment"]["bin_failures"] == 0
        assert held["swapped_assignment"]["bin_failures"] == 43
        assert held["swapped_assignment"]["worst_excluded_probability"] > 0.07260
    assert report["archived_equal_mixture"]["bin_failures"] == 0
    bank = json.loads((t.HERE / "transfer_counts.json").read_text())
    y = bank["neutral_clicks"]
    lo, hi = t.detector_weight(y, 2)
    lo_swapped, hi_swapped = t.detector_weight([y[0], y[1], y[3], y[2]], 3)
    assert (lo, hi) == (lo_swapped, hi_swapped)
    assert 0 < F(1, 4) - lo < F(6, 1000)


def test_one_field_difference_aliases_vacuum():
    control = t.run()["vacuum_with_same_neutral_difference"]
    lo, hi = control["matching_weight"]
    assert 0 < lo < hi < 0.018
    assert control["charged_challenge"]["bin_failures"] == 46
    assert control["charged_challenge"]["worst_excluded_probability"] > 0.06469


def test_two_neutral_readings_recover_three_source_fractions():
    from . import hull

    report = hull.run()
    assert report["coverage_marginals"] == 938
    lo, hi = report["determinant"]
    assert 2.8e-5 < lo <= hi < 2.82e-5
    truth = [(0.75, 0.25, 0), (0.25, 0.75, 0), (0.25, 0.125, 0.625)]
    for row, weights, rejected in zip(report["held"], truth, (47, 46, 44), strict=True):
        for (lower, upper), true in zip(row["weights"], weights, strict=True):
            assert lower <= true <= upper
            assert upper - lower < 0.012
        assert row["charged"]["excluded"] == 0
        assert row["cyclic"]["excluded"] == rejected


def test_complementary_neutral_rates_identify_three_weight_triangle():
    from . import hull

    report = hull.run()
    assert report["determinant"][0] > 2.80476e-5
    assert report["infinity_condition_upper"] < 10.286
    weights = ((0.75, 0.25, 0.0), (0.25, 0.75, 0.0), (0.25, 0.125, 0.625))
    for row, truth, missed in zip(report["held"], weights, (47, 46, 44), strict=True):
        assert all(
            lo <= p <= hi for (lo, hi), p in zip(row["weights"], truth, strict=True)
        )
        assert row["charged"]["excluded"] == 0
        assert row["cyclic"]["excluded"] == missed
    exact = [
        ((F(1), F(1)), (F(0), F(0))),
        ((F(0), F(0)), (F(1), F(1))),
        ((F(0), F(0)), (F(0), F(0))),
    ]
    assert hull.vertices(exact, ((F(1, 4), F(1, 4)), (F(1, 8), F(1, 8)))) == (
        (F(1, 4), F(1, 8), F(5, 8)),
    )
