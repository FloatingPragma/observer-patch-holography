"""Check the common clock, free matter and all-center source comparisons."""

from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def receipt(name):
    z = json.loads((HERE / f"{name}.json").read_bytes())
    for file, digest in z["inputs"].items():
        assert sha256((HERE / file).read_bytes()).hexdigest() == digest
    return z


def overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def test_ranges_and_two_species():
    z = receipt("word_routes")
    code = z["code"]
    assert (
        len(set(tuple(map(tuple, m)) for m in code)) == 10
        and max(z["common_clock_ranges"]) <= 3
    )
    for i, m in enumerate(code):
        assert sorted(v for e in m for v in e) == list(range(12))
        expected = [int(e not in m) for e in z["internal_edges"]]
        assert (
            expected
            == z["sharp_swap_arrivals"]["1"][i]
            == z["sharp_swap_arrivals"]["2"][i]
        )
    shape = np.array(z["ranging_matrix_exact_numerator"])
    assert (
        np.array_equal(shape, shape.T)
        and np.count_nonzero(shape) == 2
        and shape[2, 11] == -12
    )
    assert len(z["tables"]["radar"]) == 190
    for group, stations in z["tables"].items():
        for key, tags in stations.items():
            for data in tags.values():
                for tape in [*data["counts"].values(), *data["reversed_code_counts"]]:
                    assert sum(tape) == z["starts"] and min(tape) >= 0
                for name, metric in data["metrics"].items():
                    assert overlap(metric["held"], metric["from_source"])
                    p = metric["population"]
                    assert metric["held"][0] <= p <= metric["held"][1]
                    r = data["reversed_code_held"][name]
                    assert r[0] <= -p <= r[1]
                    if group == "radar":
                        i, j = map(int, key.split("/"))
                        assert (
                            abs(p - shape[i, j] / z["ranging_matrix_denominator"])
                            < 4e-15
                        )
    for tag in ("1", "2"):
        for m in z["tables"]["flight"]["4/6"][tag]["metrics"].values():
            assert m["held"][0] > 0
    assert (
        z["changed_blank_swap_arrivals"]["1"] > 0.999
        and z["changed_blank_swap_arrivals"]["2"] < 1e-20
    )


def test_exact_positive_source_bound_all_centers():
    from openwave.interval_echo.row_19.selected_collar.overlapping_shells import kernel

    z = receipt("overlapping_shells")
    for row in z["levels"]:
        profiles = row["center_profile_multiplicities"]
        assert sum(n for _, n in profiles) == row["cells"]
        R = row["inner_hops"]
        for offset in (0, 8):
            upper = max(
                F(kernel(h, 2 * R, offset, False), kernel(h, R, offset, True))
                for h, _ in profiles
            )
            v = row["offsets_in_sixteenths"][str(offset)]
            assert upper == F(v["max_larger_over_smaller"])
            assert v["excluded_by_counts"] == (row["measured_moment_ratio"][0] > upper)
        assert row["offsets_in_sixteenths"]["0"]["excluded_by_counts"]
    assert z["levels"][-1]["offsets_in_sixteenths"]["13"]["excluded_by_counts"]


def test_regional_travel_time_sign_at_all_resolutions():
    z = receipt("regional_travel_times")
    assert F(z["quarter_survival_lower_bound_same_code"]) > F(9, 32)
    assert F(z["limiting_shortening_ratio_half_over_quarter"]) == F(32, 9)
    q = np.array(z["directed_word_score"])
    costs = np.array(z["sharp_return_appointments"][1:])
    for data in z["levels"]:
        second = costs if data["layout"] == "same" else costs.reshape(3, 3).T.ravel()
        assert F(2 * int(q @ costs) * int(q @ second), 1296) == F(
            data["conditional_variance_gain"]
        )
        p = F(data["variance_change"])
        assert data["from_modular_read"][0] <= p <= data["from_modular_read"][1]
        for a in data["probes"].values():
            assert all(sum(t) == z["starts"] and min(t) >= 0 for t in a["counts"])
            assert a["mean_change"][0] <= 0 <= a["mean_change"][1]
            for rate, timer in a["deadlines"].items():
                forecast = data["deadline_prediction"][rate]
                change = F(forecast["population"])
                assert (
                    timer["source_corrected_duration"][0]
                    <= 0
                    <= timer["source_corrected_duration"][1]
                )
                assert (
                    timer["source_corrected_extra_round_trip"][0]
                    <= 2
                    <= timer["source_corrected_extra_round_trip"][1]
                )
                assert sum(timer["extra_round_trip_failed_passed"]) == z["starts"]
                assert (
                    timer["survival_change"][0] <= change <= timer["survival_change"][1]
                )
                assert (
                    forecast["from_modular_read"][0]
                    <= change
                    <= forecast["from_modular_read"][1]
                )
                assert all(sum(t) == z["starts"] for t in timer["failed_passed"])
            if data["level"] in (1, 3):
                t1, t2 = [
                    a["deadlines"][k]["equivalent_duration_change"]
                    for k in ("1/2", "1/4")
                ]
                assert not overlap(t1, t2) and (
                    max(t1 + t2) < 0 if p > 0 else min(t1 + t2) > 0
                )
            lo, hi = a["variance_change"]
            assert lo <= p <= hi and (lo > 0 if p > 0 else hi < 0)


def test_unit_deposit_and_disjoint_source_accelerations():
    from openwave.interval_echo.row_19.selected_collar.initial_acceleration import (
        contact_terms,
    )

    z = receipt("initial_acceleration")
    std = z["phase_energy_standards"]["species"]
    for q in (1, 2):
        reads = std[str(q)]
        for key, v in (("rest", 40 * q * q), ("localized", 40 * q * q + 2)):
            assert reads[key]["energy"][0] <= v <= reads[key]["energy"][1]
        assert (
            reads["first_gap"]["gap_per_laplacian_eigenvalue"][0]
            <= 0.5
            <= reads["first_gap"]["gap_per_laplacian_eigenvalue"][1]
        )
        d = z["unit_and_double"][str(q)]
        assert (
            d["expected"] == -contact_terms()[q] / 2
            and 0 < d["measured"][0] <= d["expected"] <= d["measured"][1]
        )
        for arm in ("empty", "ordinary"):
            for p in d[arm]["phases"]:
                assert (
                    p["held_acceleration"][0]
                    <= p["exact_acceleration"]
                    <= p["held_acceleration"][1]
                )
    rng = np.random.default_rng(28120855)

    def sym(n):
        m = rng.normal(size=(n, n))
        return m + m.T

    h = sym(6)
    a, b = np.zeros_like(h), np.zeros_like(h)
    a[:3, :3] = sym(3)
    b[3:, 3:] = sym(3)
    x = np.diag(rng.normal(size=6))

    def acc(k):
        return 2 * k @ x @ k - k @ k @ x - x @ k @ k

    np.testing.assert_allclose(
        acc(h + a + b) - acc(h + a) - acc(h + b) + acc(h), 0, atol=1e-13
    )
