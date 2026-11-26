from . import measure


def test_observer_records_under_native_one_interval_relaxation():
    r = measure.run()
    assert r["rank"] == 36
    assert r["right_inverse_defect_upper"] < 0.001
    assert r["step"] < r["reachable_halfwidth_lower"]
    assert r["native_exchange_probability_lower"] > 0
    assert r["stay_probability_lower"] > 0.5
    assert r["fixed_record_changes_on_deposited_stock"]["pair_pixel_upper"] < 1e-14
    assert r["interior_entropy_change_exact_stock_lower_nats"] > 0.0003
    expected = r["interior_entropy_change_deposited_nats"]
    lo, hi = r["simultaneous_99_counted_entropy_change_nats"]
    assert 0 < lo < expected[0] <= expected[1] < hi
    lag = r["after_one_further_interval_first_ray_change"]
    lo, hi = r["simultaneous_99_late_herald_interval"]
    assert 0 < lo < lag[0] <= lag[1] < hi
    assert all(sum(b["bulk"]) == r["shots_per_setting"] for b in r["counts"])
