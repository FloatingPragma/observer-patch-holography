import json
from fractions import Fraction

from .measure import run


def test_native_nonlinear_remainder_and_stationary_pole():
    data = run()
    native = data["written"]
    endpoints = [json.loads(value) for value in native["entropy_slopes"]]
    assert (
        endpoints[0][0] > 0.00330
        and endpoints[1] == [0.0, 0.0]
        and endpoints[2][1] < -0.00294
    )
    assert (
        Fraction(native["fisher_lower"])
        < Fraction(1, 324)
        < Fraction(native["fisher_upper"])
    )
    transfer = native["reversed_to_shared"]
    actual = json.loads(transfer["entropy_change"])[0]
    predicted = json.loads(transfer["linear_entropy_prediction"])[0]
    error = json.loads(transfer["relative_entropy_remainder"])[0]
    assert 0 < actual < 0.02 * predicted and abs(actual - (predicted - error)) < 4e-18
    for row in native["resolved_annuli"]:
        for arm in row["information_cost_reverse_fresh_shared"]:
            assert (
                abs(arm["mean_nats"] - arm["measured_nats"])
                <= 5 * arm["standard_error"]
            )
    cell = data["angular"]
    assert cell["mode_meter"]["rotations"] == 30
    assert cell["mode_meter"]["worst_rotation_error"] < 5e-13
    for arm in cell["sampled_classification"]:
        assert sum(arm["axial_transverse_other"]) == arm["trials"]
    regions = cell["nested_regions"]
    fractions = [r["full_information_fraction_home"] for r in regions]
    assert fractions == sorted(fractions) and abs(fractions[-1] - 1) < 1e-11
    assert (
        regions[0]["home"]["linear_entropy_slope"]
        > 0
        > regions[1]["home"]["linear_entropy_slope"]
    )
    assert abs(regions[-1]["home"]["linear_entropy_slope"]) < 1e-17
