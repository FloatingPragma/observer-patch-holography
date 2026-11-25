import sympy as sp
import pytest

from . import jet_lottery as jet


def test_five_matching_source_reads_and_physical_lotteries():
    record = jet.witness(*jet.ward_rows())
    assert record["fixed_readings"] == [
        "z",
        "dz",
        "ddz",
        "outward_number_radius_1",
        "outward_number_radius_2",
    ]
    assert record["rank"] == 5 and record["inverse_defect_upper"] < 0.1
    assert len(record["correction_arms"]) == 5
    assert len(record["fixed_free_weights"]) == 4
    for name, (n, d) in record["fixed_free_weights"].items():
        value = abs(n / d)
        side = "plus" if n >= 0 else "minus"
        bounds = record["lottery_probability_intervals"][side][name]
        assert bounds[0] <= value <= bounds[1]
    for name, band in record["signed_responses"].items():
        assert 0 < band[0] <= band[1], name
    for arm in record["lottery_probability_intervals"].values():
        assert min(v[0] for v in arm.values()) >= 0
        assert sum(v[0] for v in arm.values()) <= 1 <= sum(v[1] for v in arm.values())
    f = record["dyadic_control"]
    for key in ("positive", "negative"):
        assert all(isinstance(x, int) and x >= 0 for x in f[key])
        assert sum(f[key]) + f[key + "_reference"] == f["denominator"]
    assert all(v[0] > 0 for v in f["signed_responses"].values())
    for quantity in record["fixed_readings"]:
        assert f["fixed_response_absolute_upper"][quantity] < 1e-12
    regions = list(record["finite_regional_entropy"].values())
    assert regions[0]["common_admixture"] == regions[1]["common_admixture"]
    for region in regions:
        assert 0 < region["common_admixture"] < 1
        assert region["finite_entropy_gap_nats"][0] > 0


def test_incompatible_input_jet_is_rejected():
    mid, width = jet.ward_rows()
    mid = mid.copy()
    mid[1, 0] += sp.Rational(1, 100)
    with pytest.raises(ValueError, match="Ward enclosures disagree"):
        jet.intersection(mid, width)
