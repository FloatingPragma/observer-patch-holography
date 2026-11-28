import json

import numpy as np
import pytest

from openwave.interval_echo.row_19.selected_collar import reference_read as meter


def test_decoding_and_source_reversal():
    data = json.loads((meter.measure.HERE / "reference_read.json").read_bytes())
    for n in meter.SHOTS:
        banks = data["cases"][str(n)]["tapes"]
        with pytest.raises(ValueError, match="stated number"):
            meter.decode(banks, n * 2)
        result = meter.decode(banks, n)
        bounds = result["equilibrium"]["interval"]
        assert bounds[0] < 0.00022946 < bounds[1]
        assert (bounds[0] > 0) == (n == meter.SHOTS[-1])
        control = result["decorrelated"]["interval"]
        assert control[0] < 0 < control[1]
        reverse = meter.decode(
            banks | {"plus": banks["minus"], "minus": banks["plus"]}, n
        )
        np.testing.assert_allclose(
            reverse["equilibrium"]["interval"], -np.array(bounds[::-1])
        )


def test_source_implies_area_sign_with_an_independent_check():
    from openwave.interval_echo.row_19.selected_collar import source_area

    result = source_area.run()
    low, high = result["area_derivative_from_source"]
    for key in ("area_derivative_from_radar", "finite_endpoint_area_change"):
        a, b = result[key]
        assert low < a <= b < high < 0
    assert max(abs(x) for x in result["volume_derivative"]) < 3e-15
    for chart in result["flat_realizations_of_symmetric_ranges"]:
        assert chart["smallest_gram_eigenvalue"] > 0
        assert chart["maximum_symmetric_range_residual"] < 1e-12


def test_release_clocks_and_six_station_count_replay():
    from mpmath import iv
    from openwave.interval_echo.row_19.selected_collar import (
        release_clock,
        six_station_read,
    )

    clock = release_clock.run()
    assert len(clock["native_exterior_clocks"]["clock_cells"]) == 24
    assert clock["native_exterior_clocks"]["signed_endpoint_gap_change"] == 0
    for outward, zero, inward in clock["endpoints"]:
        assert outward["predicted_release_acceleration_from_added_clocks"] > 0
        assert zero["predicted_release_acceleration_from_added_clocks"] == 0
        assert inward["predicted_release_acceleration_from_added_clocks"] < 0

    stored = json.loads((meter.measure.HERE / "six_station_read.json").read_bytes())
    original_precision = iv.dps
    try:
        iv.dps = 35
        replay = six_station_read.decode(stored["counts"])
        measured = six_station_read.read_regional_source(
            stored["counts"], stored["regional_echo_reference"]
        )
        for row in stored["local_calibrations"]:
            assert (
                six_station_read.bounds(measured[tuple(row["stations"])][1])
                == row["regional_echo_charge_from_counts"]
            )
    finally:
        iv.dps = original_precision
    assert replay["changes"] == stored["changes"]
    assert all(
        lo > 0 for minors in replay["leading_gram_minors"].values() for lo, _ in minors
    )
    exact = json.loads((meter.measure.HERE / "source_area.json").read_bytes())
    a, b = replay["changes"][0]["area"]
    assert (
        a
        < exact["finite_endpoint_area_change"][0]
        <= exact["finite_endpoint_area_change"][1]
        < b
        < 0
    )
    assert any(row["area"][0] > 0 for row in replay["changes"])
