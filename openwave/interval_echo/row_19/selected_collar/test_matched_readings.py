"""Check coverage decisions and the equal-preparation control from observer counts."""

import json

import pytest

from openwave.interval_echo.row_19.selected_collar import (
    matched_readings,
    observer_records,
)


def endpoints(interval):
    return tuple(float(x) for x in interval.strip("[]").split(","))


def test_both_banks_and_sham():
    primary = json.loads(matched_readings.RECEIPT.read_text())
    precision = json.loads(
        matched_readings.RECEIPT.with_name("matched_precision.json").read_text()
    )
    for data in (primary, precision):
        assert observer_records.decode(data["counts"]) == data["decoded"]
    assert not primary["decoded"]["distinct"]
    assert endpoints(precision["decoded"]["predicted_odd_focusing"])[0] < 0
    assert endpoints(precision["decoded"]["predicted_odd_focusing"])[1] > 0
    assert endpoints(precision["decoded"]["measured_minus_predicted"])[1] < -0.64
    counts = precision["counts"]
    counts["stocks"]["collar_plus"] = counts["stocks"]["collar_minus"]
    assert not observer_records.decode(counts)["distinct"]


def test_precision_bank_requires_all_senders():
    data = json.loads(
        matched_readings.RECEIPT.with_name("matched_precision.json").read_text()
    )
    data["counts"]["stocks"]["collar_minus"]["forward"][0].pop()
    with pytest.raises(ValueError, match="Incomplete camera sender set"):
        observer_records.decode(data["counts"])
