"""Exercise unresolved and invalid nonzero clock banks in one acquisition."""

from tools.strict_json import load

from . import donor_release_checked as acquisition
from .transfer_clock import SHOTS


def test_nonzero_unresolved_clock_is_distinct_from_an_invalid_bank():
    receipt = load(acquisition.HERE / "clock.json")
    receipt["runs"][0]["counts"] = [0] * len(acquisition.clock.T)
    receipt["runs"][1]["counts"][0] = SHOTS + 1

    first, second, held, *_, zero = acquisition.run(receipt)["launches"]
    assert first["strength_in_simulator"] == 0.5
    assert first["disposition"] == "inconclusive_clock"
    assert second["strength_in_simulator"] == -0.5
    assert second["disposition"] == "invalid_clock_bank"
    for skipped, original in zip((first, second), receipt["runs"][:2], strict=True):
        assert skipped["clock_counts"] == original["counts"]
        assert skipped["clock_calibration"] is None
        assert "launches" not in skipped
    assert held["disposition"] == "passed"
    assert [row["column"] for row in held["launches"]] == [14, 15]
    assert zero["disposition"] == "inconclusive_clock"
