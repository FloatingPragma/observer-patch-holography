import json
from pathlib import Path

from ..receipt_replay import pinned_receipt, verify_recorded_band
from . import pulse_bath


def test_native_sector_relaxation_and_fixed_volume_slope():
    observed = pulse_bath.run()
    saved = pinned_receipt(
        Path(pulse_bath.__file__).with_suffix(".json"),
        "e268d53dd7768837c44f871f49397575e0d48d0d89ddc535e196532845418369",
    )
    # The gain uses a certified matrix-log eigensolver proposal. Bound its
    # signed objective derivative for both proposals in the same printed band.
    intervals = []
    for result in (saved, observed):
        interval = json.loads(
            result["matched_measurement"].pop("fixed_volume_objective_slope")
        )
        low, high = interval
        assert -0.004815774 < low < high < -0.004815759
        intervals.append(interval)
    assert max(x[0] for x in intervals) < min(x[1] for x in intervals)
    assert saved == observed
    collision = observed["native_collision"]
    assert collision["squared_norm_loss_floor"] == "1/61440"
    assert all(
        sum(sector["commutator_gram"].values()) == 15 for sector in collision["sectors"]
    )


def test_reused_radar_has_a_distinct_matched_volume_equilibrium():
    from . import reuse_equilibrium

    observed = reuse_equilibrium.run()
    saved = pinned_receipt(
        Path(reuse_equilibrium.__file__).with_suffix(".json"),
        "5a6fbdf22798896e69512add151c5fd9933e0746ff9ad8e88877c486935c7cf0",
    )
    assert json.loads(observed["fixed_volume_correlation_only_slope"])[1] < 0
    assert json.loads(observed["fixed_area_and_volume_entropy_slope"])[1] < -0.00033
    # Only terms involving a certified local matrix log depend on the proposal.
    for key in (
        "fixed_volume_marginal_coin_slope",
        "fixed_volume_correlation_only_slope",
        "fixed_area_and_volume_entropy_slope",
    ):
        verify_recorded_band(saved.pop(key), observed.pop(key))
    for key in ("entropy_slopes", "would_require_gains"):
        for old, new in zip(saved.pop(key), observed.pop(key), strict=True):
            verify_recorded_band(old, new)
    assert observed == saved
