"""Independent controls on the native clock decoder and the mass audit."""

import numpy as np
import pytest

from .unchanged_clock import QUERIES, TOLERANCE, _clock_reader, run


def test_clock_reader_recovers_a_line_and_rejects_unresolved_records():
    clicks = 0.4 + 0.2 * np.cos(0.73 * QUERIES) + 0.1 * np.sin(0.73 * QUERIES)
    rate, error = _clock_reader(clicks)
    assert rate == pytest.approx(0.73, abs=TOLERANCE)
    assert error < TOLERANCE
    with pytest.raises(ValueError, match="no resolved"):
        _clock_reader(np.full_like(clicks, 0.4))
    clicks[-1] += 0.02
    with pytest.raises(ValueError, match="held-out"):
        _clock_reader(clicks)


def test_existing_motion_leaves_both_clock_rates_as_source_time():
    audit = run()
    clock = audit["native_cage_clock"]
    assert clock["traveling_minus_parked_charge_at_final_receiver"] > 0.1
    assert clock["traveling_over_parked_clock_rate"] == pytest.approx(1, abs=TOLERANCE)
    assert max(arm["ideal_joint_tensor_residual"] for arm in clock["arms"].values()) < TOLERANCE
    tagged = audit["mass_offset_ramsey"]
    assert tagged["ratio"] == pytest.approx(1, abs=TOLERANCE)
    assert 0 < tagged["matched_spectral_mode_square_root_ratio"] < 0.99
