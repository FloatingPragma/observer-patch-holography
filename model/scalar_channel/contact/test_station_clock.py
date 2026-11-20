import math

import numpy as np
import pytest

from . import contact_response, station_clock


def test_early_pointer_and_clock_readout_sign():
    times = np.array([0, math.pi / 2, math.pi, 3 * math.pi / 2])
    p = station_clock.local_pointers({"kappa": 8}, 1, times)
    np.testing.assert_allclose(2*p[:, 0]-1, np.cos(times), atol=1e-13)
    np.testing.assert_allclose(2*p[:, 1]-1, np.sin(times), atol=1e-13)


def test_timeout_clock_is_unresolved_and_negative_counts_rejected():
    e = 6
    spec = {"endpoints": list(range(e)), "shots": 100, "clock_alpha": 1e-6,
            "clock_rates": [1, 1.25], "clock_slots": [0, 16, 17], "dt": .5}
    counts = np.zeros((2, 2, e, e, 3, 2, 3), dtype=np.int64)
    counts[:, :, ~np.eye(e, dtype=bool), ..., -1] = 100
    emission = np.zeros((e, 3, 2, 2), dtype=np.int64)
    emission[..., 0] = 100
    result = station_clock.score(spec, counts, emission)
    assert result[0][0]["resolved_first_slot_fringes"] == 0
    counts[0, 0, 0, 1, 0, 0] = [-1, 1, 100]
    with pytest.raises(ValueError):
        station_clock.score(spec, counts, emission)


def test_degenerate_levels_and_disconnected_mark():
    chart = np.array([[2, 2, 0, 0], [2, 2, 0, 0], [0, 0, 2, 2], [0, 0, 2, 2]])
    rank, moments = contact_response.marked_space(chart, np.ones(4), 0, range(4))
    assert rank == 2
    assert moments == {0: (0, 1), 1: (1, 1)}
