import json

import numpy as np
import pytest

from .ground import certificate, negative_inertia
from .measure import Assembly, HERE, bright, decode, distances


def test_graph_and_clock_commute():
    pair = Assembly(2, 0, 8)
    h = pair.h.toarray()
    rho, _, _ = bright(pair, 84.72475, 131072)
    assert (distances() == distances().T).all()
    np.testing.assert_allclose(h @ rho, rho @ h, atol=5e-14, rtol=0)
    assert 0 < rho.trace() < 1


def test_starts_and_controls():
    data = json.loads((HERE / "counts.json").read_text())
    native, half, zero = decode(data)
    t = native["terminal"][0]
    assert t["unread"]["far_limits"][0] > native["geometry"]["far_limits"][1]
    assert t["unread"]["far_limits"][1] < native["independent"]["far_limits"][0]
    assert (
        t["unread"]["mean_square_limits"][1]
        < native["independent"]["mean_square_limits"][0]
    )
    for row in (native, half, zero):
        t0 = row["terminal"][0]
        assert t0["coarse_from_near"]["far"] == 0 and t0["coarse_from_far"]["far"] == 1
    assert native["passages"][0]["limits"][0] > native["passages"][1]["limits"][1]
    data["records"][0]["passages"][0]["counts"][-1] -= 1
    with pytest.raises(ValueError):
        decode(data)


def test_exact_ground_comparison():
    assert negative_inertia([[2, 1], [1, -2]]) == 1
    assert negative_inertia([[-3, 1, 0], [1, -2, 0], [0, 0, 4]]) == 2
    result = certificate()
    assert 0.974274 < result["joined_far"][0] < result["joined_far"][1] < 0.974276
    assert (
        0.991866
        < result["independent_far"][0]
        < result["independent_far"][1]
        < 0.991867
    )
