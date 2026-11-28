"""Check conserved observables and the regression separating radial and angular rules."""

from fractions import Fraction
from math import acos, sqrt

import pytest

from .forecast import calibrated_mu, forecast
from .measure import RADIUS, TICK, observe, paired_observation


def test_flat_endpoints_from_euclidean_geometry():
    # An independent off-protocol direction exercises the local boost and reflection.
    s = 0.625
    shot = observe(0., s)
    assert shot["event"] == "returned"
    assert abs(shot["mirror_address_ticks"]*TICK - 2*acos(s)) < 1e-11
    assert abs(shot["echo_clock_ticks"]*TICK - 4*RADIUS*sqrt(1-s*s)) < 1e-11
    assert abs(shot["echo_address_ticks"]*TICK) < 1e-11


def test_bending_and_delay_require_the_angular_rule():
    gain = [0.25, 0.25]
    mu = calibrated_mu(gain, Fraction(19, 33))
    predicted = forecast(mu, 0.5)
    direct = paired_observation(19/132, 0.5)
    altered = paired_observation(19/132, 0.5, change=True)
    for key in ("delay", "bending"):
        lo, hi = predicted[key]
        assert lo <= direct[key] <= hi
        assert 0 < altered[key] < lo
    assert abs(direct["matched_blank_mirror_address_ticks"]-direct["mirror_address_ticks"])*TICK < 1e-10


def test_unresolved_exterior_forecast_fails():
    with pytest.raises(ValueError, match="turning point"):
        forecast([0., 0.01], 0.125)
