"""Independent sky projection and planar matrix checks for the merger constraint."""

from fractions import Fraction as F
import json
import math

from . import measure


def test_published_positions_and_projection():
    data = json.loads((measure.HERE / "results.json").read_text())
    scale = 15 * math.cos(math.radians(55 + 56 / 60 + 30 / 3600))
    west_peak = ((16.0 - 21.2) * scale + 7.1, 30.0 - 35.1 + 6.5)
    east_gas = (9 * scale, -5.9)
    east_peak = ((35.3 - 21.2) * scale - 21.2, 30.0 - 56.3 + 17.7)
    saved_west = data["west_peak_east_north_arcsec"]
    assert math.dist(west_peak, saved_west) < 1e-8
    margin = data["minimum_combined_position_allowance_arcsec"]
    assert abs(math.dist(east_peak, east_gas) - margin) < 1e-8
    assert 36.47 < data["western_peak_beyond_west_arcsec"] < 36.48
    assert margin < data["eastern_peak_beyond_east_arcsec"]["northern"]
    assert 96.24 < margin * 4.413 < 96.26
    assert json.loads(json.dumps(measure.build())) == data


def test_independent_shapes_can_create_one_exterior_maximum():
    # det A = det B = 1. A and B have different principal axes.
    # At t=1 the sum of the two negative quadratic forms peaks beyond west.
    # Each clipped profile can have positive extent at that point.
    for t in (F(1, 8), F(1), F(4)):
        # A=[[1,3],[3,10]], B=[[1,4],[4,17]], centres (0,0),(1,0).
        a, b, d = 1 + t, 3 + 4*t, 10 + 17*t
        determinant = a*d - b*b
        direct = ((d*t - b*4*t)/determinant, (a*4*t - b*t)/determinant)
        denominator = 1 + 3*t + t*t
        rational = ((-2*t + t*t)/denominator, t/denominator)
        assert direct == rational
        assert direct[0] < 1
        if t == 1:
            assert direct == (F(-1, 5), F(1, 5))
        assert 0 < t/(1 + t) < 1  # round-shape control


def test_single_peak_limiting_shapes_cancel():
    r = measure.build()
    p = r["west_peak_east_north_arcsec"]
    rays = [measure.difference(c, p) for c in r["plasma_centres_east_north_arcsec"]]
    ordered = sorted(math.atan2(y, x) for x, y in rays)
    q = r["west_only_upper_bound_on_smaller_axis_ratio"]

    def normal(angle, sign):
        # In the major/minor basis the input ray has angle +/- atan(q)
        # and the quadratic has eigenvalues 1 and 1/q^2.
        tilt = sign * math.atan(q)
        axis = angle - tilt
        a, b = math.cos(tilt), math.sin(tilt)/(q*q)
        return measure.unit((a*math.cos(axis) - b*math.sin(axis),
                             a*math.sin(axis) + b*math.cos(axis)))

    left, right = normal(ordered[0], -1), normal(ordered[1], 1)
    assert math.dist(left, tuple(-x for x in right)) < 1e-12
