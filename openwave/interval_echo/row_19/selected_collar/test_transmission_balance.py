"""Checks on signs, nonlinear geometry and the recorded regional evidence."""

import json
from itertools import permutations

import numpy as np
from scipy.linalg import logm

from openwave.interval_echo.row_19.selected_collar import erasure_balls as balls
from openwave.interval_echo.row_19.selected_collar import transmission_balance as meter


def entropy(X, P):
    return np.trace(X @ (logm(X) - logm(P))).real


def test_erasure_slope_including_timeout():
    X = np.array([[0.6, 0.15], [0.15, 0.4]])
    P = np.array([[0.45, -0.03], [-0.03, 0.55]])
    success = np.diag([0.8, 0.25])

    def channel(Z):
        out = np.zeros((3, 3))
        out[:2, :2] = success @ Z @ success
        out[2, 2] = np.trace(Z) - np.trace(out)
        return out

    def loss(Z):
        return entropy(Z, P) - entropy(channel(Z), channel(P))

    h = 1e-5
    numerical = (loss(X - h * (P - X)) - loss(X + h * (P - X))) / (2 * h)
    predicted = loss(X) + entropy(P, X) - entropy(channel(P), channel(X))
    assert predicted > 0
    np.testing.assert_allclose(numerical, predicted, rtol=1e-8)


def test_all_timing_partitions():
    ids = meter.reg.partitions()
    assert len(set(map(tuple, ids))) == 115975
    for did in (0, 1, 8, 511, 10000, len(ids) - 1):
        bits = ids[did]
        assert np.array_equal(np.sum(meter.reg.MASK[bits], axis=0), np.ones(10))
    assert any(np.count_nonzero(x) == 1 for x in ids)
    assert any(np.count_nonzero(x) == 10 for x in ids)


def test_fixed_volume_solid_angle_slope():
    from matched_einstein_limit.source_implication.full_algebra.enclosure import Jet

    iv = balls.iv
    iv.dps = 35
    # Tetrahedron with three orthogonal vectors from the origin. Perturb
    # one direction towards another with cosine t, leaving all lengths one.
    square = np.zeros((4, 4), dtype=object)
    for i, j in permutations(range(4), 2):
        v = 1 if 0 in (i, j) else 2
        d = -2 if {i, j} == {1, 2} else 0
        square[i, j] = Jet(iv.mpf(v), np.array([iv.mpf(d)]), np.array([[iv.mpf(0)]]))
    omega = balls.solid_angle_jets(square, 0, (1, 2, 3))
    assert omega.v.a <= iv.pi / 2 <= omega.v.b
    assert omega.g[0].a <= -1 <= omega.g[0].b


def test_port_38_has_an_independent_positive_residual():
    r = json.loads((meter.HERE / "erasure_balls.json").read_bytes())
    witnesses = [q for q in r["balls"] if q["center_port"] == 38]
    assert len(witnesses) == 6 and len(r["balls"]) == 36
    for q in witnesses:
        assert q["area_slopes"][0][1] < 0 < q["area_slopes"][3][0]
        assert q["ordinary_area_gain"][0] > 0
        assert q["additional_score_needed_any_receiver"][1] < -0.18
        assert (
            q["moving_aperture"]["any_receiver_ordinary_calibrated_erasure_residual"][
                "hidden"
            ][0]
            > 0.17
        )
