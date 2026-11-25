from fractions import Fraction
import json

import numpy as np
import pytest

from . import measure, multipoles


def test_against_exact_rational_shell_solver_and_axis_permutation():
    _, q = measure.measured_steps(None)
    phi = multipoles.field(q, 0, 0)
    exact = measure.rails.steady()["per_edge"]
    for tangent, tag in ((5, "1"), (4, "-1")):
        observed = phi[7, tangent, tangent] - phi[8, tangent, tangent]
        assert observed == pytest.approx(float(Fraction(exact[tag])), abs=2e-15)
    x = np.array(measure.predictions(measure.measured_steps(0)[1]))
    y = np.array(measure.predictions(measure.measured_steps(1)[1]))
    np.testing.assert_allclose(y, x[[1, 0, 2, 4, 3, 6, 5]], atol=2e-15, rtol=0)


def test_corner_release_and_reflection():
    with pytest.raises(ValueError, match="negative source"):
        multipoles.source(0.6, 0.5)
    weights = multipoles.source(0, 0.5)
    assert np.sum(weights) == 1
    assert weights[5, 5, 5] == 3 / 16
    # With a stop on the starting corner the first hit has time zero.
    all_corners, _ = measure.acquire(
        np.array(measure.DIRECTIONS + ((0, 0, 0),) * 2), 801, measure.CORNERS, 0.5
    )
    assert min(all_corners[5:]) > 0.2 * measure.SHOTS
    cube = multipoles.field(np.full(3, 0.125), 0, 1) - multipoles.field(
        np.full(3, 0.125), 0, 0
    )
    for size in range(1, 6):
        faces = (slice(6 - size, 6 + size),) * 3
        assert abs(cube[faces].sum()) < 2e-14


def test_held_records_and_rejected_detector_change():
    baseline = json.loads((measure.HERE / "results.json").read_text())
    measure.check(measure.run(), baseline)
    saved = json.loads((measure.HERE / "multipoles.json").read_text())
    fresh = multipoles.run()
    multipoles.verify(fresh, saved)
    saved["panels"][0]["hits"][0] += 1
    with pytest.raises(ValueError, match="hits differs"):
        multipoles.verify(fresh, saved)
