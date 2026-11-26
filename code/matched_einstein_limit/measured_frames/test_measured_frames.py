import json
from fractions import Fraction as F

import pytest

from . import measure


def test_known_high_boosts_and_slower_calibrations():
    from ..frame import FRAME

    original = [[(x / u[0], x / u[0]) for x in u[1:]] for u in FRAME]
    assert measure.certify(original) == (0, F(17, 8))
    stored = json.loads((measure.HERE / "examples.json").read_text())
    assert measure.examples() == stored
    assert all(F(r["measurable_gain_upper"]) >= F(r["exact_gain"]) for r in stored)


def test_unresolved_and_unphysical_enrollments():
    rows = measure.symmetric_menu(F(1, 4))
    rows[1] = rows[0]
    with pytest.raises(ValueError, match="distinguish"):
        measure.certify(rows)
    rows[1] = [(F(0), F(1))] * 3
    with pytest.raises(ValueError, match="light cone"):
        measure.certify(rows)


def test_noise_can_prevent_certification():
    rows = measure.symmetric_menu(F(1, 4), F(1, 7))
    with pytest.raises(ValueError, match="too wide"):
        measure.certify(rows)


def test_arbitrary_unequal_point_velocities_bound_an_independent_tensor():
    from ..frame import PAIRS, apply, inverse

    rows = measure.symmetric_menu(F(1, 4))
    # Each instrument reports its own unequal velocities, in the same tetrad.
    points = [
        [lo + F((r + 1) * (i + 2), 2048) for i, (lo, _) in enumerate(v)]
        for r, v in enumerate(rows)
    ]
    boxes = [[(x - F(1, 8192), x + F(1, 8192)) for x in v] for v in points]
    beta, gain = measure.certify(boxes)
    assert 0 < beta < 1
    directions = []
    for v in points:
        u = [1, *v]
        dilation = 1 - sum(x * x for x in v)
        directions.append(
            [u[i] * u[j] * (1 if i == j else 2) / dilation for i, j in PAIRS]
        )
    direct = inverse(directions)
    assert max(sum(abs(c) for c in row) for row in direct) <= gain
    tensor = [F((-1) ** i, i + 1) for i in range(10)]
    assert list(apply(direct, apply(directions, tensor))) == tensor


def test_rational_velocity_strings_use_numerical_order():
    rows = measure.symmetric_menu(F(1, 4))
    rows[1][0] = ("1/4", "1/3")
    expected = [[(F(lo), F(hi)) for lo, hi in row] for row in rows]
    assert measure.certify(rows) == measure.certify(expected)
    rows[1][0] = ("1/3", "1/4")
    with pytest.raises(ValueError, match="reversed"):
        measure.certify(rows)
