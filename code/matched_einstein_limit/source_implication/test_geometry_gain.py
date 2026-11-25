import numpy as np
from .geometry_gain import area_gradient, equal_volume_area, run


def test_fixed_volume_shape_and_independent_directional_secant():
    base = np.array(
        [
            [0.0, 2, 2.1, 2.2],
            [1.95, 0, 2.2, 2.1],
            [2.03, 2.19, 0, 1.9],
            [2.21, 1.99, 1.96, 0],
        ]
    )
    delta = np.array(
        [
            [0.0, 0.3, -0.1, 0.12],
            [-0.2, 0.0, 0.1, 0.0],
            [-0.03, 0.2, 0.0, -0.1],
            [0.2, 0.0, 0.16, 0.0],
        ]
    )
    np.testing.assert_allclose(
        equal_volume_area(base), equal_volume_area(1.7 * base), atol=2e-14
    )
    np.testing.assert_allclose(
        area_gradient(base, delta),
        (
            equal_volume_area(base + 1e-4 * delta)
            - equal_volume_area(base - 1e-4 * delta)
        )
        / 2e-4,
        atol=3e-10,
    )
    regular = np.ones((4, 4)) - np.eye(4)
    np.testing.assert_allclose(
        equal_volume_area(regular), np.sqrt(3) * (6 * np.sqrt(2)) ** (2 / 3)
    )


def test_ordinary_calibration_held_on_every_clock_menu():
    records = run()
    assert len(records["panels"]) == 6
    np.testing.assert_allclose(
        records["pairing_completion"], -0.01437036554, atol=1e-11
    )
    assert all(
        p["sign_of_ordinary_area_gain"] == 1 and p["ordinary_area_derivative"] > 0
        for p in records["panels"]
    )
    ratios = [p["relative_transfer"] for p in records["panels"]]
    assert 0.0009 < min(ratios) < 0.001 and 0.006 < max(ratios) < 0.0061
