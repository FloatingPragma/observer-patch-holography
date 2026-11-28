import json

import numpy as np
from scipy.integrate import quad

from . import source_profile as profile


def test_annulus_enclosure_from_independent_direction_average():
    for r, h, a, b in (
        (1.0, 2.0, 0.0, 1.0),
        (4.0, 2.0, 1.0, 2.0),
        (3.0, 2.0, 0.0, 5.0),
        (0.2, 3.0, 2.0, 4.0),
    ):

        def integrand(s):
            f = np.clip((r * r - (s - h) ** 2) / (4 * s * h), 0, 1)
            return 2 * s * f / (b * b - a * a)

        exact = quad(
            integrand, a, b, points=[x for x in (abs(r - h), r + h) if a < x < b]
        )[0]
        np.testing.assert_allclose(
            profile.annulus_fraction(
                np.array(r), np.array(h), np.array(a), np.array(b)
            ),
            exact,
            atol=1e-12,
        )


def test_held_speed_is_unused_by_calibration_and_all_rings_are_scored():
    train = {"a", "b"}
    loss = {
        "a": np.array([0.0, 1.0]),
        "b": np.array([0.0, 3.0]),
        "held": np.array([100.0, 0.0]),
    }
    assert (
        profile.choose(loss, train)
        == profile.choose(loss | {"held": np.array([0.0, 1e40])}, train)
        == 0
    )
    result = json.loads((profile.bao.HERE / "source_profile.json").read_bytes())
    for arm, count in (("bulgeless", 1621), ("resolved_maps", 369)):
        sample = result[arm]
        for name, rows in sample["forecasts"].items():
            residuals = [
                abs(np.array(v["speed_kms"]) - sample["inputs"][k]["observed_kms"])
                for k, v in rows.items()
            ]
            assert sum(map(len, residuals)) == count
            np.testing.assert_allclose(
                np.median(np.concatenate(residuals)),
                sample["scores"][name]["pooled_median_abs_kms"],
            )
