from copy import deepcopy
from math import log
from pathlib import Path
import numpy as np
from .receipt_replay import pinned_receipt
from .intrinsic_area import full_entropy_slopes, head_source, run
from model.experiments.interval_echo_einstein_station.collar_pairing.measure import (
    compare,
)


def test_complete_density_has_separate_ordinary_and_correlation_tangents():
    p, delta = head_source()
    s = 0.5

    def entropy(strength):
        w = np.linalg.eigvalsh(p + strength * delta)
        return -w @ np.log(w)

    diff = 1e-5
    ordinary, collar = full_entropy_slopes(p, delta, 0.25, s)
    np.testing.assert_allclose(ordinary, log(3), atol=1e-13)
    np.testing.assert_allclose(
        collar, (entropy(s + diff) - entropy(s - diff)) / (2 * diff), atol=2e-10
    )


def test_equal_volume_and_two_observable_matches_on_admissible_sources():
    r = run()
    frozen = pinned_receipt(
        Path(__file__).with_name("intrinsic_area.json"),
        "a0cca8487ec7227a2994cd0ed7c393fcd7f78d82692ea9852d22678708ace075",
    )
    live, saved = deepcopy(r), deepcopy(frozen)
    central = next(
        p
        for p in live["panels"]
        if (p["tick"], p["horizon"], p["odd"], p["pairing"]) == (0.5, 6, 0.25, 0.5)
    )
    volume_per_adjuster = abs(
        central["volume_slopes"][0] / central["direct_volume_adjusters"][0]
    )
    for label in ("ordinary", "pairing"):
        panels = [
            output["finite_pairs"]["direct_calibration"][label]
            for output in (live, saved)
        ]
        # Directly certify each frozen and live finite difference against its two
        # measured endpoint records, before comparing the records themselves.
        for panel in panels:
            lo, hi = panel["arms"]
            step = panel["step"]
            np.testing.assert_array_equal(
                panel.pop("centered_slopes"),
                np.subtract(hi["records"], lo["records"]) / (2 * step),
            )
        for new, old in zip(panels[0]["arms"], panels[1]["arms"], strict=True):
            adjusters = []
            for arm in (new, old):
                c = arm.pop("adjuster")
                adjusters.append(c)
                u, s = arm["odd"], arm["mean_pairing"]
                np.testing.assert_allclose(
                    arm["dephase_probabilities"],
                    2 * np.array([u * (1 - s) + c, u * s - c]) / np.array([1 - s, s]),
                    atol=1e-15,
                    rtol=0,
                )
                assert abs(arm["records"][1] - central["volume"]) <= 8 * np.spacing(
                    central["volume"]
                )
            # The volume changes only about 4.7e-5 per knob unit. A few
            # floating-point ulps in that observable move the fitted knob.
            assert abs(
                adjusters[1] - adjusters[0]
            ) * volume_per_adjuster <= 8 * np.spacing(central["volume"])
            np.testing.assert_array_max_ulp(new["records"], old["records"], maxulp=16)
    compare(live, saved)
    assert len(r["panels"]) == 54
    assert all(
        p["entropy"][0] > 0 and p["shape"][0] > 0 and p["volume_slopes"][0] < 0
        for p in r["panels"]
    )
    a, b = r["finite_pairs"]["equal_volume"]
    da, dv, ds = np.subtract(b["records"], a["records"])
    assert da > 1e-5 and abs(dv) < 1e-13 and ds > 4e-5
    held = r["finite_pairs"]["equal_area_and_volume"]
    changes = np.ptp([arm["records"] for arm in held], axis=0)
    assert max(changes[:2]) < 1e-12 and changes[2] > 0.0003
    for arm in held:
        assert all(0 < u < 1 for u in arm["dephase_probabilities"])
    assert all(
        p["direct_volume_ordinary_eta"] < 0 and p["direct_volume_transfer"] < 0
        for p in r["panels"]
    )
    central = next(
        p
        for p in r["panels"]
        if (p["tick"], p["horizon"], p["odd"], p["pairing"]) == (0.5, 6, 0.25, 0.5)
    )
    for i, name in enumerate(("ordinary", "pairing")):
        slope = r["finite_pairs"]["direct_calibration"][name]["centered_slopes"]
        np.testing.assert_allclose(
            slope[0], central["direct_volume_area_derivatives"][i], rtol=1e-4, atol=1e-8
        )
        assert abs(slope[1]) < 2e-8
        np.testing.assert_allclose(slope[2], central["entropy"][i], rtol=0.001)
