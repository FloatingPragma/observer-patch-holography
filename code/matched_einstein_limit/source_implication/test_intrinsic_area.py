from copy import deepcopy
from itertools import combinations
from math import log
from pathlib import Path
import mpmath as mp
import numpy as np
from .receipt_replay import pinned_receipt
from .intrinsic_area import full_entropy_slopes, head_source, radar_banks, run, survey_corners
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


@mp.workdps(50)
def independent_endpoint(corners, product, delta, odd, pairing, adjuster):
    """Reconstruct one recorded preparation with Cayley–Menger volumes and 50-digit head entropy."""
    u, s, c = map(mp.mpf, (odd, pairing, adjuster))
    strengths = (s + c / (1 - u), s - c / u)
    parities = (1 - u, u)
    ranges = mp.matrix(4)
    entropy = mp.mpf(0)
    for k, (mass, strength) in enumerate(zip(parities, strengths, strict=True)):
        ranges += mass * (
            (1 - strength) * mp.matrix(corners[k, 0]) + strength * mp.matrix(corners[k, 1])
        )
        eigenvalues = mp.eighe(mp.matrix(product) + strength * mp.matrix(delta), eigvals_only=True)
        entropy -= mp.fsum(mass * x * mp.log(mass * x) for x in eigenvalues)

    gains = [
        mp.root(mp.fprod(ranges[i, j] / ranges[j, i] for j in range(4) if j != i), 4)
        for i in range(4)
    ]
    squared = mp.matrix(4)
    for i, j in combinations(range(4), 2):
        squared[i, j] = squared[j, i] = ranges[i, j] * ranges[j, i] / (gains[i] * gains[j])

    def simplex(sites):
        dimension = len(sites) - 1
        cm = mp.ones(dimension + 2)
        cm[0, 0] = 0
        for i, a in enumerate(sites, start=1):
            for j, b in enumerate(sites, start=1):
                cm[i, j] = squared[a, b]
        factor = (-1) ** (dimension + 1) / (2 ** dimension * mp.factorial(dimension) ** 2)
        return mp.sqrt(factor * mp.det(cm))

    area = mp.fsum(simplex(face) for face in combinations(range(4), 3))
    return np.array([float(area), float(simplex(range(4))), float(entropy)])


def test_volume_resolution_can_hide_more_than_16_entropy_ulps():
    product, delta = head_source()
    corners = survey_corners(radar_banks(), 1, 6)
    receipt = pinned_receipt(
        Path(__file__).with_name("intrinsic_area.json"),
        "a0cca8487ec7227a2994cd0ed7c393fcd7f78d82692ea9852d22678708ace075",
    )
    arm = receipt["finite_pairs"]["direct_calibration"]["pairing"]["arms"][0]
    args = corners, product, delta, arm["odd"], arm["mean_pairing"]
    base = independent_endpoint(*args, arm["adjuster"])
    nearby = independent_endpoint(*args, arm["adjuster"] + 2e-11)
    assert abs(nearby[1] - base[1]) <= 8 * np.spacing(base[1])
    assert abs(nearby[2] - base[2]) > 16 * np.spacing(base[2])


def test_equal_volume_and_two_observable_matches_on_admissible_sources():
    r = run()
    frozen = pinned_receipt(
        Path(__file__).with_name("intrinsic_area.json"),
        "a0cca8487ec7227a2994cd0ed7c393fcd7f78d82692ea9852d22678708ace075",
    )
    live, saved = deepcopy(r), deepcopy(frozen)
    product, delta = head_source()
    corners = survey_corners(radar_banks(), 1, 6)
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
                # Calibrated knobs can differ between platforms. Compare each
                # endpoint to the same preparation at higher precision.
                expected = independent_endpoint(corners, product, delta, u, s, c)
                np.testing.assert_array_max_ulp(arm["records"], expected, maxulp=16)
            # The volume changes only about 4.7e-5 per knob unit. A few
            # floating-point ulps in that observable move the fitted knob.
            assert abs(
                adjusters[1] - adjusters[0]
            ) * volume_per_adjuster <= 8 * np.spacing(central["volume"])
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
