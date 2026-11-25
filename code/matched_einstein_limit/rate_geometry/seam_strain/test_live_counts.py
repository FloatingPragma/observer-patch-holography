"""Clock-count data integrity, native calibration and population cross-checks."""

from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest

from . import live_counts as counts
from . import live_radar, live_survey

ARCHIVE = Path(__file__).with_name("live_counts.json")
ARCHIVE_HASH = "77ed3f214c5f4904031284915335e7b3a4164ddac73c92db646072bfddf4eed8"


def archive():
    payload = ARCHIVE.read_bytes()
    assert sha256(payload).hexdigest() == ARCHIVE_HASH
    return json.loads(payload)


def test_preparations_realize_signed_tangents_with_positive_density_operators():
    prep = counts.preparations()
    ref, changes, _ = live_survey.retained(True)
    np.testing.assert_allclose(prep["reference"], ref, atol=2e-14, rtol=0)
    for density in prep.values():
        np.testing.assert_allclose(density, density.conj().T, atol=2e-14, rtol=0)
        assert abs(density.trace() - 1) < 2e-13
        assert np.linalg.eigvalsh(density)[0] > -2e-14
    for arm, tangent in changes.items():
        if arm in ("ordinary", "decorrelated"):
            observed = prep[arm] - prep["reference"]
        else:
            scale = np.sqrt(2) if arm == "hidden" else 2
            observed = (prep[f"{arm}_plus"] - prep[f"{arm}_minus"]) / scale
        np.testing.assert_allclose(observed, tangent, rtol=0, atol=2e-14)


def test_every_frozen_bin_and_timeout_score_matches_live_native_detector():
    record = archive()
    actual = counts.native_probabilities()
    ranges = (
        sum(live_radar.range_effects(1, q, True, counts.STATIONS) for q in (-1, 1)) / 2
    )
    assert record["outcomes"] == list(counts.OUTCOMES)
    assert record["score_times_four"] == list(counts.TICKS)
    for name, rho in counts.preparations().items():
        nums = np.array(record["probability_numerators"][name])
        assert np.all(nums >= 0) and np.all(nums.sum(axis=1) == counts.DEN)
        p = nums / counts.DEN
        np.testing.assert_allclose(p, actual[name], atol=1e-11, rtol=0)
        expected_range = live_survey.read(ranges, rho)
        np.testing.assert_allclose(
            p @ np.array(counts.TICKS) / 4,
            expected_range,
            atol=float(counts.MEAN_ALLOWANCE),
            rtol=0,
        )
        for bank in counts.BANKS:
            for e, row in enumerate(record["counts"][bank][name]):
                # Independent check of every drawn mean against its nominal law.
                mean = counts.score(row)
                radius = record["inference"]["hoeffding_half_width"][1]
                assert abs(float(mean) - p[e] @ np.array(counts.TICKS) / 4) < radius
                assert row[-1] > 0 and row[-2] > 0
    # Exercise the actual acquisition as independent data on each platform.
    # Exact archived counts and their native calibration are checked above.
    p = record["probability_numerators"]["hidden_plus"][3]
    fresh = counts.draw_row(
        p, 1, counts.PREPS.index("hidden_plus"), 3, seed=counts.SEED + 1
    )
    amplitude = counts.score(fresh)
    expected = np.dot(np.asarray(p) / counts.DEN, counts.TICKS) / 4
    assert (
        abs(float(amplitude) - expected)
        < record["inference"]["hoeffding_half_width"][1]
    )
    assert fresh != record["counts"]["boundary"]["hidden_plus"][3]
    bad = record["counts"]["volume"]["ordinary"][0].copy()
    bad[-1] -= 1
    with pytest.raises(ValueError, match="every|all registered"):
        counts.score(bad)


def test_independent_count_banks_recover_signed_population_with_joint_coverage():
    record = archive()
    results = counts.infer(record["counts"])
    assert results == record["inference"]
    assert record["shots_per_setting"] == counts.N
    assert record["total_shots"] == counts.SETTINGS * counts.N
    assert record["edges"] == [
        [counts.STATIONS[a], counts.STATIONS[b]] for a, b in counts.EDGES
    ]
    for name in counts.PREPS:
        assert all(
            x != y
            for x, y in zip(
                record["counts"]["volume"][name],
                record["counts"]["boundary"][name],
                strict=True,
            )
        )
    survey = json.loads(ARCHIVE.with_name("live_survey.json").read_text())
    (panel,) = (
        r
        for r in survey["panels"]
        if r["coupled"] and r["read"] == "pooled" and r["deadline"] == 12
    )
    (center,) = (c for c in panel["centers"] if c["apex"] == 5)
    for arm in counts.ARMS:
        got = results["matched"][arm]
        nominal = center["matched_boundary"]["right"][arm]
        for key, value in (
            ("right_rotation", nominal["rotation"]),
            (
                "left_control_fixed_volume_area",
                center["matched_boundary"]["left"][arm]["area_slope"],
            ),
            ("fixed_volume_boundary_area", nominal["area_slope"]),
            (
                "fixed_volume_cone_area_per_radius_squared",
                center["slopes"][arm]["angle"] / 3,
            ),
            (
                "rotated_cone_area_per_radius_squared",
                (
                    center["slopes"][arm]["angle"]
                    + nominal["rotation"] * center["slopes"]["right"]["angle"]
                )
                / 3,
            ),
        ):
            lo, hi = got[key]
            assert lo < value < hi
        assert got["fixed_volume_boundary_area"][1] < 0
        if arm != "ordinary":
            charge = panel["source_slopes"]["incoming"]
            ratio = charge[arm] / charge["ordinary"]
            assert results["right_matched_area_over_ordinary"][arm][1] < ratio
        cone = got["fixed_volume_cone_area_per_radius_squared"]
        assert cone[0] > 0 if arm == "ordinary" else cone[1] < 0
