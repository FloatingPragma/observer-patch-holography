"""Protect the adjacent-pair event from accidental replacement by shared ownership."""

import json

import numpy as np

from .measure import HERE, SHOTS, TIMES


def test_strict_position_event_and_denominator():
    packet = json.loads((HERE / "results.json").read_text())
    launches = packet["flights"]
    assert [(x["owner"], bool(x["appointments"])) for x in launches] == [
        (0, True),
        (6, False),
        (5, True),
        (15, False),
        (10, True),
        (15, False),
    ]
    for launch in launches:
        appointments = launch["appointments"]
        if not appointments:
            assert launch["admission"] < 0.001
            continue
        assert tuple(a["t"] for a in appointments) == TIMES
        for a in appointments:
            bound = a["bound"]
            strict = np.array(bound["upper_address_distance_then_other"])
            owners = np.array(a["owner_distance_then_separated"])
            assert len(strict) == len(owners)
            assert np.all(strict[:-1] <= owners[:-1] + 1e-14)
            assert strict[-1] > owners[-1]
            np.testing.assert_allclose(strict.sum(), 1.0, atol=2e-13, rtol=0)
            assert sum(bound["clicks"]) == SHOTS
            assert sum(a["clicks"]) == SHOTS
            assert abs(bound["two_or_more"] - strict[2:-1].sum()) < 1e-13
        assert appointments[0]["bound"]["two_or_more"] < 1.6e-12
        assert all(a["bound"]["resolved"] for a in appointments[1:])
        assert (
            appointments[-1]["bound"]["two_or_more"]
            < appointments[-1]["shared_owner_two_or_more"]
        )


def test_full_generator_and_spectral_counts():
    packet = json.loads((HERE / "results.json").read_text())
    levels = np.array(packet["lines"])
    thresholds = packet["threshold_counts"]
    assert [s["upper_lines_below"] for s in thresholds] == [0, 24, 24, 24, 32]
    assert (
        thresholds[0]["energy"] < levels.min() <= levels.max() < thresholds[1]["energy"]
    )
    assert len(levels) == 24
    assert thresholds[3]["energy"] - levels.max() > 0.00151
    assert packet["full_generator_residual_frobenius"] < 5e-12
    assert packet["high_address_binding_floor"] > 0.9936999
    for s in thresholds:
        assert s["schur_error_bound"] < s["clearance"]
