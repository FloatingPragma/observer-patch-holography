"""Replicate the sign clock and distinguish waveform transfer from rate transfer."""

import json

import numpy as np

from . import flight
from . import signal as s


def test_true_sharp_clock_and_two_independent_material_standards():
    run = s.replay(json.loads((s.HERE / "results.json").read_text()))
    assert abs(run["energy"]["admission"] - 2 / 33) < 1e-15
    assert run["energy"]["full_ordinary_cost"] == 561 / 16
    assert run["direct_check_max_error"] < 1e-9
    for name in s.NAMES:
        for hand in s.HANDS:
            times = run["register_crossings"][name][hand]
            assert len(times) == 3
            for i, t in enumerate(times):
                low, high = run["pointer_counts"]["records"][f"{name}/{hand}/{i + 1}"][
                    "bracket"
                ]
                assert low < t < high
    half = np.pi * np.array([0.5, 1.5, 2.5])
    np.testing.assert_allclose(
        run["register_crossings"]["empty"]["near"], half, rtol=0, atol=1e-11
    )
    d = run["delay_contrasts"]
    assert d["mass"][0]["interval"][0] > 0 > d["erasure"][0]["interval"][1]
    assert d["phase_feedback"][0]["interval"][1] < 0
    assert d["erasure"][1]["interval"][1] < 0 < d["phase_feedback"][1]["interval"][0]
    for name in ("erasure", "phase_feedback"):
        assert max(x["ratio_interval"][0] for x in d[name]) > min(
            x["ratio_interval"][1] for x in d[name]
        )
    vals = np.array(run["scan"]["signed_return_by_hand"])
    cols = list(s.NAMES)
    np.testing.assert_allclose(
        vals[:, :, cols.index("dephased")] - vals[:, :, cols.index("orbit")],
        vals[:, :, cols.index("feedback")] - vals[:, :, cols.index("paired")],
        atol=1e-11,
    )


def test_same_clock_equipped_material_flight_at_held_time():
    result = flight.run()
    s.m.fp.check(
        json.loads((s.HERE / "flight.json").read_text()), json.loads(json.dumps(result))
    )
    assert result["full_dimensions"] == [64, 1536, 17664]
    assert result["radius"]["mass"]["mean"] < -0.02
    assert result["radius"]["erasure"]["mean"] > 0.086
    assert result["radius"]["phase_feedback"]["mean"] > 0.029
    assert result["radius"]["erasure"]["mirror_launches"][0] < 0
    assert result["radius"]["erasure"]["mirror_launches"][1] > 0
