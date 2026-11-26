"""Check the physical energy meter, causal controls and all recorded original starts."""

import json
import numpy as np
from scipy.linalg import expm
from . import calorimeter, power
from ..material import route
from ..measure import rebuild


def records():
    return json.loads((calorimeter.HERE / "records.json").read_text())


def test_detector_and_upstream_cuts_factorize():
    board, dimers, heads, emitter, distances = calorimeter.apparatus()
    assert dimers[-1] == (77, 81) and distances == [4, 5]
    original = rebuild(board, dimers[:-1], route((3,), True)[2], False, True)
    arms = calorimeter.generators(board, dimers, heads, emitter, original)
    for name, _, h in arms:
        if name == "native":
            continue
        n = 4 if name == "upstream_cut" else 2
        # For every exterior operation the missing bond separates both generators.
        inner = h[::n, ::n]
        outer = h[:n, :n] - np.eye(n) * h[0, 0]
        np.testing.assert_allclose(
            h, np.kron(inner, np.eye(n)) + np.kron(np.eye(len(inner)), outer)
        )
        np.testing.assert_allclose(
            expm(-0.125j * h),
            np.kron(expm(-0.125j * inner), expm(-0.125j * outer)),
            atol=1e-13,
        )


def test_phase_receiver_measures_the_native_gap():
    meter = calorimeter.native_meter()
    clock = records()["meter_clock"]
    assert abs(clock["gap"] - np.ptp(np.linalg.eigvalsh(meter))) < clock["half_width"]
    # This is the same detector operating on an altered state evolution.
    changed = calorimeter.gap_calibration(np.random.default_rng(17), 1.5 * meter)
    assert abs(changed["gap"] - 1.5) < changed["half_width"]
    labels = records()["meter_charge"][0]["word_00_to_22"]
    assert [labels[i] for i in (0, 1, 3, 4)] == [7, 0, 0, 1]


def test_enrolled_heating_and_joint_comparisons():
    data = records()
    radius = data["familywise_99_half_width"]
    shots = data["shots_per_panel"]
    hot = lambda s, a, t: calorimeter.decode(data, s, a, t)[0] / shots
    for t in calorimeter.WAITS:
        for arm in ("pulse", "idle"):
            assert hot("detector_cut", arm, t) == hot("no_reversal", arm, t) == 0
            if t == 0:
                assert hot("native", arm, t) == 0
        assert (
            abs(hot("upstream_cut", "pulse", t) - hot("upstream_cut", "idle", t))
            < 2 * radius
        )
    assert hot("native", "pulse", 32) - hot("native", "idle", 32) > 0.01 + 2 * radius
    joint = calorimeter.decode(data, "native", "pulse", 32)[1] / shots
    controls = [
        hot(r["apparatus"], r["source"], 32)
        for r in data["rows"]
        if r["wait"] == 32 and (r["apparatus"], r["source"]) != ("native", "pulse")
    ]
    assert joint - max(controls) > 2 * radius


def test_power_counts_recover_endpoint_energy():
    data = json.loads((calorimeter.HERE / "power_counts.json").read_text())
    endpoints = records()
    for p in data["panels"]:
        hot, _ = calorimeter.decode(endpoints, p["apparatus"], p["source"], 32)
        error = (
            data["sampling_half_width"]
            + p["quadrature_bound"]
            + endpoints["familywise_99_half_width"]
        )
        assert abs(power.decode(data, p) - hot / endpoints["shots_per_panel"]) < error
    rebuilt = power.acquire()
    for actual, saved in zip(rebuilt["panels"], data["panels"], strict=True):
        assert (actual["apparatus"], actual["source"], actual["counts"]) == (
            saved["apparatus"],
            saved["source"],
            saved["counts"],
        )
        np.testing.assert_allclose(
            actual["quadrature_bound"], saved["quadrature_bound"], atol=1e-12
        )
    baseline = calorimeter.acquire()
    # Receipts contain measured decimal energies. Compare discrete readouts exactly.
    assert baseline["rows"] == endpoints["rows"]
    assert baseline["meter_clock"]["holds"] == endpoints["meter_clock"]["holds"]
