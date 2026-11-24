"""Independent branch history, exterior support, and held-receipt checks."""

import json

import numpy as np
from scipy.linalg import expm

from . import experiment as base
from . import exterior as e


def data():
    return json.loads((e.HERE / "exterior.json").read_text())


def test_native_separation_and_exact_static_support():
    grid = e.chart()
    assert set(grid["radar_belt"]).isdisjoint(grid["mobile_belt"] + [grid["pole"]])
    assert min(grid["native_distances_pole_to_radar"]) == 3
    assert e.local_occupation_cost(8)[1, 1] == -32
    assert abs(e.local_occupation_cost(8)[1, 2]) < 1e-12
    for arm in ("empty_write", "coin_even"):
        for tag in (1, 2):
            assert e.guards(arm, 8, tag) < 1e-12


def path_sum(u, prepared, last=16):
    # Enumerate pairs of first-click times with a coherent unmeasured mediator.
    # This calculation uses amplitudes for each history, independently of both
    # the stored mixture and the continuously updated two-stage density matrix.
    n = len(prepared)
    far = np.arange(5, n, 10)
    home = np.arange(0, n, 10)
    total = 0.0
    for j in range(1, last):
        prepared = u @ prepared
        relayed = np.zeros_like(prepared)
        relayed[far] = prepared[far]
        prepared[far] = 0
        for _ in range(last - j):
            relayed = u @ relayed
            total += np.vdot(relayed[home], relayed[home]).real
            relayed[home] = 0
    return total


def test_independent_full_joint_held_and_live_histories():
    r = data()
    adj = -(np.roll(np.eye(10), 1, axis=1) + np.roll(np.eye(10), -1, axis=1)) / 2
    for live in (False, True):
        for tag in (1, 2):
            for occ in (0, 1):
                med_h = adj.copy()
                med_h[5, 5] -= 32 * occ
                med = expm(-12j * med_h)[:, 0]
                joint = np.kron(np.eye(10), adj)
                if live:
                    joint += np.kron(med_h, np.eye(10))
                joint[2, 2] -= 32 * (tag == 1)
                prob = path_sum(expm(-1j * joint), np.kron(med, np.eye(10)[0]))
                row = next(
                    b
                    for b in r["continuous" if live else "stored"]
                    if b["q"] == tag
                    and b["source"] == occ
                    and b["launch"] == 0
                    and b["k"] == 8
                    and b.get("hold", 12) == 12
                )
                assert abs(prob - row["cdf" if live else "measured"][1]) < e.EPS


def test_both_launches_controls_and_held_exterior_forecast():
    r = data()
    for row in r["stored"]:
        assert (
            max(row["clock_forecast_max_error"], row["radar_forecast_max_error"])
            < 3 * e.EPS
        )
    for row in r["continuous"]:
        if row["source"] == 0:
            continue
        blank = next(
            b
            for b in r["continuous"]
            if not b["source"] and all(row[k] == b[k] for k in ("launch", "q", "k"))
        )
        change = np.array(row["cdf"]) - blank["cdf"]
        if row["k"] != 0 and row["q"] == 1:
            if row["launch"] == 0:
                assert change[-1] < -0.03
            else:
                assert change[-1] > 0.10
        else:
            assert abs(change).max() < 2 * e.EPS
    with np.load(base.HERE / "exterior.npz") as bank:
        for tag in (1, 2):
            a, b = (bank[f"k8.0_q{tag}_empty_mediator_occ{i}"] for i in (0, 1))
            assert max(abs(a.cumsum() - b.cumsum())) / e.SHOTS < 2 * e.EPS
