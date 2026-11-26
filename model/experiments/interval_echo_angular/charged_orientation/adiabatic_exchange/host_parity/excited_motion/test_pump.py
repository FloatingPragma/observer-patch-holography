import gzip
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import pump
from .measure import fixed_guide


def saved():
    return json.loads(
        gzip.decompress(Path(__file__).with_name("records.json.gz").read_bytes())
    )


def test_native_pair_gate_and_exact_same_source_control():
    assert pump.route_certificate()["pair_gate_columns"] == 16
    for row in fixed_guide().values():
        assert row["dark"][0] < 1e-13 and row["dark"][1] < 0.0013


def test_reflection_and_certified_exchange_sign():
    words = pump.stocks(2)
    reflected = [words.index((w[0], w[2], w[1], w[3])) for w in words]
    for sign, s in ((16, -1), (-16, 1)):
        for theta, want in ((0, 1), (np.pi, s)):
            f, q = pump.clocks(2, sign, theta, 1 / 16)
            line = q[:, 3 if sign > 0 else 0]
            np.testing.assert_allclose(line[reflected], want * line, atol=1e-11)
        case = next(row for row in saved()["held"] if row["detuning"] == sign)
        assert case["min_sampled_gap"] - 4 * np.pi / case["appointments"] > 0.0229
        assert case["polygon_sign"] == s


def test_deposited_probabilities_and_shot_denominators():
    packet = saved()
    for category in ("early", "held", "guarded", "parked"):
        for case in packet[category]:
            for row in case["records"]:
                for name in ("once", "twice"):
                    assert all(
                        sum(bank) == case["shots"]
                        for bank in row[name]["joint_launch_counts"]
                    )
                for arm in ("forward", "reference"):
                    for frame in row[arm + "_frame_counts"]:
                        assert all(sum(bank) == case["shots"] for bank in frame)
                if category == "held" and row["duration"] == 393216:
                    assert min(row["min_pinned"]) > 0.98
                    assert min(row["min_crossing_shared_two_pairs"]) > 0.98
                    assert max(row["twice"]["dark"]) < 0.02
                    if case["detuning"] > 0:
                        assert min(row["once"]["dark"]) > 0.99
                    else:
                        assert max(row["once"]["dark"]) < 0.01
                if category == "parked":
                    assert all(
                        sum(bank) == case["shots"] for bank in row["joint_route_counts"]
                    )


def test_step_power_independent_of_spectral_acquisition():
    g = 16
    dt = 1 / 16
    h, pairs, centre = pump.source_geometry(2)
    for theta in (0.23, 1.14, 3.9):
        d = g * centre - pairs @ pump.aperture(theta)
        kick = np.diag(np.exp(-0.5j * dt * d))
        u = kick @ expm(-1j * dt * h) @ kick
        e, q = pump.clocks(2, g, theta, dt)
        np.testing.assert_allclose(
            (q * np.exp(-7j * dt * e)) @ q.T, np.linalg.matrix_power(u, 7), atol=3e-14
        )
