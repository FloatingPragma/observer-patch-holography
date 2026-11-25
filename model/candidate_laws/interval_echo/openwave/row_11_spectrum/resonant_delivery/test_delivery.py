import json
from fractions import Fraction as Q
from itertools import product

import numpy as np
from scipy.linalg import eigh

from . import measure as delivery


def test_entire_native_two_dimer_restriction_and_opposite_current():
    h, _ = delivery.reduced((0, 1), "native", 8)
    expected = np.array(
        [[50, -1, -1, 0], [-1, 130, 0, -1], [-1, 0, 130, -1], [0, -1, -1, 146]]
    )
    np.testing.assert_array_equal(2 * h, expected)
    swap = np.eye(4)[[0, 2, 1, 3]]
    np.testing.assert_array_equal(h @ swap, swap @ h)
    odd = np.array([0, 1, -1, 0])
    np.testing.assert_array_equal(h @ odd, 65 * odd)
    hc, _ = delivery.reduced((0,), "native", 8)
    hn, _ = delivery.reduced((1,), "native", 8)
    h0 = np.kron(hc, np.eye(2)) + np.kron(np.eye(2), hn)
    levels, a = eigh(hc)
    _, b = eigh(hn)
    spectral = np.kron(a, b)
    s, d = spectral[:, 2], spectral[:, 1]
    np.testing.assert_allclose(
        h0 @ np.stack((s, d)).T, 65 * np.stack((s, d)).T, atol=1e-12
    )
    np.testing.assert_allclose(abs(d @ (h - h0) @ s), 8 / 1601, atol=1e-12)
    local = np.kron(np.eye(2), hn)
    j = 1j * (h @ local - local @ h)
    np.testing.assert_allclose(abs(d @ j @ s), 8 / np.sqrt(1601), atol=1e-12)
    for sign in (1, -1):
        packet = (s + sign * 1j * d) / np.sqrt(2)
        np.testing.assert_allclose(
            abs(np.vdot(packet, j @ packet)), 8 / np.sqrt(1601), atol=2e-13
        )
    assert np.isclose(levels[1] - levels[0], np.sqrt(1601))


def test_receipt_and_full_energy_spectral_replay():
    archive = json.loads((delivery.HERE / "results.json").read_text())
    assert archive == delivery.run()
    for label, k in (
        ("native", 8),
        ("detuned", 8),
        ("matched_cut", 8),
        ("half_contrast", 4),
    ):
        setting = archive["settings"][label]
        reg = label if label in ("detuned", "matched_cut") else "native"
        h, _ = delivery.reduced((0, 1), reg, k)
        e, v = eigh(h)
        a = eigh(delivery.reduced((0,), reg, k)[0])[1]
        b = eigh(delivery.reduced((1,), reg, k)[0])[1]
        frame = np.kron(a, b)
        for row in setting["joint_times"]:
            i = 2 if row["initial"] == "charged_excited" else 1
            amplitudes = (frame.T @ v) @ (
                np.exp(-1j * e * row["exposure"]) * (v.T @ frame[:, i])
            )
            p = abs(amplitudes) ** 2
            np.testing.assert_allclose(
                p, row["energy_probability_gg_ge_eg_ee"], atol=2e-10
            )
            assert (
                sum(row["energy_counts_gg_ge_eg_ee"])
                == sum(row["both_gateway_vacant_other"])
                == delivery.SHOTS
            )
            if label == "native" and i == 2 and row["exposure"] == 64:
                cert = setting["certified_delivery"]
                lo, hi = [
                    Q(x, cert["probability_tick_denominator"])
                    for x in cert["joint_probability_interval_ticks"]
                ]
                assert (
                    float(lo) <= p[1] <= float(hi)
                    and lo > Q(99, 100)
                    and hi - lo <= Q(2, 10**12)
                )
        for name, levels in (
            ("charged", "local_charged_levels"),
            ("neutral", "local_neutral_levels"),
        ):
            clock = setting[name + "_clock"]
            assert all(
                sum(bank) == delivery.CAL_SHOTS for bank in clock["bright_other"]
            )
            assert abs(clock["gap_from_counts"] - np.diff(setting[levels])[0]) < 0.035
    b, arms, _, park = delivery.fixture()
    for i, j in product(range(2), repeat=2):
        x = np.ones(b.n, dtype=int)
        x[[arms[0][i], arms[1][j]]] = 0
        x[park] = 2
        assert sum(x[:12] - 1) == -1 and sum(x[72:84] - 1) == 0
