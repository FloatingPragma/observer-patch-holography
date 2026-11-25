import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from . import measure
from .certify import DENOM, below


DATA = json.loads(Path(__file__).with_name("results.json").read_text())


def test_physical_head_columns_against_open_laplacian():
    assert measure.SPACINGS == (1, 2, 3, 4, 6, 8, 12, 16)
    for s in measure.SPACINGS:
        r = DATA["settings"][str(s)]
        for size, row in r["finite"].items():
            n = int(size)
            x = np.arange(-1, n + 1)
            # Each opening ends at a blank contact whose neighboring active
            # contact is blocked. Positive endpoint correction versus bulk is 4.
            d = 4.5 + 32 * ((x - 1) % (3 * s) == 0).astype(float)
            d[[0, -1]] += 4
            e = np.where(x[:-1] % 2, -0.5, -4.0)
            for sign in (-1, 1):
                actual = measure.finite(s, n, sign)
                np.testing.assert_array_equal(actual[0], d)
                np.testing.assert_array_equal(actual[1], e)
            lo, hi = row["ground_interval_numerator"]
            a, b = (2 * d).astype(int).tolist(), (2 * e).astype(int).tolist()
            assert below(a, b, 0, lo) == 0 and below(a, b, 0, hi) >= 1
        # Rebuild at interior Floquet angles; check all ordered eigenvalues
        # against the certified bands at the extremal angles.
        p = r["period"]
        for theta in np.pi * np.array([0.17, 0.43, 0.79]):
            h = np.zeros((p, p), complex)
            for i in range(p):
                j = (i + 1) % p
                g = np.zeros(p, complex)
                g[i] = 1
                g[j] = -np.exp(1j * theta) if j == 0 else -1
                h += (0.5 if i % 2 else 4) * np.outer(g, g.conjugate())
                if (i - 1) % (3 * s) == 0:
                    h[i, i] += 32
            vals = np.linalg.eigvalsh(h)
            assert np.all(vals >= np.array(r["band_minima"]) - 5.0e-11)
            assert np.all(vals <= np.array(r["band_tops"]) + 5.0e-11)
    assert DATA["settings"]["0"]["finite"]["144"]["minimum"] < 0.00042
    assert DATA["settings"]["0"]["band_minima"][0] == 0


def test_both_laboratory_ratios_and_filter_flight():
    records = [v for k, v in DATA["settings"].items() if k != "0"]
    assert 0.292 < min(v["any_bulk_triple"]["joint_error"] for v in records) < 0.293
    assert 0.292 < min(v["finite_family"]["joint_error"] for v in records) < 0.294
    r = records[-1]
    flight = r["sharp_flight"]
    clock = flight["clock"]
    assert clock["first_band_weight"] > 0.9997
    for row in clock["banks"]:
        assert row["arrived"] + row["elsewhere"] + row["refused"] == measure.SHOTS
    # Multiply the actual controlled binary-wait filters by repeated squaring,
    # then independently obtain the free propagator at the held exposure.
    d, e = measure.finite(16, 144, -1)
    h = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
    delta = h - clock["center"] * np.eye(len(d))
    u = expm(-1j * delta * clock["step"])
    psi = np.eye(len(d), dtype=complex)[:, flight["launch_ring_index"] + 1]
    for _ in range(clock["register_size"].bit_length() - 1):
        psi = (psi + u @ psi) / 2
        u = u @ u
    assert abs(np.vdot(psi, psi).real - clock["admitted"]) < 1.0e-9
    for _ in range(8):
        u = u @ u
    mask = np.isin(np.arange(-1, 145) % 144, flight["receiver_ring_indices"])
    assert (
        abs(sum(abs((u @ psi)[mask]) ** 2) - clock["banks"][-2]["joint_probability"])
        < 2e-8
    )
