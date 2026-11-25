"""Independent energy, angular-support and line-count checks on raw clock exits."""

import json
from math import asin, log, pi, sqrt

import numpy as np
from scipy.linalg import toeplitz

from .measure import DELTA, HERE, SHOTS, carrier, chart, geometry

# Fewer than 5000 binomial categories/clock quadratures, including the
# second receiver experiment. The probability allowance is half the sign's.
SIGN_ERROR = sqrt(2 * log(2 * 5000 / 0.01) / SHOTS) + 2 / 2**25
COMPLEX_ERROR = sqrt(2) * SIGN_ERROR


def sign(bins):
    assert len(bins) == 2 and sum(bins) == SHOTS
    return (bins[0] - bins[1]) / SHOTS


def response(row):
    return sign(row["cos"]) - 1j * sign(row["sin"])


def energy(row):
    return (8 * sign(row["reads"][0]["sin"]) - sign(row["reads"][1]["sin"])) / (
        6 * DELTA
    )


def symbol(axis, gradient, support):
    pts, axes = chart()
    if support == "edges":
        i, j = np.array(carrier.CELL_SEAMS).T
        rays = pts[i] - pts[j]
    else:
        rays = geometry()
    return 2.5 * np.mean(1 - np.cos(gradient * (rays @ axes[axis])))


def test_native_exchange_scale_and_dirichlet_measure():
    d = json.loads((HERE / "records.json").read_text())
    estimates = []
    for t, bins in zip([pi / 12, pi / 3], d["calibration"], strict=True):
        p = bins[1] / SHOTS
        assert sum(bins) == SHOTS
        estimates.append(
            (
                2 * asin(sqrt(p - SIGN_ERROR / 2)) / t,
                2 * asin(sqrt(p + SIGN_ERROR / 2)) / t,
            )
        )
    assert max(x for x, _ in estimates) < 1 < min(y for _, y in estimates)
    for row in d["local"]:
        assert len(row["edges"]) == 30
        assert all(sum(e) == SHOTS for e in row["edges"])
        observed = sum(e[1] for e in row["edges"]) / SHOTS
        expected = symbol(row["axis"], row["q"], "edges")
        assert abs(observed - expected) <= 30 * SIGN_ERROR / 2


def test_independent_dynamical_energy_and_rejected_angular_support():
    d = json.loads((HERE / "records.json").read_text())
    bu = 9 * SIGN_ERROR / (6 * DELTA)
    energies = {
        (r["board"], r["owner"], r["scale"], r["axis"], r["q"]): energy(r)
        for r in d["energy"]
    }
    for (board, owner, scale, axis, q), raw in energies.items():
        actual = raw - energies[board, owner, scale, axis, 0]
        radius = 2 * (bu + (6 * scale) ** 5 * DELTA**4 / 30)
        expected = scale * symbol(axis, q, "edges")
        assert abs(actual - expected) < radius
        if q == 4 and scale == 1:
            assert abs(actual - symbol(axis, q, "ports")) > radius
        if q == 4 and scale == 3:
            assert abs(actual - symbol(axis, q, "edges")) > radius


def test_rank_from_free_flight_quadratures():
    d = json.loads((HERE / "records.json").read_text())
    # Every row of the 4x4 Hermitian error has three off-diagonals.
    radius = 3 * COMPLEX_ERROR
    for b in d["frequency"]:
        vals = np.linalg.eigvalsh(toeplitz([1, *(response(x) for x in b["reads"])]))
        assert vals[0] >= -radius
        if b["board"] == "cell" and b["q"] == 0:
            assert vals[-2] <= radius
        if b["axis"] == "port" and b["q"] == 3:
            if b["board"] == "cell" and b["transverse"]:
                assert vals[-3] <= radius < vals[-2]
            else:
                assert vals[0] > radius


def test_transverse_template_lines_and_weights():
    d = json.loads((HERE / "transverse.json").read_text())
    points = geometry()
    i, j = carrier.CELL_SEAMS[0]
    latitude = 1 / (sqrt(5) * np.linalg.norm(points[i] - points[j]))
    for bank in d["banks"]:
        k, scale = bank["q"], bank["scale"]
        assert bank["holds"] == [0, 1 / 6, 1 / 2]
        for receiver, (reads, weight, level) in enumerate(
            zip(
                bank["receivers"],
                [np.cos(k * latitude) ** 2, np.sin(k * latitude) ** 2],
                [(5 - sqrt(5)) / 2, 3],
                strict=True,
            )
        ):
            amps = [response(r) for r in reads]
            assert abs(abs(amps[0]) - sqrt(weight)) < COMPLEX_ERROR
            if k == 0 and receiver == 1:
                assert all(abs(z) < COMPLEX_ERROR for z in amps)
                continue
            minimum = min(abs(z) for z in amps)
            # Angle to any exact point in the complex confidence disk.
            angle = 2 * asin(COMPLEX_ERROR / minimum)
            pilot = -np.angle(amps[1] / amps[0])
            final = -np.angle(amps[2] / amps[0])
            assert (6 * scale) / 6 + angle < pi
            assert 4 * angle < pi
            final += 2 * pi * round((3 * pilot - final) / (2 * pi))
            assert abs(2 * final - scale * level) < 2 * angle
