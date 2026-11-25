"""Independent finite-chain predictions and confidence checks on the stored records."""

from fractions import Fraction
import json
from math import log, sqrt

import numpy as np

from .measure import HERE


def uncertainty(n):
    # Fewer than 1024 scalar expectations enter the simultaneous comparisons.
    return sqrt(log(2 * 1024 / 0.01) / (2 * n))


def contact(size, left, right, cut):
    menu = 2 ** max(2, size.bit_length())
    rates = np.ones(size + 1) / menu
    rates[0] *= (1 + left) / 2
    rates[-1] *= 0 if cut else (1 + right) / 2
    middle = np.diag(rates[1:-1], k=1)
    laplace = np.diag(rates[:-1] + rates[1:]) - middle - middle.T
    trans = np.eye(size) - laplace
    conductance = 0 if cut else 1 / sum(1 / rates)
    return rates, trans, conductance


def readings():
    return json.loads((HERE / "records.json").read_text())


def test_complete_records_and_separate_thermometers():
    data = readings()
    n = data["calibration_proposals"]
    for setting, bins in data["thermometer"].items():
        w = float(Fraction(setting))
        assert sum(bins) == n
        np.testing.assert_allclose(
            np.array(bins) / n, [0.5, w / 2, (1 - w) / 2], atol=uncertainty(n), rtol=0
        )
    shots = data["paths"]
    for bank in data["driven"] + data["equilibrium"]:
        count = data["burn_in"] + 1 if "current" in bank else data["lags"] + 1
        assert sum(bank["all_offers"].values()) == count * shots
        assert len(bank["initial_refusals"]) == bank["size"]
        assert all(r > 0 for r in bank["initial_refusals"])
        assert bank["all_offers"]["left_refused"] > 0
        assert bank["all_offers"]["right_refused"] > 0
        if "current" in bank:
            counts = np.array(bank["final_offers_by_word"])
            np.testing.assert_array_equal(
                counts.sum(axis=1), bank["counts_before_final"]
            )
            # An accepted upper source gains the old endpoint sign; observer
            # energy lost is +1 exactly on old lower / accepted upper.
            plus, minus = counts[::2, 2].sum(), counts[1::2, 1].sum()
            assert bank["current"] == [minus, shots - plus - minus, plus]
        else:
            census = np.array(bank["joints"])
            np.testing.assert_array_equal(
                census.sum(axis=(1, 2)), np.full(count, shots)
            )
            np.testing.assert_array_equal(
                census.sum(axis=2), np.tile(census[0].sum(axis=1), (count, 1))
            )


def test_exact_green_kubo_and_equilibrium_current_samples():
    d = readings()
    signs = np.arange(-1, 2)
    for bank in d["equilibrium"]:
        n = bank["size"]
        rates, step, conductance = contact(n, 0.5, 0.5, bank["cut"])
        c, v = rates[0], 2 / 9
        inv = np.linalg.inv(np.eye(n) - step)
        assert abs(c - c**2 * inv[0, 0] - conductance) < 1e-15
        k = d["lags"]
        true = np.array(
            [2 * c * v]
            + [-(c**2) * v * np.linalg.matrix_power(step, t)[0, 0] for t in range(k)]
        )
        sample = signs @ np.array(bank["joints"]) @ signs / d["paths"]
        radius = uncertainty(d["paths"])
        np.testing.assert_allclose(sample, true, atol=2 * radius, rtol=0)
        assert sample[0] > radius and sample[1] < -2 * radius
        tail = c**2 * v * (np.linalg.matrix_power(step, k) @ inv)[0, 0]
        summed = sample[0] / 2 + sample[1:].sum() - tail
        assert abs(summed - conductance * v) < (2 * k + 0.5) * radius


def test_independent_imbalance_and_finite_bias_comparison():
    d = readings()
    signs = np.arange(-1, 2)
    eps = uncertainty(d["calibration_proposals"])
    for plus, minus in zip(d["driven"][::2], d["driven"][1::2], strict=True):
        assert (plus["left"], plus["right"]) == (minus["right"], minus["left"])
        measured = []
        exact = []
        transient_bounds = []
        for bank in (plus, minus):
            left, right = [float(Fraction(bank[key])) for key in ("left", "right")]
            rates, step, g = contact(bank["size"], left, right, bank["cut"])
            rho = max(abs(np.linalg.eigvalsh(step)))
            transience = rates[0] * sqrt(bank["size"]) * rho ** d["burn_in"]
            transient_bounds.append(transience)
            measured.append(signs @ bank["current"] / d["paths"])
            exact.append(g * (left / (1 + left) - right / (1 + right)))
            assert (
                abs(measured[-1] - exact[-1]) < 2 * uncertainty(d["paths"]) + transience
            )
        odds = []
        for key in (plus["left"], plus["right"]):
            low, high, _ = np.array(d["thermometer"][key]) / d["calibration_proposals"]
            odds.append(((high - eps) / (low + eps), (high + eps) / (low - eps)))
        aff_min, aff_max = log(odds[0][0] / odds[1][1]), log(odds[0][1] / odds[1][0])
        affinity = log(float(Fraction(plus["left"]) / Fraction(plus["right"])))
        assert 0 < aff_min < affinity < aff_max
        # Infer the secant from two independent heat readings and independent
        # thermometers; subtract only the predicted, stated finite-bias shift.
        transfer = (measured[0] - measured[1]) / 2
        delta = 2 * uncertainty(d["paths"]) + max(transient_bounds)
        interval = [
            (transfer + sign * delta) / a
            for sign in (-1, 1)
            for a in (aff_min, aff_max)
        ]
        equilibrium = contact(plus["size"], 0.5, 0.5, plus["cut"])[2] * 2 / 9
        correction = (exact[0] - exact[1]) / (2 * affinity) - equilibrium
        assert min(interval) <= equilibrium + correction <= max(interval)
        if not plus["cut"]:
            assert measured[1] < -delta < delta < measured[0]
            assert min(interval) - correction > 0


def test_diffusion_and_independent_four_register_conductance():
    r = json.loads((HERE / "ring.json").read_text())
    n = r["trials"]
    eps = uncertainty(n)
    signs = np.arange(-1, 2)
    for bank in r["banks"]:
        assert all(
            bins[int(site == 0)] == n
            for site, bins in enumerate(bank["original_proposals"])
        )
        predicted_d = 0 if bank["parked"] else 1 / 8
        for t, row in enumerate(bank["steps"], start=1):
            joint = np.array(row["first_and_current"])
            bins = np.array(row["displacements"])
            x = np.arange(-t, t + 1)
            assert sum(row["offers"]) == sum(bins) == joint.sum() == n
            expected = 2 * predicted_d if t == 1 else 0
            assert abs(signs @ joint @ signs / n - expected) < 2 * eps
            assert abs(x @ bins / n) <= 2 * t * eps
            measured_d = (x**2 @ bins) / (2 * t * n)
            assert abs(measured_d - predicted_d) < t * eps / 2
            if bank["parked"]:
                assert bins[t] == n
        if bank["parked"]:
            continue
        d = readings()
        low, high, _ = np.array(d["thermometer"]["1/2"]) / d["calibration_proposals"]
        err = uncertainty(d["calibration_proposals"])
        predictions = []
        for sign in (-1, 1):
            p = (high + sign * err) / (low - sign * err + high + sign * err)
            admitted = (low + high + sign * 2 * err) / 8
            diffusion = measured_d + sign * r["ticks"] * eps / 2
            predictions.append(p * (1 - p) / (2 / admitted + 3 / diffusion))
        # Independent drive pair at the smaller thermal bias.
        plus, minus = [
            x
            for x in d["driven"]
            if x["size"] == 4 and not x["cut"] and "9/16" in (x["left"], x["right"])
        ]
        wp, wm = [float(Fraction(plus[key])) for key in ("left", "right")]
        lo, hi = [], []
        for b in (plus["left"], plus["right"]):
            lower, upper, _ = np.array(d["thermometer"][b]) / d["calibration_proposals"]
            lo.append((upper - err) / (lower + err))
            hi.append((upper + err) / (lower - err))
        aff_range = (log(lo[0] / hi[1]), log(hi[0] / lo[1]))
        j = signs @ (np.array(plus["current"]) - minus["current"]) / (2 * d["paths"])
        rates, step, g = contact(4, wp, wm, False)
        delta = (
            2 * uncertainty(d["paths"])
            + 2 * rates[0] * max(abs(np.linalg.eigvalsh(step))) ** d["burn_in"]
        )
        bias = g * (wp / (1 + wp) - wm / (1 + wm)) / log(wp / wm) - 1 / 204
        observed = (
            (j - delta) / aff_range[1] - bias,
            (j + delta) / aff_range[0] - bias,
        )
        assert max(observed[0], predictions[0]) <= min(observed[1], predictions[1])
