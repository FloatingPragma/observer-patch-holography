"""Predict regional modular readings from real-time ordinary clocks and native load counts."""

import math

import numpy as np

from . import sector_thermal as thermal


def calibrate(banks):
    w, _, _, blocks, _ = thermal.apparatus()
    # Independent of the source. Each recorded start has its own visibility
    # and Fourier tail enclosure, giving its diagonal Gibbs weight.
    sharp = {k: np.zeros((len(w), 2)) for k in ("retained", "exchanged")}
    for row, g in zip(banks, blocks, strict=True):
        n = np.array(row["sharp_starts"])
        z = np.array(row["sharp_signed_counts"]) / n
        err = (
            np.sqrt(2 * math.log(2 * len(w) / 0.001) / n)
            + row["omitted_visibility_bound"]
            + 1e-9
        )
        lo, hi = z - err, z + err
        assert min(lo) > 0
        # Normalize each sector for retained records and once for a bath
        # which exchanges the sector label.
        sharp["retained"][g] = (
            np.column_stack((lo / hi.sum(), hi / lo.sum())) * len(g) / len(w)
        )
        sharp["exchanged"][g] = np.column_stack((lo, hi)) * math.exp(
            -thermal.BETA * row["energy_window"][0]
        )
    full = sharp["exchanged"].copy()
    sharp["exchanged"] /= full.sum(0)[::-1]
    return sharp


def measure(members, banks):
    w = thermal.apparatus()[0]
    reference = calibrate(banks)
    out = {}
    for j in range(4):
        var = (w[:, j] - 1) ** 2
        observed = var[members].mean(1)
        delta = observed.mean() - var.mean()
        h = math.log(2 * 4 / 0.001)
        err = math.sqrt(2 * observed.var(ddof=1) * h / len(members)) + 7 * h / (
            3 * (len(members) - 1)
        )
        receivers = {}
        for name, weights in reference.items():
            probs = np.array([weights[w[:, j] == n].sum(0) for n in range(3)])
            # At fixed mean load, each of zero and double frequency
            # changes by half the increment in the observed squared load.
            beta = 0.5 * np.log(probs[1] ** 2 / (probs[0, ::-1] * probs[2, ::-1]))
            extrema = np.outer(beta, [delta - err, delta + err])
            receivers[name] = dict(
                ordinary_reference_probability=probs.tolist(),
                modular_per_squared_reading=beta.tolist(),
                inferred_nats=[extrema.min(), extrema.max()],
            )
        out[str(j)] = dict(
            squared_load_count=int(var[members].sum()),
            variance_increment=[delta - err, delta + err],
            receivers=receivers,
        )
    return out
