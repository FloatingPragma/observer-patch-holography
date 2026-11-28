"""Read the same joint occupation with the retained and exchanged thermal references."""

import json

import numpy as np

from . import sector_ball_source as source
from .sector_passage_source import mul


def compare(words, prepared, held, lower):
    banks = json.loads((source.balls.HERE / "sector_calorimetry.json").read_text())[
        "ordinary_clock_calibration"
    ]
    calibrated = source.sector_regional.calibrate(banks)
    bath = source.thermal.bath()
    cell = 2 * (words[:, 0] > 0) + (words[:, 3] > 0)
    values = {}
    probs = {}
    intervals = {}
    for tag in ("retained", "exchanged"):
        z = np.array([calibrated[tag][cell == i].sum(0) for i in range(4)])
        probs[tag] = np.bincount(cell, weights=bath[tag].diagonal())
        odds = np.log(z[0] * z[3] / (z[1] * z[2])[::-1]) + np.array([-1e-11, 1e-11])
        intervals[tag] = odds
        values[tag] = {"probabilities": z.tolist(), "log_odds": odds.tolist()}
    coeff = intervals["exchanged"] - intervals["retained"][::-1]
    assert coeff[1] < 0
    K = np.log(probs["exchanged"] / probs["retained"])
    exact_K = float(K[3] - K[2] - K[1] + K[0])
    assert coeff[0] < exact_K < coeff[1]
    rates = {}
    for k, (d, eps) in prepared.items():
        gain = np.bincount(cell, weights=d, minlength=4)
        base = np.bincount(cell) / len(words)

        def cov(v):
            return v[3] - (v[2] + v[3]) * (v[1] + v[3])

        c = cov(base + gain) - cov(base)
        rates[str(k)] = {
            "total_relative_modular_nats": float(gain @ K),
            "connected_relative_modular_nats": float(exact_K * c),
            "distribution_error_l1": eps,
        }
    return {
        "references": values,
        "mixed_price": coeff.tolist(),
        "mixed_price_at_native_spectra": exact_K,
        "control_connected_relative_nats": (coeff / 2).tolist(),
        "full_peer_at_eighth": mul(coeff, held["covariance_interval"]).tolist(),
        "negative_connected_nats_per_rate_at_least": (-coeff[1]) * lower - 1e-10,
        "rates": rates,
    }
