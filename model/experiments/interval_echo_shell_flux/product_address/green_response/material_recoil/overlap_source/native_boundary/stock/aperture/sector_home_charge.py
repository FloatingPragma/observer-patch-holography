"""Infer ordinary rest energy using only two coarse modular readings and bath calibrations."""

from itertools import product

import numpy as np


def reconstruct(receipt):
    witness = receipt["fixed_volume_test"]["bracket_ball"]
    calibration = receipt["occupied_home_radar"]["six_marker_binary_source"]
    c = np.array(calibration["log_price_difference_from_exchanged_reference"])

    def odds(key):
        z = np.asarray(witness[key])
        return np.log(z[1] / z[[0, 2], ::-1])

    a, b = (
        odds("conditional_clock_probabilities")
        - odds("exchanged_clock_probabilities")[:, ::-1]
    )
    assert min(c) > 0 and min(b) > 0 and min(a) > 0
    rates = {}
    for k, row in receipt["held_native"].items():
        binary = calibration["measured_at_rates"][k]
        ternary = row["home_relative_modular_nats"]
        # The expression is multi-affine in the readings, a, 1/b and 1/c.
        # Corner evaluation bounds correlations between instruments as well.
        inferred = []
        for qb, qh, ca, pa, pb in product(binary, ternary, c, a, b):
            x = qb / ca
            y = (qh - pa * x) / pb
            inferred.append([48 * (-x + 3 * y), -2 * x])
        v = np.array(inferred)
        rates[k] = {
            "binary_charge": binary,
            "ternary_charge": ternary,
            "rest_energy": [v[:, 0].min() - 1e-10, v[:, 0].max() + 1e-10],
            "echo_ticks": [v[:, 1].min() - 1e-10, v[:, 1].max() + 1e-10],
        }
    return {
        "binary_coefficient": c.tolist(),
        "ternary_coefficients": [a.tolist(), b.tolist()],
        "inferred": rates,
    }
