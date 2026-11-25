"""Verify both reciprocal-camera count banks independently of the adaptive reader."""

import json

import numpy as np
from scipy.stats import beta

from . import optics


def previous_banks(charge, kernels):
    c = optics.camera
    saved = json.loads((c.HERE / "results.json").read_text())
    x = c.preparation()
    u = np.kron(np.diag([1.0, 1j]), np.eye(2))
    y = u @ x @ u.conj().T

    def tr(blocks):
        return blocks.transpose(0, 2, 1)

    cases = (
        ("plain", x, np.broadcast_to(np.eye(4), x.shape)),
        ("matched", x, tr(x)),
        ("held_fixed", y, tr(x)),
        ("held_recalibrated", y, tr(y)),
    )
    minimum = 1.0
    gains = []
    mismatch = None
    for tag, source, effect in cases:
        forward, reverse = optics.directional(kernels, source, [effect, effect])
        previous = saved["charges"][str(charge)][tag]
        for side, values in (("forward", forward), ("reverse", reverse)):
            n = np.array(previous[side]["counts"])[:, :5].T
            a = 0.01 / (2 * 201)
            low = beta.ppf(a, n, c.SHOTS - n + 1)
            high = beta.ppf(1 - a, n + 1, c.SHOTS - n)
            margin = min(np.min(values - low), np.min(high - values)) - 1e-14
            assert margin > 0
            minimum = min(minimum, margin)
            if tag == "matched":
                g = c.gain_bounds(
                    (values - 1e-14) * (1 - c.RELATIVE_STATE_ERROR) ** 2,
                    (values + 1e-14) * (1 + c.RELATIVE_STATE_ERROR) ** 2,
                )
                other = previous[side]["equilibrium_area_gain_interval"]
                assert other[0] < g[0] < g[1] < other[1]
                gains.append(g)
        if tag == "held_fixed":
            factors = (
                (1 - c.RELATIVE_STATE_ERROR) ** 2,
                (1 + c.RELATIVE_STATE_ERROR) ** 2,
            )
            vals = [
                np.max((a - 1e-14) * factors[0] - (b + 1e-14) * factors[1])
                for a, b in ((forward, reverse.T), (reverse.T, forward))
            ]
            mismatch = float(max(vals))
            assert mismatch > previous["equilibrium_pair_asymmetry_lower"]
    return dict(
        probabilities_checked=200,
        minimum_margin_inside_count_interval=float(minimum),
        equilibrium_two_way_gain=gains,
        frozen_receiver_mismatch_lower=mismatch,
    )
