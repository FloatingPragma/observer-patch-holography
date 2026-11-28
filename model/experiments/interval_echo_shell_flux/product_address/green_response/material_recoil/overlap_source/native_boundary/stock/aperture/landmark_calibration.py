"""Calibrate the energy-matched overlap repair against uniformly added ordinary singles."""

import numpy as np
from . import collective_source as source


def measured_equivalence(z):
    words, probs, stats = source.records()
    n = (words[:, :4] == 1).sum(1)
    blank = words[:, -1] == 0
    partial = np.isin(words[:, -1], source.PARTIAL)
    matched = tuple(z @ f for f in (blank, partial, n, partial * n))
    moments = stats | {"matched": matched}
    ensemble = words, probs | {"matched": z}, moments
    a, b, N, C = stats["original"]
    A, B, M, D = matched
    assert b == B and D == C
    coefficients = (4 * (A - a) / b, M - N, N - M)
    banks = {}
    for volume in (240, 960, 3840):
        for m in (1, 4):
            nr = coefficients[0] * volume + (m - 1) * coefficients[1]
            ref = ["original"] * m
            act = ["matched"] * m
            actual = source.chances(volume, act, moments=moments)
            ordinary = source.chances(volume, ref, nr)
            assert actual == ordinary and nr > 0
            # Two independent native gate banks. The second actually mixes
            # adjacent integer donor counts on fresh attempts.
            readings = []
            for j, (arms, donors) in enumerate(((act, 0), (ref, nr))):
                readings.append(
                    source.pooled_draw(
                        volume,
                        arms,
                        0,
                        donors,
                        71190 + volume + 5 * m + j,
                        ensemble=ensemble,
                    )
                )
            banks[f"{volume}/{m}"] = dict(
                ordinary_singles=str(nr),
                native_probability=str(actual[0]),
                records=readings,
            )
    return dict(
        ordinary_coefficients=list(map(str, coefficients)),
        ordinary_moments=list(map(str, stats["original"])),
        matched_moments=list(map(str, matched)),
        native_gate_acquisitions=banks,
    )
