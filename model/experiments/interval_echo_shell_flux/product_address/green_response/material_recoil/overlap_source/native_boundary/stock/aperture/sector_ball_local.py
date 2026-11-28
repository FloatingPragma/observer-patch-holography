"""Uniform sign of a local entropy tangent over the complete reporting-rate interval."""

import math
import numpy as np

from . import sector_cadence as cadence
from . import sector_retention as experiment


def certify(words, delta, subset, mats, rates, witness):
    # Fix the selected survey and its ordinary prices before the continuum
    # test. The measured home is equivalent to its three diagonal outcomes.
    j = witness["home"]
    probabilities = np.array(witness["conditional_clock_probabilities"])
    prices = np.log(probabilities[1, 0] / probabilities[[0, 2], 1]) - 1e-12
    exchanged = np.array(witness["exchanged_clock_probabilities"])
    relative = prices - np.log(exchanged[1, 1] / exchanged[[0, 2], 0]) - 1e-12
    assert min(prices) > 0 and min(relative) > 0
    clock = experiment.ReturnClock((mats[2] @ mats[1] @ mats[0]).T.tocsr(), subset)
    tick = clock.grid * np.column_stack((words[:, j] == 0, words[:, j] == 2)).astype(
        np.int64
    )
    d = np.array(delta, float)
    size = sum(abs(d))
    err = np.zeros(2)
    series = []
    for _ in range(rates["last_age"]):
        tick, step = clock.advance(tick)
        err += step
        val = tick / clock.grid
        mean = d @ val
        width = size * err + 1e-11
        series.append(np.stack((mean - width, mean + width), axis=1))
    history = np.array(series)
    osc = np.ptp(val, axis=0) + 2 * err
    power, eta = cadence.mixing(clock, len(words))
    assert (
        power == rates["continuous_cadence"]["power"]
        and eta == rates["continuous_cadence"]["oscillation_contraction_upper"]
    )
    reply = []
    for i in range(2):
        tail = (size * osc[i] / 2 + 1e-11) * power / (1 - eta) + 1e-10
        lower, upper = history[:, i].T
        # Establish strict positivity of Delta(p)/p on the closed
        # complement variable q=1-p, including its sparse limit.
        pending = [(0.0, 1.0)]
        covered = []
        while pending:
            left, right = pending.pop()
            b = cadence.evaluate(left, right, lower, upper, tail)
            if b[0] > 0:
                covered.append([left, right, b[0]])
            elif right - left <= 2**-14:
                return dict(
                    proved=False, home=j, unresolved=[left, right], interval=list(b)
                )
            else:
                mid = (left + right) / 2
                pending.extend(((left, mid), (mid, right)))
        reply.append(
            dict(
                outcome=[0, 2][i],
                probability_gain_per_rate_lower=min(b[2] for b in covered),
                certified_intervals=covered,
                sparse_limit=cadence.evaluate(1, 1, lower, upper, tail),
            )
        )
    lower = (
        sum(
            c["probability_gain_per_rate_lower"] * v
            for c, v in zip(reply, prices, strict=True)
        )
        - 1e-12
    )
    return dict(
        proved=True,
        home=j,
        home_relative_modular_per_unit_rate_lower=sum(
            c["probability_gain_per_rate_lower"] * v
            for c, v in zip(reply, relative, strict=True)
        )
        - 1e-12,
        outcomes=reply,
        relative_ordinary_prices_lower=prices.tolist(),
        retained_entropy_per_unit_rate_lower=math.nextafter(lower, -math.inf),
    )
