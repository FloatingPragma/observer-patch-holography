"""Direct native trajectories stopped by an observer's independent return clock."""

import math

import numpy as np


def report_from(tables, subset, native, rng, divisor, label):
    carrier = subset[label]
    final = np.full(len(label), -1)
    pending = np.arange(len(label))
    elapsed = 0
    for elapsed in range(1, 65537):
        for table in tables:
            carrier = table[carrier, rng.integers(128, size=len(carrier))]
        word = native[carrier]
        eligible = word >= 0
        observed = eligible & (rng.integers(divisor, size=len(carrier)) == 0)
        final[pending[observed]] = word[observed]
        carrier = carrier[~observed]
        pending = pending[~observed]
        if not len(pending):
            break
    assert not len(pending)
    return final, elapsed


def trajectories(tables, subset, native, rng, divisor, copies):
    label = np.repeat(np.arange(len(subset)), copies)
    final, elapsed = report_from(tables, subset, native, rng, divisor, label)
    return label, final, elapsed


def acquire(tables, subset, native, delta, calibration, seed=451637):
    rng = np.random.default_rng(seed)
    result = {}
    # Independent sharp start strata are recombined with the exact signed
    # source frequencies. Calibration has one interval for each ending word.
    for divisor, copies in ((8, 1024), (64, 2048)):
        label, final, elapsed = trajectories(
            tables, subset, native, rng, divisor, copies
        )
        stop = np.bincount(final, weights=delta[label] / copies, minlength=len(subset))
        middle = calibration.mean(1)
        read = float(stop @ middle)
        price_err = abs(stop) @ (np.ptp(calibration, axis=1) / 2)
        # Conditional on fixed strata, summands are independent regardless
        # of the number of unobserved returns on a trajectory.
        spread = calibration[:, 1].max() - calibration[:, 0].min()
        # Disjoint differences of two independent trips estimate the
        # variance within each start stratum. Allow the full thermal
        # calibration error on each pair before bounding the variance.
        pairs = final.reshape(len(subset), copies // 2, 2)
        half = np.ptp(calibration, axis=1) / 2
        diff = abs(middle[pairs[:, :, 0]] - middle[pairs[:, :, 1]]) + half[pairs].sum(2)
        pair_mean = (diff**2 / 2).mean(1)
        variance_upper = (delta**2) @ pair_mean + (spread**2 / 2) * math.sqrt(
            math.log(1 / 0.0005) * sum(delta**4) / copies
        )
        ell = math.log(2 / 0.002)
        err = math.sqrt(2 * variance_upper * ell / copies) + 2 * spread * max(
            abs(delta)
        ) * ell / (3 * copies)
        result[f"1/{divisor}"] = dict(
            trials=len(label),
            copies_per_start=copies,
            unreturned=0,
            maximum_sweeps=elapsed,
            heat_subtracted_read=read,
            sampling_width=err,
            weighted_variance_upper=variance_upper,
            calibration_width=price_err,
            bounds=[read - err - price_err, read + err + price_err],
            terminal_counts=np.bincount(final, minlength=len(subset)).tolist(),
        )
    return result
