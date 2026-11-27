"""Read energies, static work and local charge gains from terminal count banks."""

from functools import cache
from math import asin, log, pi, sqrt

import numpy as np
from model.experiments.count_intervals import clopper_pearson

SHOTS = 2**24
GRID = 2**40
TIMES = 2.0 ** np.arange(-8, 8)
TAIL = 0.01 / 2_000_000
BIN_ERROR = sqrt(log(2 / TAIL) / (2 * SHOTS))
LINE_WIDTH = 2**-30
ANGLE_ERROR = asin(2 * sqrt(2) * (BIN_ERROR + 1 / GRID)) + TIMES[-1] * LINE_WIDTH
ENERGY_ERROR = ANGLE_ERROR / TIMES[-1] + 2**-32


def clock_energy(counts, shots=SHOTS):
    """Lift successive measured angles using only the preceding clock reading."""
    real, imag = (2 * np.asarray(counts) / shots - 1).T
    phases = np.arctan2(imag, real)
    angle = phases[0]
    for current in phases[1:]:
        angle = current + 2 * pi * round((2 * angle - current) / (2 * pi))
    return float(angle / TIMES[-1])


@cache
def amplitude_interval(count, dimension):
    lo, hi = clopper_pearson(count, SHOTS, TAIL)
    rounding = dimension / GRID
    return sqrt(max(0, lo - rounding)), sqrt(min(1, hi + rounding))


def radius(counts, distances):
    d = np.asarray(distances)
    r2 = float(d**2 @ counts / SHOTS)
    # E[d²] is the sum of (2j-1) P(d >= j), including every unsampled far tail.
    low = high = 0.0
    rounding = len(d) / GRID
    for j in range(1, int(d.max()) + 1):
        lo, hi = clopper_pearson(int(np.asarray(counts)[d >= j].sum()), SHOTS, TAIL)
        low += (2 * j - 1) * max(0, lo - rounding)
        high += (2 * j - 1) * min(1, hi + rounding)
    return dict(rms=sqrt(r2), interval=[sqrt(low), sqrt(high)])


def contact_interval(neutral_h, aperture, bank):
    """Sum local Perron equations, using only the occupied probe neighborhood."""
    if bank["multiplicity"] != 1:
        return None
    c = np.asarray(bank["address"])
    energy = clock_energy(bank["clock"])
    # Coefficients on sqrt(p_i) in sum_{i in A} [(E-H0) sqrt(p)]_i.
    weights = (
        energy * aperture.astype(float)
        - np.asarray(neutral_h[aperture].sum(axis=0)).ravel()
    )
    ix = np.flatnonzero((weights != 0) | aperture)
    lo, hi = np.array([amplitude_interval(int(c[i]), len(c)) for i in ix]).T
    lower_den, upper_den = lo[aperture[ix]].sum(), hi[aperture[ix]].sum()
    if lower_den == 0:
        return None
    weights = weights[ix]
    lower_num = np.minimum(weights * lo, weights * hi).sum()
    upper_num = np.maximum(weights * lo, weights * hi).sum()
    ratios = [v / w for v in (lower_num, upper_num) for w in (lower_den, upper_den)]
    return [float(min(ratios) - ENERGY_ERROR), float(max(ratios) + ENERGY_ERROR)]


def signed(values):
    return sum(s * t * values[f"{s},{t}"] for s in (-1, 1) for t in (-1, 1)) / 4


def report(full, partition, masks, gain, direct, neutral_h):
    aperture = masks["neighbors"]
    energies = {key: clock_energy(bank["clock"]) for key, bank in full.items()}
    cut = {
        key: min(clock_energy(v["clock"]) for v in sides)
        for key, sides in partition.items()
    }
    two_baths = {
        key: sum(clock_energy(v["clock"]) for v in sides)
        for key, sides in partition.items()
    }
    sudden, reverse, vertex = {}, {}, {}
    for s in (-1, 1):
        sham = full[f"{s},0"]
        prob = aperture @ np.asarray(sham["address"]) / SHOTS
        for t in (-1, 1):
            key = f"{s},{t}"
            prob_after = aperture @ np.asarray(full[key]["address"]) / SHOTS
            delta = energies[key] - energies[f"{s},0"]
            coefficient = gain[str(t)]
            before, after = coefficient * prob, coefficient * prob_after
            tol = (
                2 * ENERGY_ERROR
                + abs(coefficient) * (BIN_ERROR + len(aperture) / GRID)
                + LINE_WIDTH
            )
            if not after - tol <= delta <= before + tol:
                raise ArithmeticError(
                    "Static work violates the two independent profile bounds"
                )
            sudden[key], reverse[key] = before, after
            vertex[key] = contact_interval(neutral_h[s], aperture, full[key])
    reconstruction = signed(sudden) + direct
    full_work = signed(energies) + direct
    profile_error = max(abs(v) for v in gain.values()) * (
        BIN_ERROR + len(aperture) / GRID
    )
    discrepancy = abs(full_work - reconstruction)
    return dict(
        relaxed=full_work,
        frozen_profile=reconstruction,
        reverse_profile=signed(reverse) + direct,
        relaxation=full_work - reconstruction,
        side_guard=signed(cut) + direct,
        distinct_baths=signed(two_baths) + direct,
        transmission_difference=full_work - (signed(cut) + direct),
        energy_error=ENERGY_ERROR,
        profile_error=profile_error,
        exceeds_ten_percent=bool(
            discrepancy - ENERGY_ERROR - profile_error
            > 0.1 * (abs(reconstruction) + profile_error)
        ),
        local_contact_intervals=vertex,
    )
