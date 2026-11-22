"""Predict the registered two-load home probability from independent one-load samples."""
import json
from math import asin, expm1, log, pi, sin, sqrt

import numpy as np

from model.candidate_laws.interval_echo import candidate
from .coupling_budget import native_triangle


def infer_floor(quarter, irrational, full):
    """A conservative 99% simultaneous inference, or an inconclusive reading."""
    arms = np.asarray([*quarter, *irrational, full])
    if (arms.shape != (5, 3) or not np.issubdtype(arms.dtype, np.integer)
            or (arms < 0).any() or (arms > 2**24).any() or not np.all(arms.sum(axis=1) == 2**24)):
        raise ValueError("Expected five complete three-port count vectors at the fixed exposure count")
    shots = 2**24
    if full[2] != 0:
        return dict(resolved=False, reason="full pulse populated the spectator output")
    if full[1] != 0:
        return dict(resolved=False, reason="dark port occupied")
    # Exact one-sided binomial upper bound for zero hits.
    dark_upper = -expm1(log(1/300) / shots)
    delta = asin(sqrt(dark_upper)) / pi
    # Hoeffding for a sum of two independent frequencies, one bound for each duration.
    radius = sqrt(log(300)/shots)
    # A quiet spectator port after the last beam splitter provides an apparatus
    # check independent of strength. Grossly inconsistent histograms abstain.
    for j, t in enumerate((pi/4, pi/(4*sqrt(2)))):
        spectator = np.array([arms[2*j, 2], arms[2*j+1, 1]]) / shots
        target = np.array([sin(t/2)**2, 1-sin(t/2)**2])/2
        if np.max(abs(spectator-target)) > radius:
            return dict(resolved=False, reason="spectator counts disagree with the instrument")
    coarse_upper = 1.25 - (int(arms[0, 0]) + int(arms[1, 0])) / shots + radius
    zero_test_lower = 1.25 - (int(arms[2, 0]) + int(arms[3, 0])) / shots - radius
    if coarse_upper < 0 or zero_test_lower > 1:
        return dict(resolved=False, reason="no physical fringe deficit")
    if coarse_upper >= sin((1-delta)*pi/8)**2:
        return dict(resolved=False, reason="another integer class is admissible")
    if zero_test_lower <= sin(delta*pi/(8*sqrt(2)))**2:
        return dict(resolved=False, reason="the zero class is admissible")
    # The alias_window certificate covers every nonzero signed alias on [0,37].
    # Removing the scalar part of the perturbation bounds its norm by delta.
    leakage = sqrt(1-.94) + 37*delta
    return dict(resolved=True, distance_to_nonzero_alias_upper=delta,
                quarter_deficit_upper=coarse_upper, irrational_deficit_lower=zero_test_lower,
                home_probability_lower=1-min(1., leakage)**2)


def acquire(hamiltonian):
    readings = []
    for t, seeds in ((pi/4, (74504, 74505)), (pi/(4*sqrt(2)), (74506, 74507)), (2*pi, (74508,))):
        distributions = native_triangle(hamiltonian, t)
        for row, seed in enumerate(seeds):
            prob = distributions[row]
            np.testing.assert_allclose(prob.sum(), 1, atol=4e-12, rtol=0)
            counts = np.random.default_rng(seed).multinomial(2**24, prob/prob.sum())
            readings.append(counts.tolist())
    return dict(counts=readings, forecast=infer_floor(readings[:2], readings[2:4], readings[4]))


def run():
    # Only acquisition receives the physical Hamiltonian. The inference has no
    # parameter slot, pair dataset, model-labelled outcome, or trial-level truth.
    exchange, contrast = candidate.native_law.local_term_parts()
    result = {}
    for name, k, resolved in (("installed", candidate.KAPPA, True), ("other_alias", 16, True),
                              ("reversed", -8, True), ("removed", 0, False),
                              ("half_alias", 4, False), ("detuned", 8.25, False)):
        acquired = acquire(exchange + k*contrast)
        if acquired["forecast"]["resolved"] != resolved:
            raise ArithmeticError(f"Unexpected experimental decision for {name}")
        result[name] = acquired
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
