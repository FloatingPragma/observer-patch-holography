"""Inspection cadence on the retained 80-cell reciprocal-strain source.

The material and the deadline T=6 are fixed. A trial pays t/2 on a return
and T/2 otherwise. Completion at T counts as a return and costs the cap.
The shortest distinct-station completed range is h (two inspections).

Sharp inspections have a uniform obstruction. For target projector P, let
Q=1-P and J=||PHQ||. Duhamel about PHP+QHQ gives ||P exp(-ihH)Q|| <= hJ.
Every miss leaves the live state in Q. For relay tick k and return tick j>k,
the two terminal click amplitudes thus have product norm <= h² J_out J_in;
all intervening factors contract. Summing the n(n-1)/2 disjoint histories
bounds P(return) by n(n-1) h⁴ J_out² J_in²/2 = O(h²) at fixed T.
This also bounds the distance from the timeout range by (T/2)P(return),
uniformly in the retained state. More sharp inspections erase range contrast.

At rate gamma=1 the update on a miss is exp(-gamma h P/2)exp(-ihH).
The finite-dimensional product formula gives the limit exp[t(-iH-gamma P/2)].
A relay at s and return a duration t later have joint density gamma² times
the squared amplitude of P_a exp(G_a t)P_b exp(G_b s)P_a, with s+t<=T.
Near s=t=0 its amplitude divided by st converges to -C* C, C=P_b H P_a.
The direct 0--1 contact and the independently faithful source therefore give
strictly positive continuous weak returns. The code bounds a small rectangle
in this integral. Neither detection rule identifies a geodesic or proper time.
"""

import argparse
from fractions import Fraction
import json
from math import comb, expm1, factorial
from pathlib import Path

import numpy as np
from scipy.special import jv

from . import live_radar, live_source


HERE = Path(__file__).resolve().parent
STATIONS = (0, 1, 5, 40)
EDGES = ((0, 1), (2, 3))
DIVISORS = (2, 4, 8, 16, 32, 64)
DEADLINE = 6


def kernel_grid(integer_h, rows, times):
    """Unmonitored station blocks on a grid, via the Hermitian Chebyshev series."""
    scale = 24
    degree = 260
    h = integer_h / (40 * scale)
    assert abs(h).sum(axis=1).max() <= 1
    argument = scale * np.asarray(times)
    max_x = Fraction(scale * DEADLINE, 2)
    # |J_k(x)| <= (|x|/2)^k/k! for real x (Poisson integral).
    tail = 2*max_x**(degree+1)/factorial(degree+1)/(1-max_x/(degree+2))
    assert max(argument) <= 2*max_x and tail < Fraction(1, 10**28)
    sites = rows.ravel()
    old = np.zeros((h.shape[0], len(sites)))
    old[sites] = np.eye(len(sites))
    new = h @ old
    moments = np.empty((degree+1, len(sites), len(sites)))
    moments[0], moments[1] = old[sites], new[sites]
    for n in range(2, degree+1):
        old, new = new, 2*(h @ new)-old
        moments[n] = new[sites]
    n = np.arange(degree+1)
    coefficients = jv(n[None, :], argument[:, None]) * ((-1j)**n * np.where(n == 0, 1, 2))
    result = coefficients @ moments.reshape(degree+1, -1)
    sites_count, width = rows.shape
    return result.reshape(-1, sites_count, width, sites_count, width).transpose(0, 1, 3, 2, 4)


def first_hits(travel, resident, probability):
    """Retain all coin/strain amplitudes after the first registered click."""
    if not 0 < probability <= 1:
        raise ValueError("Expected a positive single-inspection detection probability")
    weight = probability / (1 + np.sqrt(1 - probability))
    pre_click = np.zeros_like(travel)
    for t in range(1, len(travel)):
        pre_click[t] = travel[t] - weight * (resident[t-1:0:-1] @ pre_click[1:t]).sum(axis=0)
    return np.sqrt(probability) * pre_click


def effects(kernels, left, right, probability):
    """Return completion-time effects for a fair orientation, including the deadline."""
    f = first_hits(kernels[:, right, left], kernels[:, right, right], probability)
    g = first_hits(kernels[:, left, right], kernels[:, left, left], probability)
    arrived = np.zeros_like(f)
    for t in range(2, len(f)):
        forward = g[t-1:0:-1] @ f[1:t]
        reverse = f[t-1:0:-1] @ g[1:t]
        arrived[t] = (live_radar.gram(forward) + live_radar.gram(reverse)).sum(axis=0)/2
    return arrived


def leakage_bound(integer_h, rows):
    # Integer Gershgorin upper bound on ||PHQ||², with H stored at scale 40.
    strip = integer_h[rows]
    inside = strip[:, rows].toarray()
    gram = (strip @ strip.T).toarray() - inside @ inside.T
    return Fraction(int(abs(gram).sum(axis=1).max()), 1600)


def weak_limit_bound(integer_h, rows):
    cross = integer_h[rows[1]][:, rows[0]].toarray()
    gram = cross.T @ cross
    assert np.array_equal(gram @ gram, (16*1600)*gram) and np.trace(gram) == 4*16*1600
    assert abs(integer_h).sum(axis=1).max()/40 + .5 <= 24
    # Since exp(tG) contracts, Taylor's integral remainder costs t²||G||²/2.
    # For 0<s,t<=a the normalized product differs from -C*C by <=8e+e².
    a = Fraction(1, 256)
    e = a * 24**2 / 2
    singular = 16 - 8*e - e**2
    assert singular > 0
    source_floor = Fraction(14, 10000)
    modular = json.loads((HERE/"live_modular.json").read_text())
    assert Fraction(str(modular["reference_minimum_eigenvalue_lower"])) >= source_floor
    return source_floor * singular**2 * a**6 / 9


def summaries(returns, step, states):
    t = np.arange(len(returns))*step
    data = {}
    for name, state in states.items():
        p = np.einsum("tij,ji->t", returns, state).real
        trace = 0 if name == "hidden" else 1
        data[name] = dict(complete=float(p.sum()), at_deadline=float(p[-1]),
                          saturated=float(trace-p[:-1].sum()),
                          capped_range=float(DEADLINE*trace/2-(DEADLINE-t)@p/2))
    return data


def run():
    source = live_radar.source
    board = source.formation.federation(1)
    route = source.formation.region(board)
    states = live_source.postpump(True)
    rows = live_radar.station_rows(board, STATIONS)
    finest = max(DIVISORS)
    times = np.arange(DEADLINE*finest+1)/finest
    records, weak_bounds, charge_leakage = {}, {}, {}
    for charge in (-1, 1):
        h, _, _ = source.system(board, route, charge, coupled=True)
        kernels = kernel_grid(h, rows, times)
        squared = [leakage_bound(h, rr) for rr in rows]
        charge_leakage[str(charge)] = list(map(str, squared))
        weak_bounds[str(charge)] = str(weak_limit_bound(h, rows))
        for divisor in DIVISORS:
            step = 1/divisor
            samples = kernels[::finest//divisor]
            for left, right in EDGES:
                for mode, probability in (("sharp", 1.), ("rate_one", -expm1(-step))):
                    key = (divisor, left, right, mode)
                    reading = summaries(effects(samples, left, right, probability), step, states)
                    if key not in records:
                        bound = comb(DEADLINE*divisor, 2)*squared[left]*squared[right]/divisor**4
                        records[key] = dict(step=step, ports=[STATIONS[left], STATIONS[right]],
                                            detector=mode, sharp_uniform_complete_upper=float(min(1, bound)),
                                            by_charge={})
                    records[key]["by_charge"][str(charge)] = reading
    for record in records.values():
        a, b = record["by_charge"].values()
        record["pooled"] = {name: {k:(v + b[name][k])/2 for k, v in values.items()}
                            for name, values in a.items()}
    return dict(deadline=DEADLINE, gamma=1, cells=board.cells, port_order=list(STATIONS),
                source="four-tick vacancy pump, retained strain and both heads",
                leakage_squared_upper=charge_leakage, continuous_near_return_lower=weak_bounds,
                panels=list(records.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    report = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (HERE/"live_cadence.json").write_text(report)
    else:
        print(report, end="")
