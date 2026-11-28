"""Can one elapsed-time change reproduce three native exterior detection curves?"""

import argparse
from fractions import Fraction as F
import json
from math import comb, factorial
from pathlib import Path

import numpy as np
from scipy import sparse

from .measure import CHANNELS, HERE, clocks

# 5**24 < 2**56. Every matrix-power entry and its one-step products fit int64.
MAX_POWER = 24
ARMS = ("fresh", "plus", "minus", "fourth_order", "shared")


def probability_polynomials(adjacency, degree, ports):
    """Numerators of even powers in three individual detection probabilities."""
    assert max(degree) <= 5 and 2*5**MAX_POWER < 2**63
    centered = sparse.diags(degree-5, dtype=np.int64, format="csr") - adjacency.astype(np.int64)
    block = np.zeros((len(degree), 2), dtype=np.int64)
    block[ports, np.arange(2)] = 1
    powers = []
    for n in range(MAX_POWER+1):
        a, b = block[ports]
        assert a[1] == b[0]
        powers.append([int(a[0]), int(b[0]), int(b[1])])
        if n < MAX_POWER:
            block = centered @ block
    rows = []
    for n in range(0, MAX_POWER+1, 2):
        rows.append([(-1)**(n//2)*sum(
            (-1)**j*comb(n, j)*seq[j]*seq[n-j] for j in range(n+1)
        ) for seq in zip(*powers, strict=True)])
    return np.array(rows, dtype=object)


def time_identity(baseline, observed):
    """Leading coefficients forced by reparametrizing near-return elapsed time."""
    slope = baseline[1, :]
    delta = observed - baseline
    assert np.all(delta[0] == 0) and np.all(slope != 0)
    for k in range(1, len(baseline)):
        if any(delta[k]):
            shifted = [d/b for d, b in zip(delta[k], slope, strict=True)]
            residual = [d - b*shifted[0] for d, b in zip(delta[k], slope, strict=True)]
            return {"first_changed_time_power": 2*k, "blank_slopes_in_squared_time": list(map(str, slope)),
                    "first_changes": list(map(str, delta[k])),
                    "implied_squared_time_shift_coefficients": list(map(str, shifted)),
                    "residual_after_matching_near_return": list(map(str, residual))}
    raise ValueError("No change through the investigated order")


def run():
    source = json.loads((clocks.HERE/"parked_energy_charge.json").read_text())
    # Sources were observed in lexicographic order of their two sharp energy-meter
    # reports. Map the archived counts back to the calibrated nine geometry words.
    rank = np.argsort(source["meter_order"])
    stocks = {arm: np.array(source["counts"][arm]).reshape(9, 9)[np.ix_(rank, rank)] for arm in ARMS}
    np.testing.assert_array_equal(stocks["fresh"], source["fresh_counts"])
    for m in stocks.values():
        for axis in (0, 1):
            np.testing.assert_array_equal(m.sum(axis=axis), stocks["fresh"].sum(axis=axis))
    site = json.loads((clocks.HERE/"settled_midpoint_force.json").read_text())
    ports = site["potential_edges"][0]["ports"]
    pairs = site["pairs"]
    routes = clocks.force.routes
    gates, _ = routes.calibrate(routes.codebook()[0])
    active = routes.graph(clocks.force.federation(0), gates)
    totals = {arm: np.zeros((MAX_POWER//2+1, len(CHANNELS)), dtype=object) for arm in ARMS}
    for left, right in pairs:
        for a, b in np.ndindex(9, 9):
            raw = probability_polynomials(*active(left, right, a+1, b+1), ports)
            for arm in ARMS:
                totals[arm] += int(stocks[arm][a, b])*raw
    trials = len(pairs)*source["denominator"]
    for k in range(MAX_POWER//2+1):
        denom = trials * 2**(2*k) * factorial(2*k)
        for arm in ARMS:
            totals[arm][k] = [F(n, denom) for n in totals[arm][k]]
    tests = {arm: time_identity(totals["fresh"], totals[arm]) for arm in ARMS[1:]}
    for record in tests.values():
        residual = list(map(F, record["residual_after_matching_near_return"]))
        assert residual[0] == 0 and any(residual[1:])
    return {"ports": ports, "order": list(CHANNELS), "sharp_configurations": 81*len(pairs),
            "word_order_from": "parked_energy_charge.json:meter_order", "tests": tests}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE/"common_delay.json"
    actual = run()
    if parser.parse_args().write:
        path.write_text(json.dumps(actual, indent=2)+"\n")
    elif actual != json.loads(path.read_text()):
        raise AssertionError("Native common-delay certificate differs from archive")
    print("Verified exact leading obstruction to common retiming in all four parked source arms.")
