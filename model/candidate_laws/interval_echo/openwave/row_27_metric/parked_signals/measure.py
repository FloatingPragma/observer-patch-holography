"""Exact finite-exposure differences of exterior signal records on parked source pairs."""

import argparse
from fractions import Fraction
from itertools import product
import json
from math import factorial, inf, nextafter
from pathlib import Path

import numpy as np

from openwave.interval_echo.row_19.selected_collar import midpoint_clocks as clocks
from openwave.interval_echo.row_19.selected_collar import parked_energy_charge as source
from .ordinary import compare as ordinary_signals

HERE = Path(__file__).resolve().parent
TIMES = (1, 4, 16)
ORDER = 160
CHANNELS = ("near_return", "near_to_far", "far_return")


def bounds(center, error):
    return [nextafter(float(center - error), -inf), nextafter(float(center + error), inf)]


def weights(time):
    # Every Taylor monomial of exp[-it(L-5I)/2] has this common denominator.
    den = 2**ORDER * factorial(ORDER)
    return [(-1)**(n//2) * time**n * (den // (2**n*factorial(n))) for n in range(ORDER+1)]


def numerators(adjacency, degree, ports, coefficients):
    n = len(degree)
    neighbors = np.full((n, 5), n)
    for r in range(n):
        col = adjacency.indices[adjacency.indptr[r]:adjacency.indptr[r+1]]
        neighbors[r, :len(col)] = col
    # Add a zero row to pad missing local neighbors.
    power = np.zeros((n+1, 2), dtype=object)
    power[ports, np.arange(2)] = 1
    diag = np.asarray(degree-5, dtype=object)[:, None]
    moments = []
    for step in range(ORDER + 1):
        a, b = power[ports]
        if a[1] != b[0]:
            raise AssertionError("The free signal ceased to be reciprocal")
        moments.append((a[0], b[0], b[1]))
        if step != ORDER:
            power[:n] = diag*power[:n] - power[neighbors].sum(axis=1)
    result = []
    for coeff in coefficients:
        entries = []
        for seq in zip(*moments):
            terms = [w*m for w, m in zip(coeff, seq, strict=True)]
            entries.append(sum(terms[::2])**2 + sum(terms[1::2])**2)
        result.append(entries)
    return np.array(result, dtype=object)


def run():
    site = json.loads((clocks.HERE/"settled_midpoint_force.json").read_text())
    pairs = site["pairs"]
    ports = site["potential_edges"][0]["ports"]
    owners = site["potential_edges"][0]["owners"]
    assert ports == [34, 116] and owners == [2, 9]
    board = clocks.force.federation(0)
    routes = clocks.force.routes
    # This also checks the action of all local 99-state native columns on both loads.
    gates, changed = routes.calibrate(routes.codebook()[0])
    active = routes.graph(board, gates)
    baseline = active(-1, -1, 0, 0)
    lap = clocks.laplacian(baseline)
    hands = np.column_stack([clocks.clock(lap, owner)[0] for owner in owners])
    ordinary = clocks.clock(lap, 0)[0][:, 0]
    for load, interaction in clocks.force.meter.contact_terms().items():
        h, x, y = clocks.force.meter.two_loads(lap * 0.5, baseline[0], interaction, load)
        waves = ordinary[x, None]*hands[y, :]
        if load == 1:
            waves += ordinary[y, None]*hands[x, :]
        np.testing.assert_array_equal((2*h).astype(int) @ waves, waves*[8, 10, 8, 10])

    ref, fresh, shared = source.reference.local.survey.standards(board)[0][:3]
    _, _, stocks = source.controls(fresh)
    w = np.array([1, 0, -1, -1, 0, 1, 0, 0, 0])
    stocks["fourth_order"] = stocks["fresh"] + 32*np.outer(w, w)
    stocks["shared"] = 8*shared
    source_record = json.loads((clocks.HERE/"parked_energy_charge.json").read_text())
    # The first integer null vector of the stacked class, energy-bin and bit maps.
    z = np.array([0, -1, -1, 1, 0, 0, 1, 0, 0])
    null = 16*(np.outer(w, z) + np.outer(z, w))
    stocks["same_charge"] = stocks["plus"] + null
    for price in np.unique(ref):
        assert null[ref == price].sum() == 0
    seam = source.station_records()[1]
    for features in (seam, source.categorical(np.array(source_record["station_energy"])),
                     source.categorical(np.array(source_record["arrival_bit"]))):
        assert not np.any(features.T @ null @ features)
    reference = stocks["fresh"]
    for stock in stocks.values():
        assert stock.min() > 0
        for axis in (0, 1):
            np.testing.assert_array_equal(stock.sum(axis=axis), reference.sum(axis=axis))

    coefficients = [weights(t) for t in TIMES]
    totals = {key: np.zeros((len(TIMES), 3), dtype=object) for key in stocks}
    for left, right in pairs:
        for a, b in product(range(9), repeat=2):
            adjacency, deg = active(left, right, a+1, b+1)
            if max(deg) > 5:
                raise ValueError("The registered spectral error bound would fail")
            observed_lap = clocks.laplacian((adjacency, deg))
            np.testing.assert_array_equal(observed_lap @ hands, hands*[4, 6, 4, 6])
            raw = numerators(adjacency, deg, ports, coefficients)
            for tag, table in stocks.items():
                totals[tag] += int(table[a, b])*raw
    polyd = 2**ORDER*factorial(ORDER)
    normal = polyd**2*len(pairs)*int(reference.sum())
    readings = {}
    for tag, table in stocks.items():
        l1 = Fraction(int(abs(table-reference).sum()), int(reference.sum()))
        rows = {}
        for i, t in enumerate(TIMES):
            amp_err = Fraction(5*t, 2)**(ORDER+1)/factorial(ORDER+1)
            prob_err = amp_err*(2+amp_err)
            detector = {}
            for j, channel in enumerate(CHANNELS):
                value = Fraction(totals[tag][i, j], normal)
                delta = Fraction(totals[tag][i, j] - totals["fresh"][i, j], normal)
                detector[channel] = {"probability": bounds(value, prob_err),
                                     "change_from_fresh": bounds(delta, prob_err*l1)}
            rows[str(t)] = detector
        charge_tag = "plus" if tag == "same_charge" else tag
        charge = source_record["charge_change_nats"].get(charge_tag, {"exact": [0., 0.]})
        readings[tag] = {"times": rows, "charge_change_nats": charge["exact"]}
    return {"times": list(TIMES), "order": ORDER, "ports": ports, "owners": owners,
            "source_pairs": pairs, "source_words_checked": 81*len(pairs),
            "exterior_clock_gap_all_words_and_ordinary": 1,
            "one_and_two_load_identical": gates["1"] == gates["2"],
            "changed_local_rule_arrivals": changed, "stocks": readings,
            "ordinary_and_empty_signals": ordinary_signals(baseline[0], lap, ports, TIMES)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    path = HERE/"results.json"
    output = run()
    stocks = output["stocks"]
    plus = stocks["plus"]["times"]["16"]
    matched = stocks["same_charge"]["times"]["16"]
    assert stocks["same_charge"]["charge_change_nats"] == stocks["plus"]["charge_change_nats"]
    assert all(plus[k]["change_from_fresh"][1] < 0 for k in CHANNELS)
    for name in CHANNELS:
        left = plus[name]["probability"]
        right = matched[name]["probability"]
        assert left[1] < right[0] or right[1] < left[0], (name, left, right)
    if parser.parse_args().write:
        path.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    else:
        saved = json.loads(path.read_text())
        numerical = saved.pop("ordinary_and_empty_signals")
        rebuilt = output.pop("ordinary_and_empty_signals")
        for arm in ("empty", "ordinary"):
            for t in TIMES:
                np.testing.assert_allclose(numerical[arm][str(t)], rebuilt[arm][str(t)],
                                           atol=2e-11, rtol=0)
        old_gates = saved.pop("changed_local_rule_arrivals")
        new_gates = output.pop("changed_local_rule_arrivals")
        for load in ("1", "2"):
            np.testing.assert_allclose(new_gates[load], old_gates[load], atol=2e-13, rtol=0)
        if saved != output:
            raise AssertionError("The exact parked-signal replay differs from its archive")
    print("Checked all 486 sharp source configurations, exterior hands and reciprocal signals.")
