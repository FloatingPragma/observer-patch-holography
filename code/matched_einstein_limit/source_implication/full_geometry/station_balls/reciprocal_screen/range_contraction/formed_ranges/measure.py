"""Two-charge, gain-independent shape response to the same formed heads."""

import argparse
from itertools import combinations
import json
from math import ceil, floor
from pathlib import Path

from mpmath import iv
import numpy as np
import sympy as sp

from ..held_stations import native, QUARTETS, ROUTES
from ......full_algebra import certificate as cert

RECORD = Path(__file__).with_name("results.json")
SCALE = 10**12


def scaled_bounds(value):
    # One integer unit pads conversion from the 35-digit endpoint to binary64.
    assert max(abs(float(value.a)), abs(float(value.b))) * SCALE < 2**49
    return [floor(float((value * SCALE).a)) - 1, ceil(float((value * SCALE).b)) + 1]


def interval(bounds):
    return iv.mpf(bounds) / SCALE


def preparations():
    rho = native.formed_state()
    left = sp.Matrix(
        2, 2, lambda i, j: sum(rho[2 * i + k, 2 * j + k] for k in range(2))
    )
    right = sp.Matrix(
        2, 2, lambda i, j: sum(rho[2 * k + i, 2 * k + j] for k in range(2))
    )
    assert rho == rho.H and rho.trace() == 1
    return {"formed": rho, "independent": sp.kronecker_product(left, right)}


def native_ranges():
    board = native.iso.pairing.formation.federation(1)
    anchors = native.iso.pairing.formation.region(board)
    ports = (*anchors, 5, 6)
    ids = np.array([4 * p + j for p in ports for j in range(4)])
    co = native.prior.weights()
    states = preparations()
    # The existing formed-state enclosure is <1e-18 in trace norm; tensoring
    # the two separately traced marginals at most doubles this error.
    state_error = sp.Rational(1, 10**18)
    error = sp.Rational(native.prior.error_effect()) + 6 * state_error
    banks = {}
    for q, name in ((-1, "vacancy"), (1, "excess")):
        h = native.generator(board, anchors, q, False)
        kernels = native.kernel_panels(h, ids, co)
        nums = np.array(native.prior.direction_effects(kernels, 12), dtype=object)
        effects = cert.rational_effects(nums[:, None], 4 * native.prior.FRAC)
        banks[name] = {
            stock: [
                scaled_bounds(cert.ri_band(cert.inner(e[0], state), error))
                for e in effects
            ]
            for stock, state in states.items()
        }
    return ports, banks


def log_responses(banks):
    responses = {}
    for charge, stocks in banks.items():
        log_directed = {}
        for route, formed, product in zip(
            ROUTES, stocks["formed"], stocks["independent"], strict=True
        ):
            log_directed[route] = iv.ln(interval(formed) / interval(product))
        responses[charge] = {
            (i, j): (log_directed[i, j] + log_directed[j, i]) / 2
            for i, j in combinations(range(6), 2)
        }
    return responses


def compare(banks):
    responses = log_responses(banks)
    rows = []
    for i, j, k, last in QUARTETS:
        matches = [((i, j), (k, last)), ((i, k), (j, last)), ((i, last), (j, k))]
        for a, b in combinations(range(3), 2):
            per_charge = {}
            for q, v in responses.items():
                plus, minus = matches[a], matches[b]
                per_charge[q] = v[plus[0]] + v[plus[1]] - v[minus[0]] - v[minus[1]]
            gap = per_charge["excess"] - per_charge["vacancy"]
            rows.append(
                {
                    "quartet": [i, j, k, last],
                    "matching_pair": [a, b],
                    "log_cross_response": {
                        q: scaled_bounds(y) for q, y in per_charge.items()
                    },
                    "log_excess_over_vacancy": scaled_bounds(gap),
                    "excess_over_vacancy": scaled_bounds(iv.exp(gap)),
                }
            )
    return rows


def run():
    iv.dps = 35
    ports, banks = native_ranges()
    rows = compare(banks)
    return {
        "station_ports": list(ports),
        "directed_routes": list(map(list, ROUTES)),
        "integer_scale": SCALE,
        "mean_time_cap": 3,
        "stock_contrast": "formed / independent",
        "matching_order": ["ij,kl", "ik,jl", "il,jk"],
        "directed_mean_bands": banks,
        "comparisons": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    result = run()
    if parser.parse_args().write:
        RECORD.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert result == json.loads(RECORD.read_text())
        print("Replayed both species, both source stocks and all 45 quartet contrasts.")
