"""Count all 30 return clocks and compare every four-station volume on three stocks."""

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path

from mpmath import iv
import numpy as np

from ... import geometry
from . import measure as survey

history = survey.clock_history
native = history.native
ROUTES = tuple(permutations(range(6), 2))
QUARTETS = tuple(combinations(range(6), 4))
STARTS = 2**50
SEED = 261127
RECORD = Path(__file__).with_suffix(".json")


def return_bins():
    board = native.iso.pairing.formation.federation(1)
    region = native.iso.pairing.formation.region(board)
    ports = (*region, 5, 6)
    columns = np.array([4 * p + v for p in ports for v in range(4)])
    weights = native.prior.weights()
    bins = np.zeros((30, 10, 2, 4, 4), dtype=complex)
    for parity in (0, 1):
        for charge in (-1, 1):
            h = native.generator(board, region, charge, bool(parity))
            kernels = native.kernel_panels(h, columns, weights)
            hits = {route: history.first_hits(kernels, *route) for route in ROUTES}
            for k, (i, j) in enumerate(ROUTES):
                for tick in range(2, 12):
                    effect = np.zeros((2, 4, 4), dtype=object)
                    for a in range(1, tick):
                        amp = native.prior.multiply(hits[j, i][tick - a], hits[i, j][a])
                        effect += native.prior.multiply(native.prior.adjoint(amp), amp)
                    bins[k, tick - 2, parity] += (
                        effect[0].astype(float) + 1j * effect[1].astype(float)
                    ) / (2 * history.F)
    return ports, bins


def acquire(state, bins, rng):
    success = np.einsum("rtpij,pji->rt", bins, state).real
    # Some exact-zero bins acquire a tiny negative value during evaluation.
    assert success.min() > -1e-13 and success.sum(axis=1).max() < 1 - 1e-12
    success = np.maximum(success, 0)
    probs = np.column_stack([success, 1 - success.sum(axis=1)])
    return [rng.multinomial(STARTS, p).tolist() for p in probs]


def volumes(counts):
    # Bonferroni-Hoeffding across 90 capped means with scores in [1/2, 3].
    # The second summand bounds bin-effect, dyadic stock and binary64 errors.
    assert float(history.timing_error()[1]) < 1e-12
    halfwidth = iv.mpf("2.5") * iv.sqrt(iv.ln(18000) / (2 * STARTS)) + iv.mpf("1e-10")
    means = np.full((6, 6), iv.mpf(0), dtype=object)
    for route, bins in zip(ROUTES, counts, strict=True):
        assert len(bins) == 11 and min(bins) >= 0 and sum(bins) == STARTS
        score = sum(tick * n for tick, n in enumerate(bins, 2))
        mean = iv.mpf(score) / (4 * STARTS)
        means[route] = mean + iv.mpf([-1, 1]) * halfwidth
    squared = geometry.observed_squared(means)

    def gram(i, rest):
        return np.array(
            [[geometry.dot(squared, i, j, k) for k in rest] for j in rest],
            dtype=object,
        )

    # Positive leading principal minors imply that the displacement Gram
    # has rank five. Its principal subconfigurations are genuine tetrahedra.
    full = gram(0, range(1, 6))
    minors = [geometry.curved.det_interval(full[:k, :k]) for k in range(1, 6)]
    assert all(det.a > 0 for det in minors)
    values = []
    for i, *others in QUARTETS:
        det = geometry.curved.det_interval(gram(i, others))
        assert det.a > 0
        values.append(iv.sqrt(det) / 6)
    return values, [geometry.interval_numbers(det) for det in minors]


def decode(banks):
    readings = {}
    intervals = {}
    for name, counts in banks.items():
        values, minors = volumes(counts)
        intervals[name] = values
        readings[name] = {
            "volume_bands": [geometry.interval_numbers(v) for v in values],
            "five_displacement_gram_leading_minors": minors,
        }
    changes = {}
    for name in ("ordinary_write", "independent_right"):
        rows = []
        for before, after in zip(intervals["reference"], intervals[name], strict=True):
            lo, hi = geometry.interval_numbers(after - before)
            sign = "shrinks" if hi < 0 else "expands" if lo > 0 else "undecided"
            rows.append({"difference": [lo, hi], "sign": sign})
        changes[name] = rows
    return {"preparations": readings, "volume_changes": changes}


def run():
    iv.dps = 35
    rng = np.random.default_rng(SEED)
    ports, bins = return_bins()
    banks = {k: acquire(state, bins, rng) for k, state in survey.preparations().items()}
    return {
        "station_ports": list(ports),
        "quartets": list(map(list, QUARTETS)),
        "directed_routes": list(map(list, ROUTES)),
        "return_or_timeout_ticks": list(range(2, 13)),
        "starts_per_route": STARTS,
        "count_seed": SEED,
        "simultaneous_coverage": 0.99,
        "counts": banks,
        **decode(banks),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        RECORD.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    else:
        assert result == json.loads(RECORD.read_text())
        print("Replayed all 90 count banks and the 15 matched tetrahedra.")
