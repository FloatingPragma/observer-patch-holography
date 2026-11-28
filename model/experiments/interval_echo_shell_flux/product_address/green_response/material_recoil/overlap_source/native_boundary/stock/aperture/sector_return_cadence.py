"""Enclose the accepted-departure rates of local clocks across all native return cadences."""

import argparse
import itertools
import json
from hashlib import sha256

import numpy as np

from . import sector_register_pairs as pairs
from . import sector_register_source as source

HERE = source.HERE


def run():
    (ports, _, adj, _, cfg, ids, _), _, _, delta, eps = pairs.apparatus()
    degree = adj[np.searchsorted(ports, [8, 20, 32, 44])].sum(1)
    weights = degree / degree.sum()
    size, bank = pairs.get_bank()
    count, gamma = pairs.collective.sector_fronts_counts.moments(bank)
    cat = pairs.collective.peer_population.fibers.apparatus()[3]
    native, M = pairs.type_weights(cat, count, gamma, size)
    M /= size - 1
    _, _, sub, _, mats = pairs.chart.source.sources((8,))
    flags = cfg[ids]
    series, errors, tail, pos, neg, _, _ = source.conditional_history(
        sub, mats, flags, cat
    )
    base = flags.mean(0)
    ages = np.arange(len(series))
    blocks = json.loads((HERE / "sector_fronts_collective.json").read_text())[
        "continuous_cadence"
    ]["bins"]
    bins = []
    for block in blocks:
        left, right = block["rate"]
        points = np.linspace(left, right, 9)
        for a, b in itertools.pairwise(points):
            lo, hi = (1 - b) ** ages, (1 - a) ** ages
            e = hi @ errors + (1 - a) ** len(series) * tail + 1e-9
            low = (lo @ pos + hi @ neg).reshape(series.shape[1:]) - e
            high = (hi @ pos + lo @ neg).reshape(series.shape[1:]) + e
            x, dx = (low + high) / 2, (high - low) / 2
            d = (x * (M @ x)).sum(0)
            de = (dx * (abs(M) @ (2 * abs(x) + dx))).sum(0) + 1e-10
            # Relative to the blank departure frequency. Independently
            # paired members have P(two occupied)=mu^2 and E[sum]=2mu.
            mlo = np.maximum(
                0, base + np.minimum(a * (native @ low), b * (native @ low))
            )
            mhi = np.minimum(
                1, base + np.maximum(a * (native @ high), b * (native @ high))
            )
            assert np.all(mlo <= mhi)

            def independent(mu):
                return (1 - 2 * mu / 3 + mu**2 / 6) @ weights

            bins.append(
                {
                    "rate": [a, b],
                    "independent_frequency_over_blank": [
                        independent(mhi) - 1e-12,
                        independent(mlo) + 1e-12,
                    ],
                    "distinct_minus_independent_per_squared_rate": [
                        (d - de) @ weights / 6 - 1e-13,
                        (d + de) @ weights / 6 + 1e-13,
                    ],
                }
            )
    known = (np.diag(cfg.T @ delta["cross"] @ cfg) @ weights) / 6
    for row in bins:
        if row["rate"][0] <= 1 / 8 <= row["rate"][1]:
            low, high = (
                np.array(row["distinct_minus_independent_per_squared_rate"]) / 64
            )
            assert low <= known - eps["cross"] / 6 <= known + eps["cross"] / 6 <= high
    return {
        "input_sha256": {
            p: sha256((HERE / p).read_bytes()).hexdigest()
            for p in ("SECTOR_RETURN_LAPSE.md", "sector_fronts_collective.json")
        },
        "home_degrees": degree.tolist(),
        "response_at_eighth": float(known),
        "age_terms": len(series),
        "bins": bins,
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    z = run()
    dst = HERE / "sector_return_cadence.json"
    if cli.parse_args().write:
        dst.write_text(json.dumps(z, indent=2, allow_nan=False) + "\n")
    else:
        pairs.chart.source.balls.radar.lab.contact.timing.lab.lab.cold.apparatus.check(
            json.loads(dst.read_text()), z
        )
    print(z["home_degrees"], z["response_at_eighth"], len(z["bins"]))
