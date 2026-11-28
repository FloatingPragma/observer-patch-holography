"""Read native acknowledgements and collect arrival histograms on a registered beam."""

import math
from fractions import Fraction
from itertools import pairwise, product

import numpy as np
from scipy.sparse.csgraph import shortest_path
from sympy import Matrix


def enumerate_routes(adj, d, s, t):
    if s == t:
        yield (t,)
        return
    for v in np.flatnonzero(adj[s] & (d[:, t] == d[s, t] - 1)):
        for rest in enumerate_routes(adj, d, v, t):
            yield (s, *rest)


def edge_times(adj, homes, loads, edges, probe, swap):
    residents = np.zeros(len(adj), int)
    residents[homes] = loads
    occupied = set(np.flatnonzero(residents))
    register = residents.copy()
    periods = {}
    for x, y in edges:
        block = sorted({x, y} & occupied)
        rails = next(
            rails
            for rails in product(
                *(
                    sorted(set(np.flatnonzero(adj[p])) - occupied - {x, y})
                    for p in block
                )
            )
            if len(set(rails)) == len(rails)
        )
        ticks = 0

        def pulse(i, j):
            nonlocal ticks
            register[i], register[j] = divmod(
                int(swap[3 * register[i] + register[j]]), 3
            )
            ticks += 1

        for p, r in zip(block, rails, strict=True):
            pulse(p, r)
        assert register[x] == register[y] == 0
        register[x] = probe
        pulse(x, y)
        assert register[y] == probe
        pulse(y, x)
        assert register[x] == probe
        register[x] = 0
        for p, r in zip(block, rails, strict=True):
            pulse(p, r)
        np.testing.assert_array_equal(register, residents)
        periods[x, y] = ticks // 2
    return periods


def measure(adj, ports, homes, loads, arms, d, witnesses, swap):
    reports = []
    rng = np.random.default_rng(8701491)
    n = 2**22
    # Per arm there are three timing probabilities and one port probability.
    # This union bound also covers the two witnesses, both loads and both beams.
    comparisons = len(witnesses) * 2 * 2 * 2 * 4
    error = 2 * math.sqrt(math.log(2 * comparisons / 0.005) / (2 * n))
    assert np.all(swap.sum(0) == 1)
    permutation = swap.argmax(0)
    for s, t, port in witnesses:
        routes = list(enumerate_routes(adj, d, s, t))
        edges = sorted({(x, y) for route in routes for x, y in pairwise(route)})
        final = np.array([p[-2] for p in routes])
        if port is None:
            # The port most responsive to the correlation in the prompt beam.
            k = np.array([homes[0] in path and homes[3] in path for path in routes])
            port = int(np.argmax(np.bincount(final, weights=k, minlength=len(ports))))
        readings = {}
        for probe in (1, 2):
            for policy in ("fixed", "fastest"):
                raw = []
                mean = []
                for arm in arms:
                    cells = np.zeros((2, 3), dtype=int)
                    expectation = np.zeros((2, 3))
                    num_by_word = rng.multinomial(n, [0.5, 0.5])
                    for word, nw in zip(arm, num_by_word, strict=True):
                        pulse = edge_times(
                            adj, homes, loads[word], edges, probe, permutation
                        )
                        cost = np.array(
                            [sum(pulse[u, v] for u, v in pairwise(p)) for p in routes]
                        )
                        shortest = cost.min()
                        assert shortest == d[s, t]
                        delay = (cost - shortest) // 2
                        assert np.all(shortest + 2 * delay == cost) and max(delay) <= 2
                        eligible = (
                            (cost == shortest)
                            if policy == "fastest"
                            else np.ones(len(routes), dtype=bool)
                        )
                        category = (
                            3 * (final[eligible] == port).astype(int) + delay[eligible]
                        )
                        prob = np.full(len(category), 1 / len(category))
                        observations = rng.multinomial(nw, prob)
                        np.add.at(cells.ravel(), category, observations)
                        np.add.at(expectation.ravel(), category, prob / 2)
                    assert cells.sum() == n
                    raw.append(cells)
                    mean.append(expectation)
                contrast = (raw[0] - raw[1]) / n
                forecast = mean[0] - mean[1]
                bounds = np.r_[contrast.sum(0), contrast[1].sum()]
                theory = np.r_[forecast.sum(0), forecast[1].sum()]
                assert np.max(abs(bounds - theory)) < error
                readings[f"{probe}/{policy}"] = {
                    "counts_by_arm_port_time": np.array(raw).tolist(),
                    "forecast": theory.tolist(),
                    "simultaneous_intervals": np.stack(
                        (bounds - error, bounds + error), axis=1
                    ).tolist(),
                }
        reports.append(
            {
                "ends": ports[[s, t]].tolist(),
                "port": int(ports[port]),
                "routes": len(routes),
                "trials_per_arm": n,
                "records": readings,
            }
        )
    return reports


def two_delay_channels(chart):
    """Read the full register from two exterior timing channels."""
    words, _, ports, _, adj, _, loads, _ = chart.apparatus()
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    distance = shortest_path(adj, unweighted=True).astype(int)
    s, t = np.searchsorted(ports, [72, 97])
    routes = list(enumerate_routes(adj, distance, s, t))
    selected = [
        next(p for p in routes if set(p) & set(homes) == {homes[k]}) for k in (0, 3)
    ]
    edges = sorted({edge for path in selected for edge in pairwise(path)})
    swap = chart.boundary.exchange_geometry.contact_swap().argmax(0)
    stocks, lookup = np.unique(loads, axis=0, return_inverse=True)
    tables = []
    for probe in (1, 2):
        rows = []
        for stock in stocks:
            ticks = edge_times(adj, homes, stock, edges, probe, swap)
            rows.append(
                [sum(ticks[edge] for edge in pairwise(path)) for path in selected]
            )
        tables.append(np.array(rows)[lookup])
    np.testing.assert_array_equal(tables[0], tables[1])
    decoded = (tables[0] - distance[s, t]) // 2
    np.testing.assert_array_equal(tables[0], distance[s, t] + 2 * decoded)
    np.testing.assert_array_equal(decoded, words[:, [0, 3]] > 0)
    return decoded, {
        "routes": [ports[list(p)].tolist() for p in selected],
        "blank_duration": int(distance[s, t]),
        "distinct_physical_loads": len(stocks),
        "verified_return_words": len(words),
    }


def invert(calibration, records, tide):
    """Infer both source arms from exterior arrival-port counts alone."""
    inv = np.asarray(
        Matrix(
            [[Fraction(z) for z in row] for row in calibration["ordinary_mixed_matrix"]]
        ).inv(),
        float,
    )
    base = np.array([float(Fraction(v)) for v in calibration["ordinary_intercepts"]])
    result = {}
    for probe in (1, 2):
        freq = []
        for coords, r in zip(calibration["rows"], records, strict=True):
            assert coords == [r["ends"][0], r["port"], r["ends"][1]]
            counts = np.array(
                r["records"][f"{probe}/fastest"]["counts_by_arm_port_time"]
            )
            freq.append(counts[:, 1, :].sum(1) / r["trials_per_arm"])
        means = (np.array(freq).T - base) @ inv.T
        # Twelve port frequencies (three baselines, two loads and two arms).
        width = (
            np.abs(inv).sum(1)
            * math.sqrt(math.log(2 * 12 / 0.005) / (2 * records[0]["trials_per_arm"]))
            + 1e-12
        )
        lo = np.maximum(0.0, means - width)
        hi = np.minimum(1.0, means + width)
        cov_lo = lo[:, 2] - hi[:, 0] * hi[:, 1]
        cov_hi = hi[:, 2] - lo[:, 0] * lo[:, 1]
        delta = np.array([cov_lo[0] - cov_hi[1], cov_hi[0] - cov_lo[1]])
        target = np.multiply.outer(delta, tide)
        assert 0 < delta[0] < 0.5 < delta[1]
        result[str(probe)] = {
            "source_means_xy_joint": means.tolist(),
            "simultaneous_halfwidths": width.tolist(),
            "covariance_change_interval": delta.tolist(),
            "predicted_connected_commute_acceleration_at_48_82": [
                float(target.min()),
                float(target.max()),
            ],
        }
    return result
