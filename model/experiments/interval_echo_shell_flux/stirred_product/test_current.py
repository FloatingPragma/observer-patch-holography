"""Independent checks on complete cycles, routed shell flux and clock units."""

from fractions import Fraction
import json
from math import log, sqrt

import numpy as np

from .measure import HERE, N, clocks, shell, survey, tree


def load():
    return json.loads((HERE / "counts.json").read_text())


def observations(bank, tags):
    layer, offered, _, _ = survey()
    depth, bias = bank["depth"], bank["biased"]
    directed = tags[:, :6].astype(np.int64).copy()
    if bias:
        directed[:, :2] += tags[:, 6:]
    for n in range(1, depth + 1):
        moves = []
        for p in np.flatnonzero(layer == n):
            for k, q in enumerate(offered[p]):
                if layer[q] == n + 1:
                    moves.append((k, int(directed[p, k]), int(directed[q, k ^ 1])))
        yield np.array(moves)


def test_cuts_and_every_edge_resolution():
    r = load()
    assert r["cut_survey"] == [8, 24, 96, 216, 384, 600]
    with np.load(HERE / "edges.npz") as a:
        for key, b in r["banks"].items():
            samples = 0
            refusals = 0
            max_log = log(2 * 32768 / 0.005)
            for n, records in enumerate(observations(b, a[key + "_tags"]), 1):
                axis, out, inc = records.T
                c = len(out)
                assert c == 24 * n * n
                net = out - inc
                assert net.sum() == b["cycles"]
                # Condition on the numbers of moves up and down in every
                # layer. Each destination edge is drawn by independent uniform
                # choice among that category's actual (target,tag) words.
                weights = np.where((axis < 2) & b["biased"], 2.0, 1.0)
                prob = weights / sum(weights)
                variance = (sum(out) + sum(inc)) * prob * (1 - prob)
                bound = np.sqrt(2 * variance * max_log) + 2 * max_log / 3
                assert np.all(abs(net - b["cycles"] * prob) <= bound)
                assert np.all(net > 0)
                if b["biased"] and n == b["depth"]:
                    uniform = 1 / c
                    sigma = (sum(out) + sum(inc)) * uniform * (1 - uniform)
                    error = sqrt(2 * sigma * max_log) + 2 * max_log / 3
                    assert np.any(abs(net - b["cycles"] * uniform) > error)
                sites, _ = shell(n, b["biased"])
                samples += a[key + "_tags"][sites].sum()
                refusals += b["target_retries"][str(n)]
            assert samples + 2 * b["cycles"] == b["appointments"]
            assert refusals > 0


def test_vertex_by_vertex_with_tangential_and_aperture_moves():
    layer, offers, corners, _ = survey()
    r = load()
    with np.load(HERE / "edges.npz") as a:
        for k, b in r["banks"].items():
            tags = a[k + "_tags"]
            to_root = a[k + "_to_root"]
            from_root = a[k + "_from_root"]
            balance = np.zeros(N, dtype=np.int64)
            balance[corners] -= b["source_corners"]
            assert sum(b["source_corners"]) == sum(b["collectors"]) == b["cycles"]
            for n in range(1, b["depth"] + 1):
                parents = tree(n)
                idx, dst = shell(n, b["biased"])
                for p, par in parents.items():
                    if p == par:
                        assert to_root[p] == from_root[p] == 0
                    else:
                        assert par in offers[p] and layer[p] == layer[par]
                        f = int(to_root[p]) - int(from_root[p])
                        balance[p] += f
                        balance[par] -= f
                balance[idx] += tags[idx].sum(axis=1)
                np.add.at(balance, dst.ravel(), -tags[idx].ravel().astype(np.int64))
            absorb = layer == b["depth"] + 1
            np.testing.assert_array_equal(balance[absorb], -np.array(b["collectors"]))
            assert np.all(balance[~absorb] == 0)
            assert b["tangential_pulses"] == int(to_root.sum() + from_root.sum())
            assert np.count_nonzero(tags[absorb]) == 0
            assert b["direction_idles"] == (0 if b["biased"] else tags[:, 6:].sum())


def test_independent_source_and_renewal_clock():
    r = load()
    for key, b in r["banks"].items():
        D = b["depth"]
        mean, var = map(Fraction, clocks(D, b["biased"]))
        assert r["renewal_moments"][f"{D}_{int(b['biased'])}"] == [str(mean), str(var)]
        if not b["biased"]:
            assert mean == 2 + Fraction(4 * D * (D + 1), 3)
        observed = b["appointments"] / b["cycles"]
        # Eight independent banks. Cumulative clock error at most .005;
        # per-edge bounds use the remaining .005 of the simultaneous budget.
        delta = sqrt(float(var) / (b["cycles"] * (0.005 / 8)))
        assert abs(observed - float(mean)) <= delta
        other = r["banks"][key.rsplit("_", 1)[0] + "_source"]
        assert abs(observed - other["appointments"] / other["cycles"]) <= 2 * delta
    direct = r["literal"]
    assert direct["net_by_shell"] == [direct["episodes"]] * 5
    assert sum(direct["tag_counts"]) + 2 * direct["episodes"] == direct["appointments"]
    assert direct["stir_pulses"] > direct["move_pulses"] > 0 and direct["retries"] > 0
    c = r["parked"]
    assert sum(c["proposed_direction_tags"]) == c["attempted_ticks"]
    assert c["flag_fired"] == c["successful_coordinate_changes"] == c["drains"] == 0
