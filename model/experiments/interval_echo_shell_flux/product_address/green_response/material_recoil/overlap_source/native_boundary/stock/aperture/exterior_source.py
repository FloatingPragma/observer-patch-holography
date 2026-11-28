"""Local premise and held observer ranges for a resolved exterior tidal chord."""

from fractions import Fraction as Q
import json

import numpy as np
from scipy.stats import beta

from . import tidal_release as tidal

lab = tidal.lab


def source_bank():
    _, words, cells, _, codes, h, ids, _ = lab.energy.pair.setup()
    ref = np.full(len(words), Q(1, len(words)), object)
    repaired = lab.energy.least_change(cells) / len(words)
    witness = json.loads((tidal.HERE / "moving_records.json").read_bytes())[
        "composition"
    ]["same_charge_and_energy_witnesses"]
    W = {"original": ref, "source": repaired}
    controls = np.vstack(
        (
            [cells == i for i in range(4)],
            words[:, :4].T,
            h.diagonal()[codes],
            lab.energy.returned(words)[-1],
        )
    )
    for k, x in witness.items():
        weights = np.array(list(map(Q, x["group_weights"])), object)[ids] / len(words)
        assert all(controls @ weights == controls @ repaired)
        W[k] = weights
    rng = np.random.default_rng(821470)
    n = 2**29
    out = {}
    heads = np.arange(99) // 9
    partial = np.isin(heads, (2, 4, 7, 9))
    for name, weights in W.items():
        dist = tidal.flight.records.ensemble(words, weights)
        for (pop, groups), _ in dist.items():
            assert sum(pop) == 240 and pop[1] + 2 * pop[2] == 4
            assert not set(groups[-1]) & {2, 4, 7, 9} or set(groups[-1]) == {2, 4, 7, 9}
        pre = lab.before(dist)[3]
        a, b, c = [
            pre[mask].sum()
            for mask in (heads == 0, partial, partial & (np.arange(99) % 3 == 1))
        ]
        success = pre @ tidal.fastest_relay().T / 32
        assert all(success == a + (b + c) / 4)
        lo, hi = a + b / 4, a + (1 + Q(4, 240)) * b / 4
        assert lo <= success[0] <= hi
        shots = rng.multinomial(n, np.asarray(pre, float))
        ma, mb = shots[heads == 0].sum() / n, shots[partial].sum() / n
        error = np.sqrt(np.log(2 * 8 / 0.005) / (2 * n))
        lower = ma + mb / 4 - error * 1.25
        upper = ma + mb * (1 + 4 / 240) / 4 + error * (1 + (1 + 4 / 240) / 4)
        out[name] = dict(
            blank=str(a),
            partial=str(b),
            partial_with_single=str(c),
            success=str(success[0]),
            possible_without_coincidence=[str(lo), str(hi)],
            counts=shots.tolist(),
            measured_without_coincidence=[float(lower), float(upper)],
        )
    assert Q(out["source"]["possible_without_coincidence"][0]) > Q(
        out["original"]["possible_without_coincidence"][1]
    )
    assert (
        out["source"]["measured_without_coincidence"][0]
        > out["original"]["measured_without_coincidence"][1]
    )
    assert (
        Q(out["low"]["success"])
        < Q(out["original"]["success"])
        < Q(out["high"]["success"])
    )
    gap = Q(out["source"]["possible_without_coincidence"][0]) - Q(
        out["original"]["possible_without_coincidence"][1]
    )
    # Each trial returns one of three head classes. The two bounded scores
    # are in [0,1]. Reserve one Hoeffding margin for significance and one
    # for power to detect this specified population difference.
    required = int(np.ceil(4 * np.log(100) / float(gap) ** 2))
    return dict(
        shots_per_preparation=n,
        guaranteed_acceptance_gap=str(gap),
        ideal_trials_per_arm_99_significance_and_power=required,
        preparations=out,
    )


def held_witness(h, ends, adj, pair, banks, costs):
    U, V = pair
    A, B = ends
    rng = np.random.default_rng(8362722)
    n = 2**32
    N0 = list(map(int, adj[U].nonzero()[1]))
    N1 = list(map(int, adj[V].nonzero()[1]))
    paths = [(V, k) for k in N0] + [(U, k) for k in N1]
    D = int(h[U, V])
    result = {}
    for name in ("source", "original"):
        p = Q(banks[name]["success"])
        c = 1 / p
        assert c == costs[name]
        possible = banks[name]["measured_without_coincidence"]
        costs_band = [1 / possible[1], 1 / possible[0]]

        def via(i, j):
            return min(int(h[i, A] + h[B, j]), int(h[i, B] + h[A, j]))

        assert via(U, V) + costs_band[0] >= D

        def accel(x):
            return (
                sum(min(h[i, j], via(i, j) + x) ** 2 for i, j in paths)
                - len(paths) * (D**2 + 1)
            ) / D**2

        forecast = list(map(accel, costs_band))
        # Separate route-time records. Waiting for N accepted bits gives
        # a negative-binomial erasure count. Invert its two exact tails.
        low, high = [], []
        receipt = []
        for i, j in paths:
            z = via(i, j)
            if int(h[i, j]) <= z + costs_band[0]:
                low.append(int(h[i, j]))
                high.append(int(h[i, j]))
                receipt.append([i, j, 0])
                continue
            assert z + costs_band[1] < int(h[i, j])
            failures = int(rng.negative_binomial(n, float(p)))
            alpha = 0.005 / (4 * len(paths))
            pl = beta.ppf(alpha, n, failures + 1)
            ph = beta.ppf(1 - alpha, n, failures) if failures else 1.0
            low.append(z + 1 / ph)
            high.append(z + 1 / pl)
            receipt.append([i, j, failures])
        interval = [
            (sum(v * v for v in data) - len(paths) * (D**2 + 1)) / D**2
            for data in (low, high)
        ]
        truth = accel(c)
        assert (
            forecast[0] <= truth <= forecast[1] and interval[0] <= truth <= interval[1]
        )
        result[name] = dict(
            local_forecast=list(map(float, forecast)),
            held_interval=interval,
            held_paths=receipt,
            bits_per_path=n,
            initial_range=D,
        )
    # Difference between the repair and its reference; all signs matter.
    for key in ("local_forecast", "held_interval"):
        a, b = result["source"][key], result["original"][key]
        result[key + "_difference"] = [a[0] - b[1], a[1] - b[0]]
        assert result[key + "_difference"][1] < 0
    return dict(stations=pair, acquisition=result)
