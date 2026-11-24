"""Trace an epoch-stamped material relay through a settling collar population."""

import argparse
import json

import numpy as np
from scipy import sparse

from .experiment import HERE, NSEEDS, layout, geometry, history, calibrate


def adjacency(n, edges):
    a, b = edges.T
    return sparse.csr_matrix(
        (np.ones(2 * len(a), dtype=np.int8), (np.r_[a, b], np.r_[b, a])), shape=(n, n)
    )


def arrival(graphs, stations, emit):
    n = graphs[0].shape[0]
    hits = np.full((n, len(stations)), -1, dtype=np.int32)
    hits[stations, np.arange(len(stations))] = emit
    t = emit
    while True:
        front = hits == t
        reached = (graphs[min(t, len(graphs) - 1)] @ front.astype(np.int8) > 0) & (
            hits < 0
        )
        if not np.any(reached):
            break
        t += 1
        hits[reached] = t
    return hits.T


def receiver_forecasts(scales, first, start, end):
    limit = (
        int(first.max()) + int(np.ceil((end - start) * scales.max() / scales.min())) + 3
    )
    a = np.pad(scales, (0, max(0, limit - len(scales))), mode="edge")
    eta = np.r_[0, np.cumsum(1 / a)]
    target = eta[first] + eta[end] - eta[start]
    return np.interp(target, eta, np.arange(len(eta)))


def measure(level, seed, dilute):
    n, edges, ranks, stations, baseline = layout(level)
    random = np.random.default_rng(seed)
    initial = 2 * random.integers(0, 2, n, dtype=np.int8)
    if dilute:
        initial[random.random(n) < 0.5] = 1
    commit, _, last = history(n, edges, ranks, initial, 0.25, random)
    graphs = [
        adjacency(n, edges[(commit == 0) | (commit > t)]) for t in range(last + 1)
    ]
    final = geometry(n, edges, stations, commit == 0)
    dist = [
        geometry(n, edges, stations, (commit == 0) | (commit > t)) for t in range(last)
    ] + [final]
    far = baseline.max(axis=1)[:, None]
    masks = [
        (baseline > 0) & ((i + 1) * far >= 3 * baseline) & (3 * baseline > i * far)
        for i in range(3)
    ]
    emit = (0, 4, 8, 12, 16, 48, 52)
    reads = {t: arrival(graphs, stations, t) for t in emit}
    metrics = []
    for i, mask in enumerate(masks):
        cohort = mask & (final >= 0)
        scales = (
            None
            if not cohort.any()
            else np.array([ds[cohort].sum() / baseline[cohort].sum() for ds in dist])
        )
        ratios = []
        for e0, e1 in ((0, 4), (4, 8), (8, 12), (12, 16), (48, 52)):
            first, second = reads[e0], reads[e1]
            both = mask & (first >= 0) & (second >= 0)
            assert not np.any(mask & (first < 0) & (second >= 0))
            delays = (second - first)[both]
            assert np.all(delays >= e1 - e0)
            eligible = cohort
            x, y = first[eligible], second[eligible]
            # Every final-visible destination has a viable path at all epochs.
            assert np.all(x >= 0) and np.all(y >= 0)
            pred = receiver_forecasts(scales, x, e0, e1) if scales is not None else None
            if pred is None:
                outside = 0
                mean_error = None
            else:
                outside = int(np.sum(abs(y - pred) > 1 + 1e-8))
                mean_error = float(np.mean(y - pred))
            ratios.append(
                dict(
                    emission=[e0, e1],
                    both=int(both.sum()),
                    lost_second=int(np.sum(mask & (first >= 0) & (second < 0))),
                    neither=int(np.sum(mask & (first < 0))),
                    period_histogram=np.bincount(delays, minlength=5).tolist(),
                    tested=int(len(x)),
                    missed_beyond_one_round=outside,
                    mean_actual_minus_predicted=mean_error,
                )
            )
        metrics.append(
            dict(
                bin=i, scales=None if scales is None else scales.tolist(), pairs=ratios
            )
        )
    rescue = int(np.sum((reads[0] >= 0) & (final < 0)))
    lost = int(np.sum((reads[0] >= 0) & (reads[4] < 0)))
    return (
        dict(
            level=level,
            seed=seed,
            dilute=dilute,
            settled=last,
            metrics=metrics,
            transiently_accessible=rescue,
            first_only=lost,
        ),
        reads,
        commit,
    )


def run():
    records = {}
    cases = []
    for level in range(3):
        for dilute, ns in ((False, NSEEDS), (True, 2)):
            for j in range(ns):
                seed = 280210 + 10 * level + j
                case, reads, commit = measure(level, seed, dilute)
                cases.append(case)
                key = f"{seed}_{int(dilute)}"
                records[key] = np.stack(list(reads.values())).astype(np.int16)
                if dilute:
                    records[key + "_commits"] = commit
    return dict(native_calibration=calibrate(), cases=cases), records


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    report, arrays = run()
    data = json.dumps(report, indent=1) + "\n"
    if args.write:
        (HERE / "signals.json").write_text(data)
        np.savez_compressed(HERE / "signals.npz", **arrays)
    else:
        assert json.loads((HERE / "signals.json").read_text()) == json.loads(data)
        with np.load(HERE / "signals.npz") as saved:
            assert set(arrays) == set(saved)
            for name, arr in arrays.items():
                np.testing.assert_array_equal(saved[name], arr)
        print(
            "Verified every pulse arrival and missing return on 30 settling histories."
        )


if __name__ == "__main__":
    main()
