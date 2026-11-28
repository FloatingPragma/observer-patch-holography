"""Heralded binary relays measure transverse ranges on the material-release shell."""

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json

import numpy as np
from scipy.sparse.csgraph import shortest_path

from . import material_release as release

lab = release.lab
HERE = release.HERE
N = 2**38


def success_tickets():
    full, _ = lab.gates()
    table = np.zeros((2, 2, 99), int)
    for col in range(99):
        m, ab = divmod(col, 9)
        a, b = divmod(ab, 3)
        for t in (1, 2):
            for direction, (inputs, axis) in enumerate(
                [((b, 3 * t + b), 0), ((3 * a, 3 * a + t), 1)]
            ):
                probs = (
                    full[:, 9 * m + np.array(inputs)]
                    .reshape(11, 3, 3, 2)
                    .sum(0)
                    .sum(axis=axis)
                )
                first = np.outer(probs[:, 0], probs[:, 1])
                other = first.T
                assert first.sum() == other.sum() == 16
                assert first[other == 0].sum() == other[first == 0].sum()
                table[t - 1, direction, col] = first[other == 0].sum()
    return table


def material_forecast(pre, tickets, t, direction):
    """Predict relay success from transferred load and the target's start."""
    _, half = lab.gates()
    row = np.arange(99)
    head, ab = np.divmod(row, 9)
    receiver = ab % 3 if direction == 0 else ab // 3
    prepared = 9 * head + (3 * t if direction == 0 else t)
    # With the receiving endpoint emptied first, the measured positive load
    # divided by t is a one-sided transfer probability. Express twice this
    # probability in units of 1/16, the denominator of the relay experiment.
    twice_transfer = (4 // t) * (receiver @ half[:, prepared])
    plus = np.isin(head, (2, 4))
    minus = np.isin(head, (7, 9))
    if t == 1:
        sign = plus.astype(int) - minus.astype(int)
        residual = 8 * (receiver > 0) * (sign if direction == 0 else -sign)
    else:
        residual = 8 * (receiver == 1) * (plus | minus)
    np.testing.assert_array_equal(tickets[t - 1, direction], twice_transfer + residual)
    return pre @ twice_transfer / 16, pre @ residual / 16


def balanced_variance(pre, t, rng):
    # A separate experiment measures the signed current on one half pulse
    # after a fair choice of the occupied endpoint. Background stays parked.
    head = pre.reshape(11, 9).sum(1)
    half = lab.gates()[1]
    law = np.full(3, Q(0), object)
    counts = np.zeros(3, dtype=np.int64)
    start_counts = rng.multinomial(2**28, np.asarray(head, float))
    for m, prob in enumerate(head):
        dirs = rng.multinomial(start_counts[m], [0.5, 0.5])
        for d, ab in enumerate((3 * t, t)):
            col = 9 * m + ab
            outcomes = rng.multinomial(dirs[d], half[:, col] / 8)
            for k in np.flatnonzero(half[:, col]):
                change = int(k % 3) - int(ab % 3)
                assert change in (-t, 0, t)
                idx = change // t + 1
                law[idx] += prob * Q(int(half[k, col]), 16)
                counts[idx] += outcomes[k]
    assert sum(law) == 1 and law[0] == law[2]
    variance = law[0] + law[2]
    width = np.sqrt(np.log(2 * 12 / 0.01) / (2 * counts.sum()))
    assert abs((counts[0] + counts[2]) / counts.sum() - float(variance)) < width
    return dict(
        signed_counts=counts.tolist(), variance=str(variance), halfwidth=float(width)
    )


def range_survey(dist, tickets, board, marks, hops, center, rng):
    pre = lab.before(dist)
    A, B = marks[-1]
    # This shell and both detector alphabets were fixed before computing
    # the range changes. Count both directions of every pair with equal weight.
    shell = np.flatnonzero(hops[center] == hops[center, A])
    source, dest = np.where(~np.eye(len(shell), dtype=bool))
    i, j = shell[source], shell[dest]
    out = []
    for t in range(2):
        success = pre @ tickets[t].T / 16
        assert np.all(success[0] == 1) and np.all(success[1:3] == 0)
        p, q = success[3]
        assert min(p, q) > 0
        source = [material_forecast(pre[3], tickets, t + 1, d) for d in (0, 1)]
        assert all(a + b == value for (a, b), value in zip(source, (p, q), strict=True))
        costs = [1 / p, 1 / q]
        variants = np.vstack(
            (
                hops[i, j],
                hops[i, A] + costs[0] + hops[B, j],
                hops[i, B] + costs[1] + hops[A, j],
            )
        )
        choice = variants.argmin(axis=0)
        lengths = variants.min(axis=0)
        exact = sum(lengths) / len(i)
        # Independent shortest-path implementation on the complete relay graph.
        adjacency = board.adjacency().tolil().astype(float)
        for u, v in marks[:2]:
            adjacency[u, v] = adjacency[v, u] = 0
        adjacency[A, B], adjacency[B, A] = map(float, costs)
        graph = shortest_path(adjacency.tocsr(), directed=True)
        np.testing.assert_allclose(graph[center], hops[center], atol=1e-12, rtol=0)
        np.testing.assert_allclose(
            graph[i, j], np.asarray(lengths, float), atol=1e-12, rtol=0
        )
        local = []
        for direction in (0, 1):
            starts = rng.multinomial(N, np.asarray(pre[3], float))
            counts = rng.binomial(starts, tickets[t, direction] / 16).sum()
            local.append(int(counts))
        # Infer ranges from the acquired local counts alone; recalibrate
        # best routes at both ends of the simultaneous probability band.
        error = np.sqrt(np.log(2 * 24 / 0.005) / (2 * N))
        probs = np.asarray(local) / N

        def range_from(probs):
            return float(
                np.minimum.reduce(
                    [
                        hops[i, j],
                        hops[i, A] + 1 / probs[0] + hops[B, j],
                        hops[i, B] + 1 / probs[1] + hops[A, j],
                    ]
                ).mean()
            )

        inferred = [range_from(np.minimum(1, probs + error)), range_from(probs - error)]
        assert np.all(probs > error) and inferred[0] <= float(exact) <= inferred[1]
        # Replay all shortest routes, using fresh renewals on every crossing.
        # A sum of independent geometric waiting times has negative-binomial
        # failure count; deterministic blank legs need no generated draws.
        total = float(exact) * len(i)
        variance = 0
        crossings = []
        for k, prob in enumerate((p, q), 1):
            count = int(np.sum(choice == k))
            crossings.append(count)
            if count:
                successes = N * count
                failures = rng.negative_binomial(successes, float(prob))
                total += (failures - float(successes * (1 - prob) / prob)) / N
                variance += float(count * (1 - prob) / prob**2) / N
        held = total / len(i)
        sigma = np.sqrt(variance) / len(i)
        # Six preparations times two receivers, Chebyshev and a union bound.
        width = np.sqrt(12 / 0.005) * sigma
        assert abs(held - float(exact)) <= width + 1e-14
        variance_read = balanced_variance(pre[3], t + 1, rng)
        assert 4 * Q(variance_read["variance"]) == sum(a for a, b in source)
        out.append(
            dict(
                balanced_first_current=variance_read,
                third_acceptance=list(map(str, (p, q))),
                third_attempts=list(map(str, costs)),
                twice_empty_target_transfer=[str(a) for a, b in source],
                occupied_target_correction=[str(b) for a, b in source],
                local_successes=local,
                local_range_interval=inferred,
                ordered_shell_pairs=len(i),
                crossing_pairs=crossings,
                mean_range=str(exact),
                held_range=held,
                simultaneous_99_halfwidth=width,
                shell_stations=len(shell),
            )
        )
    return out


def run():
    _, words, cells, *_ = lab.energy.pair.setup()
    board, _, marks, *_ = lab.energy.returned(words)
    _, hops, _ = release.records.thermal_radar.chart(board, marks)
    w = np.full(len(words), Q(1, len(words)), object)
    arms = dict(
        original=release.records.ensemble(words, w),
        source=release.records.ensemble(
            words, lab.energy.least_change(cells) / len(words)
        ),
        single=release.records.ensemble(words, w, single=1),
        ordinary=release.records.ensemble(words, w, extra=1),
    )
    releases = json.loads((HERE / "material_release.json").read_bytes())
    arms["matched_positive"], arms["matched_reference"] = [
        release.mix(arms, weights) for weights in releases["matched_mixtures"]
    ]
    tickets = success_tickets()
    rng = np.random.default_rng(5820364)
    data = {
        name: range_survey(dist, tickets, board, marks, hops, 219, rng)
        for name, dist in arms.items()
    }

    def difference(a, b):
        return [
            Q(x["mean_range"]) - Q(y["mean_range"])
            for x, y in zip(data[a], data[b], strict=True)
        ]

    contrasts = {
        name: difference(name, "original") for name in ("source", "single", "ordinary")
    }
    contrasts["matched"] = difference("matched_positive", "matched_reference")
    for tag, a, b in [
        ("repair", "source", "original"),
        ("matched", "matched_positive", "matched_reference"),
        ("ordinary_single", "single", "original"),
        ("ordinary_double", "ordinary", "original"),
    ]:
        for t in (1, 2):
            for d, launch in enumerate((f"{t},0", f"0,{t}")):
                p1 = Q(data[a][t - 1]["twice_empty_target_transfer"][d])
                p0 = Q(data[b][t - 1]["twice_empty_target_transfer"][d])
                j = Q(
                    releases["schedules"]["forward_reverse"]["contrasts"][tag][launch][
                        "current_per_unit"
                    ][0]
                )
                assert p1 - p0 == (2 if d == 0 else -2) * j
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (HERE / "RELAY_RANGE.md", HERE / "material_release.json")
        },
        shell_radius=7,
        ordinary_blank_attempts=1,
        shots_per_link_or_ordered_pair=N,
        arms=data,
        mean_range_changes={k: list(map(str, v)) for k, v in contrasts.items()},
        repair_per_ordinary_double=[
            str(x / y) if y else None
            for x, y in zip(contrasts["source"], contrasts["ordinary"], strict=True)
        ],
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    path = HERE / "relay_range.json"
    if args.write:
        path.write_text(raw)
    else:
        assert path.read_text() == raw
    print(
        {k: list(map(float, map(Q, v))) for k, v in out["mean_range_changes"].items()}
    )
