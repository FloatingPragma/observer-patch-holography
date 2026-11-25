"""Exact release histories and an independent original-denominator count experiment."""

import argparse
import json
import math
from fractions import Fraction as Q
from functools import cache
from itertools import combinations, product
from pathlib import Path

import numpy as np

from contrast import measure_contrasts
from relay import HERE, gates, propagate, topology


def past_cones(groups, size, age):
    # Propagate sets of initial cell labels through every scheduled matching.
    cones = [1 << i for i in range(size)]
    for _ in range(age):
        for edges in groups:
            for a, b in edges:
                cones[a] = cones[b] = cones[a] | cones[b]
    return cones


def trigger_histories(age):
    for t in range(1, age + 1):
        yield t, Q(3, 4) ** (t - 1) / 4
    yield age + 1, Q(3, 4) ** age


@cache
def upper(n, warm=False):
    return Q(4, 15) + (1 if warm else -1) * Q(1, 15 * 2**n)


def histories(bank, age, di, dj):
    # Each output gives probability, contact counts, initial swap flags and hits.
    if bank in ("release", "earlier_release"):
        for t, p in trigger_histories(age):
            arrivals = (t + di, t + dj)
            yield (
                p,
                tuple(min(age, v) for v in arrivals),
                (False, False),
                tuple(v <= age for v in arrivals),
            )
    else:
        bits = {
            "shared_bit": [(i, i) for i in range(2)],
            "independent_bits": list(product(range(2), repeat=2)),
            "cut_preparation": [(0, 0)],
        }[bank]
        for b in bits:
            yield Q(1, len(bits)), (age, age), b, (False, False)


def exact_table(bank, age, di, dj):
    table = [[Q(0) for _ in range(4)] for _ in range(4)]
    for p, counts, flags, hits in histories(bank, age, di, dj):
        lo, hi = (upper(n, b) for n, b in zip(counts, flags, strict=True))
        for i, j in product(range(2), repeat=2):
            table[2 * hits[0] + i][2 * hits[1] + j] += (
                p * (lo if i else 1 - lo) * (hi if j else 1 - hi)
            )
    assert sum(map(sum, table)) == 1
    return table


def moments(table):
    a = sum(sum(row) for i, row in enumerate(table) if i % 2)
    b = sum(sum(row[1::2]) for row in table)
    ab = sum(sum(row[1::2]) for i, row in enumerate(table) if i % 2)
    return a, b, ab, ab - a * b


def centered_power(radii, age):
    # Exact mean squared spatial contrast of the conditional means after
    # centering each site's ensemble mean. Excludes intrinsic detector noise.
    distances, counts = np.unique(radii, return_counts=True)
    weights = [Q(int(n), len(radii)) for n in counts]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b, strict=True))

    def spatial_variance(v):
        return dot(weights, [x * x for x in v]) - dot(weights, v) ** 2

    mean = [Q(0)] * len(weights)
    second = Q(0)
    for t, p in trigger_histories(age):
        v = [upper(min(t + int(d), age)) for d in distances]
        mean = [x + p * y for x, y in zip(mean, v, strict=True)]
        second += p * spatial_variance(v)
    return second - spatial_variance(mean)


# State bit order: cold probe, warm partner, cold probe, warm partner.
WORDS = np.arange(16)
SWAPS = [WORDS ^ ((((WORDS >> i) ^ (WORDS >> (i + 1))) & 1) * (3 << i)) for i in (0, 2)]
INITIAL = np.array(
    [
        math.prod(
            p if w & (1 << j) else 1 - p
            for j, p in enumerate((1 / 5, 1 / 3, 1 / 5, 1 / 3))
        )
        for w in WORDS
    ]
)


def count_table(bank, age, di, dj, rng, size):
    table = np.zeros((4, 4), dtype=np.int64)
    remaining = size
    remaining_prob = Q(1)
    for p, contacts, flags, hits in histories(bank, age, di, dj):
        # Sequential conditional binomials sample the root clock (or prior bits).
        count = (
            remaining
            if p == remaining_prob
            else rng.binomial(remaining, float(p / remaining_prob))
        )
        remaining -= count
        remaining_prob -= p
        words = rng.multinomial(count, INITIAL)
        for k, perm in enumerate(SWAPS):
            if flags[k]:
                words = words[perm]
            for _ in range(contacts[k]):
                accepted = rng.binomial(words, 0.25)
                words = words - accepted + accepted[perm]
        np.add.at(
            table, (2 * hits[0] + (WORDS & 1), 2 * hits[1] + ((WORDS >> 2) & 1)), words
        )
    assert remaining == 0 and table.sum() == size
    return table


def selected_pairs(dist, cones, root):
    pairs = [p for p in combinations(range(len(dist)), 2) if root not in p]
    disjoint = [(i, j) for i, j in pairs if not cones[i] & cones[j]]
    near = min(pairs, key=lambda p: (sum(dist[root, k] for k in p), p))
    far = min(
        disjoint or pairs, key=lambda p: (-dist[p], sum(dist[root, k] for k in p), p)
    )
    return len(disjoint), dict(near=near, far=far)


def measure(protocol):
    rng = np.random.default_rng(protocol["sample_seed"])
    contrast_rng = np.random.default_rng(protocol["sample_seed"] + 1)
    root, size = protocol["root"], protocol["precision_originals_per_pair"]
    banks = (
        "release",
        "earlier_release",
        "shared_bit",
        "independent_bits",
        "cut_preparation",
    )
    ntests = (
        len(protocol["levels"]) * len(protocol["deadline_sweeps"]) * 2 * len(banks) * 3
    )
    bound = math.sqrt(math.log(2 * ntests / 0.01) / (2 * size))
    result = dict(
        registers=gates(),
        nominal_trials_per_pair_bank=size,
        simultaneous_99pct_moment_half_width=bound,
        boards=[],
    )
    for level in protocol["levels"]:
        board, groups, dist = topology(level)
        radii = dist[root]
        depth = int(radii.max())
        arrivals = propagate(groups, board.cells, root, depth)
        assert np.array_equal(arrivals, radii)
        cut = propagate(groups, board.cells, root, depth, cut=True)
        assert set(np.flatnonzero(cut >= 0)) == {root}
        eager = propagate(groups, board.cells, root, depth, stagger=False)
        assert np.all(eager <= radii) and np.any(eager < radii)
        nonroot = np.delete(radii, root)
        record = dict(
            cells=board.cells,
            root=root,
            eccentricity=depth,
            matches=[len(g) for g in groups],
            distance_by_cell=radii.tolist(),
            eager_first_sweep=int(sum(eager == 1)),
            deadlines=[],
        )
        result["boards"].append(record)
        for h in protocol["deadline_sweeps"]:
            cones = past_cones(groups, board.cells, h)
            disjoint, pairs = selected_pairs(dist, cones, root)
            deadline = dict(
                sweeps=h,
                disjoint_pairs=disjoint,
                total_pairs=math.comb(board.cells - 1, 2),
                centered_release_power=str(centered_power(nonroot, h)),
                centered_earlier_release_power=str(centered_power(nonroot, depth + h)),
                centered_shared_bit_power="0",
                pairs=[],
            )
            deadline["local_differences"] = measure_contrasts(
                protocol,
                groups,
                radii,
                cones,
                h,
                root,
                contrast_rng,
                upper,
                trigger_histories,
            )
            record["deadlines"].append(deadline)
            for label, (i, j) in pairs.items():
                di, dj = int(radii[i]), int(radii[j])
                row = dict(
                    selection=label,
                    cells=[i, j],
                    radii=[di, dj],
                    separation=int(dist[i, j]),
                    common_initial_cells=(cones[i] & cones[j]).bit_count(),
                    banks={},
                )
                deadline["pairs"].append(row)
                for bank in banks:
                    age = depth + h if bank == "earlier_release" else h
                    exact = exact_table(bank, age, di, dj)
                    means = moments(exact)
                    if bank == "shared_bit":
                        assert means[-1] == Q(1, 225 * 4**h)
                    if bank in ("independent_bits", "cut_preparation") or (
                        bank == "release" and max(di, dj) >= h - 1
                    ):
                        assert means[-1] == 0
                    counts = count_table(bank, age, di, dj, rng, size)
                    observed = moments(counts / size)
                    assert (
                        max(
                            abs(float(a) - b)
                            for a, b in zip(means[:3], observed[:3], strict=True)
                        )
                        <= bound
                    )
                    row["banks"][bank] = dict(
                        probabilities=[[str(x) for x in v] for v in exact],
                        counts=counts.tolist(),
                        covariance=str(means[-1]),
                        sampled_covariance=float(observed[-1]),
                        covariance_half_width=bound * (1 + observed[0] + observed[1])
                        + bound**2,
                    )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out", type=Path, help="destination JSON (default checks the archived record)"
    )
    args = parser.parse_args()
    result = measure(json.loads((HERE / "protocol.json").read_text()))
    data = json.dumps(result, indent=1) + "\n"
    if args.out:
        args.out.write_text(data)
    elif data != (HERE / "results.json").read_text():
        raise SystemExit("archived run differs")
