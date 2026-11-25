"""A four-probe measurement using two disjoint local energy differences."""

from fractions import Fraction as Q
from itertools import combinations
import math

import numpy as np


def select_edges(groups, radii, cones, root):
    edges = sorted(
        tuple(sorted((a, b), key=lambda u: radii[u]))
        for group in groups
        for a, b in group
        if root not in (a, b) and radii[a] != radii[b]
    )
    best = None
    total = 0
    for u, v in combinations(edges, 2):
        if (cones[u[0]] | cones[u[1]]) & (cones[v[0]] | cones[v[1]]):
            continue
        total += 1
        rank = (radii[u[0]] + radii[v[0]], u, v)
        if best is None or rank < best:
            best = rank
    return total, None if best is None else (*best[1], *best[2])


def measure_contrasts(
    protocol, groups, radii, cones, age, root, rng, upper, trigger_histories
):
    pairs, cells = select_edges(groups, radii, cones, root)
    record = dict(eligible_edge_pairs=pairs)
    if cells is None:
        return record
    record.update(cells=list(cells), radii=[int(radii[i]) for i in cells], banks={})
    words = np.arange(256)
    swaps = [
        words ^ ((((words >> i) ^ (words >> (i + 1))) & 1) * (3 << i))
        for i in (0, 2, 4, 6)
    ]
    p_init = np.array(
        [
            math.prod(
                (p if int(w) & (1 << k) else 1 - p)
                for k, p in enumerate((1 / 5, 1 / 3) * 4)
            )
            for w in words
        ]
    )
    terminal = sum(((words >> (2 * i)) & 1) << i for i in range(4))
    indices = np.arange(16)
    a = (indices & 1) - ((indices >> 1) & 1)
    b = ((indices >> 2) & 1) - ((indices >> 3) & 1)
    size = protocol["contrast_trials"]
    # Bonferroni includes even deadlines where the disjoint-edge menu is empty.
    m = len(protocol["levels"]) * len(protocol["deadline_sweeps"]) * 3 * 3
    bound = math.sqrt(2 * math.log(2 * m / 0.01) / size)
    for mode in ("release", "earlier_release", "shared_bit"):
        elapsed = age + int(max(radii)) if mode == "earlier_release" else age
        history = [
            (p, [min(t + r, elapsed) for r in record["radii"]], False)
            for t, p in trigger_histories(elapsed)
        ]
        if mode == "shared_bit":
            history = [(Q(1, 2), [age] * 4, flag) for flag in (False, True)]
        remaining, mass = size, Q(1)
        exact = [Q(0)] * 3
        output = np.zeros(16, dtype=np.int64)
        for p, contacts, swap_initial in history:
            count = remaining if p == mass else rng.binomial(remaining, float(p / mass))
            remaining -= count
            mass -= p
            means = [upper(n, swap_initial) for n in contacts]
            u, v = means[0] - means[1], means[2] - means[3]
            exact = [x + p * y for x, y in zip(exact, (u, v, u * v), strict=True)]
            state = rng.multinomial(count, p_init)
            for k, perm in enumerate(swaps):
                if swap_initial:
                    state = state[perm]
                for _ in range(contacts[k]):
                    accepted = rng.binomial(state, 0.25)
                    state = state - accepted + accepted[perm]
            np.add.at(output, terminal, state)
        assert remaining == 0 and output.sum() == size
        observed = [float(output @ d) / size for d in (a, b, a * b)]
        assert all(
            abs(float(x) - y) <= bound for x, y in zip(exact, observed, strict=True)
        )
        covariance = exact[2] - exact[0] * exact[1]
        if mode != "earlier_release":
            assert covariance == 0
        record["banks"][mode] = dict(
            counts=output.tolist(),
            covariance=str(covariance),
            sampled_covariance=observed[2] - observed[0] * observed[1],
            covariance_half_width=bound * (1 + abs(observed[0]) + abs(observed[1]))
            + bound**2,
        )
    return record
