"""Resolve what extra source readings the exterior material comparison requires."""

from itertools import combinations, product

import networkx as nx
import numpy as np
from sympy import Matrix


def calibrate(ports, adj, cfg, field):
    homes = np.searchsorted(ports, [8, 20, 32, 44])
    blocks = [
        list(np.flatnonzero(np.isin(homes, list(block))))
        for block in nx.biconnected_components(nx.from_numpy_array(adj))
    ]
    assert sorted(b for b in blocks if b) == [[0, 3], [1], [2]]
    s = (cfg[:, None, :] + cfg[None, :, :]).reshape(-1, 4)
    features = [np.ones(len(s), int)]
    names = ["intercept"]
    for h, k in product(range(4), (1, 2)):
        features.append((s[:, h] == k).astype(int))
        names.append(f"home_{ports[homes[h]]}_total_{k}")
    for x, y in product((1, 2), repeat=2):
        features.append(((s[:, 0] == x) & (s[:, 3] == y)).astype(int))
        names.append(f"joint_8_44_totals_{x}_{y}")
    F = np.array(features).T
    # A unit terminal current visits just the unique chain of articulation
    # blocks. Effective resistance is the sum of the per-block resistances.
    # This also holds after applying the linear ordinary release operator.
    columns = list(Matrix(F).rref()[1])
    F = F[:, columns]
    enrolled = list(Matrix(F.T).rref()[1])
    inv = np.asarray(Matrix(F[enrolled]).inv(), float)
    field = field.reshape(len(s), -1)
    coef = inv @ field[enrolled]
    residual = float(abs(F @ coef - field).max())
    assert residual < 2e-9
    return {
        "features": [names[i] for i in columns],
        "calibration_pair_classes": [list(divmod(i, len(cfg))) for i in enrolled],
        "additional_ordered_pairs_checked": len(s) - len(enrolled),
        "maximum_held_residual": residual,
    }


def compare_swaps(cfg, qlo, qhi, field, field_error):
    # Distinct choices per member have disjoint sets of changed home flags.
    # All complete single-member laws AND all four pooled single-home
    # distributions coincide on the two sides of each switch.
    groups = list(combinations(range(len(cfg)), 2))
    choices = []
    for (a, b), (c, d) in combinations(groups, 2):
        if np.any((cfg[a] != cfg[b]) & (cfg[c] != cfg[d])):
            continue
        plus = ((a, c), (b, d))
        minus = ((a, d), (b, c))

        def contrast(x, y, plus=plus, minus=minus):
            return (sum(x[i, j] for i, j in plus) - sum(y[i, j] for i, j in minus)) / 2

        l, h = contrast(qlo, qhi), contrast(qhi, qlo)
        f = contrast(field, field)
        if l <= 0 <= h:
            continue
        direction = 1 if l > 0 else -1
        if direction < 0:
            l, h = -h, -l
        flow = np.array(
            [direction * f - 2 * field_error, direction * f + 2 * field_error]
        )
        ratio = np.divide.outer(flow, [l, h])
        choices.append(
            {
                "class_switch": [a, b, c, d],
                "direction": direction,
                "positive_charge": [l, h],
                "release_per_nat": [ratio.min(), ratio.max()],
            }
        )
    if len(choices) < 2:
        return {"resolved_switches": len(choices)}
    # The comparison is exhaustive and exploratory. Witnesses maximize the
    # certified gap between two release-per-charge intervals at this baseline.
    low = min(choices, key=lambda z: z["release_per_nat"][1])
    high = max(choices, key=lambda z: z["release_per_nat"][0])
    target = min(low["positive_charge"][0], high["positive_charge"][0]) / 4
    return {
        "resolved_switches": len(choices),
        "lower": low,
        "upper": high,
        "shared_charge_nats": target,
        "matched_acceleration_intervals": [
            (np.array(z["release_per_nat"]) * target).tolist() for z in (low, high)
        ],
        "certified_slope_gap": high["release_per_nat"][0] - low["release_per_nat"][1],
    }
