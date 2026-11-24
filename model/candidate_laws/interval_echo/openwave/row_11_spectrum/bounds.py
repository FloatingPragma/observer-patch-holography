"""Integer composition bounds on entire, unshuttered native federations."""

import argparse
import json

import networkx as nx
import numpy as np

from model.executable_model import cell, federation, rotations
from .measure import HERE, word_table


def layout_type(gate):
    return min(tuple(sorted(p[i] for i in gate)) for p in rotations())


def local_certificate(layouts):
    board = cell()
    _, _, words = word_table(board.n)
    s = words.astype(np.int16) - 1
    q, d = s.sum(axis=1), (s**2).sum(axis=1)
    internal = sum((s[:, i] - s[:, j]) ** 2 for i, j in board.seams).astype(np.int16)
    holes = np.sum(s == -1, axis=1)
    gains = np.sum(s == 1, axis=1)
    cls = holes * 13 + gains
    local_min = np.full((13, 13), 10000, dtype=int)
    descriptions = []
    for gate in sorted(layouts):
        b = np.ones(12, dtype=int)
        b[list(gate)] = 0
        cost = internal + (s**2 + 2 * s) @ b
        slack = int(np.min(4 * cost - 6 * q - 3 * d))
        if slack < 0:
            raise ValueError("Proposed owner inequality has a counterexample")
        groups = np.full(169, 10000, dtype=int)
        np.minimum.at(groups, cls, cost)
        local_min = np.minimum(local_min, groups.reshape(13, 13))
        descriptions.append(
            dict(gate=gate, words=len(words), min_slack_times_four=slack)
        )
    # The three largest single-port corrections give the most favorable gateways of any geometry.
    largest_three = 3 * np.minimum(3, gains) - np.maximum(0, holes - 9)
    best_possible_cost = internal + 3 * gains - holes - largest_three
    assert np.min(4 * best_possible_cost - 6 * q - 3 * d) == 0
    best_groups = np.full(169, 10000, dtype=int)
    np.minimum.at(best_groups, cls, best_possible_cost)
    return (
        descriptions,
        small_bound(local_min),
        small_bound(best_groups.reshape(13, 13)),
    )


def small_bound(local_min):
    # Permit any layout in each occupied owner and omit all external contacts.
    lower = {(0, 0): 0}
    for d in range(1, 8):
        for z in range(d + 1):
            g = d - z
            lower[z, g] = min(
                lower[z - i, g - j] + int(local_min[i, j])
                for i in range(z + 1)
                for j in range(g + 1)
                if i + j
            )
    return {d: lower[(d + 1) // 2, (d - 1) // 2] for d in (1, 3, 5, 7)}


def competitors(board):
    nb = nx.Graph()
    nb.add_nodes_from(range(board.cells))
    nb.add_edges_from((a // 12, b // 12) for a, b in board.external)
    depth = nx.single_source_shortest_path_length(nb, 0)
    owners = sorted(nb, key=lambda j: (depth[j], j))
    order = np.array([12 * j + i for j in owners for i in range(12)])
    a, b = np.array(board.seams).T
    idle = 6 - board.degree
    s = np.zeros(board.n, dtype=np.int16)
    s[order[0]] = -1
    costs = []
    for p in range(board.n // 2):
        if p:
            s[order[p]] = -1
            s[order[-p]] = 1
        cost = int(np.sum((s[a] - s[b]) ** 2) + np.dot(idle, s**2 + 2 * s))
        twice_kinetic = int(np.sum(s[a] != s[b]))
        costs.append(16 * cost + twice_kinetic)
    return costs


def run():
    boards = [federation(i) for i in range(3)]
    layouts = {
        layout_type(tuple(np.flatnonzero(board.degree[i : i + 12] == 6)))
        for board in boards
        for i in range(0, board.n, 12)
    }
    desc, small, any_small = local_certificate(layouts)
    out = []
    for board in boards:
        costs = competitors(board)
        n = board.n
        floor = min(small.values())
        max_d = n - 1
        ds = np.arange(1, max_d + 1, 2)
        lb = [small[d] if d < 9 else 0.75 * d - 1.5 for d in ds]
        lower = 8 * min(lb)
        assert lower == 8 * floor
        uniform = 8 * (6 * ds - 1.5 + (21 / 4) * (ds - 1) / (n - 1))
        best = np.minimum(uniform, np.array(costs) / 2)
        min_d = {}
        for target in (207, 3477):
            match = ds[best >= target * lower]
            min_d[str(target)] = None if not len(match) else int(match[0])
        out.append(
            dict(
                cells=board.cells,
                ports=n,
                fillings=len(ds),
                lowest_floor=lower,
                every_minimum_ceiling=float(best.max()),
                mass_ratio_ceiling=float(max(best) / lower),
                uniform_only_ceiling=float(max(uniform)),
                minimum_defect_count=min_d,
                twice_sharp_energies=costs,
            )
        )
    any_floor = 8 * min(min(any_small.values()), 0.75 * 9 - 1.5)
    needed = {
        str(t): next(
            d for d in range(1, 10001, 2) if 8 * (6 * d + 3.75) > t * any_floor
        )
        for t in (207, 3477)
    }
    return dict(
        owner_types=desc,
        few_defect_potential_floors=small,
        boards=out,
        arbitrary_three_gateways=dict(
            few_defect_potential_floors=any_small,
            universal_energy_floor=any_floor,
            minimum_defects_for_ratio=needed,
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    d = json.loads(json.dumps(run()))
    p = HERE / "bounds.json"
    if parser.parse_args().write:
        p.write_text(json.dumps(d, indent=1) + "\n")
    elif d != json.loads(p.read_text()):
        raise ValueError("Exact federation bounds differ")
    print(
        "Verified the full-board variational inequalities in every negative-unit composition."
    )


if __name__ == "__main__":
    main()
