"""Read wall and separated-cluster energies on the same open three-belt guide."""

import argparse
import json
from collections import Counter
from itertools import permutations, product
from math import comb, sqrt
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.candidate_laws.interval_echo.openwave.quark_winding.belt_walls import (
    belt,
    composite,
)

ROOT = Path(__file__).parent
SHOTS = 16384


def is_wall(word):
    sites = [j for j in range(10) if word & (1 << j)]
    if len(sites) not in (3, 4):
        return False
    gaps = sorted((b - a) % 10 for a, b in zip(sites, sites[1:] + sites[:1]))
    return gaps in ([3, 3, 4], [1, 3, 3, 3])


def occupied_pairs(word, edges):
    return sum(bool(word & (1 << a)) and bool(word & (1 << b)) for a, b in edges)


def boundary_groups(cell, ends, wall_only):
    """Words with these boundary bits and count have identical bridge interactions."""
    groups = {}
    seams = belt.ring(10).seams
    for word in range(1 << 10):
        if wall_only and not is_wall(word):
            continue
        global_word = word << (10 * cell)
        key = word.bit_count(), global_word & ends
        cost = occupied_pairs(word, seams)
        count, best, witness = groups.get(key, (0, -1, 0))
        if cost > best:
            best, witness = cost, global_word
        groups[key] = count + 1, best, witness
    return groups


def maxima(bridges, wall_only):
    ends = sum(1 << i for i in {p for edge in bridges for p in edge})
    groups = [boundary_groups(cell, ends, wall_only) for cell in range(3)]
    stats = {n: {"count": 0, "maximum": -1, "witness": 0} for n in (10, 11)}
    for entries in product(*(g.items() for g in groups)):
        keys, values = zip(*entries)
        total = sum(n for n, _ in keys)
        if total not in stats:
            continue
        info = stats[total]
        (c1, c2, c3), costs, witnesses = zip(*values)
        info["count"] += c1 * c2 * c3
        bits = sum(mask for _, mask in keys)
        pairs = sum(costs) + occupied_pairs(bits, bridges)
        if pairs > info["maximum"]:
            info.update(maximum=pairs, witness=sum(witnesses))
    return stats


def low_band_bound(loads, edges, wall_max, threshold):
    wall_words = {
        n: [w for w in range(1 << 10) if w.bit_count() == n and is_wall(w)]
        for n in (3, 4)
    }
    walls = {
        a | (b << 10) | (c << 20)
        for order in set(permutations(loads))
        for a, b, c in product(*(wall_words[n] for n in order))
    }
    moves = composite.hop_list(edges)
    parents = Counter()
    max_exits = 0
    for wall in walls:
        exits = [
            word for word in composite.neighbours(wall, moves) if word not in walls
        ]
        parents.update(exits)
        max_exits = max(max_exits, len(exits))
    max_parents = max(parents.values())
    norm_bound = sqrt(max_exits * max_parents) / 2
    lower = -2 * candidate.KAPPA * wall_max
    if lower <= threshold:
        raise ValueError("Wall-block energy bounds overlap the trial threshold")
    return {
        "wall_block_energy_floor": lower,
        "distinct_complement_words_on_cut": len(parents),
        "maximum_exits_per_wall": max_exits,
        "maximum_wall_parents": max_parents,
        "off_block_norm_upper": norm_bound,
        "any_energy_filtered_wall_probability_upper": min(
            1.0, (norm_bound / (lower - threshold)) ** 2
        ),
    }


def comparison_word(loads, bridges, wall):
    result = 0
    for cell, n in enumerate(loads):
        blocked = {p % 10 for edge in bridges for p in edge if p // 10 == cell}
        if wall:
            words = (w for w in range(1 << 10) if w.bit_count() == n and is_wall(w))
        else:
            words = (
                sum(1 << ((start + j) % 10) for j in range(n)) for start in range(10)
            )
        word = next(w for w in words if all(not (w & (1 << p)) for p in blocked))
        result |= word << (10 * cell)
    return result


def energy_records(word, edges, rng):
    pairs = occupied_pairs(word, edges)
    records = []
    boundaries = 0
    for a, b in edges:
        left, right = (word >> a) & 1, (word >> b) & 1
        boundaries += left != right
        odd = int(rng.binomial(SHOTS, 0.5)) if left != right else 0
        both = SHOTS if left and right else 0
        records.append([both, odd, SHOTS - both - odd])
    sums = np.sum(records, axis=0)
    reads = {}
    for coupling in (candidate.KAPPA, 0, -candidate.KAPPA):
        reads[str(coupling)] = {
            "mean": boundaries / 2 - 2 * coupling * pairs,
            "observed": float((sums[1] - 2 * coupling * sums[0]) / SHOTS),
        }
    return {
        "word": word,
        "pair_count": pairs,
        "cut_size": boundaries,
        "energies": reads,
        "counts": records,
    }


def run():
    geometry, edges = composite.path()
    bridges = geometry["local_bridges"]
    all_states, wall_states = (maxima(bridges, only) for only in (False, True))
    trials = []
    rng = np.random.default_rng(131305)
    for loads in ((3, 3, 4), (3, 4, 4)):
        n = sum(loads)
        sector, wall = all_states[n], wall_states[n]
        if sector["count"] != comb(30, n) or wall["count"] != 3000:
            raise ArithmeticError("Incomplete wall or full-sector count")
        for group in (sector, wall):
            assert occupied_pairs(group["witness"], edges) == group["maximum"]
        trial = energy_records(comparison_word(loads, bridges, False), edges, rng)
        source = energy_records(comparison_word(loads, bridges, True), edges, rng)
        kappa = candidate.KAPPA
        threshold = trial["energies"][str(kappa)]["mean"]
        ceiling = (sector["maximum"] + threshold / (2 * kappa)) / (
            sector["maximum"] - wall["maximum"]
        )
        trials.append(
            {
                "loads": loads,
                "sector": sector,
                "walls": wall,
                "separated": trial,
                "three_walls": source,
                "maximum_wall_probability_at_or_below_separated_mean": ceiling,
                "energy_filter": low_band_bound(
                    loads, edges, wall["maximum"], threshold
                ),
            }
        )
    return {
        "shots_per_seam": SHOTS,
        "seed": 131305,
        "edges": edges,
        "geometry": geometry,
        "sectors": trials,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.output:
        with args.output.open("x") as f:
            f.write(result)
    elif json.loads(result) != json.loads((ROOT / "results.json").read_text()):
        raise ValueError("Wall-energy records differ")
    else:
        print(
            "Reproduced the full-sector energy bounds and the two local energy comparisons."
        )


if __name__ == "__main__":
    main()
