"""Calibrated current parity and wall position on the native ten-port belt."""

import argparse
import json
from collections import deque
from math import factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, expm
from scipy.stats import beta

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo.candidate import KAPPA, native_law
from ...quark_winding.belt_walls import belt

HERE = Path(__file__).parent
WAITS = [0, 1, 2, 4, 8, 16, 32, 64]
SHOTS = 2**30
SEED = 122013
MODES = ["native", "wall_only", "exchange_only", "parked", "unwired"]


def crossing(a, b, cut):
    """Signed crossing of the cut just before site cut along the short arc."""
    x, y = (a - cut) % 10, (b - cut) % 10
    delta = (y - x + 5) % 10 - 5
    return (x + delta - y) // 10


def edges(sector, cut):
    left, right = np.nonzero(np.triu(sector.K, 1))
    movements = []
    for i, j in zip(left, right, strict=True):
        src, dst = map(set, (sector.states[i], sector.states[j]))
        (a,) = src - dst
        (b,) = dst - src
        movements.append(crossing(a, b, cut))
    return left, right, np.array(movements)


def counter_signs(sec, cut):
    i, j, current = edges(sec, cut)
    signs = np.ones_like(sec.K)
    signs[i, j] = signs[j, i] = 1 - 2 * abs(current)
    return signs


def bits(sec):
    return [sum(1 << v for v in w) for w in sec.states]


def tour(sec, cut):
    walls = belt.selected(sec)
    start = int(walls[0])
    pos = {int(i): belt.wall(sec.states[i], 10)[2] for i in walls}
    i = start
    trace = []
    while True:
        (j,) = [j for j in pos if sec.K[i, j] != 0 and (pos[j] - pos[i]) % 10 == 3]
        src, dst = map(set, (sec.states[i], sec.states[j]))
        (a,) = src - dst
        (b,) = dst - src
        trace.append(
            [int(i), int(j), crossing(a, b, cut), int((pos[i] - cut) % 10 + 3 >= 10)]
        )
        i = j
        if i == start:
            break
        if len(trace) >= 10:
            raise ValueError("The directed wall cycle did not close")
    return trace


def shortest(sec, cut, allowed):
    i, j, current = edges(sec, cut)
    graph = [[] for _ in sec.states]
    for a, b, s in zip(i, j, current, strict=True):
        if a in allowed and b in allowed:
            graph[a].append((int(b), int(abs(s))))
            graph[b].append((int(a), int(abs(s))))
    start = int(belt.selected(sec)[0])
    queue = deque([(start, 0)])
    parent = {(start, 0): None}
    dist = {(start, 0): 0}
    count = {(start, 0): 1}
    while queue:
        v, k = queue.popleft()
        for w, c in graph[v]:
            target = (w, k ^ c)
            step = dist[v, k] + 1
            if target not in parent:
                parent[target] = (v, k)
                dist[target] = step
                count[target] = 0
                queue.append(target)
            if dist[target] == step:
                count[target] += count[v, k]
    end = (start, 1)
    if end not in dist:
        return None
    path = []
    at = end
    while at is not None:
        path.append(list(at))
        at = parent[at]
    length = dist[end]
    return {
        "steps": length,
        "paths": count[end],
        "coefficient_numerator": count[end],
        "coefficient_denominator": 2**length * factorial(length),
        "word_and_pointer_path": path[::-1],
    }


def flight(h, source):
    vals, vecs = eigh(h)
    return (np.exp(-1j * np.outer(WAITS, vals)) * vecs[source]) @ vecs.T


def probabilities(sec, mode, cut):
    walls = belt.selected(sec)
    start = int(walls[0])
    selected = walls if mode == "wall_only" else np.arange(len(sec.states))
    h = sec.hamiltonian(0 if mode == "exchange_only" else KAPPA)
    if mode == "parked":
        h = np.diag(np.diag(h))
    h = h[np.ix_(selected, selected)]
    sign = counter_signs(sec, cut)[np.ix_(selected, selected)]
    source = int(np.where(selected == start)[0][0])
    plus = flight(h, source)
    minus = plus if mode == "unwired" else flight(h * sign, source)
    prob = np.zeros((len(WAITS), len(sec.states), 2))
    prob[:, selected, 0] = abs((plus + minus) / 2) ** 2
    prob[:, selected, 1] = abs((plus - minus) / 2) ** 2
    np.testing.assert_allclose(prob.sum(axis=(1, 2)), 1, atol=2e-13, rtol=0)
    return prob


def native_counter():
    # The installed parity bit couples by two controlled flips to one
    # endpoint's occupancy. Apply the full native 99-state contact to the
    # blank-memory input to calibrate it on an independently held ordinary unit.
    exchange, contrast = native_law.local_term_parts()
    h = (exchange + KAPPA * contrast)[:9, :9]
    block = expm(-1j * np.pi * h)
    checked = []
    for input_word in [(1, 0), (0, 1), (0, 0)]:
        start = native_law.SRC.index(input_word)
        tagged = np.zeros((9, 2), complex)
        for idx, word in enumerate(native_law.SRC):
            tagged[idx, (input_word[0] - word[0]) % 2] = block[idx, start]
        checked.append([list(input_word), float((abs(tagged[:, 1]) ** 2).sum())])
    return checked


def cp(k, n, alpha):
    return [
        float(beta.ppf(alpha / 2, k, n + 1 - k)) if k else 0.0,
        float(beta.ppf(1 - alpha / 2, k + 1, n - k)) if k < n else 1.0,
    ]


def acquire():
    rng = np.random.default_rng(SEED)
    output = {}
    receipt = {
        "native_port_order": list(next(iter(belt.belts().values()))),
        "waits": WAITS,
        "modes": MODES,
        "shots": SHOTS,
        "seed": SEED,
        "calibration": native_counter(),
        "sectors": {},
    }
    for n in [3, 4]:
        sec = BinarySector(belt.ring(10), n)
        walls = belt.selected(sec)
        source = int(walls[0])
        banks = np.zeros(
            (len(MODES), 10, len(WAITS), len(sec.states), 2), dtype="int64"
        )
        outcomes = {}
        for m, mode in enumerate(MODES):
            # Moving the cut across a set of ports changes the parity by the
            # XOR of their initial and final occupations. The state gauge has
            # no time, width or wall-charge parameters.
            base = probabilities(sec, mode, 0)
            words = np.array(bits(sec))
            for cut in range(10):
                region = (1 << cut) - 1
                flip = np.array(
                    [
                        int((int(w ^ words[source]) & region).bit_count() % 2)
                        for w in words
                    ],
                    bool,
                )
                prob = base.copy()
                if mode != "unwired":
                    prob[:, flip] = prob[:, flip, ::-1]
                for t, p in enumerate(prob):
                    p = np.round(p.ravel(), 14)
                    banks[m, cut, t] = rng.multinomial(SHOTS, p / p.sum()).reshape(
                        -1, 2
                    )
            odd = banks[m, :, :, source, 1]
            outcomes[mode] = {
                "odd_source_cut0": odd[0].tolist(),
                "min_across_cuts": odd.min(axis=0).tolist(),
                "max_across_cuts": odd.max(axis=0).tolist(),
            }
        output[str(n)] = banks
        trace = tour(sec, 0)
        alpha = 0.01 / (2 * 10 * len(WAITS))
        held = []
        native, disabled = (
            banks[0, :, :, source, 1],
            banks[MODES.index("unwired"), :, :, source, 1],
        )
        for t in range(len(WAITS)):
            bounds = [cp(int(x), SHOTS, alpha) for x in native[:, t]]
            empty = [cp(int(x), SHOTS, alpha) for x in disabled[:, t]]
            held.append(
                {
                    "time": WAITS[t],
                    "lowest_lower": min(x[0] for x in bounds),
                    "highest_upper": max(x[1] for x in bounds),
                    "resolved_on_all_cuts": all(
                        lit[0] > dark[1]
                        for lit, dark in zip(bounds, empty, strict=True)
                    ),
                }
            )
        receipt["sectors"][str(n)] = {
            "source": list(sec.states[source]),
            "word_masks": bits(sec),
            "cut0_wall_tour": trace,
            "min_full": shortest(sec, 0, set(range(len(sec.states)))),
            "min_only_walls": shortest(sec, 0, set(walls)),
            "same_word": outcomes,
            "native_intervals": held,
        }
    return receipt, output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    r, c = acquire()
    js = HERE / "results.json"
    ar = HERE / "counts.npz"
    if args.write:
        js.write_text(json.dumps(r, indent=2, allow_nan=False) + "\n")
        np.savez_compressed(ar, **c)
    else:
        if json.loads(js.read_bytes()) != r:
            raise ValueError("Summary replay changed")
        with np.load(ar) as old:
            for key, counts in c.items():
                np.testing.assert_array_equal(counts, old[key])
    for key, value in r["sectors"].items():
        print(
            key,
            value["min_full"]["steps"],
            value["min_full"]["paths"],
            value["same_word"]["native"]["odd_source_cut0"],
        )


if __name__ == "__main__":
    main()
