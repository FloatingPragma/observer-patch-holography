"""Acquire and score terminal repair flags on a closed cycle of owners."""

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import product
from math import log, sqrt
from pathlib import Path

import numpy as np

from .apparatus import board_layout, load_table, local_calibration

HERE = Path(__file__).resolve().parent
SHOTS = 65536
SEED = 1282611
MENUS = (
    ("half", "shared", "repair"),
    ("half", "split", "repair"),
    ("third", "shared", "repair"),
    ("third", "split", "repair"),
    ("half", "shared", "unchanged_loads"),
    ("half", "shared", "held_parked"),
)


def tapes(layout, private, sharing, mode):
    """Only preparation sees the source bits. The last dimension is two distinct meter bits."""
    _, writers = local_calibration()
    pairs = 2 * private[:, :, 0] + private[:, :, 1]
    endpoints = np.stack(
        [
            load_table(s, sharing == "shared")[pairs[:, i]]
            for i, s in enumerate(layout["stations"])
        ],
        axis=1,
    )
    n = endpoints.shape[1]
    reads = []
    for i in range(n):
        table = writers[
            "parked"
            if mode == "held_parked" and i == n - 1
            else "repair"
            if mode == "held_parked"
            else mode
        ]
        reads.append(table[endpoints[:, i, 1], endpoints[:, (i + 1) % n, 0]])
    powers = 1 << np.arange(n)
    return np.stack(reads, axis=1).transpose(0, 2, 1) @ powers


def acquire():
    rng = np.random.Generator(np.random.PCG64(SEED))
    grid, _ = local_calibration()
    panels = []
    for level in range(3):
        layout = board_layout(level)
        n = len(layout["owners"])
        for split, sharing, mode in MENUS:
            private = rng.integers(4, size=(SHOTS, n, 2), dtype=np.uint8) < grid[split]
            words = tapes(layout, private.astype(np.int8), sharing, mode)
            joint = Counter((int(a), int(b)) for a, b in words)
            panels.append(
                dict(
                    level=level,
                    menu=[split, sharing, mode],
                    source_hits=private.sum(axis=0).tolist(),
                    counts=[[a, b, n] for (a, b), n in sorted(joint.items())],
                )
            )
    return dict(
        seed=SEED,
        starts=SHOTS,
        geometry=[board_layout(i) for i in range(3)],
        panels=panels,
    )


def oracle(n, split, sharing, mode):
    """Independent cut law from counting vertex assignments or independent edge endpoints."""
    p = Fraction(1, 2 if split == "half" else 4)
    if mode == "unchanged_loads":
        return {0: Fraction(1)}
    if sharing == "split":
        q = 2 * p * (1 - p)
        return {
            w: q ** w.bit_count() * (1 - q) ** (n - w.bit_count())
            for w in range(1 << n)
        }
    weights = Counter()
    for bits in product((0, 1), repeat=n):
        w = sum((bits[i] ^ bits[(i + 1) % n]) << i for i in range(n))
        if mode == "held_parked":
            w &= (1 << (n - 1)) - 1
        k = sum(bits)
        weights[w] += p**k * (1 - p) ** (n - k)
    return weights


def score(data):
    keys = {(i, *menu) for i in range(3) for menu in MENUS}
    if data["starts"] != SHOTS or data["seed"] != SEED:
        raise ValueError("Unexpected exposure or seed")
    tests = sum(
        len(MENUS) * ((1 << len(board_layout(i)["owners"])) + 2) for i in range(3)
    )
    radius = sqrt(log(2 * tests / 0.01) / (2 * SHOTS))
    result = {}
    for panel in data["panels"]:
        idx = panel["level"]
        menu = tuple(panel["menu"])
        key = (idx, *menu)
        if key not in keys:
            raise ValueError("Unexpected or duplicate panel")
        keys.remove(key)
        n = len(board_layout(idx)["owners"])
        expected = oracle(n, *menu)
        seen = {}
        words = Counter()
        for a, b, count in panel["counts"]:
            if (
                any(type(x) is not int for x in (a, b, count))
                or not 0 <= a < 1 << n
                or not 0 <= b < 1 << n
                or count <= 0
                or (a, b) in seen
            ):
                raise ValueError("Malformed endpoint receipt")
            seen[a, b] = count
            if a != b:
                raise ValueError("Endpoint meters disagree")
            words[a] += count
        if sum(words.values()) != SHOTS:
            raise ValueError("Missing starts")
        if set(words) - set(expected):
            raise ValueError("Observed a forbidden tape")
        even = sum(c for w, c in words.items() if not w.bit_count() % 2)
        for measured, p in [(words[w], expected.get(w, 0)) for w in range(1 << n)] + [
            (words[0], expected[0]),
            (even, sum(p for w, p in expected.items() if not w.bit_count() % 2)),
        ]:
            if abs(measured / SHOTS - float(p)) > radius:
                raise ValueError("Forecast rejected")
        result[key] = (words[0], even, len(words))
    if keys:
        raise ValueError("Missing panels")
    return result, radius


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = json.loads(json.dumps(acquire()))
    scores, width = score(data)
    path = HERE / "records.json"
    if args.write:
        with path.open("x") as f:
            json.dump(data, f, indent=1)
            f.write("\n")
    elif json.loads(path.read_text()) != data:
        raise ValueError("Archive differs from new acquisition")
    print("simultaneous radius", width)
    for key, v in scores.items():
        print(*key, *v)


if __name__ == "__main__":
    main()
