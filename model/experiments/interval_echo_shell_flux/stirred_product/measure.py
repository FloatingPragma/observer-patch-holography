"""Renewal currents on three independently ranged Interval Echo addresses."""

import argparse
from collections import deque
from fractions import Fraction
from functools import cache
from itertools import product
import json
from math import pi
from pathlib import Path

import numpy as np
from sympy import Matrix, Rational, ones

from ..product_address import measure as address
from model.experiments.interval_echo_consensus.native_path_reference.occupation_action.measure import (
    native_kernel,
)

HERE = Path(__file__).resolve().parent
CYCLES = 2**18
N = 12**3
SHAPE = (12, 12, 12)
IDLE = (0, 0)


def index(pos):
    return int(np.ravel_multi_index(pos, SHAPE))


@cache
def survey():
    address.calibrate()
    b, axes, _ = address.layout()
    ranges, _ = address.timing()
    census = address.measured_cuts(ranges, None)["outward"]
    corners = [
        index(address.release_corner(c)) for c in product(address.CENTRES, repeat=3)
    ]
    if len(set(corners)) != 8:
        raise ValueError("Source release failed")
    if native_kernel((pi / 2, pi / 2)) != ((Fraction(1, 2),) * 2,) * 2:
        raise ValueError("Material tag is biased")
    layer = np.array([ranges[pos] for pos in product(range(12), repeat=3)])
    offers = np.full((N, 6), -1, dtype=int)
    for p, pos in enumerate(product(range(12), repeat=3)):
        for axis, path in enumerate(axes):
            for sign in (-1, 1):
                j = pos[axis] + sign
                if 0 <= j < 12:
                    a, bp = path[pos[axis]], path[j]
                    if bp not in b.neighbours[a]:
                        raise ValueError("Off-board contact")
                    q = list(pos)
                    q[axis] = path.index(address.pulse((a, bp), a))
                    offers[p, 2 * axis + int(sign == 1)] = index(tuple(q))
    return layer, offers, corners, census


@cache
def shell(n, bias):
    layer, offers, _, _ = survey()
    sites = np.flatnonzero(layer == n)
    extras = offers[sites, :2] if bias else np.repeat(sites[:, None], 2, axis=1)
    dest = np.column_stack((offers[sites], extras))
    if dest.min() < 0:
        raise ValueError("Active source touches external edge of address tape")
    return sites, dest


@cache
def tree(n):
    layer, offers, _, _ = survey()
    root = int(np.flatnonzero(layer == n)[0])
    parent = {root: root}
    order = deque([root])
    while order:
        p = order.popleft()
        for q in sorted(offers[p]):
            if q >= 0 and layer[q] == n and q not in parent:
                parent[int(q)] = p
                order.append(int(q))
    if len(parent) != sum(layer == n):
        raise ValueError("Disconnected stirred shell")
    return parent


def distribute(rng, k, bins, counts):
    counts[bins] = rng.multinomial(k, np.full(len(bins), 1 / len(bins)))


def clocks(depth, bias):
    """Independent mean and variance for the finite first absorption time."""
    layer, _, _, _ = survey()
    p = [[Rational(0) for _ in range(depth)] for _ in range(depth)]
    for n in range(1, depth + 1):
        _, dest = shell(n, bias)
        for target in range(1, depth + 1):
            p[n - 1][target - 1] = Rational(
                int((layer[dest] == target).sum()), dest.size
            )
    f = (Matrix.eye(depth) - Matrix(p)).inv()
    mean = f * ones(depth, 1)
    second = f * (2 * mean - ones(depth, 1))
    return str(2 + mean[0]), str(second[0] - mean[0] ** 2)


def complete(rng, depth, bias, cycles):
    layer, _, corners, _ = survey()
    source = rng.multinomial(cycles, [1 / 8] * 8)
    old = np.zeros(N, dtype=np.int64)
    old[corners] = source
    new = np.zeros_like(old)
    tags = np.zeros((N, 8), dtype=np.int64)
    inward = np.zeros_like(old)
    outward = np.zeros_like(old)
    up = cycles
    refused = {}
    for n in range(depth, 0, -1):
        sites, dest = shell(n, bias)
        move = layer[dest].ravel() - n
        success, failure = np.flatnonzero(move > 0), np.flatnonzero(move <= 0)
        failures = rng.negative_binomial(up, len(success) / len(move))
        counts = np.zeros(len(move), dtype=np.int64)
        distribute(rng, up, success, counts)
        distribute(rng, failures, failure, counts)
        tags[sites] = counts.reshape(dest.shape)
        np.add.at(old, dest.ravel(), counts)
        new[sites] = tags[sites].sum(axis=1)
        up = cycles + int(counts[move < 0].sum())
        bits = (len(sites) - 1).bit_length()
        refused[str(n)] = int(
            rng.negative_binomial(int(counts.sum()), len(sites) / 2**bits)
        )
    # Old arrival word and independently drawn target are both transported
    # through the tree root. Every traversed leg is observable and counted.
    for n in range(1, depth + 1):
        parents = tree(n)
        for p, parent in reversed(parents.items()):
            if p != parent:
                inward[p] += old[p]
                outward[p] += new[p]
                inward[parent] += inward[p]
                outward[parent] += outward[p]
        root = next(iter(parents))
        if inward[root] + old[root] != outward[root] + new[root]:
            raise ValueError("Missing visit or stir transfer")
        inward[root] = outward[root] = 0
    if sum(old[layer == depth + 1]) != cycles:
        raise ValueError("Missing absorption")
    appointments = int(tags.sum()) + 2 * cycles
    return (
        dict(
            depth=depth,
            biased=bias,
            cycles=cycles,
            source_corners=source.tolist(),
            collectors=old[layer == depth + 1].tolist(),
            appointments=appointments,
            direction_idles=int(tags[:, 6:].sum()) if not bias else 0,
            tangential_pulses=int(inward.sum() + outward.sum()),
            target_retries=refused,
        ),
        (tags, inward, outward),
    )


def sample_literal(rng, depth, episodes):
    """Uncompressed check: actual offered word/bit choices, then local route legs."""
    layer, _, corners, _ = survey()
    crossing = np.zeros(depth, dtype=int)
    clock = 0
    hops = 0
    refused = 0
    cut = np.zeros(8, dtype=int)
    full_steps = 0
    for _ in range(episodes):
        p = corners[int(rng.integers(8))]
        clock += 1
        while layer[p] <= depth:
            n = int(layer[p])
            sites, dest = shell(n, False)
            bits = (len(sites) - 1).bit_length()
            while True:
                x = int(rng.integers(2**bits))
                if x < len(sites):
                    break
                refused += 1
            y = int(sites[x])
            parents = tree(n)
            a = p
            while parents[a] != a:
                a = parents[a]
                hops += 1
            a = y
            while parents[a] != a:
                a = parents[a]
                hops += 1
            tag = int(rng.integers(8))
            q = int(dest[x, tag])
            if tag < 6:
                full_steps += 1
            cut[tag] += 1
            if layer[q] != n:
                crossing[min(layer[q], n) - 1] += int(layer[q] - n)
            p = q
            clock += 1
        clock += 1
    return dict(
        episodes=episodes,
        appointments=clock,
        tag_counts=cut.tolist(),
        net_by_shell=crossing.tolist(),
        retries=refused,
        stir_pulses=hops,
        move_pulses=full_steps,
    )


def parked(rng):
    _, axes, flag = address.layout()
    clicks = [
        address.pulse((path[5], path[6]), path[5], True) == path[6] for path in axes
    ]
    fired = address.pulse(flag, flag[0], True) == flag[1]
    if any(clicks) or fired:
        raise ValueError("Changed kinetic transfers unexpectedly")
    # Independent challenge starts with prepared triples after the source.
    # A failed stirring relay is read and that appointment ends on the old
    # word. Even an accepted same-word target cannot leave on its step.
    counts = rng.multinomial(10 * CYCLES, np.full(8, 1 / 8))
    return dict(
        source_trials=CYCLES,
        flag_fired=int(fired) * CYCLES,
        preloaded_trials=CYCLES,
        attempted_ticks=10 * CYCLES,
        proposed_direction_tags=counts.tolist(),
        successful_coordinate_changes=0,
        drains=0,
    )


def run():
    rng = np.random.default_rng(1396123)
    survey_result = survey()
    banks = {}
    archive = {}
    for depth, bias in ((3, False), (4, False), (5, False), (5, True)):
        desc = f"{depth}_{'x' if bias else 'regular'}"
        for meter in ("source", "edge"):
            bank, arrays = complete(rng, depth, bias, CYCLES)
            k = f"{desc}_{meter}"
            banks[k] = bank
            for name, x in zip(("tags", "to_root", "from_root"), arrays, strict=True):
                archive[f"{k}_{name}"] = x.astype(np.int32)
    return dict(
        cut_survey=survey_result[3],
        banks=banks,
        renewal_moments={
            f"{d}_{int(b)}": clocks(d, b)
            for d, b in ((3, False), (4, False), (5, False), (5, True))
        },
        literal=sample_literal(rng, 5, 512),
        parked=parked(rng),
    ), archive


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    record, bins = run()
    encoded = json.dumps(record, separators=(",", ":"), allow_nan=False) + "\n"
    if args.write:
        (HERE / "counts.json").write_text(encoded)
        np.savez_compressed(HERE / "edges.npz", **bins)
    else:
        if (HERE / "counts.json").read_text() != encoded:
            raise SystemExit("Source, drain or clock census changed")
        with np.load(HERE / "edges.npz", allow_pickle=False) as expected:
            if set(expected.files) != set(bins):
                raise SystemExit("Wrong aperture enrollment")
            for k, value in bins.items():
                np.testing.assert_array_equal(expected[k], value)
    for key, b in record["banks"].items():
        print(key, b["cycles"] / b["appointments"], b["tangential_pulses"])
