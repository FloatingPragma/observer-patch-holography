"""Certificates of which parked ring configurations the two material probes see."""

import json
from functools import cache
from itertools import product
from pathlib import Path

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation
from ..ring import geometry

BOARDS = ((0, 14), (1, 14), (2, 14), (0, 16))
HERE = Path(__file__).resolve().parent


@cache
def apertures(level, size):
    board = federation(level)
    ring = tuple(geometry(level, size))
    on_ring = set(ring)

    def probes(offsets):
        out = []
        for j in offsets:
            ends = {ring[j], ring[(j + 1) % size]}
            for a in sorted(
                set.intersection(*(set(board.neighbours[x]) for x in ends)) - on_ring
            ):
                if set(board.neighbours[a]) & on_ring != ends:
                    continue
                for b in sorted(set(board.neighbours[a]) - on_ring):
                    if a // 12 == b // 12 and not on_ring.intersection(
                        board.neighbours[b]
                    ):
                        out.append((j, a, b))
        return sorted(out)

    for left, right in product(probes(range(size - 5, size - 2)), probes(range(2, 5))):
        if (left[0] - right[0]) % 2 == 0 or left[1] // 12 == right[1] // 12:
            continue
        if any(y in board.neighbours[x] for x in left[1:] for y in right[1:]):
            continue
        return ring, (left, right)
    raise ValueError("No circuit meets the registered local probe criteria")


def wall_pairs(word, n):
    return tuple(
        j for j in range(n) if ((word >> j) ^ (word >> ((j + 1) % n))) & 1 == 0
    )


def source_census(size, offsets, sign):
    """Enumerate every same-sign wall record and every clean local swap witness."""
    tips = {i for j in offsets for i in (j, (j + 1) % size)}
    assert len(tips) == 4 and sign in (-1, 1)
    fixed = sum(1 << x for x in tips)
    target = fixed if sign > 0 else 0
    same = []
    opposite = []
    originals = []
    for word in range(1 << size):
        if word & fixed != target:
            continue
        delta = word.bit_count() - size // 2
        if delta == sign:
            same.append(word)
            if set(wall_pairs(word, size)) == set(offsets):
                originals.append(word)
        if delta == -sign:
            opposite.append(word)
    witnesses = []
    for word in originals:
        moves = []
        for j in range(size):
            k = (j + 1) % size
            mask = (1 << j) | (1 << k)
            if {j, k} & tips or (word & mask) in (0, mask):
                continue
            moved = word ^ mask
            if len(wall_pairs(moved, size)) == 4:
                assert moved in same
                moves.append([j, k, moved])
        witnesses.append({"source_word": word, "swap_bonds_and_resulting_words": moves})
    return dict(
        same_inventory_words=len(same),
        opposite_inventory_words=len(opposite),
        two_wall_words=len(originals),
        opposite_inventory_example=opposite[0],
        witnesses=witnesses,
    )


@cache
def blank_terms():
    """Obtain local matrix elements, including leakage checks, from the full 99-state law."""
    x, c = candidate.native_law.local_term_parts()
    for h in (x, c):
        np.testing.assert_allclose(h[9:, :9], 0, atol=4e-15)
    return x[:9, :9], c[:9, :9]


def rail_matrix(board, sensor, installed, sign, kappa=8):
    """Sum the physical contact and private-zero terms in each probe branch."""
    a, b = sensor
    if not {a, b}.isdisjoint(installed):
        raise ValueError("Probe rails already occupied by source")
    stock = 1 - sign
    paths = [(1, stock), (stock, 1)]
    src = candidate.native_law.SRC
    exchange, contrast = blank_terms()
    idx = [src.index(w) for w in paths]
    h = (exchange + kappa * contrast)[np.ix_(idx, idx)].copy()
    for branch, terms in enumerate(paths):
        for p, other, load in zip(sensor, (b, a), terms, strict=True):
            neighbours = [v for v in board.neighbours[p] if v != other]
            exterior = [installed.get(v, 0) for v in neighbours] + [0] * (
                6 - board.degree[p]
            )
            for static in exterior:
                i = src.index((load, static))
                h[branch, branch] += kappa * contrast[i, i]
    return h - np.trace(h) * np.eye(2) / 2


def check_operators(level, n, ring, settings, source, amount, sign):
    words = source["witnesses"]
    board = federation(level)
    checks = 0

    def instrument(word):
        loads = {p: amount * ((word >> j) & 1) for j, p in enumerate(ring)}
        return np.array(
            [rail_matrix(board, rails[1:], loads, sign) for rails in settings]
        )

    for row in words:
        expected = instrument(row["source_word"])
        for w in [
            source["opposite_inventory_example"],
            *(s[2] for s in row["swap_bonds_and_resulting_words"]),
        ]:
            np.testing.assert_allclose(expected, instrument(w), atol=2e-13, rtol=0)
            checks += 2
    return checks


def empirical(level, n, offsets, amount, sign, label):
    """Count the earlier sharp word records; excitation times and words are public metadata."""
    tips = {i for j in offsets for i in (j, (j + 1) % n)}
    fixed = sum(1 << x for x in tips)
    target = fixed if sign > 0 else 0
    with np.load(HERE.parent / "words.npz") as data:
        words = data[f"words_{n}"].astype(int)
        if sign == -1:
            words ^= (1 << n) - 1
        records = data[f"counts_{level}_{n}_{amount}_{sign}_{label}"][4].astype(
            np.int64
        )
    local = (words & fixed) == target
    pure = np.array([set(wall_pairs(w, n)) == set(offsets) for w in words])
    return dict(
        starts=int(sum(records)),
        local_matches=int(sum(records[local])),
        two_wall_matches=int(sum(records[local & pure])),
    )


def main():
    runs = []
    for level, n in BOARDS:
        ring, sensors = apertures(level, n)
        offsets = tuple(s[0] for s in sensors)
        for sign in (-1, 1):
            census = source_census(n, offsets, sign)
            checks = sum(
                check_operators(level, n, ring, sensors, census, amount, probe_sign)
                for amount, probe_sign in product((1, 2), (-1, 1))
            )
            prior = {
                f"{amount}_{law}": empirical(level, n, offsets, amount, sign, law)
                for amount, law in product((1, 2), ("native", "contrast_off"))
            }
            runs.append(
                dict(
                    level=level,
                    size=n,
                    sign=sign,
                    sensors=sensors,
                    **census,
                    same_hamiltonian_checks=checks,
                    prior_release_counts=prior,
                )
            )
    return runs


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
