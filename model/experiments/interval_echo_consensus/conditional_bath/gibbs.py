"""Fixed-appointment Gibbs projection using measured native partial exchanges."""

import argparse
import json
from functools import cache
from itertools import combinations
from math import asin, pi, sin, sqrt

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.candidate import KAPPA, native_law

from .measure import HERE, INDEX, N, WORDS, physical_path

PACKED = np.asarray(WORDS, dtype=np.int64)
MOMENT = np.array([(w & 3).bit_count() % 2 for w in WORDS])
PAIR_SLOTS = tuple(combinations(range(8), 2))
EXPOSURES = (0, pi / 4, pi / 3, pi / 2, 2 * pi / 3, 3 * pi / 4, pi)


def adjacent(word, edge):
    return word ^ (3 << edge) if ((word >> edge) & 3) in (1, 2) else word


def shuffle(start, length):
    return [
        (start + i, i + 1, i + 2) for j in range(1, length) for i in reversed(range(j))
    ]


def transposition(left, right):
    return (*range(left, right), *reversed(range(left, right - 1)))


def matching_flip(word):
    occupied = [i for i in range(8) if word & (1 << i)]
    empty = [i for i in range(8) if not word & (1 << i)]
    pairs = {tuple(sorted(pair)) for pair in zip(occupied, empty, strict=True)}
    for pair in PAIR_SLOTS:
        if pair in pairs:
            for edge in transposition(*pair):
                word = adjacent(word, edge)
    return word


FLIP = np.array([INDEX[matching_flip(w)] for w in WORDS])
SWAPS = np.array([[INDEX[adjacent(w, e)] for w in WORDS] for e in range(7)])


@cache
def native_exchange(tau, disabled=False):
    exchange, closure = native_law.local_term_parts()
    h = KAPPA * closure if disabled else exchange + KAPPA * closure
    gate = np.abs(expm(-1j * tau * h)) ** 2
    probs = []
    for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
        col = native_law.SRC.index((a, b))
        swap = native_law.SRC.index((b, a))
        allowed = {col, swap}
        if np.max(np.delete(gate[:, col], sorted(allowed))) > 2e-13:
            raise ValueError("Native contact escaped its binary rail")
        if a == b:
            np.testing.assert_allclose(gate[col, col], 1.0, atol=2e-13, rtol=0)
        else:
            p = gate[swap, col]
            np.testing.assert_allclose(gate[col, col] + p, 1.0, atol=2e-13, rtol=0)
            probs.append(p)
    np.testing.assert_allclose(probs[0], probs[1], atol=2e-13, rtol=0)
    return float(probs[0])


@cache
def swap_rate(num, den, disabled=False):
    return native_exchange(2 * asin(sqrt(num / den)), disabled)


def exact_and_native(schedule, complement=False, disabled=False, fair_only=False):
    n = len(WORDS)
    integer = np.eye(n, dtype=int).astype(object)
    real = np.eye(n)
    den = 1
    if complement:
        integer += integer[:, FLIP]
        den = 2
        p = native_exchange(pi / 2)
        real = (1 - p) * real + p * real[:, FLIP]
    for edge, a, b in schedule:
        integer = (b - a) * integer + a * integer[:, SWAPS[edge]]
        den *= b
        p = swap_rate(1, 2) if fair_only else swap_rate(a, b, disabled)
        real = (1 - p) * real + p * real[:, SWAPS[edge]]
    return integer, den, real


def datum(label, word):
    if label == "reference":
        return 0
    if label == "source":
        return (word & 3).bit_count() % 2
    if label == "cut":
        return (word & 3).bit_count()
    port = int(label[-1])
    return word & (15 ^ (1 << port))


def certify(label, schedule, complement=False):
    k, den, real = exact_and_native(schedule, complement)
    labels = np.asarray([datum(label, w) for w in WORDS])
    same = labels[:, None] == labels
    fibre_sizes = same.sum(axis=1)
    np.testing.assert_array_equal(k * fibre_sizes[:, None], same.astype(int) * den)
    target = same / fibre_sizes[:, None]
    np.testing.assert_allclose(real, target, atol=2e-13, rtol=0)
    return {
        "complete_rows": len(k),
        "entries": int(k.size),
        "fibre_sizes": sorted(set(map(int, fibre_sizes))),
        "native_max_abs_error": float(np.max(abs(real - target))),
        "shuffle_contacts": len(schedule),
    }


SOURCE = shuffle(0, 2) + shuffle(2, 6)
REFERENCE = shuffle(0, 8)
LOCAL = shuffle(3, 5)
ROUTE = [(2, 1, 1)]


def run_shots(rng, start, schedule, complement=False):
    indices = np.array(start, dtype=int, copy=True)
    if complement:
        toss = rng.random(len(indices)) < native_exchange(pi / 2)
        indices[toss] = FLIP[indices[toss]]
    for e, a, b in schedule:
        hit = rng.random(len(indices)) < swap_rate(a, b)
        indices[hit] = SWAPS[e, indices[hit]]
    return indices


def run():
    np.testing.assert_array_equal(PACKED[FLIP], 255 ^ PACKED)
    np.testing.assert_allclose(native_exchange(pi), 1.0, atol=2e-13, rtol=0)
    checks = {
        name: certify(name, schedule, orbit)
        for name, schedule, orbit in (
            ("reference", REFERENCE, False),
            ("source", SOURCE, True),
            ("cut", SOURCE, False),
            ("port3", LOCAL, False),
            ("port2", ROUTE + LOCAL + ROUTE, False),
        )
    }
    _, _, half = exact_and_native(LOCAL, fair_only=True)
    _, _, full = exact_and_native(LOCAL)
    _, _, disabled = exact_and_native(LOCAL, disabled=True)
    controls = {
        "half_only_max_fibre_error": float(np.max(abs(half - full))),
        "disabled_max_identity_error": float(
            np.max(abs(disabled - np.eye(len(WORDS))))
        ),
        "disabled_max_fibre_error": float(np.max(abs(disabled - full))),
    }
    banks = []
    for level in range(3):
        rng = np.random.default_rng(51484 + level)
        for i, tau in enumerate(EXPOSURES):
            selected = native_exchange(tau)
            np.testing.assert_allclose(selected, sin(tau / 2) ** 2, atol=2e-13, rtol=0)
            cal = int((rng.random(N) < selected).sum())
            paired = []
            for even, odd in [(60, 29), (195, 226)]:
                flags = rng.random(N) < selected
                starts = np.where(flags, INDEX[odd], INDEX[even])
                np.testing.assert_array_equal(MOMENT[starts], flags)
                indices = run_shots(rng, starts, SOURCE, complement=True)
                np.testing.assert_array_equal(MOMENT[indices], flags)
                counts = np.bincount(PACKED[indices] & 15, minlength=16)
                paired.append(
                    {
                        "even": even,
                        "odd": odd,
                        "odd_starts": int(flags.sum()),
                        "word_counts": counts.tolist(),
                    }
                )
            banks.append(
                dict(level=level, tau_index=i, calibration_ones=cal, runs=paired)
            )
    return {
        "shots_per_bank": N,
        "exposures_in_pi": ["0", "1/4", "1/3", "1/2", "2/3", "3/4", "1"],
        "paths": [physical_path(level) for level in range(3)],
        "fixed_matching_slots": sum(len(transposition(*pair)) for pair in PAIR_SLOTS),
        "checks": checks,
        "controls": controls,
        "banks": banks,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    receipt = run()
    dest = HERE / "gibbs.json"
    data = json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        dest.write_text(data)
    elif data != dest.read_text():
        raise ValueError("Replay mismatch")
    print("Five exact all-input kernels, three controls and 42 forward source banks")
