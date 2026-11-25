"""Read the bottom energy from an initially maximally mixed three-dimer ring."""

import argparse
from functools import lru_cache
from itertools import product
import json

import numpy as np
from scipy.stats import binom

from model.executable_model.carrier import federation
from model.candidate_laws.interval_echo import candidate
from .geometry import layouts
from .measure import HERE

BINS, PERIOD, REPEATS, LAUNCHES = 2**18, 1024, 25, 256
GRID = 2**25


def native_matrix(layout, charges, cut=-1, kinetic=True):
    """Project the microscopic qutrit seam columns directly at these physical ports."""
    board = federation(0)
    pairs = layout["pairs"]
    n = len(pairs)
    raw, contrast = candidate.native_law.local_term_parts()
    local = (raw if kinetic else 0) + candidate.KAPPA * contrast
    index = candidate.native_law.SRC.index
    out = np.zeros((2**n, 2**n))
    for word in range(2**n):
        occ = {pin: 1 + q for pin, q in zip(layout["pins"], charges, strict=True)}
        for k, (a, b) in enumerate(pairs):
            occ[a] = 2 * ((word >> k) & 1)
            occ[b] = 2 - occ[a]
        for p, v in occ.items():
            out[word, word] += candidate.KAPPA * (6 - board.degree[p]) * v**2
        for k, (a, b) in enumerate(pairs):
            old = index((occ[a], occ[b]))
            out[word, word] += local[old, old]
            out[word ^ (1 << k), word] += local[index((occ[b], occ[a])), old]
        links = [
            (b, pairs[(k + 1) % n][0]) for k, (_, b) in enumerate(pairs) if k != cut
        ]
        links += [
            (p, pairs[site][0])
            for p, site in zip(layout["pins"], layout["marks"], strict=True)
        ]
        for a, b in links:
            i = index((occ[a], occ[b]))
            out[word, word] += local[i, i]
    return out


@lru_cache(maxsize=8)
def clock_distribution(tick):
    x = tick / GRID * BINS / PERIOD
    bins = np.arange(BINS)
    # Rotate each spectrum into the public window [-128,896).
    d = (x - bins + BINS / 2) % BINS - BINS / 2
    weights = (np.sinc(d) / np.sinc(d / BINS)) ** 2
    if abs(weights.sum() - 1) > 1e-9:
        raise ValueError("Clock probability loss")
    cdf = np.cumsum(weights / weights.sum())
    cdf[-1] = 1
    return cdf


def acquire():
    rng = np.random.default_rng(520118)
    layout = layouts()[0]
    panels = []
    for q in product((-1, 0, 1), repeat=3):
        h = native_matrix(layout, q)
        energies = np.linalg.eigvalsh(h)
        ticks = np.rint(energies * GRID).astype(int)
        if min(energies) < 0 or max(energies) > 768:
            raise ValueError("Populated line outside the interior of the clock window")
        # Identity/8 is the equally weighted mixture in every eigenbasis.
        # Repeated Fourier energy instruments are diagonal in this same basis.
        branch = rng.integers(0, 8, LAUNCHES)
        u = rng.random((LAUNCHES, REPEATS))
        counts = np.empty(u.shape, dtype=np.int32)
        for k, t in enumerate(ticks):
            mask = branch == k
            counts[mask] = np.searchsorted(clock_distribution(t), u[mask])
        panels.append(counts)
    return np.stack(panels)


def ground_reads(raw):
    decoded = np.where(raw >= 7 * BINS // 8, raw - BINS, raw) * (PERIOD / BINS)
    return np.median(decoded, axis=-1).min(axis=-1).reshape((3, 3, 3))


def failure_bound():
    # Two one-sided tails outside 8 bins each have mass at most
    # (1/4) sum_{j>=8} 1/j² <= 1/28.
    return 27 * (
        (7 / 8) ** LAUNCHES + LAUNCHES * binom.sf(REPEATS // 2, REPEATS, 1 / 14)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    raw = acquire()
    path = HERE / "clock.npz"
    if args.write:
        np.savez_compressed(path, exits=raw)
    else:
        with np.load(path, allow_pickle=False) as f:
            np.testing.assert_array_equal(f["exits"], raw)
    e = ground_reads(raw)
    out = [
        (s[0, 0] + s[2, 2] - s[0, 2] - s[2, 0]) / 4
        for s in (e[:, :, 1], e[:, 1, :], e[1, :, :])
    ]
    print(
        json.dumps(
            dict(pair_readings=out, whole_failure_bound=failure_bound()), indent=2
        )
    )
