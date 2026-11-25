"""Integer tests of the spectator detunings used to isolate individual material transitions."""

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from .circuit import K, local_h

HERE = Path(__file__).resolve().parent
PAIRS = tuple(combinations(range(3), 2))
BASE = 256


def diagonal(a, b):
    """Twice the diagonal of one blank native contact."""
    return int(2 * K) * (a - b) ** 2 + (a != b)


def symbols(a, b):
    return [diagonal(b, c) - diagonal(a, c) for c in range(3)]


def detuning(a, b, spectators):
    weights = [BASE**i for i in range(len(spectators))]
    offset = 131 * sum(weights) + 1
    digits = symbols(a, b)
    return (b * b - a * a) * offset + sum(
        w * digits[c] for w, c in zip(weights, spectators, strict=True)
    )


def certify():
    local = np.zeros((9, 9), int)
    for a, b in product(range(3), repeat=2):
        local[3 * a + b, 3 * a + b] = diagonal(a, b)
        if a != b:
            local[3 * b + a, 3 * a + b] = -1
    np.testing.assert_allclose(2 * local_h(K), local, rtol=0, atol=1e-12)
    sq = [0] + [(b * b - a * a) ** 2 for a, b in PAIRS]
    interpolation = {}
    for pair, key in zip(PAIRS, sq[1:], strict=True):

        def value(t):
            v = Fraction(1)
            for other in sq:
                if key != other:
                    v *= Fraction(t - other, key - other)
            return v

        if [value(t) for t in sq] != [int(t == key) for t in sq]:
            raise ValueError("Species projection failed")
        interpolation[str(pair)] = {
            "spectator_digits": symbols(*pair),
            "double_commutator_gap": key,
        }
    counts = []
    for size in range(2, 9):
        gaps_per_pair = []
        for a, b in PAIRS:
            ws = list(product(range(3), repeat=size - 2))
            ds = {detuning(a, b, w) for w in ws}
            if len(ds) != len(ws) or min(ds) <= 0:
                raise ValueError("Spectators collided or orientation ambiguous")
            gaps_per_pair.append(len(ds))
        counts.append([size, gaps_per_pair])
    return {
        "native_twice_blank": local.tolist(),
        "species": interpolation,
        "base": BASE,
        "checked_distinct_positive_gaps_per_species": counts,
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()
    result = certify()
    path = HERE / "certificate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_bytes()):
        raise ValueError("Changed inventory controllability certificate")
    print(
        "Three native exchange channels separated with every spectator placement through eight data ports."
    )
