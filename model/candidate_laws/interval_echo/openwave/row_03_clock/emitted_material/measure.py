"""Isolate exact Bohr frequencies for the existing neutral-emission materials."""

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import product
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

from ...row_09_neutral.extended_decay.material import route
from ...row_09_neutral.extended_decay.measure import rebuild
from ...row_09_neutral.passive_decay.material import edge, hardware
from ...row_09_neutral.passive_decay.measure import physical

HERE = Path(__file__).resolve().parent
EPSILON = sp.Rational(1, 10**14)


def twice_generator(diagonal):
    """Integer 2H on one orientation bit per native blank dimer."""
    h = np.diag(2 * np.asarray(diagonal, dtype=int))
    words = np.arange(len(h))
    for bit in range((len(h) - 1).bit_length()):
        h[words, words ^ (1 << bit)] = -1
    return h


@lru_cache(maxsize=None)
def spectrum(diagonal):
    polynomial = sp.Matrix(twice_generator(diagonal)).charpoly().as_poly()
    intervals = polynomial.intervals(eps=EPSILON)
    assert sum(m for _, m in intervals) == len(diagonal)
    assert intervals[0][1] == 1
    levels = [(a / 2, b / 2) for (a, b), _ in intervals]
    # Multiplicity can add a zero-frequency term and repeats no positive difference.
    low, first = levels[:2]
    rest = (first[0] - low[1], first[1] - low[0])
    assert rest[0] > 0
    bounds, overlaps = [], []
    for i, (lo, hi) in enumerate(levels[1:], 1):
        for j, (alo, ahi) in enumerate(levels[1:i], 1):
            delta = (lo - ahi - rest[1], hi - alo - rest[0])
            if delta[0] <= 0 <= delta[1]:
                overlaps.append([j, i])
            bounds.append(max(sp.Integer(0), delta[0], -delta[1]))
    # Keep a rigorous downward decimal truncation in exchange units.
    lower = str(sp.floor(min(bounds) * 10**12) / sp.Integer(10**12))
    return {
        "rest_interval": [str(x) for x in rest],
        "excited_rest_candidates": overlaps,
        "detuning_lower_rational": lower,
        "polynomial_2H": [int(x) for x in polynomial.all_coeffs()],
    }


def acquire():
    board, dimers, head_template = hardware()
    addresses = sorted(head_template)
    grouped = defaultdict(Counter)
    for signs in product((32, 33), repeat=len(addresses)):
        heads = dict(zip(addresses, signs, strict=True))
        native = physical(board, dimers, heads)
        # A=B=0 eliminates the core-receiver ac term; compare with a fixed source.
        diag = tuple(map(int, native[:4] - native[0]))
        left = (
            heads[edge(dimers[0][0], dimers[2][0])]
            != heads[edge(dimers[0][0], dimers[2][1])]
        )
        right = (
            heads[edge(dimers[2][0], dimers[3][0])]
            != heads[edge(dimers[2][0], dimers[3][1])]
        )
        mixed = diag[3] - diag[1] - diag[2] + diag[0]
        assert bool(mixed) == right
        grouped[diag]["intact" if left and right else "other"] += 1
    banks = []
    for diag, counts in sorted(grouped.items()):
        result = spectrum(diag)
        factorizes = diag[3] == diag[1] + diag[2] - diag[0]
        if factorizes:
            assert result["excited_rest_candidates"]
        else:
            assert not result["excited_rest_candidates"], (diag, result)
        banks.append(
            {
                "diagonal": diag,
                "both_links_conduct": counts["intact"],
                "other": counts["other"],
                "receiver_factorizes": factorizes,
                **result,
            }
        )
    board, dimers, heads, _, _ = route((3,), True)
    all_contacts = rebuild(board, dimers, heads, False, True)
    diag = tuple(map(int, np.diag(all_contacts)))
    np.testing.assert_array_equal(twice_generator(diag), 2 * all_contacts)
    crossing = spectrum(diag)
    return dict(histories=banks, crossing={"dimers": dimers, **crossing})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(acquire(), indent=2) + "\n"
    if args.output:
        args.output.write_text(result)
    else:
        assert (HERE / "frequencies.json").read_text() == result
        print(
            "Verified all sharp material stocks and the complete joined cross-owner spectrum."
        )
