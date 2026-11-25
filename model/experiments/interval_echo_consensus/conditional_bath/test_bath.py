from fractions import Fraction
from itertools import product
import json
from math import comb

import numpy as np

from .measure import (
    ALL,
    CUT,
    HERE,
    N,
    RADIUS,
    WORDS,
    forward_kernel,
    grouped,
    law,
    physical_path,
    sweep_maps,
)


def test_native_kernel_equals_complete_fair_tape_and_preserves_inventory():
    for edges in (ALL, CUT):
        tape = sweep_maps(edges)
        empirical = np.array(
            [np.bincount(tape[:, i], minlength=70) / len(tape) for i in range(70)]
        )
        np.testing.assert_allclose(forward_kernel(edges), empirical, atol=4e-14, rtol=0)
        assert np.all(empirical.diagonal() > 0)
        for i, w in enumerate(WORDS):
            if edges == CUT:
                assert all(
                    (w & 3).bit_count() == (WORDS[j] & 3).bit_count()
                    for j in np.flatnonzero(empirical[i])
                )
        np.testing.assert_array_equal(empirical.sum(axis=0), 1)
        # The same deterministic all-left sweep repeated eight times sorts all words
        # on each active connected component; this proves almost sure finite capture.
        image = np.arange(len(WORDS))
        for _ in range(8):
            image = tape[0, image]
        assert len(set(image)) == (1 if edges == ALL else 3)


def test_correlated_reference_from_physical_words_and_every_fibre():
    counts = np.bincount(np.array(WORDS) & 15, minlength=16)
    r = law(True)
    p = law()
    assert list(Fraction(int(n), 70) for n in counts) == r
    assert sum(p) == sum(r) == 1 and min(p) > 0

    def mark(x):
        return int((x & 3).bit_count() == 1)

    assert sum(p[x] * mark(x) for x in range(16)) == Fraction(1, 4)
    for x in range(16):
        assert p[x] / r[x] == Fraction(7, 4) * Fraction(1, 4) ** mark(x)
        for port in (2, 3):
            y = x ^ (1 << port)
            one = x | (1 << port)
            others = (x & ~(1 << port)).bit_count()
            assert (
                r[one] / (r[x] + r[y])
                == p[one] / (p[x] + p[y])
                == Fraction(4 - others, 5)
            )
    for view in ("pair_totals", "bath_ports", "left_window", "right_window", "full"):
        rr, pp = grouped(r, view), grouped(p, view)
        assert all(pp[k] / rr[k] in (Fraction(7, 4), Fraction(7, 16)) for k in pp)


def test_coupled_and_forward_observer_banks_including_controls():
    d = json.loads((HERE / "results.json").read_bytes())
    p = np.array([float(x) for x in law()])
    r = np.array([float(x) for x in law(True)])
    eq = np.array([w >> 2 in (0, 3) for w in range(16)])
    for level, bank in enumerate(d["records"]):
        assert bank["path"] == physical_path(level)
        for x, expected in [("active", p), ("reference", r)]:
            observed = np.array(bank[x]) / N
            assert sum(bank[x]) == N and max(abs(observed - expected)) < RADIUS
        assert sum(bank["reference_depth"]) == N
        assert [sum(bank["active_depth"][str(k)]) for k in range(3)] == bank[
            "source_sectors"
        ]
        rate = np.dot(bank["active"], eq) / N
        assert (
            abs(rate - 0.45) < RADIUS
            and abs(rate - 13 / 24) > RADIUS
            and abs(rate - 0.5) > RADIUS
        )
    for bank in d["forward"]:
        frequencies = np.array(bank["counts"]) / N
        assert max(abs(frequencies - bank["expected"])) < RADIUS
        if bank["sweeps"] == 64:
            if bank["mode"] == "native":
                assert max(abs(frequencies - p)) < RADIUS
            else:
                assert max(abs(frequencies - p)) > RADIUS


def test_supplemental_fibres_and_observing_the_stopping_time():
    d = json.loads((HERE / "transitions.json").read_bytes())
    assert {(r["word"], r["old"], r["port"]) for r in d["rows"]} == set(
        product(range(8), range(2), (2, 3))
    )
    for row in d["rows"]:
        assert np.sum(row["depths"], axis=0).tolist() == row["outcome"]
        target = Fraction(4 - row["word"].bit_count(), 5)
        assert row["expected"] == str(target)
        assert sum(row["outcome"]) == N
        assert abs(Fraction(row["outcome"][1], N) - target) < RADIUS
    # For a single particle enumerate one sweep independently in physical coordinates.
    forced = []
    for bits in product((False, True), repeat=4):
        possible = set(range(5))
        for site, go_right in enumerate(bits):
            hit = bool(possible & {site, site + 1})
            possible -= {site, site + 1}
            if hit:
                possible.add(site + int(go_right))
        if len(possible) == 1:
            forced.append(next(iter(possible)))
    assert len(forced) == 2 and 0 not in forced
    assert d["first_sweep"]["conditional_one"] == "0"


def test_dyadic_fixed_wait_cannot_realize_exact_reference():
    assert law(True)[0] == Fraction(1, 70)
    assert law()[12] == Fraction(3, 20)
    assert all((2**bits) % 5 != 0 for bits in range(100))
    assert comb(8, 4) == 70


def test_composed_public_transactions_have_the_same_unique_fibre_attractor():
    from scipy.sparse.csgraph import connected_components

    r, p = law(True), law()
    moves = []
    for mode in (0, 2, 3):
        mat = np.zeros((16, 16), dtype=object)
        for x in range(16):
            if mode == 0:
                flip = (x & 12) | ((x & 1) << 1) | ((x & 2) >> 1)
                mat[x, x] = Fraction(1, 2)
                mat[x, flip] += Fraction(1, 2)
            else:
                held = x & ~(1 << mode)
                f = Fraction(4 - held.bit_count(), 5)
                mat[x, held] = 1 - f
                mat[x, held | (1 << mode)] = f
        for stationary in (p, r):
            flows = np.array(stationary, dtype=object)[:, None] * mat
            np.testing.assert_array_equal(flows, flows.T)
        moves.append(mat)
    k = moves[0] @ moves[1] @ moves[2]
    assert np.all(k.diagonal() > 0)
    n, labels = connected_components(
        np.asarray(k, float), directed=True, connection="strong"
    )
    assert n == 3
    for x, y in product(range(16), repeat=2):
        assert (labels[x] == labels[y]) == ((x & 3).bit_count() == (y & 3).bit_count())
