from fractions import Fraction
from itertools import product
import json
from math import comb, sin

import numpy as np

from .gibbs import (
    EXPOSURES,
    FLIP,
    LOCAL,
    PACKED,
    REFERENCE,
    ROUTE,
    SOURCE,
    datum,
    exact_and_native,
    shuffle,
)
from .measure import HERE, N, RADIUS, WORDS, physical_path


def test_one_appointment_performs_the_complete_conditional_on_every_input():
    for label, schedule, flip in (
        ("reference", REFERENCE, False),
        ("source", SOURCE, True),
        ("cut", SOURCE, False),
        ("port2", ROUTE + LOCAL + ROUTE, False),
        ("port3", LOCAL, False),
    ):
        raw, den, measured = exact_and_native(schedule, flip)
        labels = [datum(label, w) for w in WORDS]
        for i, w in enumerate(WORDS):
            allowed = [z for z in WORDS if datum(label, z) == labels[i]]
            expected = [Fraction(int(z in allowed), len(allowed)) for z in WORDS]
            assert [Fraction(v, den) for v in raw[i]] == expected
            np.testing.assert_allclose(
                measured[i], list(map(float, expected)), atol=2e-13, rtol=0
            )
    np.testing.assert_array_equal(PACKED[FLIP], 255 ^ PACKED)
    assert len(shuffle(0, 8)) == 28


def test_projection_fibres_and_one_state_family_on_overlapping_reads():
    micro = [w & 15 for w in WORDS]
    r = [Fraction(micro.count(x), 70) for x in range(16)]
    for s in (
        Fraction(0),
        Fraction(1, 8),
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(4, 7),
        Fraction(13, 14),
        Fraction(1),
    ):
        p = [
            Fraction(micro.count(x)) * (s / 40 if datum("source", x) else (1 - s) / 30)
            for x in range(16)
        ]
        assert sum(p) == 1 and sum(p[x] * datum("source", x) for x in range(16)) == s
        if 0 < s < 1:
            assert all(
                p[x] / r[x]
                == Fraction(7, 3)
                * (1 - s)
                * (3 * s / (4 * (1 - s))) ** datum("source", x)
                for x in range(16)
            )
        for x, port in product(range(16), (2, 3)):
            held = x & ~(1 << port)
            one = held | (1 << port)
            expected = Fraction(4 - held.bit_count(), 5)
            assert r[one] == expected * (r[held] + r[one])
            assert p[one] == expected * (p[held] + p[one])
        for view in (
            lambda x: (datum("source", x), x & 7),
            lambda x: (datum("source", x), x >> 1),
        ):
            for v in set(map(view, range(16))):
                members = [x for x in range(16) if view(x) == v]
                assert all(
                    p[x] * sum(r[y] for y in members)
                    == r[x] * sum(p[y] for y in members)
                    for x in members
                )


def test_forward_banks_and_discriminating_controls():
    d = json.loads((HERE / "gibbs.json").read_text())
    assert d["paths"] == [physical_path(i) for i in range(3)]
    assert all(
        check["entries"] == 4900 and check["native_max_abs_error"] < 2e-13
        for check in d["checks"].values()
    )
    assert d["fixed_matching_slots"] == 140
    assert d["controls"]["half_only_max_fibre_error"] > 0.29
    assert d["controls"]["disabled_max_fibre_error"] > 0.89
    assert d["controls"]["disabled_max_identity_error"] < 2e-13
    assert {(b["level"], b["tau_index"]) for b in d["banks"]} == set(
        product(range(3), range(7))
    )
    for b in d["banks"]:
        s = sin(EXPOSURES[b["tau_index"]] / 2) ** 2
        assert abs(b["calibration_ones"] / N - s) < RADIUS
        assert {(u["even"] & 3).bit_count() for u in b["runs"]} == {0, 2}
        for u in b["runs"]:
            counts = u["word_counts"]
            assert sum(counts) == N
            assert (
                sum(n for x, n in enumerate(counts) if datum("source", x))
                == u["odd_starts"]
            )
            target = [
                comb(4, 4 - x.bit_count())
                * (s / 40 if datum("source", x) else (1 - s) / 30)
                for x in range(16)
            ]
            assert max(abs(np.array(counts) / N - target)) < RADIUS
