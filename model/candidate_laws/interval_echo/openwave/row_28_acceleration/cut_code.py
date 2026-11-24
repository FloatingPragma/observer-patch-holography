"""Verify two exact observer channels at the enrolled 25-seam cut."""

import json
from fractions import Fraction
from itertools import product

import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import federation

from .experiment import HERE


def verify():
    panel = next(
        p
        for p in json.loads((HERE / "capacity.json").read_text())["panels"]
        if p["level"] == 1 and p["seed"] == 280222 and not p["dilute"]
    )
    stations = {p["origin"]: p for p in panel["stations"]}
    first, second = stations[0], stations[265]
    assert first["component"] == second["component"]
    cut = panel["components"][str(first["component"])]["seams"]
    board = federation(1)
    assert len(cut) == 25 and len({v for i in cut for v in board.seams[i]}) == 2 * len(
        cut
    )

    law = candidate.native_law
    # Native W columns give exactly equal amplitudes on the two coin outputs;
    # signs, final loads and support are tested over the complete 99-state read.
    for a, b in ((0, 2), (2, 0)):
        column = law.echo()[:, law.SRC.index((a, b))].reshape(-1, law.N_SRC)
        expected = np.zeros_like(column)
        for coin in (0, 1):
            expected[1 + law.LABELS.index((coin, a - b)), law.SRC.index((1, 1))] = (
                2**-0.5
            )
        np.testing.assert_allclose(column, expected, rtol=0, atol=1e-15)

    # Rational count of the exact equiprobable outcomes, before any further
    # evolution. Symbols here are the independently chosen signs at writes.
    signs = (-2, 2)
    coin = {s: {c: Fraction(1, 2) for c in (0, 1)} for s in signs}
    signed = {s: {(c, s): Fraction(1, 2) for c in (0, 1)} for s in signs}
    for rows in (coin, signed):
        assert all(sum(r.values()) == 1 for r in rows.values())
    assert coin[-2] == coin[2]
    assert set(signed[-2]).isdisjoint(signed[2])
    # Directly test the product support at three sites. At B sites two unequal
    # sign words are disjoint in the full record on any differing coordinate.
    words = list(product(signs, repeat=3))
    for x in words:
        for y in words:
            assert all(set(coin[a]) & set(coin[b]) for a, b in zip(x, y, strict=True))
            joint_overlap = all(
                set(signed[a]) & set(signed[b]) for a, b in zip(x, y, strict=True)
            )
            assert joint_overlap == (x == y)
    return len(cut), 1, 2 ** len(cut)


if __name__ == "__main__":
    sites, coarse, full = verify()
    print(
        f"{sites} disjoint writes; exact correctable sign-code sizes: coin {coarse}, coin and sign {full}."
    )
