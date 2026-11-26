from itertools import combinations

import numpy as np
import pytest

from . import apparatus as lab, measure, readback


def direct_offer(source, probe, rate, time):
    """Propagate all 325 two-unit words of the 26 physical rails."""
    n = lab.RAILS
    states = list(combinations(range(2 * n), 2))
    labels = dict(zip(states, range(len(states)), strict=True))
    wave = np.zeros(len(states), complex)
    for a, s in enumerate(source):
        for b, p in enumerate(probe):
            wave[labels[a, n + b]] = s * p
    gate = lab.contact(rate, time)
    for r in range(n):
        out = np.zeros_like(wave)
        for j, word in enumerate(states):
            present = r in word
            other = r + n in word
            if present and other:
                out[j] += gate[3, 3] * wave[j]
            elif not (present or other):
                out[j] += gate[0, 0] * wave[j]
            else:
                here = r if present else r + n
                there = r + n if present else r
                row = 1 if present else 2
                rest = word[1] if word[0] == here else word[0]
                out[j] += gate[row, row] * wave[j]
                out[labels[tuple(sorted((there, rest)))]] += (
                    gate[3 - row, row] * wave[j]
                )
        wave = out
    both = [
        sum(
            abs(wave[j]) ** 2
            for j, pair in enumerate(states)
            if sum(x < n for x in pair) == k
        )
        for k in (0, 2)
    ]
    mixed = np.array([[wave[labels[a, n + b]] for b in range(n)] for a in range(n)])
    kept = source.conj() @ mixed
    return kept, [*both, np.vdot(mixed, mixed).real - np.vdot(kept, kept).real]


def test_independent_325_state_contact_and_winding_return():
    rng = np.random.default_rng(226)
    p = rng.normal(size=13) + 1j * rng.normal(size=13)
    p /= np.linalg.norm(p)
    for q, sign, rate, time in [
        (1, 1, 1, np.pi / 2),
        (3, -1, 0.5, np.pi / 2),
        (0, -1, 1, 0.73),
    ]:
        s = lab.source(q, sign, rate)
        target, flags = direct_offer(s, p, rate, time)
        found, counts = lab.offer(p[:, None], q, sign, rate, time)
        np.testing.assert_allclose(found[:, 0], target, atol=2e-13)
        np.testing.assert_allclose(counts[0], flags, atol=2e-13)
        np.testing.assert_allclose(
            np.vdot(target, target).real + sum(flags), 1, atol=2e-13
        )


def test_legal_hardware_contacts_and_signed_pointer_calibration():
    layout = lab.layout()
    rails = set(layout["sources"][0] + layout["sources"][1] + layout["probe"])
    assert len(rails) == 39 and layout["double_standard"] not in rails
    assert len(layout["contacts"]) == 26
    for rate in (1, 0.5):
        np.testing.assert_allclose(
            lab.standards(rate), [[0, 1], [1, 0]] * 2, atol=1e-13
        )


def test_bilinear_contrast_uses_failures_and_survives_every_orientation():
    obs = readback.ideal_geometry()
    for index, center in enumerate(measure.CENTERS):
        p = np.zeros((2, 2, 2, 4, len(measure.COORDS), 9))
        anti = lab.geometry()[1]
        for order, seq in enumerate(((0, 1), (1, 0))):
            for n, (s, t) in enumerate(((-1, -1), (-1, 1), (1, -1), (1, 1))):
                for j, (x, y) in enumerate(
                    (
                        (0, center),
                        (0, anti[center]),
                        (anti[0], center),
                        (anti[0], anti[center]),
                    )
                ):
                    p[order, n // 2, n % 2, j] = lab.banks(
                        seq, (x, y), (s, t), 1, np.pi / 2, 0, measure.COORDS
                    )
        observed = readback.order_contrast(p * measure.PAIR_SHOTS)
        np.testing.assert_allclose(observed, obs[:, 12 + index], atol=1e-15)


def test_archived_rank_bound_and_invalid_or_missing_acquisitions():
    with np.load(measure.HERE / "counts.npz") as file:
        data = dict(file)
    result = readback.decode(data)
    assert result["lower_99"] > 0 and result["total_simulated_starts"] == 2993592205312
    assert result["failed_encounters"] > result["total_simulated_starts"] / 2
    with pytest.raises(ValueError, match="Missing held"):
        readback.decode({k: v for k, v in data.items() if k != "held"})
    bad = data["single"].copy()
    bad[0, 0, 0, 0, 0] += 1
    with pytest.raises(ValueError, match="Incomplete single"):
        readback.decode(dict(data, single=bad))
    doubled = data["double"].copy()
    doubled[0, 0] = doubled[0, 1]
    with pytest.raises(ValueError, match="misses forecast"):
        readback.decode(dict(data, double=doubled))


def test_integer_cdf_near_large_conditional_quantiles():
    from fractions import Fraction
    from math import comb
    from .sampling import quantiles, resolving_cdf

    for n, w, mass in ((41, 1, 3), (80, 7, 15), (36, 8, 9)):
        p = Fraction(w, mass)
        for k in (0, n // 3, n - 1):
            exact = sum(comb(n, x) * p**x * (1 - p) ** (n - x) for x in range(k + 1))
            assert abs(float(exact) - resolving_cdf(k, n, w, mass)) < 2e-14
    n = np.array([402646717, 402641923])
    levels = np.array([0.2815121297271622, 0.16446213629296447])
    estimated = quantiles(levels, n, np.ones(2, dtype=np.int64), np.full(2, 3))
    np.testing.assert_array_equal(estimated, [134210101, 134204740])
