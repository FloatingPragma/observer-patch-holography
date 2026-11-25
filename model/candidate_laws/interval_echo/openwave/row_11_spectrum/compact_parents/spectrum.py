"""All finite-volume product energies for a compact single defect per occupied owner."""

from fractions import Fraction
from functools import cache
from math import ceil, isqrt, lcm
from pathlib import Path
import json

import numpy as np

from model.experiments.interval_echo_encounter import confinement as compact

HERE = Path(__file__).resolve().parent
TICKS = 10**12
ROOT = isqrt(5 * TICKS**2)
AUDIT = (
    HERE.parents[1] / "electron_field/lepton_family/arrival_spin/anchor_spectrum.json"
)


@cache
def modes(gates):
    _, h, exits = compact.blocks(gates, "single", 0)
    dark, polynomial = compact.restriction(h, exits)
    mult = lcm(*(v.denominator for v in dark.flat))
    w = np.array(dark * mult, dtype=int)
    assert not np.any(w[list(gates)])
    identity = np.eye(12, dtype=int)
    q = h @ h - 10 * h + 20 * identity
    base = q @ w
    column = np.flatnonzero(np.any(base, axis=0))[0]
    v = base[:, column]
    np.testing.assert_array_equal(h @ v, 6 * v)
    result = {0: (v, 0 * v)}
    w = (h - 6 * identity) @ w
    if w.any():
        v = w[:, np.flatnonzero(np.any(w, axis=0))[0]]
        u = (h - 5 * identity) @ v
        np.testing.assert_array_equal(h @ u, 5 * u + 5 * v)
        np.testing.assert_array_equal(h @ v, u + 5 * v)
        result[1] = (u, v)
        result[-1] = (u, -v)
    characteristic = [1, -6]
    expected = np.polymul(characteristic, characteristic)
    if len(result) == 3:
        expected = np.polymul(np.polymul(expected, characteristic), [1, -10, 20])
    np.testing.assert_array_equal(expected, polynomial)
    for a, b in result.values():
        assert not a[list(gates)].any() and not b[list(gates)].any()
    return result


def orb(gates, tag):
    a, b = modes(gates)[tag]
    z = a + np.sqrt(5) * b
    return z / np.linalg.norm(z)


def terms(n, special):
    """Each pair (c,d) identifies exactly one energy (c + d sqrt 5)/2."""
    assert 0 <= special < 204
    pairs = []
    exotic = []
    signed = []
    for p in range((n + 1) // 2):
        for h in range(min(2 * p + 1, special) + 1):
            d = np.arange(-h, h + 1, 2, dtype=np.int16)
            signed.append(d)
            pairs.append(np.full_like(d, p))
            exotic.append(np.full_like(d, h))
    p, h, d = (np.concatenate(x).astype(np.int64) for x in (pairs, exotic, signed))
    return p, h, d, 70 + 204 * p - h


def recipe(p, h, d):
    return dict(
        pairs=int(p),
        vacancies=int(p) + 1,
        plus_branch=int((h + d) // 2),
        minus_branch=int((h - d) // 2),
        ordinary_branch=int(2 * p + 1 - h),
    )


def scan(n, special, ground, targets):
    p, h, d, c = terms(n, special)
    # Each integer below encloses 2 E in units 1/TICKS to within special ticks.
    ticks = c * TICKS + d * ROOT
    assert ROOT**2 < 5 * TICKS**2 < (ROOT + 1) ** 2

    def entry(i):
        return dict(
            **recipe(p[i], h[i], d[i]),
            twice_energy=[int(c[i]), int(d[i])],
            energy=(int(c[i]) + int(d[i]) * np.sqrt(5)) / 2,
        )

    denominator = {
        "ordinary_compact": (35, 35),
        "ground": tuple(Fraction(v, 10**8) for v in ground),
    }
    if special:
        for name, sign in (("lower_compact", -1), ("upper_compact", 1)):
            denominator[name] = sorted(
                Fraction(69 * TICKS + sign * w, 2 * TICKS) for w in [ROOT, ROOT + 1]
            )
    matches = {}
    for light, (lo, hi) in denominator.items():
        matches[light] = {}
        value = float((lo + hi) / 2)
        for name, decimal in targets.items():
            ratio = Fraction(decimal)
            limits = [2 * TICKS * ratio * x for x in (lo, hi)]
            target = round(sum(limits) / 2)
            error = special + ceil(max(abs(a - target) for a in limits))
            distances = abs(ticks - target)
            i, j = np.argpartition(distances, 1)[:2]
            if distances[i] > distances[j]:
                i, j = j, i
            separation = int(distances[j]) - int(distances[i])
            assert separation > 2 * error
            best = entry(i)
            matches[light][name] = dict(
                **best,
                ratio_to_light=best["energy"] / value,
                ratio_offset=best["energy"] / value - float(ratio),
                certified_runner_up_separation_in_native_units=(separation - 2 * error)
                // 2
                / TICKS,
            )
    return dict(
        cells=n,
        exotic_cells=special,
        distinct_energies=len(c),
        maximum=entry(np.argmax(ticks)),
        comparisons=matches,
    )


def datasets():
    audit = json.loads(AUDIT.read_text())["scale_audit"]
    ground = json.loads((HERE.parent / "lightest/results.json").read_text())["boards"]
    targets = dict(
        measured_muon=str(audit["codata_eq_187_ratio"]),
        register_middle=audit["closed_form_ratios"][1],
        register_heavy=audit["closed_form_ratios"][2],
    )
    return ground, targets, audit["ratio_standard_uncertainty"]
