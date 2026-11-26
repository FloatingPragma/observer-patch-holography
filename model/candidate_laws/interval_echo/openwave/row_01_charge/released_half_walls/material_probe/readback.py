"""Score original-start words and conditional probes against independent material standards."""

import json
from itertools import product
from math import log, sqrt

import numpy as np

from .flight import ANALYZERS, BOARDS, DWELLS, mount, probability
from .. import ring
from .measure import HERE, LAWS, SEED, SHOTS


def enrolled(rows, expected, width):
    banks = {}
    for key, counts in rows:
        key = tuple(key)
        if (
            key in banks
            or len(counts) != width
            or any(type(x) is not int or x < 0 for x in counts)
            or sum(counts) != SHOTS
        ):
            raise ValueError("Invalid counts or duplicate bank")
        banks[key] = np.asarray(counts)
    if set(expected) != banks.keys():
        raise ValueError("Missing or unregistered banks")
    return banks


def fraction(k, n):
    if n == 0:
        return None
    p = k / n
    r = sqrt(log(2 * 32768 / 0.01) / (2 * n))
    return np.array([max(0, p - r), min(1, p + r)])


def overlap(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


def analyze(raw):
    if raw["seed"] != SEED or raw["shots"] != SHOTS:
        raise ValueError("Unregistered acquisition")
    enrolled_mounts = []
    for lvl, n in BOARDS:
        contacts, distant = mount(lvl, n)
        enrolled_mounts.append(
            dict(
                level=lvl,
                length=n,
                contacts=[list(v) for v in contacts],
                distant=[list(v) for v in distant],
                ring=list(ring.geometry(lvl, n)),
            )
        )
    if raw["mounts"] != enrolled_mounts:
        raise ValueError("Unregistered physical probe mount")
    standard_keys = [
        (*g, k, side, x, y, s, t, r)
        for g, k, side, x, y, s, t, r in product(
            BOARDS, (4, 8), range(2), range(3), range(3), (-1, 1), range(2), ANALYZERS
        )
    ]
    held_keys = [
        (*g, law, amount, source, probe, t, r)
        for g, law, amount, source, probe, t, r in product(
            BOARDS, LAWS, (1, 2), (-1, 1), (-1, 1), range(2), ANALYZERS
        )
    ]
    standards = enrolled(raw["standard"], standard_keys, 2)
    banks = enrolled(raw["measurements"], held_keys, 8)
    cal_fails = 0
    for (*_, k, side, x, y, s, t, r), v in standards.items():
        pred = probability(-2 * k * s * (x + y), DWELLS[t], r)
        cal_fails += int(not overlap(fraction(v[1], SHOTS), [pred, pred]))
    tally = {}
    examples = []
    for (lvl, n, law, amt, src, s, t, r), v in banks.items():
        _, k, distant = LAWS[law]
        yes = v.reshape(2, 2, 2)[1]
        for side, hits in enumerate((yes[1].sum(), yes[:, 1].sum())):
            obs = fraction(hits, yes.sum())
            stats = tally.setdefault(
                law,
                dict(
                    arms=0,
                    empty=0,
                    full_match=0,
                    mixture_excluded=0,
                    half_detuning_excluded=0,
                    zero_match=0,
                ),
            )
            stats["arms"] += 1
            if obs is None:
                stats["empty"] += 1
                continue
            if distant or k == 0:
                stats["zero_match"] += int(overlap(obs, [0.5, 0.5]))
                continue
            endpoint = amt if src > 0 else 0

            def cal(x, y):
                return fraction(standards[lvl, n, k, side, x, y, s, t, r][1], SHOTS)

            full = cal(endpoint, endpoint)
            base = cal(0, amt)
            stats["full_match"] += int(overlap(obs, full))
            if amt == 1:
                half_mixture = (full + base) / 2
            else:
                half_mixture = cal(1, endpoint)
            stats["mixture_excluded"] += int(not overlap(obs, half_mixture))
            pred = probability(-2 * k * s * amt * (1 + 0.5 * src), DWELLS[t], r)
            stats["half_detuning_excluded"] += int(not overlap(obs, [pred, pred]))
            if (lvl, n, law, side, t, r) == (
                0,
                14,
                "native",
                0,
                0,
                "phase",
            ) and amt == 1:
                examples.append(
                    dict(
                        source=src,
                        probe=s,
                        accepted=int(yes.sum()),
                        first_rail=int(hits),
                        fraction=hits / yes.sum(),
                        interval=obs.tolist(),
                        full_standard=full.tolist(),
                        half_mixture=half_mixture.tolist(),
                        half_detuning=pred,
                    )
                )
    return dict(
        sharp_calibrations=len(standards),
        calibration_failures=cal_fails,
        held_banks=len(banks),
        counts=tally,
        examples=examples,
    )


if __name__ == "__main__":
    print(json.dumps(analyze(json.loads((HERE / "counts.json").read_text())), indent=2))
