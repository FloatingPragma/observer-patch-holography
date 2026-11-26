"""Infer concentric marker volumes using disjoint forward, reverse and source-monitor cohorts."""

from fractions import Fraction
from itertools import combinations
from math import comb, log, sqrt
from pathlib import Path
import json
import numpy as np

from .. import clock_ball
from ..tracer_flux import measure as tracer
from . import flight, shutters

SHOTS = 2**34
TIMES = (Fraction(1, 8), Fraction(1, 4))
STOCKS = dict(
    zero=tracer.ZERO, one=tracer.ONE, split=clock_ball.SPLIT, spread=clock_ball.SPREAD
)
WORDS = sorted({w for words in STOCKS.values() for w in words})
EPS = sqrt(log(2 * 2000 / 0.01) / (2 * SHOTS)) + 1e-7


def band(n):
    return (max(0, n / SHOTS - EPS), min(1, n / SHOTS + EPS))


def difference(left, right):
    return (left[0] - right[1], left[1] - right[0])


def quotient(num, den):
    if den[0] <= 0 <= den[1]:
        raise ValueError("A registered gain or reciprocal count is unresolved")
    terms = [x / y for x in num for y in den]
    return (min(terms), max(terms))


def product(a, b):
    terms = [x * y for x in a for y in b]
    return (min(terms), max(terms))


def arm(radii, radius, states, starts, kernel, stock, j, backwards):
    inside = np.array([all(radii[t] <= radius for t in tags) for _, tags in states])
    rows = np.array([w in stock for w, _ in states])
    shares = rows.astype(float) / len(stock)
    n_sites = sum(d <= radius for d in radii)
    source_count = comb(n_sites if backwards else 12, j)
    keep = (inside[starts] if backwards else np.ones(len(starts), dtype=bool)) & rows[
        starts
    ]
    valid = keep.sum() / (source_count * len(stock))
    # Receiver independently draws w' from the same stock. Resident equality is
    # an ordinary binary event, implemented after the unbiased native flight.
    targets = np.ones(len(inside), dtype=bool) if backwards else inside
    land = kernel[:, keep].sum(1) / (source_count * len(stock))
    hit = float(land @ (shares * targets))
    match_miss = float(land @ (shares * ~targets))
    unchecked = float(land @ targets)
    categories = np.array(
        [hit, match_miss, unchecked - hit, valid - unchecked - match_miss, 1 - valid]
    )
    assert categories.min() > -2e-12
    categories = categories.clip(min=0)
    categories = categories / categories.sum()
    # Among authentic hits, count the exchanges where the terminal resident
    # word differs from the incoming one. Equality refers to the independent stock.
    changed = 0.0
    selected = starts[keep]
    for w in stock:
        columns = [
            c
            for c, source_ix in zip(np.flatnonzero(keep), selected, strict=True)
            if states[source_ix][0] == w
        ]
        others = np.array([state[0] != w for state in states])
        changed += (shares * targets * others) @ kernel[:, columns].sum(1)
    return categories, float(changed / (source_count * len(stock)))


def counts(rng, categories, changed):
    raw = rng.multinomial(SHOTS, categories).tolist()
    admitted = SHOTS - raw[-1]
    return dict(
        outcomes=raw,
        expected=categories.tolist(),
        matched_after_source_changed=changed,
        conditional_rate=categories[0] / (1 - categories[-1]),
        measured_conditional_interval=quotient(band(raw[0]), band(admitted)),
    )


def assays(rng, times, j):
    h, states, rad = flight.hardware(j)
    starts = np.array([i for i, (w, _) in enumerate(states) if w in WORDS])
    transitions, cert = flight.propagators(h, starts, times)
    assert max(abs(k.sum(0) - 1).max() for k in transitions) < 2e-12
    reader = None
    if j == 1:
        basis = [(tags[0], w) for w, tags in states]
        reader = tracer.certify_acquisition(
            basis, [i for i, r in enumerate(rad) if r == 2], rad.index(3)
        )
    all_times = []
    for time, kernel in zip(times, transitions, strict=True):
        readings = {}
        for name, stock in STOCKS.items():
            shells = []
            for radius in (1, 2):
                pair = [
                    counts(
                        rng,
                        *arm(rad, radius, states, starts, kernel, stock, j, backwards),
                    )
                    for backwards in (False, True)
                ]
                expected = sum(
                    comb(
                        sum(d <= radius for d in rad)
                        - sum(rad[t] <= radius for t in w),
                        j,
                    )
                    for w in stock
                ) / len(stock)
                unit = comb(10, j)
                # Reference admission is counted too: its number of placements
                # cancels when multiplying the conditional forward rate.
                reference_bound = quotient(
                    band(pair[0]["outcomes"][0]),
                    pair[1]["measured_conditional_interval"],
                )
                raw_reference = comb(12, j)
                lo, hi = (raw_reference * x for x in reference_bound)
                q = unit * pair[0]["conditional_rate"] / pair[1]["conditional_rate"]
                assert abs(q - expected) < 1e-10
                original = pair[0]["outcomes"]
                unchecked_bound = [
                    raw_reference * x for x in band(original[0] + original[2])
                ]
                blind = (
                    unit
                    * (pair[0]["expected"][0] + pair[0]["expected"][2])
                    / (1 - pair[0]["expected"][-1])
                )
                shells.append(
                    dict(
                        radius=radius,
                        forward=pair[0],
                        reverse=pair[1],
                        recovered_placements=q,
                        independently_counted_99=[lo, hi],
                        without_source_match=blind,
                        without_source_match_99=unchecked_bound,
                    )
                )
            readings[name] = shells
        changes = [
            difference(
                readings[b][0]["independently_counted_99"],
                readings[a][0]["independently_counted_99"],
            )
            for a, b in [("zero", "one"), ("split", "spread")]
        ]
        all_times.append(
            dict(
                time=float(time),
                arms=readings,
                ordinary_to_held_first_variation_nats=product(
                    (log(5), log(5)), quotient(changes[1], changes[0])
                ),
            )
        )
    return (
        dict(
            tags=j,
            states=len(states),
            reference_placements=comb(10, j),
            dynamics_bound=cert,
            times=all_times,
        ),
        reader,
        rad,
    )


def bulk(rng, rad):
    interior = [p for p, d in enumerate(rad) if d < 2]
    admitted = [p for p, d in enumerate(rad) if d <= 2]
    labels = [tuple(p for p in w if p in interior) for w in combinations(admitted, 2)]
    keys = sorted(set(labels))
    probs = np.array([labels.count(w) / len(labels) for w in keys])
    base = rng.multinomial(SHOTS, probs).tolist()
    scores = {
        w: tuple(log(x) for x in band(n)) for w, n in zip(keys, base, strict=True)
    }
    # The signed boundary capacity derivative equals minus the interior score
    # derivative at the counted equilibrium. Its additive log(55) cancels.
    arms = {}
    for name, stock in STOCKS.items():
        samples = rng.multinomial(SHOTS, [1 / len(stock)] * len(stock)).tolist()
        intervals = [
            product(band(n), scores[tuple(p for p in w if p in interior)])
            for w, n in zip(stock, samples, strict=True)
        ]
        arms[name] = dict(
            counts=samples,
            log_probability=[sum(v[k] for v in intervals) for k in (0, 1)],
        )
    return dict(
        uniform_interior_counts=base,
        labels=keys,
        held_arms=arms,
        measured_held_capacity_nats=difference(
            arms["spread"]["log_probability"], arms["split"]["log_probability"]
        ),
        measured_ordinary_capacity_nats=difference(
            arms["one"]["log_probability"], arms["zero"]["log_probability"]
        ),
    )


def run():
    rng = np.random.default_rng(803217)
    one, reader, rad = assays(rng, TIMES, 1)
    two, _, _ = assays(rng, TIMES, 2)
    interior = bulk(rng, rad)
    for bank in (one, two):
        for q in bank["times"]:
            arms = q["arms"]

            def rawchange(a, b):
                return difference(
                    arms[b][0]["independently_counted_99"],
                    arms[a][0]["independently_counted_99"],
                )

            q["entirely_counted_prediction_nats"] = product(
                interior["measured_ordinary_capacity_nats"],
                quotient(rawchange("split", "spread"), rawchange("zero", "one")),
            )

            def decoded(q):
                if bank["tags"] == 1:
                    return q
                return [(1 + sqrt(1 + 8 * x)) / 2 for x in q]

            shells = {
                name: difference(
                    decoded(s[1]["independently_counted_99"]),
                    decoded(s[0]["independently_counted_99"]),
                )
                for name, s in arms.items()
            }
            q["measured_shell_volume_per_step"] = shells
            q["finite_shell_secant_prediction_nats"] = product(
                interior["measured_ordinary_capacity_nats"],
                quotient(
                    difference(shells["spread"], shells["split"]),
                    difference(shells["one"], shells["zero"]),
                ),
            )
    return dict(
        shots_per_setting=SHOTS,
        simultaneous_coverage=0.99,
        probability_radius=EPS,
        dwell_hands=[tracer.material_hand(float(t), rng) for t in TIMES],
        shutters=shutters.run(rng, SHOTS, STOCKS),
        receiver=reader,
        clock_radii=rad,
        independent_interior=interior,
        tagged_volumeters=[one, two],
    )


if __name__ == "__main__":
    Path(__file__).with_name("results.json").write_text(
        json.dumps(run(), indent=2, allow_nan=False) + "\n"
    )
