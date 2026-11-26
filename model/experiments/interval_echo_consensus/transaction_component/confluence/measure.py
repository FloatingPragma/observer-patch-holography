"""Complete finite rewrite census and held, routed overlapping critical pairs."""

import argparse
from collections import Counter
from functools import cache
import gzip
from itertools import combinations, product
import json
from pathlib import Path

from .apparatus import Experiment, MASKS, calibration, layout, mismatch
from ..contract import batch, launch, require, validate

OFFERS = tuple((side, m) for side in range(2) for m in MASKS)
WORDS = tuple(product((0, 1), repeat=3))


def move(b, x, o):
    side, mask = o
    y = list(x)
    y[side] = tuple(b[side] if j in mask else val for j, val in enumerate(x[side]))
    return tuple(y)


def graph():
    total = Counter()
    paths = Counter()
    fixed = set()
    rows = []
    for b in product((0, 1), repeat=2):

        @cache
        def end_paths(x):
            succ = [move(b, x, o) for o in OFFERS]
            n = sum(yy != x for yy in succ)
            if not n:
                return {x: Counter({0: 1})}
            ends = {}
            for yy in succ:
                if yy == x:
                    continue
                for v, counts in end_paths(yy).items():
                    add = ends.setdefault(v, Counter())
                    for length, count in counts.items():
                        add[length + 1] += count
            return ends

        for x in product(WORDS, repeat=2):
            rank = sum(mismatch(bb, w) for bb, w in zip(b, x, strict=True))
            targets = {}
            for i, o in enumerate(OFFERS):
                y = move(b, x, o)
                if y == x:
                    continue
                require(
                    sum(mismatch(bb, w) for bb, w in zip(b, y, strict=True)) < rank,
                    "Rank failed",
                )
                side, mask = o
                require(
                    all(y[side][j] == b[side] for j in mask),
                    "Covered mismatch survived",
                )
                targets[i] = y
            total["states"] += 1
            total["enabled_offers"] += len(targets)
            total["distinct_arrows"] += len(set(targets.values()))
            total["normal"] += int(not targets)
            total["stuttering_offers"] += len(OFFERS) - len(targets)
            for i, j in combinations(targets, 2):
                a, c = targets[i], targets[j]
                u = move(b, a, OFFERS[j])
                require(u == move(b, c, OFFERS[i]), "Unjoined one-step fork")
                total["peaks"] += 1
                if a == c:
                    total["coincident"] += 1
                elif u != a and u != c:
                    total["two_proper_legs"] += 1
                    if OFFERS[i][0] == OFFERS[j][0]:
                        total["overlapping_two_proper_legs"] += 1
                else:
                    total["one_zero_leg"] += 1
            endpoints = end_paths(x)
            expected = tuple((bb, bb, bb) for bb in b)
            require(set(endpoints) == {expected}, "Nonunique normal endpoint")
            fixed.add((b, expected))
            counts = endpoints[expected]
            paths.update(counts)
            rows.append(
                [
                    list(b),
                    x,
                    rank,
                    list(targets),
                    [counts.get(length, 0) for length in range(7)],
                ]
            )
        end_paths.cache_clear()
    return dict(
        total=dict(total),
        paths={str(k): v for k, v in sorted(paths.items())},
        fixed=sorted(fixed),
        rows=rows,
    )


def prior_program():
    """Compute the closure of the old aggregate operators if batch priority is relaxed."""
    starts = 0
    schedules = Counter()
    endpoints = set()
    for a, b, x, y in product((0, 1), (0, 1), WORDS, WORDS):
        if not mismatch(a, x) or not mismatch(b, y):
            continue
        starts += 1

        def visit(states, ids):
            enabled = batch(states)
            if not enabled:
                expected = ((a, 0, a, a, a, 1, 1), (b, 0, b, b, b, 1, 1))
                require(tuple(states) == expected, "Prior aggregate endpoint changed")
                endpoints.add(expected)
                schedules["".join(map(str, ids))] += 1
                return
            for tx in enabled:
                require(
                    validate(tx, states[tx.side]) == "ok", "Prior aggregate refused"
                )
                changed = list(states)
                changed[tx.side] = tx.result()
                visit(changed, (*ids, tx.identifier))

        visit(launch(a, b, x, y), ())
    return dict(
        sources=starts,
        terminal_count=len(endpoints),
        orders=dict(sorted(schedules.items())),
    )


def run():
    def seq(*settings):
        return [(s, m, False) for s, m in settings]

    one = seq((0, 0), (1, 1), (0, 2), (1, 3), (0, 3))
    other = seq((0, 2), (0, 0), (1, 3), (0, 3), (1, 1))
    source = ((1, 0), ((0, 0, 0), (1, 0, 1)))
    panels = {
        "left_first": Experiment(*source).exercise(one),
        "right_first": Experiment(*source).exercise(other),
        "different_data": Experiment((1, 0), ((0, 1, 0), (0, 1, 1))).exercise(one),
        "hidden_ones": Experiment(*source, hidden=1).exercise(one),
        "whole_first": Experiment(*source).exercise(seq((0, 3), (1, 3))),
        "parked": Experiment(*source).exercise([(0, 0, True)]),
        "coarse_zero": Experiment((0, 0), ((1, 1, 1), (1, 1, 1))).exercise(
            seq((1, 3), (0, 3))
        ),
        "coarse_one": Experiment((1, 1), ((0, 0, 0), (0, 0, 0))).exercise(
            seq((0, 3), (1, 3))
        ),
        "consistent": Experiment((1, 0), ((1, 1, 1), (0, 0, 0))).exercise(
            seq((0, 0), (1, 3))
        ),
    }
    for key in (
        "left_first",
        "right_first",
        "whole_first",
        "different_data",
        "hidden_ones",
    ):
        require(
            panels[key]["end"] == panels["whole_first"]["end"],
            "Held normal endpoints disagree",
        )
    for key, middle in (("left_first", [1, 1, 0]), ("right_first", [0, 1, 1])):
        hits = [
            r for r in panels[key]["observations"] if r[0] == "published" and r[1] == 0
        ]
        require(
            [r[-1] for r in hits] == [1, 0] and hits[0][-2] == middle,
            "Critical-pair alias",
        )
    z, v = panels["coarse_zero"], panels["coarse_one"]
    require(
        z["end"][1] == v["end"][1] and z["end"][0] != v["end"][0],
        "Parity fixture failed",
    )
    require(panels["parked"]["initial"] == panels["parked"]["end"], "Parked arm moved")
    return dict(
        schema=1,
        offers=OFFERS,
        allocation=[layout(i) for i in range(2)],
        calibration=calibration(),
        panels=panels,
        graph=graph(),
        prior_program=prior_program(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", action="store_true")
    args = parser.parse_args()
    result = run()
    serialized = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    archive = Path(__file__).with_name("records.json.gz")
    if args.acquire:
        with archive.open("xb") as f:
            f.write(gzip.compress(serialized, mtime=0))
    else:
        require(
            serialized == gzip.decompress(archive.read_bytes()),
            "Public receipt differs",
        )
    print(
        "Verified complete overlapping-repair graph, two native branch orders and protected parity control."
    )


if __name__ == "__main__":
    main()
