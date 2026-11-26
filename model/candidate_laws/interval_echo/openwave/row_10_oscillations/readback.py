"""Interpret original-start receiver counts and the six local ternary standards."""

import json
from itertools import product

from .measure import BASES, HERE, OWNERS, PANELS, SEED, SIZE, TIMES, WORDS
from .mixture import checked, interval


def incompatible(x, y):
    return x[0] > y[1] or y[0] > x[1]


def report(raw):
    n = raw["starts_per_bank"]
    if (
        type(n) is not int
        or n != SIZE
        or raw["seed"] != SEED
        or tuple(raw["source_bases"]) != BASES
        or tuple(raw["final_appointments"]) != TIMES
    ):
        raise ValueError("Unregistered acquisition size")
    banks = {tuple(row["key"]): checked(row["counts"], 9, n) for row in raw["banks"]}
    designs = list(
        product(
            ((8, False), (4, False), (8, True)), OWNERS, TIMES, BASES, range(3), PANELS
        )
    )
    expected = {(*law, o, t, b, s, p) for law, o, t, b, s, p in designs}
    if banks.keys() != expected or len(raw["banks"]) != len(expected):
        raise ValueError("Unregistered or missing flight banks")
    cal = {
        (r["panel"], r["word"]): checked(r["counts"], 6, n)
        for r in raw["local_word_standards"]
    }
    if (
        set(cal) != set(product(PANELS, range(6)))
        or len(raw["local_word_standards"]) != 18
    ):
        raise ValueError("Incomplete local standards")
    for panel, s in product(PANELS, range(6)):
        clicks = cal[panel, s]
        target = {s}
        if panel != "sharp":
            k = PANELS.index(panel) - 1
            moved = list(WORDS[s])
            moved[k], moved[k + 1] = moved[k + 1], moved[k]
            target.add(WORDS.index(tuple(moved)))
        for i, x in enumerate(clicks):
            if incompatible(interval(x, n), [float(i in target) / len(target)] * 2):
                raise ValueError("Native isolated contact failed a sharp calibration")
    coherent, changed, conjugates = [0, 0], [0, 0], [0, 0]
    no_arrivals = 0
    for law, o, t, b, s, p in designs:
        counts = banks[*law, o, t, b, s, p]
        if law[1]:
            no_arrivals += sum(counts[:6])
            continue
        if p != "sharp":
            q = PANELS.index(p) - 1
            sharp = banks[*law, o, t, b, s, "sharp"]
            for i, w in enumerate(WORDS):
                sw = list(w)
                sw[q], sw[q + 1] = sw[q + 1], sw[q]
                classical = interval(sharp[i] + sharp[WORDS.index(tuple(sw))], n) / 2
                coherent[0] += int(incompatible(classical, interval(counts[i], n)))
                coherent[1] += 1
        if law[0] == 8:
            alt = banks[4, False, o, t, b, s, p]
            changed[0] += int(
                incompatible(interval(sum(counts[:6]), n), interval(sum(alt[:6]), n))
            )
            changed[1] += 1
        if b == "character" and s == 1:
            partner = banks[*law, o, t, b, 2, p]
            for x, y in zip(counts, partner, strict=True):
                conjugates[0] += int(incompatible(interval(x, n), interval(y, n)))
                conjugates[1] += 1
    if no_arrivals:
        raise ValueError("Pair reached a graph-disconnected receiver")
    examples = []
    for s in range(3):
        row = banks[8, False, 17, 2048, "arm", s, "sharp"]
        pair, k = sum(row[:6]), row[WORDS.index((2, 1, 0))]
        examples.append(
            dict(
                arm=s,
                pairs=pair,
                word_210=k,
                conditional_interval=interval(k, pair).tolist(),
            )
        )
    return dict(
        arm_standards=examples,
        blocked_clicks=no_arrivals,
        total_banks=len(banks),
        resolved_against_incoherent_half_swap=coherent,
        resolved_free_contrast_pair_rates=changed,
        incompatible_conjugate_bin_pairs=conjugates,
    )


if __name__ == "__main__":
    print(json.dumps(report(json.loads((HERE / "records.json").read_text())), indent=2))
