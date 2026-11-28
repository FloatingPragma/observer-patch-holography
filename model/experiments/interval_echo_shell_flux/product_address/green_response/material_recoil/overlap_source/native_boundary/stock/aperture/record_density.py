"""Source abundance and recovery from finite records of a four-owner repair service."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import product
import json
from math import comb, log

import numpy as np
from scipy.sparse.csgraph import connected_components

from . import pair_contact as contact

HERE = contact.HERE
stock = contact.stock
MENUS = ("none", "third_head", "ordinary", "all")
LEVELS = (Q(1, 16), Q(1, 4), Q(1, 2))


def classes():
    write, exchange, _, _ = stock.kernels()
    _, seams, _, edges, powers, index, words, codes = stock.carrier()
    graph = stock.graph(write, edges, powers, index, codes) + stock.graph(
        exchange, edges, powers, index, codes
    )
    _, labels = connected_components(graph, directed=False)
    source = []
    for k in range(5):
        starts = (
            (words[:, 4:] == 0).all(1)
            & (words[:, :4] != 1).all(1)
            & ((words[:, :4] == 2).sum(1) == k)
        )
        assert starts.sum() == comb(4, k)
        unique = np.unique(labels[starts])
        assert len(unique) == 1
        source.append(words[labels == unique[0]])
    return seams[:2], source


def counts(states):
    # Preserve all 36 accepted fine labels until the magnitude sufficient statistic is checked.
    fine = Counter()
    tables = {menu: defaultdict(lambda: [0] * 4) for menu in MENUS}
    fine_labels = np.flatnonzero((np.arange(11) != 0) & (stock.MAG != 1))
    for word in states:
        a, b, c = map(int, word[4:])
        if a not in fine_labels or b not in fine_labels:
            continue
        x, y = int(stock.MAG[a] == 2), int(stock.MAG[b] == 2)
        fine[a, b] += 1
        ordinary = tuple(map(int, word[:4]))
        for menu, r in zip(MENUS, [(), (c,), ordinary, (*ordinary, c)], strict=True):
            tables[menu][r][2 * x + y] += 1
    total = np.array(tables["none"][()])
    bits = stock.MAG[fine_labels] // 2
    multiplicity = np.bincount(bits, minlength=2)
    for i, j in product(range(6), repeat=2):
        a, b = int(bits[i]), int(bits[j])
        assert (
            fine[int(fine_labels[i]), int(fine_labels[j])]
            * multiplicity[a]
            * multiplicity[b]
            == total[2 * a + b]
        )
    return {menu: dict(table) for menu, table in tables.items()}


def input_weights(p):
    return [comb(4, k) * p**k * (1 - p) ** (4 - k) for k in range(5)]


def combine(panel, weights, menu):
    result = defaultdict(lambda: [Q(0)] * 4)
    for (size, tables), weight in zip(panel, weights, strict=True):
        if not weight:
            continue
        for r, values in tables[menu].items():
            target = result[r]
            for i, x in enumerate(values):
                target[i] += weight * x / size
    return dict(result)


def recover(tables):
    original = [Q(0)] * 4
    after = [Q(0)] * 4
    outside = Q(0)
    for row in tables.values():
        n = sum(row)
        assert n > 0
        a, b, c, d = row
        product_row = [u * v / n for u in (a + b, c + d) for v in (a + c, b + d)]
        for j in range(4):
            original[j] += row[j]
            after[j] += product_row[j]
            if row[j] == 0:
                outside += product_row[j]
    difference = [y - x for x, y in zip(original, after, strict=True)]
    delta = difference[3]
    assert difference == [delta, -delta, -delta, delta]
    # The comparison receiver sees A and B, with the witness omitted on both arms.
    a, b, c, d = original
    if min(original) > 0:
        charge = -float(delta) * log(a * d / (b * c))
    elif delta == 0:
        charge = 0.0
    else:
        charge = None
    return (
        original,
        after,
        dict(
            accepted=str(sum(original)),
            delta_both=str(delta),
            receiver_nats_per_attempt=charge,
            weight_into_witness_reference_null=str(outside),
        ),
    )


def dilute_terms(panel, menu):
    """Leading term of each witness determinant divided by its acceptance."""
    strata = set().union(*(table[menu] for _, table in panel))
    leading = defaultdict(Q)
    for r in strata:
        # Expand each cell in p exactly; the binomial inventory has degree four.
        cells = [[Q(0)] * 5 for _ in range(4)]
        for k, (n, tables) in enumerate(panel):
            values = tables[menu].get(r, [0] * 4)
            for j in range(5 - k):
                weight = Q(comb(4, k) * comb(4 - k, j) * (-1) ** j, n)
                for cell, value in zip(cells, values, strict=True):
                    cell[k + j] += weight * value
        a, b, c, d = cells
        numerator = [Q(0)] * 9
        for i, j in product(range(5), repeat=2):
            numerator[i + j] += b[i] * c[j] - a[i] * d[j]
        if any(numerator):
            denominator = [sum(t) for t in zip(*cells, strict=True)]
            top = next(i for i, v in enumerate(numerator) if v)
            bottom = next(i for i, v in enumerate(denominator) if v)
            leading[top - bottom] += numerator[top] / denominator[bottom]
    order = min(leading)
    # Higher-order stratum terms cannot affect this first nonzero order.
    assert leading[order] != 0
    return dict(order=order, coefficient=str(leading[order]))


def polynomial(panel):
    import sympy as sp

    p = sp.Symbol("p")
    base = [
        sum(
            w * Q(t["none"][()][i], n)
            for (n, t), w in zip(panel, input_weights(p), strict=True)
        )
        for i in range(4)
    ]
    x = sp.Symbol("x")

    def symmetric(expr):
        # x=p(1-p); power-reduce the exact polynomial modulo p²-p+x.
        reduced = sp.rem(sp.expand(expr), p * p - p + x, p)
        assert p not in reduced.free_symbols
        return sp.factor(reduced)

    rows = [symmetric(expr) for expr in base]
    total = sum(rows)
    delta = sp.factor((rows[1] * rows[2] - rows[0] * rows[3]) / total)
    num, den = sp.fraction(delta)

    def bernstein(expr):
        poly = sp.Poly(expr, x)
        # Bernstein coefficients on 0 <= x <= 1/4 enclose the polynomial range.
        d = poly.degree()
        return [
            sum(
                poly.nth(j) * sp.Rational(comb(i, j), comb(d, j) * 4**j)
                for j in range(i + 1)
            )
            for i in range(d + 1)
        ]

    order = sp.Poly(num, x).terms()[-1][0][0]
    assert order >= 0
    numerator = bernstein(sp.cancel(num / x**order))
    denominator = bernstein(den)
    # Return isolated sign changes too if coefficients do not have a common sign.
    roots = sp.polys.polytools.intervals(
        sp.Poly(num / x**order, x), eps=sp.Rational(1, 10**13)
    )
    roots = [
        (str(a), str(b)) for (a, b), _ in roots if b >= 0 and a <= sp.Rational(1, 4)
    ]
    return dict(
        variable="x=p(1-p)",
        joint_rows=[str(q) for q in rows],
        delta_both=str(delta),
        dilute_order=order,
        dilute_coefficient=str(sp.limit(delta / x**order, x, 0)),
        odds_at_dilution=str(sp.limit(rows[0] * rows[3] / (rows[1] * rows[2]), x, 0)),
        numerator_bernstein=list(map(str, numerator)),
        denominator_bernstein=list(map(str, denominator)),
        interior_zeros=roots,
    )


def held(old, new, rng):
    attempts = 2**42

    def bank(row):
        p = np.array([*map(float, row), float(1 - sum(row))])
        return rng.multinomial(attempts, p / p.sum()).tolist()

    first, last = bank(old), bank(new)
    e = np.sqrt(contact.RISK / (2 * attempts))

    def prob(count):
        return count / attempts

    d = prob(last[3] - first[3])
    band = d + np.array([-2, 2]) * e
    p = np.array(first[:4]) / attempts
    lo, hi = p - e, p + e
    assert min(lo) > 0
    odds = np.log([lo[0] * lo[3] / (hi[1] * hi[2]), hi[0] * hi[3] / (lo[1] * lo[2])])
    box = -np.outer(band, odds)
    assert band[0] <= float(new[3] - old[3]) <= band[1]
    return dict(
        attempts=attempts,
        native=first,
        recovered=last,
        delta_both=band.tolist(),
        reference_log_odds=odds.tolist(),
        source_nats_per_attempt=[float(box.min()), float(box.max())],
    )


def run():
    seams, source = classes()
    panel = [(len(states), counts(states)) for states in source]
    _, table, rejected, size = contact.source_counts()
    assert panel[2][0] == size
    assert panel[2][1]["none"][()] == table.ravel().tolist()
    rng = np.random.default_rng(621074)
    archive = {}
    for name, weights in [("fixed_2200", [Q(i == 2) for i in range(5)])] + [
        (str(p), input_weights(p)) for p in LEVELS
    ]:
        family = {}
        for menu in MENUS:
            old, new, row = recover(combine(panel, weights, menu))
            pooled = combine(panel, weights, "none")[()]
            assert old == pooled
            row["changed_passive_receiver"] = str(
                sum(abs(a - b) for a, b in zip(old, pooled, strict=True)) / 2
            )
            if name == str(LEVELS[0]):
                row["counts"] = held(old, new, rng)
                lo, hi = row["counts"]["source_nats_per_attempt"]
                assert lo <= row["receiver_nats_per_attempt"] <= hi
            family[menu] = row
        archive[name] = family
    return dict(
        sha256={
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in (
                HERE / "RECORD_DENSITY.md",
                HERE / "pair_contact.json",
                stock.HERE / "results.json",
            )
        },
        seams=seams,
        classes=[dict(words=n, magnitude_counts=t["none"][()]) for n, t in panel],
        third_magnitude_two_words=[
            int(np.count_nonzero(stock.MAG[states[:, 6]] == 2)) for states in source
        ],
        third_head_counts=[
            {str(r[0]): row for r, row in t["third_head"].items()} for _, t in panel
        ],
        dilute={menu: dilute_terms(panel, menu) for menu in MENUS},
        polynomial=polynomial(panel),
        cases=archive,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    path = HERE / "record_density.json"
    if args.write:
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        contact.timing.lab.lab.cold.apparatus.check(json.loads(path.read_bytes()), out)
    print(out["polynomial"])
    print(
        {
            n: {
                m: (float(Q(c["delta_both"])), c["receiver_nats_per_attempt"])
                for m, c in panel.items()
            }
            for n, panel in out["cases"].items()
        }
    )
