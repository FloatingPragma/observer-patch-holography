"""Exact incidence and material marginals for equality of two native register tallies."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
import json
import math

import numpy as np

from . import peer_fibers as fibers

HERE = fibers.source.tidal.HERE


def histogram(k, collapse=True, draw_bank=None):
    words, cells, groups, cat, types, counts, protected, unit, B = fibers.apparatus()
    base = k * counts // unit
    hold = k * protected // unit
    nfree = int(base.sum())
    size = int(hold.sum() + nfree)
    companions = {}
    by_cell = np.zeros(4, int)
    for g in np.unique(groups[cat < 0]):
        mask = groups == g
        assert np.all(cat[mask] < 0) and len(set(cells[mask])) == 1
        quota, remainder = divmod(k * int(mask.sum()), unit)
        assert remainder == 0
        companions[str(g)] = quota
        by_cell[cells[mask][0]] += quota
    assert np.array_equal(by_cell, hold)
    codes = np.array([s[0] for s in types])
    ref_table = hold + np.bincount(codes, weights=base, minlength=4).astype(int)
    r0 = int(ref_table[:2].sum())
    c0 = int(ref_table[::2].sum())
    total_shuffles = math.comb(size, c0)
    # The binomial peer count follows by placing c0 right-zero reports
    # uniformly into the two fixed left blocks. No target table is supplied.
    shuffles = {
        x: math.comb(r0, x) * math.comb(size - r0, c0 - x)
        for x in range(max(0, c0 + r0 - size), min(c0, r0) + 1)
    }
    assert sum(shuffles.values()) == total_shuffles
    multiplicities = counts // unit
    batches = fibers.batches(B, base, multiplicities, types, collapse)
    banks = defaultdict(lambda: [0, np.zeros(len(types), object)])
    total = 0
    for v, number, sums, count in batches:
        x = int(v[codes == 0].sum() + hold[0])
        if draw_bank is not None:
            draw_bank.append((number * shuffles[x], v.copy()))
        row = banks[x]
        row[0] += number
        row[1] += sums
        total += count
    denom = sum(row[0] for row in banks.values())
    numerator = sum(shuffles[x] * row[0] for x, row in banks.items())
    probability = Q(numerator, denom * total_shuffles)
    uniform = fibers.source.records()[1]["original"]
    out = {}
    ensembles = {}
    for tag, weights, norm in (
        ("ordinary_only", {x: 1 for x in banks}, denom),
        ("agreed", shuffles, numerator),
        (
            "disagreed",
            {x: total_shuffles - p for x, p in shuffles.items()},
            denom * total_shuffles - numerator,
        ),
    ):
        sums = sum(weights[x] * row[1] for x, row in banks.items())
        z = uniform.copy()
        for i, (amount, number) in enumerate(zip(sums, counts, strict=True)):
            z[cat == i] = Q(amount, norm * size * int(number))
        assert sum(z) == 1 and np.all(z[cat < 0] == uniform[cat < 0])
        ensembles[tag] = z
        out[tag] = dict(
            blank_change=str((z - uniform) @ (words[:, -1] == 0)),
            joint=[str(z @ (cells == i)) for i in range(4)],
            signed_variance=str(z @ (words[:, 1] - 1) ** 2),
            type_probabilities=[str(Q(t, norm * size)) for t in sums],
        )
    assert np.all(
        ensembles["agreed"] * probability + ensembles["disagreed"] * (1 - probability)
        == ensembles["ordinary_only"]
    )
    p_guard = Q(denom, int(multiplicities.sum()) ** nfree)
    info = dict(
        size=size,
        free_draws=nfree,
        protected_by_cell=hold.tolist(),
        protected_by_component=companions,
        histograms=total,
        ordinary_gate=str(p_guard),
        agreement_given_ordinary=str(probability),
        acceptance_given_free_register=str(p_guard * probability),
        material=out,
        head_tables={
            str(x): dict(
                ordinary_probability=str(Q(row[0], denom)),
                peer_agreement=str(Q(shuffles[x], total_shuffles)),
            )
            for x, row in sorted(banks.items())
        },
    )
    return info, ensembles


def run():
    result = {}
    for k in range(1, 17):
        panel, _ = histogram(k)
        if k > 4:
            del panel["head_tables"]
            for arm in panel["material"].values():
                del arm["type_probabilities"]
        result[str(panel["size"])] = panel
    return dict(
        sha256={
            "PEER_RETURN.md": sha256((HERE / "PEER_RETURN.md").read_bytes()).hexdigest()
        },
        registers=result,
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    out = run()
    raw = json.dumps(out, indent=2, allow_nan=False) + "\n"
    dst = HERE / "peer_population.json"
    if args.write:
        dst.write_text(raw)
    else:
        assert dst.read_text() == raw
    for key, row in out["registers"].items():
        print(
            key,
            row["histograms"],
            float(Q(row["ordinary_gate"])),
            float(Q(row["agreement_given_ordinary"])),
            [(k, float(Q(v["blank_change"]))) for k, v in row["material"].items()],
        )
