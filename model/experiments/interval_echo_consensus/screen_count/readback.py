"""Independent component-charge prediction and simultaneous histogram checks."""

from argparse import ArgumentParser
import gzip
import json
from math import log, log2
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
MAX_CHECKS = 262144
LOG = log(2 * MAX_CHECKS / 0.01)
N = 262144


def components(size, pairs):
    groups = [{u} for u in range(size)]
    for a, b in pairs:
        left = next(g for g in groups if a in g)
        right = next(g for g in groups if b in g)
        if left is not right:
            left.update(right)
            groups.remove(right)
    return [sum(1 << i for i in g) for g in groups]


def summarize(archive):
    count = 0
    worst = 0

    def test(measured, probs, n):
        nonlocal count, worst
        p = np.array(probs, dtype=float)
        observed = np.array(measured)
        if (
            observed.dtype.kind not in "iu"
            or np.any(observed < 0)
            or observed.sum() != n
        ):
            raise ValueError("Incomplete comparison bank")
        if np.any((p == 0) & (observed != 0)) or p.sum() != 1:
            raise ValueError("Out-of-support comparison or wrong probability sum")
        allowance = np.sqrt(2 * n * p * (1 - p) * LOG) + 2 * LOG / 3
        rat = float(np.max(np.abs(observed - n * p) / allowance))
        if rat > 1:
            raise ValueError("Simultaneous count allowance exceeded")
        worst = max(worst, rat)
        count += len(p)

    levels = []
    for row in archive:
        c = row["config"]
        owners = c["owners"]
        n = len(owners)
        es = [(owners.index(u // 12), owners.index(v // 12)) for u, v in c["inside"]]
        site = owners.index(c["boundary"][0] // 12)
        w0 = 0
        for a, _ in es:
            w0 ^= 1 << a
        level = {}
        for mode, bank in row["panels"].items():
            q = {
                "parked": 0,
                "third": 0.25,
                "half": 0.5,
                "full": 1,
                "diagonal": 0,
                "cut": 0,
            }[mode]
            active = [
                e
                for i, e in enumerate(es)
                if mode != "diagonal" and not (mode == "cut" and i in c["skip"])
            ]
            groups = components(n, active)
            multiplicity = 2 ** (n - len(groups))
            words = np.arange(2**n)
            expected = np.zeros((2, len(words)))
            for exterior in (0, 1):
                shift = w0 ^ ((1 - exterior) << site)
                admissible = [
                    all(((int(w) ^ shift) & g).bit_count() % 2 == 0 for g in groups)
                    for w in words
                ]
                expected[exterior, admissible] = (
                    1 - q if exterior else q
                ) / multiplicity
            measured = np.array(bank["joint"]).reshape(expected.shape)
            if bank["refused"] != 0:
                raise ValueError("A primary station refused")
            test(measured.ravel(), expected.ravel(), N)
            marginal = measured.sum(axis=0)
            law = expected.sum(axis=0)
            for omitted in range(n):
                projection = (words & ((1 << omitted) - 1)) | (
                    (words >> (omitted + 1)) << omitted
                )
                accum = np.zeros(2 ** (n - 1), dtype=int)
                np.add.at(accum, projection, marginal)
                probs = np.zeros_like(accum, dtype=float)
                np.add.at(probs, projection, law)
                if len(groups) == 1 and np.any(probs != 1 / len(probs)):
                    raise ValueError("Missing product marginal")
                test(accum, probs, N)
            for g in groups:
                mask = np.array([(int(w) & g).bit_count() % 2 for w in words]).astype(
                    bool
                )
                test(
                    [sum(marginal[~mask]), sum(marginal[mask])],
                    [sum(law[~mask]), sum(law[mask])],
                    N,
                )
            odd = np.array([int(w).bit_count() % 2 for w in words], dtype=bool)
            level[mode] = dict(
                component_sizes=sorted(g.bit_count() for g in groups),
                support=int(np.count_nonzero(law)),
                observed=int(np.count_nonzero(marginal)),
                entropy_bits=round(-sum(p * log2(p) for p in law if p), 8),
                odd_reports=int(sum(marginal[odd])),
                external_occupied=int(measured[1].sum()),
                minimum_nonempty_bin=int(marginal[law > 0].min()),
            )
        for mode, counts in row["clocks"].items():
            q = {"parked": 0, "third": 0.25, "half": 0.5, "full": 1, "diagonal": 0}[
                mode
            ]
            test(counts, [1 - q, q], 65536)
        levels.append(level)
    assert count <= MAX_CHECKS
    return dict(
        checked_bins=count, maximum_fraction_of_allowance=round(worst, 8), levels=levels
    )


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    with gzip.open(HERE / "records.json.gz", "rt") as src:
        out = summarize(json.load(src))
    data = json.dumps(out, indent=2) + "\n"
    target = HERE / "summary.json"
    if args.write:
        target.write_text(data)
    elif target.read_text() != data:
        raise ValueError("Screen summary has changed")
