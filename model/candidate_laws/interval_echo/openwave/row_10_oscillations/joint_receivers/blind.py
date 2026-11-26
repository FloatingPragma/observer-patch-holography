"""Train a material click classifier and score independent three-source cohorts."""

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np

from .. import eraser
from ..measure import OWNERS, TIMES, sector_fields
from .apparatus import access, pointer
from .measure import STARTS, bins, sample
from .readback import interval

SEED = 26112642
PANELS = (-1, 1)


def acquire():
    access()
    source = eraser.launch()
    waves = sector_fields(8, False)
    rows = []
    for owner, (i, time), sign, panel in product(
        OWNERS, enumerate(TIMES), range(3), PANELS
    ):
        field = waves[owner][:, i, :]
        joint = source[sign] * (field @ field.conj().T)
        far = np.trace(joint).real
        stock = np.trace(source[sign]).real
        local = bins(pointer(panel), joint)
        p = [*local, stock - far, 0, 1 / 3 - stock, 2 / 3]
        for role in ("training", "validation"):
            key = [owner, time, sign, panel, role]
            rows.append(dict(key=key, counts=sample(key, p, STARTS, seed=SEED)))
    return dict(seed=SEED, starts=STARTS, banks=rows)


def decode(record):
    b = {tuple(r["key"]): r["counts"] for r in record["banks"]}
    if len(b) != len(OWNERS) * len(TIMES) * 3 * len(PANELS) * 2 or not all(
        sum(v) == record["starts"] and min(v) >= 0 for v in b.values()
    ):
        raise ValueError("incomplete three-recipe record")
    rows = []
    for owner, time in product(OWNERS, TIMES):
        labels = []
        conditional = np.zeros((3, 3))
        correct = []
        pair_counts = []
        separated = set()
        for panel in PANELS:
            train = np.array(
                [b[owner, time, j, panel, "training"][:3] for j in range(3)]
            )
            val = np.array(
                [b[owner, time, j, panel, "validation"][:3] for j in range(3)]
            )
            guess = (train / train.sum(axis=1)[:, None]).argmax(axis=0)
            labels.append(guess.tolist())
            for left, right in combinations(range(3), 2):
                for outcome in range(3):
                    a, bnd = [
                        interval(int(val[s, outcome]), int(val[s].sum()))
                        for s in (left, right)
                    ]
                    if a[0] > bnd[1] or bnd[0] > a[1]:
                        separated.add((left, right))
            for j in range(3):
                n = int(val[j].sum())
                pair_counts.append(n)
                grouped = np.bincount(guess, weights=val[j], minlength=3)
                k = int(grouped[j])
                correct.append([k, n])
                conditional[j] += grouped / (2 * n)
        np.testing.assert_allclose(conditional.sum(axis=1), 1, atol=1e-15)
        limits = np.mean([interval(k, n) for k, n in correct], axis=0)
        rows.append(
            dict(
                owner=owner,
                time=time,
                decision=labels,
                conditional_confusion=conditional.tolist(),
                panel_source_correct=correct,
                validation_far_counts=pair_counts,
                balanced_success=np.trace(conditional) / 3,
                limits=limits.tolist(),
                raw_pairwise_separations=sorted(separated),
                label_count=len(set(sum(labels, []))),
                equally_likely_given_pair_herald_only=1 / 3,
            )
        )
    return rows


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--acquire", type=Path)
    args = p.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        print(
            json.dumps(
                decode(
                    json.loads(
                        Path(__file__).with_name("blind_counts.json").read_text()
                    )
                ),
                indent=2,
            )
        )
