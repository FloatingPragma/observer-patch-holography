"""Transfer a single calibrated source-reference read to held flight appointments."""

import argparse
import json
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np

from .. import eraser
from ..joint_receivers.apparatus import access, pointer
from ..joint_receivers.measure import STARTS, bins, sample
from ..joint_receivers.readback import interval
from ..measure import OWNERS, sector_fields

HERE = Path(__file__).resolve().parent
TRAINING = HERE.parent / "joint_receivers" / "blind_counts.json"
TIMES = (384, 768)
PANELS = (-1, 1)
MODES = ("native", "which_head", "contrast_four")
SEED = 26112673


def acquire(times=TIMES, seed=SEED):
    access()
    sources = eraser.launch()
    banks = []
    preserve = np.equal.outer(np.arange(4) % 2, np.arange(4) % 2)
    for strength, modes in [(8, MODES[:2]), (4, MODES[2:])]:
        waves = sector_fields(strength, False, times=times)
        for owner, (t, time), recipe, panel in product(
            OWNERS, enumerate(times), range(3), PANELS
        ):
            field = waves[owner][:, t, :]
            coincidence = sources[recipe] * (field @ field.conj().T)
            at_source = np.trace(sources[recipe]).real
            pair = np.trace(coincidence).real
            for mode in modes:
                joint = coincidence * preserve if mode == "which_head" else coincidence
                outcomes = bins(pointer(panel), joint)
                np.testing.assert_allclose(sum(outcomes), pair, atol=1e-25, rtol=1e-12)
                key = [mode, owner, time, recipe, panel]
                counts = sample(
                    key,
                    [*outcomes, at_source - pair, 0, 1 / 3 - at_source, 2 / 3],
                    STARTS,
                    seed=seed,
                )
                banks.append(dict(key=key, counts=counts))
    return dict(starts=STARTS, seed=seed, banks=banks)


def calibrate(training, recipes):
    selected = {
        (panel, source): row["counts"][:3]
        for row in training["banks"]
        for owner, time, source, panel, role in [row["key"]]
        if (owner, time, role) == (12, 256, "training") and source in recipes
    }
    decoder = {}
    for panel in PANELS:
        counts = [selected[panel, source] for source in recipes]
        if any(sum(row) == 0 for row in counts):
            raise ValueError("training cohort has no intercepted pairs")
        decoder[panel] = [
            max(
                recipes,
                key=lambda j: Fraction(
                    selected[panel, j][exit_id], sum(selected[panel, j])
                ),
            )
            for exit_id in range(3)
        ]
    return decoder


def score(banks, decoder, recipes, starts):
    gains, intervals, per_bank = [], [], []
    for recipe, panel in product(recipes, PANELS):
        counts = banks[recipe, panel]
        n = sum(counts[:3])
        k = sum(
            count
            for count, guess in zip(counts[:3], decoder[panel], strict=True)
            if guess == recipe
        )
        gains.append(Fraction(len(recipes) * k - n, len(recipes) * starts))
        per_bank.append(dict(recipe=recipe, panel=panel, correct=k, pairs=n))
        if n:
            intervals.append(interval(k, n))
    complete = len(intervals) == len(per_bank)
    return dict(
        cohorts=per_bank,
        conditional_accuracy=(
            sum(row["correct"] / row["pairs"] for row in per_bank) / len(per_bank)
        )
        if complete
        else None,
        limits=np.mean(intervals, axis=0).tolist() if complete else None,
        pairs_per_trillion=sum(row["pairs"] for row in per_bank)
        * 1e12
        / (starts * len(per_bank)),
        extra_correct_per_trillion=float(sum(gains) / len(gains) * 10**12),
    )


def decode(data, training):
    labels = (0, 1, 2)
    signs = (1, 2)
    full, signed = calibrate(training, labels), calibrate(training, signs)
    keys = set(product(MODES, OWNERS, TIMES, labels, PANELS))
    bank = {tuple(r["key"]): r["counts"] for r in data["banks"]}
    if (
        set(bank) != keys
        or len(data["banks"]) != len(keys)
        or data["starts"] != STARTS
        or any(
            len(c) != 7
            or sum(c) != STARTS
            or any(type(n) is not int or n < 0 for n in c)
            for c in bank.values()
        )
    ):
        raise ValueError("incomplete or invalid original-start cohorts")
    rows = []
    for mode, owner, time in product(MODES, OWNERS, TIMES):
        b = {(j, p): bank[mode, owner, time, j, p] for j, p in product(labels, PANELS)}
        rows.append(
            dict(
                mode=mode,
                owner=owner,
                time=time,
                three=score(b, full, labels, STARTS),
                signs_with_three_labels=score(b, full, signs, STARTS),
                two=score(b, signed, signs, STARTS),
            )
        )
    return dict(three_label_map=full, two_label_map=signed, rows=rows)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--acquire", type=Path)
    args = cli.parse_args()
    if args.acquire:
        args.acquire.write_text(
            json.dumps(acquire(), separators=(",", ":"), allow_nan=False) + "\n"
        )
    else:
        data = json.loads((HERE / "counts.json").read_text())
        print(
            json.dumps(
                decode(data, json.loads(TRAINING.read_text())),
                indent=2,
                allow_nan=False,
            )
        )
